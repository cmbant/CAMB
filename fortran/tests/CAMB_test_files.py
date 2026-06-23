import argparse
import bisect
import copy
import filecmp
import fnmatch
import math
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

import numpy as np

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
FORTRAN_DIR = os.path.abspath(os.path.join(TESTS_DIR, ".."))
REPO_ROOT = os.path.abspath(os.path.join(FORTRAN_DIR, ".."))
DEFAULT_BASE_SETTINGS = os.path.join(REPO_ROOT, "inifiles", "params.ini")
DEFAULT_TOLERANCE_SCALE = 1.0
LPA_AUTO_LMAX_SETTINGS = {
    4000: ["l_max_scalar = 4200", "k_eta_max_scalar = 90000", "lens_output_margin = 200", "do_nonlinear = 2"],
    6000: ["l_max_scalar = 6200", "k_eta_max_scalar = 162000", "lens_output_margin = 200", "do_nonlinear = 2"],
}
THETA_CASES = {
    "thetastar_lmax4000_lpa_auto": ("thetastar", "0.01039"),
    "cosmomc_theta_lmax4000_lpa_auto": ("cosmomc_theta", "0.01043"),
    "cosmomc_theta_ppf_w-0p82_wa0p35_lmax4000": ("cosmomc_theta", "0.0104085"),
    "thetastar_ppf_w-1p2_wa0p5_lmax4000": ("thetastar", "0.0104085"),
    "thetastar_ppf_w-0p7_wa-0p3_lmax4000": ("thetastar", "0.0104085"),
}
SUPPRESSED_CMB_OUTPUT_SETTINGS = ["scalar_output_file =", "total_output_file ="]
TRANSFER_FILE_PATTERNS = [
    "*transfer_out.dat",
    "*transfer_out2.dat",
    "*transfer_out_z0.dat",
    "*transfer_2.dat",
    "*transfer_3.dat",
]
TRANSFER_CHECKED_COLUMNS = ["CDM", "baryon", "v_CDM", "v_b", "Weyl"]
TRANSFER_POWER_VAR_COLUMNS = {
    2: "CDM",
    3: "baryon",
    7: "total",
    8: "no_nu",
    9: "total_de",
    10: "Weyl",
    11: "v_CDM",
    12: "v_b",
}


def parse_override_assignment(text):
    if "=" not in text:
        raise argparse.ArgumentTypeError("Overrides must use KEY=VALUE syntax")
    key, value = text.split("=", 1)
    key = key.strip()
    if not key:
        raise argparse.ArgumentTypeError("Override key cannot be empty")
    return key, value.strip()


def build_parser():
    parser = argparse.ArgumentParser(description="Run CAMB tests")
    parser.add_argument("ini_dir", help="ini file directory")
    parser.add_argument("--make_ini", "--make-ini", action="store_true", help="If set, output ini files to ini_dir")
    parser.add_argument("--out_files_dir", "--out-files-dir", default="test_outputs", help="output files directory")
    parser.add_argument(
        "--base_settings",
        "--base-settings",
        default=DEFAULT_BASE_SETTINGS,
        help="settings to include as defaults for all combinations",
    )
    parser.add_argument("--no_run_test", "--no-run-test", action="store_true", help="Don't run tests on files")
    parser.add_argument(
        "--runner",
        choices=("module", "command"),
        default="module",
        help="Use the compiled Python camb module or the legacy command-line executable",
    )
    parser.add_argument("--prog", default="./camb", help="executable to run when --runner=command")
    parser.add_argument("--no_validate", "--no-validate", action="store_true", help="Skip ini validation")
    parser.add_argument("--clean", action="store_true", help="delete output dir before run")
    parser.add_argument("--diff_to", "--diff-to", help="output directory to compare to, e.g. test_outputs2")
    parser.add_argument(
        "--diff_tolerance",
        "--diff-tolerance",
        type=float,
        help=(
            "scale factor applied to numerical comparison tolerances; "
            "1 uses the default per-output tolerance table [default: 1]"
        ),
        default=DEFAULT_TOLERANCE_SCALE,
    )
    parser.add_argument(
        "--verbose_diff_output",
        "--verbose-diff-output",
        "--verbose",
        action="store_true",
        help="during diff_to print more error messages",
    )
    parser.add_argument("--num_diff", "--num-diff", action="store_true", help="during diff_to use absolute diffs")
    parser.add_argument("--no_sources", "--no-sources", action="store_true", help="turn off CAMB sources tests")
    parser.add_argument("--no_de", "--no-de", action="store_true", help="Don't run dark energy tests")
    parser.add_argument("--max_tests", "--max-tests", type=int, help="maximum tests to run")
    parser.add_argument(
        "--boosted_reference",
        "--boosted-reference",
        action="store_true",
        help="generate/run ini files with check_accuracy-style boosted reference accuracy settings",
    )
    parser.add_argument(
        "--strict_reference",
        "--strict-reference",
        action="store_true",
        help="raise the boosted reference preset to the strict check_accuracy settings",
    )
    parser.add_argument("--accuracy_boost", "--accuracy-boost", type=float, help="boosted reference AccuracyBoost")
    parser.add_argument("--l_sample_boost", "--l-sample-boost", type=float, help="boosted reference lSampleBoost")
    parser.add_argument("--l_accuracy_boost", "--l-accuracy-boost", type=float, help="boosted reference lAccuracyBoost")
    parser.add_argument(
        "--do_late_rad_truncation",
        "--do-late-rad-truncation",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="boosted reference DoLateRadTruncation",
    )
    parser.add_argument(
        "--min_l_logl_sampling",
        "--min-l-logl-sampling",
        type=int,
        help="boosted reference min_l_logl_sampling",
    )
    parser.add_argument(
        "--override",
        "--set",
        action="append",
        type=parse_override_assignment,
        default=[],
        metavar="KEY=VALUE",
        help="override an ini parameter for every test; can be repeated",
    )
    return parser


args = argparse.Namespace()
prog = ""
out_files_dir = ""

logfile = None


def printlog(text):
    global logfile
    print(text)
    sys.stdout.flush()
    if logfile is None:
        logfile = open(os.path.join(args.ini_dir, "test_results.log"), "a", encoding="utf-8")
    logfile.write(text + "\n")


def close_logfile():
    global logfile
    if logfile is not None:
        logfile.close()
        logfile = None


def repo_inifile(name):
    return os.path.join(REPO_ROOT, "inifiles", name)


def resolve_from_ini_dir(path):
    if os.path.isabs(path):
        return path
    return os.path.join(args.ini_dir, path)


def ensure_camb_on_path():
    if REPO_ROOT not in sys.path:
        sys.path.insert(0, REPO_ROOT)


_check_accuracy_module = None
_recombination_module = None


def get_check_accuracy_module():
    global _check_accuracy_module
    if _check_accuracy_module is None:
        ensure_camb_on_path()
        from camb import check_accuracy

        _check_accuracy_module = check_accuracy
    return _check_accuracy_module


def get_recombination_module():
    global _recombination_module
    if _recombination_module is None:
        ensure_camb_on_path()
        from camb import recombination

        _recombination_module = recombination
    return _recombination_module


def ini_value(value):
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def recfast_approx_settings(recfast_approx_model):
    recombination = get_recombination_module()
    params = recombination.recfast_approx_model_params[recfast_approx_model]
    settings = ["recombination_model = Recfast"]
    for key, value in params.items():
        ini_key = "RECFAST_H_fudge" if key == "RECFAST_fudge" else key
        settings.append(f"{ini_key} = {ini_value(value)}")
    return settings


def recfast_cosmorec_settings():
    recombination = get_recombination_module()
    return recfast_approx_settings(recombination.recfast_cosmorec)


def check_accuracy_reference_defaults():
    check_accuracy = get_check_accuracy_module()
    return {
        "accuracy_boost": check_accuracy.DEFAULT_ACCURACY_SETTINGS["AccuracyBoost"],
        "l_sample_boost": check_accuracy.DEFAULT_ACCURACY_SETTINGS["lSampleBoost"],
        "l_accuracy_boost": check_accuracy.DEFAULT_ACCURACY_SETTINGS["lAccuracyBoost"],
        "do_late_rad_truncation": check_accuracy.DEFAULT_ACCURACY_SETTINGS["DoLateRadTruncation"],
        "min_l_logl_sampling": check_accuracy.DEFAULT_ACCURACY_SETTINGS["min_l_logl_sampling"],
        "reference_min_transfer_k_per_logint": check_accuracy.REFERENCE_MIN_TRANSFER_K_PER_LOGINT,
    }


def check_accuracy_strict_reference_defaults():
    check_accuracy = get_check_accuracy_module()
    return {
        "accuracy_boost": check_accuracy.STRICT_REFERENCE_SETTINGS["AccuracyBoost"],
        "l_sample_boost": check_accuracy.STRICT_REFERENCE_SETTINGS["lSampleBoost"],
        "l_accuracy_boost": check_accuracy.STRICT_REFERENCE_SETTINGS["lAccuracyBoost"],
        "min_l_logl_sampling": check_accuracy.STRICT_REFERENCE_SETTINGS["min_l_logl_sampling"],
    }


def make_ini_file(*args, **kwargs):
    ensure_camb_on_path()
    from camb.inifile import IniFile

    return IniFile(*args, **kwargs)


def apply_ini_overrides(filename):
    if not args.override:
        return filename
    ini = make_ini_file(filename)
    for key, value in args.override:
        if key not in ini.params:
            ini.readOrder.append(key)
        ini.params[key] = value
    ini.saveFile(filename)
    return filename


def set_ini_param(ini, key, value):
    if key not in ini.params:
        ini.readOrder.append(key)
    ini.params[key] = value


def boosted_reference_settings():
    settings = check_accuracy_reference_defaults()
    if args.strict_reference:
        for key, value in check_accuracy_strict_reference_defaults().items():
            if key in settings and isinstance(settings[key], (float, int)):
                settings[key] = max(settings[key], value)
            else:
                settings[key] = value

    for key in [
        "accuracy_boost",
        "l_sample_boost",
        "l_accuracy_boost",
        "do_late_rad_truncation",
        "min_l_logl_sampling",
    ]:
        value = getattr(args, key)
        if value is not None:
            settings[key] = value
    return settings


def apply_boosted_reference_settings(filename):
    if not args.boosted_reference:
        return filename

    ini = make_ini_file(filename)
    settings = boosted_reference_settings()
    for key in ["accuracy_boost", "l_sample_boost", "l_accuracy_boost"]:
        value = float(settings[key])
        if key in ini.params:
            value = max(ini.float(key), value)
        set_ini_param(ini, key, ini_value(value))
    min_l_logl_sampling = int(settings["min_l_logl_sampling"])
    if "min_l_logl_sampling" in ini.params:
        min_l_logl_sampling = max(ini.int("min_l_logl_sampling"), min_l_logl_sampling)
    set_ini_param(ini, "min_l_logl_sampling", ini_value(min_l_logl_sampling))
    if "do_late_rad_truncation" not in ini.params:
        set_ini_param(ini, "do_late_rad_truncation", ini_value(settings["do_late_rad_truncation"]))
    if ini.bool("get_transfer", False):
        for key, value in check_accuracy_strict_reference_defaults().items():
            if key in {"accuracy_boost", "l_sample_boost", "l_accuracy_boost"}:
                set_ini_param(ini, key, ini_value(max(ini.float(key, 1.0), float(value))))
        transfer_k_per_logint = max(
            ini.int("transfer_k_per_logint", 0),
            int(settings["reference_min_transfer_k_per_logint"]),
        )
        set_ini_param(ini, "transfer_k_per_logint", ini_value(transfer_k_per_logint))
    ini.saveFile(filename)
    return filename


def postprocess_test_ini(test_name, filename):
    if test_name not in THETA_CASES:
        return filename

    key, value = THETA_CASES[test_name]
    ini = make_ini_file(filename)
    ini.delete_keys(["hubble"])
    set_ini_param(ini, key, value)
    ini.saveFile(filename)
    return filename


def write_flat_ini_file(filename, parameter_lines):
    filename = os.path.abspath(filename)
    with tempfile.NamedTemporaryFile(
        "w", suffix=".ini", dir=os.path.dirname(filename), delete=False, encoding="utf-8"
    ) as handle:
        handle.write("\n".join(parameter_lines))
        handle.write("\n")
        temp_filename = handle.name

    try:
        make_ini_file(temp_filename).saveFile(filename)
    finally:
        os.remove(temp_filename)

    return filename


def remove_tree(path):
    def onerror(func, target, exc_info):
        for _ in range(10):
            try:
                os.chmod(target, stat.S_IWRITE | stat.S_IREAD)
            except OSError:
                pass
            try:
                func(target)
                return
            except PermissionError:
                time.sleep(0.2)
        raise exc_info[1]

    shutil.rmtree(path, onerror=onerror)


# The tolerance matrix gives the tolerances for comparing two values of the actual results with
# results given in a diff_to. Filename globbing is supported by fnmatch. The first glob matching
# is the winner and its tolerances will be used. To implement a first match order, a regular array
# had to be used instead of a dictionary for the filetolmatrix. The first match is then implemented
# in the routine getToleranceVector().


class ColTol(dict):
    """
    Specify the column tolerances for all columns in a file.
     This class is inherited from the dict class and overrides the missing
     method to return the tolerance of the asterisk key, which is denoting
     the tolerances for all not explicitly specified columns. The
     tolerance for a column can be Ignore()d be using the <b>Ignore()</b> value or
     has to be a tupple, where the first item tells whether the second is to
     be evaluated. The first item can be a bool (value does not matter), to always
     select the second item for evaluation, or a function accepting the
     dictionary of inifile setting. The function then has to return true,
     when the second item of the tupple has to be evaluated.
     A tolerance for the |old-new| < tol can be specified by giving the
     scalar tolerance or a function of two vectors for the second item of
     the tupple. The first vector contains all columns of the old values
     the second all values of the new values. The values are addressed by
     the columns names taken from the newer file. The function has to return
     true, when the new value is ok, false else.
     Additionally ranges of tolerances or functions can be given as a sorted
     list of 2-tupples. The first value is the lower bound of the column "L" for
     which the second value/function is applicable. The lists is traversed as
     long as "L" is smaller then the first value or the list ends. The second
     value is then taken for the comparison. That value also can be an Ignore()
     object, denoting that the value is to always accepted.
    """

    def __missing__(self, item):
        return self["*"]


class Ignore:
    """
    Ignore() files of this class completely.
    """


def diffnsqrt(old, new, tol, c1, c2):
    """
    Implement |C1'x'C2_{new} - C1'x'C2_{old}| / sqrt(C1'x'C1_{old} * C2'x'C2_{old}) < tol.
    :param old: The row of the old values.
    :param new: The row of the new values.
    :param tol: The tolerance to match.
    :param c1: The name of the first component.
    :param c2: The name of the second component.
    :return: True, when |C1'x'C2_{new} - C1'x'C2_{old}| / sqrt(C1'x'C1_{old} * C2'x'C2_{old}) < tol, false else.
    :rtype : bool
    """
    oc1c1 = old[c1 + "x" + c1]
    oc2c2 = old[c2 + "x" + c2]
    # Skip the test when exactly one variable is negative, but not both.
    if (oc1c1 < 0 or oc2c2 < 0) and (oc1c1 >= 0 or oc2c2 >= 0):
        return True
    res = math.fabs(new[c1 + "x" + c2] - old[c1 + "x" + c2]) / math.sqrt(oc1c1 * oc2c2) < tol
    if args.verbose_diff_output and not res:
        printlog(
            "diffnsqrt: |{:g} - {:g}|/sqrt({:g} * {:g}) = {:g} > {:g}".format(
                new[c1 + "x" + c2],
                old[c1 + "x" + c2],
                oc1c1,
                oc2c2,
                math.fabs(new[c1 + "x" + c2] - old[c1 + "x" + c2]) / math.sqrt(oc1c1 * oc2c2),
                tol,
            )
        )
    return res


def normabs(o, n, tol):
    """
    Compute |o - n| / |o| < tol
    :param o: The old value
    :param n: The new value
    :param tol: toleranace
    :return: True when |o - n| / |o| < tol, false else
    """
    res = (math.fabs(o - n) / math.fabs(o) if o != 0.0 else math.fabs(o - n)) < tol
    if args.verbose_diff_output and not res:
        printlog(
            f"normabs: |{o:g} - {n:g}| / |{o:g}| = {math.fabs(o - n) / math.fabs(o) if o != 0.0 else math.fabs(o - n):g} > {tol:g}"
        )
    return res


# The filetolmatrix as described above.
filetolmatrix = [
    ["*scalCls.dat", Ignore()],  # Ignore() all scalCls.dat files.
    [
        "*tensCls.dat",
        ColTol(
            {
                "TE": (True, lambda o, n: diffnsqrt(o, n, scale_tolerance(1e-2), "T", "E")),
                "*": (True, [(0, 1e-2), (600, Ignore())]),
            }
        ),
    ],
    ["*sharp_cl_*.dat", ColTol({"CL": (True, 1e-3), "P": (True, 1e-3), "P_vv": (True, 1e-3), "*": Ignore()})],
    ["*", ColTol({"*": (True, 1e-4)})],
]


def runScript(fname):
    now = time.time()
    try:
        if args.runner == "module":
            ensure_camb_on_path()
            import camb

            camb.run_ini(fname, no_validate=args.no_validate)
            res = ""
        else:
            res = subprocess.check_output([prog, fname], stderr=subprocess.STDOUT, text=True)
        code = 0
    except subprocess.CalledProcessError as error:
        res = error.output
        code = error.returncode
    except Exception as error:
        res = str(error)
        code = 1
    return time.time() - now, res, code


def getInis(ini_dir):
    ini_files = []
    for fname in sorted(os.listdir(ini_dir)):
        if fnmatch.fnmatch(fname, "params_*.ini"):
            ini_files.append(os.path.join(args.ini_dir, fname))
    return ini_files


def getTestParams():
    params = [["base"]]

    for lmax in [2500, 3000, 4500, 6000]:
        params.append([f"lmax{lmax}", f"l_max_scalar = {lmax}", "k_eta_max_scalar  = %s" % (lmax * 2.5)])

    for lmax in [2500, 3000, 4500]:
        params.append(
            [
                f"nonlin_lmax{lmax}",
                "do_nonlinear =2",
                "get_transfer= T",
                f"l_max_scalar = {lmax}",
                "k_eta_max_scalar  = %s" % (lmax * 2.5),
            ]
        )

    for lmax in [400, 600, 1000]:
        params.append(
            [
                f"tensor_lmax{lmax}",
                "get_tensor_cls = T",
                f"l_max_tensor = {lmax}",
                "k_eta_max_tensor  = %s" % (lmax * 2),
            ]
        )

    params.append(["tensoronly", "get_scalar_cls=F", "get_tensor_cls = T"])
    params.append(
        ["tensor_transfer", "get_scalar_cls=F", "get_tensor_cls = T", "get_transfer= T", "transfer_high_precision = T"]
    )
    params.append(["transfer_only", "get_scalar_cls=F", "get_transfer= T", "transfer_high_precision = F"])
    params.append(["transfer_highprec", "get_scalar_cls=F", "get_transfer= T", "transfer_high_precision = T"])

    params.append(["all", "get_scalar_cls=T", "get_tensor_cls = T", "get_transfer= T"])
    params.append(["all_nonlin1", "get_scalar_cls=T", "get_tensor_cls = T", "get_transfer= T", "do_nonlinear=1"])
    params.append(["all_nonlin2", "get_scalar_cls=T", "get_tensor_cls = T", "get_transfer= T", "do_nonlinear=2"])
    params.append(
        [
            "all_nonlinhigh",
            "get_scalar_cls=T",
            "get_tensor_cls = T",
            "get_transfer= T",
            "do_nonlinear=2",
            "transfer_high_precision = T",
        ]
    )
    params.append(
        [
            "transfer_delta10",
            "get_scalar_cls=F",
            "get_transfer= T",
            "transfer_high_precision = T",
            "transfer_k_per_logint =10",
            "transfer_matterpower(1) =",
        ]
    )
    params.append(
        ["transfer_redshifts", "get_scalar_cls=F", "get_transfer= T", "transfer_num_redshifts=2"]
        + [
            "transfer_redshift(1)=1",
            "transfer_redshift(2)=0.7",
            "transfer_filename(2)=transfer_out2.dat",
            "transfer_matterpower(2)=matterpower2.dat",
        ]
    )
    params.append(
        ["transfer_redshifts2", "get_scalar_cls=F", "get_transfer= T", "transfer_num_redshifts=2"]
        + [
            "transfer_redshift(1)=0.7",
            "transfer_redshift(2)=0",
            "transfer_filename(2)=transfer_out2.dat",
            "transfer_matterpower(2)=matterpower2.dat",
        ]
    )

    params.append(["transfer_nonu", "get_scalar_cls=F", "get_transfer= T", "transfer_power_var = 8"])
    params.append(
        [
            "raw_transfer_two_redshifts",
            "get_scalar_cls = F",
            "get_transfer = T",
            "transfer_high_precision = T",
            "transfer_kmax = 5",
            "transfer_k_per_logint = 50",
            "transfer_num_redshifts = 2",
            "transfer_redshift(1) = 1.0",
            "transfer_redshift(2) = 0.0",
            "transfer_filename(2) = transfer_out_z0.dat",
            "transfer_matterpower(2) = matterpower_z0.dat",
        ]
    )

    # AM - Added HMcode and halomodel tests (halofit_version=5,6)
    params.append(["HMcode", "transfer_kmax=100", "halofit_version=5", "do_nonlinear=1", "get_transfer= T"])
    params.append(["halomodel", "transfer_kmax=100", "halofit_version=6", "do_nonlinear=1", "get_transfer= T"])
    # AM - End of edits

    params.append(["zre", "re_use_optical_depth = F", "re_redshift  = 8.5"])
    params.append(["noderived", "derived_parameters = F"])
    params.append(["no_rad_trunc", "do_late_rad_truncation   = F"])

    recfast_cosmorec = recfast_cosmorec_settings()
    params.append(
        [
            "lmax4000_lpa_auto",
            *LPA_AUTO_LMAX_SETTINGS[4000],
            "ombh2 = 0.0222",
            "omch2 = 0.115",
            "scalar_spectral_index(1) = 0.965",
        ]
    )
    params.append(
        [
            "lmax6000_lpa_auto",
            *LPA_AUTO_LMAX_SETTINGS[6000],
            "ombh2 = 0.0229",
            "omch2 = 0.109",
            "re_optical_depth = 0.06",
        ]
    )
    params.append(
        ["thetastar_lmax4000_lpa_auto"]
        + LPA_AUTO_LMAX_SETTINGS[4000]
        + [
            "ombh2 = 0.0219",
            "omch2 = 0.118",
            "re_optical_depth = 0.06",
        ]
        + recfast_cosmorec
    )
    params.append(
        ["cosmomc_theta_lmax4000_lpa_auto"]
        + LPA_AUTO_LMAX_SETTINGS[4000]
        + [
            "ombh2 = 0.0231",
            "omch2 = 0.108",
            "helium_fraction = 0.245",
        ]
        + recfast_cosmorec
    )

    for acc in [1.1, 1.5, 2.2]:
        params.append([f"accuracy_boost{acc}", f"accuracy_boost = {acc}"])

    for acc in [1, 1.5, 2]:
        params.append([f"l_accuracy_boost{acc}", f"l_accuracy_boost = {acc}"])

    params.append(["nu_massless", "omnuh2 =0"])

    for mnu in [0, 0.01, 0.03, 0.1]:
        omnu = mnu / 100.0
        params.append([f"nu_mass{mnu}", f"omnuh2 ={omnu}", "massive_neutrinos  = 3"])
    params.append(
        [
            "nu_masssplit",
            "omnuh2 =0.03",
            "massive_neutrinos = 1 1",
            "nu_mass_fractions=0.2 0.8",
            "nu_mass_degeneracies = 1 1",
            "nu_mass_eigenstates = 2",
            "massless_neutrinos = 1.046",
        ]
    )
    for omnuh2 in [0.0032, 0.00325, 0.0134, 0.02, 0.023, 0.024, 0.025, 0.03]:
        case_settings = [
            f"omnuh2 = {omnuh2}",
            "massless_neutrinos = 0.046",
            "massive_neutrinos = 1 1 1",
            "nu_mass_fractions = 0.3333333333333333 0.3333333333333333 0.3333333333333333",
            "nu_mass_degeneracies = 1.0 1.0 1.0",
            "nu_mass_eigenstates = 3",
            "share_delta_neff = F",
            "l_max_scalar = 4200",
            "k_eta_max_scalar = 108000",
            "lens_output_margin = 200",
            "do_nonlinear = 2",
        ]
        if omnuh2 == 0.0134:
            case_settings.append("accurate_BB = T")
        params.append(
            [
                f"three_deg_omnuh2_{omnuh2:.5f}",
                *case_settings,
            ]
        )

    pars = {
        "ombh2": [0.0219, 0.0226, 0.0253],
        "omch2": [0.1, 0.08, 0.15],
        "omk": [0, -0.03, 0.04, 0.001, -0.001],
        "hubble": [62, 67, 71, 78],
        "w": [-1.2, -1, -0.98, -0.75],
        "helium_fraction": [0.21, 0.23, 0.27],
        "scalar_spectral_index(1)": [0.94, 0.98],
        "scalar_nrun(1)": [-0.015, 0, 0.03],
        "re_optical_depth": [0.03, 0.05, 0.08, 0.11],
    }

    for par, vals in pars.items():
        for val in vals:
            params.append(
                [
                    f"{par}_{val:.3f}",
                    "get_transfer= T",
                    "do_nonlinear=1",
                    "transfer_high_precision = T",
                    f"{par} = {val}",
                ]
            )

    for omk in [-0.03, -0.02, -0.01, -0.001, -1e-5, 0.0, 1e-5, 0.001, 0.01, 0.02, 0.03]:
        label = f"nearflat_omk_{omk:+.5g}".replace("+", "p").replace("-", "m").replace(".", "p")
        params.append(
            [
                label,
                f"omk = {omk}",
                "l_max_scalar = 4200",
                "k_eta_max_scalar = 10500",
                "lens_output_margin = 200",
            ]
        )

    for omk in [-1e-6, -5.1e-7, 5.1e-7, 1e-6]:
        label = f"pp_continuity_omk_{omk:+.2g}".replace("+", "p").replace("-", "m").replace(".", "p")
        params.append(
            [
                label,
                f"omk = {omk}",
                "l_max_scalar = 4200",
                "k_eta_max_scalar = 10500",
                "lens_output_margin = 200",
            ]
        )

    if not args.no_de and not os.environ.get("CAMB_TESTS_NO_DE"):
        for wa in [-0.3, -0.01, 0.5]:
            for w in [-1.2, -0.998, -0.7]:
                params.append(
                    [
                        f"ppf_w{w}_wa{wa}",
                        f"w = {w}",
                        f"wa ={wa}",
                        "do_nonlinear = 2",
                        "get_transfer= T",
                        "dark_energy_model=PPF",
                    ]
                )

        ppf_theta_lmax4000 = [
            "l_max_scalar = 4200",
            "k_eta_max_scalar = 18000",
            "lens_output_margin = 200",
            "do_nonlinear = 3",
            "get_transfer = T",
            "transfer_high_precision = T",
            "transfer_kmax = 5",
            "transfer_num_redshifts = 3",
            "transfer_redshift(1) = 1.0",
            "transfer_redshift(2) = 0.5",
            "transfer_redshift(3) = 0.0",
            "dark_energy_model = PPF",
            "ombh2 = 0.02237",
            "omch2 = 0.12",
            "omnuh2 = 0.000644866570625114",
            "re_optical_depth = 0.0544",
        ]
        params.extend(
            [
                [
                    "cosmomc_theta_ppf_w-0p82_wa0p35_lmax4000",
                    "w = -0.82",
                    "wa = 0.35",
                    *ppf_theta_lmax4000,
                ],
                [
                    "thetastar_ppf_w-1p2_wa0p5_lmax4000",
                    "w = -1.2",
                    "wa = 0.5",
                    *ppf_theta_lmax4000,
                ],
                [
                    "thetastar_ppf_w-0p7_wa-0p3_lmax4000",
                    "w = -0.7",
                    "wa = -0.3",
                    *ppf_theta_lmax4000,
                ],
            ]
        )

        params.append(
            [
                "ppf_w-1.000_wa0.000",
                "w = -1.0",
                "wa = 0.0",
                "do_nonlinear = 1",
                "get_transfer= T",
                "transfer_high_precision = T",
                "dark_energy_model=PPF",
            ]
        )

    if not args.no_sources and not os.environ.get("CAMB_TESTS_NO_SOURCES"):
        # ##CAMB sources options and new outputs
        params.append(
            ["delta_xe", "evolve_delta_xe =T", "get_transfer= T", "do_nonlinear=2", "transfer_high_precision = T"]
        )

        def make_win(i, z, kind, bias, sigma, s):
            return [
                f"redshift({i}) = {z}",
                f"redshift_kind({i}) = {kind}",
                f"redshift_bias({i}) = {bias}",
                f"redshift_sigma({i}) = {sigma}",
                f"redshift_dlog10Ndm({i}) = {s}",
            ]

        source_lmax = ["l_max_scalar = 2500", "k_eta_max_scalar = 20000"]
        counts_def = [f"DEFAULT({repo_inifile('params_counts.ini')})"]
        source_counts = (
            source_lmax
            + ["num_redshiftwindows = 2"]
            + make_win(1, 0.3, "counts", 1.5, 0.06, 0.42)
            + make_win(2, 1, "counts", 2, 0.3, 0)
        )
        bool_options = ["counts_evolve", "DoRedshiftLensing", "counts_redshift", "evolve_delta_xe"]
        for b1 in ["T", "F"]:
            for b2 in ["T", "F"]:
                for b3 in ["T", "F"]:
                    for b4 in ["T", "F"]:
                        bs = [b1, b2, b3, b4]
                        pars = copy.copy(source_counts)
                        for opt, b in zip(bool_options, bs):
                            pars += [opt + " = " + b]
                        params.append(["counts_opts_" + "_".join(bs)] + counts_def + pars)
        params.append(
            ["counts_1bin"]
            + counts_def
            + source_lmax
            + ["num_redshiftwindows = 1"]
            + make_win(1, 0.15, "counts", 1.2, 0.04, -0.2)
        )
        params.append(["counts_no_cmb_lmax2500", "want_CMB = F"] + counts_def + source_counts)

        params.append(
            ["counts_overlap"]
            + counts_def
            + source_lmax
            + ["num_redshiftwindows = 2"]
            + make_win(1, 0.17, "counts", 1.2, 0.04, -0.2)
            + make_win(2, 0.2, "counts", 1.2, 0.04, -0.2)
        )
        params.append(["lensing_base", f"DEFAULT({repo_inifile('params_lensing.ini')})"] + source_lmax)
        params.append(["21cm_base", f"DEFAULT({repo_inifile('params_21cm.ini')})", "l_max_scalar = 6000"])
        params.append(
            ["21cm_base2", f"DEFAULT({repo_inifile('params_21cm.ini')})", "get_transfer = T", "l_max_scalar = 6000"]
        )
        params.append(
            ["counts_lens", f"DEFAULT({repo_inifile('params_counts.ini')})"]
            + source_lmax
            + ["num_redshiftwindows = 2"]
            + make_win(1, 0.17, "counts", 1.2, 0.04, -0.2)
            + make_win(2, 0.5, "lensing", 0, 0.07, 0.2)
        )

    max_tests = args.max_tests or os.environ.get("CAMB_TESTS_MAX")
    if max_tests:
        params = params[: int(max_tests)]
    return params


def list_files(file_dir):
    return [f for f in os.listdir(file_dir) if ".ini" not in f]


def output_file_num(file_dir):
    return len(list_files(file_dir))


def remove_stale_generated_inis(expected_ini_files):
    expected = {os.path.basename(filename) for filename in expected_ini_files}
    for fname in os.listdir(args.ini_dir):
        if fnmatch.fnmatch(fname, "params_*.ini") and fname not in expected:
            os.remove(os.path.join(args.ini_dir, fname))


def makeIniFiles():
    printlog("Making test ini files...")
    params = getTestParams()
    ini_files = []
    base_settings = os.path.abspath(args.base_settings)
    for pars in params:
        name = "params_" + pars[0]
        fname = os.path.join(args.ini_dir, name + ".ini")
        ini_files.append(fname)
    remove_stale_generated_inis(ini_files)
    for pars, fname in zip(params, ini_files):
        name = "params_" + pars[0]
        write_flat_ini_file(
            fname,
            [
                "output_root=" + os.path.join(out_files_dir, name),
                *pars[1:],
                *SUPPRESSED_CMB_OUTPUT_SETTINGS,
                f"DEFAULT({base_settings})",
            ],
        )
        postprocess_test_ini(pars[0], fname)
        apply_ini_overrides(fname)
        apply_boosted_reference_settings(fname)
    printlog("Made test ini files.")
    return ini_files


def getPreparedIniFiles():
    if args.make_ini:
        return makeIniFiles(), None

    ini_files = getInis(args.ini_dir)
    if not args.override and not args.boosted_reference:
        return ini_files, None

    override_dir = tempfile.mkdtemp(prefix="test_ini_overrides_", dir=args.ini_dir)
    prepared = []
    for ini in ini_files:
        copied_ini = os.path.join(override_dir, os.path.basename(ini))
        shutil.copy(ini, copied_ini)
        apply_ini_overrides(copied_ini)
        prepared.append(apply_boosted_reference_settings(copied_ini))
    return prepared, override_dir


def get_tolerance_vector(filename, cols):
    """
    Get the tolerances for the given filename.
    :param filename: The name of the file to retrieve the tolerances for.
    :param cols: Gives the column names.
    :returns: False, when the file is to be Ignored completely;
    the vector of tolerances when a pattern in the filetolmatrix matched;
    an empty ColTol when no match was found.
    """
    for key, val in filetolmatrix:
        if fnmatch.fnmatch(filename, key):
            if isinstance(val, Ignore):
                return False
            else:
                return [val.get(t, (False, False)) for t in cols]
    return ColTol()


def parse_float(value):
    try:
        return float(value)
    except ValueError:
        sp = customsplit(value)
        return float(sp[0] + "E" + sp[1])


def numeric_row(row):
    return [parse_float(value) for value in row]


def align_rows_by_key(orig_rows, new_rows, key_func):
    try:
        orig_by_key = {key_func(row): row for row in orig_rows}
        new_by_key = {key_func(row): row for row in new_rows}
    except (IndexError, ValueError):
        return None
    common_keys = sorted(set(orig_by_key).intersection(new_by_key))
    if not common_keys:
        return None
    return common_keys, orig_by_key, new_by_key


def align_rows_by_l(orig_rows, new_rows):
    aligned = align_rows_by_key(orig_rows, new_rows, lambda row: int(numeric_row(row)[0]))
    if aligned is None:
        return None
    common_l, orig_by_l, new_by_l = aligned
    return [orig_by_l[ell] for ell in common_l], [new_by_l[ell] for ell in common_l]


def interpolate_row(x, orig_x, orig_coord, orig_rows):
    coord = math.log(x)
    index = bisect.bisect_left(orig_coord, coord)
    if index < len(orig_coord) and orig_x[index] == x:
        return orig_rows[index]
    if index == 0 or index == len(orig_coord):
        return None
    low_coord = orig_coord[index - 1]
    high_coord = orig_coord[index]
    weight = (coord - low_coord) / (high_coord - low_coord)
    low_row = orig_rows[index - 1]
    high_row = orig_rows[index]
    return [low + weight * (high - low) for low, high in zip(low_row, high_row)]


def align_rows_by_interpolated_k(orig_rows, new_rows):
    orig_numeric = [numeric_row(row) for row in orig_rows]
    new_numeric = [numeric_row(row) for row in new_rows]
    if not orig_numeric or not new_numeric:
        return None
    if any(row[0] <= 0 for row in orig_numeric + new_numeric):
        return None
    orig_x = [row[0] for row in orig_numeric]
    if any(x2 <= x1 for x1, x2 in zip(orig_x, orig_x[1:])):
        return None
    orig_coord = [math.log(x) for x in orig_x]
    aligned_orig = []
    aligned_new = []
    for row in new_numeric:
        x = row[0]
        if x < orig_x[0] or x > orig_x[-1]:
            continue
        interpolated = interpolate_row(x, orig_x, orig_coord, orig_numeric)
        if interpolated is None:
            continue
        aligned_orig.append(interpolated)
        aligned_new.append(row)
    if not aligned_orig:
        return None
    return aligned_orig, aligned_new


def align_mismatched_rows(filename, orig_rows, new_rows, cols):
    if not cols:
        return None
    first_col = cols[0]
    if first_col == "L":
        return align_rows_by_l(orig_rows, new_rows)
    if first_col == "k/h":
        return align_rows_by_interpolated_k(orig_rows, new_rows)
    return None


def params_ini_for_output(filename):
    suffixes = [
        "_lenspotentialCls.dat",
        "_lensedtotCls.dat",
        "_scalarCovCls.dat",
        "_transfer_out_z0.dat",
        "_matterpower_z0.dat",
        "_transfer_out2.dat",
        "_matterpower2.dat",
        "_transfer_out.dat",
        "_matterpower.dat",
        "_lensedCls.dat",
        "_sharp_cl_1.dat",
        "_sharp_cl_2.dat",
        "_sharp_cl_3.dat",
        "_scalCls.dat",
        "_tensCls.dat",
        "_totCls.dat",
    ]
    for suffix in suffixes:
        if filename.endswith(suffix):
            return filename[: -len(suffix)] + "_params.ini"
    inifilenameparts = filename.rsplit("_", 2)
    inifilename = (
        "_".join(inifilenameparts[0:2])
        if inifilenameparts[1] not in {"transfer", "matterpower"}
        else inifilenameparts[0]
    )
    return inifilename + "_params.ini"


_diff_ini_cache = {}
_diff_reference_ini_cache = {}


def get_diff_ini(filename):
    inifilename = os.path.join(args.ini_dir, args.out_files_dir, params_ini_for_output(filename))
    if inifilename not in _diff_ini_cache:
        inifile = make_ini_file()
        if os.path.exists(inifilename):
            inifile.readFile(inifilename)
        _diff_ini_cache[inifilename] = inifile
    return _diff_ini_cache[inifilename]


def get_diff_reference_ini(filename):
    inifilename = os.path.join(resolve_from_ini_dir(args.diff_to), params_ini_for_output(filename))
    if inifilename not in _diff_reference_ini_cache:
        inifile = make_ini_file()
        if os.path.exists(inifilename):
            inifile.readFile(inifilename)
        _diff_reference_ini_cache[inifilename] = inifile
    return _diff_reference_ini_cache[inifilename]


def is_boosted_reference_diff(filename):
    inifile = get_diff_ini(filename)
    reference = get_diff_reference_ini(filename)
    for key in ["accuracy_boost", "l_sample_boost", "l_accuracy_boost"]:
        if reference.float(key, 1.0) > inifile.float(key, 1.0):
            return True
    return reference.int("min_l_logl_sampling", 0) > inifile.int("min_l_logl_sampling", 0)


def read_numeric_output_table(filename):
    with open(filename, encoding="utf-8") as f:
        matrix = [row for row in (line.split() for line in f) if row]
    if len(matrix) < 2:
        raise ValueError(f"{filename} has fewer than two non-empty rows")
    if len(matrix[0]) == len(matrix[1]) + 1:
        matrix[0] = matrix[0][1:]
        base = 1
    else:
        base = 0
    cols = [s[0] + "x" + s[1] if len(s) == 2 and s != "nu" else s for s in matrix[0]]
    return cols, [numeric_row(row) for row in matrix[base:]]


def align_numeric_rows_by_l(orig_rows, new_rows):
    aligned = align_rows_by_key(orig_rows, new_rows, lambda row: int(row[0]))
    if aligned is None:
        return None
    common_l, orig_by_l, new_by_l = aligned
    return [(ell, orig_by_l[ell], new_by_l[ell]) for ell in common_l]


def finite_fractional_delta(value, reference):
    denominator = math.fabs(reference)
    if denominator == 0:
        return 0.0 if value == reference else math.inf
    return (value - reference) / denominator


def finite_normalized_cross_delta(value, reference, reference_diag1, reference_diag2):
    if (reference_diag1 < 0 or reference_diag2 < 0) and (reference_diag1 >= 0 or reference_diag2 >= 0):
        return 0.0
    denominator = math.sqrt(max(reference_diag1 * reference_diag2, 0.0))
    if denominator == 0:
        return 0.0 if value == reference else math.inf
    return (value - reference) / denominator


def tolerance_scale():
    return getattr(args, "diff_tolerance", DEFAULT_TOLERANCE_SCALE)


def scale_tolerance(tolerance):
    return tolerance * tolerance_scale()


def check_accuracy_rows_failed(filename, rows):
    scale = tolerance_scale()
    failed = False
    for row in rows:
        tolerance = row.tolerance * scale
        if row.max_abs <= tolerance:
            continue
        failed = True
        if args.verbose_diff_output:
            location = ""
            if row.location:
                location_prefix = "at" if "=" in str(row.location) else "at L="
                location = f" {location_prefix}{row.location}"
            scale_note = f" (base {row.tolerance:g}, scale {scale:g})" if scale != 1 else ""
            printlog(
                f"check_accuracy-style mismatch in {filename}: {row.quantity} {row.range_label} "
                f"max {row.max_abs:g} > {tolerance:g}{scale_note}{location}"
            )
    return failed


def compare_error_ranges(filename, quantity, ell_errors, ranges, min_ell=2):
    check_accuracy = get_check_accuracy_module()
    ell = np.array([ell for ell, _ in ell_errors], dtype=int)
    errors = np.array([error for _, error in ell_errors], dtype=float)
    return check_accuracy_rows_failed(
        filename,
        check_accuracy.compare_l_ranges(quantity, errors, ell, ranges, min_ell=min_ell),
    )


def column_error_series(pairs, cols, column):
    col = cols.index(column)
    return [(ell, finite_fractional_delta(new_row[col], orig_row[col])) for ell, orig_row, new_row in pairs]


def cap_pairs(pairs, max_l):
    if max_l is None:
        return pairs
    capped = [pair for pair in pairs if pair[0] <= max_l]
    return capped or pairs


def check_accuracy_cmb_pairs(pairs, inifile):
    # These caps preserve the historical text-output comparison range; the tolerances
    # and pass/fail statistics themselves are delegated to camb.check_accuracy below.
    if inifile.bool("Do21cm") and pairs[-1][0] > 6000:
        return cap_pairs(pairs, 6000)
    if inifile.int("l_max_scalar", 0) <= 3000 and pairs[-1][0] > 2500:
        return cap_pairs(pairs, 2500)
    return pairs


def check_accuracy_lensing_pairs(pairs, inifile):
    if inifile.bool("Do21cm") and pairs[-1][0] > 6000:
        return cap_pairs(pairs, 6000)
    return pairs


def check_accuracy_params_from_ini(inifile):
    want_scalars = inifile.bool("get_scalar_cls", True)
    want_tensors = inifile.bool("get_tensor_cls")
    want_vectors = inifile.bool("get_vector_cls")
    return SimpleNamespace(
        WantCls=want_scalars or want_tensors or want_vectors,
        WantScalars=want_scalars,
        WantTensors=want_tensors,
        DoLensing=inifile.bool("do_lensing", True),
        Want_CMB=inifile.bool("want_CMB", True),
        Accuracy=SimpleNamespace(AccurateBB=inifile.bool("accurate_BB")),
    )


def check_accuracy_run_output(label, params, *, lensed_cls=None, lens_potential_cls=None, matter_power=None):
    check_accuracy = get_check_accuracy_module()
    return check_accuracy.RunOutput(
        label,
        params,
        None,
        0.0,
        0.0,
        {},
        lensed_cls,
        lens_potential_cls,
        matter_power,
    )


def cls_array_from_pairs(pairs, cols, *, use_new):
    check_accuracy = get_check_accuracy_module()
    cls = np.zeros((max(ell for ell, _, _ in pairs) + 1, len(check_accuracy.CL_COLUMNS)))
    columns = {"TxT": "TT", "ExE": "EE", "BxB": "BB", "TxE": "TE"}
    for text_column, check_accuracy_column in columns.items():
        col = cols.index(text_column)
        target_col = check_accuracy.CL_COLUMNS[check_accuracy_column]
        for ell, orig_row, new_row in pairs:
            row = new_row if use_new else orig_row
            cls[ell, target_col] = row[col]
    return cls


def lensing_array_from_pairs(pairs, cols, *, use_new):
    pp_col = cols.index("PxP")
    lensing_cls = np.zeros((max(ell for ell, _, _ in pairs) + 1, 1))
    for ell, orig_row, new_row in pairs:
        row = new_row if use_new else orig_row
        lensing_cls[ell, 0] = row[pp_col]
    return lensing_cls


def matter_power_data_from_rows(rows, cols, inifile, filename):
    if "k/h" not in cols or "P" not in cols:
        raise ValueError(f"{filename} does not have k/h and P columns")
    k_col = cols.index("k/h")
    p_col = cols.index("P")
    k = np.array([row[k_col] for row in rows], dtype=float)
    if len(k) == 0:
        raise ValueError(f"{filename} has no matter power rows")
    if not np.all(k > 0):
        raise ValueError(f"{filename} has non-positive k/h values")
    pk = np.array([[row[p_col] for row in rows]], dtype=float)
    z = np.array([matter_power_redshift_for_output(filename, inifile)], dtype=float)
    return get_check_accuracy_module().MatterPowerData(
        k=k,
        z=z,
        pk=pk,
        requested_kmax=matter_power_requested_kmax(inifile, float(k[-1])),
        npoints=len(k),
        compare_at_input_nodes=inifile.int("transfer_k_per_logint", 0) != 0,
    )


def output_basename(inifile):
    return os.path.basename(inifile.string("output_root", ""))


def indexed_ini_string(inifile, name, index, default=""):
    return inifile.string(f"{name}({index})", default)


def matter_power_redshift_for_output(filename, inifile):
    root = output_basename(inifile)
    for index in range(1, inifile.int("transfer_num_redshifts", 1) + 1):
        output_name = indexed_ini_string(inifile, "transfer_matterpower", index)
        if output_name and filename == f"{root}_{os.path.basename(output_name)}":
            return inifile.float(f"transfer_redshift({index})", 0.0)
    return inifile.float("transfer_redshift(1)", 0.0)


def matter_power_requested_kmax(inifile, output_kmax):
    if not inifile.hasKey("transfer_kmax") or not inifile.hasKey("hubble"):
        return output_kmax
    return min(output_kmax, inifile.float("transfer_kmax") / (inifile.float("hubble") / 100.0))


def check_accuracy_matter_power_params_from_ini(inifile):
    nonlinear = {
        0: "NonLinear_none",
        1: "NonLinear_pk",
        2: "NonLinear_lens",
        3: "NonLinear_both",
    }.get(inifile.int("do_nonlinear", 0), "NonLinear_none")
    return SimpleNamespace(
        NonLinear=nonlinear,
        Transfer=SimpleNamespace(high_precision=inifile.bool("transfer_high_precision", False)),
    )


def compare_check_accuracy_matter_power(filename, orig_rows, new_rows, cols, inifile):
    check_accuracy = get_check_accuracy_module()
    try:
        standard_matter_power = matter_power_data_from_rows(new_rows, cols, inifile, filename)
        reference_matter_power = matter_power_data_from_rows(orig_rows, cols, inifile, filename)
        params = check_accuracy_matter_power_params_from_ini(inifile)
        tolerance = check_accuracy.default_mpk_tolerance(params)
        standard = check_accuracy_run_output("standard", params, matter_power=standard_matter_power)
        reference = check_accuracy_run_output("reference", params, matter_power=reference_matter_power)
        rows = check_accuracy.compare_matter_power(standard, reference, tolerance)
    except ValueError as error:
        if args.verbose_diff_output:
            printlog(f"could not compare matter power in {filename}: {error}")
        return True
    return check_accuracy_rows_failed(filename, rows)


def is_transfer_output_file(filename):
    return any(fnmatch.fnmatch(filename, pattern) for pattern in TRANSFER_FILE_PATTERNS)


def tolerance_ranges(tolerance):
    if isinstance(tolerance, float):
        return [(None, None, tolerance)]
    return tolerance


def signed_log_interpolate_to_grid(source_k, source_values, target_k):
    log_source_k = np.log(source_k)
    log_target_k = np.log(target_k)
    interpolated = np.empty_like(target_k)
    for index, k_value in enumerate(target_k):
        exact = np.searchsorted(source_k, k_value)
        if exact < len(source_k) and source_k[exact] == k_value:
            interpolated[index] = source_values[exact]
            continue
        upper = np.searchsorted(source_k, k_value)
        if upper == 0 or upper == len(source_k):
            interpolated[index] = np.nan
            continue
        values = source_values[upper - 1 : upper + 1]
        if values[0] != 0 and values[1] != 0 and np.sign(values[0]) == np.sign(values[1]):
            interpolated[index] = np.sign(values[0]) * np.exp(
                np.interp(log_target_k[index], log_source_k[upper - 1 : upper + 1], np.log(np.abs(values)))
            )
        else:
            interpolated[index] = np.interp(log_target_k[index], log_source_k[upper - 1 : upper + 1], values)
    return interpolated


def transfer_reference_on_standard_grid(reference_rows, standard_rows, reference_col, standard_col):
    reference_k = np.array([row[0] for row in reference_rows])
    standard_k = np.array([row[0] for row in standard_rows])
    if not len(reference_k) or not len(standard_k):
        raise ValueError("empty transfer table")
    if np.any(reference_k <= 0) or np.any(standard_k <= 0):
        raise ValueError("transfer k/h values must be positive")
    if np.any(reference_k[1:] <= reference_k[:-1]):
        raise ValueError("reference transfer k/h values are not strictly increasing")
    in_range = (standard_k >= reference_k[0]) & (standard_k <= reference_k[-1])
    if not np.any(in_range):
        raise ValueError("transfer k/h grids do not overlap")
    reference_values = np.array([row[reference_col] for row in reference_rows])
    standard_values = np.array([row[standard_col] for row in standard_rows])[in_range]
    common_k = standard_k[in_range]
    interpolated_reference = signed_log_interpolate_to_grid(reference_k, reference_values, common_k)
    return common_k, standard_values, interpolated_reference


def transfer_stat_rows(column, common_k, standard_values, reference_values, mpk_tolerance):
    check_accuracy = get_check_accuracy_module()
    errors = check_accuracy.fractional_delta(standard_values, reference_values)
    locations = np.array([f"k/h={k:.6g}" for k in common_k])
    rows = []
    for start, stop, tolerance in tolerance_ranges(mpk_tolerance):
        mask = np.ones_like(common_k, dtype=bool)
        if start is not None:
            mask &= common_k >= start
        if stop is not None:
            mask &= common_k < stop
        if not np.any(mask):
            continue
        rows.append(
            check_accuracy.finite_stats(
                f"transfer {column}",
                check_accuracy.matter_power_range_label(start, stop),
                errors[mask],
                tolerance / 2,
                locations=locations[mask],
            )
        )
    return rows


def compare_check_accuracy_transfer(filename, reference_rows, standard_rows, cols, inifile):
    if "k/h" not in cols:
        return False
    check_accuracy = get_check_accuracy_module()
    mpk_tolerance = check_accuracy.default_mpk_tolerance(check_accuracy_matter_power_params_from_ini(inifile))
    columns = list(TRANSFER_CHECKED_COLUMNS)
    power_column = TRANSFER_POWER_VAR_COLUMNS.get(inifile.int("transfer_power_var", 7), "total")
    if power_column not in columns:
        columns.append(power_column)
    failed = False
    for column in columns:
        if column not in cols:
            continue
        col_index = cols.index(column)
        try:
            common_k, standard_values, reference_values = transfer_reference_on_standard_grid(
                reference_rows,
                standard_rows,
                col_index,
                col_index,
            )
        except ValueError as error:
            if args.verbose_diff_output:
                printlog(f"could not compare transfer column {column} in {filename}: {error}")
            return True
        rows = transfer_stat_rows(column, common_k, standard_values, reference_values, mpk_tolerance)
        failed |= check_accuracy_rows_failed(filename, rows)
    return failed


def should_compare_lensed_cmb_file(filename, inifile):
    if not fnmatch.fnmatch(filename, "*lensedCls.dat"):
        return True
    return not inifile.bool("get_tensor_cls") or not inifile.string("lensed_total_output_file", "")


def compare_check_accuracy_cmb(filename, pairs, cols, inifile):
    if not inifile.bool("want_CMB", True):
        return False
    if any(column not in cols for column in ["TxT", "ExE", "BxB", "TxE"]):
        return False
    check_accuracy = get_check_accuracy_module()
    params = check_accuracy_params_from_ini(inifile)
    pairs = check_accuracy_cmb_pairs(pairs, inifile)
    standard = check_accuracy_run_output(
        "standard",
        params,
        lensed_cls=cls_array_from_pairs(pairs, cols, use_new=True),
    )
    reference = check_accuracy_run_output(
        "reference",
        params,
        lensed_cls=cls_array_from_pairs(pairs, cols, use_new=False),
    )
    return check_accuracy_rows_failed(filename, check_accuracy.compare_cls(standard, reference))


def compare_check_accuracy_lensing_pp(filename, pairs, cols, inifile):
    if "PxP" not in cols:
        return False
    check_accuracy = get_check_accuracy_module()
    params = SimpleNamespace()
    pairs = check_accuracy_lensing_pairs(pairs, inifile)
    standard = check_accuracy_run_output(
        "standard",
        params,
        lens_potential_cls=lensing_array_from_pairs(pairs, cols, use_new=True),
    )
    reference = check_accuracy_run_output(
        "reference",
        params,
        lens_potential_cls=lensing_array_from_pairs(pairs, cols, use_new=False),
    )
    return check_accuracy_rows_failed(filename, check_accuracy.compare_lensing(standard, reference))


def components_for_column(column):
    if "x" not in str(column):
        return None
    left, right = column.split("x", 1)
    if not left or not right:
        return None
    return left, right


def auto_column(component):
    return f"{component}x{component}"


def compare_cross_column(filename, pairs, cols, column, tolerance, *, max_l=None):
    components = components_for_column(column)
    if components is None:
        return False
    left, right = components
    if left == right:
        errors = column_error_series(pairs, cols, column)
    else:
        left_auto = auto_column(left)
        right_auto = auto_column(right)
        if left_auto not in cols or right_auto not in cols:
            return False
        col = cols.index(column)
        left_col = cols.index(left_auto)
        right_col = cols.index(right_auto)
        errors = [
            (
                ell,
                finite_normalized_cross_delta(new_row[col], orig_row[col], orig_row[left_col], orig_row[right_col]),
            )
            for ell, orig_row, new_row in pairs
        ]
    stop = None if max_l is None else max_l + 1
    return compare_error_ranges(
        filename,
        column,
        errors,
        [get_check_accuracy_module().RangeTolerance(0, stop, tolerance)],
        min_ell=2,
    )


def compare_lensing_cross_columns(filename, pairs, cols, inifile):
    if not inifile.bool("want_CMB", True):
        return False
    txp_tolerance = 1.4e-2 if is_boosted_reference_diff(filename) else 1e-2
    failed = False
    if "TxP" in cols:
        failed |= compare_cross_column(filename, pairs, cols, "TxP", txp_tolerance, max_l=100)
    if "ExP" in cols:
        failed |= compare_cross_column(filename, pairs, cols, "ExP", 2e-2, max_l=60)
    return failed


def should_compare_source_column(column, inifile):
    components = components_for_column(column)
    if components is None:
        return False
    left, right = components
    if not (left.startswith("W") or right.startswith("W")):
        return False
    if not inifile.bool("want_CMB", True) and (left in {"T", "E"} or right in {"T", "E"}):
        return False
    return True


def compare_source_window_columns(filename, pairs, cols, inifile):
    if is_boosted_reference_diff(filename):
        tolerance = 6e-2 if inifile.bool("Do21cm") else 7e-3
    else:
        tolerance = 5e-3
    failed = False
    for column in cols:
        if should_compare_source_column(column, inifile):
            failed |= compare_cross_column(filename, pairs, cols, column, tolerance)
    return failed


def check_accuracy_output_num_unequal(filename):
    is_lensed_cmb = fnmatch.fnmatch(filename, "*lensedCls.dat") or fnmatch.fnmatch(filename, "*lensedtotCls.dat")
    is_lensing = fnmatch.fnmatch(filename, "*lenspotentialCls.dat")
    is_scalar_cov = fnmatch.fnmatch(filename, "*scalarCovCls.dat")
    is_matter_power = fnmatch.fnmatch(filename, "*matterpower*.dat")
    is_transfer = is_transfer_output_file(filename)
    if not (is_lensed_cmb or is_lensing or is_scalar_cov or is_matter_power or is_transfer):
        return None

    orig_name = os.path.join(resolve_from_ini_dir(args.diff_to), filename)
    new_name = os.path.join(args.ini_dir, args.out_files_dir, filename)
    try:
        orig_cols, orig_rows = read_numeric_output_table(orig_name)
        new_cols, new_rows = read_numeric_output_table(new_name)
    except ValueError as error:
        if args.verbose_diff_output:
            printlog(f"could not parse numeric output in {filename}: {error}")
        return True
    if orig_cols != new_cols:
        if args.verbose_diff_output:
            printlog(f"num columns do not match in {filename}: reference and test headers differ")
        return True

    inifile = get_diff_ini(filename)
    if is_matter_power:
        return compare_check_accuracy_matter_power(filename, orig_rows, new_rows, new_cols, inifile)
    if is_transfer:
        return compare_check_accuracy_transfer(filename, orig_rows, new_rows, new_cols, inifile)

    pairs = align_numeric_rows_by_l(orig_rows, new_rows)
    if pairs is None:
        if args.verbose_diff_output:
            printlog(f"num rows do not overlap in {filename}")
        return True

    failed = False
    if is_lensed_cmb and should_compare_lensed_cmb_file(filename, inifile):
        failed |= compare_check_accuracy_cmb(filename, pairs, new_cols, inifile)
    if is_lensing:
        failed |= compare_check_accuracy_lensing_pp(filename, pairs, new_cols, inifile)
    if is_lensing or is_scalar_cov:
        failed |= compare_lensing_cross_columns(filename, pairs, new_cols, inifile)
    if is_scalar_cov:
        failed |= compare_source_window_columns(filename, pairs, new_cols, inifile)
    return failed


def num_unequal(filename, cmpFcn):
    """
    Check whether two files are numerically unequal for the given compare function.
    :param filename: The base name of the files to check.
    :param cmpFcn: The default comparison function. Can be overriden by the filetolmatrix.
    :return: True, when the files do not match, false else.
    """
    check_accuracy_mismatch = check_accuracy_output_num_unequal(filename)
    if check_accuracy_mismatch is not None:
        return check_accuracy_mismatch

    orig_name = os.path.join(resolve_from_ini_dir(args.diff_to), filename)
    with open(orig_name) as f:
        origMat = [[_x for _x in ln.split()] for ln in f]
        # Check if the first row has one more column, which is the #
        if len(origMat[0]) == len(origMat[1]) + 1:
            origBase = 1
            origMat[0] = origMat[0][1:]
        else:
            origBase = 0
    new_name = os.path.join(args.ini_dir, args.out_files_dir, filename)
    with open(new_name) as f:
        newMat = [[_x for _x in ln.split()] for ln in f]
        if len(newMat[0]) == len(newMat[1]) + 1:
            newBase = 1
            newMat[0] = newMat[0][1:]
        else:
            newBase = 0
    if newBase == 1:
        cols = [s[0] + "x" + s[1] if len(s) == 2 and s != "nu" else s for s in newMat[0]]
    else:
        cols = range(len(newMat[0]))

    orig_rows = origMat[origBase:]
    new_rows = newMat[newBase:]
    if len(orig_rows) != len(new_rows):
        aligned = align_mismatched_rows(filename, orig_rows, new_rows, cols)
        if aligned is None:
            if args.verbose_diff_output:
                printlog("num rows do not match in %s: %d != %d" % (filename, len(orig_rows), len(new_rows)))
            return True
        orig_rows, new_rows = aligned

    tolerances = get_tolerance_vector(filename, cols)
    row = 0
    col = 0
    try:
        if tolerances:
            inifilename = params_ini_for_output(filename)
            inifilename = os.path.join(args.ini_dir, args.out_files_dir, inifilename)
            if not os.path.exists(inifilename):
                if "sharp_cl_params" in inifilename:
                    inifile = make_ini_file()
                else:
                    printlog(f"ini filename does not exist: {inifilename}")
            else:
                try:
                    # The following split fails for *_transfer_out.* files where it not needed anyway.
                    inifile = make_ini_file()
                    inifile.readFile(inifilename)
                except OSError:
                    printlog(f"Could not open ini filename: {inifilename}")
            for o_row, n_row in zip(orig_rows, new_rows):
                row += 1
                if len(o_row) != len(n_row):
                    if args.verbose_diff_output:
                        printlog("num columns do not match in %s: %d != %d" % (filename, len(o_row), len(n_row)))
                    return True
                col = 0
                of_row = numeric_row(o_row)
                nf_row = numeric_row(n_row)
                oldrowdict = False
                newrowdict = False
                for o, n in zip(of_row, nf_row):
                    if isinstance(tolerances[col], Ignore):
                        pass
                    else:
                        cond, tols = tolerances[col]
                        compare_column = cond if isinstance(cond, bool) else cond(inifile)
                        if compare_column:
                            if isinstance(tols, float):
                                if not cmpFcn(o, n, scale_tolerance(tols)):
                                    if args.verbose_diff_output:
                                        printlog(
                                            'value mismatch at %d, %d ("%s") of %s: %s != %s'
                                            % (row, col + 1, cols[col], filename, o, n)
                                        )
                                    return True
                            elif not isinstance(tols, Ignore):
                                if not oldrowdict:
                                    oldrowdict = dict(zip(cols, of_row))
                                    newrowdict = dict(zip(cols, nf_row))
                                if isinstance(tols, list):
                                    cand = False
                                    for lim, rhs in tols:
                                        if lim < newrowdict["L"]:
                                            cand = rhs
                                        else:
                                            break
                                    if isinstance(cand, float):
                                        if not cmpFcn(o, n, scale_tolerance(cand)):
                                            if args.verbose_diff_output:
                                                printlog(
                                                    'value mismatch at %d, %d ("%s") of %s: %s != %s'
                                                    % (row, col + 1, cols[col], filename, o, n)
                                                )
                                            return True
                                    elif not isinstance(cand, (bool, Ignore)):
                                        if not cand(oldrowdict, newrowdict):
                                            if args.verbose_diff_output:
                                                printlog(
                                                    'value mismatch at %d, %d ("%s") of %s: %s != %s'
                                                    % (row, col + 1, cols[col], filename, o, n)
                                                )
                                            return True
                                else:
                                    if not tols(oldrowdict, newrowdict):
                                        if args.verbose_diff_output:
                                            printlog(
                                                'value mismatch at %d, %d ("%s") of %s: %s != %s'
                                                % (row, col + 1, cols[col], filename, o, n)
                                            )
                                        return True
                    col += 1
            return False
        else:
            #            if args.verbose_diff_output:
            #                printlog("Skipped file %s" % (filename))
            return False
    except ValueError as e:
        printlog("ValueError: '%s' at %d, %d in file: %s" % (e, row, col + 1, filename))
        return True


def customsplit(s):
    """
    Need to implement our own split, because for exponents of three digits
    the 'E' marking the exponent is dropped, which is not supported by python.
    :param s: The string to split.
    :return: An array containing the mantissa and the exponent, or the value, when no split was possible.
    """
    n = len(s)
    i = n - 1
    # Split the exponent from the string by looking for ['E']('+'|'-')D+
    while i > 4:
        if s[i] == "+" or s[i] == "-":
            return [s[0 : i - 1], s[i:n]]
        i -= 1
    return [s]


def textualcmp(o, n, tolerance):
    """
    Do a textual comparison for numbers whose exponent is zero or greater.
    The fortran code writes floating point values, with 5 significant digits
    after the comma and an exponent. I.e., for numbers with a positive
    exponent the usual comparison against a delta fails.
    :param o: The old value.
    :param n: The new value.
    :param tolerance: The allowed tolerance.
    :return: True, when |o - n| is greater then the tolerance allows, false else.
    """
    o_s = customsplit(o)
    n_s = customsplit(n)
    if len(o_s) > 1 and len(n_s) > 1:
        o_mantise = float(o_s[0])
        o_exp = int(o_s[1])
        n_mantise = float(n_s[0])
        n_exp = int(n_s[1])
        # Check without respect of the exponent, when that is greater zero.
        if 0 <= o_exp:
            if o_exp != n_exp:
                # Quit when exponent difference is significantly larger
                if abs(o_exp - n_exp) > 1:
                    return True
                if o_exp > n_exp:
                    o_mantise *= 10.0
                else:
                    n_mantise *= 10.0
            return math.fabs(float(o_mantise) - float(n_mantise)) >= tolerance
        return math.fabs(float(o_s[0] + "E" + o_s[1]) - float(n_s[0] + "E" + n_s[1])) >= tolerance
    # In all other cases do a numerical check
    return math.fabs(float(o) - float(n)) >= tolerance


def run_diff():
    printlog("Running diff_to...")
    if args.num_diff:
        defCmpFcn = lambda o, n, t: math.fabs(float(o) - float(n)) >= t
    else:
        defCmpFcn = normabs
    out_files_dir2 = resolve_from_ini_dir(args.diff_to)
    _, mismatch, errors = filecmp.cmpfiles(
        out_files_dir, out_files_dir2, list(set(list_files(out_files_dir)) | set(list_files(out_files_dir2)))
    )
    len_errors = len(errors)
    if len_errors and len_errors != 1 and errors[0] != args.diff_to:
        printlog("Missing/Extra files:")
        for err in errors:
            if err != args.diff_to:
                printlog("  " + err)
    if mismatch:
        numerical_mismatch = [f for f in mismatch if num_unequal(f, defCmpFcn)]
        if numerical_mismatch:
            printlog("Files do not match:")
            for err in numerical_mismatch:
                printlog("  " + err)
        len_num_mismatch = len(numerical_mismatch)
    else:
        len_num_mismatch = 0

    printlog("Done with %d numerical accuracy mismatches and %d extra/missing files" % (len_num_mismatch, len_errors))
    return 1 if len_errors > 0 or len_num_mismatch > 0 else 0


def main(argv=None):
    global args, prog, out_files_dir

    args = build_parser().parse_args(argv)
    args.ini_dir = os.path.abspath(args.ini_dir)
    args.base_settings = os.path.abspath(args.base_settings)
    prog = os.path.abspath(args.prog)

    os.makedirs(args.ini_dir, exist_ok=True)
    out_files_dir = os.path.join(args.ini_dir, args.out_files_dir)

    if args.clean and os.path.exists(out_files_dir):
        remove_tree(out_files_dir)
    os.makedirs(out_files_dir, exist_ok=True)

    override_dir = None
    try:
        if args.diff_to:
            return run_diff()

        inis, override_dir = getPreparedIniFiles()
        if args.no_run_test:
            return 0

        errors = 0
        files = output_file_num(out_files_dir)
        if files:
            printlog(f"Output directory is not empty (run with --clean to force delete): {out_files_dir}")
            return 1
        start = time.time()
        error_list = []
        for ini in inis:
            printlog(os.path.basename(ini) + "...")
            timing, output, return_code = runScript(ini)
            if return_code:
                printlog(f"error {return_code}")
                if output:
                    printlog(str(output).strip())
            nfiles = output_file_num(out_files_dir)
            if nfiles > files:
                msg = f"..OK, produced {nfiles - files} files in {timing:.2f}s"
            else:
                errors += 1
                error_list.append(os.path.basename(ini))
                msg = f"..no files in {timing:.2f}s"
            printlog(msg)
            files = nfiles
        printlog(f"Done, {errors} errors in {time.time() - start:.2f}s (outputs not checked yet)")
        if errors:
            printlog(f"Fails in : {error_list}")
        return 1 if errors else 0
    finally:
        if override_dir and os.path.exists(override_dir):
            remove_tree(override_dir)
        close_logfile()


if __name__ == "__main__":
    raise SystemExit(main())
