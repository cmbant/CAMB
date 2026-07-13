"""Make the lensed-CMB validation plot for the CAMB2 paper draft.

The plotted residual metric is the pointwise envelope

    max(|Delta TT|/TT, |Delta EE|/EE, |Delta TE|/sqrt(TT EE)).

For the paper figure the high-ell tail uses |Delta EE|/EE only for ell>3000,
since this is where foregrounds strongly dominate.

Running with no arguments reproduces the published figure
(fig:lensed-cl-validation in docs/CAMB2.tex):

    python docs/changelog/make_lensed_cls_validation_plot.py

That is, a CAMB reference using ``STRICT_REFERENCE_SETTINGS``
(AccuracyBoost=lSampleBoost=lAccuracyBoost=3 and min_l_logl_sampling=100000),
the plain CLASS ``cl_permille`` curve, and a CLASS ``cl_ref`` curve plotted over
the full L<=lmax range. By default both CLASS curves are computed to
L=lmax+4000 and then truncated;
all curves use Planck-era RECFAST, fixed Y_He, HMCode 2020 (no feedback), and
matched physical CMB lensing source k_max. The output is written to
docs/changelog/lensed_cls_validation/lensed_cls_validation_residuals_planck_recfast_cl_ref_ee_gt3000.pdf,
the file referenced by the paper.

The higher-accuracy CLASS curve is very expensive over the full range (tens of
minutes). For a quick plot, lower its plotted range or skip it, e.g.

    python docs/changelog/make_lensed_cls_validation_plot.py --class-ref-lmax 500
    python docs/changelog/make_lensed_cls_validation_plot.py --skip-class-ref

To change where the plot switches to EE-only residuals, use e.g.

    python docs/changelog/make_lensed_cls_validation_plot.py --ee-only-above 3500
    python docs/changelog/make_lensed_cls_validation_plot.py --no-ee-only-above
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

import camb
from camb import model
from camb.check_accuracy import (
    DEFAULT_NOISE_CONFIGS,
    STRICT_REFERENCE_SETTINGS,
    NoiseConfig,
    apply_accuracy_settings,
    make_cmb_noise_spectra,
)

try:
    from classy import Class
except ImportError as exc:  # pragma: no cover - exercised only when CLASS is absent
    raise SystemExit(
        "Could not import classy. Set PYTHONPATH to the local CLASS wrapper, e.g. "
        "PYTHONPATH=/tmp/classy-comparisons/classy-local:$PYTHONPATH"
    ) from exc


DEFAULT_OUTDIR = Path("docs/changelog/lensed_cls_validation")
DEFAULT_CLASS_WORKDIR = Path(os.environ.get("CLASS_WORKDIR", os.environ.get("TMPDIR", "/tmp") + "/classy-comparisons"))
DEFAULT_OUTPUT_STEM = "lensed_cls_validation_residuals_planck_recfast_cl_ref_ee_gt3000"
CLASS_PRECISION_IGNORED_KEYS = {"recfast_Nz0", "l_max_ur_ten"}
SPECTRUM_COLUMNS = {"TT": 0, "EE": 1, "BB": 2, "TE": 3}
CLASS_PRECISION_PROFILES = ("cl_permille", "cl_targeted", "cl_targeted_tight", "cl_ref")
RECOMBINATION_TAG = "planck_recfast"
CAMB_RECFAST_APPROX_MODEL = "planck"
CLASS_PLANCK_RECFAST_PARAMS = {
    "recombination": "RECFAST",
    "recfast_photoion_dependence": "Tmat",
    "recfast_Heswitch": 6,
    "recfast_fudge_He": 0.86,
    "recfast_Hswitch": 1,
    "recfast_fudge_H": 1.14,
    "recfast_delta_fudge_H": -0.015,
    "recfast_AGauss1": -0.14,
    "recfast_AGauss2": 0.079,
    "recfast_zGauss1": 7.28,
    "recfast_zGauss2": 6.73,
    "recfast_wGauss1": 0.18,
    "recfast_wGauss2": 0.33,
}
COSMOLOGY = {
    "H0": 67.5,
    "ombh2": 0.02237,
    "omch2": 0.12,
    "As": 2.1e-9,
    "ns": 0.965,
    "tau": 0.0544,
    "YHe": 0.245,
    "TCMB": 2.7255,
    "nnu": 3.044,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--lmax", type=int, default=4000)
    parser.add_argument(
        "--class-ref-lmax",
        type=int,
        default=None,
        help=(
            "Plotted-multipole limit for the higher-accuracy CLASS curve; it is masked above this. "
            "Defaults to --lmax (full range, matching the published figure). The higher-accuracy CLASS "
            "run is expensive, so pass a smaller value (e.g. 500) or --skip-class-ref for a quick plot."
        ),
    )
    parser.add_argument(
        "--class-ref-profile",
        choices=CLASS_PRECISION_PROFILES,
        default="cl_ref",
        help=(
            "CLASS precision profile for the higher-accuracy comparison curve. "
            "cl_targeted and cl_targeted_tight are moderate CMB/lensing profiles between "
            "cl_permille.pre and cl_ref.pre (default cl_ref, used for the published figure)."
        ),
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--class-workdir", type=Path, default=DEFAULT_CLASS_WORKDIR)
    parser.add_argument("--reference-lens-potential-accuracy", type=float, default=8.0)
    parser.add_argument("--reference-lens-output-margin", type=int, default=1000)
    parser.add_argument(
        "--class-lensing-margin",
        type=int,
        default=2000,
        help=(
            "Fallback CLASS lensed-spectrum margin above the plotted/reference maximum. "
            "Use the more specific CLASS margin options below for the plotted curves."
        ),
    )
    parser.add_argument(
        "--class-permille-lensing-margin",
        type=int,
        default=4000,
        help=(
            "Compute the plotted CLASS cl_permille curve to this many multipoles above the plotted maximum "
            "before truncating."
        ),
    )
    parser.add_argument(
        "--class-ref-lensing-margin",
        type=int,
        default=4000,
        help=(
            "Compute the higher-accuracy CLASS reference curve to this many multipoles above --class-ref-lmax. "
            "Default 4000, matching the published figure."
        ),
    )
    parser.add_argument(
        "--class-kmax-mode",
        choices=("boosted-camb", "precision-file"),
        default="boosted-camb",
        help=(
            "For CLASS runs, either override k_max_limber_over_l_max_scalars to match the boosted CAMB "
            "physical CMB lensing source kmax, or leave the precision file value unchanged"
        ),
    )
    parser.add_argument(
        "--class-k-max-limber-over-l-max-scalars",
        type=float,
        help="Explicit CLASS k_max_limber_over_l_max_scalars override; takes precedence over --class-kmax-mode",
    )
    parser.add_argument(
        "--class-k-max-tau0-over-l-max",
        type=float,
        help=(
            "Deprecated explicit CLASS primary-transfer k_max_tau0_over_l_max override. "
            "Only use for non-full-Limber checks; takes precedence over --class-kmax-mode."
        ),
    )
    parser.add_argument(
        "--mnu",
        type=float,
        default=0.0,
        help=(
            "Summed neutrino mass in eV. Default 0 for the baseline paper plot; "
            "use --mnu 0.06 for the standard minimal-mass variant."
        ),
    )
    parser.add_argument(
        "--num-massive-neutrinos",
        type=int,
        help="Number of massive neutrinos in CAMB. Defaults to 0 for mnu=0 and 1 for mnu>0.",
    )
    parser.add_argument(
        "--neutrino-hierarchy",
        choices=("degenerate", "normal", "inverted"),
        default="degenerate",
        help="CAMB neutrino hierarchy setting for mnu>0.",
    )
    parser.add_argument("--feedback-logt-agn", type=float, default=7.8)
    parser.add_argument(
        "--strict-camb-reference",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "use camb.check_accuracy.STRICT_REFERENCE_SETTINGS (AccuracyBoost=3) for the boosted CAMB "
            "reference (default; pass --no-strict-camb-reference for the lighter AccuracyBoost=2 reference)"
        ),
    )
    parser.add_argument(
        "--linear-lensing-only",
        action="store_true",
        help="Use linear matter power for all lensed spectra instead of HMCode nonlinear corrections",
    )
    parser.add_argument("--skip-class-ref", action="store_true", help="omit the higher-accuracy CLASS curve")
    parser.add_argument(
        "--output-stem",
        default=DEFAULT_OUTPUT_STEM,
        help="base filename for outputs (default writes the file referenced by docs/CAMB2.tex)",
    )
    parser.add_argument("--ymax", type=float, default=1.0e-2, help="upper limit for the residual-envelope axis")
    parser.add_argument(
        "--ee-only-above",
        type=int,
        default=3000,
        help=("Use only |Delta EE|/EE in the residual curve above this multipole. Default 3000 for the paper figure."),
    )
    parser.add_argument(
        "--no-ee-only-above",
        dest="ee_only_above",
        action="store_const",
        const=None,
        help="use the full TT/EE/TE residual envelope at all plotted multipoles",
    )
    parser.add_argument(
        "--historical-camb-cache",
        type=Path,
        help="optional cached CAMB case from another checkout to compare against the boosted reference",
    )
    parser.add_argument(
        "--historical-camb-label",
        default="CAMB origin/master default / boosted CAMB",
        help="legend label for --historical-camb-cache",
    )
    parser.add_argument("--force", action="store_true", help="rerun cases even if cached data exist")
    return parser


def find_class_root(class_workdir: Path) -> Path:
    roots = [path for path in class_workdir.glob("class_public-[0-9]*") if path.is_dir()]
    if not roots:
        raise FileNotFoundError(f"No CLASS source tree found under {class_workdir}")
    return sorted(roots)[-1]


def load_class_precision_file(path: Path) -> dict[str, float | int | str]:
    overrides: dict[str, float | int | str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if not value or key in CLASS_PRECISION_IGNORED_KEYS:
            continue
        try:
            parsed: float | int | str = float(value)
            if parsed.is_integer() and "." not in value and "e" not in value.lower():
                parsed = int(parsed)
        except ValueError:
            parsed = value
        overrides[key] = parsed
    return overrides


def class_precision_profiles(class_root: Path) -> dict[str, dict[str, float | int | str]]:
    cl_permille = load_class_precision_file(class_root / "cl_permille.pre")
    cl_ref = load_class_precision_file(class_root / "cl_ref.pre")
    cl_targeted = {
        **cl_permille,
        # Tighter k and perturbation sampling, but much lighter than cl_ref.pre.
        "k_step_sub": 0.025,
        "k_step_super": 0.0005,
        "tol_perturbations_integration": 3.0e-6,
        "perturbations_sampling_stepsize": 0.03,
        # Denser transfer/Bessel sampling at high l.
        "l_logstep": 1.05,
        "l_linstep": 25,
        # Use the more accurate lensed-Cl angular quadrature without the full cl_ref source settings.
        "l_switch_limber": 40.0,
        "accurate_lensing": 1,
        "num_mu_minus_lmax": 400,
        "delta_l_max": 1000,
        # Keep scalar sources to higher k before truncating tails.
        "transfer_neglect_delta_k_S_t0": 5.0,
        "transfer_neglect_delta_k_S_t1": 5.0,
        "transfer_neglect_delta_k_S_t2": 5.0,
        "transfer_neglect_delta_k_S_e": 5.0,
    }
    cl_targeted_tight = {
        **cl_targeted,
        # Ref-like scalar sampling/integration, while keeping the all-purpose cl_ref source tail lighter.
        "k_step_sub": 0.015,
        "k_step_super": 0.0002,
        "tol_perturbations_integration": 1.0e-6,
        "perturbations_sampling_stepsize": 0.015,
        "l_logstep": 1.035,
        "l_linstep": 15,
        "l_max_g": 35,
        "l_max_pol_g": 22,
        "l_max_ur": 40,
        "radiation_streaming_approximation": 2,
        "radiation_streaming_trigger_tau_over_tau_k": 160.0,
        "radiation_streaming_trigger_tau_c_over_tau": 60.0,
        "ur_fluid_approximation": 2,
        "ur_fluid_trigger_tau_over_tau_k": 45.0,
        "hyper_sampling_flat": 10.0,
        "hyper_flat_approximation_nu": 1.0e6,
        "num_mu_minus_lmax": 800,
        "transfer_neglect_delta_k_S_t0": 20.0,
        "transfer_neglect_delta_k_S_t1": 20.0,
        "transfer_neglect_delta_k_S_t2": 20.0,
        "transfer_neglect_delta_k_S_e": 20.0,
        "neglect_CMB_sources_below_visibility": 1.0e-20,
        "transfer_neglect_late_source": 1500.0,
    }
    return {
        "cl_permille": cl_permille,
        "cl_targeted": cl_targeted,
        "cl_targeted_tight": cl_targeted_tight,
        "cl_ref": cl_ref,
    }


def class_profile_label(profile: str) -> str:
    labels = {
        "cl_permille": "CLASS cl_permille",
        "cl_targeted": "CLASS cl_permille tuned",
        "cl_targeted_tight": "CLASS cl_permille boosted",
        "cl_ref": "CLASS cl_ref",
    }
    return labels[profile]


def make_camb_params(
    *,
    lmax: int,
    lens_potential_accuracy: float | None,
    lens_output_margin: int,
    boosted: bool,
    mnu: float,
    num_massive_neutrinos: int,
    neutrino_hierarchy: str,
    feedback_logt_agn: float | None = None,
    use_nonlinear_lensing: bool = True,
    strict_reference: bool = False,
) -> camb.CAMBparams:
    params = camb.CAMBparams()
    params.set_classes(recombination_model="Recfast")
    params.Recomb.set_params(recfast_approx_model=CAMB_RECFAST_APPROX_MODEL)
    params.set_cosmology(
        H0=COSMOLOGY["H0"],
        ombh2=COSMOLOGY["ombh2"],
        omch2=COSMOLOGY["omch2"],
        mnu=mnu,
        num_massive_neutrinos=num_massive_neutrinos,
        neutrino_hierarchy=neutrino_hierarchy,
        nnu=COSMOLOGY["nnu"],
        YHe=COSMOLOGY["YHe"],
        tau=COSMOLOGY["tau"],
        TCMB=COSMOLOGY["TCMB"],
    )
    params.InitPower.set_params(As=COSMOLOGY["As"], ns=COSMOLOGY["ns"])
    params.NonLinear = model.NonLinear_both if use_nonlinear_lensing else model.NonLinear_none
    if not use_nonlinear_lensing:
        params.NonLinearModel.set_params(halofit_version="mead2020")
    elif feedback_logt_agn is None:
        params.NonLinearModel.set_params(halofit_version="mead2020")
    else:
        params.NonLinearModel.set_params(halofit_version="mead2020_feedback", HMCode_logT_AGN=feedback_logt_agn)
    params.set_for_lmax(
        lmax,
        lens_potential_accuracy=lens_potential_accuracy,
        lens_output_margin=lens_output_margin,
        nonlinear=use_nonlinear_lensing,
    )
    if boosted:
        if strict_reference:
            apply_accuracy_settings(params, STRICT_REFERENCE_SETTINGS, boost_from_raw=True)
        else:
            params.set_accuracy(
                AccuracyBoost=2.0,
                lSampleBoost=2.0,
                lAccuracyBoost=2.0,
                min_l_logl_sampling=100000,
            )
            params.Accuracy.IntTolBoost = 2.0
    params.WantCls = True
    params.Want_CMB = True
    params.WantTransfer = False
    params.WantTensors = False
    params.WantDerivedParameters = False
    return params


def run_camb_case(params: camb.CAMBparams, lmax: int) -> dict[str, Any]:
    start_cpu = time.process_time()
    start_wall = time.perf_counter()
    results = camb.get_results(params)
    cpu_s = time.process_time() - start_cpu
    wall_s = time.perf_counter() - start_wall
    eta0 = float(results.conformal_time(0.0))
    return {
        "dls": results.get_lensed_scalar_cls(lmax=lmax, CMB_unit="muK"),
        "raw_cls": results.get_lensed_scalar_cls(lmax=lmax, CMB_unit="muK", raw_cl=True),
        "cpu_s": cpu_s,
        "wall_s": wall_s,
        "max_l": int(params.max_l),
        "max_eta_k": float(params.max_eta_k),
        "eta0": eta0,
        "k_max_1/Mpc": float(params.max_eta_k / eta0),
    }


def ensure_camb_lensing_metadata(payload: dict[str, Any], params: camb.CAMBparams) -> dict[str, Any]:
    if "eta0" not in payload or "k_max_1/Mpc" not in payload:
        background = camb.get_background(params)
        eta0 = float(background.conformal_time(0.0))
        payload["eta0"] = eta0
        payload["k_max_1/Mpc"] = float(payload["max_eta_k"] / eta0)
    return payload


def class_lensed_dls(params: dict[str, float | int | str], lmax: int) -> tuple[np.ndarray, np.ndarray, float, float]:
    cosmo = Class()
    cosmo.set(params)
    start_cpu = time.process_time()
    start_wall = time.perf_counter()
    try:
        cosmo.compute()
        cls = cosmo.lensed_cl(lmax)
    finally:
        cpu_s = time.process_time() - start_cpu
        wall_s = time.perf_counter() - start_wall
        cosmo.struct_cleanup()
        cosmo.empty()

    ell = np.arange(lmax + 1, dtype=np.float64)
    factor = ell * (ell + 1.0) / (2.0 * math.pi) * (COSMOLOGY["TCMB"] * 1.0e6) ** 2
    dls = np.zeros((lmax + 1, 4), dtype=np.float64)
    raw_cls = np.zeros((lmax + 1, 4), dtype=np.float64)
    for name, column in SPECTRUM_COLUMNS.items():
        key = name.lower()
        if key in cls:
            raw_cls[:, column] = cls[key][: lmax + 1] * (COSMOLOGY["TCMB"] * 1.0e6) ** 2
            dls[:, column] = cls[key][: lmax + 1] * factor
    return dls, raw_cls, cpu_s, wall_s


def make_class_params(
    *,
    lmax: int,
    precision_overrides: dict[str, float | int | str],
    neutrino_params: dict[str, float | int | str],
    use_nonlinear_lensing: bool = True,
) -> dict[str, float | int | str]:
    params: dict[str, float | int | str] = {
        "output": "tCl,pCl,lCl",
        "lensing": "yes",
        "non linear": "hmcode" if use_nonlinear_lensing else "none",
        "l_max_scalars": lmax,
        "H0": COSMOLOGY["H0"],
        "omega_b": COSMOLOGY["ombh2"],
        "omega_cdm": COSMOLOGY["omch2"],
        "Omega_k": 0.0,
        "A_s": COSMOLOGY["As"],
        "n_s": COSMOLOGY["ns"],
        "tau_reio": COSMOLOGY["tau"],
        "YHe": COSMOLOGY["YHe"],
        "T_cmb": COSMOLOGY["TCMB"],
    }
    params.update(neutrino_params)
    if use_nonlinear_lensing:
        params["hmcode_version"] = "2020"
    params.update(CLASS_PLANCK_RECFAST_PARAMS)
    params.update(precision_overrides)
    params.update(CLASS_PLANCK_RECFAST_PARAMS)
    params["l_max_scalars"] = lmax
    return params


def run_class_case(
    *,
    lmax: int,
    precision_overrides: dict[str, float | int | str],
    k_max_limber_over_l_max_scalars: float | None,
    k_max_tau0_over_l_max: float | None = None,
    neutrino_params: dict[str, float | int | str],
    use_nonlinear_lensing: bool = True,
) -> dict[str, Any]:
    params = make_class_params(
        lmax=lmax,
        precision_overrides=precision_overrides,
        neutrino_params=neutrino_params,
        use_nonlinear_lensing=use_nonlinear_lensing,
    )
    if k_max_limber_over_l_max_scalars is not None:
        params["k_max_limber_over_l_max_scalars"] = k_max_limber_over_l_max_scalars
    if k_max_tau0_over_l_max is not None:
        params["k_max_tau0_over_l_max"] = k_max_tau0_over_l_max
    dls, raw_cls, cpu_s, wall_s = class_lensed_dls(params, lmax)
    return {
        "dls": dls,
        "raw_cls": raw_cls,
        "cpu_s": cpu_s,
        "wall_s": wall_s,
        "class_params": params,
    }


def save_case(path: Path, payload: dict[str, Any]) -> None:
    np.savez_compressed(
        path,
        dls=payload["dls"],
        raw_cls=payload["raw_cls"],
        metadata=np.array(json.dumps({k: v for k, v in payload.items() if k not in {"dls", "raw_cls"}})),
    )


def load_case(path: Path) -> dict[str, Any]:
    data = np.load(path, allow_pickle=False)
    metadata = json.loads(str(data["metadata"]))
    return {"dls": data["dls"], "raw_cls": data["raw_cls"], **metadata}


def cached_case(path: Path, force: bool, runner: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    if path.exists() and not force:
        print(f"Using cached {path}")
        return load_case(path)
    print(f"Running {path.stem}...")
    payload = runner()
    save_case(path, payload)
    return payload


def kmax_tag(k_max_limber_over_l_max_scalars: float | None, k_max_tau0_over_l_max: float | None = None) -> str:
    if k_max_tau0_over_l_max is not None:
        return f"kmaxtau{k_max_tau0_over_l_max:.6g}".replace(".", "p")
    if k_max_limber_over_l_max_scalars is not None:
        return f"kmaxlimber{k_max_limber_over_l_max_scalars:.6g}".replace(".", "p")
    return "precision"


def float_tag(value: float) -> str:
    return f"{value:.6g}".replace("-", "m").replace(".", "p")


def mnu_tag(mnu: float) -> str:
    return "" if mnu == 0 else f"_mnu{float_tag(mnu)}"


def class_neutrino_params_from_camb(params: camb.CAMBparams) -> dict[str, float | int | str]:
    if params.nu_mass_eigenstates == 0:
        return {"N_ur": params.num_nu_massless, "N_ncdm": 0}
    t_std = (4 / 11) ** (1 / 3)
    n = params.nu_mass_eigenstates
    return {
        "N_ur": params.num_nu_massless,
        "N_ncdm": n,
        "deg_ncdm": ",".join(str(params.nu_mass_numbers[i]) for i in range(n)),
        "T_ncdm": ",".join(
            str(t_std * (params.nu_mass_degeneracies[i] / params.nu_mass_numbers[i]) ** 0.25) for i in range(n)
        ),
        "omega_ncdm": ",".join(str(params.omnuh2 * params.nu_mass_fractions[i]) for i in range(n)),
    }


def residual_components(candidate: np.ndarray, reference: np.ndarray) -> dict[str, np.ndarray]:
    with np.errstate(divide="ignore", invalid="ignore"):
        tt = (candidate[:, 0] - reference[:, 0]) / reference[:, 0]
        ee = (candidate[:, 1] - reference[:, 1]) / reference[:, 1]
        te = (candidate[:, 3] - reference[:, 3]) / np.sqrt(np.maximum(reference[:, 0] * reference[:, 1], 0.0))
    return {"TT": tt, "EE": ee, "TE": te}


def residual_envelope(
    candidate: np.ndarray, reference: np.ndarray, lmax: int, ee_only_above: int | None = None
) -> np.ndarray:
    envelope = np.full(lmax + 1, np.nan)
    max_l = min(lmax, candidate.shape[0] - 1, reference.shape[0] - 1)
    components = residual_components(candidate[: max_l + 1], reference[: max_l + 1])
    abs_components = {name: np.abs(values) for name, values in components.items()}
    if ee_only_above is not None:
        start = max(0, ee_only_above + 1)
        abs_components["TT"][start:] = np.nan
        abs_components["TE"][start:] = np.nan
    stacked = np.vstack(list(abs_components.values()))
    finite_columns = np.any(np.isfinite(stacked), axis=0)
    envelope[: max_l + 1][finite_columns] = np.nanmax(stacked[:, finite_columns], axis=0)
    envelope[:2] = np.nan
    return envelope


def mask_above(values: np.ndarray, lmax: int) -> np.ndarray:
    masked = values.copy()
    if masked.shape[0] > lmax + 1:
        masked[lmax + 1 :] = np.nan
    return masked


def amplitude_fisher_per_ell(dls: np.ndarray, config: NoiseConfig, lmax: int, fields: tuple[str, ...]) -> np.ndarray:
    ell = np.arange(lmax + 1)
    noise = make_cmb_noise_spectra(ell, config)
    fisher = np.zeros(lmax + 1)
    for multipole in range(max(2, config.lmin), lmax + 1):
        if fields == ("E",):
            signal = dls[multipole, 1]
            total = signal + noise["E"][multipole]
            if total <= 0:
                continue
            fisher[multipole] = config.fsky * (2 * multipole + 1) * (signal / total) ** 2 / 2.0
        elif fields == ("T", "E"):
            signal = np.array(
                [
                    [dls[multipole, 0], dls[multipole, 3]],
                    [dls[multipole, 3], dls[multipole, 1]],
                ]
            )
            total = signal + np.diag([noise["T"][multipole], noise["E"][multipole]])
            try:
                inv_total = np.linalg.inv(total)
            except np.linalg.LinAlgError:
                continue
            fisher[multipole] = (
                config.fsky * (2 * multipole + 1) * np.trace(inv_total @ signal @ inv_total @ signal) / 2.0
            )
        else:
            raise ValueError(f"Unsupported fields for amplitude Fisher: {fields}")
    return fisher


def binned_amplitude_sigma(
    dls: np.ndarray,
    config: NoiseConfig,
    lmax: int,
    fractional_bin_width: float = 0.5,
    fields: tuple[str, ...] = ("T", "E"),
    ee_only_above: int | None = None,
) -> np.ndarray:
    fisher_per_ell = amplitude_fisher_per_ell(dls, config, lmax, fields)
    if ee_only_above is not None:
        if fields != ("T", "E"):
            raise ValueError("--ee-only-above only applies to the TT+EE amplitude band")
        fisher_ee = amplitude_fisher_per_ell(dls, config, lmax, ("E",))
        fisher_per_ell[ee_only_above + 1 :] = fisher_ee[ee_only_above + 1 :]
    cumulative = np.concatenate(([0.0], np.cumsum(fisher_per_ell)))
    fisher = np.zeros(lmax + 1)
    for multipole in range(max(2, config.lmin), lmax + 1):
        half_width = fractional_bin_width * multipole / 2.0
        lo = max(config.lmin, int(math.ceil(multipole - half_width)))
        hi = min(lmax, int(math.floor(multipole + half_width)))
        if hi >= lo:
            fisher[multipole] = cumulative[hi + 1] - cumulative[lo]
    sigma = np.full(lmax + 1, np.nan)
    mask = fisher > 0
    sigma[mask] = 1.0 / np.sqrt(fisher[mask])
    return sigma


def summarize_curve(values: np.ndarray, ranges: dict[str, tuple[int, int]]) -> dict[str, dict[str, float | int | None]]:
    summary = {}
    for name, (lo, hi) in ranges.items():
        hi = min(hi, values.shape[0] - 1)
        segment = values[lo : hi + 1]
        if not np.any(np.isfinite(segment)):
            summary[name] = {"max": None, "ell": -1}
            continue
        offset = int(np.nanargmax(segment))
        summary[name] = {"max": float(segment[offset]), "ell": lo + offset}
    return summary


def make_plot(
    path: Path,
    ell: np.ndarray,
    curves: dict[str, np.ndarray],
    sigma_a: np.ndarray,
    summaries: dict[str, dict[str, dict[str, float | int]]],
    y_max: float = 1.0e-2,
    ee_only_above: int | None = None,
) -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    def color_for(label: str) -> str:
        if label.startswith(("CAMB origin/master", "CAMB v1")):
            return "#8a5a00"
        if label.startswith("CAMB default"):
            return "#111111"
        if label.startswith("CLASS cl_permille"):
            return "#2f8f46" if "boosted /" in label or "tuned /" in label else "#1b6ca8"
        if label.startswith("CLASS cl_ref"):
            return "#2f8f46"
        if label.startswith("HMCode feedback"):
            return "#b23a48"
        if "boosted" in label or "tuned" in label:
            return "#2f8f46"
        return "0.2"

    def linestyle_for(label: str):
        if label.startswith(("CAMB origin/master", "CAMB v1")):
            return (0, (1.5, 1.5))
        if label.startswith("CAMB default"):
            return "-"
        if label.startswith("CLASS cl_permille") and "boosted /" not in label and "tuned /" not in label:
            return (0, (4, 2))
        if "boosted" in label or "tuned" in label or label.startswith("CLASS cl_ref"):
            return (0, (4, 2))
        if label.startswith("HMCode feedback"):
            return "-."
        return "-"

    fig, ax = plt.subplots(figsize=(7.15, 4.15))
    ymin = 0.7e-4
    ax.fill_between(
        ell,
        ymin,
        sigma_a,
        where=np.isfinite(sigma_a) & (sigma_a > ymin),
        color="0.86",
        alpha=0.85,
        linewidth=0,
        zorder=0,
    )
    ax.axhline(1.0e-3, color="0.35", linewidth=0.9, linestyle=":")
    if ee_only_above is not None:
        ax.axvline(ee_only_above, color="0.55", linewidth=0.8, linestyle="--", zorder=0.5)
    for label, values in curves.items():
        ax.plot(
            ell,
            values,
            color=color_for(label),
            linestyle=linestyle_for(label),
            linewidth=1.25,
            label=label,
        )
    ax.set_xscale("function", functions=(np.sqrt, np.square))
    ax.set_yscale("log")
    ax.set_xlim(2, ell[-1])
    ax.set_ylim(ymin, y_max)
    xticks = np.array([2, 10, 30, 100, 300, 1000, 2000, 3000, 4000, 6000])
    xticks = xticks[(xticks >= 2) & (xticks <= ell[-1])]
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(int(tick)) for tick in xticks])
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"lensed $C_\ell$ fractional difference")
    ax.yaxis.grid(True, which="major", color="0.88", linewidth=0.6)
    ax.yaxis.grid(True, which="minor", color="0.93", linewidth=0.35)
    for tick in xticks:
        if tick not in (1000, 2000):
            ax.axvline(tick, color="0.88", linewidth=0.6, zorder=0.4)
    legend = ax.legend(
        loc="upper left",
        frameon=True,
        framealpha=0.96,
        borderpad=0.6,
    )
    fig.canvas.draw()
    legend_bbox = legend.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transAxes.inverted())
    ax.text(
        0.018,
        max(0.05, legend_bbox.y0 - 0.006),
        (
            r"$\max(TT,EE,TE)$ residual"
            if ee_only_above is None
            else rf"$\max(TT,EE,TE)$ for $\ell\leq {ee_only_above}$; EE only at $\ell>{ee_only_above}$"
        ),
        transform=ax.transAxes,
        color="0.25",
        fontsize=8,
        va="top",
    )
    fig.tight_layout(pad=0.6)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)


def main() -> None:
    args = build_parser().parse_args()
    if args.lmax <= 2:
        raise ValueError("--lmax must be greater than 2")
    if args.ee_only_above is not None and args.ee_only_above < 2:
        raise ValueError("--ee-only-above must be at least 2")
    if args.mnu < 0:
        raise ValueError("--mnu must be non-negative")
    if args.num_massive_neutrinos is None:
        args.num_massive_neutrinos = 0 if args.mnu == 0 else 1
    if args.mnu == 0 and args.num_massive_neutrinos != 0:
        raise ValueError("--num-massive-neutrinos must be 0 when --mnu is 0")
    if args.mnu > 0 and args.num_massive_neutrinos <= 0:
        raise ValueError("--num-massive-neutrinos must be positive when --mnu is non-zero")
    if args.class_ref_lmax is None:
        args.class_ref_lmax = args.lmax
    if args.class_ref_lmax > args.lmax:
        raise ValueError("--class-ref-lmax cannot exceed --lmax")
    if args.class_k_max_limber_over_l_max_scalars is not None and args.class_k_max_tau0_over_l_max is not None:
        raise ValueError(
            "--class-k-max-limber-over-l-max-scalars and --class-k-max-tau0-over-l-max are mutually exclusive"
        )
    if args.output_stem == DEFAULT_OUTPUT_STEM and args.mnu != 0:
        args.output_stem = f"{args.output_stem}{mnu_tag(args.mnu)}"

    args.outdir.mkdir(parents=True, exist_ok=True)
    class_root = find_class_root(args.class_workdir)
    class_version = class_root.name.removeprefix("class_public-")
    precision = class_precision_profiles(class_root)
    class_ref_lensing_margin = (
        args.class_lensing_margin if args.class_ref_lensing_margin is None else args.class_ref_lensing_margin
    )
    class_permille_compute_lmax = args.lmax + args.class_permille_lensing_margin
    class_ref_compute_lmax = args.class_ref_lmax + class_ref_lensing_margin
    lensing_tag = "linear_lensing" if args.linear_lensing_only else "nonlinear_lensing"
    cache_tag = RECOMBINATION_TAG if not args.linear_lensing_only else f"{RECOMBINATION_TAG}_{lensing_tag}"
    cache_tag = f"{cache_tag}{mnu_tag(args.mnu)}"
    reference_tag = "strictref" if args.strict_camb_reference else "boost2ref"
    camb_reference_settings = (
        STRICT_REFERENCE_SETTINGS
        if args.strict_camb_reference
        else {
            "AccuracyBoost": 2.0,
            "lSampleBoost": 2.0,
            "lAccuracyBoost": 2.0,
            "IntTolBoost": 2.0,
            "min_l_logl_sampling": 100000,
        }
    )
    baseline_nonlinear = (
        "linear matter power only" if args.linear_lensing_only else "CAMB/CLASS standard HMCode 2020 with no feedback"
    )
    neutrino_model = {
        "mnu": args.mnu,
        "num_massive_neutrinos": args.num_massive_neutrinos,
        "neutrino_hierarchy": args.neutrino_hierarchy,
    }

    print(f"CAMB {camb.__version__}: {Path(camb.__file__).resolve()}")
    print(f"CLASS {class_version}: {class_root}")

    default_params = make_camb_params(
        lmax=args.lmax,
        lens_potential_accuracy=None,
        lens_output_margin=200,
        boosted=False,
        **neutrino_model,
        use_nonlinear_lensing=not args.linear_lensing_only,
    )
    default = cached_case(
        args.outdir / f"camb_default_{cache_tag}.npz",
        args.force,
        lambda: run_camb_case(default_params, args.lmax),
    )
    ensure_camb_lensing_metadata(default, default_params)
    boosted_params = make_camb_params(
        lmax=args.lmax,
        lens_potential_accuracy=args.reference_lens_potential_accuracy,
        lens_output_margin=args.reference_lens_output_margin,
        boosted=True,
        **neutrino_model,
        use_nonlinear_lensing=not args.linear_lensing_only,
        strict_reference=args.strict_camb_reference,
    )
    boosted = cached_case(
        args.outdir / f"camb_{reference_tag}_{cache_tag}.npz",
        args.force,
        lambda: run_camb_case(boosted_params, args.lmax),
    )
    ensure_camb_lensing_metadata(boosted, boosted_params)
    class_neutrino_params = class_neutrino_params_from_camb(default_params)
    class_permille_kmax_tau0 = args.class_k_max_tau0_over_l_max
    class_ref_kmax_tau0 = args.class_k_max_tau0_over_l_max
    if args.class_k_max_tau0_over_l_max is not None:
        class_permille_kmax_limber = None
        class_ref_kmax_limber = None
    elif args.class_k_max_limber_over_l_max_scalars is not None:
        class_permille_kmax_limber = args.class_k_max_limber_over_l_max_scalars
        class_ref_kmax_limber = args.class_k_max_limber_over_l_max_scalars
    elif args.class_kmax_mode == "boosted-camb":
        class_permille_kmax_limber = float(boosted["k_max_1/Mpc"]) / float(class_permille_compute_lmax)
        class_ref_kmax_limber = float(boosted["k_max_1/Mpc"]) / float(class_ref_compute_lmax)
    else:
        class_permille_kmax_limber = None
        class_ref_kmax_limber = None
    class_permille_kmax_tag = kmax_tag(class_permille_kmax_limber, class_permille_kmax_tau0)
    class_ref_kmax_tag = kmax_tag(class_ref_kmax_limber, class_ref_kmax_tau0)
    feedback = None
    if not args.linear_lensing_only:
        feedback = cached_case(
            args.outdir / f"camb_hmcode_feedback_{cache_tag}_logt{float_tag(args.feedback_logt_agn)}.npz",
            args.force,
            lambda: run_camb_case(
                make_camb_params(
                    lmax=args.lmax,
                    lens_potential_accuracy=args.reference_lens_potential_accuracy,
                    lens_output_margin=args.reference_lens_output_margin,
                    boosted=True,
                    **neutrino_model,
                    feedback_logt_agn=args.feedback_logt_agn,
                ),
                args.lmax,
            ),
        )
    class_permille = cached_case(
        args.outdir / f"class_cl_permille_{cache_tag}_lmax{class_permille_compute_lmax}_{class_permille_kmax_tag}.npz",
        args.force,
        lambda: run_class_case(
            lmax=class_permille_compute_lmax,
            precision_overrides=precision["cl_permille"],
            k_max_limber_over_l_max_scalars=class_permille_kmax_limber,
            k_max_tau0_over_l_max=class_permille_kmax_tau0,
            neutrino_params=class_neutrino_params,
            use_nonlinear_lensing=not args.linear_lensing_only,
        ),
    )
    class_ref = None
    if not args.skip_class_ref:
        class_ref = cached_case(
            args.outdir
            / (
                f"class_{args.class_ref_profile}_{cache_tag}_compute_lmax{class_ref_compute_lmax}"
                f"_plot_lmax{args.class_ref_lmax}_{class_ref_kmax_tag}.npz"
            ),
            args.force,
            lambda: run_class_case(
                lmax=class_ref_compute_lmax,
                precision_overrides=precision[args.class_ref_profile],
                k_max_limber_over_l_max_scalars=class_ref_kmax_limber,
                k_max_tau0_over_l_max=class_ref_kmax_tau0,
                neutrino_params=class_neutrino_params,
                use_nonlinear_lensing=not args.linear_lensing_only,
            ),
        )

    ell = np.arange(args.lmax + 1)
    curves = {
        "CAMB default / boosted": residual_envelope(default["dls"], boosted["dls"], args.lmax, args.ee_only_above),
        "CLASS cl_permille / boosted CAMB": residual_envelope(
            class_permille["dls"], boosted["dls"], args.lmax, args.ee_only_above
        ),
    }
    historical_camb = None
    if args.historical_camb_cache is not None:
        historical_camb = load_case(args.historical_camb_cache)
        curves[args.historical_camb_label] = residual_envelope(
            historical_camb["dls"], boosted["dls"], args.lmax, args.ee_only_above
        )
    if class_ref is not None:
        curves[
            (
                f"{class_profile_label(args.class_ref_profile)} / boosted CAMB"
                if args.class_ref_lmax == args.lmax
                else f"{class_profile_label(args.class_ref_profile)} (L<={args.class_ref_lmax}) / boosted CAMB"
            )
        ] = mask_above(
            residual_envelope(class_ref["dls"], boosted["dls"], args.lmax, args.ee_only_above), args.class_ref_lmax
        )
    if feedback is not None:
        curves["HMCode feedback / no feedback"] = residual_envelope(
            feedback["dls"], boosted["dls"], args.lmax, args.ee_only_above
        )
    ranges = {
        "low": (2, 599),
        "interior": (600, 3499),
        "tail": (3500, args.lmax),
        "all": (2, args.lmax),
    }
    summaries = {label: summarize_curve(values, ranges) for label, values in curves.items()}
    so_defaults = DEFAULT_NOISE_CONFIGS["so"]
    noise_config = NoiseConfig(
        name="so",
        noise_muK_arcmin_T=so_defaults["noise_muK_arcmin_T"],
        noise_muK_arcmin_P=so_defaults["noise_muK_arcmin_P"],
        fwhm_arcmin=so_defaults["fwhm_arcmin"],
        fsky=so_defaults["fsky"],
        lmin=2,
        lmax=args.lmax,
        fields=("T", "E"),
    )
    amplitude_error_bin_width = 0.5
    sigma_a = binned_amplitude_sigma(
        boosted["dls"], noise_config, args.lmax, amplitude_error_bin_width, ee_only_above=args.ee_only_above
    )

    plot_path = args.outdir / f"{args.output_stem}.pdf"
    make_plot(plot_path, ell, curves, sigma_a, summaries, args.ymax, args.ee_only_above)

    metadata = {
        "plot": str(plot_path),
        "png": str(plot_path.with_suffix(".png")),
        "lmax": args.lmax,
        "class_ref_lmax": args.class_ref_lmax,
        "camb_version": camb.__version__,
        "camb_path": str(Path(camb.__file__).resolve()),
        "class_version": class_version,
        "class_root": str(class_root),
        "class_ref_profile": args.class_ref_profile,
        "class_ref_profile_label": class_profile_label(args.class_ref_profile),
        "class_precision_overrides": precision[args.class_ref_profile],
        "lensing_power": "linear" if args.linear_lensing_only else "hmcode_nonlinear",
        "skipped_class_ref": args.skip_class_ref,
        "ymax": args.ymax,
        "ee_only_above": args.ee_only_above,
        "residual_envelope": (
            "max(|Delta TT|/TT, |Delta EE|/EE, |Delta TE|/sqrt(TT EE))"
            if args.ee_only_above is None
            else (
                "max(|Delta TT|/TT, |Delta EE|/EE, |Delta TE|/sqrt(TT EE)) for ell <= ee_only_above; "
                "|Delta EE|/EE for ell > ee_only_above"
            )
        ),
        "camb_reference_settings": camb_reference_settings,
        "class_lensing_margin": args.class_lensing_margin,
        "class_permille_lensing_margin": args.class_permille_lensing_margin,
        "class_ref_lensing_margin": class_ref_lensing_margin,
        "class_permille_compute_lmax": class_permille_compute_lmax,
        "class_ref_compute_lmax": class_ref_compute_lmax,
        "recombination": {
            "tag": RECOMBINATION_TAG,
            "camb_model": "Recfast",
            "camb_recfast_approx_model": CAMB_RECFAST_APPROX_MODEL,
            "class_params": CLASS_PLANCK_RECFAST_PARAMS,
        },
        "cosmology": {**COSMOLOGY, **neutrino_model},
        "class_neutrino_params": class_neutrino_params,
        "baseline_nonlinear": baseline_nonlinear,
        "boosted_reference": {
            **camb_reference_settings,
            "lens_potential_accuracy": args.reference_lens_potential_accuracy,
            "lens_output_margin": args.reference_lens_output_margin,
            "max_eta_k": boosted["max_eta_k"],
            "eta0": boosted["eta0"],
            "k_max_1/Mpc": boosted["k_max_1/Mpc"],
            "max_l": boosted["max_l"],
        },
        "class_permille_k_max_limber_over_l_max_scalars": class_permille_kmax_limber,
        "class_ref_k_max_limber_over_l_max_scalars": class_ref_kmax_limber,
        "class_permille_k_max_tau0_over_l_max": class_permille_kmax_tau0,
        "class_ref_k_max_tau0_over_l_max": class_ref_kmax_tau0,
        "class_matched_max_eta_k": boosted["max_eta_k"] if args.class_kmax_mode == "boosted-camb" else None,
        "class_matched_k_max_1/Mpc": boosted["k_max_1/Mpc"] if args.class_kmax_mode == "boosted-camb" else None,
        "class_kmax_mode": args.class_kmax_mode,
        "standard_reference_noise": noise_config.__dict__,
        "amplitude_error_band": {
            "kind": "single-bin spectral-amplitude 1 sigma",
            "fields": (
                "TT+EE at all plotted multipoles"
                if args.ee_only_above is None
                else (
                    f"bin-summed TT+EE Fisher contributions for ell <= {args.ee_only_above}; "
                    f"EE-only contributions for ell > {args.ee_only_above}"
                )
            ),
            "fractional_bin_width": amplitude_error_bin_width,
            "bin_definition": "for each plotted L, use integer multipoles in [ceil(0.75 L), floor(1.25 L)]",
        },
        "summaries": summaries,
        "run_times": {
            "camb_default": {"cpu_s": default["cpu_s"], "wall_s": default["wall_s"]},
            "camb_boosted_reference": {"cpu_s": boosted["cpu_s"], "wall_s": boosted["wall_s"]},
            "class_cl_permille": {"cpu_s": class_permille["cpu_s"], "wall_s": class_permille["wall_s"]},
        },
        "note": "CLASS lensed spectra are computed above the plotted maximum and then truncated, to reduce high-ell "
        "lensing-convolution edge effects.",
    }
    if feedback is not None:
        metadata["feedback_curve"] = {
            "halofit_version": "mead2020_feedback",
            "HMCode_logT_AGN": args.feedback_logt_agn,
        }
        metadata["run_times"]["camb_hmcode_feedback"] = {"cpu_s": feedback["cpu_s"], "wall_s": feedback["wall_s"]}
    if historical_camb is not None:
        metadata["historical_camb_curve"] = {
            "label": args.historical_camb_label,
            "cache": str(args.historical_camb_cache),
            "metadata": {key: value for key, value in historical_camb.items() if key not in {"dls", "raw_cls"}},
        }
        metadata["run_times"]["historical_camb"] = {
            "cpu_s": historical_camb.get("cpu_s"),
            "wall_s": historical_camb.get("wall_s"),
        }
    if class_ref is not None:
        metadata["run_times"][
            f"class_{args.class_ref_profile}_compute_lmax{class_ref_compute_lmax}_{class_ref_kmax_tag}"
        ] = {
            "cpu_s": class_ref["cpu_s"],
            "wall_s": class_ref["wall_s"],
        }
        metadata["note"] += " The CLASS reference curve is masked above class_ref_lmax."
    json_path = args.outdir / f"{args.output_stem}.json"
    json_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {plot_path}")
    print(f"Wrote {plot_path.with_suffix('.png')}")
    print(f"Wrote {json_path}")


if __name__ == "__main__":
    main()
