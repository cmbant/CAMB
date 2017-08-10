from __future__ import print_function
import os
import subprocess
import argparse
import fnmatch
import time
import sys
import shutil
import copy
from iniFile import iniFile

parser = argparse.ArgumentParser(description='Run CAMB tests')
parser.add_argument('ini_dir', help='ini file directory')
parser.add_argument('--make_ini', action='store_true', help='if set, output ini files to ini_dir')
parser.add_argument('--out_files_dir', default='test_outputs', help='output files directory')
parser.add_argument('--base_settings', default='params.ini', help='settings to include as defaults for all combinations')
parser.add_argument('--no_run_test', action='store_true', help='don''t run tests on files')
parser.add_argument('--prog', default='./camb', help='executable to run')
parser.add_argument('--clean', action='store_true', help='delete output dir before run')
parser.add_argument('--diff_to', help='output directory to compare to, e.g. test_outputs2')
parser.add_argument('--diff_tolerance', type=float, help='the tolerance for the numerical diff when no explicit ' +
                                                         'diff is given [default: 1e-4]', default=1e-4)
parser.add_argument('--verbose_diff_output', action='store_true', help='during diff_to print more error messages')
parser.add_argument('--num_diff', action='store_true', help='during diff_to print more error messages')
parser.add_argument('--no_sources', action='store_true', help='turn off CAMB sources (counts/lensing/21cm) tests')
parser.add_argument('--no_de', action='store_true', help='don''t run dark energy tests')
parser.add_argument('--max_tests', type=int, help='maximum tests to run (for quick testing of pipeline)')



args = parser.parse_args()

logfile = None

def printlog(text):
    global logfile
    print(text)
    if logfile is None:
        logfile = open(os.path.join(args.ini_dir,'test_results.log'), 'a')
    logfile.write(text+"\n")

# The tolerance matrix gives the tolerances for comparing two values of the actual results with
# results given in a diff_to. Filename globing is supported by fnmatch. The first glob matching
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
    pass

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
    oc1c1 = old[c1 + 'x' + c1]
    oc2c2 = old[c2 + 'x' + c2]
    # Skip the test when exactly one variable is negative, but not both.
    if (oc1c1 < 0 or oc2c2 < 0) and (oc1c1 >= 0 or oc2c2 >= 0):
        return True
    res = math.fabs(new[c1 + 'x' + c2] - old[c1 + 'x' + c2]) / math.sqrt(oc1c1 * oc2c2) < tol
    if args.verbose_diff_output and not res:
        printlog("diffnsqrt: |%g - %g|/sqrt(%g * %g) = %g > %g" % (new[c1 + 'x' + c2], old[c1 + 'x' + c2], oc1c1, oc2c2,
                                                                math.fabs(new[c1 + 'x' + c2] - old[c1 + 'x' + c2]) / math.sqrt(oc1c1 * oc2c2),
                                                                tol))
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
        printlog("normabs: |%g - %g| / |%g| = %g > %g" % (o, n, o, math.fabs(o - n) / math.fabs(o) if o != 0.0 else math.fabs(o - n), tol))
    return res

def wantCMBTandlmaxscalarge2000(ini):
    """
    Return true when want_CMB is set in the ini file and l_max_scalar is >= 2000.
    :param ini: The dictionary all inifile settings.
    :return: True, when want_CMB and l_max_scalar >= 2000, false else.
    """
    return ini.int("l_max_scalar") >= 2000 and ini.bool("want_CMB")

def wantCMBT(ini):
    """
    Return true when want_CMB is set.
    :param ini: The dictionary all inifile settings.
    :return: True, when want_CMB is set.
    """
    return ini.bool("want_CMB")

# A short cut for lensedCls and lenspotentialCls files.
coltol1 = ColTol({"L": Ignore(),
                  "TxT": (wantCMBT ,
                          [(0, 3e-3),
                           (600, 1e-3),
                           (2500, 3e-3),
                           (6000, 0.02)]),
                  "ExE": (wantCMBT,
                          [(0, 3e-3),
                           (600, 1e-3),
                           (2500, 3e-3),
                           (6000, 0.02),
                           (8000, 0.1)]),
                  "BxB": (wantCMBTandlmaxscalarge2000,
                          [(0, 5e-3),
                           (1000, 1e-2),
                           (6000, 0.02),
                           (8000, 0.1)]),
                  "TxE": (wantCMBT,
                          [(0, lambda o, n: diffnsqrt(o, n, 3e-3, 'T', 'E')),
                           (600, lambda o, n: diffnsqrt(o, n, 1e-3, 'T', 'E')),
                           (2500, lambda o, n: diffnsqrt(o, n, 3e-3, 'T', 'E')),
                           (6000, lambda o, n: diffnsqrt(o, n, 3e-2, 'T', 'E'))]),
                  "PxP": (True,
                          [(0, 5e-3),
                           (1000, 1e-2),
                           (6000, 0.02)]),
                  "TxP": (wantCMBT,
                          [(0, lambda o, n: diffnsqrt(o, n, 0.01, 'T', 'P')),
                           (100, Ignore())]),
                  "ExP": (wantCMBT,
                          [(0, lambda o, n: diffnsqrt(o, n, 0.02, 'E', 'P')),
                           (60, Ignore())]),
                  "TxW1": (wantCMBT, lambda o, n: diffnsqrt(o, n, 5e-3, 'T', 'W1')),
                  "ExW1": (wantCMBT, lambda o, n: diffnsqrt(o, n, 5e-3, 'E', 'W1')),
                  "PxW1": (True, lambda o, n: diffnsqrt(o, n, 5e-3, 'P', 'W1')),
                  "W1xT": (wantCMBT, lambda o, n: diffnsqrt(o, n, 5e-3, 'W1', 'T')),
                  "W1xE": (wantCMBT, lambda o, n: diffnsqrt(o, n, 5e-3, 'W1', 'E')),
                  "W1xP": (True, lambda o, n: diffnsqrt(o, n, 5e-3, 'W1', 'P')),
                  "W1xW1": (True, 5e-3),
                  "PxW2": (True, lambda o, n: diffnsqrt(o, n, 5e-3, 'P', 'W2')),
                  "W1xW2": (True, lambda o, n: diffnsqrt(o, n, 5e-3, 'W1', 'W2')),
                  "W2xT": (wantCMBT, lambda o, n: diffnsqrt(o, n, 5e-3, 'W2', 'T')),
                  "W2xE": (wantCMBT, lambda o, n: diffnsqrt(o, n, 5e-3, 'W2', 'E')),
                  "W2xP": (True, lambda o, n: diffnsqrt(o, n, 5e-3, 'W2', 'P')),
                  "W2xW1": (True, lambda o, n: diffnsqrt(o, n, 5e-3, 'W2', 'W1')),
                  "W2xW2": (True, 5e-3),
                  "*" : Ignore()})

# The filetolmatrix as described above.
filetolmatrix = [["*scalCls.dat", Ignore()],  # Ignore() all scalCls.dat files.
                 ["*lensedCls.dat", coltol1],  # lensed and lenspotential files both use coltol1 given above
                 ["*lensedtotCls.dat", coltol1],
                 ["*lenspotentialCls.dat", coltol1],
                 ["*scalCovCls.dat", coltol1],
                 ["*tensCls.dat", ColTol({"TE": (True, lambda o, n: diffnsqrt(o, n, 1e-2, 'T', 'E')),
                                          "*": (True, [(0, 1e-2),
                                                (600, Ignore())])})],
                 ["*matterpower.dat", ColTol({"P": (True, lambda o, n: normabs(o["P"], n["P"], 1e-3 if n["k/h"] < 1 else 3e-3)),
                                              "*": Ignore()})],
                 ["*transfer_out.dat", ColTol({"baryon": (True, lambda o, n: normabs(o["baryon"], n["baryon"], 1e-3 if n["k/h"] < 1 else 3e-3)),
                                               "CDM": (True, lambda o, n: normabs(o["CDM"], n["CDM"], 1e-3 if n["k/h"] < 1 else 3e-3)),
                                               "v_CDM": (True, lambda o, n: normabs(o["v_CDM"], n["v_CDM"], 1e-3 if n["k/h"] < 1 else 3e-3)),
                                               "v_b": (True, lambda o, n: normabs(o["v_b"], n["v_b"], 1e-3 if n["k/h"] < 1 else 3e-3)),
                                               "*": Ignore()})],
                 ["*sharp_cl_*.dat", ColTol({"CL": (True, 1e-3),
                                             "P": (True, 1e-3),
                                             "P_vv": (True, 1e-3),
                                             "*": Ignore()})],
                 ["*", ColTol({"*": (True, args.diff_tolerance)})],
                ]


prog = os.path.abspath(args.prog)
if not os.path.exists(args.ini_dir):
    os.mkdir(args.ini_dir)

out_files_dir = os.path.join(args.ini_dir, args.out_files_dir)

if args.clean:
    if os.path.exists(out_files_dir): shutil.rmtree(out_files_dir)

if not os.path.exists(out_files_dir):
    os.mkdir(out_files_dir)

def runScript(fname):
    now = time.time()
    try:
        res = str(subprocess.check_output([prog, fname]))
        code = 0
    except subprocess.CalledProcessError as e:
        res = e.output
        code = e.returncode
    return time.time() - now, res, code

def getInis(ini_dir):
    inis = []
    for fname in os.listdir(ini_dir):
        if fnmatch.fnmatch(fname, '*.ini'):
            inis.append(os.path.join(args.ini_dir, fname))
    return inis


def getTestParams():
    params = [['base']]

    for lmax in [1000, 2000, 2500, 3000, 4500, 6000]:
        params.append(['lmax%s' % lmax, 'l_max_scalar = %s' % lmax, 'k_eta_max_scalar  = %s' % (lmax * 2.5)])

    for lmax in [1000, 2000, 2500, 3000, 4500]:
        params.append(['nonlin_lmax%s' % lmax, 'do_nonlinear =2', 'get_transfer= T', 'l_max_scalar = %s' % lmax, 'k_eta_max_scalar  = %s' % (lmax * 2.5)])

    for lmax in [400, 600, 1000]:
        params.append(['tensor_lmax%s' % lmax, 'get_tensor_cls = T', 'l_max_tensor = %s' % lmax, 'k_eta_max_tensor  = %s' % (lmax * 2)])

    params.append(['tensoronly', 'get_scalar_cls=F', 'get_tensor_cls = T'])
    params.append(['tensor_tranfer', 'get_scalar_cls=F', 'get_tensor_cls = T', 'get_transfer= T', 'transfer_high_precision = T'])
    params.append(['tranfer_only', 'get_scalar_cls=F', 'get_transfer= T', 'transfer_high_precision = F'])
    params.append(['tranfer_highprec', 'get_scalar_cls=F', 'get_transfer= T', 'transfer_high_precision = T'])

    params.append(['all', 'get_scalar_cls=T', 'get_tensor_cls = T', 'get_transfer= T'])
    params.append(['all_nonlin1', 'get_scalar_cls=T', 'get_tensor_cls = T', 'get_transfer= T', 'do_nonlinear=1'])
    params.append(['all_nonlin2', 'get_scalar_cls=T', 'get_tensor_cls = T', 'get_transfer= T', 'do_nonlinear=2'])
    params.append(['all_nonlinhigh', 'get_scalar_cls=T', 'get_tensor_cls = T', 'get_transfer= T', 'do_nonlinear=2', 'transfer_high_precision = T'])
    params.append(['tranfer_delta10', 'get_scalar_cls=F', 'get_transfer= T', 'transfer_high_precision = T', 'transfer_k_per_logint =10'])
    params.append(['tranfer_redshifts', 'get_scalar_cls=F', 'get_transfer= T', 'transfer_num_redshifts=2']
                  + ['transfer_redshift(1)=1', 'transfer_redshift(2)=0.7', 'transfer_filename(2)=transfer_out2.dat', 'transfer_matterpower(2)=matterpower2.dat']
                  )
    params.append(['tranfer_redshifts2', 'get_scalar_cls=F', 'get_transfer= T', 'transfer_num_redshifts=2']
                  + ['transfer_redshift(1)=0.7', 'transfer_redshift(2)=0', 'transfer_filename(2)=transfer_out2.dat', 'transfer_matterpower(2)=matterpower2.dat']
                  )


    params.append(['tranfer_nonu', 'get_scalar_cls=F', 'get_transfer= T', 'transfer_power_var = 8'])

    #AM - Added HMcode and halomodel tests (halofit_version=5,6)
    params.append(['HMcode','transfer_kmax=100', 'halofit_version=5', 'do_nonlinear=1', 'get_transfer= T'])
    params.append(['halomodel','transfer_kmax=100', 'halofit_version=6', 'do_nonlinear=1', 'get_transfer= T'])
    #AM - End of edits

    params.append(['zre', 're_use_optical_depth = F', 're_redshift  = 8.5'])
    params.append(['nolens', 'lensing = F'])
    params.append(['noderived', 'derived_parameters = F'])
    params.append(['no_rad_trunc', 'do_late_rad_truncation   = F'])

    for acc in [0.95, 1.1, 1.5, 2.2]:
        params.append(['accuracy_boost%s' % acc, 'accuracy_boost = %s' % acc])

    for acc in [1, 1.5, 2]:
        params.append(['l_accuracy_boost%s' % acc, 'l_accuracy_boost = %s' % acc])

    params.append(['acc', 'l_accuracy_boost =2', 'accuracy_boost=2'])
    params.append(['accsamp', 'l_accuracy_boost =2', 'accuracy_boost=2', 'l_sample_boost = 1.5'])

    params.append(['mu_massless', 'omnuh2 =0'])

    for mnu in [0, 0.01, 0.03, 0.1]:
        omnu = mnu / 100.
        params.append(['mu_mass%s' % mnu, 'omnuh2 =%s' % omnu, 'massive_neutrinos  = 3'])
    params.append(['mu_masssplit', 'omnuh2 =0.03', 'massive_neutrinos = 1 1', 'nu_mass_fractions=0.2 0.8',
                    'nu_mass_degeneracies = 1 1', 'nu_mass_eigenstates = 2', 'massless_neutrinos = 1.046'])



    for etamax in [10000, 14000, 20000, 40000]:
        params.append(['acclens_ketamax%s' % etamax, 'do_nonlinear = 2', 'l_max_scalar  = 6000', 'k_eta_max_scalar  = %s' % etamax, 'accurate_BB = F'])

    for etamax in [10000, 14000, 20000, 40000]:
        params.append(['acclensBB_ketamax%s' % etamax, 'do_nonlinear = 2', 'l_max_scalar = 2500', 'k_eta_max_scalar  = %s' % etamax, 'accurate_BB = T'])

    pars = {
        'ombh2':[ 0.0219, 0.0226, 0.0253],
        'omch2':[ 0.1, 0.08, 0.15],
        'omk':[ 0, -0.03, 0.04, 0.001, -0.001],
        'hubble':[ 62, 67, 71, 78],
        'w':[ -1.2, -1, -0.98, -0.75],
        'helium_fraction':[ 0.21, 0.23, 0.27],
        'scalar_spectral_index(1)' :[0.94, 0.98],
        'scalar_nrun(1)' :[-0.015, 0, 0.03],
        're_optical_depth': [0.03, 0.05, 0.08, 0.11],
    }

    for par, vals in pars.items():
        for val in vals:
            params.append(['%s_%.3f' % (par, val), 'get_transfer= T', 'do_nonlinear=1', 'transfer_high_precision = T',
                           '%s = %s' % (par, val)])

    if not args.no_de and not os.environ.get('CAMB_TESTS_NO_DE'):
        for wa in [-0.3, -0.01, 0.5]:
            for w in [-1.2, -0.998, -0.7]:
                params.append(['ppf_w%s_wa%s' % (w, wa), 'w = %s' % w, 'wa =%s' % wa, 'do_nonlinear = 2', 'get_transfer= T', 'dark_energy_model=PPF'])

        params.append(['ppf_w-1.000_wa0.000', 'w = -1.0', 'wa = 0.0', 'do_nonlinear = 1', 'get_transfer= T',
                       'transfer_high_precision = T', 'dark_energy_model=PPF'])

    if not args.no_sources and not os.environ.get('CAMB_TESTS_NO_SOURCES'):
        # ##CAMB sources options and new outputs
        params.append(['delta_xe', 'evolve_delta_xe =T', 'get_transfer= T', 'do_nonlinear=2', 'transfer_high_precision = T'])

        def make_win(i, z, kind, bias, sigma, s):
            return ['redshift(%s) = %s' % (i, z), 'redshift_kind(%s) = %s' % (i, kind), 'redshift_bias(%s) = %s' % (i, bias),
                    'redshift_sigma(%s) = %s' % (i, sigma), 'redshift_dlog10Ndm(%s) = %s' % (i, s)]

        counts_def = ['DEFAULT(params_counts.ini)']
        source_counts = ['num_redshiftwindows = 2'] + make_win(1, 0.3, 'counts', 1.5, 0.06, 0.42) + make_win(2, 1, 'counts', 2, 0.3, 0)
        bool_options = ['counts_evolve', 'DoRedshiftLensing', 'counts_redshift', 'evolve_delta_xe']
        for b1 in ['T', 'F']:
            for b2 in ['T', 'F']:
                for b3 in ['T', 'F']:
                    for b4 in ['T', 'F']:
                        bs = [b1, b2, b3, b4]
                        pars = copy.copy(source_counts)
                        for opt, b in zip(bool_options, bs):
                            pars += [opt + ' = ' + b]
                        params.append(['counts_opts_' + '_'.join(bs)] + counts_def + pars)
        params.append(['counts_1bin'] + counts_def
                      + ['num_redshiftwindows = 1'] + make_win(1, 0.15, 'counts', 1.2, 0.04, -0.2))
        params.append(['counts_lmax1', 'l_max_scalar = 400', 'want_CMB = F'] + counts_def + source_counts)
        params.append(['counts_lmax2', 'l_max_scalar = 1200'] + counts_def + source_counts)

        params.append(['counts_overlap'] + counts_def
                      + ['num_redshiftwindows = 2'] + make_win(1, 0.17, 'counts', 1.2, 0.04, -0.2) + make_win(2, 0.2, 'counts', 1.2, 0.04, -0.2))
        params.append(['lensing_base', 'DEFAULT(params_lensing.ini)'])
        params.append(['21cm_base', 'DEFAULT(params_21cm.ini)'])
        params.append(['21cm_base2', 'DEFAULT(params_21cm.ini)', 'get_transfer = T'])
        params.append(['counts_lens', 'DEFAULT(params_counts.ini)']
                      + ['num_redshiftwindows = 2'] + make_win(1, 0.17, 'counts', 1.2, 0.04, -0.2) + make_win(2, 0.5, 'lensing', 0, 0.07, 0.2))

    max_tests =  args.max_tests or os.environ.get('CAMB_TESTS_MAX')
    if max_tests:
        params = params[:int(max_tests)]
    return params

def list_files(file_dir):
    return [f for f in os.listdir(file_dir) if not '.ini' in f]

def output_file_num(file_dir):
    return len(list_files(file_dir))

def makeIniFiles():
    params = getTestParams()
    inis = []
    base_ini = 'inheritbase_' + os.path.basename(args.base_settings)
    shutil.copy(args.base_settings, os.path.join(args.ini_dir, base_ini))
    for pars in params:
        name = 'params_' + pars[0]
        fname = os.path.join(args.ini_dir, name + '.ini')
        inis.append(fname)
        with open(fname, 'w') as f:
            f.write('output_root=' + os.path.join(out_files_dir, name) + '\n'
                    + '\n'.join(pars[1:]) + '\nDEFAULT(' + base_ini + ')\n')
    return inis

def get_tolerance_vector(filename, cols):
    """
    Get the tolerances for the given filename.
    :param filename: The name of the file to retrieve the tolerances for.
    :param cols: Gives the column names.
    :returns: False, when the file is to be Ignore()d completely;
    the vector of tolerances when a pattern in the filetolmatrix matched;
    an empty ColTol when no match was found.
    """
    for key, val in filetolmatrix:
        if fnmatch.fnmatch(filename, key):
            if isinstance(val, Ignore):
                return False
            else:
                return [val[t] for t in cols]
    return ColTol()

def num_unequal(filename, cmpFcn):
    """
    Check whether two files are numerically unequal for the given compare function.
    :param filename: The base name of the files to check.
    :param cmpFcn: The default comparison function. Can be overriden by the filetolmatrix.
    :return: True, when the files do not match, false else.
    """
    orig_name = os.path.join(args.ini_dir, args.diff_to, filename)
    with open(orig_name) as f:
        origMat = [[x for x in ln.split()] for ln in f]
        # Check if the first row has one more column, which is the #
        if len(origMat[0]) == len(origMat[1]) + 1:
            origBase = 1
            origMat[0] = origMat[0][1:]
        else:
            origBase = 0
    new_name = os.path.join(args.ini_dir, args.out_files_dir, filename)
    with open(new_name) as f:
        newMat = [[x for x in ln.split()] for ln in f]
        if len(newMat[0]) == len(newMat[1]) + 1:
            newBase = 1
            newMat[0] = newMat[0][1:]
        else:
            newBase = 0
    if len(origMat) - origBase != len(newMat) - newBase:
        if args.verbose_diff_output:
            printlog('num rows do not match in %s: %d != %d' % (filename, len(origMat), len(newMat)))
        return True
    if newBase == 1:
        cols = [s[0] + 'x' + s[1] if len(s) == 2 and s != 'nu' else s for s in newMat[0]]
    else:
        cols = range(len(newMat[0]))

    tolerances = get_tolerance_vector(filename, cols)
    row = 0
    col = 0
    try:
        if tolerances:
            inifilenameparts = filename.rsplit('_', 2)
            inifilename = '_'.join(inifilenameparts[0:2]) if inifilenameparts[1] != 'transfer' else inifilenameparts[0]
            inifilename += "_params.ini"
            inifilename = os.path.join(args.ini_dir, args.out_files_dir, inifilename)
            try:
                # The following split fails for *_transfer_out.* files where it not needed anyway.
                inifile = iniFile()
                inifile.readFile(inifilename)
            except:
                printlog("Could not open inifilename: %s" % inifilename)
            for o_row, n_row in zip(origMat[origBase:], newMat[newBase:]):
                row += 1
                if len(o_row) != len(n_row):
                    if args.verbose_diff_output:
                        printlog('num columns do not match in %s: %d != %d' % (filename, len(o_row), len(n_row)))
                    return True
                col = 0
                of_row = [float(f) for f in o_row]
                nf_row = []
                for f in n_row:
                    try:
                        nf_row += [float(f)]
                    except ValueError:
                        sp = customsplit(f)
                        nf_row += [ float(sp[0] + 'E' + sp[1]) ]
                oldrowdict = False
                newrowdict = False
                for o, n in zip(of_row, nf_row):
                    if isinstance(tolerances[col], Ignore):
                        pass
                    else:
                        cond, tols = tolerances[col]
                        # When the column condition is bool (True or False) or a function
                        # returning False, then skip this column.
                        if isinstance(cond, bool) or not cond(inifile):
                            pass
                        else:
                            if isinstance(tols, float):
                                if not cmpFcn(o, n, tols):
                                    if args.verbose_diff_output:
                                        printlog('value mismatch at %d, %d ("%s") of %s: %s != %s' % (row, col + 1, cols[col], filename, o, n))
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
                                        if not cmpFcn(o, n, cand):
                                            if args.verbose_diff_output:
                                                printlog('value mismatch at %d, %d ("%s") of %s: %s != %s' % (row, col + 1, cols[col], filename, o, n))
                                            return True
                                    elif not isinstance(cand, Ignore):
                                        if not cand(oldrowdict, newrowdict):
                                            if args.verbose_diff_output:
                                                printlog('value mismatch at %d, %d ("%s") of %s: %s != %s' % (row, col + 1, cols[col], filename, o, n))
                                            return True
                                else:
                                    if not tols(oldrowdict, newrowdict):
                                        if args.verbose_diff_output:
                                            printlog('value mismatch at %d, %d ("%s") of %s: %s != %s' % (row, col + 1, cols[col], filename, o, n))
                                        return True
                    col += 1
            return False
        else:
#            if args.verbose_diff_output:
#                printlog("Skipped file %s" % (filename))
            return False
    except ValueError as e:
        printlog("ValueError: '%s' at %d, %d in file: %s" % (e.message, row, col + 1, filename))
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
        if s[i] == '+' or s[i] == '-':
            return [ s[0: i - 1], s[i: n] ]
        i -= 1
    return [ s ]

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
    os = customsplit(o)
    ns = customsplit(n)
    if len(os) > 1 and len(ns) > 1:
        o_mantise = float(os[0])
        o_exp = int(os[1])
        n_mantise = float(ns[0])
        n_exp = int(ns[1])
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
        return math.fabs(float(os[0] + 'E' + os[1]) - float(ns[0] + 'E' + ns[1])) >= tolerance
    # In all other cases do a numerical check
    return math.fabs(float(o) - float(n)) >= tolerance

if args.diff_to:
    import filecmp
    import math
    if args.num_diff:
        defCmpFcn = lambda o, n, t: math.fabs(float(o) - float(n)) >= t
    else:
        defCmpFcn = normabs
    out_files_dir2 = os.path.join(args.ini_dir, args.diff_to)
    match, mismatch, errors = filecmp.cmpfiles(out_files_dir, out_files_dir2,
         list(set(list_files(out_files_dir)) | set(list_files(out_files_dir2))))
    len_errors = len(errors)
    if len_errors and len_errors != 1 and errors[0] != args.diff_to:
        printlog('Missing/Extra files:')
        for err in errors:
            # Only print files that are not the diff_to
            if err != args.diff_to:
                printlog('  ' + err)
    if len(mismatch):
        numerical_mismatch = [f for f in mismatch if num_unequal(f, defCmpFcn)]
        if len(numerical_mismatch):
            printlog('Files do not match:')
            for err in numerical_mismatch:
                printlog('  ' + err)
        len_num_mismatch = len(numerical_mismatch)
    else:
        len_num_mismatch = 0

    printlog("Done with %d numerical accuracy mismatches and %d extra/missing files" % (len_num_mismatch, len_errors))
    if len_errors > 0 or len_num_mismatch > 0:
       sys.exit(1)
    else:
       sys.exit()

if args.make_ini:
    inis = makeIniFiles()
else:
    inis = getInis(args.ini_dir)
if not args.no_run_test:
    errors = 0
    files = output_file_num(out_files_dir)
    if files:
        printlog('Output directory is not empty (run with --clean to force delete): %s' % out_files_dir)
        sys.exit()
    start = time.time()
    error_list = []
    for ini in inis:
        printlog(os.path.basename(ini) + '...')
        timing, output, return_code = runScript(ini)
        if return_code:
            printlog('error %s' % return_code)
        nfiles = output_file_num(out_files_dir)
        if nfiles > files:
            msg = '..OK, produced %s files in %.2fs' % (nfiles - files, timing)
        else:
            errors += 1
            error_list.append(os.path.basename(ini))
            msg = '..no files in %.2fs' % timing
        printlog(msg)
        files = nfiles
    printlog('Done, %s errors in %.2fs (outputs not checked)' % (errors, time.time() - start))
    if errors:
        printlog('Fails in : %s' % error_list)
