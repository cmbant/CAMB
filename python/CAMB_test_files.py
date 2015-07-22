from __future__ import print_function
import os
import subprocess
import argparse
import fnmatch
import time
import sys
import shutil
import copy

parser = argparse.ArgumentParser(description='Run CAMB tests')
parser.add_argument('ini_dir', help='ini file directory')
parser.add_argument('--make_ini', action='store_true', help='if set, output ini files to ini_dir')
parser.add_argument('--out_files_dir', default='test_outputs', help='output files directory')
parser.add_argument('--base_settings', default='params.ini', help='settings to include as defaults for all combinations')
parser.add_argument('--no_run_test', action='store_true', help='don''t run tests on files')
parser.add_argument('--prog', default='./camb', help='executable to run')
parser.add_argument('--clean', action='store_true', help='delete output dir before run')
parser.add_argument('--diff_to', help='output directory to compare to, e.g. test_outputs2')
parser.add_argument('--diff_tolerance', type=float, help='the tolerance for the numerical diff [default: 1e-4]', default=1e-4)
parser.add_argument('--verbose_diff_output', action='store_true', help='during diff_to print more error messages')
parser.add_argument('--num_diff', action='store_true', help='during diff_to print more error messages')

args = parser.parse_args()

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
    params = []

    params.append(['base'])

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

    for wa in [-0.3, -0.01, 0.5]:
        for w in [-1.2, -0.998, -0.7]:
            params.append(['ppf_w%s_wa%s' % (w, wa), 'w = %s' % w, 'wa =%s' % wa, 'do_nonlinear = 2', 'get_transfer= T', 'dark_energy_model=PPF'])

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

    for par, vals in pars.iteritems():
        for val in vals:
            params.append(['%s_%.3f' % (par, val), 'get_transfer= T', 'do_nonlinear=1', 'transfer_high_precision = T',
                           '%s = %s' % (par, val)])

    params.append(['ppf_w-1.000_wa0.000', 'w = -1.0', 'wa = 0.0', 'do_nonlinear = 1', 'get_transfer= T',
                   'transfer_high_precision = T', 'dark_energy_model=PPF'])



    ###CAMB sources options and new outputs
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
    params.append(['counts_lens', 'DEFAULT(params_counts.ini)']
                  + ['num_redshiftwindows = 2'] + make_win(1, 0.17, 'counts', 1.2, 0.04, -0.2) + make_win(2, 0.5, 'lensing', 0, 0.07, 0.2))

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

def num_unequal(filename, cmpFcn):
    # Check whether two files are unequal for the given tolerance
    orig_name = os.path.join(args.ini_dir, args.diff_to, filename)
    with open(orig_name) as f:
        origMat = [[x for x in ln.split()] for ln in f]
    new_name = os.path.join(args.ini_dir, args.out_files_dir, filename)
    with open(new_name) as f:
        newMat = [[x for x in ln.split()] for ln in f]
    if len(origMat) != len(newMat):
        if args.verbose_diff_output:
            print('num rows do not match in %s: %d != %d' % (filename, len(origMat), len(newMat)))
        return True
    row = 0
    for o_row, n_row in zip(origMat, newMat):
        row = row + 1
        col = 0
        if len(o_row) != len(n_row):
            if args.verbose_diff_output:
                print('num columns do not match in %s: %d != %d' % (filename, len(o_row), len(n_row)))
            return True
        for o, n in zip(o_row, n_row):
            col = col + 1
            if cmpFcn(o, n):
                if args.verbose_diff_output:
                    print('value mismatch at %d, %d of %s: %s != %s' % (row, col, filename, o, n))
                return True
    return False

# Do a textual comparision for numbers whose exponent is zero or greater.
# The fortran code writes floating point values, with 5 significant digits
# after the comma and an exponent. I.e., for numbers with a positive
# exponent the usual comparison against a delta fails.
def textualcmp(o, n):
    os = o.split('E', 1)
    ns = n.split('E', 1)
    if len(os) > 1 and len(ns) > 1:
        o_mantise = float(os[0])
        o_exp = int(os[1])
        n_mantise = float(ns[0])
        n_exp = int(ns[1])
        if o_exp != n_exp:
            return True
        # Check without respect of the exponent, when that is greater zero.
        if 0 <= o_exp:
            return math.fabs(float(o_mantise) - float(n_mantise)) >= args.diff_tolerance
    # In all other cases do a numerical check
    return math.fabs(float(o) - float(n)) >= args.diff_tolerance

if args.diff_to:
    import filecmp
    import math
    if args.num_diff:
        defCmpFcn = lambda o, n: math.fabs(float(o) - float(n)) >= args.diff_tolerance
    else:
        defCmpFcn = textualcmp
    out_files_dir2 = os.path.join(args.ini_dir, args.diff_to)
    match, mismatch, errors = filecmp.cmpfiles(out_files_dir, out_files_dir2,
         list(set(list_files(out_files_dir)) | set(list_files(out_files_dir2))))
    if len(errors) and len(errors) != 1 and errors[0] != args.diff_to:
        len_errors = len(errors) - 1
        print('Missing/Extra files:')
        for err in errors:
            # Only print files that are not the diff_to
            if (err != args.diff_to):
                print('  ' + err)
    else:
        len_errors = 0
    if len(mismatch):
        numerical_mismatch = [f for f in mismatch if num_unequal(f, defCmpFcn)]
        if len(numerical_mismatch):
            print('Files do not match:')
            for err in numerical_mismatch:
                print('  ' + err)
        len_num_mismatch = len(numerical_mismatch)
    else:
        len_num_mismatch = 0

    print("Done with %d (%d) mismatches and %d extra/missing files" % (len_num_mismatch, len(mismatch), len_errors))
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
        print('Output directory is not empty (run with --clean to force delete): %s' % out_files_dir)
        sys.exit()
    start = time.time()
    error_list = []
    for ini in inis:
        print(os.path.basename(ini) + '...')
        timing, output, return_code = runScript(ini)
        if return_code:
            print('error %s' % return_code)
        nfiles = output_file_num(out_files_dir)
        if nfiles > files:
            msg = '..OK, produced %s files in %.2fs' % (nfiles - files, timing)
        else:
            errors += 1
            error_list.append(os.path.basename(ini))
            msg = '..no files in %.2fs' % (timing)
        print(msg)
        files = nfiles
    print('Done, %s errors in %.2fs (outputs not checked)' % (errors, time.time() - start))
    if errors:
        print('Fails in : %s' % error_list)
