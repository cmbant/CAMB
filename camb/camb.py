from .baseconfig import camblib, CAMBError, CAMBValueError, CAMBUnknownArgumentError, np
from ctypes import c_double, c_bool, POINTER, byref
import ctypes
from . import model, constants
from ._config import config
from .model import CAMBparams
from .results import CAMBdata, MatterTransferData, ClTransferData
import logging
import os
import numbers
from inspect import getfullargspec

_debug_params = False


def set_feedback_level(level=1):
    """
    Set the feedback level for internal CAMB calls

    :param level:  zero for nothing, >1 for more
    """
    config.FeedbackLevel = level


def get_results(params):
    """
    Calculate results for specified parameters and return :class:`~.results.CAMBdata` instance for getting results.

    :param params: :class:`.model.CAMBparams` instance
    :return: :class:`~.results.CAMBdata` instance
    """
    res = CAMBdata()
    if _debug_params:
        print(params)
    res.calc_power_spectra(params)
    return res


def get_transfer_functions(params, only_time_sources=False):
    """
    Calculate transfer functions for specified parameters and return :class:`~.results.CAMBdata` instance for
    getting results and subsequently calculating power spectra.

    :param params: :class:`.model.CAMBparams` instance
    :param only_time_sources: does not calculate the CMB l,k transfer functions and does not apply any non-linear
                              correction scaling. Results with only_time_sources=True can therefore be used with
                              different initial power spectra to get consistent non-linear lensed spectra.
    :return: :class:`~.results.CAMBdata` instance
    """

    res = CAMBdata()
    res.calc_transfers(params, only_transfers=True, only_time_sources=only_time_sources)
    return res


def get_background(params, no_thermo=False):
    """
    Calculate background cosmology for specified parameters and return :class:`~.results.CAMBdata`, ready to get derived
    parameters and use background functions like :func:`~results.CAMBdata.angular_diameter_distance`.

    :param params: :class:`.model.CAMBparams` instance
    :param no_thermo: set True if thermal and ionization history not required.
    :return: :class:`~.results.CAMBdata` instance
    """

    res = CAMBdata()
    if no_thermo:
        res.calc_background_no_thermo(params)
    else:
        res.calc_background(params)
    return res


def get_age(params):
    """
    Get age of universe for given set of parameters

    :param params:  :class:`.model.CAMBparams` instance
    :return: age of universe in gigayears
    """
    return CAMB_GetAge(byref(params))


def get_zre_from_tau(params, tau):
    """
    Get reionization redshift given optical depth tau

    :param params: :class:`.model.CAMBparams` instance
    :param tau: optical depth
    :return: reionization redshift (or negative number if error)
    """
    return params.Reion.get_zre(params, tau)


def set_params(cp=None, verbose=False, **params):
    """

    Set all CAMB parameters at once, including parameters which are part of the
    CAMBparams structure, as well as global parameters.

    E.g.::

      cp = camb.set_params(ns=1, H0=67, ombh2=0.022, omch2=0.1, w=-0.95, Alens=1.2, lmax=2000,
                           WantTransfer=True, dark_energy_model='DarkEnergyPPF')

    This is equivalent to::

      cp = model.CAMBparams()
      cp.DarkEnergy = DarkEnergyPPF()
      cp.DarkEnergy.set_params(w=-0.95)
      cp.set_cosmology(H0=67, omch2=0.1, ombh2=0.022, Alens=1.2)
      cp.set_for_lmax(lmax=2000)
      cp.InitPower.set_params(ns=1)
      cp.WantTransfer = True

    The wrapped functions are (in this order):

    * :meth:`.model.CAMBparams.set_accuracy`
    * :meth:`.model.CAMBparams.set_classes`
    * :meth:`.dark_energy.DarkEnergyEqnOfState.set_params` (or equivalent if a different dark energy model class used)
    * :meth:`.model.CAMBparams.set_cosmology`
    * :meth:`.model.CAMBparams.set_matter_power`
    * :meth:`.model.CAMBparams.set_for_lmax`
    * :meth:`.initialpower.InitialPowerLaw.set_params`  (or equivalent if a different initial power model class used)
    * :meth:`.nonlinear.Halofit.set_params`

    :param params: the values of the parameters
    :param cp: use this CAMBparams instead of creating a new one
    :param verbose: print out the equivalent set of commands
    :return: :class:`.model.CAMBparams` instance

    """

    if 'ALens' in params:
        raise ValueError('Use Alens not ALens')

    if cp is None:
        cp = model.CAMBparams()
    else:
        assert isinstance(cp, model.CAMBparams), "cp should be an instance of CAMBparams"

    used_params = set()

    def do_set(setter):
        kwargs = {kk: params[kk] for kk in getfullargspec(setter).args[1:] if kk in params}
        used_params.update(kwargs)
        if kwargs:
            if verbose:
                logging.warning('Calling %s(**%s)' % (setter.__name__, kwargs))
            setter(**kwargs)

    # Note order is important: must call DarkEnergy.set_params before set_cosmology if setting theta rather than H0
    # set_classes allows redefinition of the classes used, so must be called before setting class parameters
    do_set(cp.set_accuracy)
    do_set(cp.set_classes)
    do_set(cp.DarkEnergy.set_params)
    do_set(cp.set_cosmology)
    do_set(cp.set_matter_power)
    do_set(cp.set_for_lmax)
    do_set(cp.InitPower.set_params)
    do_set(cp.NonLinearModel.set_params)

    if cp.InitPower.has_tensors():
        cp.WantTensors = True

    unused_params = set(params) - used_params
    if unused_params:
        for k in unused_params:
            obj = cp
            if '.' in k:
                parts = k.split('.')
                for p in parts[:-1]:
                    obj = getattr(obj, p)
                par = parts[-1]
            else:
                par = k
            if hasattr(obj, par):
                setattr(obj, par, params[k])
            else:
                raise CAMBUnknownArgumentError("Unrecognized parameter: %s" % k)
    return cp


def get_valid_numerical_params(transfer_only=False, **class_names):
    """
    Get numerical parameter names that are valid input to :func:`set_params`

    :param transfer_only: if True, exclude parameters that affect only initial power spectrum or non-linear model
    :param class_names: class name parameters that will be used by :meth:`.model.CAMBparams.set_classes`
    :return: set of valid input parameter names for :func:`set_params`
    """
    cp = CAMBparams()
    cp.set_classes(**class_names)
    params = set()

    def extract_params(set_func):
        pars = getfullargspec(set_func)
        for arg in pars.args[1:len(pars.args) - len(pars.defaults or [])]:
            params.add(arg)
        if pars.defaults:
            for arg, v in zip(pars.args[len(pars.args) - len(pars.defaults):], pars.defaults):
                if (isinstance(v, numbers.Number) or v is None) and 'version' not in arg:
                    params.add(arg)

    extract_params(cp.DarkEnergy.set_params)
    extract_params(cp.set_cosmology)
    if not transfer_only:
        extract_params(cp.InitPower.set_params)
        extract_params(cp.NonLinearModel.set_params)
    for f, tp in cp._fields_:
        if not f.startswith('_') and tp == ctypes.c_double:
            params.add(f)
    return params - {'max_eta_k_tensor', 'max_eta_k', 'neutrino_hierarchy', 'standard_neutrino_neff',
                     'pivot_scalar', 'pivot_tensor', 'num_massive_neutrinos', 'num_nu_massless', 'bbn_predictor'}


def set_params_cosmomc(p, num_massive_neutrinos=1, neutrino_hierarchy='degenerate', halofit_version='mead',
                       dark_energy_model='ppf', lmax=2500, lens_potential_accuracy=1, inpars=None):
    """
    get CAMBParams for dictionary of cosmomc-named parameters assuming Planck 2018 defaults

    :param p: dictionary of cosmomc parameters (e.g. from getdist.types.BestFit's getParamDict() function)
    :param num_massive_neutrinos: usually 1 if fixed mnu=0.06 eV, three if mnu varying
    :param neutrino_hierarchy: hierarchy
    :param halofit_version: name of the soecific Halofit model to use for non-linear modelling
    :param dark_energy_model: ppf or fluid dark energy model
    :param lmax: lmax for accuracy settings
    :param lens_potential_accuracy: lensing accuracy parameter
    :param inpars: optional input CAMBParams to set
    :return:
    """
    pars = inpars or model.CAMBparams()
    if p.get('alpha1', 0) or p.get('Aphiphi', 1) != 1:
        raise ValueError('Parameter not currrently supported by set_params_cosmomc')
    pars.set_cosmology(H0=p['H0'], ombh2=p['omegabh2'], omch2=p['omegach2'], mnu=p.get('mnu', 0.06),
                       omk=p.get('omegak', 0), tau=p['tau'], deltazrei=p.get('deltazrei', None),
                       nnu=p.get('nnu', constants.default_nnu), Alens=p.get('Alens', 1.0),
                       YHe=p.get('yheused', None), meffsterile=p.get('meffsterile', 0),
                       num_massive_neutrinos=num_massive_neutrinos, neutrino_hierarchy=neutrino_hierarchy)
    pars.InitPower.set_params(ns=p['ns'], r=p.get('r', 0), As=p['A'] * 1e-9, nrun=p.get('nrun', 0),
                              nrunrun=p.get('nrunrun', 0))
    pars.set_dark_energy(w=p.get('w', -1), wa=p.get('wa', 0), dark_energy_model=dark_energy_model)
    pars.set_for_lmax(lmax, lens_potential_accuracy=lens_potential_accuracy)
    pars.NonLinearModel.set_params(halofit_version=halofit_version)
    pars.WantTensors = pars.InitPower.has_tensors()
    return pars


def validate_ini_file(filename):
    # Check if fortran .ini file parameters are valid; catch error stop in separate process
    import subprocess
    import sys
    try:
        err = ''
        command = '"%s" "%s" "%s" --validate' % (
            sys.executable, os.path.join(os.path.dirname(__file__), '_command_line.py'), filename)
        subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as E:
        err = E.output.decode().replace('ERROR STOP', '').strip()
    if err:
        raise CAMBValueError(err + ' (%s)' % filename)
    return True


def run_ini(ini_filename, no_validate=False):
    """
    Run the command line camb from a .ini file (producing text files as with the command line program).
    This does the same as the command line program, except global config parameters are not read and set (which does not
    change results in almost all cases).

    :param ini_filename: .ini file to use
    :param no_validate: do not pre-validate the ini file (faster, but may crash kernel if error)
    """
    if not os.path.exists(ini_filename):
        raise CAMBValueError('File not found: %s' % ini_filename)
    if not no_validate:
        validate_ini_file(ini_filename)
    run_inifile = camblib.__camb_MOD_camb_runinifile
    run_inifile.argtypes = [ctypes.c_char_p, POINTER(ctypes.c_long)]
    run_inifile.restype = c_bool
    s = ctypes.create_string_buffer(ini_filename.encode("latin-1"))
    if not run_inifile(s, ctypes.c_long(len(ini_filename))):
        config.check_global_error('run_ini')


def read_ini(ini_filename, no_validate=False):
    """
    Get a :class:`.model.CAMBparams` instance using parameter specified in a .ini parameter file.

    :param ini_filename: path of the .ini file to read
    :param no_validate: do not pre-validate the ini file (faster, but may crash kernel if error)
    :return: :class:`.model.CAMBparams` instance
    """
    if not os.path.exists(ini_filename):
        raise CAMBValueError('File not found: %s' % ini_filename)
    if not no_validate:
        validate_ini_file(ini_filename)
    cp = model.CAMBparams()
    read_inifile = camblib.__camb_MOD_camb_readparamfile
    read_inifile.argtypes = [POINTER(CAMBparams), ctypes.c_char_p, POINTER(ctypes.c_long)]
    read_inifile.restype = ctypes.c_bool
    s = ctypes.create_string_buffer(ini_filename.encode("latin-1"))
    if not read_inifile(cp, s, ctypes.c_long(len(ini_filename))):
        config.check_global_error('read_ini')
    return cp


def get_matter_power_interpolator(params, zmin=0, zmax=10, nz_step=100, zs=None, kmax=10, nonlinear=True,
                                  var1=None, var2=None, hubble_units=True, k_hunit=True,
                                  return_z_k=False, k_per_logint=None, log_interp=True, extrap_kmax=None):
    r"""
    Return a 2D spline interpolation object to evaluate matter power spectrum as function of z and k/h, e.g.

    .. code-block:: python

       from camb import get_matter_power_interpolator
       PK = get_matter_power_interpolator(params);
       print('Power spectrum at z=0.5, k/h=0.1/Mpc is %s (Mpc/h)^3 '%(PK.P(0.5, 0.1)))

    For a description of outputs for different var1, var2 see :ref:`transfer-variables`.
    If you already have a :class:`~.results.CAMBdata` result object, you can instead
    use :meth:`~.results.CAMBdata.get_matter_power_interpolator`.

    :param params: :class:`.model.CAMBparams` instance
    :param zmin: minimum z (use 0 or smaller than you want for good interpolation)
    :param zmax: maximum z (use larger than you want for good interpolation)
    :param nz_step: number of steps to sample in z (default max allowed is 100)
    :param zs: instead of zmin,zmax, nz_step, can specific explicit array of z values to spline from
    :param kmax: maximum k
    :param nonlinear: include non-linear correction from halo model
    :param var1: variable i (index, or name of variable; default delta_tot)
    :param var2: variable j (index, or name of variable; default delta_tot)
    :param hubble_units: if true, output power spectrum in :math:`({\rm Mpc}/h)^{3}` units,
                         otherwise :math:`{\rm Mpc}^{3}`
    :param k_hunit: if true, matter power is a function of k/h, if false, just k (both :math:`{\rm Mpc}^{-1}` units)
    :param return_z_k: if true, return interpolator, z, k where z, k are the grid used
    :param k_per_logint: specific uniform sampling over log k (if not set, uses optimized irregular sampling)
    :param log_interp: if true, interpolate log of power spectrum (unless any values are negative in which case ignored)
    :param extrap_kmax: if set, use power law extrapolation beyond kmax to extrap_kmax (useful for tails of integrals)
    :return: An object PK based on :class:`~scipy:scipy.interpolate.RectBivariateSpline`, that can be called
             with PK.P(z,kh) or PK(z,log(kh)) to get log matter power values.
             If return_z_k=True, instead return interpolator, z, k where z, k are the grid used.
    """

    pars = params.copy()

    if zs is None:
        zs = zmin + np.exp(np.log(zmax - zmin + 1) * np.linspace(0, 1, nz_step)) - 1
    pars.set_matter_power(redshifts=zs, kmax=kmax, k_per_logint=k_per_logint, silent=True)
    pars.NonLinear = model.NonLinear_none
    results = get_results(pars)

    return results.get_matter_power_interpolator(nonlinear=nonlinear, var1=var1, var2=var2, hubble_units=hubble_units,
                                                 k_hunit=k_hunit, return_z_k=return_z_k, log_interp=log_interp,
                                                 extrap_kmax=extrap_kmax)


CAMB_GetAge = camblib.__camb_MOD_camb_getage
CAMB_GetAge.restype = c_double
CAMB_GetAge.argtypes = [POINTER(model.CAMBparams)]
