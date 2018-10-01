from .baseconfig import camblib, CAMB_Structure, CAMBError, CAMBParamRangeError, dll_import
from ctypes import c_bool, c_int, c_double, c_float, byref, POINTER
from . import reionization as ion
from . import recombination as recomb
from . import initialpower as ipow
from . import constants
import numpy as np
from . import bbn
import logging
import six
import copy

# ---Parameters

max_nu = 5
max_transfer_redshifts = 150
nthermo_derived = 13
Transfer_kh = 1
Transfer_cdm = 2
Transfer_b = 3
Transfer_g = 4
Transfer_r = 5
Transfer_nu = 6
Transfer_tot = 7
Transfer_nonu = 8
Transfer_tot_de = 9
Transfer_Weyl = 10
Transfer_Newt_vel_cdm = 11
Transfer_Newt_vel_baryon = 12
Transfer_vel_baryon_cdm = 13
Transfer_max = Transfer_vel_baryon_cdm

NonLinear_none = 0
NonLinear_pk = 1
NonLinear_lens = 2
NonLinear_both = 3

derived_names = ['age', 'zstar', 'rstar', 'thetastar', 'DAstar', 'zdrag',
                 'rdrag', 'kd', 'thetad', 'zeq', 'keq', 'thetaeq', 'thetarseq']

transfer_names = ['k/h', 'delta_cdm', 'delta_baryon', 'delta_photon', 'delta_neutrino', 'delta_nu', 'delta_tot',
                  'delta_nonu', 'delta_tot_de', 'Weyl', 'v_newtonian_cdm', 'v_newtonian_baryon', 'v_baryon_cdm']

evolve_names = transfer_names + ['a', 'etak', 'H', 'growth', 'v_photon', 'pi_photon', \
                                 'E_2', 'v_neutrino', 'T_source', 'E_source', 'lens_potential_source']

background_names = ['x_e', 'opacity', 'visibility', 'cs2b', 'T_b']

neutrino_hierarchies = ['normal', 'inverted', 'degenerate']
neutrino_hierarchy_normal = 1
neutrino_hierarchy_inverted = 2
neutrino_hierarchy_degenerate = 3

# ---Variables in modules.f90
# To set the value please just put
# variable_name.value = new_value


# logical
_HighAccuracyDefault = dll_import(POINTER(c_bool), "modelparams", "highaccuracydefault")
_HighAccuracyDefault.value = True

_lSampleBoost = dll_import(c_double, "modelparams", "lsampleboost")
_AccuracyBoost = dll_import(c_double, "modelparams", "accuracyboost")
_lAccuracyBoost = dll_import(c_float, "modelparams", "laccuracyboost")
_DoLateRadTruncation = dll_import(c_bool, "gaugeinterface", "dolateradtruncation")

DebugParam = dll_import(c_double, "modelparams", "debugparam")
# DebugParam.value = 1000000*2

# logical
do_bispectrum = dll_import(c_int, "modelparams", "do_bispectrum")
# do_bispectrum.value = False

max_bessels_l_index = dll_import(c_int, "modelparams", "max_bessels_l_index")
# max_bessels_l_index.value = 1000000

max_bessels_etak = dll_import(c_double, "modelparams", "max_bessels_etak")
# max_bessels_etak.value = 1000000*2

# logical
call_again = dll_import(c_int, "modelparams", "call_again")
# call_again.value = False


grhom = dll_import(c_double, "modelparams", "grhom")
grhog = dll_import(c_double, "modelparams", "grhog")
grhor = dll_import(c_double, "modelparams", "grhor")
grhob = dll_import(c_double, "modelparams", "grhob")
grhoc = dll_import(c_double, "modelparams", "grhoc")
grhov = dll_import(c_double, "modelparams", "grhov")
grhornomass = dll_import(c_double, "modelparams", "grhornomass")
grhok = dll_import(c_double, "modelparams", "grhok")

taurst = dll_import(c_double, "modelparams", "taurst")
dtaurec = dll_import(c_double, "modelparams", "dtaurec")
taurend = dll_import(c_double, "modelparams", "taurend")
tau_maxvis = dll_import(c_double, "modelparams", "tau_maxvis")
adotrad = dll_import(c_double, "modelparams", "adotrad")

grhormass = dll_import(c_double * max_nu, "modelparams", "grhormass")
nu_masses = dll_import(c_double * max_nu, "modelparams", "nu_masses")

akthom = dll_import(c_double, "modelparams", "akthom")
fHe = dll_import(c_double, "modelparams", "fhe")
Nnow = dll_import(c_double, "modelparams", "nnow")

limber_phiphi = dll_import(c_int, "modelparams", "limber_phiphi")
# limber_phiphi.value = 0

num_extra_redshiftwindows = dll_import(c_int, "modelparams", "num_extra_redshiftwindows")
# num_extra_redshiftwindows.value = 0

num_redshiftwindows = dll_import(c_int, "modelparams", "num_redshiftwindows")

num_custom_sources = dll_import(c_int, "modelparams", "num_custom_sources")

# logical
use_spline_template = dll_import(c_bool, "modelparams", "use_spline_template")
# use_spline_template.value = True

ThermoDerivedParams = dll_import(c_double * nthermo_derived, "modelparams", "thermoderivedparams")
# ThermoDerivedParams.value = 1.

# logical
Log_lvalues = dll_import(c_bool, "lvalues", "log_lvalues")
# Log_lvalues.value = False

# Variables from module ModelData

# logical
has_cl_2D_array = dll_import(c_bool, "modeldata", "has_cl_2d_array")
# has_cl_2D_array.value = False


lmax_lensed = dll_import(c_int, "modeldata", "lmax_lensed")

# Variable from module Transfer
# logical
transfer_interp_matterpower = dll_import(c_bool, "transfer", "transfer_interp_matterpower")
# transfer_interp_matterpower.value = False

transfer_power_var = dll_import(c_int, "transfer", "transfer_power_var")
# transfer_power_var.value = Transfer_tot

# logical
get_growth_sigma8 = dll_import(c_bool, "transfer", "get_growth_sigma8")
# get_growth_sigma8.value = True

CAMB_validateparams = camblib.__camb_MOD_camb_validateparams
CAMB_validateparams.restype = c_bool

# args for these set below after CAMBparams defined
CAMB_setinitialpower = camblib.__handles_MOD_camb_setinitialpower
CAMB_SetNeutrinoHierarchy = camblib.__camb_MOD_camb_setneutrinohierarchy

numpy_1d = np.ctypeslib.ndpointer(c_double, flags='C_CONTIGUOUS')
CAMB_primordialpower = camblib.__handles_MOD_camb_primordialpower
CAMB_primordialpower.restype = c_bool


class TransferParams(CAMB_Structure):
    """
    Object storing parameters for the matter power spectrum calculation. PK variables are for setting main outputs.
    Other entries are used internally, e.g. for sampling to get correct non-linear corrections and lensing.

    :ivar high_precision: True for more accuracy
    :ivar kmax: k_max to output
    :ivar k_per_logint: number of points per log k interval. If zero, set an irregular optimized spacing.
    :ivar PK_num_redshifts: number of redshifts to calculate
    :ivar PK_redshifts: redshifts to output for the matter transfer and power

    """
    _fields_ = [
        ("high_precision", c_int),  # logical
        ("accurate_massive_neutrinos", c_int),  # logical
        ("num_redshifts", c_int),
        ("kmax", c_double),
        ("k_per_logint", c_int),
        ("redshifts", c_double * max_transfer_redshifts),
        ("PK_redshifts", c_double * max_transfer_redshifts),
        ("NLL_redshifts", c_double * max_transfer_redshifts),
        ("PK_redshifts_index", c_int * max_transfer_redshifts),
        ("NLL_redshifts_index", c_int * max_transfer_redshifts),
        ("PK_num_redshifts", c_int),
        ("NLL_num_redshifts", c_int)
    ]


class CAMBparams(CAMB_Structure):
    """
    Object storing the parameters for a CAMB calculation, including cosmological parameters and
    settings for what to calculate. When a new object is instantiated, default parameters are set automatically.

    You can view the set of underlying parameters used by the Fortran code by printing the CAMBparams instance.

    To add a new parameter, add it to the CAMBparams type in modules.f90, then  edit the _fields_ list in the CAMBparams
    class in model.py to add the new parameter in the corresponding location of the member list. After rebuilding the
    python version you can then access the parameter by using params.new_parameter where params is a CAMBparams instance.
    You could also modify the wrapper functions to set the field value less directly.
    """

    def __init__(self):
        getattr(camblib, '__camb_MOD_camb_setdefparams')(byref(self))
        if _HighAccuracyDefault.value:
            self.AccurateReionization = True
        self.InitPower.set_params()

    _fields_ = [
        ("WantCls", c_int),  # logical
        ("WantTransfer", c_int),  # logical
        ("WantScalars", c_int),  # logical
        ("WantTensors", c_int),  # logical
        ("WantVectors", c_int),  # logical
        ("DoLensing", c_int),  # logical
        ("want_zstar", c_int),  # logical
        ("want_zdrag", c_int),  # logical
        ("PK_WantTransfer", c_int),  # logical
        ("NonLinear", c_int),
        ("Want_CMB", c_int),  # logical
        ("max_l", c_int),
        ("max_l_tensor", c_int),
        ("max_eta_k", c_double),
        ("max_eta_k_tensor", c_double),
        ("omegab", c_double),
        ("omegac", c_double),
        ("omegav", c_double),
        ("omegan", c_double),
        ("H0", c_double),
        ("TCMB", c_double),
        ("YHe", c_double),
        ("num_nu_massless", c_double),
        ("num_nu_massive", c_int),
        ("nu_mass_eigenstates", c_int),
        ("share_delta_neff", c_int),  # logical
        ("nu_mass_degeneracies", c_double * max_nu),
        ("nu_mass_fractions", c_double * max_nu),
        ("nu_mass_numbers", c_int * max_nu),
        ("scalar_initial_condition", c_int),
        ("OutputNormalization", c_int),
        ("AccuratePolarization", c_int),  # logical
        ("AccurateBB", c_int),  # logical
        ("AccurateReionization", c_int),  # logical
        ("MassiveNuMethod", c_int),
        ("InitPower", ipow.InitialPowerParams),
        ("Reion", ion.ReionizationParams),
        ("Recomb", recomb.RecombinationParams),
        ("Transfer", TransferParams),
        ("InitialConditionVector", c_double * 10),
        ("OnlyTransfers", c_int),  # logical
        ("DerivedParameters", c_int),  # logical
        ("ReionHist", ion.ReionizationHistory),
        ("flat", c_int),  # logical
        ("closed", c_int),  # logical
        ("open", c_int),  # logical
        ("omegak", c_double),
        ("curv", c_double),
        ("r", c_double),
        ("Ksign", c_double),
        ("tau0", c_double),
        ("chi0", c_double)
    ]

    def copy(self):
        """
        Make independent copy.

         :return: copy of self
        """
        return copy.deepcopy(self)

    def validate(self):
        """
        Do some quick tests for sanity

        :return: True if OK
        """
        return CAMB_validateparams(byref(self))

    def set_accuracy(self, AccuracyBoost=1., lSampleBoost=1., lAccuracyBoost=1.,
                     HighAccuracyDefault=True, DoLateRadTruncation=True):
        """
        Set parameters determining calculation accuracy (large values may give big slow down).
        Note curently these are set globally, not just per parameter set.

        :param AccuracyBoost: increase AccuracyBoost to decrease integration step size, increase density of k sampling, etc.
        :param lSampleBoost: increase lSampleBoost to increase density of L sampling for CMB
        :param lAccuracyBoost: increase lAccuracyBoost to increase the maximum L included in the Boltzmann hierarchies
        :param HighAccuracyDefault: True for Planck-level accuracy (False is WMAP)
        :param DoLateRadTruncation: If True, use approximation to radiation perturbation evolution at late times
        :return: self
        """
        _lSampleBoost.value = lSampleBoost
        _AccuracyBoost.value = AccuracyBoost
        _lAccuracyBoost.value = lAccuracyBoost
        _HighAccuracyDefault.value = HighAccuracyDefault
        _DoLateRadTruncation.value = DoLateRadTruncation
        logging.warning('accuracy parameters are changed globally, not yet per parameter set')
        return self

    def set_initial_power(self, initial_power_params):
        """
        Set the InitialPower primordial power spectrum parameters

        :param initial_power_params: :class:`.initialpower.InitialPowerParams` instance
        :return: self
        """
        assert (isinstance(initial_power_params, ipow.InitialPowerParams))
        CAMB_setinitialpower(byref(self), byref(initial_power_params))
        return self

    def set_cosmology(self, H0=67.0, cosmomc_theta=None, ombh2=0.022, omch2=0.12, omk=0.0,
                      neutrino_hierarchy='degenerate', num_massive_neutrinos=1,
                      mnu=0.06, nnu=3.046, YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046,
                      TCMB=constants.COBE_CMBTemp, tau=None, deltazrei=None, bbn_predictor=None,
                      theta_H0_range=[10, 100]):
        """
        Sets cosmological parameters in terms of physical densities and parameters used in Planck 2015 analysis.
        Default settings give a single distinct neutrino mass eigenstate, by default one neutrino with mnu = 0.06eV.
        Set the neutrino_hierarchy parameter to normal or inverted to use a two-eigenstate model that is a good
        approximation to the known mass splittings seen in oscillation measurements.
        If you require more fine-grained control you can set the neutrino parameters directly rather than using this function.

        :param H0: Hubble parameter (in km/s/Mpc)
        :param cosmomc_theta: The CosmoMC theta parameter. You must set H0=None to solve for H0 given cosmomc_theta. Note that
                you must have already set the dark energy model, you can't use set_cosmology with cosmomc_theta and then
                change the background evolution (which would change cosmomc_theta at the calculated H0 value).
        :param ombh2: physical density in baryons
        :param omch2:  physical density in cold dark matter
        :param omk: Omega_K curvature parameter
        :param neutrino_hierarchy: 'degenerate', 'normal', or 'inverted' (1 or 2 eigenstate approximation)
        :param num_massive_neutrinos:  number of massive neutrinos (ignored unless hierarchy == 'degenerate')
        :param mnu: sum of neutrino masses (in eV)
        :param nnu: N_eff, effective relativistic degrees of freedom
        :param YHe: Helium mass fraction. If None, set from BBN consistency.
        :param meffsterile: effective mass of sterile neutrinos
        :param standard_neutrino_neff:  default value for N_eff in standard cosmology (non-integer to allow for partial
                heating of neutrinos at electron-positron annihilation and QED effects)
        :param TCMB: CMB temperature (in Kelvin)
        :param tau: optical depth; if None, current Reion settings are not changed
        :param deltazrei: redshift width of reionization; if None, uses default
        :param bbn_predictor: :class:`.bbn.BBNPredictor` instance used to get YHe from BBN consistency if YHe is None
        :param theta_H0_range: if cosmomc_theta is specified, the min, max interval of H0 values to map to; outside this range
                 it will raise an exception.

        """

        if YHe is None:
            # use BBN prediction
            self.bbn_predictor = bbn_predictor or bbn.get_default_predictor()
            YHe = self.bbn_predictor.Y_He(ombh2, nnu - standard_neutrino_neff)
        self.YHe = YHe

        if cosmomc_theta is not None:
            if not (0.001 < cosmomc_theta < 0.1):
                raise CAMBParamRangeError('cosmomc_theta looks wrong (parameter is just theta, not 100*theta)')

            kw = locals()
            [kw.pop(x) for x in ['self', 'H0', 'cosmomc_theta']]

            if H0 is not None:
                raise CAMBError('Set H0=None when setting cosmomc_theta.')

            try:
                from scipy.optimize import brentq
            except ImportError:
                raise CAMBError('You need SciPy to set cosmomc_theta.')

            from . import camb

            def f(H0):
                self.set_cosmology(H0=H0, **kw)
                return camb.get_background(self, no_thermo=True).cosmomc_theta() - cosmomc_theta

            try:
                self.H0 = brentq(f, theta_H0_range[0], theta_H0_range[1], rtol=1e-4)
            except ValueError:
                raise CAMBParamRangeError('No solution for H0 inside of theta_H0_range')
        else:
            self.H0 = H0

        self.TCMB = TCMB
        fac = (self.H0 / 100.0) ** 2
        self.omegab = ombh2 / fac
        self.omegac = omch2 / fac

        neutrino_mass_fac = 94.07
        # conversion factor for thermal with Neff=3 TCMB=2.7255

        if isinstance(neutrino_hierarchy, six.string_types):
            if not neutrino_hierarchy in neutrino_hierarchies:
                raise CAMBError('Unknown neutrino_hierarchy {0:s}'.format(neutrino_hierarchy))
            neutrino_hierarchy = neutrino_hierarchies.index(neutrino_hierarchy) + 1

        if (nnu >= standard_neutrino_neff or neutrino_hierarchy != neutrino_hierarchy_degenerate):
            omnuh2 = mnu / neutrino_mass_fac * (standard_neutrino_neff / 3) ** 0.75
        else:
            omnuh2 = mnu / neutrino_mass_fac * (nnu / 3.0) ** 0.75
        omnuh2_sterile = meffsterile / neutrino_mass_fac
        if omnuh2_sterile > 0 and nnu < standard_neutrino_neff:
            raise CAMBError('sterile neutrino mass required Neff>3.046')
        if omnuh2 and not num_massive_neutrinos:
            raise CAMBError('non-zero mnu with zero num_massive_neutrinos')

        omnuh2 = omnuh2 + omnuh2_sterile
        self.omegan = omnuh2 / fac
        self.omegam = self.omegab + self.omegac + self.omegan
        self.omegav = 1 - omk - self.omegam
        # self.share_delta_neff = False
        # self.nu_mass_eigenstates = 0
        # self.num_nu_massless = nnu
        # self.nu_mass_numbers[0] = 0
        # self.num_nu_massive = 0
        # if omnuh2 > omnuh2_sterile:
        #     neff_massive_standard = num_massive_neutrinos * standard_neutrino_neff / 3.0
        #     self.num_nu_massive = num_massive_neutrinos
        #     self.nu_mass_eigenstates = self.nu_mass_eigenstates + 1
        #     if nnu > neff_massive_standard:
        #         self.num_nu_massless = nnu - neff_massive_standard
        #     else:
        #         self.num_nu_massless = 0
        #         neff_massive_standard = nnu
        #
        #     self.nu_mass_numbers[self.nu_mass_eigenstates - 1] = num_massive_neutrinos
        #     self.nu_mass_degeneracies[self.nu_mass_eigenstates - 1] = neff_massive_standard
        #     self.nu_mass_fractions[self.nu_mass_eigenstates - 1] = (omnuh2 - omnuh2_sterile) / omnuh2
        # else:
        #     neff_massive_standard = 0
        if omnuh2_sterile > 0:
            if nnu < standard_neutrino_neff:
                raise CAMBError('nnu < 3.046 with massive sterile')
        # self.num_nu_massless = standard_neutrino_neff - neff_massive_standard
        #     self.num_nu_massive = self.num_nu_massive + 1
        #     self.nu_mass_eigenstates = self.nu_mass_eigenstates + 1
        #     self.nu_mass_numbers[self.nu_mass_eigenstates - 1] = 1
        #     self.nu_mass_degeneracies[self.nu_mass_eigenstates - 1] = max(1e-6, nnu - standard_neutrino_neff)
        #     self.nu_mass_fractions[self.nu_mass_eigenstates - 1] = omnuh2_sterile / omnuh2
        assert num_massive_neutrinos == int(num_massive_neutrinos)
        CAMB_SetNeutrinoHierarchy(byref(self), byref(c_double(omnuh2)), byref(c_double(omnuh2_sterile)),
                                  byref(c_double(nnu)), byref(c_int(neutrino_hierarchy)),
                                  byref(c_int(int(num_massive_neutrinos))))

        if tau is not None:
            self.Reion.set_tau(tau, delta_redshift=deltazrei)
        elif deltazrei:
            raise CAMBError('must set tau if setting deltazrei')

        return self

    def set_dark_energy(self, w=-1.0, sound_speed=1.0, dark_energy_model='fluid'):
        r"""
        Set dark energy parameters. Not that in this version these are not actually stored in
        the CAMBparams variable but set globally. So be careful!

        :param w: :math:`w\equiv p_{\rm de}/\rho_{\rm de}`, assumed constant
        :param sound_speed: rest-frame sound speed of dark energy fluid
        :param dark_energy_model: model to use, default is 'fluid'
        :return: self
        """
        # Variables from module LambdaGeneral
        if dark_energy_model != 'fluid':
            raise CAMBError('This version only supports the fluid energy model')
        if w != -1 or sound_speed != 1:
            logging.warning('Currently dark energy parameters are changed globally, not per parameter set')
        w_lam = dll_import(c_double, "lambdageneral", "w_lam")
        w_lam.value = w
        cs2_lam = dll_import(c_double, "lambdageneral", "cs2_lam")
        cs2_lam.value = sound_speed
        return self

    def get_omega_k(self):
        r"""
        Get curvature parameter :math:`\Omega_K`

        :return: :math:`\Omega_K`
        """
        return 1 - self.omegab - self.omegac - self.omegan - self.omegav

    def get_zre(self):
        if self.Reion.use_optical_depth:
            from . import camb
            return camb.get_zre_from_tau(self, self.Reion.optical_depth)
        else:
            return self.Reion.redshift

    def N_eff(self):
        """
        :return: Effective number of degrees of freedom in relativistic species at early times.
        """
        return sum(self.nu_mass_degeneracies[:self.nu_mass_eigenstates]) + self.num_nu_massless

    def get_Y_p(self, ombh2=None, delta_neff=None):
        r"""
        Get BBN helium nucleon fraction (NOT the same as the mass fraction Y_He) by intepolation using the
        :class:`.bbn.BBNPredictor` instance passed to :meth:`.model.CAMBparams.set_cosmology`
        (or the default one, if `Y_He` has not been set).

        :param ombh2: :math:`\Omega_b h^2` (default: value passed to :meth:`.model.CAMBparams.set_cosmology`)
        :param delta_neff:  additional :math:`N_{\rm eff}` relative to standard value (of 3.046) (default: from values passed to :meth:`.model.CAMBparams.set_cosmology`)
        :return:  :math:`Y_p^{\rm BBN}` helium nucleon fraction predicted by BBN.
        """
        try:
            ombh2 = ombh2 if ombh2 != None else self.omegab * (self.H0 / 100.) ** 2
            delta_neff = delta_neff if delta_neff != None else self.N_eff() - 3.046
            return self.bbn_predictor.Y_p(ombh2, delta_neff)
        except AttributeError:
            raise CAMBError('Not able to compute Y_p: not using an interpolation table for BBN abundances.')

    def get_DH(self, ombh2=None, delta_neff=None):
        r"""
        Get deuterium ration D/H by intepolation using the
        :class:`.bbn.BBNPredictor` instance passed to :meth:`.model.CAMBparams.set_cosmology`
        (or the default one, if `Y_He` has not been set).

        :param ombh2: :math:`\Omega_b h^2` (default: value passed to :meth:`.model.CAMBparams.set_cosmology`)
        :param delta_neff:  additional :math:`N_{\rm eff}` relative to standard value (of 3.046) (default: from values passed to :meth:`.model.CAMBparams.set_cosmology`)
        :return: BBN helium nucleon fraction D/H
        """
        try:
            ombh2 = ombh2 if ombh2 != None else self.omegab * (self.H0 / 100.) ** 2
            delta_neff = delta_neff if delta_neff != None else self.N_eff() - 3.046
            return self.bbn_predictor.DH(ombh2, delta_neff)
        except AttributeError:
            raise CAMBError('Not able to compute DH: not using an interpolation table for BBN abundances.')

    def set_matter_power(self, redshifts=[0.], kmax=1.2, k_per_logint=None, nonlinear=None,
                         accurate_massive_neutrino_transfers=False, silent=False):
        """
        Set parameters for calculating matter power spectra and transfer functions.

        :param redshifts: array of redshifts to calculate
        :param kmax: maximum k to calculate
        :param k_per_logint: number of k steps per log k. Set to zero to use default optimized spacing.
        :param nonlinear: if None, uses existing setting, otherwise boolean for whether to use non-linear matter power.
        :param accurate_massive_neutrino_transfers: if you want the massive neutrino transfers accurately
        :param silent: if True, don't give warnings about sort order
        :return: self
        """

        self.WantTransfer = True
        self.Transfer.high_precision = True
        self.Transfer.accurate_massive_neutrinos = accurate_massive_neutrino_transfers
        self.Transfer.kmax = kmax
        if nonlinear is not None:
            if nonlinear:
                if self.NonLinear in [NonLinear_lens, NonLinear_both]:
                    self.NonLinear = NonLinear_both
                else:
                    self.NonLinear = NonLinear_pk
            else:
                if self.NonLinear in [NonLinear_lens, NonLinear_both]:
                    self.NonLinear = NonLinear_lens
                else:
                    self.NonLinear = NonLinear_none
        if not k_per_logint:
            self.Transfer.k_per_logint = 0
        else:
            self.Transfer.k_per_logint = k_per_logint
        zs = sorted(redshifts, reverse=True)
        if not silent and np.any(np.array(zs) - np.array(redshifts) != 0):
            print("Note: redshifts have been re-sorted (earliest first)")
        if len(redshifts) > max_transfer_redshifts:
            raise CAMBError('You can have at most %s redshifts' % max_transfer_redshifts)
        self.Transfer.PK_num_redshifts = len(redshifts)
        for i, z in enumerate(zs):
            self.Transfer.PK_redshifts[i] = z
        return self

    def set_nonlinear_lensing(self, nonlinear):
        """
        Settings for whether or not to use non-linear corrections for the CMB lensing potential.
        Note that set_for_lmax also sets lensing to be non-linear if lens_potential_accuracy>0

        :param nonlinear: true to use non-linear corrections
        """
        if nonlinear:
            if self.NonLinear in [NonLinear_pk, NonLinear_both]:
                self.NonLinear = NonLinear_both
            else:
                self.NonLinear = NonLinear_lens
        else:
            if self.NonLinear in [NonLinear_pk, NonLinear_both]:
                self.NonLinear = NonLinear_pk
            else:
                self.NonLinear = NonLinear_none

    def set_for_lmax(self, lmax, max_eta_k=None, lens_potential_accuracy=0,
                     lens_margin=150, k_eta_fac=2.5, lens_k_eta_reference=18000.0):
        r"""
        Set parameters to get CMB power spectra accurate to specific a l_lmax.
        Note this does not fix the actual output L range, spectra may be calculated above l_max (but may not be accurate there).
        To fix the l_max for output arrays use the optional input argument to :meth:`.camb.CAMBdata.get_cmb_power_spectra` etc.

        :param lmax: :math:`\ell_{\rm max}` you want
        :param max_eta_k: maximum value of :math:`k \eta_0\approx k\chi_*` to use, which indirectly sets k_max. If None, sensible value set automatically.
        :param lens_potential_accuracy: Set to 1 or higher if you want to get the lensing potential accurate
        :param lens_margin: the :math:`\Delta \ell_{\rm max}` to use to ensure lensed :math:`C_\ell` are correct at :math:`\ell_{\rm max}`
        :param k_eta_fac:  k_eta_fac default factor for setting max_eta_k = k_eta_fac*lmax if max_eta_k=None
        :param lens_k_eta_reference:  value of max_eta_k to use when lens_potential_accuracy>0; use k_eta_max = lens_k_eta_reference*lens_potential_accuracy
        :return: self
        """
        if self.DoLensing:
            self.max_l = lmax + lens_margin
        else:
            self.max_l = lmax
        self.max_eta_k = max_eta_k or self.max_l * k_eta_fac
        if lens_potential_accuracy:
            self.set_nonlinear_lensing(True)
            self.max_eta_k = max(self.max_eta_k, lens_k_eta_reference * lens_potential_accuracy)
        return self

    def scalar_power(self, k):
        return self.primordial_power(k, 0)

    def tensor_power(self, k):
        return self.primordial_power(k, 2)

    def primordial_power(self, k, ix):
        if np.isscalar(k):
            karr = np.array([float(k)])
        else:
            karr = np.array(k)
        n = karr.shape[0]
        powers = np.empty(n)
        CAMB_primordialpower(byref(self), karr, powers, byref(c_int(n)), byref(c_int(ix)))
        if np.isscalar(k):
            return powers[0]
        else:
            return powers


def Transfer_SetForNonlinearLensing(P):
    camblib.__transfer_MOD_transfer_setfornonlinearlensing(byref(P))


def Transfer_SortAndIndexRedshifts(P):
    camblib.__transfer_MOD_transfer_sortandindexredshifts(byref(P))


CAMB_primordialpower.argtypes = [POINTER(CAMBparams), numpy_1d, numpy_1d, POINTER(c_int), POINTER(c_int)]

CAMB_SetNeutrinoHierarchy.argtypes = [POINTER(CAMBparams), POINTER(c_double), POINTER(c_double),
                                      POINTER(c_double), POINTER(c_int), POINTER(c_int)]
