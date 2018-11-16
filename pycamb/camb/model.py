from .baseconfig import camblib, CAMB_Structure, F2003Class, \
    CAMBError, CAMBValueError, CAMBParamRangeError, dll_import, f_allocatable
from ctypes import c_bool, c_int, c_double, byref, POINTER
from . import reionization as ion
from . import recombination as recomb
from . import constants
from .nonlinear import NonLinearModel
from .initialpower import InitialPower, SplinedInitialPower
import numpy as np
from . import bbn
import ctypes
import six

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
density_names = ['tot', 'K', 'cdm', 'baryon', 'photon', 'neutrino', 'nu', 'de']

dark_energy_models = ['fluid', 'ppf']

neutrino_hierarchies = ['normal', 'inverted', 'degenerate']
neutrino_hierarchy_normal = 1
neutrino_hierarchy_inverted = 2
neutrino_hierarchy_degenerate = 3

# ---Variables in modules.f90
# To set the value please just put
# variable_name.value = new_value

DebugParam = dll_import(c_double, "cambsettings", "debugparam")
# DebugParam.value = 1000000*2

# logical
do_bispectrum = dll_import(c_int, "cambsettings", "do_bispectrum")
# do_bispectrum.value = False

limber_phiphi = dll_import(c_int, "cambsettings", "limber_phiphi")
# limber_phiphi.value = 0

num_extra_redshiftwindows = dll_import(c_int, "cambsettings", "num_extra_redshiftwindows")
# num_extra_redshiftwindows.value = 0

num_redshiftwindows = dll_import(c_int, "cambsettings", "num_redshiftwindows")

num_custom_sources = dll_import(c_int, "cambsettings", "num_custom_sources")

# logical
use_spline_template = dll_import(c_bool, "cambsettings", "use_spline_template")
# use_spline_template.value = True

# Variable from module Transfer
# logical
transfer_interp_matterpower = dll_import(c_bool, "transfer", "transfer_interp_matterpower")
# transfer_interp_matterpower.value = False

transfer_power_var = dll_import(c_int, "transfer", "transfer_power_var")
# transfer_power_var.value = Transfer_tot

CAMB_validateparams = camblib.__camb_MOD_camb_validateparams
CAMB_validateparams.restype = c_bool

# args for these set below after CAMBparams defined

CAMB_SetNeutrinoHierarchy = camblib.__camb_MOD_camb_setneutrinohierarchy

numpy_1d = np.ctypeslib.ndpointer(c_double, flags='C_CONTIGUOUS')
CAMB_primordialpower = camblib.__handles_MOD_camb_primordialpower
CAMB_primordialpower.restype = c_bool

CAMBparams_SetDarkEnergy = camblib.__handles_MOD_cambparams_setdarkenergy
CAMBparams_GetAllocatables = camblib.__handles_MOD_cambparams_getallocatables

CAMBparams_DE_SetTable = camblib.__handles_MOD_cambparams_setdarkenergytable
CAMBparams_DE_GetStressEnergy = camblib.__handles_MOD_cambparams_darkenergystressenergy


class TransferParams(CAMB_Structure):
    """
    Object storing parameters for the matter power spectrum calculation. PK variables are for setting main outputs.
    Other entries are used internally, e.g. for sampling to get correct non-linear corrections and lensing.

    :ivar high_precision: True for more accuracy
    :ivar accurate_massive_neutrinos: True if you want neutrino transfer functions accurate (false by default)
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


class DarkEnergyParams(CAMB_Structure):
    """
    Object storing dark energy parameters. Set w and wa to use w(a) = w + (1-a)*wa parameterization,
    or call set_w_a_table to set another tabulated w(a). If tabulated w(a) is used, w and wa are set
    to approximate values at z=0.

    :ivar w: equation of state parameter at z=0
    :ivar wa: evolution of w, using default parameterization with  w(a) = w + (1-a)*wa
    :ivar cs2: rest-frame sound speed of dark energy perturbations
    :ivar use_tabulated_w: if using an explicitly set w(a) function rather than default w + (1-a)*wa

    """
    _fields_ = [
        ("w", c_double),
        ("wa", c_double),
        ("cs2", c_double),
        ("use_tabulated_w", c_int),
        ("no_perturbations", c_int)
    ]

    def set_w_a_table(self, a, w):
        """
        Set w(a) from numerical values (used as cublic spline)

        :param a: array of scale factors
        :param w: array of w(a)
        :return: self
        """
        if len(a) != len(w): raise ValueError('Dark energy w(a) table non-equal sized arrays')
        if not np.isclose(a[-1], 1):  raise ValueError('Dark energy w(a) arrays must end at a=1')
        CAMBparams_DE_SetTable(byref(self), a, w, byref(c_int(len(a))))
        return self


class AccuracyParams(CAMB_Structure):
    _fields_ = [
        ("lSampleBoost", c_double),
        ("AccuracyBoost", c_double),
        ("lAccuracyBoost", c_double),
        ("AccuratePolarization", c_int),  # logical
        ("AccurateBB", c_int),  # logical
        ("AccurateReionization", c_int),  # logical
        ("TimeStepBoost", c_double),
        ("IntTolBoost", c_double),
        ("SourcekAccuracyBoost", c_double),
        ("IntkAccuracyBoost", c_double),
        ("TransferkBoost", c_double),
        ("NonFlatIntAccuracyBoost", c_double),
        ("BessIntBoost", c_double),
        ("BesselBoost", c_double),
        ("LimberBoost", c_double),
        ("Kmax_Boost", c_double),
        ("neutrino_q_boost", c_double),
        ("thermo_boost", c_double)
    ]


class SourceTermParams(CAMB_Structure):
    _fields_ = [
        ("limber_windows", c_int),
        ("do_counts_lensing", c_int),
        ("line_phot_dipole", c_int),
        ("line_phot_quadrupole", c_int),
        ("line_basic", c_int),
        ("line_extra", c_int),
        ("line_distortions", c_int),
        ("line_reionization", c_int),
        ("counts_velocity", c_int),
        ("counts_radial", c_int),  # does not include time delay; subset of counts_velocity, just 1 / (chi * H) term
        ("counts_density", c_int),
        ("counts_redshift", c_int),
        ("counts_timedelay", c_int),  # time delay terms * 1 / (H * chi)
        ("counts_ISW", c_int),
        ("counts_potential", c_int),  # terms in potentials at source
        ("counts_evolve", c_int),
        ("use_21cm_mK", c_int)]


class CAMBparams(F2003Class):
    """
    Object storing the parameters for a CAMB calculation, including cosmological parameters and
    settings for what to calculate. When a new object is instantiated, default parameters are set automatically.

    You can view the set of underlying parameters used by the Fortran code by printing the CAMBparams instance.

    To add a new parameter, add it to the CAMBparams type in modules.f90, then  edit the _fields_ list in the CAMBparams
    class in model.py to add the new parameter in the corresponding location of the member list. After rebuilding the
    python version you can then access the parameter by using params.new_parameter where params is a CAMBparams instance.
    You could also modify the wrapper functions to set the field value less directly.
    """

    _fields_ = [
        ("WantCls", c_int),  # logical
        ("WantTransfer", c_int),  # logical
        ("WantScalars", c_int),  # logical
        ("WantTensors", c_int),  # logical
        ("WantVectors", c_int),  # logical
        ("want_cl_2D_array", c_int),  # logical
        ("DoLensing", c_int),  # logical
        ("Do21cm", c_int),  # logical
        ("want_zstar", c_int),  # logical
        ("want_zdrag", c_int),  # logical
        ("PK_WantTransfer", c_int),  # logical
        ("NonLinear", c_int),
        ("Want_CMB", c_int),  # logical
        ("Want_CMB_lensing", c_int),  # logical
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
        ("Alens", c_double),
        ("MassiveNuMethod", c_int),
        ("Reion", ion.ReionizationParams),
        ("Recomb", recomb.RecombinationParams),
        ("Transfer", TransferParams),
        ("Accuracy", AccuracyParams),
        ("SourceTerms", SourceTermParams),
        ("DoLateRadTruncation", c_int),
        ("Evolve_baryon_cs", c_int),
        ("Evolve_delta_xe", c_int),
        ("Evolve_delta_Ts", c_int),
        ("InitialConditionVector", c_double * 10),
        ("OnlyTransfers", c_int),  # logical
        ("transfer_21cm_cl", c_int),  # logical
        ("DerivedParameters", c_int),  # logical
        ("Log_lvalues", c_int),  # boolean
        ("_InitPower", f_allocatable),  # resolved class accessed as self.InitPower
        ("_DarkEnergy", f_allocatable),  # resolved class accessed as self.DarkEnergy
        ("_NonLinearModel", f_allocatable),  # resolved class accessed as self.NonLinearModel
    ]

    def __init__(self):
        self._init_members()
        self.InitPower.set_params()

    def _init_members(self):
        de = POINTER(DarkEnergyParams)()
        nonlin = POINTER(NonLinearModel)()
        initpower = ctypes.c_void_p()

        i = c_int(0)
        j = c_int(0)
        CAMBparams_GetAllocatables(byref(self), byref(i), byref(de), byref(nonlin), byref(j), byref(initpower))
        self.DarkEnergy = de.contents
        self.NonLinearModel = nonlin.contents
        self.InitPower = InitialPower.void_contents(initpower, j.value)

    def validate(self):
        """
        Do some quick tests for sanity

        :return: True if OK
        """
        return CAMB_validateparams(byref(self))

    def set_accuracy(self, AccuracyBoost=1., lSampleBoost=1., lAccuracyBoost=1., DoLateRadTruncation=True):
        """
        Set parameters determining calculation accuracy (large values may give big slow down).
        Note curently these are set globally, not just per parameter set.

        :param AccuracyBoost: increase AccuracyBoost to decrease integration step size, increase density of k sampling, etc.
        :param lSampleBoost: increase lSampleBoost to increase density of L sampling for CMB
        :param lAccuracyBoost: increase lAccuracyBoost to increase the maximum L included in the Boltzmann hierarchies
        :param DoLateRadTruncation: If True, use approximation to radiation perturbation evolution at late times
        :return: self
        """
        self.Accuracy.lSampleBoost = lSampleBoost
        self.Accuracy.AccuracyBoost = AccuracyBoost
        self.Accuracy.lAccuracyBoost = lAccuracyBoost
        self.DoLateRadTruncation = DoLateRadTruncation
        return self

    def set_initial_power_function(self, P_scalar, P_tensor=None, kmin=1e-6, kmax=100., N_min=200, rtol=5e-5, args=()):
        """
        Set the initial power spectrum from a function P_scalar(k, *args), and optionally also the tensor spectrum.
        The function is called to make a pre-computed array which is then interpolated inside CAMB. The sampling in k
        is set automatically so that the spline is accurate, but you may also need to increase other accuracy parameters.

        :param P_scalar: function returning normalized initial scalar curvature power as function of k (in Mpc^{-1})
        :param P_tensor: optional function returning normalized initial tensor power spectrum
        :param kmin: minimum wavenumber to compute
        :param kmax: maximum wavenumber to compute
        :param N_min: minimum number of spline points for the pre-computation
        :param rtol: relative tolerance for deciding how many points are enough
        :param args: optional list of arguments passed to P_scalar (and P_tensor)
        :return: self
        """

        from scipy.interpolate import UnivariateSpline
        assert N_min > 7
        assert kmin < kmax
        # sample function logspace, finely enough that it interpolates accurately
        N = N_min
        ktest = np.logspace(np.log10(kmin), np.log10(kmax), N // 2)
        PK_test = P_scalar(ktest, *args)
        while True:
            ks = np.logspace(np.log10(kmin), np.log10(kmax), N)
            PK_compare = UnivariateSpline(ktest, PK_test, s=0)(ks)
            PK = P_scalar(ks, *args)
            if np.allclose(PK, PK_compare, atol=np.max(PK) * 1e-6, rtol=rtol):
                break
            N *= 2
            PK_test = PK
            ktest = ks
        PK_t = None if P_tensor is None else P_tensor(ks, *args)
        self.set_initial_power_table(ks, PK, PK_t)
        return self

    def set_initial_power_table(self, k, pk=None, pk_tensor=None):
        """
        Set a general intial power spectrum from tabulated values. It's up to you to ensure the sampling
        of the k values is high enough that it can be interpolated accurately.

        :param k: array of k values (Mpc^{-1})
        :param pk: array of primordial curvature perturbation power spectrum values P(k_i)
        param pk_tensor: array of tensor spectrum values
        """
        initpower = POINTER(SplinedInitialPower)()
        if pk is None:
            pk = np.asarray([])
        elif len(k) != len(pk):
            raise CAMBValueError("k and P(k) arrays must be same size")
        if pk_tensor is None:
            pk_tensor = np.asarray([])
        elif len(k) != len(pk_tensor):
            raise CAMBValueError("k and P_tensor(k) arrays must be same size")
        self.call_func('setpktable', extra_args=[POINTER(c_int), POINTER(c_int), numpy_1d, numpy_1d, numpy_1d,
                                                 POINTER(POINTER(SplinedInitialPower))], args=[
            byref(c_int(len(pk))), byref(c_int(len(pk_tensor))), k, pk, pk_tensor, byref(initpower)])
        self.InitPower = initpower.contents
        return self

    def set_initial_power(self, initial_power_params):
        """
        Set the InitialPower primordial power spectrum parameters

        :param initial_power_params: :class:`.initialpower.InitialPowerLaw` or :class:`.initialpower.SplinedInitialPower`instance
        :return: self
        """
        assert isinstance(initial_power_params, InitialPower)
        p, class_id = initial_power_params._pointer_id()
        self.call_func(tag='setinitialpower', extra_args=[POINTER(ctypes.c_void_p), POINTER(c_int)], args=[
            byref(p), byref(c_int(class_id))])
        self._init_members()
        return self

    def set_cosmology(self, H0=67.0, cosmomc_theta=None, ombh2=0.022, omch2=0.12, omk=0.0,
                      neutrino_hierarchy='degenerate', num_massive_neutrinos=1,
                      mnu=0.06, nnu=3.046, YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046,
                      TCMB=constants.COBE_CMBTemp, tau=None, deltazrei=None, Alens=1.0,
                      bbn_predictor=None, theta_H0_range=[10, 100]):
        r"""
        Sets cosmological parameters in terms of physical densities and parameters used in Planck 2015 analysis.
        Default settings give a single distinct neutrino mass eigenstate, by default one neutrino with mnu = 0.06eV.
        Set the neutrino_hierarchy parameter to normal or inverted to use a two-eigenstate model that is a good
        approximation to the known mass splittings seen in oscillation measurements.
        If you require more fine-grained control can set the neutrino parameters directly rather than using this function.

        :param cosmomc_theta: The CosmoMC theta parameter :math:`\theta_{\rm MC}`. You must set H0=None to solve for H0 given cosmomc_theta. Note that
                you must have already set the dark energy model, you can't use set_cosmology with cosmomc_theta and then
                change the background evolution (which would change cosmomc_theta at the calculated H0 value). Likewise the dark energy model
                cannot depend on H0.
        :param ombh2: physical density in baryons
        :param omch2:  physical density in cold dark matter
        :param omk: Omega_K curvature parameter
        :param neutrino_hierarchy: 'degenerate', 'normal', or 'inverted' (1 or 2 eigenstate approximation)
        :param num_massive_neutrinos:  number of massive neutrinos
        :param mnu: sum of neutrino masses (in eV)
        :param nnu: N_eff, effective relativistic degrees of freedom
        :param YHe: Helium mass fraction. If None, set from BBN consistency.
        :param meffsterile: effective mass of sterile neutrinos
        :param standard_neutrino_neff:  default value for N_eff in standard cosmology (non-integer to allow for partial
                heating of neutrinos at electron-positron annihilation and QED effects)
        :param TCMB: CMB temperature (in Kelvin)
        :param tau: optical depth; if None, current Reion settings are not changed
        :param deltazrei: redshift width of reionization; if None, uses default
        :param Alens: (non-physical) scaling of the lensing potential compared to prediction
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
        self.Alens = Alens

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

    def set_dark_energy(self, w=-1.0, cs2=1.0, wa=0, dark_energy_model='fluid'):
        r"""
        Set dark energy parameters (use set_dark_energy_w_a to set w(a) from numerical table instead)

        :param w: :math:`w\equiv p_{\rm de}/\rho_{\rm de}`, assumed constant
        :param wa: evolution of w (for dark_energy_model=ppf)
        :param cs2: rest-frame sound speed squared of dark energy fluid
        :param dark_energy_model: model to use ('fluid' or 'ppf'), default is 'fluid'
        :return: self
        """

        if dark_energy_model == 'fluid' and wa and (w < -1 - 1e-6 or 1 + w + wa < -1 - 1e-6):
            raise CAMBError('fluid dark energy model does not support w crossing -1')

        de = POINTER(DarkEnergyParams)()
        CAMBparams_SetDarkEnergy(byref(self), byref(c_int(dark_energy_models.index(dark_energy_model))), byref(de))
        self.DarkEnergy = de.contents
        self.DarkEnergy.w = w
        self.DarkEnergy.wa = wa
        self.DarkEnergy.cs2 = cs2
        self.DarkEnergy.use_tabulated_w = False

        return self

    def set_dark_energy_w_a(self, a, w, dark_energy_model='fluid'):
        """
        Set the dark energy equation of state from tabulated values (which are cubic spline interpolated).

        :param a: array of sampled a = 1/(1+z) values
        :param w: array of w(a)
        :param dark_energy_model:  model to use ('fluid' or 'ppf'), default is 'fluid'
        :return: self
        """
        if dark_energy_model == 'fluid' and np.any(w < -1):
            raise CAMBError('fluid dark energy model does not support w crossing -1')
        de = POINTER(DarkEnergyParams)()
        CAMBparams_SetDarkEnergy(byref(self), byref(c_int(dark_energy_models.index(dark_energy_model))), byref(de))
        self.DarkEnergy = de.contents
        self.DarkEnergy.set_w_a_table(a, w)

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

    def get_dark_energy_rho_w(self, a):
        """
        get dark energy density in units of the dark energy density today, and w=P/rho
        :param a: scalar factor or array of scale factors
        :return: rho, w arrays at redshifts 1/a-1 [or scalars if a is scalar]
        """
        if np.isscalar(a):
            scales = np.array([a])
        else:
            scales = np.asarray(a)
        rho = np.zeros(scales.shape)
        w = np.zeros(scales.shape)
        CAMBparams_DE_GetStressEnergy(byref(self), scales, rho, w, byref(c_int(len(scales))))
        if np.isscalar(a):
            return rho[0], w[0]
        else:
            return rho, w

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

    def set_z_outputs(self, z_outputs):
        """
        Set the redshifts for calculating BAO parameters at

        :param z_outputs: array of redshifts
        """
        z_outputs = np.array(z_outputs, dtype=np.float64)
        self.call_func('setbackgroundoutputs_z', extra_args=[numpy_1d, POINTER(c_int)],
                       args=[z_outputs, byref(c_int(len(z_outputs)))])


CAMB_primordialpower.argtypes = [POINTER(CAMBparams), numpy_1d, numpy_1d, POINTER(c_int), POINTER(c_int)]
CAMBparams_SetDarkEnergy.argtypes = [POINTER(CAMBparams), POINTER(c_int), POINTER(POINTER(DarkEnergyParams))]
CAMBparams_GetAllocatables.argtypes = CAMBparams_SetDarkEnergy.argtypes + [POINTER(POINTER(NonLinearModel))] + \
                                      [POINTER(c_int), POINTER(ctypes.c_void_p)]

CAMBparams_DE_SetTable.argtypes = [POINTER(DarkEnergyParams), numpy_1d, numpy_1d, POINTER(c_int)]

CAMBparams_DE_GetStressEnergy.argtypes = [POINTER(CAMBparams), numpy_1d, numpy_1d, numpy_1d, POINTER(c_int)]

CAMB_SetNeutrinoHierarchy.argtypes = [POINTER(CAMBparams), POINTER(c_double), POINTER(c_double),
                                      POINTER(c_double), POINTER(c_int), POINTER(c_int)]
