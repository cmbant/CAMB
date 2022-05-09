from .baseconfig import camblib, CAMB_Structure, F2003Class, fortran_class, numpy_1d, np, \
    CAMBError, CAMBValueError, CAMBParamRangeError, AllocatableArrayDouble, AllocatableObject, \
    AllocatableObjectArray, AllocatableArrayInt, numpy_1d_int
from ctypes import c_bool, c_int, c_double, byref, POINTER, c_void_p
import ctypes
from . import reionization as reion
from . import recombination as recomb
from . import constants
from .initialpower import InitialPower, SplinedInitialPower
from .nonlinear import NonLinearModel
from .dark_energy import DarkEnergyModel, DarkEnergyEqnOfState
from .recombination import RecombinationModel
from .sources import SourceWindow
from . import bbn
import logging
from typing import Union, Optional

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

# for 21cm case
Transfer_monopole = 4
Transfer_vnewt = 5
Transfer_Tmat = 6

NonLinear_none = "NonLinear_none"
NonLinear_pk = "NonLinear_pk"
NonLinear_lens = "NonLinear_lens"
NonLinear_both = "NonLinear_both"
NonLinear_names = [NonLinear_none, NonLinear_pk, NonLinear_lens, NonLinear_both]

derived_names = ['age', 'zstar', 'rstar', 'thetastar', 'DAstar', 'zdrag',
                 'rdrag', 'kd', 'thetad', 'zeq', 'keq', 'thetaeq', 'thetarseq']

transfer_names = ['k/h', 'delta_cdm', 'delta_baryon', 'delta_photon', 'delta_neutrino', 'delta_nu', 'delta_tot',
                  'delta_nonu', 'delta_tot_de', 'Weyl', 'v_newtonian_cdm', 'v_newtonian_baryon', 'v_baryon_cdm']

evolve_names = transfer_names + ['a', 'etak', 'H', 'growth', 'v_photon', 'pi_photon',
                                 'E_2', 'v_neutrino', 'T_source', 'E_source', 'lens_potential_source']

background_names = ['x_e', 'opacity', 'visibility', 'cs2b', 'T_b', 'dopacity', 'ddopacity', 'dvisibility',
                    'ddvisibility']
density_names = ['tot', 'K', 'cdm', 'baryon', 'photon', 'neutrino', 'nu', 'de']

neutrino_hierarchy_normal = 'normal'
neutrino_hierarchy_inverted = 'inverted'
neutrino_hierarchy_degenerate = 'degenerate'
neutrino_hierarchies = [neutrino_hierarchy_normal, neutrino_hierarchy_inverted, neutrino_hierarchy_degenerate]


class TransferParams(CAMB_Structure):
    """
    Object storing parameters for the matter power spectrum calculation.

    """
    _fields_ = [
        ("high_precision", c_bool, "True for more accuracy"),
        ("accurate_massive_neutrinos", c_bool,
         "True if you want neutrino transfer functions accurate (false by default)"),
        ("kmax", c_double, "k_max to output (no h in units)"),
        ("k_per_logint", c_int, "number of points per log k interval. If zero, set an irregular optimized spacing"),
        ("PK_num_redshifts", c_int, "number of redshifts to calculate"),
        ("PK_redshifts", c_double * max_transfer_redshifts, {"size": "PK_num_redshifts"},
         "redshifts to output for the matter transfer and power"),
    ]


class AccuracyParams(CAMB_Structure):
    """
    Structure with parameters governing numerical accuracy. AccuracyBoost will also scale almost all the other
    parameters except for lSampleBoost (which is specific to the output interpolation) and lAccuracyBoost
    (which is specific to the multipole hierarchy evolution), e.g setting AccuracyBoost=2, IntTolBoost=1.5, means
    that internally the k sampling for integration will be boosed by AccuracyBoost*IntTolBoost = 3.
    """
    _fields_ = [
        ("AccuracyBoost", c_double, "general accuracy setting effecting everything related to step sizes etc. "
                                    "(including separate settings below except the next two)"),
        ("lSampleBoost", c_double,
         "accuracy for sampling in ell for interpolation for the C_l (if >=50, all ell are calculated)"),
        ("lAccuracyBoost", c_double, "Boosts number of multipoles integrated in Boltzman heirarchy"),
        ("AccuratePolarization", c_bool, "Do you care about the accuracy of the polarization Cls?"),
        ("AccurateBB", c_bool, "Do you care about BB accuracy (e.g. in lensing)"),
        ("AccurateReionization", c_bool, "Do you care about pecent level accuracy on EE signal from reionization?"),
        ("TimeStepBoost", c_double, "Sampling timesteps"),
        ("BackgroundTimeStepBoost", c_double,
         "Number of time steps for background thermal history and source window interpolation"),
        ("IntTolBoost", c_double, "Tolerances for integrating differential equations"),
        ("SourcekAccuracyBoost", c_double, "Accuracy of k sampling for source time integration"),
        ("IntkAccuracyBoost", c_double, "Accuracy of k sampling for integration"),
        ("TransferkBoost", c_double, "Accuracy of k sampling for transfer functions"),
        ("NonFlatIntAccuracyBoost", c_double, "Accuracy of non-flat time integration"),
        ("BessIntBoost", c_double, "Accuracy of bessel integration truncation"),
        ("LensingBoost", c_double, "Accuracy of the lensing of CMB power spectra"),
        ("NonlinSourceBoost", c_double, "Accuracy of steps and kmax used for the non-linear correction"),
        ("BesselBoost", c_double, "Accuracy of bessel pre-computation sampling"),
        ("LimberBoost", c_double, "Accuracy of Limber approximation use"),
        ("SourceLimberBoost", c_double, "Scales when to switch to Limber for source windows"),
        ("KmaxBoost", c_double, "Boost max k for source window functions"),
        ("neutrino_q_boost", c_double, "Number of momenta integrated for neutrino perturbations"),
    ]


class SourceTermParams(CAMB_Structure):
    """
    Structure with parameters determining how galaxy/lensing/21cm power spectra and transfer functions are calculated.
    """
    _fields_ = [
        ("limber_windows", c_bool,
         "Use Limber approximation where appropriate. CMB lensing uses Limber even if limber_window is false, " +
         "but method is changed to be consistent with other sources if limber_windows is true"),
        ("limber_phi_lmin", c_int, "When limber_windows=True, the minimum L to use Limber approximation for the "
                                   "lensing potential and other sources (which may use higher but not lower)"),
        ("counts_density", c_bool, "Include the density perturbation source"),
        ("counts_redshift", c_bool, "Include redshift distortions"),
        ("counts_lensing", c_bool, "Include magnification bias for number counts"),
        ("counts_velocity", c_bool, "Non-redshift distortion velocity terms"),
        ("counts_radial", c_bool, "Radial displacement velocity term; does not include time delay; "
                                  "subset of counts_velocity, just 1 / (chi * H) term"),
        ("counts_timedelay", c_bool, "Include time delay terms * 1 / (H * chi)"),
        ("counts_ISW", c_bool, "Include tiny ISW terms"),
        ("counts_potential", c_bool, "Include tiny terms in potentials at source"),
        ("counts_evolve", c_bool, "Accout for source evolution"),
        ("line_phot_dipole", c_bool, "Dipole sources for 21cm"),
        ("line_phot_quadrupole", c_bool, "Quadrupole sources for 21cm"),
        ("line_basic", c_bool, "Include main 21cm monopole density/spin temerature sources"),
        ("line_distortions", c_bool, "Redshift distortions for 21cm"),
        ("line_extra", c_bool, "Include other sources"),
        ("line_reionization", c_bool, "Replace the E modes with 21cm polarization"),
        ("use_21cm_mK", c_bool, "Use mK units for 21cm")]


class CustomSources(CAMB_Structure):
    """
    Structure containing symoblic-compiled custom CMB angular power spectrum source functions.
    Don't change this directly, instead call  :meth:`.model.CAMBparams.set_custom_scalar_sources`.
    """
    _fields_ = [("num_custom_sources", c_int, "number of sources set"),
                ("c_source_func", c_void_p, "Don't directly change this"),
                ("custom_source_ell_scales", AllocatableArrayInt, "scaling in L for outputs")]


@fortran_class
class CAMBparams(F2003Class):
    """
    Object storing the parameters for a CAMB calculation, including cosmological parameters and
    settings for what to calculate. When a new object is instantiated, default parameters are set automatically.

    To add a new parameter, add it to the CAMBparams type in model.f90, then  edit the _fields_ list in the CAMBparams
    class in model.py to add the new parameter in the corresponding location of the member list. After rebuilding the
    python version you can then access the parameter by using params.new_parameter_name where params is a CAMBparams
    instance. You could also modify the wrapper functions to set the field value less directly.

    You can view the set of underlying parameters used by the Fortran code by printing the CAMBparams instance.
    In python, to set cosmology parameters it is usually best to use :meth:`set_cosmology` and
    equivalent methods for most other parameters. Alternatively the convenience function :func:`.camb.set_params`
    can construct a complete instance from a dictionary of relevant parameters.

    """
    _fields_ = [
        ("WantCls", c_bool, "Calculate C_L"),
        ("WantTransfer", c_bool, "Calculate matter transfer functions and matter power spectrum"),
        ("WantScalars", c_bool, "Calculates scalar modes"),
        ("WantTensors", c_bool, "Calculate tensor modes"),
        ("WantVectors", c_bool, "Calculate vector modes"),
        ("WantDerivedParameters", c_bool, "Calculate derived parameters"),
        ("Want_cl_2D_array", c_bool, "For the C_L, include NxN matrix of all possible cross-spectra between sources"),
        ("Want_CMB", c_bool, "Calculate the temperature and polarization power spectra"),
        ("Want_CMB_lensing", c_bool, "Calculate the lensing potential power spectrum"),
        ("DoLensing", c_bool, "Include CMB lensing"),
        ("NonLinear", c_int, {"names": NonLinear_names}),
        ("Transfer", TransferParams),
        ("want_zstar", c_bool),
        ("want_zdrag", c_bool),
        ("min_l", c_int, "l_min for the scalar C_L (1 or 2, L=1 dipoles are Newtonian Gauge)"),
        ("max_l", c_int, "l_max for the scalar C_L"),
        ("max_l_tensor", c_int, "l_max for the tensor C_L"),
        ("max_eta_k", c_double, "Maximum k*eta_0 for scalar C_L, where eta_0 is the conformal time today"),
        ("max_eta_k_tensor", c_double, "Maximum k*eta_0 for tensor C_L, where eta_0 is the conformal time today"),
        ("ombh2", c_double, "Omega_baryon h^2"),
        ("omch2", c_double, "Omega_cdm h^2"),
        ("omk", c_double, "Omega_K"),
        ("omnuh2", c_double, "Omega_massive_neutrino h^2"),
        ("H0", c_double, "Hubble parameter is km/s/Mpc units"),
        ("TCMB", c_double, "CMB temperature today in Kelvin"),
        ("YHe", c_double, "Helium mass fraction"),
        ("num_nu_massless", c_double, "Effective number of massless neutrinos"),
        ("num_nu_massive", c_int, "Total physical (integer) number of massive neutrino species"),
        ("nu_mass_eigenstates", c_int, "Number of non-degenerate mass eigenstates"),
        ("share_delta_neff", c_bool, "Share the non-integer part of num_nu_massless between the eigenstates "),
        ("nu_mass_degeneracies", c_double * max_nu, {"size": "nu_mass_eigenstates"},
         "Degeneracy of each distinct eigenstate"),
        ("nu_mass_fractions", c_double * max_nu, {"size": "nu_mass_eigenstates"},
         "Mass fraction in each distinct eigenstate"),
        ("nu_mass_numbers", c_int * max_nu, {"size": "nu_mass_eigenstates"},
         "Number of physical neutrinos per distinct eigenstate"),
        ("InitPower", AllocatableObject(InitialPower)),
        ("Recomb", AllocatableObject(recomb.RecombinationModel)),
        ("Reion", AllocatableObject(reion.ReionizationModel)),
        ("DarkEnergy", AllocatableObject(DarkEnergyModel)),
        ("NonLinearModel", AllocatableObject(NonLinearModel)),
        ("Accuracy", AccuracyParams),
        ("SourceTerms", SourceTermParams),
        ("z_outputs", AllocatableArrayDouble, "redshifts to always calculate BAO output parameters"),
        ("scalar_initial_condition", c_int,
         {"names": ["initial_vector", "initial_adiabatic", "initial_iso_CDM", "initial_iso_baryon",
                    "initial_iso_neutrino", "initial_iso_neutrino_vel"]}),
        ("InitialConditionVector", AllocatableArrayDouble,
         "if scalar_initial_condition is initial_vector, the vector of initial condition amplitudes"),
        ("OutputNormalization", c_int, "If non-zero, multipole to normalize the C_L at"),
        ("Alens", c_double, "non-physical scaling amplitude for the CMB lensing spectrum power"),
        ("MassiveNuMethod", c_int, {"names": ["Nu_int", "Nu_trunc", "Nu_approx", "Nu_best"]}),
        ("DoLateRadTruncation", c_bool,
         "If true, use smooth approx to radiation perturbations after decoupling on small"
         " scales, saving evolution of irrelevant oscillatory multipole equations"),
        ("Evolve_baryon_cs", c_bool,
         "Evolve a separate equation for the baryon sound speed rather than using background approximation"),
        ("Evolve_delta_xe", c_bool, "Evolve ionization fraction perturbations"),
        ("Evolve_delta_Ts", c_bool, "Evolve the spin temperature perturbation (for 21cm)"),
        ("Do21cm", c_bool, "21cm is not yet implemented via the python wrapper"),
        ("transfer_21cm_cl", c_bool, "Get 21cm C_L at a given fixed redshift"),
        ("Log_lvalues", c_bool, "Use log spacing for sampling in L"),
        ("use_cl_spline_template", c_bool,
         "When interpolating use a fiducial spectrum shape to define ratio to spline"),
        ("SourceWindows", AllocatableObjectArray(SourceWindow)),
        ("CustomSources", CustomSources)
    ]

    _fortran_class_module_ = 'model'

    _methods_ = [('SetNeutrinoHierarchy', [POINTER(c_double), POINTER(c_double),
                                           POINTER(c_double), POINTER(c_int), POINTER(c_int)]),
                 ('Validate', None, c_int),
                 ('PrimordialPower', [numpy_1d, numpy_1d, POINTER(c_int), POINTER(c_int)]),
                 ('SetCustomSourcesFunc',
                  [POINTER(c_int), POINTER(ctypes.c_void_p), numpy_1d_int])]

    def __init__(self, **kwargs):
        set_default_params(self)
        self.InitPower.set_params()
        super().__init__(**kwargs)

    def validate(self):
        """
        Do some quick tests for sanity

        :return: True if OK
        """
        return self.f_Validate() != 0

    def set_accuracy(self, AccuracyBoost=1., lSampleBoost=1., lAccuracyBoost=1., DoLateRadTruncation=True):
        """
        Set parameters determining overall calculation accuracy (large values may give big slow down).
        For finer control you can set individual accuracy parameters by changing CAMBParams.Accuracy
        (:class:`.model.AccuracyParams`) .

        :param AccuracyBoost: increase AccuracyBoost to decrease integration step size, increase density of k
                              sampling, etc.
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

    def set_initial_power_function(self, P_scalar, P_tensor=None, kmin=1e-6, kmax=100., N_min=200, rtol=5e-5,
                                   effective_ns_for_nonlinear=None, args=()):
        r"""
        Set the initial power spectrum from a function P_scalar(k, \*args), and optionally also the tensor spectrum.
        The function is called to make a pre-computed array which is then interpolated inside CAMB. The sampling in k
        is set automatically so that the spline is accurate, but you may also need to increase other
        accuracy parameters.

        :param P_scalar: function returning normalized initial scalar curvature power as function of k (in Mpc^{-1})
        :param P_tensor: optional function returning normalized initial tensor power spectrum
        :param kmin: minimum wavenumber to compute
        :param kmax: maximum wavenumber to compute
        :param N_min: minimum number of spline points for the pre-computation
        :param rtol: relative tolerance for deciding how many points are enough
        :param effective_ns_for_nonlinear: an effective n_s for use with approximate non-linear corrections
        :param args: optional list of arguments passed to P_scalar (and P_tensor)
        :return: self
        """

        from scipy.interpolate import InterpolatedUnivariateSpline
        assert N_min > 7
        assert kmin < kmax
        # sample function logspace, finely enough that it interpolates accurately
        N = N_min
        ktest = np.logspace(np.log10(kmin), np.log10(kmax), N // 2)
        PK_test = P_scalar(ktest, *args)
        while True:
            ks = np.logspace(np.log10(kmin), np.log10(kmax), N)
            PK_compare = InterpolatedUnivariateSpline(ktest, PK_test)(ks)
            PK = P_scalar(ks, *args)
            if np.allclose(PK, PK_compare, atol=np.max(PK) * 1e-6, rtol=rtol):
                break
            N *= 2
            PK_test = PK
            ktest = ks
        PK_t = None if P_tensor is None else P_tensor(ks, *args)
        self.set_initial_power_table(ks, PK, PK_t, effective_ns_for_nonlinear)
        return self

    def set_initial_power_table(self, k, pk=None, pk_tensor=None, effective_ns_for_nonlinear=None):
        """
        Set a general intial power spectrum from tabulated values. It's up to you to ensure the sampling
        of the k values is high enough that it can be interpolated accurately.

        :param k: array of k values (Mpc^{-1})
        :param pk: array of primordial curvature perturbation power spectrum values P(k_i)
        :param pk_tensor: array of tensor spectrum values
        :param effective_ns_for_nonlinear: an effective n_s for use with approximate non-linear corrections
        """
        self.InitPower = SplinedInitialPower()
        initpower = self.InitPower
        if effective_ns_for_nonlinear is not None:
            initpower.effective_ns_for_nonlinear = effective_ns_for_nonlinear
        if pk is None:
            pk = np.empty(0)
        elif len(k) != len(pk):
            raise CAMBValueError("k and P(k) arrays must be same size")
        if pk_tensor is not None:
            if len(k) != len(pk_tensor):
                raise CAMBValueError("k and P_tensor(k) arrays must be same size")
            initpower.set_tensor_table(k, pk_tensor)
        initpower.set_scalar_table(k, pk)
        return self

    def set_initial_power(self, initial_power_params):
        """
        Set the InitialPower primordial power spectrum parameters

        :param initial_power_params: :class:`.initialpower.InitialPowerLaw`
                                     or :class:`.initialpower.SplinedInitialPower` instance
        :return: self
        """
        self.InitPower = initial_power_params
        return self

    def set_H0_for_theta(self, theta, cosmomc_approx=False, theta_H0_range=(10, 100), est_H0=67.0,
                         iteration_threshold=8):
        r"""
        Set H0 to give a specified value of the acoustic angular scale parameter theta.

        :param theta: value of :math:`r_s/D_M` at redshift :math:`z_\star`
        :param cosmomc_approx: if true, use approximate fitting formula for :math:`z_\star`,
                               if false do full numerical calculation
        :param theta_H0_range: min, max iterval to search for H0 (in km/s/Mpc)
        :param est_H0: an initial guess for H0 in km/s/Mpc, used in the case comsomc_approx=False.
        :param iteration_threshold: differnce in H0 from est_H0 for which to iterate, used for cosmomc_approx=False
        """

        if not (0.001 < theta < 0.1):
            raise CAMBParamRangeError('theta looks wrong (parameter is just theta, not 100*theta)')

        try:
            from scipy.optimize import brentq
        except ImportError:
            raise CAMBError('You need SciPy to set cosmomc_theta.')

        from . import camb

        data = camb.CAMBdata()
        if not cosmomc_approx:
            zstar = c_double()
            self.H0 = est_H0
            data.calc_background_no_thermo(self)
            # get_zstar initializes the recombination model
            zstar = data.f_get_zstar(byref(zstar))

        def f(H0):
            self.H0 = H0
            data.calc_background_no_thermo(self)
            if cosmomc_approx:
                theta_test = data.cosmomc_theta()
            else:
                rs = data.sound_horizon(zstar)
                theta_test = rs / (data.angular_diameter_distance(zstar) * (1 + zstar))
            return theta_test - theta

        try:
            # noinspection PyTypeChecker
            self.H0: float = brentq(f, theta_H0_range[0], theta_H0_range[1], rtol=5e-5)
            if not cosmomc_approx and abs(self.H0 - est_H0) > iteration_threshold:
                # iterate with recalculation of recombination and zstar
                self.set_H0_for_theta(theta, theta_H0_range=theta_H0_range, est_H0=self.H0,
                                      iteration_threshold=iteration_threshold)
        except ValueError:
            raise CAMBParamRangeError('No solution for H0 inside of theta_H0_range')

    def set_cosmology(self, H0: Optional[float] = None, ombh2=0.022, omch2=0.12, omk=0.0,
                      cosmomc_theta: Optional[float] = None, thetastar: Optional[float] = None,
                      neutrino_hierarchy: Union[str, int] = 'degenerate', num_massive_neutrinos=1,
                      mnu=0.06, nnu=constants.default_nnu, YHe: Optional[float] = None, meffsterile=0.0,
                      standard_neutrino_neff=constants.default_nnu, TCMB=constants.COBE_CMBTemp,
                      tau: Optional[float] = None, zrei: Optional[float] = None, deltazrei: Optional[float] = None,
                      Alens=1.0, bbn_predictor: Union[None, str, bbn.BBNPredictor] = None, theta_H0_range=(10, 100)):
        r"""
        Sets cosmological parameters in terms of physical densities and parameters (e.g. as used in Planck analyses).
        Default settings give a single distinct neutrino mass eigenstate, by default one neutrino with mnu = 0.06eV.
        Set the neutrino_hierarchy parameter to normal or inverted to use a two-eigenstate model that is a good
        approximation to the known mass splittings seen in oscillation measurements.
        For more fine-grained control can set the neutrino parameters directly rather than using this function.

        Instead of setting the Hubble parameter directly, you can instead set the acoustic scale parameter
        (cosmomc_theta, which is based on a fitting forumula for simple models, or thetastar, which is numerically
        calculated more generally). Note that you must have already set the dark energy model, you can't use
        set_cosmology with theta and then change the background evolution (which would change theta at the calculated
        H0 value).Likewise the dark energy model cannot depend explicitly on H0.

        :param H0: Hubble parameter today in km/s/Mpc. Can leave unset and instead set thetastar or cosmomc_theta
                  (which solves for the required H0).
        :param ombh2: physical density in baryons
        :param omch2:  physical density in cold dark matter
        :param omk: Omega_K curvature parameter
        :param cosmomc_theta: The approximate CosmoMC theta parameter :math:`\theta_{\rm MC}`. The angular
                              diamter distance is calculated numerically, but the redshift :math:`z_\star`
                              is calculated using an approximate (quite accurate but non-general) fitting formula.
                              Leave unset to use H0 or thetastar.
        :param thetastar: The angular acoustic scale parameter :math:`\theta_\star = r_s(z_*)/D_M(z_*)`, defined as
                    the ratio of the photon-baryon sound horizon :math:`r_s` to the angular diameter
                    distance :math:`D_M`, where both quantities are evaluated at :math:`z_*`, the redshift at
                    which the optical depth (excluding reionization) is unity. Leave unset to use H0 or cosmomc_theta.
        :param neutrino_hierarchy: 'degenerate', 'normal', or 'inverted' (1 or 2 eigenstate approximation)
        :param num_massive_neutrinos:  number of massive neutrinos
        :param mnu: sum of neutrino masses (in eV). Omega_nu is calculated approximately from this assuming neutrinos
               non-relativistic today; i.e. here is defined as a direct proxy for Omega_nu. Internally the actual
               physical mass is calculated from the Omega_nu accounting for small mass-dependent velocity corrections
               but neglecting spectral distortions to the neutrino distribution.
               Set the neutrino field values directly if you need finer control or more complex neutrino models.
        :param nnu: N_eff, effective relativistic degrees of freedom
        :param YHe: Helium mass fraction. If None, set from BBN consistency.
        :param meffsterile: effective mass of sterile neutrinos
        :param standard_neutrino_neff:  default value for N_eff in standard cosmology (non-integer to allow for partial
                heating of neutrinos at electron-positron annihilation and QED effects)
        :param TCMB: CMB temperature (in Kelvin)
        :param tau: optical depth; if None and zrei is None, current Reion settings are not changed
        :param zrei: reionization mid-point optical depth (set tau=None to use this)
        :param deltazrei: redshift width of reionization; if None, uses default
        :param Alens: (non-physical) scaling of the lensing potential compared to prediction
        :param bbn_predictor: :class:`.bbn.BBNPredictor` instance used to get YHe from BBN consistency if YHe is None,
         or name of a BBN predictor class, or file name of an interpolation table
        :param theta_H0_range: if thetastar or cosmomc_theta is specified, the min, max interval of H0 values to map to;
          if H0 is outside this range it will raise an exception.
        """

        if YHe is None:
            # use BBN prediction
            if isinstance(bbn_predictor, str):
                self.bbn_predictor = bbn.get_predictor(bbn_predictor)
            else:
                self.bbn_predictor = bbn_predictor or bbn.get_predictor()
            YHe = self.bbn_predictor.Y_He(ombh2 * (constants.COBE_CMBTemp / TCMB) ** 3, nnu - standard_neutrino_neff)
        self.YHe = YHe
        self.TCMB = TCMB
        self.ombh2 = ombh2
        self.omch2 = omch2
        self.Alens = Alens

        neutrino_mass_fac = constants.neutrino_mass_fac * (constants.COBE_CMBTemp / TCMB) ** 3

        if not isinstance(neutrino_hierarchy, str):
            neutrino_hierarchy = neutrino_hierarchies[neutrino_hierarchy - 1]

        if nnu >= standard_neutrino_neff or neutrino_hierarchy != neutrino_hierarchy_degenerate:
            omnuh2 = mnu / neutrino_mass_fac * (standard_neutrino_neff / 3) ** 0.75
        else:
            omnuh2 = mnu / neutrino_mass_fac * (nnu / 3.0) ** 0.75
        omnuh2_sterile = meffsterile / neutrino_mass_fac
        if omnuh2_sterile > 0 and nnu < standard_neutrino_neff:
            raise CAMBError('sterile neutrino mass required Neff> %.3g' % constants.default_nnu)
        if omnuh2 and not num_massive_neutrinos:
            raise CAMBError('non-zero mnu with zero num_massive_neutrinos')

        omnuh2 = omnuh2 + omnuh2_sterile
        self.omnuh2 = omnuh2
        self.omk = omk
        if omnuh2_sterile > 0:
            if nnu < standard_neutrino_neff:
                raise CAMBError('nnu < %.3g with massive sterile' % constants.default_nnu)
        assert num_massive_neutrinos == int(num_massive_neutrinos)
        self.f_SetNeutrinoHierarchy(byref(c_double(omnuh2)), byref(c_double(omnuh2_sterile)),
                                    byref(c_double(nnu)),
                                    byref(c_int(neutrino_hierarchies.index(neutrino_hierarchy) + 1)),
                                    byref(c_int(int(num_massive_neutrinos))))

        if cosmomc_theta or thetastar:
            if H0 is not None:
                raise CAMBError('Set H0=None when setting theta.')
            if cosmomc_theta and thetastar:
                raise CAMBError('Cannot set both cosmomc_theta and thetastar')

            self.set_H0_for_theta(cosmomc_theta or thetastar, cosmomc_approx=cosmomc_theta is not None,
                                  theta_H0_range=theta_H0_range)
        else:
            if H0 is None:
                raise CAMBError('Must set H0, cosmomc_theta or thetastar')
            if H0 < 1:
                raise CAMBValueError('H0 is the value in km/s/Mpc, your value looks very small')
            self.H0 = H0

        if tau is not None:
            if zrei is not None:
                raise CAMBError('Cannot set both tau and zrei')
            self.Reion.set_tau(tau, delta_redshift=deltazrei)
        elif zrei is not None:
            self.Reion.set_zrei(zrei, delta_redshift=deltazrei)
        elif deltazrei:
            raise CAMBError('must set tau if setting deltazrei')

        return self

    @property
    def h(self):
        return self.H0 / 100

    @h.setter
    def h(self, value):
        self.H0 = value * 100

    @property
    def omegab(self):
        return self.ombh2 / (self.H0 / 100) ** 2

    @property
    def omegac(self):
        return self.omch2 / (self.H0 / 100) ** 2

    @property
    def omeganu(self):
        return self.omnuh2 / (self.H0 / 100) ** 2

    @property
    def omegam(self):
        return (self.ombh2 + self.omch2 + self.omnuh2) / (self.H0 / 100) ** 2

    @property
    def N_eff(self):
        """
        :return: Effective number of degrees of freedom in relativistic species at early times.
        """
        if self.share_delta_neff:
            return self.num_nu_massless + self.num_nu_massive
        else:
            return sum(self.nu_mass_degeneracies[:self.nu_mass_eigenstates]) + self.num_nu_massless

    def set_classes(self, dark_energy_model=None, initial_power_model=None,
                    non_linear_model=None, recombination_model=None):
        """
        Change the classes used to implement parts of the model.

        :param dark_energy_model: 'fluid', 'ppf', or name of a DarkEnergyModel class
        :param initial_power_model: name of an InitialPower class
        :param non_linear_model: name of a NonLinearModel class
        :param recombination_model: name of recombination_model class
        """
        if dark_energy_model:
            self.DarkEnergy = self.make_class_named(dark_energy_model, DarkEnergyModel)
        if initial_power_model:
            self.InitPower = self.make_class_named(initial_power_model, InitialPower)
        if non_linear_model:
            self.NonLinear = self.make_class_named(non_linear_model, NonLinearModel)
        if recombination_model:
            self.Recomb = self.make_class_named(recombination_model, RecombinationModel)

    def set_dark_energy(self, w=-1.0, cs2=1.0, wa=0, dark_energy_model='fluid'):
        r"""
        Set dark energy parameters (use set_dark_energy_w_a to set w(a) from numerical table instead)
        To use a custom dark energy model, assign the class instance to the DarkEnergy field instead.

        :param w: :math:`w\equiv p_{\rm de}/\rho_{\rm de}`, assumed constant
        :param wa: evolution of w (for dark_energy_model=ppf)
        :param cs2: rest-frame sound speed squared of dark energy fluid
        :param dark_energy_model: model to use ('fluid' or 'ppf'), default is 'fluid'
        :return: self
        """

        de = self.make_class_named(dark_energy_model, DarkEnergyEqnOfState)
        de.set_params(w=w, wa=wa, cs2=cs2)
        self.DarkEnergy = de
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
        self.DarkEnergy = self.make_class_named(dark_energy_model, DarkEnergyEqnOfState)
        # Note that assigning to allocatable fields makes deep copies of the object
        self.DarkEnergy.set_w_a_table(a, w)
        return self

    def get_zre(self):
        return self.Reion.get_zre(self)

    # alias consistent with input parameter name
    get_zrei = get_zre

    def get_Y_p(self, ombh2=None, delta_neff=None):
        r"""
        Get BBN helium nucleon fraction (NOT the same as the mass fraction Y_He) by intepolation using the
        :class:`.bbn.BBNPredictor` instance passed to :meth:`set_cosmology`
        (or the default one, if `Y_He` has not been set).

        :param ombh2: :math:`\Omega_b h^2` (default: value passed to :meth:`set_cosmology`)
        :param delta_neff:  additional :math:`N_{\rm eff}` relative to standard value (of 3.046)
                           (default: from values passed to :meth:`set_cosmology`)
        :return:  :math:`Y_p^{\rm BBN}` helium nucleon fraction predicted by BBN.
        """
        try:
            ombh2 = ombh2 if ombh2 is not None else self.ombh2
            delta_neff = delta_neff if delta_neff is not None else self.N_eff - constants.default_nnu
            return self.bbn_predictor.Y_p(ombh2, delta_neff)
        except AttributeError:
            raise CAMBError('Not able to compute Y_p: not using an interpolation table for BBN abundances.')

    def get_DH(self, ombh2=None, delta_neff=None):
        r"""
        Get deuterium ration D/H by intepolation using the
        :class:`.bbn.BBNPredictor` instance passed to :meth:`set_cosmology`
        (or the default one, if `Y_He` has not been set).

        :param ombh2: :math:`\Omega_b h^2` (default: value passed to :meth:`set_cosmology`)
        :param delta_neff:  additional :math:`N_{\rm eff}` relative to standard value (of 3.046)
                           (default: from values passed to :meth:`set_cosmology`)
        :return: BBN helium nucleon fraction D/H
        """
        try:
            ombh2 = ombh2 if ombh2 is not None else self.ombh2
            delta_neff = delta_neff if delta_neff is not None else self.N_eff - constants.default_nnu
            return self.bbn_predictor.DH(ombh2, delta_neff)
        except AttributeError:
            raise CAMBError('Not able to compute DH: not using an interpolation table for BBN abundances.')

    def set_matter_power(self, redshifts=(0.,), kmax=1.2, k_per_logint=None, nonlinear=None,
                         accurate_massive_neutrino_transfers=False, silent=False):
        """
        Set parameters for calculating matter power spectra and transfer functions.

        :param redshifts: array of redshifts to calculate
        :param kmax: maximum k to calculate (where k is just k, not k/h)
        :param k_per_logint: minimum number of k steps per log k. Set to zero to use default optimized spacing.
        :param nonlinear: if None, uses existing setting, otherwise boolean for whether to use non-linear matter power.
        :param accurate_massive_neutrino_transfers: if you want the massive neutrino transfers accurately
        :param silent: if True, don't give warnings about sort order
        :return: self
        """
        if not len(redshifts):
            raise CAMBError('set_matter_power redshifts list is empty')

        self.WantTransfer = True
        self.Transfer.high_precision = True
        self.Transfer.accurate_massive_neutrinos = accurate_massive_neutrino_transfers
        self.Transfer.kmax = kmax
        zs = sorted(redshifts, reverse=True)
        if nonlinear is not None:
            if nonlinear:
                if self.NonLinear in [NonLinear_lens, NonLinear_both]:
                    self.NonLinear = NonLinear_both
                else:
                    self.NonLinear = NonLinear_pk
                if not silent and (kmax < 5 or kmax < 20 and np.max(zs) > 4):
                    logging.warning("Using kmax=%s with Halofit non-linear models may give inaccurate results" % kmax)
            else:
                if self.NonLinear in [NonLinear_lens, NonLinear_both]:
                    self.NonLinear = NonLinear_lens
                else:
                    self.NonLinear = NonLinear_none
        self.Transfer.k_per_logint = k_per_logint if k_per_logint else 0
        if not silent and np.any(np.array(zs) - np.array(redshifts) != 0):
            print("Note: redshifts have been re-sorted (earliest first)")
        if len(redshifts) > max_transfer_redshifts:
            raise CAMBError('You can have at most %s redshifts' % max_transfer_redshifts)
        self.Transfer.PK_redshifts = zs
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
                     lens_margin=150, k_eta_fac=2.5, lens_k_eta_reference=18000.0, nonlinear=None):
        r"""
        Set parameters to get CMB power spectra accurate to specific a l_lmax.
        Note this does not fix the actual output L range, spectra may be calculated above l_max
        (but may not be accurate there). To fix the l_max for output arrays use the optional input argument
        to :meth:`.results.CAMBdata.get_cmb_power_spectra` etc.

        :param lmax: :math:`\ell_{\rm max}` you want
        :param max_eta_k: maximum value of :math:`k \eta_0\approx k\chi_*` to use, which indirectly sets k_max.
                          If None, sensible value set automatically.
        :param lens_potential_accuracy: Set to 1 or higher if you want to get the lensing potential accurate
                                        (1 is only Planck-level accuracy)
        :param lens_margin: the :math:`\Delta \ell_{\rm max}` to use to ensure lensed :math:`C_\ell` are correct
                            at :math:`\ell_{\rm max}`
        :param k_eta_fac:  k_eta_fac default factor for setting max_eta_k = k_eta_fac*lmax if max_eta_k=None
        :param lens_k_eta_reference:  value of max_eta_k to use when lens_potential_accuracy>0; use
                                      k_eta_max = lens_k_eta_reference*lens_potential_accuracy
        :param nonlinear: use non-linear power spectrum; if None, sets nonlinear if lens_potential_accuracy>0 otherwise
                          preserves current setting
        :return: self
        """
        if self.DoLensing:
            self.max_l = lmax + lens_margin
        else:
            self.max_l = lmax
        self.max_eta_k = max_eta_k or self.max_l * k_eta_fac
        if lens_potential_accuracy:
            self.set_nonlinear_lensing(nonlinear is not False)
            self.max_eta_k = max(self.max_eta_k, lens_k_eta_reference * lens_potential_accuracy)
        elif nonlinear is not None:
            self.set_nonlinear_lensing(nonlinear)
        return self

    def scalar_power(self, k):
        r"""
        Get the primordial scalar curvature power spectrum at :math:`k`

        :param k: wavenumber :math:`k` (in :math:`{\rm Mpc}^{-1}` units)
        :return: power spectrum at :math:`k`
        """
        return self.primordial_power(k, 0)

    def tensor_power(self, k):
        r"""
        Get the primordial tensor curvature power spectrum at :math:`k`

        :param k:  wavenumber :math:`k` (in :math:`{\rm Mpc}^{-1}` units)
        :return: tensor power spectrum at :math:`k`
        """

        return self.primordial_power(k, 2)

    def primordial_power(self, k, ix):
        if np.isscalar(k):
            karr = np.array([k], dtype=np.float64)
        else:
            karr = np.array(k, dtype=np.float64)
        n = karr.shape[0]
        powers = np.empty(n)
        self.f_PrimordialPower(karr, powers, byref(c_int(n)), byref(c_int(ix)))
        if np.isscalar(k):
            return powers[0]
        else:
            return powers

    _custom_source_name_dict = {}

    def set_custom_scalar_sources(self, custom_sources, source_names=None, source_ell_scales=None,
                                  frame='CDM', code_path=None):
        r"""
        Set custom sources for angular power spectrum using camb.symbolic sympy expressions.

        :param custom_sources: list of sympy expressions for the angular power spectrum sources
        :param source_names: optional list of string naes for the sources
        :param source_ell_scales: list or dictionary of scalings for each source name, where for integer entry n,
            the source for multipole :math:`\ell` is scalled by :math:`\sqrt{(\ell+n)!/(\ell-n)!}`,
            i.e. :math:`n=2` for a new polarization-like source.
        :param frame: if the source is not gauge invariant, frame in which to interpret result
        :param code_path: optional path for output of source code for CAMB f90 source function
        """

        from . import symbolic

        if isinstance(custom_sources, dict):
            assert (not source_names)
            if source_ell_scales and not isinstance(source_ell_scales, dict):
                raise CAMBValueError('source_ell_scales must be a dictionary if custom_sources is')
            lst = []
            source_names = []
            for name in custom_sources.keys():
                source_names.append(name)
                lst.append(custom_sources[name])
            custom_sources = lst
        elif not isinstance(custom_sources, (list, tuple)):
            custom_sources = [custom_sources]
            if source_names:
                source_names = [source_names]
        custom_source_names = source_names or ["C%s" % (i + 1) for i in range(len(custom_sources))]
        if len(custom_source_names) != len(custom_sources):
            raise CAMBValueError('Number of custom source names does not match number of sources')
        scales = np.zeros(len(custom_sources), dtype=np.int32)
        if source_ell_scales:
            if isinstance(source_ell_scales, dict):
                if set(source_ell_scales.keys()) - set(custom_source_names):
                    raise CAMBValueError('scale dict key not in source names list')
                for i, name in enumerate(custom_source_names):
                    if name in source_ell_scales:
                        scales[i] = source_ell_scales[name]
            else:
                scales[:] = source_ell_scales

        _current_source_func = symbolic.compile_sympy_to_camb_source_func(custom_sources, frame=frame,
                                                                          code_path=code_path)

        custom_source_func = ctypes.cast(_current_source_func, c_void_p)
        self._custom_source_name_dict[custom_source_func.value] = custom_source_names
        self.f_SetCustomSourcesFunc(byref(c_int(len(custom_sources))), byref(custom_source_func), scales)

    def get_custom_source_names(self):
        if self.CustomSources.num_custom_sources:
            return self._custom_source_name_dict[self.CustomSources.c_source_func]
        else:
            return []

    def clear_custom_scalar_sources(self):
        self.f_SetCustomSourcesFunc(byref(c_int(0)), byref(ctypes.c_void_p(0)), np.zeros(0, dtype=np.int32))

    def diff(self, params):
        """
        Print differences between this set of parameters and params

        :param params: another CAMBparams instance
        """
        p1 = str(params)
        p2 = str(self)
        for line1, line2 in zip(p1.split('\n'), p2.split('\n')):
            if line1 != line2:
                print(line1, ' <-> ', line2)


def set_default_params(P):
    """
    Set default values for all parameters

    :param P: :class:`.model.CAMBparams`
    :return: P
    """
    assert (isinstance(P, CAMBparams))
    camblib.__camb_MOD_camb_setdefparams(byref(P))
    return P
