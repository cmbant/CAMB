from .baseconfig import camblib, CAMBError, CAMBValueError, CAMBUnknownArgumentError, CAMB_Structure, F2003Class, \
    fortran_class, numpy_1d, numpy_2d, numpy_1d_int, fortran_array, AllocatableArrayDouble, ndpointer, np, lib_import
from ctypes import c_float, c_int, c_double, c_bool, POINTER, byref
import ctypes
from . import model, constants
from ._config import config
from .model import set_default_params, CAMBparams
import logging
from scipy.interpolate import RectBivariateSpline, interp1d

int_arg = POINTER(c_int)
d_arg = POINTER(c_double)


class _MatterTransferData(CAMB_Structure):
    # contains complex types with pointers, so just set up dummy
    _fields_ = [('num_q_trans', c_int),
                ('q_trans', POINTER(c_double)),
                ('sigma_8', POINTER(c_double)),
                ('sigma2_vdelta_8', POINTER(c_double)),
                ('TransferData', POINTER(c_float)),
                ('sigma_8_size', c_int),
                ('sigma2_vdelta_8_size', c_int),
                ('TransferData_size', c_int * 3)
                ]


class _ClTransferData(CAMB_Structure):
    _fields_ = [('NumSources', c_int),
                ('q_size', c_int),
                ('q', POINTER(c_double)),
                ('delta_size', c_int * 3),
                ('delta_p_l_k', POINTER(c_double)),
                ('l_size', c_int),
                ('L', POINTER(c_int))
                ]


def save_cmb_power_array(filename, array, labels, lmin=0):
    """
    Save a zero-based 2-d array of CL to a text file, with each line starting with L.

    :param filename: filename to save
    :param array: 2D array of power spectra
    :param labels:  header names for each column in the output
    :param lmin: L to start output in file (usually 0 or 2)
    """
    lmax = array.shape[0] - 1
    ls = np.atleast_2d(np.arange(lmin, lmax + 1)).T
    ncol = array.shape[1]
    if isinstance(labels, str):
        labels = labels.split()
    # noinspection PyTypeChecker
    np.savetxt(filename, np.hstack((ls, array[lmin:, :])), fmt=['%4u'] + ['%12.7e'] * ncol,
               header=' L ' + ' '.join(['{:13s}'.format(lab) for lab in labels]))


class MatterTransferData:
    r"""
    MatterTransferData is the base class for storing matter power transfer function data for various q values.
    In a flat universe q=k, in a closed universe q is quantized.

    To get an instance of this data, call :meth:`.results.CAMBdata.get_matter_transfer_data`.

    For a description of the different Transfer_xxx outputs (and 21cm case) see :ref:`transfer-variables`; the
    array is indexed by index+1 gven by:

    - Transfer_kh = 1 (k/h)
    - Transfer_cdm = 2 (cdm)
    - Transfer_b = 3 (baryons)
    - Transfer_g = 4 (photons)
    - Transfer_r = 5 (massless neutrinos)
    - Transfer_nu = 6 (massive neutrinos)
    - Transfer_tot = 7 (total matter)
    - Transfer_nonu = 8 (total matter excluding neutrinos)
    - Transfer_tot_de = 9 (total including dark energy perturbations)
    - Transfer_Weyl = 10 (Weyl potential)
    - Transfer_Newt_vel_cdm = 11 (Newtonian CDM velocity)
    - Transfer_Newt_vel_baryon = 12 (Newtonian baryon velocity)
    - Transfer_vel_baryon_cdm = 13 (relative baryon-cdm velocity)


    :ivar nq:  number of q modes calculated
    :ivar q: array of q values calculated
    :ivar sigma_8: array of :math:`\sigma_8` values for each redshift
    :ivar sigma2_vdelta_8: array of v-delta8 correlation, so sigma2_vdelta_8/sigma_8 can define growth
    :ivar transfer_data: numpy array T[entry, q_index, z_index] storing transfer functions for each redshift and q;
                         entry+1 can be one of the Transfer_xxx variables above.
    """

    nq: int
    q: np.ndarray
    sigma_8: np.ndarray
    sigma2_vdelta_8: np.ndarray
    transfer_data: np.ndarray

    def transfer_z(self, name, z_index=0):
        """
        Get transfer function (function of q, for each q in self.q_trans) by name for given redshift index

        :param name:  parameter name
        :param z_index: which redshift
        :return: array of transfer function values for each calculated k
        """

        if name not in model.transfer_names:
            raise CAMBError('Unknown name %s; must be one of %s' % (name, model.transfer_names))
        return self.transfer_data[model.transfer_names.index(name), :, z_index]


class ClTransferData:
    r"""
    ClTransferData is the base class for storing CMB power transfer functions, as a function of q and :math:`\ell`.
    To get an instance of this data, call :meth:`.results.CAMBdata.get_cmb_transfer_data`

    :ivar NumSources:  number of sources calculated (size of p index)
    :ivar q: array of q values calculated (=k in flat universe)
    :ivar L: int array of :math:`\ell` values calculated
    :ivar delta_p_l_k: transfer functions, indexed by source, L, q
    """
    NumSources: int
    q: np.ndarray
    L: np.ndarray
    delta_p_l_k: np.ndarray

    def get_transfer(self, source=0):
        r"""
        Return :math:`C_\ell` trasfer functions as a function of :math:`\ell`
        and :math:`q` (:math:`= k` in a flat universe).

        :param source: index of source: e.g. 0 for temperature, 1 for E polarization, 2 for lensing potential
        :return: array of computed L, array of computed q, transfer functions T(L,q)
        """

        return self.L, self.q, self.delta_p_l_k[source, :, :]


@fortran_class
class CAMBdata(F2003Class):
    """
    An object for storing calculational data, parameters and transfer functions.
    Results for a set of parameters (given in a :class:`~.model.CAMBparams` instance)
    are returned by the :func:`.camb.get_background`, :func:`.camb.get_transfer_functions` or :func:`.camb.get_results`
    functions. Exactly which quantities are already calculated depends on which of these functions you use and the
    input parameters.

    To quickly make a fully calculated CAMBdata instance for a set of parameters you can call :func:`.camb.get_results`.

    """
    _fortran_class_module_ = 'results'

    _fields_ = [("Params", CAMBparams),
                ("ThermoDerivedParams", c_double * model.nthermo_derived,
                 "array of derived parameters, see :meth:`get_derived_params` to get as a dictionary"),
                ("flat", c_bool, "flat universe"),
                ("closed", c_bool, "closed universe"),
                # grho gives the contribution to the expansion rate from: (g) photons,
                # (r) one flavor of relativistic neutrino (2 degrees of freedom),
                #  grho is actually 8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
                ("grhocrit", c_double, "kappa*a^2*rho_c(0)/c^2 with units of Mpc**(-2)"),
                ("grhog", c_double, "kappa/c^2*4*sigma_B/c^3 T_CMB^4"),
                ("grhor", c_double, "7/8*(4/11)^(4/3)*grhog (per massless neutrino species)"),
                ("grhob", c_double, "baryon contribution"),
                ("grhoc", c_double, "CDM contribution"),
                ("grhov", c_double, "Dark energy contribution"),
                ("grhornomass", c_double, "grhor*number of massless neutrino species"),
                ("grhok", c_double, "curvature contribution to critical density"),
                ("taurst", c_double, "time at start of recombination"),
                ("dtaurec", c_double, "time step in recombination"),
                ("taurend", c_double, "time at end of recombination"),
                ("tau_maxvis", c_double, "time at peak visibility"),
                ("adotrad", c_double, "da/d tau in early radiation-dominated era"),
                ("omega_de", c_double, "Omega for dark energy today"),
                ("curv", c_double, "curvature K"),
                ("curvature_radius", c_double, r":math:`1/\sqrt{|K|}`"),
                ("Ksign", c_double, "Ksign = 1,0 or -1"),
                ("tau0", c_double, "conformal time today"),
                ("chi0", c_double, "comoving angular diameter distance of big bang; rofChi(tau0/curvature_radius)"),
                ("scale", c_double, "relative to flat. e.g. for scaling L sampling"),
                ("akthom", c_double, "sigma_T * (number density of protons now)"),
                ("fHe", c_double, "n_He_tot / n_H_tot"),
                ("Nnow", c_double, "number density today"),
                ("z_eq", c_double, "matter-radiation equality redshift assuming all neutrinos relativistic"),
                ("grhormass", c_double * model.max_nu),
                ("nu_masses", c_double * model.max_nu),
                ("num_transfer_redshifts", c_int,
                 "Number of calculated redshift outputs for the matter transfer (including those for CMB lensing)"),
                ("transfer_redshifts", AllocatableArrayDouble, "Calculated output redshifts"),
                ("PK_redshifts_index", c_int * model.max_transfer_redshifts, "Indices of the requested PK_redshifts"),
                ("OnlyTransfers", c_bool, "Only calculating transfer functions, not power spectra"),
                ("HasScalarTimeSources", c_bool, "calculate and save time source functions, not power spectra")]

    # Note there are many more fields in Fortran. Since F2003Class is memory-managed by Fortran, we don't need
    # need to define them all in python.
    # _methods_ refer to the imported functions in the corresponding fortran class.

    _methods_ = [('AngularDiameterDistance', [d_arg], c_double),
                 ('AngularDiameterDistanceArr', [numpy_1d, numpy_1d, int_arg]),
                 ('AngularDiameterDistance2', [d_arg, d_arg], c_double),
                 ('AngularDiameterDistance2Arr', [numpy_1d, numpy_1d, numpy_1d, int_arg]),
                 ('ComovingRadialDistance', [d_arg], c_double),
                 ('ComovingRadialDistanceArr', [numpy_1d, numpy_1d, int_arg, d_arg]),
                 ('Hofz', [d_arg], c_double),
                 ('HofzArr', [numpy_1d, numpy_1d, int_arg]),
                 ('DeltaPhysicalTimeGyr', [d_arg, d_arg, d_arg], c_double),
                 ('DeltaPhysicalTimeGyrArr', [numpy_1d, numpy_1d, numpy_1d, int_arg, d_arg]),
                 ('GetBackgroundDensities', [int_arg, numpy_1d, numpy_2d]),
                 ('DeltaTime', [d_arg, d_arg, d_arg], c_double),
                 ('DeltaTimeArr', [numpy_1d, numpy_1d, numpy_1d, int_arg, d_arg]),
                 ('TimeOfzArr', [numpy_1d, numpy_1d, int_arg, d_arg]),
                 ('sound_horizon_zArr', [numpy_1d, numpy_1d, int_arg]),
                 ('RedshiftAtTimeArr', [numpy_1d, numpy_1d, int_arg]),
                 ('CosmomcTheta', [], c_double),
                 ('DarkEnergyStressEnergy', [numpy_1d, numpy_1d, numpy_1d, int_arg]),
                 ('get_lmax_lensed', [], c_int),
                 ('get_zstar', [d_arg], c_double),
                 ('SetParams', [POINTER(CAMBparams), int_arg, int_arg, int_arg, int_arg])
                 ]

    def __init__(self):
        set_default_params(self.Params)

    def set_params(self, params):
        """
        Set parameters from params. Note that this does not recompute anything;
        you will need to call :meth:`calc_transfers` if you change any parameters affecting the
        background cosmology or the transfer function settings.

        :param params: a :class:`~.model.CAMBparams` instance
        """
        self.Params = params

    def get_derived_params(self):
        """
        :return: dictionary of derived parameter values, indexed by name ('kd', 'age', etc..)
        """
        res = {}
        for name, value in zip(model.derived_names, self.ThermoDerivedParams):
            res[name] = value
        return res

    def get_background_outputs(self):
        """
        Get BAO values for redshifts set in Params.z_outputs

        :return: rs/DV, H, DA, F_AP for each requested redshift (as 2D array)
        """
        n = len(self.Params.z_outputs)
        if not n:
            raise CAMBError('Set z_outputs with required redshifts (and then calculate transfers/results)'
                            ' before calling get_background_outputs')
        outputs = np.empty((n, 4))
        CAMB_GetBackgroundOutputs(byref(self), outputs, byref(c_int(n)))
        return outputs

    def get_BAO(self, redshifts, params):
        """
        Get BAO parameters at given redshifts, using parameters in params

        :param redshifts: list of redshifts
        :param params: optional :class:`~.model.CAMBparams` instance to use
        :return: array of rs/DV, H, DA, F_AP for each redshift as 2D array
        """
        P = params.copy()
        P.z_outputs = redshifts
        self.calc_background(P)
        return self.get_background_outputs()

    @staticmethod
    def _check_params(params):
        if not isinstance(params, CAMBparams):
            raise CAMBValueError('Must pass a CAMBparams instance')
        if not params.ombh2:
            raise CAMBValueError('Parameter values not set')

    def calc_background_no_thermo(self, params, do_reion=False):
        """
        Calculate the background evolution without calculating thermal or ionization history.
        e.g. call this if you want to just use :meth:`angular_diameter_distance` and similar background functions

        :param params:  :class:`~.model.CAMBparams` instance to use
        :param do_reion: whether to initialize the reionization model
        """
        self._check_params(params)
        self.f_SetParams(byref(params), None, byref(c_int(1 if do_reion else 0)), None, byref(c_int(1)))
        config.check_global_error('calc_background_no_thermo')

    def calc_background(self, params):
        """
        Calculate the background evolution and thermal history.
        e.g. call this if you want to get derived parameters and call background functions

        :param params:  :class:`~.model.CAMBparams` instance to use
        """
        self._check_params(params)
        if CAMBdata_CalcBackgroundTheory(byref(self), byref(params)):
            config.check_global_error('calc_background')

    def calc_transfers(self, params, only_transfers=True, only_time_sources=False):
        """
        Calculate the transfer functions (for CMB and matter power, as determined
        by params.WantCls, params.WantTransfer).

        :param params: :class:`~.model.CAMBparams` instance with parameters to use
        :param only_transfers: only calculate transfer functions, no power spectra
        :param only_time_sources: only calculate time transfer functions, no (p,l,k) transfer functions or
                                    non-linear scaling
        :return: non-zero if error, zero if OK
        """
        self._check_params(params)
        if not (only_transfers or only_time_sources):
            self._check_powers(params)
        if CAMBdata_gettransfers(byref(self), byref(params), byref(c_int(1 if only_transfers else 0)),
                                 byref(c_int(1 if only_time_sources else 0))):
            config.check_global_error('calc_transfer')

    def _check_powers(self, params=None):
        if params is None:
            params = self.Params
        if params.InitPower.has_tensors() and not params.WantTensors:
            raise CAMBError('r>0 but params.WantTensors = F')
        if params.WantScalars and params.WantCls and params.DoLensing and params.scalar_power(0.05) > 2e-8:
            raise CAMBError('Lensing requires a realistically normalized spectrum, you have P(k=0.05/Mpc) > 2e-8')

    def calc_power_spectra(self, params=None):
        """
        Calculates transfer functions and power spectra.

        :param params: optional :class:`~.model.CAMBparams` instance with parameters to use

        """
        if params is not None:
            self.calc_transfers(params, only_transfers=False)
        else:
            self._check_powers()
            CAMBdata_transferstopowers(byref(self))
            config.check_global_error()

    def power_spectra_from_transfer(self, initial_power_params=None, silent=False):
        """
        Assuming :meth:`calc_transfers` or :meth:`calc_power_spectra` have already been used, re-calculate the
        power spectra using a new set of initial power spectrum parameters with otherwise the same cosmology.
        This is typically much faster that re-calculating everything, as the transfer functions can be re-used.
        NOTE: if non-linear lensing is on, the transfer functions have the non-linear correction included when
        they are calculated, so using this function with a different initial power spectrum will not give quite the
        same results as doing a full recalculation.

        :param initial_power_params: :class:`.initialpower.InitialPowerLaw` or :class:`.initialpower.SplinedInitialPower`
               instance with new primordial power spectrum parameters, or None to use current power spectrum.
        :param silent: suppress warnings about non-linear corrections not being recalculated
        """
        if not silent and not self.HasScalarTimeSources and \
                self.Params.NonLinear in [model.NonLinear_lens, model.NonLinear_both] and \
                self.Params.WantScalars and self.Params.WantCls and not getattr(self, '_suppress_power_warn', False):
            logging.warning(
                'power_spectra_from_transfer with non-linear lensing does not recalculate the non-linear correction')
            self._suppress_power_warn = True
        if initial_power_params:
            self.Params.set_initial_power(initial_power_params)
        self._check_powers()
        CAMBdata_transferstopowers(byref(self))
        config.check_global_error()

    def _CMB_unit(self, CMB_unit):
        if isinstance(CMB_unit, str):
            if CMB_unit == 'muK':
                CMB_unit = self.Params.TCMB * 1e6
            elif CMB_unit == 'K':
                CMB_unit = self.Params.TCMB
            else:
                raise CAMBValueError('Unknown CMB_unit: %s' % CMB_unit)
        return CMB_unit

    def _scale_cls(self, cls, CMB_unit=None, raw_cl=False, lens_potential=False):
        if raw_cl:
            ls = np.arange(1, cls.shape[0])[..., np.newaxis]
            ls = np.float64(ls * (ls + 1))
            if lens_potential:
                cls[1:, 0] /= ls[:, 0] ** 2 / (2 * np.pi)
                cls[1:, 1:] /= ls ** (3. / 2) / (2 * np.pi)
            else:
                cls[1:, :] /= ls / (2 * np.pi)

        if CMB_unit is not None:
            CMB_unit = self._CMB_unit(CMB_unit)
            if lens_potential:
                cls[:, 1:] *= CMB_unit
            else:
                cls *= CMB_unit ** 2

        return cls

    def _lmax_setting(self, lmax=None, unlensed=False):
        if self.Params.DoLensing and not unlensed:
            lmax_calc = self.f_get_lmax_lensed()
            if not lmax_calc:
                raise CAMBError('lensed CL have not been calculated')
        else:
            lmax_calc = self.Params.max_l
        if lmax is None:
            lmax = lmax_calc
        elif lmax > lmax_calc:
            logging.warning('getting CMB power spectra to higher L than calculated, may be innacurate/zeroed.')
        return lmax

    def save_cmb_power_spectra(self, filename, lmax=None, CMB_unit='muK'):
        r"""
        Save CMB power to a plain text file. Output is lensed total :math:`\ell(\ell+1)C_\ell/2\pi` then
        lensing potential and cross: L TT EE BB TE PP PT PE.

        :param filename: filename to save
        :param lmax: lmax to save
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
               and :math:`\mu K` units for lensing cross.
        """
        lmax = self._lmax_setting(lmax)
        cmb = self.get_total_cls(lmax, CMB_unit=CMB_unit)
        lens = self.get_lens_potential_cls(lmax, CMB_unit=CMB_unit)
        save_cmb_power_array(filename, np.hstack((cmb, lens)), 'TT EE BB TE PP PT PE')

    def get_cmb_power_spectra(self, params=None, lmax=None,
                              spectra=('total', 'unlensed_scalar', 'unlensed_total', 'lensed_scalar', 'tensor',
                                       'lens_potential'), CMB_unit=None, raw_cl=False):
        r"""
        Get CMB power spectra, as requested by the 'spectra' argument. All power spectra are
        :math:`\ell(\ell+1)C_\ell/2\pi` self-owned numpy arrays (0..lmax, 0..3), where 0..3 index
        are TT, EE, BB, TE, unless raw_cl is True in which case return just :math:`C_\ell`.
        For the lens_potential the power spectrum returned is that of the deflection.

        Note that even if lmax is None, all spectra a returned to the same lmax, appropriate
        for lensed spectra. Use the individual functions instead if you want to the full unlensed
        and lensing potential power spectra to the higher lmax actually computed.

        :param params: optional :class:`~.model.CAMBparams` instance with parameters to use. If None, must have
          previously set parameters and called `calc_power_spectra` (e.g. if you got this instance
          using :func:`.camb.get_results`),
        :param lmax: maximum L
        :param spectra: list of names of spectra to get
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
          and :math:`\mu K` units for lensing cross.
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: dictionary of power spectrum arrays, indexed by names of requested spectra
        """
        P = {}
        if params is not None:
            self.calc_power_spectra(params)
        lmax = self._lmax_setting(lmax)
        for spectrum in spectra:
            P[spectrum] = getattr(self, 'get_' + spectrum + '_cls')(lmax, CMB_unit=CMB_unit,
                                                                    raw_cl=raw_cl)
        return P

    def get_cmb_correlation_functions(self, params=None, lmax=None, spectrum='lensed_scalar',
                                      xvals=None, sampling_factor=1):
        r"""
        Get the CMB correlation functions from the power spectra.
        By default evaluated at points :math:`\cos(\theta)` = xvals that are roots of Legendre polynomials,
        for accurate back integration with :func:`.correlations.corr2cl`.
        If xvals is explicitly given, instead calculates correlations at provided :math:`\cos(\theta)` values.

        :param params: optional :class:`~.model.CAMBparams` instance with parameters to use. If None, must have
          previously set parameters and called :meth:`calc_power_spectra` (e.g. if you got this instance
          using :func:`.camb.get_results`),
        :param lmax: optional maximum L to use from the cls arrays
        :param spectrum: type of CMB power spectrum to get; default 'lensed_scalar', one of
          ['total', 'unlensed_scalar', 'unlensed_total', 'lensed_scalar', 'tensor']
        :param xvals: optional array of :math:`\cos(\theta)` values at which to calculate correlation function.
        :param sampling_factor: multiple of lmax for the Gauss-Legendre order if xvals not given (default 1)
        :return: if xvals not given: corrs, xvals, weights; if xvals specified, just corrs.
          corrs is 2D array corrs[i, ix], where ix=0,1,2,3 are T, Q+U, Q-U and cross, and i indexes xvals
        """

        if spectrum not in ['total', 'unlensed_scalar', 'unlensed_total', 'lensed_scalar', 'tensor']:
            raise CAMBValueError('Can only get CMB correlation functions for known CMB spectrum')
        from . import correlations

        cls = self.get_cmb_power_spectra(params, lmax, spectra=[spectrum])[spectrum]
        if xvals is None:
            return correlations.gauss_legendre_correlation(cls, sampling_factor=sampling_factor)
        else:
            return correlations.cl2corr(cls, xvals, lmax=lmax)

    def get_cmb_transfer_data(self, tp='scalar'):
        r"""
        Get :math:`C_\ell` transfer functions

        :return: :class:`.ClTransferData` instance holding output arrays (copies, not pointers)
        """

        cdata = _ClTransferData()
        CAMBdata_cltransferdata(byref(self), byref(cdata),
                                byref(c_int(['scalar', 'vector', 'tensor'].index(tp))))
        data = ClTransferData()
        data.NumSources = cdata.NumSources
        data.q = fortran_array(cdata.q, cdata.q_size)
        data.L = fortran_array(cdata.L, cdata.l_size, dtype=c_int)
        data.delta_p_l_k = fortran_array(cdata.delta_p_l_k, cdata.delta_size)
        return data

    def get_time_evolution(self, q, eta, vars=model.evolve_names, lAccuracyBoost=4, frame='CDM'):
        """
        Get the mode evolution as a function of conformal time for some k values.

        Note that gravitational potentials (e.g. Weyl) are not integrated in the code and are
        calculated as derived parameters; they may be numerically unstable far outside the horizon.
        (use the series expansion result if needed far outside the horizon)

        :param q: wavenumber values to calculate (or array of k values)
        :param eta: array of requested conformal times to output
        :param vars: list of variable names or sympy symbolic expressions to output (using camb.symbolic)
        :param lAccuracyBoost: factor by which to increase l_max in hierarchies compared to default - often
          needed to get nice smooth curves of acoustic oscillations for plotting.
        :param frame: for symbolic expressions, can specify frame name if the variable is not gauge invariant.
            e.g. specifying Delta_g and frame='Newtonian' would give the Newtonian gauge photon density perturbation.
        :return: nd array, A_{qti}, size(q) x size(times) x len(vars), or 2d array if q is scalar
        """

        old_boost = self.Params.Accuracy.lAccuracyBoost
        try:
            if lAccuracyBoost:
                self.Params.Accuracy.lAccuracyBoost = lAccuracyBoost
            if not isinstance(vars, (tuple, list)):
                vars = [vars]
            import sympy
            named_vars = [var for var in vars if isinstance(var, str)]

            unknown = set(named_vars) - set(model.evolve_names)
            if unknown:
                raise CAMBError('Unknown names %s; valid names are %s' % (unknown, model.evolve_names))

            num_standard_names = len(model.evolve_names)

            custom_vars = []
            ix = np.empty(len(vars), dtype=int)
            for i, var in enumerate(vars):
                if var in model.evolve_names:
                    ix[i] = model.evolve_names.index(var)
                elif isinstance(var, sympy.Expr):
                    custom_vars.append(var)
                    ix[i] = num_standard_names + len(custom_vars) - 1
                else:
                    raise CAMBError(
                        'Variables must be variable names, or a sympy expression (using camb.symbolic variables)')

            if np.isscalar(q):
                k = np.array([q], dtype=np.float64)
            else:
                k = np.array(q, dtype=np.float64)
            times = np.array(np.atleast_1d(eta), dtype=np.float64)
            indices = np.argsort(times)  # times must be in increasing order
            ncustom = len(custom_vars)
            if ncustom:
                from . import symbolic
                funcPtr = symbolic.compile_sympy_to_camb_source_func(custom_vars, frame=frame)
                custom_source_func = ctypes.cast(funcPtr, ctypes.c_void_p)
            else:
                custom_source_func = ctypes.c_void_p(0)
            nvars = num_standard_names + ncustom
            outputs = np.empty((k.shape[0], times.shape[0], nvars))
            if CAMB_TimeEvolution(byref(self), byref(c_int(k.shape[0])), k, byref(c_int(times.shape[0])),
                                  times[indices], byref(c_int(nvars)), outputs,
                                  byref(c_int(ncustom)), byref(custom_source_func)):
                config.check_global_error('get_time_evolution')
            i_rev = np.zeros(times.shape, dtype=int)
            i_rev[indices] = np.arange(times.shape[0])
            outputs = outputs[:, i_rev, :]
        finally:
            self.Params.Accuracy.lAccuracyBoost = old_boost
        if np.isscalar(q):
            return outputs[0, :, :][:, ix]
        else:
            return outputs[:, :, ix]

    def get_redshift_evolution(self, q, z, vars=model.evolve_names, lAccuracyBoost=4):
        """
        Get the mode evolution as a function of redshift for some k values.

        :param q: wavenumber values to calculate (or array of k values)
        :param z: array of redshifts to output
        :param vars: list of variable names or camb.symbolic sympy expressions to output
        :param lAccuracyBoost: boost factor for ell accuracy (e.g. to get nice smooth curves for plotting)
        :return: nd array, A_{qti}, size(q) x size(times) x len(vars), or 2d array if q is scalar
        """
        return self.get_time_evolution(q, self.conformal_time(z), vars, lAccuracyBoost)

    def get_background_time_evolution(self, eta, vars=model.background_names, format='dict'):
        """
        Get the evolution of background variables a function of conformal time.
        For the moment a and H are rather perversely only available via :meth:`get_time_evolution`

        :param eta: array of requested conformal times to output
        :param vars: list of variable names to output
        :param format: 'dict' or 'array', for either dict of 1D arrays indexed by name, or 2D array
        :return: n_eta x len(vars) 2D numpy array of outputs or dict of 1D arrays
        """

        if isinstance(vars, str):
            vars = [vars]
        unknown = set(vars) - set(model.background_names)
        if unknown:
            raise CAMBError('Unknown names %s; valid names are %s' % (unknown, model.background_names))
        outputs = np.zeros((eta.shape[0], 9))
        CAMB_BackgroundThermalEvolution(byref(self), byref(c_int(eta.shape[0])), eta, outputs)
        indices = [model.background_names.index(var) for var in vars]
        if format == 'dict':
            res = {}
            for var, index in zip(vars, indices):
                res[var] = outputs[:, index]
            return res
        else:
            assert format == 'array', "format must be dict or array"
            return outputs[:, np.array(indices)]

    def get_background_redshift_evolution(self, z, vars=model.background_names, format='dict'):
        """
        Get the evolution of background variables a function of redshift.
        For the moment a and H are rather perversely only available via :meth:`get_time_evolution`

        :param z: array of requested redshifts to output
        :param vars: list of variable names to output
        :param format: 'dict' or 'array', for either dict of 1D arrays indexed by name, or 2D array
        :return: n_eta x len(vars) 2D numpy array of outputs or dict of 1D arrays
        """

        return self.get_background_time_evolution(self.conformal_time(z), vars, format)

    def get_background_densities(self, a, vars=model.density_names, format='dict'):
        r"""
        Get the individual densities as a function of scale factor. Returns :math:`8\pi G a^4 \rho_i` in Mpc units.
        :math:`\Omega_i` can be simply obtained by taking the ratio of the components to tot.

        :param a: scale factor or array of scale factors
        :param vars: list of variables to output (default all)
        :param format: 'dict' or 'array', for either dict of 1D arrays indexed by name, or 2D array
        :return: n_a x len(vars) 2D numpy array or dict of 1D arrays of :math:`8\pi G a^4 \rho_i` in Mpc units.
        """
        if isinstance(vars, str):
            vars = [vars]
        unknown = set(vars) - set(model.density_names)
        if unknown:
            raise CAMBError('Unknown names %s; valid names are %s' % (unknown, model.density_names))
        arr = np.atleast_1d(a)
        outputs = np.zeros((arr.shape[0], 8))
        self.f_GetBackgroundDensities(byref(c_int(arr.shape[0])), arr, outputs)
        indices = [model.density_names.index(var) for var in vars]
        if format == 'dict':
            res = {}
            for var, index in zip(vars, indices):
                res[var] = outputs[:, index]
            return res
        else:
            assert format == 'array', "format must be dict or array"
            return outputs[:, np.array(indices)]

    def get_dark_energy_rho_w(self, a):
        r"""
        Get dark energy density in units of the dark energy density today, and equation of state parameter
        :math:`w\equiv P/\rho`

        :param a: scalar factor or array of scale factors
        :return: rho, w arrays at redshifts :math:`1/a-1` [or scalars if :math:`a` is scalar]
        """
        if np.isscalar(a):
            scales = np.array([a])
        else:
            scales = np.ascontiguousarray(a)
        rho = np.zeros(scales.shape)
        w = np.zeros(scales.shape)
        self.f_DarkEnergyStressEnergy(scales, rho, w, byref(c_int(len(scales))))
        if np.isscalar(a):
            return rho[0], w[0]
        else:
            return rho, w

    def get_Omega(self, var, z=0):
        r"""
        Get density relative to critical density of variables var

        :param var: one of 'K', 'cdm', 'baryon', 'photon', 'neutrino' (massless), 'nu' (massive neutrinos), 'de'
        :param z: redshift
        :return:  :math:`\Omega_i(a)`
        """
        dic = self.get_background_densities(1. / (1 + z), ['tot', var])
        res = dic[var] / dic['tot']
        if np.isscalar(z):
            return res[0]
        else:
            return res

    def get_matter_transfer_data(self) -> MatterTransferData:
        """
        Get matter transfer function data and sigma8 for calculated results.

        :return: :class:`.MatterTransferData` instance holding output arrays (copies, not pointers)
        """
        if not self.Params.WantTransfer:
            raise CAMBError("must have Params.WantTransfer to get matter transfers and power")

        cdata = _MatterTransferData()
        CAMBdata_mattertransferdata(byref(self), byref(cdata))
        data = MatterTransferData()
        data.nq = cdata.num_q_trans
        from numpy import ctypeslib as nplib
        data.q = nplib.as_array(cdata.q_trans, shape=(data.nq,)).copy()
        if cdata.sigma_8_size:
            data.sigma_8 = nplib.as_array(cdata.sigma_8, shape=(cdata.sigma_8_size,)).copy()
        else:
            data.sigma_8 = None
        if cdata.sigma2_vdelta_8_size:
            data.sigma2_vdelta_8 = nplib.as_array(cdata.sigma2_vdelta_8, shape=(cdata.sigma2_vdelta_8_size,)).copy()
        else:
            data.sigma2_vdelta_8 = None
        data.transfer_data = fortran_array(cdata.TransferData, cdata.TransferData_size, dtype=np.float32)
        return data

    @staticmethod
    def _transfer_var(var1, var2):
        if var1 is None:
            var1 = config.transfer_power_var
        if var2 is None:
            var2 = config.transfer_power_var
        if isinstance(var1, str):
            var1 = model.transfer_names.index(var1) + 1
        if isinstance(var2, str):
            var2 = model.transfer_names.index(var2) + 1
        return c_int(var1), c_int(var2)

    def get_linear_matter_power_spectrum(self, var1=None, var2=None, hubble_units=True, k_hunit=True,
                                         have_power_spectra=True, params=None, nonlinear=False):
        r"""
        Calculates :math:`P_{xy}(k)`, where x, y are one of model.Transfer_cdm, model.Transfer_xx etc.
        The output k values are not regularly spaced, and not interpolated. They are either k or k/h depending on the
        value of k_hunit (default True gives k/h).

        For a description of outputs for different var1, var2 see :ref:`transfer-variables`.

        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param hubble_units: if true, output power spectrum in (Mpc/h) units, otherwise Mpc
        :param k_hunit: if true, matter power is a function of k/h, if false, just k (both :math:`{\rm Mpc}^{-1}` units)
        :param have_power_spectra: set to False if not already computed power spectra
        :param params: if have_power_spectra=False, optional :class:`~.model.CAMBparams` instance
                       to specify new parameters
        :param nonlinear: include non-linear correction from halo model
        :return: k/h or k, z, PK, where kz and z are arrays of k/h or k and z respectively,
                 and PK[i,j] is the value at z[i], k[j]/h or k[j]
        """
        if self.OnlyTransfers or params is not None or not have_power_spectra:
            self.calc_power_spectra(params)

        num_k = c_int(0)
        CAMBdata_mattertransferks(byref(self), byref(num_k), np.array([]))
        nk = num_k.value

        ks = np.empty(nk, dtype=np.float64)
        CAMBdata_mattertransferks(byref(self), byref(num_k), ks)

        nz = self.Params.Transfer.PK_num_redshifts
        if k_hunit:
            kh = ks / (self.Params.H0 / 100)
        else:
            kh = ks
        var1, var2 = self._transfer_var(var1, var2)

        hubble_units = c_int(hubble_units)
        PK = np.empty((nz, nk))
        if nonlinear:
            CAMBdata_GetNonLinearMatterPower(byref(self), PK, byref(var1), byref(var2), byref(hubble_units))
            config.check_global_error('get_[non]linear_matter_power_spectrum')
        else:
            CAMBdata_GetLinearMatterPower(byref(self), PK, byref(var1), byref(var2), byref(hubble_units))

        z = self.Params.Transfer.PK_redshifts[:nz]
        z.reverse()
        return kh, np.array(z), PK

    def get_nonlinear_matter_power_spectrum(self, var1=None, var2=None, hubble_units=True, k_hunit=True,
                                            have_power_spectra=True, params=None):
        r"""
        Calculates :math:`P_{xy}(k/h)`, where x, y are one of model.Transfer_cdm, model.Transfer_xx etc.
        The output k values are not regularly spaced, and not interpolated.

        For a description of outputs for different var1, var2 see :ref:`transfer-variables`.

        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param hubble_units: if true, output power spectrum in :Math:`({\rm Mpc}/h)^{3}` units,
                             otherwise :math:`{\rm Mpc}^{3}`
        :param k_hunit: if true, matter power is a function of k/h, if false, just k (both :math:`{\rm Mpc}^{-1}` units)
        :param have_power_spectra: set to False if not already computed power spectra
        :param params: if have_power_spectra=False, optional :class:`~.model.CAMBparams` instance
                       to specify new parameters
        :return: k/h or k, z, PK, where kz and z are arrays of k/h or k and z respectively,
                 and PK[i,j] is the value at z[i], k[j]/h or k[j]
        """
        return self.get_linear_matter_power_spectrum(var1=var1, var2=var2, hubble_units=hubble_units, k_hunit=k_hunit,
                                                     have_power_spectra=have_power_spectra, params=params,
                                                     nonlinear=True)

    def get_sigmaR(self, R, z_indices=None, var1=None, var2=None, hubble_units=True, return_R_z=False):
        r"""
        Calculate :math:`\sigma_R` values, the RMS linear matter fluctuation in spheres of radius R in linear theory.
        Accuracy depends on the sampling with which the matter transfer functions are computed.

        For a description of outputs for different var1, var2 see :ref:`transfer-variables`.
        Note that numerical errors are slightly different to get_sigma8 for R=8 Mpc/h.

        :param R: radius in Mpc or h^{-1} Mpc units, scalar or array
        :param z_indices: indices of redshifts in Params.Transfer.PK_redshifts to calculate
                          (default None gives all computed in order of increasing time (decreasing redshift);
                          -1 gives redshift 0; list gives all listed indices)
        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param hubble_units: if true, R is in h^{-1} Mpc, otherwise Mpc
        :param return_R_z: if true, return tuple of R, z, sigmaR
                           (where R always Mpc units not h^{-1}Mpc and R, z are arrays)
        :return: array of :math:`\sigma_R` values, or 2D array indexed by (redshift, R)
        """
        assert self.Params.WantTransfer
        var1, var2 = self._transfer_var(var1, var2)
        if hubble_units:
            if not np.isscalar(R):
                R = np.array(R, dtype=np.float64)
            R = R / (self.Params.H0 / 100)

        radii = np.atleast_1d(np.ascontiguousarray(R, dtype=np.float64)) * self.Params.H0 / 100
        valid_indices = np.array(self.PK_redshifts_index[:self.Params.Transfer.PK_num_redshifts], dtype=np.int32)
        if z_indices is None:
            z_ix = valid_indices
        else:
            z_ix = np.atleast_1d(valid_indices[z_indices])
        sigma_R = np.empty((len(z_ix), len(radii)), dtype=np.float64)
        CAMBdata_GetSigmaRArray(byref(self), sigma_R, radii, byref(c_int(len(radii))), z_ix,
                                byref(c_int(len(z_ix))), byref(var1), byref(var2))
        if z_indices is not None and np.isscalar(z_indices):
            sigma_R = sigma_R[0, :]
            if np.isscalar(R):
                sigma_R = sigma_R[0]
        elif np.isscalar(R):
            sigma_R = sigma_R[:, 0]
        if return_R_z:
            return R, np.array(self.transfer_redshifts)[z_ix - 1], sigma_R
        else:
            return sigma_R

    def get_sigma8(self):
        r"""
        Get :math:`\sigma_8` values at Params.PK_redshifts (must previously have calculated power spectra)

        :return: array of :math:`\sigma_8` values, in order of increasing time (decreasing redshift)
        """
        assert self.Params.WantTransfer
        sigma8 = np.empty(self.Params.Transfer.PK_num_redshifts, dtype=np.float64)
        CAMBdata_GetSigma8(byref(self), sigma8, byref(c_int(0)))
        return sigma8

    def get_sigma8_0(self):
        r"""
        Get :math:`\sigma_8` value today (must previously have calculated power spectra)

        :return: :math:`\sigma_8` today
        """
        nz = self.Params.Transfer.PK_num_redshifts
        if not nz or self.Params.Transfer.PK_redshifts[nz - 1] > 1e-5:
            raise CAMBError("sigma8 requested at z=0, but P(z=0) not calcaulted")
        return self.get_sigma8()[-1]

    def get_fsigma8(self):
        r"""
        Get :math:`f\sigma_8` growth values (must previously have calculated power spectra).
        For general models :math:`f\sigma_8` is defined as in the Planck 2015 parameter paper in terms of
        the velocity-density correlation: :math:`\sigma^2_{vd}/\sigma_{dd}` for :math:`8 h^{-1} {\rm Mpc}` spheres.

        :return: array of f*sigma_8 values, in order of increasing time (decreasing redshift)
        """
        assert self.Params.WantTransfer
        fsigma8 = np.empty(self.Params.Transfer.PK_num_redshifts, dtype=np.float64)
        CAMBdata_GetSigma8(byref(self), fsigma8, byref(c_int(1)))
        return fsigma8

    def get_matter_power_spectrum(self, minkh=1e-4, maxkh=1.0, npoints=100,
                                  var1=None, var2=None,
                                  have_power_spectra=False, params=None):
        """
        Calculates :math:`P_{xy}(k/h)`, where x, y are one of Transfer_cdm, Transfer_xx etc.
        The output k values are regularly log spaced and interpolated. If NonLinear is set, the result is non-linear.

        For a description of outputs for different var1, var2 see :ref:`transfer-variables`.

        :param minkh: minimum value of k/h for output grid (very low values < 1e-4 may not be calculated)
        :param maxkh: maximum value of k/h (check consistent with input params.Transfer.kmax)
        :param npoints: number of points equally spaced in log k
        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param have_power_spectra: set to True if already computed power spectra
        :param params: if have_power_spectra=False and want to specify new parameters,
                       a :class:`~.model.CAMBparams` instance
        :return: kh, z, PK, where kz and z are arrays of k/h and z respectively, and PK[i,j] is value at z[i], k/h[j]
        """

        if not have_power_spectra:
            self.calc_power_spectra(params)

        if not npoints >= 2:
            raise CAMBError('Need at least two points in get_matter_power_spectrum')

        assert self.Params.WantTransfer
        if self.Params.Transfer.kmax < maxkh * self.Params.h:
            logging.warning("get_matter_power_spectrum using larger k_max than input parameter Transfer.kmax")
        if self.Params.NonLinear != model.NonLinear_none and self.Params.Transfer.kmax < 1:
            logging.warning("get_matter_power_spectrum Transfer.kmax small to get non-linear spectrum")

        nz = self.Params.Transfer.PK_num_redshifts
        PK = np.empty((nz, npoints))
        var1, var2 = self._transfer_var(var1, var2)

        dlnkh = (np.log(maxkh) - np.log(minkh)) / (npoints - 1)
        CAMBdata_GetMatterPower(byref(self), PK, byref(c_double(minkh)),
                                byref(c_double(dlnkh)), byref(c_int(npoints)), byref(var1), byref(var2))
        z = self.Params.Transfer.PK_redshifts[:nz]
        z.reverse()
        return minkh * np.exp(np.arange(npoints) * dlnkh), z, PK

    def get_matter_power_interpolator(self, nonlinear=True, var1=None, var2=None, hubble_units=True, k_hunit=True,
                                      return_z_k=False, log_interp=True, extrap_kmax=None, silent=False):
        r"""
        Assuming transfers have been calculated, return a 2D spline interpolation object to evaluate matter
        power spectrum as function of z and k/h (or k). Uses self.Params.Transfer.PK_redshifts as the spline node
        points in z. If fewer than four redshift points are used the interpolator uses a reduced order spline in z
        (so results at intermediate z may be innaccurate), otherwise it uses bicubic.
        Usage example:

        .. code-block:: python

           PK = results.get_matter_power_interpolator();
           print('Power spectrum at z=0.5, k/h=0.1 is %s (Mpc/h)^3 '%(PK.P(0.5, 0.1)))

        For a description of outputs for different var1, var2 see :ref:`transfer-variables`.

        :param nonlinear: include non-linear correction from halo model
        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param hubble_units: if true, output power spectrum in :math:`({\rm Mpc}/h)^{3}` units,
                             otherwise :math:`{\rm Mpc}^{3}`
        :param k_hunit: if true, matter power is a function of k/h, if false, just k (both :math:`{\rm Mpc}^{-1}` units)
        :param return_z_k: if true, return interpolator, z, k where z, k are the grid used
        :param log_interp: if true, interpolate log of power spectrum
                           (unless any values cross zero in which case ignored)
        :param extrap_kmax: if set, use power law extrapolation beyond kmax to extrap_kmax
                            (useful for tails of integrals)
        :param silent: Set True to silence warnings
        :return: An object PK based on :class:`~scipy:scipy.interpolate.RectBivariateSpline`,
                 that can be called with PK.P(z,kh) or PK(z,log(kh)) to get log matter power values.
                 If return_z_k=True, instead return interpolator, z, k where z, k are the grid used.
        """

        class PKInterpolator(RectBivariateSpline):
            islog: bool
            logsign: int

            def P(self, z, kh, grid=None):
                if grid is None:
                    grid = not np.isscalar(z) and not np.isscalar(kh)
                if self.islog:
                    return self.logsign * np.exp(self(z, np.log(kh), grid=grid))
                else:
                    return self(z, np.log(kh), grid=grid)

        class PKInterpolatorSingleZ(interp1d):
            islog: bool
            logsign: int

            def __init__(self, *args, **kwargs):
                self._single_z = np.array(args[0])
                super().__init__(*(args[1:]), kind=kwargs.get("ky"))

            def check_z(self, z):
                if not np.allclose(z, self._single_z):
                    raise CAMBError(
                        "P(z,k) requested at z=%g, but only computed for z=%s. "
                        "Cannot extrapolate!" % (z, self._single_z))

            def __call__(self, *args):
                self.check_z(args[0])
                # NB returns dimensionality as the 2D one: 1 dimension if z single
                return (lambda x: x[0] if np.isscalar(args[0]) else x)(super().__call__(*(args[1:])))

            def P(self, z, kh, **_kwargs):
                # grid kwarg is ignored
                if self.islog:
                    return self.logsign * np.exp(self(z, np.log(kh)))
                else:
                    return self(z, np.log(kh))

        assert self.Params.WantTransfer
        khs, zs, pk = self.get_linear_matter_power_spectrum(var1, var2, hubble_units, nonlinear=nonlinear)
        kh_max = khs[-1]
        if not k_hunit:
            khs *= self.Params.H0 / 100
        sign = 1
        if log_interp and np.any(pk <= 0):
            if np.all(pk < 0):
                sign = -1
            else:
                log_interp = False
        p_or_log_p = np.log(sign * pk) if log_interp else pk
        logkh = np.log(khs)
        deg_z = min(len(zs) - 1, 3)
        kmax = khs[-1]
        PKInterpolator = PKInterpolator if deg_z else PKInterpolatorSingleZ
        if extrap_kmax and extrap_kmax > kmax:
            # extrapolate to ultimate power law
            # TODO: use more physical extrapolation function for linear case
            if not silent and (kh_max < 3 and extrap_kmax > 2 and nonlinear or kh_max < 0.4):
                logging.warning("Extrapolating to higher k with matter transfer functions "
                                "only to k=%.3g Mpc^{-1} may be inaccurate.\n " % (kh_max * self.Params.H0 / 100))
            if not log_interp:
                raise CAMBValueError(
                    "Cannot use extrap_kmax with log_inter=False (e.g. PK crossing zero for %s, %s.)" % (var1, var2))

            logextrap = np.log(extrap_kmax)
            log_p_new = np.empty((pk.shape[0], pk.shape[1] + 2))
            log_p_new[:, :-2] = p_or_log_p
            delta = logextrap - logkh[-1]

            dlog = (log_p_new[:, -3] - log_p_new[:, -4]) / (logkh[-1] - logkh[-2])
            log_p_new[:, -1] = log_p_new[:, -3] + dlog * delta
            log_p_new[:, -2] = log_p_new[:, -3] + dlog * delta * 0.9
            logkh = np.hstack((logkh, logextrap - delta * 0.1, logextrap))
            p_or_log_p = log_p_new

        deg_k = min(len(logkh) - 1, 3)
        res = PKInterpolator(zs, logkh, p_or_log_p, kx=deg_z, ky=deg_k)
        res.kmin = np.min(khs)
        res.kmax = kmax
        res.islog = log_interp
        res.logsign = sign
        res.zmin = np.min(zs)
        res.zmax = np.max(zs)
        if return_z_k:
            return res, zs, khs
        else:
            return res

    def get_total_cls(self, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get lensed-scalar + tensor CMB power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE
        """
        lmax = self._lmax_setting(lmax)
        res = np.empty((lmax + 1, 4))
        opt = c_int(lmax)
        CAMB_SetTotCls(byref(self), byref(opt), res)
        self._scale_cls(res, CMB_unit, raw_cl)
        return res

    def get_tensor_cls(self, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get tensor CMB power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE
        """

        if lmax is None:
            lmax = self.Params.max_l_tensor
        lmax = self._lmax_setting(lmax, unlensed=True)
        res = np.empty((lmax + 1, 4))
        opt = c_int(lmax)
        CAMB_SetTensorCls(byref(self), byref(opt), res)
        self._scale_cls(res, CMB_unit, raw_cl)
        return res

    def get_unlensed_scalar_cls(self, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get unlensed scalar CMB power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE. CL[:,2] will be zero.
        """
        lmax = self._lmax_setting(lmax, unlensed=True)
        res = np.empty((lmax + 1, 4))
        opt = c_int(lmax)
        CAMB_SetUnlensedScalCls(byref(self), byref(opt), res)
        self._scale_cls(res, CMB_unit, raw_cl)
        return res

    def get_unlensed_total_cls(self, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get unlensed CMB power spectra, including tensors if relevant. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE.
        """
        lmax = self._lmax_setting(lmax, unlensed=True)
        return self.get_unlensed_scalar_cls(lmax, CMB_unit, raw_cl) + self.get_tensor_cls(lmax, CMB_unit, raw_cl)

    def get_lensed_scalar_cls(self, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get lensed scalar CMB power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE.
        """

        lmax = self._lmax_setting(lmax)
        res = np.empty((lmax + 1, 4))
        opt = c_int(lmax)
        CAMB_SetLensedScalCls(byref(self), byref(opt), res)
        self._scale_cls(res, CMB_unit, raw_cl)
        return res

    def get_lens_potential_cls(self, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get lensing deflection angle potential power spectrum, and cross-correlation with T and E. Must have already
        calculated power spectra.
        Power spectra are :math:`[L(L+1)]^2C_L^{\phi\phi}/2\pi` and corresponding deflection cross-correlations.

        :param lmax: lmax to output to

        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K` units for lensing cross.
        :param raw_cl: return lensing potential :math:`C_L` rather than :math:`[L(L+1)]^2C_L/2\pi`
        :return: numpy array CL[0:lmax+1,0:3], where 0..2 indexes PP, PT, PE.
        """

        lmax = self._lmax_setting(lmax, unlensed=True)
        res = np.empty((lmax + 1, 3))
        opt = c_int(lmax)
        CAMB_SetLensPotentialCls(byref(self), byref(opt), res)
        self._scale_cls(res, CMB_unit, raw_cl, lens_potential=True)
        return res

    def get_unlensed_scalar_array_cls(self, lmax=None):
        r"""
        Get array of all cross power spectra. Must have already calculated power spectra.
        Results are dimensionless, and not scaled by custom_scaled_ell_fac.

        :param lmax: lmax to output to
        :return: numpy array CL[0:, 0:,0:lmax+1], where 0.. index T, E, lensing potential, source window functions
        """

        lmax = self._lmax_setting(lmax, unlensed=True)
        if not self.Params.Want_cl_2D_array:
            raise CAMBError('unlensed_scalar_array not calculated (set Want_cl_2D_array)')

        n = 3 + len(self.Params.SourceWindows) + self.Params.CustomSources.num_custom_sources
        res = np.empty((n, n, lmax + 1), order='F')
        CAMB_SetUnlensedScalarArray(byref(self), byref(c_int(lmax)), res, byref(c_int(n)))
        return res

    def get_cmb_unlensed_scalar_array_dict(self, params=None, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get all unlensed auto and cross power spectra, including any custom source functions set
        using :meth:`.model.CAMBparams.set_custom_scalar_sources`.

        :param params: optional :class:`~.model.CAMBparams` instance with parameters to use. If None, must have
          previously set parameters and called :meth:`calc_power_spectra` (e.g. if you got this instance
          using :func:`.camb.get_results`),
        :param lmax: maximum :math:`\ell`
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for
          CMB :math:`C_\ell` and :math:`\mu K` units for lensing cross.
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: dictionary of power spectrum arrays, index as TxT, TxE, PxW1, W1xW2, custom_name_1xT... etc.
                 Note that P is the lensing deflection, lensing windows Wx give convergence.
        """

        old_val = None
        try:
            if params is not None:
                old_val = params.Want_cl_2D_array
                params.Want_cl_2D_array = True
                self.calc_power_spectra(params)
            elif not self.Params.Want_cl_2D_array:
                raise CAMBValueError('Want_cl_2D_array must be true to have array C_L')
            nwindows = len(self.Params.SourceWindows)
            lmax = lmax or self.Params.max_l
            arr = self.get_unlensed_scalar_array_cls(lmax)
            custom_source_names = self.Params.get_custom_source_names()
            names = ['T', 'E', 'P'] + ["W%s" % (i + 1) for
                                       i in range(nwindows)] + custom_source_names
            CMB_unit = self._CMB_unit(CMB_unit) or 1
            CMB_units = [CMB_unit, CMB_unit, 1] + [1] * nwindows + [CMB_unit] * len(custom_source_names)

            result = {}
            for i, name in enumerate(names):
                for j, name2 in enumerate(names):
                    tag = name + 'x' + name2
                    if j < i:
                        result[tag] = result[name2 + 'x' + name]
                    else:
                        cls = arr[i, j, :]
                        if raw_cl:
                            ls = np.arange(1, cls.shape[0])
                            fac = np.float64(ls * (ls + 1))
                            if i == 2 and j == 2:
                                fac *= fac
                            elif i == 2 or j == 2:
                                fac *= np.sqrt(fac)
                            cls[1:] /= (fac / (2 * np.pi))

                        if CMB_unit is not None:
                            cls *= CMB_units[i] * CMB_units[j]
                        result[tag] = cls
        finally:
            if old_val is not None:
                params.Want_cl_2D_array = old_val
        return result

    def get_source_cls_dict(self, params=None, lmax=None, raw_cl=False):
        r"""
        Get all source window function and CMB lensing and cross power spectra. Does not include CMB spectra.
        Note that P is the deflection angle, but lensing windows return the kappa power.

        :param params: optional :class:`~.model.CAMBparams` instance with parameters to use. If None, must have
          previously set parameters and called :meth:`calc_power_spectra`
          (e.g. if you got this instance using :func:`.camb.get_results`),
        :param lmax: maximum :math:`\ell`
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: dictionary of power spectrum arrays, index as PXP, PxW1, W1xW2, ... etc.
        """

        try:
            if params is not None:
                old_val = params.Want_cl_2D_array
                params.Want_cl_2D_array = True
                self.calc_power_spectra(params)
            elif not self.Params.Want_cl_2D_array:
                raise CAMBValueError('Want_cl_2D_array must be true to have array C_L')
            nwindows = len(self.Params.SourceWindows)
            lmax = lmax or self.Params.max_l
            arr = self.get_unlensed_scalar_array_cls(lmax)
            names = ['P'] + ["W%s" % (i + 1) for i in range(nwindows)]

            result = {}
            for i, name in enumerate(names):
                for j, name2 in enumerate(names):
                    tag = name + 'x' + name2
                    if j < i:
                        result[tag] = result[name2 + 'x' + name]
                    else:
                        cls = arr[i + 2, j + 2, :]
                        if raw_cl:
                            ls = np.arange(1, cls.shape[0])
                            fac = np.float64(ls * (ls + 1))
                            if i == 0 and j == 0:
                                fac *= fac
                            elif i == 0 or j == 0:
                                fac *= np.sqrt(fac)
                            cls[1:] /= (fac / (2 * np.pi))
                        result[tag] = cls
        finally:
            if params is not None:
                params.Want_cl_2D_array = old_val
        return result

    def get_lensed_gradient_cls(self, lmax=None, CMB_unit=None, raw_cl=False, clpp=None):
        r"""
        Get lensed gradient scalar CMB power spectra in flat sky approximation
        (`arXiv:1101.2234 <https://arxiv.org/abs/1101.2234>`_).
        Note that lmax used to calculate results may need to be substantially larger than the lmax
        output from this function (there is no extrapolation as in the main lensing routines).
        Lensed power spectra must be already calculated.

        :param lmax: lmax to output to
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :param clpp: custom array of :math:`[L(L+1)]^2 C_L^{\phi\phi}/2\pi` lensing potential power spectrum
             to use (zero based), rather than calculated specturm from this model
        :return: numpy array CL[0:lmax+1,0:8], where CL[:,i] are :math:`T\nabla T`, :math:`E\nabla E`,
                 :math:`B\nabla B`, :math:`PP_\perp`, :math:`T\nabla E`, :math:`TP_\perp`, :math:`(\nabla T)^2`,
                 :math:`\nabla T\nabla T` where the first six are as defined in appendix C
                 of `1101.2234 <https://arxiv.org/abs/1101.2234>`_.
        """
        assert self.Params.DoLensing
        lmax = self._lmax_setting(lmax)
        res = np.empty((lmax + 1, 8))
        opt = c_int(lmax)
        if clpp is not None:
            if clpp.shape[0] < self.Params.max_l + 1:
                raise CAMBValueError('clpp must go to at least Params.max_l (zero based)')
            clpp = np.array(clpp, dtype=np.float64)
            GetFlatSkyCgrads = lib_import('lensing', '', 'getflatskycgradswithspectrum')
            GetFlatSkyCgrads.argtypes = [POINTER(CAMBdata), numpy_1d, int_arg, numpy_1d]
            GetFlatSkyCgrads(byref(self), clpp, byref(opt), res)
        else:
            GetFlatSkyCgrads = lib_import('lensing', '', 'getflatskycgrads')
            GetFlatSkyCgrads.argtypes = [POINTER(CAMBdata), int_arg, numpy_1d]
            GetFlatSkyCgrads(byref(self), byref(opt), res)
        self._scale_cls(res, CMB_unit, raw_cl)
        return res

    def get_lensed_cls_with_spectrum(self, clpp, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get lensed CMB power spectra using curved-sky correlation function method, using
        cpp as the lensing spectrum. Useful for e.g. getting partially-delensed spectra.

        :param clpp: array of :math:`[L(L+1)]^2 C_L^{\phi\phi}/2\pi` lensing potential power spectrum (zero based)
        :param lmax: lmax to output to
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE.
        """
        assert self.Params.DoLensing
        lmax_unlens = self.Params.max_l
        if clpp.shape[0] < lmax_unlens + 1:
            raise CAMBValueError('clpp must go to at least Params.max_l (zero based)')
        res = np.zeros((lmax_unlens + 1, 4), dtype=np.float64)
        lmax_lensed = c_int(0)
        lensClsWithSpectrum = lib_import('lensing', '', 'lensclswithspectrum')
        lensClsWithSpectrum.argtypes = [POINTER(CAMBdata), numpy_1d, numpy_2d, int_arg]
        clpp = np.array(clpp, dtype=np.float64)
        lensClsWithSpectrum(byref(self), clpp, res, byref(lmax_lensed))
        res = res[:min(lmax_lensed.value, lmax or lmax_lensed.value) + 1, :]
        self._scale_cls(res, CMB_unit, raw_cl)
        return res

    def get_partially_lensed_cls(self, Alens, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
           Get lensed CMB power spectra using curved-sky correlation function method, using
           true lensing spectrum scaled by Alens. Alens can be an array in L for realistic delensing estimates.
           Note that if Params.Alens is also set, the result is scaled by the product of both

           :param Alens: scaling of the lensing relative to true, with Alens=1 being the standard result.
              Can be a scalar in which case all L are scaled, or a zero-based array with the L by L scaling
              (with L larger than the size of the array having Alens_L=1).
           :param lmax: lmax to output to
           :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell`
           :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
           :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE.
           """
        clpp = self.get_lens_potential_cls()[:, 0]
        if np.isscalar(Alens):
            clpp *= Alens
        else:
            clpp[:Alens.size] *= Alens
        return self.get_lensed_cls_with_spectrum(clpp, lmax, CMB_unit, raw_cl)

    def angular_diameter_distance(self, z):
        """
        Get (non-comoving) angular diameter distance to redshift z.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer
        functions or power spectra.

        :param z: redshift or array of redshifts
        :return: angular diameter distances, matching rank of z
        """
        if np.isscalar(z):
            return self.f_AngularDiameterDistance(byref(c_double(z)))
        else:
            z = np.asarray(z)
            arr = np.empty(z.shape)
            indices = np.argsort(z)
            redshifts = np.array(z[indices], dtype=np.float64)
            self.f_AngularDiameterDistanceArr(arr, redshifts, byref(c_int(z.shape[0])))
            arr[indices] = arr.copy()
            return arr

    def _make_scalar_or_arrays(self, z1, z2):
        if np.isscalar(z1):
            if np.isscalar(z2):
                return z1, z2
            else:
                z1 = np.ones(len(z2)) * z1
        else:
            z1 = np.ascontiguousarray(z1, dtype=np.float64)
            if np.isscalar(z2):
                z2 = np.ones(len(z1)) * z2
            else:
                z2 = np.ascontiguousarray(z2, dtype=np.float64)
                if len(z1) != len(z2):
                    raise CAMBError('z1 nand z2 must be scalar or same-length 1D arrays')
        return z1, z2

    def angular_diameter_distance2(self, z1, z2):
        r"""
        Get angular diameter distance between two redshifts
        :math:`\frac{r}{1+z_2}\text{sin}_K\left(\frac{\chi(z_2) - \chi(z_1)}{r}\right)`
        where :math:`r` is curvature radius and :math:`\chi` is the comoving radial distance.
        If z_1 >= z_2 returns zero.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer
        functions or power spectra.

        :param z1: redshift 1, or orray of redshifts
        :param z2: redshift 2, or orray of redshifts
        :return: result (scalar or array of distances between pairs of z1, z2)
        """
        z1, z2 = self._make_scalar_or_arrays(z1, z2)
        if not np.isscalar(z1):
            dists = np.empty(z1.shape)
            self.f_AngularDiameterDistance2Arr(dists, z1, z2, byref(c_int(dists.shape[0])))
            return dists
        else:
            return self.f_AngularDiameterDistance2(byref(c_double(z1)), byref(c_double(z2)))

    def comoving_radial_distance(self, z, tol=1e-4):
        """
        Get comoving radial distance from us to redshift z in Mpc. This is efficient for arrays.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer
        functions or power spectra.

        :param z: redshift
        :param tol: numerical tolerance parameter
        :return: comoving radial distance (Mpc)
        """
        if not np.isscalar(z):
            indices = np.argsort(z)
            redshifts = np.array(z[indices], dtype=np.float64)
            chis = np.empty(redshifts.shape)
            self.f_ComovingRadialDistanceArr(chis, redshifts, byref(c_int(chis.shape[0])), byref(c_double(tol)))
            chis[indices] = chis.copy()
            return chis
        else:
            return self.f_ComovingRadialDistance(byref(c_double(z)))

    def redshift_at_comoving_radial_distance(self, chi):
        """
        Convert comoving radial distance array to redshift array.

        :param chi: comoving radial distance (in Mpc), scalar or array
        :return: redshift at chi, scalar or array
        """
        return self.redshift_at_conformal_time(self.tau0 - chi)

    def redshift_at_conformal_time(self, eta):
        """
        Convert conformal time array to redshift array.
        Note that this function requires the transfers or background to have been calculated with no_thermo=False
        (the default).

        :param eta: conformal time from bing bang (in Mpc), scalar or array
        :return: redshift at eta, scalar or array
        """

        if np.isscalar(eta):
            times = np.array([eta], dtype=np.float64)
        else:
            times = np.ascontiguousarray(eta, dtype=np.float64)
        redshifts = np.empty(times.shape)
        self.f_RedshiftAtTimeArr(redshifts, times, byref(c_int(times.shape[0])))
        config.check_global_error('redshift_at_conformal_time')
        if np.isscalar(eta):
            return redshifts[0]
        else:
            return redshifts

    def luminosity_distance(self, z):
        """
        Get luminosity distance from to redshift z.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions or
        power spectra.

        :param z: redshift or array of redshifts
        :return: luminosity distance (matches rank of z)
        """

        if not np.isscalar(z):
            z = np.asarray(z)
        return self.angular_diameter_distance(z) * (1.0 + z) ** 2

    def h_of_z(self, z):
        r"""
        Get Hubble rate at redshift z, in :math:`{\rm Mpc}^{-1}` units, scalar or array

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions
        or power spectra.

        Use hubble_parameter instead if you want in [km/s/Mpc] units.

        :param z: redshift
        :return: H(z)
        """
        if not np.isscalar(z):
            z = np.array(z, dtype=np.float64)
            arr = np.empty(z.shape)
            self.f_HofzArr(arr, z, byref(c_int(z.shape[0])))
            return arr
        else:
            return self.f_Hofz(byref(c_double(z)))

    def hubble_parameter(self, z):
        """
        Get Hubble rate at redshift z, in km/s/Mpc units. Scalar or array.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions
        or power spectra.

        :param z: redshift
        :return: H(z)/[km/s/Mpc]
        """
        return (constants.c / 1e3) * self.h_of_z(z)

    def physical_time_a1_a2(self, a1, a2):
        """
        Get physical time between two scalar factors in Julian Gigayears

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions
        or power spectra.

        :param a1: scale factor 1
        :param a2: scale factor 2
        :return: (age(a2)-age(a1))/Gigayear
        """
        a1, a2 = self._make_scalar_or_arrays(a1, a2)
        if not np.isscalar(a1):
            times = np.empty(a1.shape)
            self.f_DeltaPhysicalTimeGyrArr(times, a1, a2, byref(c_int(times.shape[0])), None)
            return times
        else:
            return self.f_DeltaPhysicalTimeGyr(byref(c_double(a1)), byref(c_double(a2)), None)

    def physical_time(self, z):
        """
        Get physical time from hot big bang to redshift z in Julian Gigayears.

        :param z:  redshift
        :return: t(z)/Gigayear
        """
        if not np.isscalar(z):
            z = np.asarray(z, dtype=np.float64)
        return self.physical_time_a1_a2(0, 1.0 / (1 + z))

    def conformal_time_a1_a2(self, a1, a2):
        """
        Get conformal time between two scale factors (=comoving radial distance travelled by light on light cone)

        :param a1: scale factor 1
        :param a2: scale factor 2
        :return: eta(a2)-eta(a1) = chi(a1)-chi(a2) in Megaparsec
        """
        a1, a2 = self._make_scalar_or_arrays(a1, a2)
        if not np.isscalar(a1):
            times = np.empty(a1.shape)
            self.f_DeltaTimeArr(times, a1, a2, byref(c_int(times.shape[0])), None)
            return times
        else:
            return self.f_DeltaTime(byref(c_double(a1)), byref(c_double(a2)), None)

    def conformal_time(self, z, presorted=None, tol=None):
        """
        Conformal time from hot big bang to redshift z in Megaparsec.

        :param z: redshift or array of redshifts
        :param presorted: if True, redshifts already sorted to be monotonically increasing, if False decreasing,
           or if None unsorted.  If presorted is True or False no checks are done.
        :param tol: integration tolerance
        :return: eta(z)/Mpc
        """
        if np.isscalar(z):
            redshifts = np.array([z], dtype=np.float64)
        else:
            redshifts = np.array(z, dtype=np.float64)
            if presorted is True:
                redshifts = redshifts[::-1].copy()
            elif presorted is None:
                indices = np.argsort(redshifts)[::-1]
                redshifts = redshifts[indices]

        eta = np.empty(redshifts.shape)
        if tol:
            tol = byref(c_double(tol))

        self.f_TimeOfzArr(eta, redshifts, byref(c_int(eta.shape[0])), tol)
        if np.isscalar(z):
            return eta[0]
        else:
            if presorted is False:
                return eta
            elif presorted is True:
                return eta[::-1]
            else:
                eta[indices] = eta.copy()
                return eta

    def sound_horizon(self, z):
        """
        Get comoving sound horizon as function of redshift in Megaparsecs, the integral of the sound speed
        up to given redshift.

        :param z: redshift or array of redshifts
        :return: r_s(z)
        """
        if np.isscalar(z):
            redshifts = np.array([z], dtype=np.float64)
        else:
            redshifts = np.array(z, dtype=np.float64)
        rs = np.empty(redshifts.shape)
        self.f_sound_horizon_zArr(rs, redshifts, byref(c_int(redshifts.shape[0])))
        if np.isscalar(z):
            return rs[0]
        else:
            return rs

    def cosmomc_theta(self):
        r"""
        Get :math:`\theta_{\rm MC}`, an approximation of the ratio of the sound horizon to the angular diameter
        distance at recombination.

        :return: :math:`\theta_{\rm MC}`
        """
        return self.f_CosmomcTheta()


CAMBdata_gettransfers = camblib.__handles_MOD_cambdata_gettransfers
CAMBdata_gettransfers.argtypes = [POINTER(CAMBdata), POINTER(model.CAMBparams),
                                  POINTER(c_int), POINTER(c_int)]
CAMBdata_gettransfers.restype = c_int

CAMBdata_transferstopowers = camblib.__camb_MOD_camb_transferstopowers
CAMBdata_transferstopowers.argtypes = [POINTER(CAMBdata)]

CAMBdata_mattertransferdata = camblib.__handles_MOD_cambdata_mattertransferdata
CAMBdata_mattertransferdata.argtypes = [POINTER(CAMBdata), POINTER(_MatterTransferData)]

CAMBdata_mattertransferks = camblib.__handles_MOD_cambdata_getmattertransferks
CAMBdata_mattertransferks.argtypes = [POINTER(CAMBdata), POINTER(c_int), numpy_1d]

CAMBdata_cltransferdata = camblib.__handles_MOD_cambdata_cltransferdata
CAMBdata_cltransferdata.argtypes = [POINTER(CAMBdata), POINTER(_ClTransferData), int_arg]
CAMBdata_GetLinearMatterPower = camblib.__handles_MOD_cambdata_getlinearmatterpower
CAMBdata_GetLinearMatterPower.argtypes = [POINTER(CAMBdata), numpy_2d, int_arg, int_arg, int_arg]

CAMBdata_GetNonLinearMatterPower = camblib.__handles_MOD_cambdata_getnonlinearmatterpower
CAMBdata_GetNonLinearMatterPower.argtypes = [POINTER(CAMBdata), numpy_2d, int_arg, int_arg, int_arg]

CAMBdata_GetMatterPower = camblib.__handles_MOD_cambdata_getmatterpower
CAMBdata_GetMatterPower.argtypes = [POINTER(CAMBdata), numpy_2d,
                                    d_arg, d_arg, int_arg, int_arg, int_arg]

CAMBdata_GetSigma8 = camblib.__handles_MOD_cambdata_getsigma8
CAMBdata_GetSigma8.argtypes = [POINTER(CAMBdata), numpy_1d, int_arg]

CAMBdata_GetSigmaRArray = camblib.__handles_MOD_cambdata_getsigmararray
CAMBdata_GetSigmaRArray.argtypes = [POINTER(CAMBdata), numpy_2d, numpy_1d, int_arg, numpy_1d_int, int_arg, int_arg,
                                    int_arg]

CAMBdata_CalcBackgroundTheory = camblib.__handles_MOD_cambdata_calcbackgroundtheory
CAMBdata_CalcBackgroundTheory.argtypes = [POINTER(CAMBdata), POINTER(model.CAMBparams)]
CAMBdata_CalcBackgroundTheory.restype = c_int

CAMB_SetTotCls = camblib.__handles_MOD_camb_settotcls
CAMB_SetUnlensedCls = camblib.__handles_MOD_camb_setunlensedcls
CAMB_SetLensPotentialCls = camblib.__handles_MOD_camb_setlenspotentialcls
CAMB_SetUnlensedScalCls = camblib.__handles_MOD_camb_setunlensedscalcls
CAMB_SetLensedScalCls = camblib.__handles_MOD_camb_setlensedscalcls
CAMB_SetTensorCls = camblib.__handles_MOD_camb_settensorcls

_set_cl_args = [POINTER(CAMBdata), int_arg, numpy_1d]

CAMB_SetTotCls.argtypes = _set_cl_args
CAMB_SetUnlensedCls.argtypes = _set_cl_args
CAMB_SetLensPotentialCls.argtypes = _set_cl_args
CAMB_SetUnlensedScalCls.argtypes = _set_cl_args
CAMB_SetTensorCls.argtypes = _set_cl_args
CAMB_SetLensedScalCls.argtypes = _set_cl_args

CAMB_SetUnlensedScalarArray = camblib.__handles_MOD_camb_setunlensedscalararray
CAMB_SetUnlensedScalarArray.argtypes = [POINTER(CAMBdata), int_arg, ndpointer(c_double, flags='F_CONTIGUOUS', ndim=3),
                                        int_arg]

del _set_cl_args

CAMB_TimeEvolution = camblib.__handles_MOD_camb_timeevolution
CAMB_TimeEvolution.restype = c_bool
CAMB_TimeEvolution.argtypes = [POINTER(CAMBdata), int_arg, numpy_1d, int_arg, numpy_1d,
                               int_arg, ndpointer(c_double, flags='C_CONTIGUOUS', ndim=3),
                               int_arg, POINTER(ctypes.c_void_p)]

CAMB_BackgroundThermalEvolution = camblib.__handles_MOD_getbackgroundthermalevolution
CAMB_BackgroundThermalEvolution.argtypes = [POINTER(CAMBdata), int_arg, numpy_1d, numpy_2d]

CAMB_GetBackgroundOutputs = camblib.__handles_MOD_camb_getbackgroundoutputs
CAMB_GetBackgroundOutputs.argtypes = [POINTER(CAMBdata), numpy_1d, int_arg]
