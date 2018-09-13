from .baseconfig import camblib, CAMBError, CAMBValueError, CAMBUnknownArgumentError, CAMB_Structure, dll_import, \
    mock_load
import ctypes
from ctypes import c_float, c_int, c_double, c_bool, POINTER, byref
from . import model, constants, initialpower, lensing
import numpy as np
from numpy import ctypeslib as nplib
from numpy.ctypeslib import ndpointer
import logging
import sys
import six
import copy
from inspect import ismethod

if not mock_load:
    from scipy.interpolate import UnivariateSpline, RectBivariateSpline

if six.PY3:
    from inspect import getfullargspec as getargspec
else:
    from inspect import getargspec

_debug_params = False


class _CAMBdata(CAMB_Structure):
    # contains complex types with pointers, so just set up dummy
    _fields_ = []


class _MatterTransferData(CAMB_Structure):
    # contains complex types with pointers, so just set up dummy
    _fields_ = [('num_q_trans', c_int),
                ('q_trans', POINTER(c_double)),
                ('sigma_8', POINTER(c_double)),
                ('sigma2_vdelta_8', POINTER(c_double)),
                ('TransferData', POINTER(c_float)),
                ('sigma_8_size', c_int * 2),
                ('sigma2_vdelta_8_size', c_int * 2),
                ('TransferData_size', c_int * 3)
                ]


class _ClTransferData(CAMB_Structure):
    _fields_ = [('NumSources', c_int),
                ('q_size', c_int),
                ('q', POINTER(c_double)),
                ('delta_size', c_int * 3),
                ('delta_p_l_k', POINTER(c_double)),
                ('l_size', c_int),
                ('l', POINTER(c_int))
                ]


# Use FeedbackLevel.value to read and set
FeedbackLevel = dll_import(c_int, "modelparams", "feedbacklevel")

model.has_cl_2D_array.value = True

int_arg = POINTER(c_int)
d_arg = POINTER(c_double)

# for the case where CAMB wrapper functions do the F-C conversion, so use C here
numpy_2d = ndpointer(c_double, flags='C_CONTIGUOUS', ndim=2)
numpy_1d = ndpointer(c_double, flags='C_CONTIGUOUS')

CAMBdata_new = camblib.__handles_MOD_cambdata_new
CAMBdata_new.argtypes = [POINTER(POINTER(_CAMBdata))]

CAMBdata_free = camblib.__handles_MOD_cambdata_free
CAMBdata_free.argtypes = [POINTER(POINTER(_CAMBdata))]

CAMBdata_getparams = camblib.__handles_MOD_cambdata_getparams
CAMBdata_getparams.argtypes = [POINTER(_CAMBdata), POINTER(POINTER(model.CAMBparams))]

CAMBdata_setparams = camblib.__handles_MOD_cambdata_setparams
CAMBdata_setparams.argtypes = [POINTER(_CAMBdata), POINTER(model.CAMBparams)]

CAMBdata_gettransfers = camblib.__handles_MOD_cambdata_gettransfers
CAMBdata_gettransfers.argtypes = [POINTER(_CAMBdata), POINTER(model.CAMBparams),
                                  POINTER(c_bool)]
CAMBdata_gettransfers.restype = c_int

CAMBdata_transferstopowers = camblib.__camb_MOD_camb_transferstopowers
CAMBdata_transferstopowers.argtypes = [POINTER(_CAMBdata)]

CAMBdata_mattertransferdata = camblib.__handles_MOD_cambdata_mattertransferdata
CAMBdata_mattertransferdata.argtypes = [POINTER(_CAMBdata), POINTER(_MatterTransferData)]

CAMBdata_cltransferdata = camblib.__handles_MOD_cambdata_cltransferdata
CAMBdata_cltransferdata.argtypes = [POINTER(_CAMBdata), POINTER(_ClTransferData), int_arg]
CAMB_SetTotCls = camblib.__handles_MOD_camb_settotcls
CAMB_SetUnlensedCls = camblib.__handles_MOD_camb_setunlensedcls
CAMB_SetLensPotentialCls = camblib.__handles_MOD_camb_setlenspotentialcls
CAMB_SetUnlensedScalCls = camblib.__handles_MOD_camb_setunlensedscalcls
CAMB_SetLensedScalCls = camblib.__handles_MOD_camb_setlensedscalcls
CAMB_SetTensorCls = camblib.__handles_MOD_camb_settensorcls

_set_cl_args = [int_arg, numpy_1d, int_arg]

CAMB_SetTotCls.argtypes = _set_cl_args
CAMB_SetUnlensedCls.argtypes = _set_cl_args
CAMB_SetLensPotentialCls.argtypes = _set_cl_args
CAMB_SetUnlensedScalCls.argtypes = _set_cl_args
CAMB_SetTensorCls.argtypes = _set_cl_args
CAMB_SetLensedScalCls.argtypes = _set_cl_args

CAMB_SetUnlensedScalarArray = camblib.__handles_MOD_camb_setunlensedscalararray
CAMB_SetUnlensedScalarArray.argtypes = [int_arg, ndpointer(c_double, flags='F_CONTIGUOUS', ndim=3),
                                        int_arg, int_arg]

del _set_cl_args

CAMB_SetBackgroundOutputs_z = camblib.__handles_MOD_camb_setbackgroundoutputs_z
CAMB_SetBackgroundOutputs_z.argtypes = [numpy_1d, int_arg]
CAMB_GetBackgroundOutputs = camblib.__handles_MOD_camb_getbackgroundoutputs
CAMB_GetBackgroundOutputs.argtypes = [numpy_1d, int_arg]
CAMB_GetNumBackgroundOutputs = camblib.__handles_MOD_camb_getnumbackgroundoutputs
CAMB_GetNumBackgroundOutputs.restype = c_int

CAMB_GetAge = camblib.__camb_MOD_camb_getage
CAMB_GetAge.restype = c_double
CAMB_GetAge.argtypes = [POINTER(model.CAMBparams)]

CAMB_GetZreFromTau = camblib.__camb_MOD_camb_getzrefromtau
CAMB_GetZreFromTau.restype = c_double
CAMB_GetZreFromTau.argtypes = [POINTER(model.CAMBparams), d_arg]

CAMB_SetParamsForBackground = camblib.__handles_MOD_cambdata_setparamsforbackground
CAMB_SetParamsForBackground.argtypes = [POINTER(_CAMBdata), POINTER(model.CAMBparams)]

CAMB_CalcBackgroundTheory = camblib.__handles_MOD_cambdata_calcbackgroundtheory
CAMB_CalcBackgroundTheory.argtypes = [POINTER(_CAMBdata), POINTER(model.CAMBparams)]
CAMB_CalcBackgroundTheory.restype = c_int

CAMBdata_GetLinearMatterPower = camblib.__handles_MOD_cambdata_getlinearmatterpower
CAMBdata_GetLinearMatterPower.argtypes = [POINTER(_CAMBdata), numpy_2d, int_arg, int_arg, int_arg]

CAMBdata_GetNonLinearMatterPower = camblib.__handles_MOD_cambdata_getnonlinearmatterpower
CAMBdata_GetNonLinearMatterPower.argtypes = [POINTER(_CAMBdata), numpy_2d, int_arg, int_arg, int_arg]

CAMBdata_GetMatterPower = camblib.__handles_MOD_cambdata_getmatterpower
CAMBdata_GetMatterPower.argtypes = [POINTER(_CAMBdata), numpy_2d,
                                    d_arg, d_arg, int_arg, int_arg, int_arg]

AngularDiameterDistance = camblib.__modelparams_MOD_angulardiameterdistance
AngularDiameterDistance.argtyes = [d_arg]
AngularDiameterDistance.restype = c_double

AngularDiameterDistanceArr = camblib.__modelparams_MOD_angulardiameterdistancearr
AngularDiameterDistanceArr.argtypes = [numpy_1d, numpy_1d, int_arg]

AngularDiameterDistance2 = camblib.__modelparams_MOD_angulardiameterdistance2
AngularDiameterDistance2.argtyes = [d_arg]
AngularDiameterDistance2.restype = c_double

ComovingRadialDistance = camblib.__modelparams_MOD_comovingradialdistance
ComovingRadialDistance.argtyes = [d_arg]
ComovingRadialDistance.restype = c_double

ComovingRadialDistanceArr = camblib.__modelparams_MOD_comovingradialdistancearr
ComovingRadialDistanceArr.argtypes = [numpy_1d, numpy_1d, int_arg, d_arg]

TimeOfzArr = camblib.__modelparams_MOD_timeofzarr
TimeOfzArr.argtypes = [int_arg, numpy_1d, numpy_1d]

Hofz = camblib.__modelparams_MOD_hofz
Hofz.argtyes = [d_arg]
Hofz.restype = c_double

HofzArr = camblib.__modelparams_MOD_hofzarr
HofzArr.argtypes = [numpy_1d, numpy_1d, int_arg]

DeltaPhysicalTimeGyr = camblib.__modelparams_MOD_deltaphysicaltimegyr
DeltaPhysicalTimeGyr.argtypes = [d_arg, d_arg, d_arg]
DeltaPhysicalTimeGyr.restype = c_double

DeltaTime = camblib.__modelparams_MOD_deltatime
DeltaTime.argtypes = [d_arg, d_arg, d_arg]
DeltaTime.restype = c_double

CosmomcTheta = camblib.__modelparams_MOD_cosmomctheta
CosmomcTheta.restype = c_double

CAMB_TimeEvolution = camblib.__handles_MOD_camb_timeevolution
CAMB_TimeEvolution.restype = c_bool
CAMB_TimeEvolution.argtypes = [int_arg, numpy_1d, int_arg, numpy_1d,
                               int_arg, ndpointer(c_double, flags='C_CONTIGUOUS', ndim=3),
                               int_arg, POINTER(ctypes.c_void_p)]

CAMB_SetCustomSourcesFunc = camblib.__handles_MOD_camb_setcustomsourcesfunc
CAMB_SetCustomSourcesFunc.argtypes = [int_arg, POINTER(ctypes.c_void_p), ndpointer(c_int, flags='C_CONTIGUOUS')]

CAMB_BackgroundEvolution = camblib.__thermodata_MOD_getbackgroundevolution
CAMB_BackgroundEvolution.argtypes = [int_arg, numpy_1d, numpy_2d]


class MatterTransferData(object):
    r"""
    MatterTransferData is the base class for storing matter power transfer function data for various q values.
    In a flat universe q=k, in a closed universe q is quantized.

    To get an instance of this data, call :meth:`camb.CAMBdata.get_matter_transfer_data`

    :ivar nq:  number of q modes calculated
    :ivar q: array of q values calculated
    :ivar sigma_8: array of :math:`\sigma_8` values for each redshift for each power spectrum
    :ivar sigma2_vdelta_8: array of v-delta8 correlation, so sigma2_vdelta_8/sigma_8 can define growth
    :ivar transfer_data: numpy array T[entry, q_index, z_index] storing transfer functions for each redshift and q; entry+1 can be

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
    """

    def transfer_z(self, name, z_index=0):
        """
        Get transfer function (function of q, for each q in self.q_trans) by name for given redshift index

        :param name:  parameter name
        :param z_index: which redshift
        :return: array of transfer function values for each calculated k
        """

        if not name in model.transfer_names:
            raise CAMBError('Unknown name %s; must be one of %s' % (name, model.transfer_names))
        return self.transfer_data[model.transfer_names.index(name), :, z_index]


class ClTransferData(object):
    r"""
    ClTransferData is the base class for storing CMB power transfer functions, as a function of q and :math:`\ell`.
    To get an instance of this data, call :meth:`camb.CAMBdata.get_cmb_transfer_data`

    :ivar NumSources:  number of sources calculated (size of p index)
    :ivar q: array of q values calculated (=k in flat universe)
    :ivar l: int array of :math:`\ell` values calculated
    :ivar delta_p_l_k: transfer functions, indexed by source, l, q
    """

    def get_transfer(self, source=0):
        r"""
        Return :math:`C_\ell` trasfer functions as a function of :math:`\ell`
        and :math:`q` (:math:`= k` in a flat universe).

        :param source: index of source: e.g. 0 for temperature, 1 for E polarization, 2 for lensing potential
        :return: array of computed l, array of computed q, transfer functions T(l,q)
        """

        return self.l, self.q, self.delta_p_l_k[source, :, :]


def set_feedback_level(level=1):
    """
    Set the feedback level for internal CAMB calls
    :param level:  zero for nothing, >1 for more
    """
    FeedbackLevel.value = level


def set_default_params(P):
    """
    Set default values for all parameters
    :param P: :class:`.model.CAMBparams`
    :return: P
    """
    assert (isinstance(P, model.CAMBparams))
    camblib.__camb_MOD_camb_setdefparams(byref(P))
    return P


def fortran_array(c_pointer, shape, dtype=np.float64, order='F', own_data=True):
    if not hasattr(shape, '__len__'):
        shape = np.atleast_1d(shape)
    arr_size = np.prod(shape[:]) * np.dtype(dtype).itemsize
    if sys.version_info.major >= 3:
        buf_from_mem = ctypes.pythonapi.PyMemoryView_FromMemory
        buf_from_mem.restype = ctypes.py_object
        buf_from_mem.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        buffer = buf_from_mem(c_pointer, arr_size, 0x100)
    else:
        buffer_from_memory = ctypes.pythonapi.PyBuffer_FromMemory
        buffer_from_memory.restype = ctypes.py_object
        buffer = buffer_from_memory(c_pointer, arr_size)
    arr = np.ndarray(tuple(shape[:]), dtype, buffer, order=order)
    if own_data and not arr.flags.owndata:
        return arr.copy()
    else:
        return arr


def set_z_outputs(z_outputs):
    """
    Set the redshifts for calculating BAO parameters at

    :param z_outputs: array of redshifts
    """
    z_outputs = np.array(z_outputs)
    CAMB_SetBackgroundOutputs_z(z_outputs, byref(c_int(len(z_outputs))))


class CAMBdata(object):
    """
    An object for storing transfer function data and parameters for CAMB.
    Not that it *only* stores transfer functions. If you want to get power spectra or background functions,
    you must have called one of the calculation functions for the parameters of interest more recently than
    any other call to these functions. You can can make multiple instances of CAMBdata and then later call
    :func:`~camb.CAMBdata.power_spectra_from_transfer` to calculate other quantities.

    To quickly make a fully calculated CAMBdata instance for a set of parameters you can call :func:`get_results`.

    :ivar Params: the :class:`.model.CAMBparams` parameters being used

    """

    def __init__(self):
        self._key = POINTER(_CAMBdata)()
        CAMBdata_new(byref(self._key))
        self.Params = self.get_params()
        self._one = c_int(1)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def __del__(self):
        self.free()

    def free(self):
        if self._key:
            CAMBdata_free(byref(self._key))
            self._key = None

    def set_params(self, params):
        """
        Set parameters from params. Note that this does not recompute anything;
        you will need to call :meth:`calc_transfers` if you change any parameters affecting the
        background cosmology or the transfer function settings.

        :param params: a :class:`.model.CAMBparams` instance
        """
        assert (isinstance(params, model.CAMBparams))
        CAMBdata_setparams(self._key, byref(params))

    def get_params(self):
        """
        Get the parameters currently set. Returned object references stored data, so elements can be modified without
        calling set_params again.

        :return: :class:`.model.CAMBparams` instance pointing to the underlying parameters used by CAMB.
        """

        p = POINTER(model.CAMBparams)()
        CAMBdata_getparams(self._key, byref(p))
        return p.contents

    def get_derived_params(self):
        """
        :return: dictionary of derived parameter values, indexed by name ('kd', 'age', etc..)
        """
        res = {}
        for name, value in zip(model.derived_names, model.ThermoDerivedParams):
            res[name] = value
        return res

    def get_background_outputs(self):
        """
        Get BAO values for redshifts (set redshifts using :func:`set_z_outputs`)

        :return: rs/DV, H, DA, F_AP for each requested redshift (as 2D array)
        """
        n = CAMB_GetNumBackgroundOutputs()
        if not n:
            raise CAMBError(
                'Call camb.set_z_outputs with required redshifts (and then calculate transfers/results) before calling get_background_outputs')
        outputs = np.empty((n, 4))
        CAMB_GetBackgroundOutputs(outputs, byref(c_int(n)))
        return outputs

    def get_BAO(self, redshifts, params):
        """
        Get BAO parameters at given redshifts, using parameters in params

        :param redshifts: list of redshifts
        :param params: optional :class:`.model.CAMBparams` instance to use
        :return: array of rs/DV, H, DA, F_AP for each redshift as 2D array
        """
        set_z_outputs(redshifts)
        self.calc_background(params)
        res = self.get_background_outputs()
        set_z_outputs([])
        return res

    def calc_background_no_thermo(self, params):
        """
        Calculate the background evolution without calculating thermal history.
        e.g. call this if you want to just use :meth:`angular_diameter_distance` and similar background functions

        :param params:  :class:`.model.CAMBparams` instance to use
        """
        CAMB_SetParamsForBackground(self._key, byref(params))

    def calc_background(self, params):
        """
        Calculate the background evolution and thermal history.
        e.g. call this if you want to get derived parameters and call background functions
        :param params:  :class:`.model.CAMBparams` instance to use
        """
        res = CAMB_CalcBackgroundTheory(self._key, byref(params))
        if res:
            raise CAMBError('Error %s in calc_background' % res)

    def calc_transfers(self, params, only_transfers=True):
        """
        Calculate the transfer functions (for CMB and matter power, as determined by params.WantCls, params.WantTransfer)

        :param params: :class:`.model.CAMBparams` instance with parameters to use
        :param only_transfers: only calculate transfer functions, no power spectra
        :return: non-zero if error, zero if OK
        """
        opt = c_bool()
        opt.value = only_transfers
        if not only_transfers: self._check_powers(params)
        return CAMBdata_gettransfers(self._key, byref(params), byref(opt))

    def _check_powers(self, params=None):
        if params is None: params = self.Params
        if params.InitPower.has_tensors() and not params.WantTensors:
            raise CAMBError('r>0 but params.WantTensors = F')
        if params.WantScalars and params.WantCls and params.DoLensing and params.scalar_power(0.05) > 2e-8:
            raise CAMBError('Lensing requires a realistically normalized spectrum, you have P(k=0.05/Mpc) > 2e-8')

    def calc_power_spectra(self, params=None):
        """
        Calculates transfer functions and power spectra.

        :param params: optional :class:`.model.CAMBparams` instance with parameters to use

        """
        if params is not None:
            result = self.calc_transfers(params, only_transfers=False)
            if result != 0:
                raise CAMBError('Error getting transfer functions: %u' % result)
        else:
            self._check_powers()
            CAMBdata_transferstopowers(self._key)

    def power_spectra_from_transfer(self, initial_power_params):
        """
        Assuming :meth:`calc_transfers` or :meth:`calc_power_spectra` have already been used, re-calculate the power spectra
        using a new set of initial power spectrum parameters with otherwise the same cosmology.
        This is typically much faster that re-calculating everything, as the transfer functions can be re-used.

        :param initial_power_params: :class:`.initialpower.InitialPowerParams` instance with new primordial power spectrum parameters
        """
        self.Params.set_initial_power(initial_power_params)
        self._check_powers()
        CAMBdata_transferstopowers(self._key)

    def _CMB_unit(self, CMB_unit):
        if isinstance(CMB_unit, six.string_types):
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
                cls[1:, 0:] /= ls ** 2 / (2 * np.pi)
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
            lmax_calc = model.lmax_lensed.value
        else:
            lmax_calc = self.Params.max_l
        if lmax is None:
            lmax = lmax_calc
        elif lmax > lmax_calc:
            logging.warning('getting CMB power spectra to higher L than calculated, may be innacurate/zeroed.')
        return lmax

    def save_cmb_power_spectra(self, filename, lmax, CMB_unit='muK'):
        r"""
        Save CMB power to a plain text file. Output is lensed total :math:`\ell(\ell+1)C_\ell/2\pi` then lensing potential and cross: L TT EE BB TE PP PT PE.

        :param filename: filename to save
        :param lmax: lmax to save
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell` and :math:`\mu K` units for lensing cross.
        """
        cmb = self.get_total_cls(lmax, CMB_unit=CMB_unit)
        lens = self.get_lens_potential_cls(lmax, CMB_unit=CMB_unit)
        ls = np.atleast_2d(np.arange(lmax + 1)).T
        np.savetxt(filename, np.hstack((ls, cmb, lens)), fmt=['%4u'] + ['%12.7e'] * 7,
                   header=' L ' + 'TT EE BB TE PP PT PE'.replace(' ', ' ' * 12))

    def get_cmb_power_spectra(self, params=None, lmax=None,
                              spectra=['total', 'unlensed_scalar', 'unlensed_total', 'lensed_scalar', 'tensor',
                                       'lens_potential'], CMB_unit=None, raw_cl=False):
        r"""
        Get CMB power spectra, as requested by the 'spectra' argument. All power spectra are :math:`\ell(\ell+1)C_\ell/2\pi` self-owned
        numpy arrays (0..lmax, 0..3), where 0..3 index are TT, EE, BB TT, unless raw_cl is True in which case return just :math:`C_\ell`.
        For the lens_potential the power spectrum returned is that of the deflection.

        :param params: optional :class:`.model.CAMBparams` instance with parameters to use. If None, must have
          previously set parameters and called `calc_power_spectra` (e.g. if you got this instance using :func:`get_results`),
        :param lmax: maximum l
        :param spectra: list of names of spectra to get
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell` and :math:`\mu K` units for lensing cross.
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

        :param params: optional :class:`.model.CAMBparams` instance with parameters to use. If None, must have
          previously set parameters and called :meth:`calc_power_spectra` (e.g. if you got this instance using :func:`get_results`),
        :param lmax: optional maximum L to use from the cls arrays
        :param spectrum: type of CMB power spectrum to get; default 'lensed_scalar', one of
          ['total', 'unlensed_scalar', 'unlensed_total', 'lensed_scalar', 'tensor']
        :param xvals: optional array of :math:`\cos(\theta)` values at which to calculate correlation function.
        :param sampling_factor: multiple of lmax for the Gauss-Legendre order if xvals not given (default 1)
        :return: if xvals not given: corrs, xvals, weights; if xvals specified, just corrs.
          corrs is 2D array corrs[i, ix], where ix=0,1,2,3 are T, Q+U, Q-U and cross, and i indexes xvals
        """

        if not spectrum in ['total', 'unlensed_scalar', 'unlensed_total', 'lensed_scalar', 'tensor']:
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
        CAMBdata_cltransferdata(self._key, byref(cdata), byref(c_int(['scalar', 'vector', 'tensor'].index(tp))))
        data = ClTransferData()
        data.NumSources = cdata.NumSources
        data.q = fortran_array(cdata.q, cdata.q_size)
        data.l = fortran_array(cdata.l, cdata.l_size, dtype=c_int)
        data.delta_p_l_k = fortran_array(cdata.delta_p_l_k, cdata.delta_size)
        return data

    def get_time_evolution(self, q, eta, vars=model.evolve_names, lAccuracyBoost=4, frame='CDM'):
        """
        Get the mode evolution as a function of conformal time for some k values.

        :param q: wavenumber values to calculate (or array of k values)
        :param eta: array of requested conformal times to output
        :param vars: list of variable names or sympy symbolic expressions to output (using camb.symbolic)
        :param lAccuracyBoost: factor by which to increase l_max in hierarchies compared to default - often
          needed to get nice smooth curves of acoustic oscillations for plotting.
        :param frame: for symbolic expressions, can specify frame name if the variable is not gauge invariant.
            e.g. specifying Delta_g and frame='Newtonian' would give the Newtonian gauge photon density perturbation.
        :return: nd array, A_{qti}, size(q) x size(times) x len(vars), or 2d array if q is scalar
        """

        old_boost = model._lAccuracyBoost.value
        try:
            if lAccuracyBoost: model._lAccuracyBoost.value = lAccuracyBoost
            if not isinstance(vars, (tuple, list)):
                vars = [vars]
            import sympy
            named_vars = [var for var in vars if isinstance(var, six.string_types)]

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
                custom_source_func = ctypes.cast(funcPtr, ctypes.c_voidp)
            else:
                custom_source_func = ctypes.c_voidp(0)
            nvars = num_standard_names + ncustom
            outputs = np.empty((k.shape[0], times.shape[0], nvars))
            if CAMB_TimeEvolution(byref(c_int(k.shape[0])), k, byref(c_int(times.shape[0])), times[indices],
                                  byref(c_int(nvars)), outputs,
                                  byref(c_int(ncustom)), byref(custom_source_func)):
                raise CAMBError('Error in evolution')
            i_rev = np.zeros(times.shape, dtype=int)
            i_rev[indices] = np.arange(times.shape[0])
            outputs = outputs[:, i_rev, :]
        finally:
            model._lAccuracyBoost.value = old_boost
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

        if isinstance(vars, six.string_types): vars = [vars]
        unknown = set(vars) - set(model.background_names)
        if unknown:
            raise CAMBError('Unknown names %s; valid names are %s' % (unknown, model.background_names))
        outputs = np.zeros((eta.shape[0], 5))
        CAMB_BackgroundEvolution(byref(c_int(eta.shape[0])), eta, outputs)
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

    def get_matter_transfer_data(self):
        """
        Get matter transfer function data and sigma8 for calculated results.

        :return: :class:`.MatterTransferData` instance holding output arrays (copies, not pointers)
        """
        if not self.Params.WantTransfer:
            raise CAMBError("must have Params.WantTransfer to get matter transfers and power")

        cdata = _MatterTransferData()
        CAMBdata_mattertransferdata(self._key, byref(cdata))
        data = MatterTransferData()
        data.nq = cdata.num_q_trans
        data.q = nplib.as_array(cdata.q_trans, shape=(data.nq,))
        data.sigma_8 = fortran_array(cdata.sigma_8, cdata.sigma_8_size)
        data.sigma2_vdelta_8 = fortran_array(cdata.sigma2_vdelta_8, cdata.sigma2_vdelta_8_size)
        data.transfer_data = fortran_array(cdata.TransferData, cdata.TransferData_size, dtype=np.float32)
        return data

    def _transfer_var(self, var1, var2):
        if var1 is None: var1 = model.transfer_power_var.value
        if var2 is None: var2 = model.transfer_power_var.value
        if isinstance(var1, six.string_types): var1 = model.transfer_names.index(var1) + 1
        if isinstance(var2, six.string_types): var2 = model.transfer_names.index(var2) + 1
        return c_int(var1), c_int(var2)

    def get_linear_matter_power_spectrum(self, var1=None, var2=None,
                                         hubble_units=True, have_power_spectra=False, params=None, nonlinear=False):
        r"""
        Calculates :math:`P_{xy}(k/h)`, where x, y are one of model.Transfer_cdm, model.Transfer_xx etc.
        The output k values are not regularly spaced, and not interpolated.

        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param hubble_units: if true, output power spectrum in (Mpc/h) units, otherwise Mpc
        :param have_power_spectra: set to True if already computed power spectra
        :param params: if have_power_spectra=False, optional :class:`.model.CAMBparams` instance to specify new parameters
        :param nonlinear: include non-linear correction from halo model
        :return: kh, z, PK, where kz an z are arrays of k/h and z respectively, and PK[i,j] is value at z[i], k/h[j]
        """
        if not have_power_spectra:
            self.calc_power_spectra(params)
        data = self.get_matter_transfer_data()

        nk = data.nq
        nz = self.Params.Transfer.PK_num_redshifts
        kh = data.transfer_data[model.Transfer_kh - 1, :, 0]

        var1, var2 = self._transfer_var(var1, var2)

        hubble_units = c_int(hubble_units)
        PK = np.empty((nz, nk))
        if nonlinear:
            CAMBdata_GetNonLinearMatterPower(self._key, PK, byref(var1), byref(var2), byref(hubble_units))
        else:
            CAMBdata_GetLinearMatterPower(self._key, PK, byref(var1), byref(var2), byref(hubble_units))

        z = self.Params.Transfer.PK_redshifts[:nz]
        z.reverse()
        return np.array(kh), np.array(z), PK

    def get_nonlinear_matter_power_spectrum(self, **kwargs):
        r"""
        Calculates :math:`P_{xy}(k/h)`, where x, y are one of model.Transfer_cdm, model.Transfer_xx etc.
        The output k values are not regularly spaced, and not interpolated.

        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param hubble_units: if true, output power spectrum in :Math:`({\rm Mpc}/h)^{3}` units, otherwise :math:`{\rm Mpc}^{3}`
        :param have_power_spectra: set to True if already computed power spectra
        :param params: if have_power_spectra=False, optional :class:`.model.CAMBparams` instance to specify new parameters
        :return: kh, z, PK, where kz an z are arrays of k/h and z respectively, and PK[i,j] is value at z[i], k/h[j]
        """
        kwargs['nonlinear'] = True
        return self.get_linear_matter_power_spectrum(**kwargs)

    def get_sigma8(self):
        r"""
        Get :math:`\sigma_8` values (must previously have calculated power spectra)

        :return: array of :math:`\sigma_8` values, in order of increasing time (decreasing redshift)
        """
        mtrans = self.get_matter_transfer_data()
        return mtrans.sigma_8[:, 0]

    def get_fsigma8(self):
        r"""
        Get :math:`f\sigma_8` growth values (must previously have calculated power spectra).
        For general models :math:`f\sigma_8` is defined as in the Planck 2015 parameter paper in terms of
        the velocity-density correlation: :math:`\sigma^2_{vd}/\sigma_{dd}` for :math:`8 h^{-1} {\rm Mpc}` spheres.

        :return: array of f*sigma_8 values, in order of increasing time (decreasing redshift)
        """
        mtrans = self.get_matter_transfer_data()
        return mtrans.sigma2_vdelta_8[:, 0] / mtrans.sigma_8[:, 0]

    def get_matter_power_spectrum(self, minkh=1e-4, maxkh=1.0, npoints=100,
                                  var1=None, var2=None,
                                  have_power_spectra=False, params=None):
        """
        Calculates :math:`P_{xy}(k/h)`, where x, y are one of Transfer_cdm, Transfer_xx etc defined in ModelParams.
        The output k values are regularly log spaced and interpolated. If NonLinear is set, the result is non-linear.

        :param minkh: minimum value of k/h for output grid (very low values < 1e-4 may not be calculated)
        :param maxkh: maximum value of k/h (check consistent with input params.Transfer.kmax)
        :param npoints: number of points equally spaced in log k
        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param have_power_spectra: set to True if already computed power spectra
        :param params: if have_power_spectra=False and want to specify new parameters, a :class:`.model.CAMBparams` instance
        :return: kh, z, PK, where kz an z are arrays of k/h and z respectively, and PK[i,j] is value at z[i], k/h[j]
        """

        if not have_power_spectra:
            self.calc_power_spectra(params)

        assert self.Params.WantTransfer
        if self.Params.Transfer.kmax < maxkh:
            logging.warning("get_matter_power_spectrum using larger k_max than input parameter Transfer.kmax")
        if self.Params.NonLinear == model.NonLinear_none and self.Params.Transfer.kmax < 1:
            logging.warning("get_matter_power_spectrum Transfer.kmax small to get non-linear spectrum")

        nz = self.Params.Transfer.PK_num_redshifts
        PK = np.empty((nz, npoints))
        var1, var2 = self._transfer_var(var1, var2)

        dlnkh = (np.log(maxkh) - np.log(minkh)) / (npoints - 1)
        CAMBdata_GetMatterPower(self._key, PK, byref(c_double(minkh)),
                                byref(c_double(dlnkh)), byref(c_int(npoints)), byref(var1), byref(var2))
        z = self.Params.Transfer.PK_redshifts[:nz]
        z.reverse()
        return minkh * np.exp(np.arange(npoints) * dlnkh), z, PK

    def get_matter_power_interpolator(self, nonlinear=True, var1=None, var2=None, hubble_units=True, k_hunit=True,
                                      return_z_k=False, log_interp=True, extrap_kmax=None):
        r"""
        Assuming transfers have been calculated, return a 2D spline interpolation object to evaluate matter
        power spectrum as function of z and k/h (or k)
        e.g::
          PK = results.get_matter_power_interpolator();
          print 'Power spectrum at z=0.5, k/h=0.1 is %s (Mpc/h)^3 '%(PK.P(0.5, 0.1))

        :param nonlinear: include non-linear correction from halo model
        :param var1: variable i (index, or name of variable; default delta_tot)
        :param var2: variable j (index, or name of variable; default delta_tot)
        :param hubble_units: if true, output power spectrum in :math:`({\rm Mpc}/h)^{3}` units, otherwise :math:`{\rm Mpc}^{3}`
        :param k_hunit: if true, matter power is a function of k/h, if false, just k (both :math:`{\rm Mpc}^{-1}` units)
        :param return_z_k: if true, return interpolator, z, k where z, k are the grid used
        :param log_interp: if true, interpolate log of power spectrum (unless any values are negative in which case ignored)
        :param extrap_kmax: if set, use power law extrapolation beyond kmax to extrap_kmax (useful for tails of integrals)
        :return: RectBivariateSpline object PK, that can be called with PK(z,log(kh)) to get log matter power values.
            if return_z_k=True, instead return interpolator, z, k where z, k are the grid used
        """

        class PKInterpolator(RectBivariateSpline):

            def P(self, z, kh, grid=None):
                if grid is None:
                    grid = not np.isscalar(z) and not np.isscalar(kh)
                if self.islog:
                    return np.exp(self(z, np.log(kh), grid=grid))
                else:
                    return self(z, np.log(kh), grid=grid)

        assert self.Params.WantTransfer
        kh, z, pk = self.get_linear_matter_power_spectrum(var1, var2, hubble_units, nonlinear=nonlinear)
        if not k_hunit:
            kh *= self.Params.H0 / 100
        if log_interp and np.any(pk <= 0):
            log_interp = False
        logkh = np.log(kh)
        if extrap_kmax and extrap_kmax > kh[-1]:
            logextrap = np.log(extrap_kmax)
            logpknew = np.empty((pk.shape[0], pk.shape[1] + 1))
            logpknew[:, :-1] = np.log(pk)
            logpknew[:, -1] = logpknew[:, -2] + (logpknew[:, -2] - logpknew[:, -3]) / (logkh[-2] - logkh[-3]) * (
                    logextrap - logkh[-1])
            logkhnew = np.hstack((logkh, logextrap))
            if log_interp:
                res = PKInterpolator(z, logkhnew, logpknew)
            else:
                res = PKInterpolator(z, logkhnew, np.exp(logpknew))
            res.kmax = extrap_kmax
        else:
            if log_interp:
                res = PKInterpolator(z, logkh, np.log(pk))
            else:
                res = PKInterpolator(z, logkh, pk)
            res.kmax = np.max(kh)

        res.kmin = np.min(kh)
        res.islog = log_interp
        res.zmin = np.min(z)
        res.zmax = np.max(z)
        if return_z_k:
            return res, z, kh
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
        CAMB_SetTotCls(byref(opt), res, byref(self._one))
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
        CAMB_SetTensorCls(byref(opt), res, byref(self._one))
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
        CAMB_SetUnlensedScalCls(byref(opt), res, byref(self._one))
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
        return self.get_unlensed_scalar_cls(lmax, CMB_unit, raw_cl) + \
               self.get_tensor_cls(lmax, CMB_unit, raw_cl)

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
        CAMB_SetLensedScalCls(byref(opt), res, byref(self._one))
        self._scale_cls(res, CMB_unit, raw_cl)
        return res

    def get_lens_potential_cls(self, lmax=None, CMB_unit=None, raw_cl=False):
        r"""
        Get lensing deflection angle potential power spectrum, and cross-correlation with T and E. Must have already calculated power spectra.
        Power spectra are :math:`[L(L+1)]^2C_L^{\phi\phi}/2\pi` and corresponding deflection cross-correlations.

        :param lmax: lmax to output to

        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K` units for lensing cross.
        :param raw_cl: return lensing potential :math:`C_L` rather than :math:`[L(L+1)]^2C_L/2\pi`
        :return: numpy array CL[0:lmax+1,0:3], where 0..2 indexes PP, PT, PE.
        """

        lmax = self._lmax_setting(lmax, unlensed=True)
        res = np.empty((lmax + 1, 3))
        opt = c_int(lmax)
        CAMB_SetLensPotentialCls(byref(opt), res, byref(self._one))
        self._scale_cls(res, CMB_unit, raw_cl, lens_potential=True)
        return res

    def get_unlensed_scalar_array_cls(self, lmax=None):
        """
        Get array of all cross power spectra. Must have already calculated power spectra.
        Results are dimensionless, and not scaled by custom_scaled_ell_fac.

        :param lmax: lmax to output to
        :return: numpy array CL[0:, 0:,0:lmax+1], where 0.. index T, E, deflection angle, source window functions
        """

        lmax = self._lmax_setting(lmax, unlensed=True)
        if not model.has_cl_2D_array.value:
            raise CAMBError('unlensed_scalar_array not calculated (set model.has_cl_2D_array)')
        n = 3 + model.num_redshiftwindows.value + len(custom_source_names)
        res = np.empty((n, n, lmax + 1), order='F')
        CAMB_SetUnlensedScalarArray(byref(c_int(lmax)), res, byref(self._one), byref(c_int(n)))
        return res

    def get_cmb_unlensed_scalar_array_dict(self, params=None, lmax=None, CMB_unit=None, raw_cl=False):
        """
        Get all unlensed auto and cross power spectra, including any custom source functions set using :func:`set_custom_scalar_sources`.

        :param params: optional :class:`.model.CAMBparams` instance with parameters to use. If None, must have
          previously set parameters and called :meth:`calc_power_spectra` (e.g. if you got this instance using :func:`get_results`),
        :param lmax: maximum :math:`\ell`
        :param CMB_unit: scale results from dimensionless. Use 'muK' for :math:`\mu K^2` units for CMB :math:`C_\ell` and :math:`\mu K` units for lensing cross.
        :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
        :return: dictionary of power spectrum arrays, index as TxT, TxE, W1xW2, custom_name_1xT... etc.
        """

        old_val = model.has_cl_2D_array.value
        try:
            nwindows = model.num_redshiftwindows.value
            if params is not None:
                model.has_cl_2D_array.value = True
                self.calc_power_spectra(params)
            elif not model.has_cl_2D_array:
                raise CAMBValueError('model.has_cl_2D_array must be true to have array C_L')
            lmax = lmax or self.Params.max_l
            arr = self.get_unlensed_scalar_array_cls(lmax)
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
                            ls = np.arange(1, cls.shape[0])[..., np.newaxis]
                            fac = np.float64(ls * (ls + 1))
                            if i == 3 and j == 3:
                                fac *= fac
                            elif i == 3 or j == 3:
                                fac *= np.sqrt(fac)
                            cls[1:] /= (fac / (2 * np.pi))

                        if CMB_unit is not None:
                            cls *= CMB_units[i] * CMB_units[j]
                        result[tag] = cls
        finally:
            model.has_cl_2D_array.value = old_val
        return result

    def angular_diameter_distance(self, z):
        """
        Get (non-comoving) angular diameter distance to redshift z.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions or power spectra.

        :param z: redshift or array of redshifts
        :return: angular diameter distances, matching rank of z
        """
        if np.isscalar(z):
            return AngularDiameterDistance(byref(c_double(z)))
        else:
            z = np.asarray(z)
            arr = np.empty(z.shape)
            indices = np.argsort(z)
            redshifts = np.array(z[indices], dtype=np.float64)
            AngularDiameterDistanceArr(arr, redshifts, byref(c_int(z.shape[0])))
            arr[indices] = arr.copy()
            return arr

    def angular_diameter_distance2(self, z1, z2):
        r"""
        Get angular diameter distance between two redshifts
        :math:`\frac{r}{1+z_2}\text{sin}_K\left(\frac{\chi(z_2) - \chi(z_1)}{r}\right)`
        where :math:`r` is curvature radius and :math:`\chi` is the comoving radial distance.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions or power spectra.

        :param z1: redshift 1
        :param z2: redshift 2
        :return: result
        """
        if not np.isscalar(z1) or not np.isscalar(z2):
            raise CAMBError('vector z not supported yet')
        return AngularDiameterDistance2(byref(c_double(z1)), byref(c_double(z2)))

    def comoving_radial_distance(self, z, tol=1e-4):
        """
        Get comoving radial distance from us to redshift z in Mpc. This is efficient for arrays.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions or power spectra.

        :param z: redshift
        :return: comoving radial distance (Mpc)
        """
        if not np.isscalar(z):
            indices = np.argsort(z)
            redshifts = np.array(z[indices], dtype=np.float64)
            chis = np.empty(redshifts.shape)
            ComovingRadialDistanceArr(chis, redshifts, byref(c_int(chis.shape[0])), byref(c_double(tol)))
            chis[indices] = chis.copy()
            return chis
        else:
            return ComovingRadialDistance(byref(c_double(z)))

    def redshift_at_comoving_radial_distance(self, chi, nz_step=150, zmax=10000):
        """
        Convert comoving radial distance array to redshift array.
        This is not calculated directly, but fit from a spline to a forward calculation of chi from z.
        This is a utility routine, not optimized to be fast (though it can work on a large vector efficiently)

        :param chi: comoving radial distance (in Mpc), scalar or array
        :param nz_step: number of redshifts calculated internally for solving grid
        :param zmax: maximum redshift in internal solving grid
        :return: redshift at chi, scalar or array
        """

        zs = np.exp(np.log(zmax + 1) * np.linspace(0, 1, nz_step)) - 1
        chis = self.comoving_radial_distance(zs, tol=1e-5)
        f = UnivariateSpline(chis, zs, s=0)
        if np.isscalar(chi):
            return np.asscalar(f(chi))
        else:
            return f(chi)

    def luminosity_distance(self, z):
        """
        Get luminosity distance from to redshift z.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions or power spectra.

        :param z: redshift or array of redshifts
        :return: luminosity distance (matches rank of z)
        """

        if not np.isscalar(z): z = np.asarray(z)
        return self.angular_diameter_distance(z) * (1.0 + z) ** 2

    def h_of_z(self, z):
        r"""
        Get Hubble rate at redshift z, in :math:`{\rm Mpc}^{-1}` units, scalar or array

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions or power spectra.

        Use hubble_parameter instead if you want in [km/s/Mpc] units.

        :param z: redshift
        :return: H(z)
        """
        if not np.isscalar(z):
            z = np.array(z, dtype=np.float64)
            arr = np.empty(z.shape)
            HofzArr(arr, z, byref(c_int(z.shape[0])))
            return arr
        else:
            return Hofz(byref(c_double(z)))

    def hubble_parameter(self, z):
        """
        Get Hubble rate at redshift z, in km/s/Mpc units. Scalar or array.

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions or power spectra.

        :param z: redshift
        :return: H(z)/[km/s/Mpc]
        """
        return constants.c * self.h_of_z(z) / 1e3

    def physical_time_a1_a2(self, a1, a2):
        """
        Get physical time between two scalar factors in Gigayears

        Must have called :meth:`calc_background`, :meth:`calc_background_no_thermo` or calculated transfer functions or power spectra.

        :param a1: scale factor 1
        :param a2: scale factor 2
        :return: (age(a2)-age(a1))/Gigayear
        """
        if not np.isscalar(a1) or not np.isscalar(a2):
            raise CAMBError('vector inputs not supported yet')
        return DeltaPhysicalTimeGyr(byref(c_double(a1)), byref(c_double(a2)), None)

    def physical_time(self, z):
        """
        Get physical time from hot big bang to redshift z in Gigayears.

        :param z:  redshift
        :return: t(z)/Gigayear
        """
        return self.physical_time_a1_a2(0, 1.0 / (1 + z))

    def conformal_time_a1_a2(self, a1, a2):
        """
        Get conformal time between two scale factors (=comoving radial distance travelled by light on light cone)

        :param a1: scale factor 1
        :param a2: scale factor 2
        :return: eta(a2)-eta(a1) = chi(a1)-chi(a2) in Megaparsec
        """

        if not np.isscalar(a1) or not np.isscalar(a2):
            raise CAMBError('vector inputs not supported yet')
        return DeltaTime(byref(c_double(a1)), byref(c_double(a2)), None)

    def conformal_time(self, z):
        """
        Conformal time from hot big bang to redshift z in Megaparsec.
        Use comoving_radial_distance for faster result for arrays.

        :param z: redshift or array of redshifts
        :return: eta(z)/Mpc
        """
        if np.isscalar(z):
            redshifts = np.array([z], dtype=np.float64)
        else:
            redshifts = np.array(z, dtype=np.float64)
        eta = np.empty(redshifts.shape)
        TimeOfzArr(byref(c_int(eta.shape[0])), redshifts, eta)
        if np.isscalar(z):
            return eta[0]
        else:
            return eta

    def cosmomc_theta(self):
        r"""
        Get :math:`\theta_{\rm MC}`, an approximation of the ratio of the sound horizon to the angular diameter distance at recombination.

        :return: :math:`\theta_{\rm MC}`
        """
        return CosmomcTheta()


def get_results(params):
    """
    Calculate results for specified parameters and return :class:`CAMBdata` instance for getting results.

    :param params: :class:`.model.CAMBparams` instance
    :return: :class:`CAMBdata` instance
    """
    res = CAMBdata()
    if _debug_params: print(params)
    res.calc_power_spectra(params)
    return res


def get_transfer_functions(params):
    """
    Calculate transfer functions for specified parameters and return :class:`CAMBdata` instance for getting results
    and subsequently calculating power spectra.

    :param params: :class:`.model.CAMBparams` instance
    :return: :class:`CAMBdata` instance
    """

    res = CAMBdata()
    res.calc_transfers(params)
    return res


def get_background(params, no_thermo=False):
    """
    Calculate background cosmology for specified parameters and return :class:`CAMBdata`, ready to get derived
     parameters and use background functions like :func:`~camb.CAMBdata.angular_diameter_distance`.

    :param params: :class:`.model.CAMBparams` instance
    :param no_thermo: Use calc_background_no_thermo instead
    :return: :class:`CAMBdata` instance
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
    cTau = c_double(tau)
    return CAMB_GetZreFromTau(byref(params), byref(cTau))


def cleanup():
    camblib.__camb_MOD_camb_cleanup()


def set_params(cp=None, verbose=False, **params):
    """

    Set all CAMB parameters at once, including parameters which are part of the
    CAMBparams structure, as well as global parameters.

    E.g.::

      cp = camb.set_params(ns=1, omch2=0.1, w=-0.95, ALens=1.2, lmax=2000)

    This is equivalent to::

      cp = model.CAMBparams()
      cp.set_dark_energy(w=-0.95)
      cp.set_cosmology(omch2=0.1)
      cp.set_for_lmax(lmax=2000)
      cp.InitPower.set_params(ns=1)
      lensing.ALens.value = 1.2


    :param **params: the values of the parameters
    :param cp: use this CAMBparams instead of creating a new one
    :param verbose: print out the equivalent set of commands

    """

    if cp is None:
        cp = model.CAMBparams()
    else:
        assert isinstance(cp, model.CAMBparams), "cp should be an instance of CAMBparams"

    _used_params = set()
    if 'ALens' in params:
        _used_params.add('ALens')
        lensing.ALens.value = params['ALens']

    #Note order is important: must call set_dark_energy before set_cosmology if setting cosmomc_theta
    setters = [cp.set_accuracy, cp.set_dark_energy, cp.set_cosmology, cp.set_initial_power,
               cp.set_matter_power, cp.set_for_lmax, cp.InitPower.set_params]

    for setter in setters:
        kwargs = {k: params[k] for k in getargspec(setter).args[1:] if k in params}
        _used_params.update(kwargs.keys())
        if kwargs:
            if verbose:
                logging.warning('Calling %s(**%s)' % (setter.__name__, kwargs))
            setter(**kwargs)

    if cp.InitPower.has_tensors():
        cp.WantTensors = True

    unused_params = set(params) - set(_used_params)
    if unused_params:
        raise CAMBUnknownArgumentError("Unrecognized parameters: %s" % unused_params)
    return cp


def get_matter_power_interpolator(params, zmin=0, zmax=10, nz_step=100, zs=None, kmax=10, nonlinear=True,
                                  var1=None, var2=None, hubble_units=True, k_hunit=True,
                                  return_z_k=False, k_per_logint=None, log_interp=True, extrap_kmax=None):
    r"""
    Return a 2D spline interpolation object to evaluate matter power spectrum as function of z and k/h
    e.g::
      from camb import get_matter_power_interpolator
      PK = get_matter_power_interpolator(params);
      print('Power spectrum at z=0.5, k/h=0.1/Mpc is %s (Mpc/h)^3 '%(PK.P(0.5, 0.1)))

    :param params: :class:`.model.CAMBparams` instance
    :param zmin: minimum z (use 0 or smaller than you want for good interpolation)
    :param zmax: maximum z (use larger than you want for good interpolation)
    :param nz_step: number of steps to sample in z (default max allowed is 100)
    :param zs: instead of zmin,zmax, nz_step, can specific explicit array of z values to spline from
    :param kmax: maximum k
    :param nonlinear: include non-linear correction from halo model
    :param var1: variable i (index, or name of variable; default delta_tot)
    :param var2: variable j (index, or name of variable; default delta_tot)
    :param hubble_units: if true, output power spectrum in :math:`({\rm Mpc}/h)^{3}` units, otherwise :math:`{\rm Mpc}^{3}`
    :param k_hunit: if true, matter power is a function of k/h, if false, just k (both :math:`{\rm Mpc}^{-1}` units)
    :param return_z_k: if true, return interpolator, z, k where z, k are the grid used
    :param k_per_logint: specific uniform sampling over log k (if not set, uses optimized irregular sampling)
    :param log_interp: if true, interpolate log of power spectrum (unless any values are negative in which case ignored)
    :param extrap_kmax: if set, use power law extrapolation beyond kmax to extrap_kmax (useful for tails of integrals)
    :return: RectBivariateSpline object PK, that can be called with PK(z,log(kh)) to get log matter power values.
        if return_z_k=True, instead return interpolator, z, k where z, k are the grid used
    """

    pars = copy.deepcopy(params)
    if zs is None:
        zs = zmin + np.exp(np.log(zmax - zmin + 1) * np.linspace(0, 1, nz_step)) - 1
    pars.set_matter_power(redshifts=zs, kmax=kmax, k_per_logint=k_per_logint, silent=True)
    pars.NonLinear = model.NonLinear_none
    results = get_results(pars)

    return results.get_matter_power_interpolator(nonlinear=nonlinear, var1=var1, var2=var2, hubble_units=hubble_units,
                                                 k_hunit=k_hunit, return_z_k=return_z_k, log_interp=log_interp,
                                                 extrap_kmax=extrap_kmax)


custom_source_names = []
_current_source_func = None


def set_custom_scalar_sources(custom_sources, source_names=None, source_ell_scales=None,
                              frame='CDM', code_path=None):
    r"""
    Set custom sources for angular power spectrum using camb.symbolic sympy expressions.

    :param custom_sources: list of sympy expressions for the angular power spectrum sources
    :param source_names: optional list of string naes for the sources
    :param source_ell_scales: list or dictionary of scalings for each source name, where for integer entry n, the source for
     multipole :math:`\ell` is scalled by :math:`\sqrt{(\ell+n)!/(\ell-n)!}`, i.e. :math:`n=2` for a new polarization-like source.
    :param frame: if the source is not gauge invariant, frame in which to interpret result
    :param code_path: optional path for output of source code for CAMB f90 source function
    """

    from . import symbolic
    global _current_source_func

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
        if source_names: source_names = [source_names]
    custom_source_names[:] = source_names or ["C%s" % (i + 1) for i in range(len(custom_sources))]
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

    _current_source_func = symbolic.compile_sympy_to_camb_source_func(custom_sources, frame=frame, code_path=code_path)
    custom_source_func = ctypes.cast(_current_source_func, ctypes.c_voidp)
    CAMB_SetCustomSourcesFunc(byref(c_int(len(custom_sources))), byref(custom_source_func), scales)


def clear_custom_scalar_sources():
    global _current_source_func
    custom_source_names[:] = []
    CAMB_SetCustomSourcesFunc(byref(c_int(0)), byref(ctypes.c_void_p(0)), np.zeros(0, dtype=np.int32))
    _current_source_func = None
