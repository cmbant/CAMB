import ctypes
from baseconfig import CAMBL, camblib, CAMBError
import model as model
import os
import os.path as osp
import numpy as np
from numpy import ctypeslib as nplib
from numpy.ctypeslib import ndpointer
import constants
import logging


# Use FeedbackLevel.value to read and set
FeedbackLevel = ctypes.c_int.in_dll(camblib, "__modelparams_MOD_feedbacklevel")

model.has_cl_2D_array.value = True


class _CAMBdata(ctypes.Structure):
    # contains complex types with pointers, so just set up dummy
    _fields_ = []


class _MatterTransferData(ctypes.Structure):
    # contains complex types with pointers, so just set up dummy
    _fields_ = [('num_q_trans', ctypes.c_int),
                ('q_trans', ctypes.POINTER(ctypes.c_double)),
                ('sigma_8', ctypes.POINTER(ctypes.c_double)),
                ('sigma2_vdelta_8', ctypes.POINTER(ctypes.c_double)),
                ('TransferData', ctypes.POINTER(ctypes.c_float)),
                ('sigma_8_size', ctypes.c_int * 2),
                ('sigma2_vdelta_8_size', ctypes.c_int * 2),
                ('TransferData_size', ctypes.c_int * 3)
                ]


int_arg = ctypes.POINTER(ctypes.c_int)
d_arg = ctypes.POINTER(ctypes.c_double)

# for the case where CAMB wrapper functions do the F-C conversion, so use C here
numpy_2d = ndpointer(ctypes.c_double, flags='C_CONTIGUOUS', ndim=2)
numpy_1d = ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')

CAMBdata_new = camblib.__handles_MOD_cambdata_new
CAMBdata_new.restype = ctypes.POINTER(_CAMBdata)

CAMBdata_free = camblib.__handles_MOD_cambdata_free
CAMBdata_free.argtypes = [ctypes.POINTER(ctypes.POINTER(_CAMBdata))]

CAMBdata_getparams = camblib.__handles_MOD_cambdata_getparams
CAMBdata_getparams.restype = ctypes.POINTER(model.CAMBparams)
CAMBdata_getparams.argtypes = [ctypes.POINTER(_CAMBdata)]

CAMBdata_setparams = camblib.__handles_MOD_cambdata_setparams
CAMBdata_setparams.argtypes = [ctypes.POINTER(_CAMBdata), ctypes.POINTER(model.CAMBparams)]

CAMBdata_gettransfers = camblib.__handles_MOD_cambdata_gettransfers
CAMBdata_gettransfers.argtypes = [ctypes.POINTER(_CAMBdata), ctypes.POINTER(model.CAMBparams),
                                  ctypes.POINTER(ctypes.c_bool)]
CAMBdata_gettransfers.restype = ctypes.c_int

CAMBdata_transferstopowers = camblib.__camb_MOD_camb_transferstopowers
CAMBdata_transferstopowers.argtypes = [ctypes.POINTER(_CAMBdata)]

CAMBdata_mattertransferdata = camblib.__handles_MOD_cambdata_mattertransferdata
CAMBdata_mattertransferdata.argtypes = [ctypes.POINTER(_CAMBdata), ctypes.POINTER(_MatterTransferData)]

CAMB_SetTotCls = camblib.__handles_MOD_camb_settotcls
CAMB_SetUnlensedCls = camblib.__handles_MOD_camb_setunlensedcls
CAMB_SetLensPotentialCls = camblib.__handles_MOD_camb_setlenspotentialcls
CAMB_SetUnlensedScalCls = camblib.__handles_MOD_camb_setunlensedscalcls
CAMB_SetLensedScalCls = camblib.__handles_MOD_camb_setlensedscalcls
CAMB_SetTensorCls = camblib.__handles_MOD_camb_settensorcls
CAMB_SetUnlensedScalarArray = camblib.__handles_MOD_camb_setunlensedscalararray

_set_cl_args = [ctypes.POINTER(ctypes.c_int), numpy_1d, int_arg]

CAMB_SetTotCls.argtypes = _set_cl_args
CAMB_SetUnlensedCls.argtypes = _set_cl_args
CAMB_SetLensPotentialCls.argtypes = _set_cl_args
CAMB_SetUnlensedScalCls.argtypes = _set_cl_args
CAMB_SetTensorCls.argtypes = _set_cl_args
CAMB_SetLensedScalCls.argtypes = _set_cl_args
CAMB_SetUnlensedScalarArray.argtypes = _set_cl_args + [int_arg]

del _set_cl_args

CAMB_SetBackgroundOutputs_z = camblib.__handles_MOD_camb_setbackgroundoutputs_z
CAMB_SetBackgroundOutputs_z.argtypes = [numpy_1d, int_arg]
CAMB_GetBackgroundOutputs = camblib.__handles_MOD_camb_getbackgroundoutputs
CAMB_GetBackgroundOutputs.argtypes = [numpy_1d, int_arg]
CAMB_GetNumBackgroundOutputs = camblib.__handles_MOD_camb_getnumbackgroundoutputs
CAMB_GetNumBackgroundOutputs.restype = ctypes.c_int

CAMB_GetAge = camblib.__camb_MOD_camb_getage
CAMB_GetAge.restype = ctypes.c_double
CAMB_GetAge.argtypes = [ctypes.POINTER(model.CAMBparams)]

CAMB_GetZreFromTau = camblib.__camb_MOD_camb_getzrefromtau
CAMB_GetZreFromTau.restype = ctypes.c_double
CAMB_GetZreFromTau.argtypes = [ctypes.POINTER(model.CAMBparams), d_arg]

CAMB_SetParamsForBackground = camblib.__handles_MOD_cambdata_setparamsforbackground
CAMB_SetParamsForBackground.argtypes = [ctypes.POINTER(_CAMBdata), ctypes.POINTER(model.CAMBparams)]

CAMB_CalcBackgroundTheory = camblib.__handles_MOD_cambdata_calcbackgroundtheory
CAMB_CalcBackgroundTheory.argtypes = [ctypes.POINTER(_CAMBdata), ctypes.POINTER(model.CAMBparams)]
CAMB_CalcBackgroundTheory.restype = ctypes.c_int

CAMBdata_GetLinearMatterPower = camblib.__handles_MOD_cambdata_getlinearmatterpower
CAMBdata_GetLinearMatterPower.argtypes = [ctypes.POINTER(_CAMBdata), numpy_2d, int_arg, int_arg, int_arg]

CAMBdata_GetMatterPower = camblib.__handles_MOD_cambdata_getmatterpower
CAMBdata_GetMatterPower.argtypes = [ctypes.POINTER(_CAMBdata), numpy_2d,
                                    d_arg, d_arg, int_arg, int_arg, int_arg]

AngularDiameterDistance = camblib.__modelparams_MOD_angulardiameterdistance
AngularDiameterDistance.argtyes = [d_arg]
AngularDiameterDistance.restype = ctypes.c_double

AngularDiameterDistanceArr = camblib.__modelparams_MOD_angulardiameterdistancearr
AngularDiameterDistanceArr.argtypes = [numpy_1d, numpy_1d, int_arg]

AngularDiameterDistance2 = camblib.__modelparams_MOD_angulardiameterdistance2
AngularDiameterDistance2.argtyes = [d_arg]
AngularDiameterDistance2.restype = ctypes.c_double

ComovingRadialDistance = camblib.__modelparams_MOD_comovingradialdistance
ComovingRadialDistance.argtyes = [d_arg]
ComovingRadialDistance.restype = ctypes.c_double


Hofz = camblib.__modelparams_MOD_hofz
Hofz.argtyes = [d_arg]
Hofz.restype = ctypes.c_double

DeltaPhysicalTimeGyr = camblib.__modelparams_MOD_deltaphysicaltimegyr
DeltaPhysicalTimeGyr.argtypes = [d_arg, d_arg, d_arg]
DeltaPhysicalTimeGyr.restype = ctypes.c_double

DeltaTime = camblib.__modelparams_MOD_deltatime
DeltaTime.argtypes = [d_arg, d_arg, d_arg]
DeltaTime.restype = ctypes.c_double

CosmomcTheta = camblib.__modelparams_MOD_cosmomctheta
CosmomcTheta.restype = ctypes.c_double


class MatterTransferData(object):
    """
    MatterTransferData is the base class for storing matter power transfer function data for various q values.
    In a flat universe q=k, in a closed universe q is quantised.

    Got get an instance of this data, call `CAMBdata.get_matter_power_spectrum`

    :ivar num_q_trans:  number of q modes calcualted
    :ivar q_trans: array of q values calculated
    :ivar sigma_8: array of sigma8 values for each redshift for each power spectrum
    :ivar sigma2_vdelta_8: array of v-delta8 correlation, so sigma2_vdelta_8/sigma_8 can define growth
    :ivar TransferData: numpy array T[entry, q_index, z_index] storing transfer functions for each redshift and q; entry+1 can be

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
    pass


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
    camblib.__camb_MOD_camb_setdefparams(ctypes.byref(P))
    return P


def fortran_array(cP, shape, dtype=np.float64):
    # this often seems to make a copy anway, so enforce that for consistency
    # for non-copy as_array().reshape() might work, but doesn't allow non-default types...
    # return nplib.as_array(cP, tuple(shape[::-1])).reshape(shape, order='F')
    arr= np.ndarray(tuple(shape), dtype, np.core.multiarray.int_asbuffer(
        ctypes.addressof(cP.contents), np.prod(shape) * np.dtype(dtype).itemsize), order='F')
    if not arr.flags.owndata:
        arr = arr.copy()
    return arr


def set_z_outputs(z_outputs):
    """
    Set the redshifts for calculating BAO parameters at
    :param z_outputs: array of redshifts
    """
    z_outputs = np.array(z_outputs)
    CAMB_SetBackgroundOutputs_z(z_outputs, ctypes.byref(ctypes.c_int(len(z_outputs))))


class CAMBdata(object):
    def __init__(self):
        """
        Makes an object for storing transfer function data and parameters for CAMB.
        Not that it *only* stores transfer functions. If you want to get power spectra or background functions,
        you must have called one of the calculation functions for the parametes of interest more recently than
        any other call to these functions. You can can make multiple instances of CAMBdata and then later call
        `power_spectra_from_transfer` to calculate other quantities.

        To quickly make a fully calculated CAMBdata instance for a set of parameters you can call `get_results`.

        """
        self._key = None
        self._key = CAMBdata_new()
        self.Params = self.get_params()
        self._one = ctypes.c_int(1)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def __del__(self):
        self.free()

    def free(self):
        if self._key:
            CAMBdata_free(ctypes.byref(self._key))
            self._key = None

    def set_params(self, params):
        """
        Set parameters from params

        :param params: a :class:`.model.CAMBparams` instance
        """
        CAMBdata_setparams(self._key, ctypes.byref(params))

    def get_params(self):
        """
        :return: :class:`.model.CAMBparams` instance pointing to the underlying parameters used by CAMB.
        Elemenets can be modified without calling set_params again.
        """
        return CAMBdata_getparams(self._key).contents

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
        Get BAO values for redshifts (set redshifts using 'set_z_outputs')

        :return: rs/DV, H, DA, F_AP for each requested redshift (as 2D array)
        """
        n = CAMB_GetNumBackgroundOutputs()
        if not n:
            raise CAMBError(
                'Call camb.set_z_outputs with required redshifts (and then calculate transfers/results) before calling get_background_outputs')
        outputs = np.empty((n, 4))
        CAMB_GetBackgroundOutputs(outputs, ctypes.byref(ctypes.c_int(n)))
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
        e.g. call this if you want to just use `angular_diameter_distance` and similar background functions
        :param params:  :class:`.model.CAMBparams` instance to use
        """
        CAMB_SetParamsForBackground(self._key, ctypes.byref(params))

    def calc_background(self, params):
        """
        Calculate the background evolution and thermal history.
        e.g. call this if you want to get derived parameters and call background functions
        :param params:  :class:`.model.CAMBparams` instance to use
        """
        res = CAMB_CalcBackgroundTheory(self._key, ctypes.byref(params))
        if res:
            raise CAMBError('Error %s in calc_background' % res)

    def calc_transfers(self, params, only_transfers=True):
        """
        Calculate the transfer functions (for CMB and matter power, as determined by params.WantCls, params.WantTransfer)

        :param params: :class:`.model.CAMBparams` instance with parameters to use
        :param only_transfers: only calcualte transfer functions, no power spectra
        :return: non-zero if error, zero if OK
        """
        opt = ctypes.c_bool()
        opt.value = only_transfers
        return CAMBdata_gettransfers(self._key, ctypes.byref(params), ctypes.byref(opt))

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
            CAMBdata_transferstopowers(self._key)

    def power_spectra_from_transfer(self, initial_power_params):
        """
        Assuming `calc_transfers` or `calc_power_spectra` have already been used, re-calcualte the power spectra
        using a new set of initial power spectrum parameters with otherwise the same cosmologiy.
        This is typically much faster that re-calculating everything, as the transfer functions can be re-used.

        :param initial_power_params: :class:`initialpower.InitialPowerParams` instance with new primordial power spectrum parameters
        """
        if initial_power_params.has_tensors() and not self.Params.WantTensors:
            raise CAMBError('r>0 but params.WantTensors = F')
        self.get_params().set_initial_power(initial_power_params)
        CAMBdata_transferstopowers(self._key)

    def get_cmb_power_spectra(self, params=None, lmax=None,
                              spectra=['total', 'unlensed_scalar', 'lensed_scalar', 'tensor', 'lens_potential']):
        """
        return a dictionary of CMB power spectra, as requested by the 'spectra' argument. All power spectra are self-owned
        numpy arrays (0..lmax, 0..3), where 0..3 index are TT, EE, BB TT.

        :param params: optional :class:`.model.CAMBparams` instance with parameters to use. If None, must have
        previously set parameters and called `calc_power_spectra` (e.g. if you got this instance using `camb.get_results`),
        :param lmax: maximum l
        :param spectra:
        :return:
        """
        P = {}
        if params is not None:
            self.calc_power_spectra(params)
        if self.Params.InitPower.has_tensors() and not self.Params.WantTensors:
            raise CAMBError('r>0 but params.WantTensors = F')

        if self.Params.DoLensing:
            lmax_calc = model.lmax_lensed.value
        else:
            lmax_calc = self.Params.Max_l
        if lmax is None:
            lmax = lmax_calc
        elif lmax > lmax_calc:
            logging.warning('getting CMB power spectra to higher L than calculated, may be innacurate/zeroed.')
        for spectrum in spectra:
            P[spectrum] = getattr(self, 'get_' + spectrum + '_cls')(lmax)
        return P

    def get_matter_transfer_data(self):
        """
        Get matter transfer function data and sigma8 for calculated results.

        :return: :class:`.MatterTransferData` instance holding output arrays (copies, not pointers)
        """
        if not self.Params.WantTransfer:
            raise CAMBError("must have Params.WantTransfer to get matter transfers and power")

        cdata = _MatterTransferData()
        CAMBdata_mattertransferdata(self._key, ctypes.byref(cdata))
        data = MatterTransferData()
        data.num_q_trans = cdata.num_q_trans
        data.q_trans = nplib.as_array(cdata.q_trans, shape=(data.num_q_trans,))
        data.sigma_8 = fortran_array(cdata.sigma_8, cdata.sigma_8_size)
        data.sigma2_vdelta_8 = fortran_array(cdata.sigma2_vdelta_8, cdata.sigma2_vdelta_8_size)
        data.TransferData = fortran_array(cdata.TransferData, cdata.TransferData_size, dtype=np.float32)
        return data

    def get_linear_matter_power_spectrum(self, var1=model.transfer_power_var.value, var2=model.transfer_power_var.value,
                                         hubble_units=True, have_power_spectra=False, params=None):
        """
        Calculates P_{xy}(k/h), where x, y are one of model.Transfer_cdm, model.Transfer_xx etc defined in ModelParams.
        The output k values are not regularly spaced, and not interpolated.

        :param var1: variable i
        :param var2: variable j
        :param hubble_units: if true, output power spectrum in (Mpc/h)^{-3} units, otherwise Mpc^{-3}
        :param have_power_spectra: set to True if already computed power spectra
        :param params: if have_power_spectra=False, optional :class:`.model.CAMBparams` instance to specify new parameters
        :return: kh, z, PK, where kz an z are arrays of k/h and z respectively, and PK[i,j] is value at z[i], k/h[j]
        """
        if not have_power_spectra:
            self.calc_power_spectra(params)
        data = self.get_matter_transfer_data()

        nk = data.num_q_trans
        nz = self.Params.Transfer.PK_num_redshifts
        kh = data.TransferData[model.Transfer_kh - 1, :, 0]

        var1 = ctypes.c_int(var1)
        var2 = ctypes.c_int(var2)
        hubble_units = ctypes.c_int(hubble_units)
        PK = np.empty((nz, nk))
        CAMBdata_GetLinearMatterPower(self._key, PK, ctypes.byref(var1), ctypes.byref(var2), ctypes.byref(hubble_units))

        z = self.Params.Transfer.PK_redshifts[:nz]
        z.reverse()
        return np.array(kh), np.array(z), PK

    def get_sigma8(self):
        """
        Get sigma_8 values (must previous have calculated power spectra)

        :return: array of sigma_8 values, in order of increasing time (decreasing redshift)
        """
        mtrans = self.get_matter_transfer_data()
        return mtrans.sigma_8[:, 0]

    def get_matter_power_spectrum(self, minkh=1e-4, maxkh=1.0, npoints=100,
                                  var1=model.transfer_power_var.value, var2=model.transfer_power_var.value,
                                  have_power_spectra=False, params=None):
        """
        Calculates P_{xy}(k/h), where x, y are one of Transfer_cdm, Transfer_xx etc defined in ModelParams.
        The output k values are regularly log spaced and interpolated. If NonLinear is set, the result is non-linear.
        :param minkh: minimum value of k/h for output grid (very low values < 1e-4 may not be calculated)
        :param maxkh: maximum value of k/h (check consistent with input params.Transfer.kmax)
        :param npoints: number of points equally spaced in log k
        :param var1: variable i
        :param var2: variable j
        :param have_power_spectra: set to True if already computed power spectra
        :param params: if have_power_spectra=False and want to specify new parameters, a :class:`.model.CAMBparams` instance
        :return: kh, z, PK, where kz an z are arrays of k/h and z respectively, and PK[i,j] is value at z[i], k/h[j]
        """

        if not have_power_spectra:
            self.calc_power_spectra(params)

        assert (self.Params.WantTransfer)
        if self.Params.Transfer.kmax < maxkh:
            logging.warning("get_matter_power_spectrum using larger k_max than input parameter Transfer.kmax")
        if self.Params.NonLinear == model.NonLinear_none and self.Params.Transfer.kmax < 1:
            logging.warning("get_matter_power_spectrum Transfer.kmax small to get non-linear spectrum")

        nz = self.Params.Transfer.PK_num_redshifts
        PK = np.empty((nz, npoints))
        dlnkh = (np.log(maxkh) - np.log(minkh)) / (npoints - 1)
        CAMBdata_GetMatterPower(self._key, PK, ctypes.byref(ctypes.c_double(minkh)),
                                ctypes.byref(ctypes.c_double(dlnkh)),
                                ctypes.byref(ctypes.c_int(npoints)), ctypes.byref(ctypes.c_int(var1)),
                                ctypes.byref(ctypes.c_int(var2)))
        z = self.Params.Transfer.PK_redshifts[:nz]
        z.reverse()
        return minkh * np.exp(np.arange(npoints) * dlnkh), z, PK

    def get_total_cls(self, lmax):
        """
        Get lensed-scalar + tensor CMB power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE
        """
        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetTotCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_tensor_cls(self, lmax):
        """
        Get tensor CMB power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE
        """

        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetTensorCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_unlensed_scalar_cls(self, lmax):
        """
        Get unlensed scalar CMB power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE. CL[:,2] will be zero.
        """

        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetUnlensedScalCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_lensed_scalar_cls(self, lmax):
        """
        Get lensed scalar CMB power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :return: numpy array CL[0:lmax+1,0:4], where 0..3 indexes TT, EE, BB, TE.
        """

        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetLensedScalCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_lens_potential_cls(self, lmax):
        """
        Get lensing deflection angle potential power spectrum, and cross-correlation with T and E. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :return: numpy array CL[0:lmax+1,0:3], where 0..2 indexes PP, PT, PE.
        """

        res = np.empty((lmax + 1, 3))
        opt = ctypes.c_int(lmax)
        CAMB_SetLensPotentialCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_unlensed_scalar_array_cls(self, lmax):
        """
        Get array of all cross power spectra. Must have already calculated power spectra.

        :param lmax: lmax to output to
        :return: numpy array CL[0:lmax+1,0:, 0:], where 0.. index T, E, deflection angle, source window functions
        """

        if not model.has_cl_2D_array.value:
            raise CAMBError('unlensed_scalar_array not calculated')
        n = 3 + model.num_redshiftwindows.value
        res = np.empty((lmax + 1, n, n))
        opt = ctypes.c_int(lmax)
        CAMB_SetUnlensedScalarArray(ctypes.byref(opt), res, ctypes.byref(self._one), ctypes.byref(ctypes.c_int(n)))
        return res

    def angular_diameter_distance(self, z):
        """
        Get (non-comoving) angular diameter distance to redshift z.

        Must have called calc_background, calc_background_no_thermo or calculated transfer functions or power spectra.

        :param z: redshift or array of redshifts
        :return: angular diameter distances, matching rank of z
        """
        if np.isscalar(z):
            return AngularDiameterDistance(ctypes.byref(ctypes.c_double(z)))
        else:
            z=np.asarray(z)
            arr = np.empty(z.shape)
            AngularDiameterDistanceArr(arr, z, ctypes.byref(ctypes.c_int(z.shape[0])))
            return arr

    def angular_diameter_distance2(self, z1, z2):
        """
        Get angular diameter distance between two redshifts

        [r/(1+z2)]*sin_K([chi(z_1) - chi(z_1)]/r) where r is curvature radius and chi is the comoving radial distance

        Must have called calc_background, calc_background_no_thermo or calculated transfer functions or power spectra.

        :param z1: redshift 1
        :param z2: redshift 2
        :return: result
        """
        return AngularDiameterDistance2(ctypes.byref(ctypes.c_double(z1)), ctypes.byref(ctypes.c_double(z2)))

    def comoving_radial_distance(self, z):
        """
        Get comoving radial distance from us to redshift z.

        Must have called calc_background, calc_background_no_thermo or calculated transfer functions or power spectra.

        :param z: redshift
        :return: comoving radial distance
        """
        return ComovingRadialDistance(ctypes.byref(ctypes.c_double(z)))

    def luminosity_distance(self, z):
        """
        Get luminosity distance from to redshift z.

        Must have called calc_background, calc_background_no_thermo or calculated transfer functions or power spectra.

        :param z: redshift or array of redshifts
        :return: luminosity distance (matches rank of z)
        """

        return self.angular_diameter_distance(z)*(1.0+np.asarray(z))**2

    def h_of_z(self, z):
        """
        Get Hubble rate at redshift z, in Mpc^{-1} units.

        Must have called calc_background, calc_background_no_thermo or calculated transfer functions or power spectra.

        :param z: reshift
        :return: H(z)
        """
        if not np.isscalar(z):
            raise CAMBError('vector z not supported yet')
        return Hofz(ctypes.byref(ctypes.c_double(z)))

    def h_of_z_hunit(self, z):
        """
        Get Huuble rate at redshift z, in km/s/Mpc units.

        Must have called calc_background, calc_background_no_thermo or calculated transfer functions or power spectra.

        :param z: reshift
        :return: H(z)/[km/s/Mpc]
        """
        return constants.c * self.h_of_z(z) / 1e3

    def physical_time_gyr(self, a1, a2):
        """
        Get physical time between two scalar factors in Gigayears

        Must have called calc_background, calc_background_no_thermo or calculated transfer functions or power spectra.

        :param a1: scale factor 1
        :param a2: scale factor 2
        :return: (age(a2)-age(a1))/Gigayear
        """
        return DeltaPhysicalTimeGyr(ctypes.byref(ctypes.c_double(a1)), ctypes.byref(ctypes.c_double(a2)), None)

    def physical_time_z(self, z):
        """
        Get physical time from hot big bang to redshift z in Gigayears.
        :param z:  redshift
        :return: t(z)/Gigayear
        """
        if not np.isscalar(z):
            raise CAMBError('vector z not supported yet')
        return self.physical_time_gyr(0, 1.0 / (1 + z))

    def conformal_time(self, a1, a2):
        """
        Get conformal time between two scale facors (=comoving radial distance travelled by light on light cone)
        :param a1: scale factor 1
        :param a2: scale factor 2
        :return: eta(a2)-eta(a1) = chi(a1)-chi(a2) in Megaparsec
        """

        if not np.isscalar(a1) or not np.isscalar(a2):
            raise CAMBError('vector z not supported yet')
        return DeltaTime(ctypes.byref(ctypes.c_double(a1)), ctypes.byref(ctypes.c_double(a2)), None)

    def conformal_time_z(self, z):
        """
        Conformal time from hot big bang to redshift z in Megaparsec.

        :param z: redshift
        :return: eta(z)/Mpc
        """
        return self.conformal_time(0, 1.0 / (1 + z))

    def cosmomc_theta(self):
        """
        Get theta_MC, an approximation of the radio of the sound horizon to the angular diameter distance at recombination
        :return: theta_MC
        """
        return CosmomcTheta()


def get_results(params):
    """
    Calculate results for specified parameters and return :class:`CAMBdata` instance for getting results.

    :param params: :class:`.model.CAMBparams` instance
    :return: :class:`CAMBdata` instance
    """
    res = CAMBdata()
    res.calc_power_spectra(params)
    return res

def get_transfer_functions(params):
    """
    Calculate transfer functions for specified parameters and return :class:`CAMBdata` instance for getting results
    and subsequently calcualating power spectra.

    :param params: :class:`.model.CAMBparams` instance
    :return: :class:`CAMBdata` instance
    """

    res = CAMBdata()
    res.calc_transfers(params)
    return res


def get_age(params):
    """
    Get age of universe for given set of parameters

    :param params:  :class:`.model.CAMBparams` instance
    :return: age of universe in gigayears
    """
    return CAMB_GetAge(ctypes.byref(params))


def get_zre_from_tau(params, tau):
    """
    Get reionization redshift given optical depth tau

    :param params: :class:`.model.CAMBparams` instance
    :param tau: optical depth
    :return: reionization redshift
    """
    cTau = ctypes.c_double(tau)
    return CAMB_GetZreFromTau(ctypes.byref(params), ctypes.byref(cTau))


def cleanup():
    camblib.__camb_MOD_camb_cleanup()
