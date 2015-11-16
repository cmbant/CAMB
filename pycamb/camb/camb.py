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


class MatterTransferData(object):
    pass


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
    pass


def set_feedback_level(level=1):
    FeedbackLevel.value = level


def set_default_params(P):
    assert (isinstance(P, model.CAMBparams))
    camblib.__camb_MOD_camb_setdefparams(ctypes.byref(P))


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
    z_outputs = np.array(z_outputs)
    CAMB_SetBackgroundOutputs_z(z_outputs, ctypes.byref(ctypes.c_int(len(z_outputs))))


class CAMBdata(object):
    def __init__(self):
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

    def set_params(self, Params):
        CAMBdata_setparams(self._key, ctypes.byref(Params))

    def get_params(self):
        """
        :return: a model.CAMBparams structure, pointing to the underlying component of CAMBdata
        """
        return CAMBdata_getparams(self._key).contents

    def get_derived_params(self):
        """
        :return: dictionary of derived parameter values
        """
        res = {}
        for name, value in zip(model.derived_names, model.ThermoDerivedParams):
            res[name] = value
        return res

    def get_background_outputs(self):
        """
        Get BAO values for redshifts (set using camb.set_z_outputs)

        :return: rs/DV, H, DA, F_AP for each requested redshift (as 2D array)
        """
        n = CAMB_GetNumBackgroundOutputs()
        if not n:
            raise CAMBError(
                'Call camb.set_z_outputs with required redshifts (and then calculate transfers/results) before calling get_background_outputs')
        outputs = np.empty((n, 4))
        CAMB_GetBackgroundOutputs(outputs, ctypes.byref(ctypes.c_int(n)))
        return outputs

    def get_BAO(self, redshifts, Params):
        """
        Get BAO parameters at given redshifts, using Params

        :param redshifts: list of redshifts
        :param Params: optional model.CAMBparams
        :return: array of rs/DV, H, DA, F_AP for each redshift as 2D array
        """
        set_z_outputs(redshifts)
        self.calc_background(Params)
        res = self.get_background_outputs()
        set_z_outputs([])
        return res

    def calc_background_no_thermo(self, Params):
        CAMB_SetParamsForBackground(self._key, ctypes.byref(Params))

    def calc_background(self, Params):
        res = CAMB_CalcBackgroundTheory(self._key, ctypes.byref(Params))
        if res:
            raise CAMBError('Error %s in calc_background' % res)

    def calc_transfers(self, Params, only_transfers=True):
        opt = ctypes.c_bool()
        opt.value = only_transfers
        return CAMBdata_gettransfers(self._key, ctypes.byref(Params), ctypes.byref(opt))

    def calc_power_spectra(self, Params=None):
        if Params is not None:
            result = self.calc_transfers(Params, only_transfers=False)
            if result != 0:
                raise CAMBError('Error getting transfer functions: %u' % result)
        else:
            CAMBdata_transferstopowers(self._key)

    def power_spectra_from_transfer(self, initial_power_params):
        if initial_power_params.has_tensors() and not self.Params.WantTensors:
            raise CAMBError('r>0 but Params.WantTensors = F')
        self.get_params().set_initial_power(initial_power_params)
        CAMBdata_transferstopowers(self._key)

    def get_cmb_power_spectra(self, Params=None, lmax=None,
                              spectra=['total', 'unlensed_scalar', 'lensed_scalar', 'tensor', 'lens_potential']):
        P = {}
        if Params is not None:
            self.calc_power_spectra(Params)
        if self.Params.InitPower.has_tensors() and not self.Params.WantTensors:
            raise CAMBError('r>0 but Params.WantTensors = F')

        if lmax is None:
            if self.Params.DoLensing:
                lmax = model.lmax_lensed.value
            else:
                lmax = self.Params.Max_l
        for spectrum in spectra:
            P[spectrum] = getattr(self, 'get_' + spectrum + '_cls')(lmax)
        return P

    def get_matter_transfer_data(self):
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
                                         hubble_units=True, have_power_spectra=False, Params=None):
        """
        Calculates P_{xy}(k/h), where x, y are one of Transfer_cdm, Transfer_xx etc defined in ModelParams.
        The output k values are not regularly spaced, and not interpolated.

        :param var1: variable i
        :param var2: variable j
        :param hubble_units: if true, output power spectrum in (Mpc/h)^{-3} units, otherwise Mpc^{-3}
        :param have_power_spectra: if already computed power spectra
        :param Params: if have_power_spectra=False and want to specify new parameters
        :return: kh, z, PK, where kz an z are arrays of k/h and z respectively, and PK[i,j] is value at z[i], k/h[j]
        """
        if not have_power_spectra:
            self.calc_power_spectra(Params)
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
        mtrans = self.get_matter_transfer_data()
        return np.array(mtrans.sigma_8[:, 0])

    def get_matter_power_spectrum(self, minkh=1e-4, maxkh=1.0, npoints=100,
                                  var1=model.transfer_power_var.value, var2=model.transfer_power_var.value,
                                  have_power_spectra=False, Params=None):
        """
        Calculates P_{xy}(k/h), where x, y are one of Transfer_cdm, Transfer_xx etc defined in ModelParams.
        The output k values are regularly log spaced and interpolated. If NonLinear is set, the result is non-linear.
        :param minkh: minimum value of k/h for output grid (very low values < 1e-4 may not be calculated)
        :param maxkh: maximum value of k/h (check consistent with input Params.Transfer.kmax)
        :param npoints: number of points equally spaced in log k
        :param var1: variable i
        :param var2: variable j
        :param have_power_spectra: if already computed power spectra
        :param Params: if have_power_spectra=False and want to specify new parameters
        :return: kh, z, PK, where kz an z are arrays of k/h and z respectively, and PK[i,j] is value at z[i], k/h[j]
        """

        if not have_power_spectra:
            self.calc_power_spectra(Params)

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
        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetTotCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_tensor_cls(self, lmax):
        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetTensorCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_unlensed_scalar_cls(self, lmax):
        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetUnlensedScalCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_lensed_scalar_cls(self, lmax):
        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetLensedScalCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_tensor_cls(self, lmax):
        res = np.empty((lmax + 1, 4))
        opt = ctypes.c_int(lmax)
        CAMB_SetTensorCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_lens_potential_cls(self, lmax):
        res = np.empty((lmax + 1, 3))
        opt = ctypes.c_int(lmax)
        CAMB_SetLensPotentialCls(ctypes.byref(opt), res, ctypes.byref(self._one))
        return res

    def get_unlensed_scalar_array_cls(self, lmax):
        if not model.has_cl_2D_array.value:
            raise CAMBError('unlensed_scalar_array not calculated')
        n = 3 + model.num_redshiftwindows.value
        res = np.empty((lmax + 1, n, n))
        opt = ctypes.c_int(lmax)
        CAMB_SetUnlensedScalarArray(ctypes.byref(opt), res, ctypes.byref(self._one), ctypes.byref(ctypes.c_int(n)))
        return res

    def angular_diameter_distance(self, z):
        if np.isscalar(z):
            return AngularDiameterDistance(ctypes.byref(ctypes.c_double(z)))
        else:
            z=np.asarray(z)
            arr = np.empty(z.shape)
            AngularDiameterDistanceArr(arr, z, ctypes.byref(ctypes.c_int(z.shape[0])))
            return arr

    def angular_diameter_distance2(self, z1, z2):
        return AngularDiameterDistance2(ctypes.byref(ctypes.c_double(z1)), ctypes.byref(ctypes.c_double(z2)))

    def comoving_radial_distance(self, z):
        return ComovingRadialDistance(ctypes.byref(ctypes.c_double(z)))

    def luminosity_distance(self, z):
        return self.angular_diameter_distance(z)*(1.0+np.asarray(z))**2

    def h_of_z(self, z):
        if not np.isscalar(z):
            raise CAMBError('vector z not supported yet')
        return Hofz(ctypes.byref(ctypes.c_double(z)))

    def h_of_z_hunit(self, z):
        return constants.c * self.h_of_z(z) / 1e3

    def physical_time_gyr(self, a1, a2):
        return DeltaPhysicalTimeGyr(ctypes.byref(ctypes.c_double(a1)), ctypes.byref(ctypes.c_double(a2)), None)

    def physical_time_z(self, z):
        if not np.isscalar(z):
            raise CAMBError('vector z not supported yet')
        return self.physical_time_gyr(0, 1.0 / (1 + z))

    def conformal_time(self, a1, a2):
        if not np.isscalar(a1) or not np.isscalar(a2):
            raise CAMBError('vector z not supported yet')
        return DeltaTime(ctypes.byref(ctypes.c_double(a1)), ctypes.byref(ctypes.c_double(a2)), None)

    def conformal_time_z(self, z):
        return self.conformal_time(0, 1.0 / (1 + z))

    def cosmomc_theta(self):
        return CosmomcTheta()


def get_results(Params):
    res = CAMBdata()
    res.calc_power_spectra(Params)
    return res

def get_transfer_functions(Params):
    res = CAMBdata()
    res.calc_transfers(Params)
    return res

def get_age(P):
    return CAMB_GetAge(ctypes.byref(P))


def get_zre_from_tau(P, tau):
    cTau = ctypes.c_double(tau)
    return CAMB_GetZreFromTau(ctypes.byref(P), ctypes.byref(cTau))


def cleanup():
    camblib.__camb_MOD_camb_cleanup()
