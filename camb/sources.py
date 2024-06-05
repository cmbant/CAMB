from .baseconfig import F2003Class, fortran_class, numpy_1d, numpy_2d, np, numpy_1d_or_null
from ctypes import POINTER, c_int, c_double, byref


class SourceWindow(F2003Class):
    """
    Abstract base class for a number count/lensing/21cm source window function.
    A list of instances of these classes can be assigned to the SourceWindows field of :class:`.model.CAMBparams`.

    Note that source windows can currently only be used in flat models.
    """
    _fields_ = [("source_type", c_int, {"names": ["21cm", "counts", "lensing"], "start": 1}),
                ("bias", c_double),
                ("dlog10Ndm", c_double)]

    _fortran_class_module_ = "SourceWindows"
    _fortran_class_name_ = "TSourceWindow"


@fortran_class
class GaussianSourceWindow(SourceWindow):
    """
    A Gaussian W(z) source window function.
    """
    _fields_ = [("redshift", c_double),
                ("sigma", c_double)]

    _fortran_class_name_ = "TGaussianSourceWindow"


@fortran_class
class SplinedSourceWindow(SourceWindow):
    """
    A numerical W(z) source window function constructed by interpolation from a numerical table.
    """
    _fortran_class_name_ = "TSplinedSourceWindow"

    _methods_ = [("SetTable", [POINTER(c_int), numpy_1d, numpy_1d, numpy_1d_or_null]),
                 ("SetTable2DBias", [POINTER(c_int), POINTER(c_int), numpy_1d, numpy_1d, numpy_1d, numpy_2d])
                 ]

    def __init__(self, **kwargs):
        z = kwargs.pop('z', None)
        if z is not None:
            self.set_table(z, kwargs.pop('W'), kwargs.pop('bias_z', None),
                           kwargs.pop('k_bias', None), kwargs.pop('bias_kz', None))
        super().__init__(**kwargs)

    def __getstate__(self):
        raise TypeError("Cannot save class with splines")

    def set_table(self, z, W, bias_z=None, k_bias=None, bias_kz=None):
        """
        Set arrays of z and W(z) for cubic spline interpolation. Note that W(z) is the total count distribution
        observed, not a fractional selection function on an underlying distribution.

        :param z: array of redshift values (monotonically increasing)
        :param W: array of window function values. It must be well enough sampled to smoothly cubic-spline interpolate
        :param bias_z: optional array of bias values at each z for scale-independent bias
        :param k_bias: optional array of k values for bias (Mpc^-1)
        :param bias_kz: optional 2D contiguous array for space-dependent bias(k, z).
                        Must ensure range of k is large enough to cover required values.

        """
        if len(W) != len(z) or z[-1] < z[1] or len(z) < 5:
            raise ValueError(
                "Redshifts must be well sampled and in ascending order, with window function the same length as z")
        z = np.ascontiguousarray(z, dtype=np.float64)
        W = np.ascontiguousarray(W, dtype=np.float64)
        if bias_z is not None:
            if bias_kz is not None:
                raise ValueError("set bias_k or bias_zk")
            bias_z = np.ascontiguousarray(bias_z, dtype=np.float64)
            if len(bias_z) != len(z):
                raise ValueError("bias array must be same size as the redshift array")
        if bias_kz is not None:
            k = np.ascontiguousarray(k_bias, dtype=np.float64)
            if bias_kz.shape[0] != len(k) or bias_kz.shape[1] != len(z):
                raise ValueError('Bias array does not match shape of k,z arrays')
            self.f_SetTable2DBias(byref(c_int(len(z))), byref(c_int(len(k))), z, k, W, bias_kz)
        else:
            self.f_SetTable(byref(c_int(len(z))), z, W, bias_z)
