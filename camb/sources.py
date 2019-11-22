from .baseconfig import F2003Class, fortran_class, numpy_1d, np
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

    _methods_ = [("SetTable", [POINTER(c_int), numpy_1d, numpy_1d])]

    def __init__(self, **kwargs):
        z = kwargs.pop('z', None)
        if z is not None:
            self.set_table(z, kwargs.pop('W'))
        super(F2003Class, self).__init__(**kwargs)

    def set_table(self, z, W):
        """
        Set arrays of z and W(z) for cublic spline interpolation. Note that W(z) is the total count distribution
        observed, not a fractional selection function on an underlying distribution.

        :param z: array of redshift values (monotonically increasing)
        :param W: array of window function values. It must be well enough sampled to smoothly cubic-spline interpolate
        """
        if len(W) != len(z) or z[-1] < z[1] or len(z) < 5:
            raise ValueError(
                "Redshifts must be well sampled and in ascending order, with window function the same length as z")
        self.f_SetTable(byref(c_int(len(z))), np.asarray(z, dtype=np.float64), np.asarray(W, dtype=np.float64))
