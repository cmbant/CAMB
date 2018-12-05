from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np
from ctypes import c_int, c_double, byref, POINTER, c_bool


class DarkEnergyModel(F2003Class):
    """
    Abstract base class for dark energy model implementations.
    """
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int)]


class DarkEnergyEqnOfState(DarkEnergyModel):
    """
    Abstract base class for models using w and wa parameterization with use w(a) = w + (1-a)*wa parameterization,
    or call set_w_a_table to set another tabulated w(a). If tabulated w(a) is used, w and wa are set
    to approximate values at z=0.

    See :meth:`.model.CAMBparams.set_initial_power_function` for a convenience constructor function to
    set a general interpolated P(k) model from a python function.

    """
    _fortran_class_module_ = 'DarkEnergyInterface'
    _fortran_class_name_ = 'TDarkEnergyEqnOfState'

    _fields_ = [
        ("w", c_double),
        ("wa", c_double),
        ("cs2", c_double),
        ("use_tabulated_w", c_bool),
        ("no_perturbations", c_bool)
    ]

    def set_params(self, w=-1.0, wa=0, cs2=1.0):
        """
         Set the parameters so that P(a)/rho(a) = w(a) = w + (1-a)*wa

        :param w: w(0)
        :param wa: -dw/da(0)
        :param cs2: fluid rest-frame sound speed
        """
        self.w = w
        self.wa = wa
        self.cs2 = cs2

    def set_w_a_table(self, a, w):
        """
        Set w(a) from numerical values (used as cublic spline). Note this is quite slow.

        :param a: array of scale factors
        :param w: array of w(a)
        :return: self
        """
        if len(a) != len(w): raise ValueError('Dark energy w(a) table non-equal sized arrays')
        if not np.isclose(a[-1], 1):  raise ValueError('Dark energy w(a) arrays must end at a=1')

        self.call_method('setwtable', extra_args=[numpy_1d, numpy_1d, POINTER(c_int)],
                         args=[a, w, byref(c_int(len(a)))])
        return self


@fortran_class
class DarkEnergyFluid(DarkEnergyEqnOfState):
    """
    Class implementing the w, wa or splined w(a) parameterization using the constant sound-speed signle fluid model (as for single-field quintessense).

    """

    _fortran_class_module_ = 'DarkEnergyFluid'
    _fortran_class_name_ = 'TDarkEnergyFluid'

    def set_params(self, w=-1.0, wa=0, cs2=1.0):
        """
        Set the parameters so that P(a)/rho(a) = w(a) = w + (1-a)*wa.
        The fluid model must have w(a) <= -1 at all times.

        :param w: w(0)
        :param wa: -dw/da(0)
        :param cs2: fluid rest-frame sound speed
        """

        if wa and (w < -1 - 1e-6 or 1 + w + wa < -1 - 1e-6):
            raise CAMBError('fluid dark energy model does not support w crossing -1')
        self.w = w
        self.wa = wa
        self.cs2 = cs2


@fortran_class
class DarkEnergyPPF(DarkEnergyEqnOfState):
    """
    Class implementating the w, wa or splined w(a) parameterization in the PPF perturbation approximation (`arXiv:0808.3125 <http://arxiv.org/abs/0808.3125>`_)
    Use inherited methods to set parameters or interpolation table.

    """
    _fields_ = [("c_Gamma_ppf", c_double)]
    _fortran_class_module_ = 'DarkEnergyPPF'
    _fortran_class_name_ = 'TDarkEnergyPPF'


@fortran_class
class AxionEffectiveFluid(DarkEnergyModel):
    """
    Example implementation of a specifc (early) dark energy fluid model (`arXiv:1806.10608 <http://arxiv.org/abs/1806.10608>`_).
    Not well tested, but should serve to demonstrate how to make your own custom classes.
    """
    _fields_ = [
        ("w_n", c_double),
        ("om", c_double),
        ("a_c", c_double),
        ("theta_i", c_double)]
    _fortran_class_name_ = 'TAxionEffectiveFluid'
    _fortran_class_module_ = 'DarkEnergyFluid'

    def set_params(self, w_n, om, a_c):
        self.w_n = w_n
        self.om = om
        self.a_c = a_c


# short names for models that support w/wa
F2003Class._class_names.update({'fluid': DarkEnergyFluid, 'ppf': DarkEnergyPPF})
