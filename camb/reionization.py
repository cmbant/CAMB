from ctypes import POINTER, byref, c_bool, c_double, c_void_p

from .baseconfig import F2003Class, f_pointer, fortran_class


class ReionizationModel(F2003Class):
    """
    Abstract base class for reionization models.

    The optional ``include_heating`` switch adds an approximate model of the photo-heating of the
    intergalactic medium during reionization. It is a very crude instantaneous-Jeans approximation for
    an order of magnitude estimate. When enabled, the baryon temperature is smoothly raised
    towards ``reion_temperature`` (default :math:`10^4` K) following the ionization-fraction shape, and the
    baryon sound speed is mapped towards the corresponding ideal-gas value. This raises the baryon pressure,
    suppressing small-scale baryon (and hence total matter) power. It is off by default and only affects the
    low-redshift matter power spectrum.

    """

    _fields_ = (
        ("Reionization", c_bool, "Is there reionization? (can be off for matter power, which is independent of it)"),
        (
            "include_heating",
            c_bool,
            "Whether to include approximate smooth reionization heating following the x_e shape (default False)",
        ),
        (
            "reion_temperature",
            c_double,
            "Asymptotic gas temperature in K reached during reionization heating (default 1e4 K)",
        ),
    )

    def write_ini(self, state) -> None:
        state.set("reionization", self.Reionization)
        state.set("reion_include_heating", self.include_heating)
        state.set("reion_temperature", self.reion_temperature)


class BaseTauWithHeReionization(ReionizationModel):
    """
    Abstract class for models that map z_re to tau, and include second reionization of Helium
    """

    _fields_ = (
        ("use_optical_depth", c_bool, "Whether to use the optical depth or redshift parameters"),
        ("redshift", c_double, "Reionization redshift (xe=0.5) if use_optical_depth=False"),
        ("optical_depth", c_double, "Optical depth if use_optical_depth=True"),
        (
            "fraction",
            c_double,
            "Reionization fraction when complete, "
            "or -1 for full ionization of hydrogen and first ionization of helium.",
        ),
        ("include_helium_fullreion", c_bool, "Whether to include second reionization of helium"),
        ("helium_redshift", c_double, "Redshift for second reionization of helium"),
        ("helium_delta_redshift", c_double, "Width in redshift for second reionization of helium"),
        ("helium_redshiftstart", c_double, "Include second helium reionization below this redshift"),
        ("tau_solve_accuracy_boost", c_double, "Accuracy boosting parameter for solving for z_re from tau"),
        (
            "timestep_boost",
            c_double,
            "Accuracy boosting parameter for the minimum number of time sampling steps through reionization",
        ),
        ("max_redshift", c_double, "Maximum redshift allowed when mapping tau into reionization redshift"),
        ("__min_redshift", c_double, "Minimum redshift allowed when mapping tau into reionization redshift"),
        ("__fHe", c_double, "Helium fraction"),
        ("__state", f_pointer),
    )

    _fortran_class_module_ = "Reionization"
    _fortran_class_name_ = "TBaseTauWithHeReionization"

    _methods_ = (("GetZreFromTau", [c_void_p, POINTER(c_double)], c_double, {"nopass": True}),)

    def get_zre(self, params, tau=None):
        """
        Get the midpoint redshift of reionization.

        :param params: :class:`.model.CAMBparams` instance with cosmological parameters
        :param tau: if set, calculate the redshift for optical depth tau, otherwise uses currently set parameters
        :return: reionization mid-point redshift
        """
        if self.use_optical_depth or tau:
            from .camb import CAMBparams

            assert isinstance(params, CAMBparams)
            return self.f_GetZreFromTau(byref(params), c_double(tau or self.optical_depth))
        else:
            return self.redshift

    def set_zrei(self, zrei):
        """
        Set the mid-point reionization redshift

        :param zrei: mid-point redshift
        :return: self
        """
        self.use_optical_depth = False
        self.redshift = zrei
        return self

    def set_tau(self, tau):
        """
        Set the optical depth

        :param tau: optical depth
        :return: self
        """
        self.use_optical_depth = True
        self.optical_depth = tau
        return self

    def set_extra_params(self, max_zrei=None):
        """
        Set extra parameters (not tau, or zrei)

        :param max_zrei: maximum redshift allowed when mapping tau into reionization redshift
        """
        if max_zrei is not None:
            self.max_redshift = max_zrei

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.set("re_use_optical_depth", self.use_optical_depth)
        if self.use_optical_depth:
            state.set("re_optical_depth", self.optical_depth)
        else:
            state.set("re_redshift", self.redshift)
        state.set("re_ionization_frac", self.fraction)
        state.set("include_helium_fullreion", self.include_helium_fullreion)
        state.set("re_helium_redshift", self.helium_redshift)
        state.set("re_helium_delta_redshift", self.helium_delta_redshift)
        state.set("re_helium_redshiftstart", self.helium_redshiftstart)
        state.set("max_zrei", self.max_redshift)


@fortran_class
class TanhReionization(BaseTauWithHeReionization):
    """
    This default (unphysical) tanh x_e parameterization is described in
    Appendix B of `arXiv:0804.3865 <https://arxiv.org/abs/0804.3865>`_
    """

    _fields_ = (("delta_redshift", c_double, "Duration of reionization"),)

    _fortran_class_name_ = "TTanhReionization"

    def set_extra_params(self, deltazrei=None, max_zrei=None) -> None:
        """
        Set extra parameters (not tau, or zrei)

        :param deltazrei: delta z for reionization
        :param max_zrei: maximum redshift allowed when mapping tau into reionization
        """
        super().set_extra_params(max_zrei)
        if deltazrei is not None:
            self.delta_redshift = deltazrei

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.set("re_delta_redshift", self.delta_redshift)


@fortran_class
class ExpReionization(BaseTauWithHeReionization):
    """
    An ionization fraction that decreases exponentially at high z, saturating to fully ionized at fixed redshift.
    This model has a minimum non-zero tau around 0.04 for reion_redshift_complete=6.1.
    Similar to e.g. arXiv:1509.02785, arXiv:2006.16828, but not attempting to fit shape near x_e~1 at z<6.1
    """

    _fields_ = (
        ("reion_redshift_complete", c_double, "end of reionization"),
        ("reion_exp_smooth_width", c_double, "redshift scale to smooth exponential"),
        ("reion_exp_power", c_double, "power in exponential, exp(-lambda(z-redshift_complete)^exp_power)"),
    )

    _fortran_class_name_ = "TExpReionization"

    def set_extra_params(
        self, reion_redshift_complete=None, reion_exp_power=None, reion_exp_smooth_width=None, max_zrei=None
    ) -> None:
        """
        Set extra parameters (not tau, or zrei)

        :param reion_redshift_complete: redshift at which reionization complete (e.g. around 6)
        :param reion_exp_power: power in exponential decay with redshift
        :param reion_exp_smooth_width: smoothing parameter to keep derivative smooth
        :param max_zrei: maximum redshift allowed when mapping tau into reionization
        """
        super().set_extra_params(max_zrei)
        if reion_redshift_complete is not None:
            self.reion_redshift_complete = reion_redshift_complete
        if reion_exp_power is not None:
            self.reion_exp_power = reion_exp_power
        if reion_exp_smooth_width is not None:
            self.reion_exp_smooth_width = reion_exp_smooth_width

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.write_fields(self, names=("reion_redshift_complete", "reion_exp_power", "reion_exp_smooth_width"))
