from .baseconfig import F2003Class, fortran_class, f_pointer
from ctypes import c_bool, c_double, POINTER, byref, c_void_p


class ReionizationModel(F2003Class):
    """
    Abstract base class for reionization models.
    """

    _fields_ = [
        ('Reionization', c_bool, 'Is there reionization? (can be off for matter power, which is independent of it)')
    ]


class BaseTauWithHeReionization(ReionizationModel):
    """
    Abstract class for models that map z_re to tau, and include second reionization of Helium
    """

    _fields_ = [
        ('use_optical_depth', c_bool, 'Whether to use the optical depth or redshift parameters'),
        ('redshift', c_double, 'Reionization redshift (xe=0.5) if use_optical_depth=False'),
        ('optical_depth', c_double, 'Optical depth if use_optical_depth=True'),
        (
            'fraction',
            c_double,
            'Reionization fraction when complete, or -1 for full ionization of hydrogen and first ionization of helium.',
        ),
        ('include_helium_fullreion', c_bool, 'Whether to include second reionization of helium'),
        ('helium_redshift', c_double, 'Redshift for second reionization of helium'),
        ('helium_delta_redshift', c_double, 'Width in redshift for second reionization of helium'),
        ('helium_redshiftstart', c_double, 'Include second helium reionization below this redshift'),
        ('tau_solve_accuracy_boost', c_double, 'Accuracy boosting parameter for solving for z_re from tau'),
        (
            'timestep_boost',
            c_double,
            'Accuracy boosting parameter for the minimum number of time sampling steps through reionization',
        ),
        ('max_redshift', c_double, 'Maximum redshift allowed when mapping tau into reionization redshift'),
        ('__min_redshift', c_double, 'Minimum redshift allowed when mapping tau into reionization redshift'),
        ('__fHe', c_double, 'Helium fraction'),
        ('__state', f_pointer),
    ]

    _fortran_class_module_ = 'Reionization'
    _fortran_class_name_ = 'TBaseTauWithHeReionization'

    _methods_ = [('GetZreFromTau', [c_void_p, POINTER(c_double)], c_double, {'nopass': True})]

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


@fortran_class
class TanhReionization(BaseTauWithHeReionization):
    """
    This default (unphysical) tanh x_e parameterization is described in
    Appendix B of `arXiv:0804.3865 <https://arxiv.org/abs/0804.3865>`_
    """

    _fields_ = [('delta_redshift', c_double, 'Duration of reionization')]

    _fortran_class_name_ = 'TTanhReionization'

    def set_extra_params(self, deltazrei=None, max_zrei=None) -> None:
        """
        Set extra parameters (not tau, or zrei)

        :param deltazrei: delta z for reionization
        :param max_zrei: maximum redshift allowed when mapping tau into reionization
        """
        super().set_extra_params(max_zrei)
        if deltazrei is not None:
            self.delta_redshift = deltazrei


@fortran_class
class ExpReionization(BaseTauWithHeReionization):
    """
    An ionization fraction that decreases exponentially at high z, saturating to fully ionized at fixed redshift.
    This model has a minimum non-zero tau around 0.04 for reion_redshift_complete=6.1.
    Similar to e.g. arXiv:1509.02785, arXiv:2006.16828, but not attempting to fit shape near x_e~1 at z<6.1
    """

    _fields_ = [
        ('reion_redshift_complete', c_double, 'end of reionization'),
        ('reion_exp_smooth_width', c_double, 'redshift scale to smooth exponential'),
        ('reion_exp_power', c_double, 'power in exponential, exp(-lambda(z-redshift_complete)^exp_power)'),
    ]

    _fortran_class_name_ = 'TExpReionization'

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


@fortran_class
class WeibullReionization(BaseTauWithHeReionization):
    """
    Weibull function reionization model following:
    - Trac et al. 2022, ApJ 927, 186 (arXiv:2109.10375) Eqs. 10-15
    - Cain et al. 2025, arXiv:2505.15899v1 Eq. 1

    This model uses a Weibull distribution to parameterize the reionization history
    with parameters for completion redshift, duration, and asymmetry.
    """

    _fields_ = [
        ('reion_redshift_complete', c_double, 'redshift at which reionization is complete (5% neutral)'),
        ('reion_duration', c_double, 'duration parameter (Delta z_90)'),
        ('reion_asymmetry', c_double, 'asymmetry parameter (A_z)'),
        ('__weibull_lambda', c_double, 'Weibull scale parameter (computed)'),
        ('__weibull_k', c_double, 'Weibull shape parameter (computed)'),
    ]

    _fortran_class_name_ = 'TWeibullReionization'

    def set_extra_params(
        self, reion_redshift_complete=None, reion_duration=None, reion_asymmetry=None, max_zrei=None
    ) -> None:
        """
        Set extra parameters (not tau, or zrei)

        :param reion_redshift_complete: redshift at which reionization is complete (5% neutral)
        :param reion_duration: duration parameter (Delta z_90)
        :param reion_asymmetry: asymmetry parameter (A_z)
        :param max_zrei: maximum redshift allowed when mapping tau into reionization
        """
        super().set_extra_params(max_zrei)
        if reion_redshift_complete is not None:
            self.reion_redshift_complete = reion_redshift_complete
        if reion_duration is not None:
            self.reion_duration = reion_duration
        if reion_asymmetry is not None:
            self.reion_asymmetry = reion_asymmetry
        # Update Weibull parameters when any parameter changes
        if hasattr(self, '_fortran_pointer') and self._fortran_pointer:
            self.f_SetParamsForZre()
