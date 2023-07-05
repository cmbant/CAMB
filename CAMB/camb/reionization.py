from .baseconfig import F2003Class, fortran_class
from ctypes import c_bool, c_double, POINTER, byref, c_void_p


class ReionizationModel(F2003Class):
    """
    Abstract base class for reionization models.
    """
    _fields_ = [
        ("Reionization", c_bool, "Is there reionization? (can be off for matter power which is independent of it)")]


@fortran_class
class TanhReionization(ReionizationModel):
    """
    This default (unphysical) tanh x_e parameterization is described in
    Appendix B of `arXiv:0804.3865 <https://arxiv.org/abs/0804.3865>`_
    """
    _fields_ = [
        ("use_optical_depth", c_bool, "Whether to use the optical depth or redshift paramters"),
        ("redshift", c_double, "Reionization redshift if use_optical_depth-False"),
        ("optical_depth", c_double, "Optical depth if use_optical_depth=True"),
        ("delta_redshift", c_double, "Duration of reionization"),
        ("fraction", c_double,
         "Reionization fraction when complete, or -1 for full ionization of hydrogen and first ionization of helium."),
        ("include_helium_fullreion", c_bool, "Whether to include second reionization of helium"),
        ("helium_redshift", c_double, "Redshift for second reionization of helium"),
        ("helium_delta_redshift", c_double, "Width in redshift for second reionization of helium"),
        ("helium_redshiftstart", c_double, "Include second helium reionizatio below this redshift"),
        ("tau_solve_accuracy_boost", c_double, "Accuracy boosting parameter for solving for z_re from tau"),
        ("timestep_boost", c_double,
         "Accuracy boosting parameter for the minimum number of time sampling steps through reionization"),
        ("max_redshift", c_double, "Maxmimum redshift allowed when mapping tau into reionization redshift")]

    _fortran_class_module_ = 'Reionization'
    _fortran_class_name_ = 'TTanhReionization'

    _methods_ = [('GetZreFromTau', [c_void_p, POINTER(c_double)], c_double, {"nopass": True})]

    def set_zrei(self, zrei, delta_redshift=None):
        """
        Set the mid-point reionization redshift

        :param zrei: mid-point redshift
        :param delta_redshift:  delta z for reionization
        :return:  self
        """
        self.use_optical_depth = False
        self.redshift = zrei
        if delta_redshift is not None:
            self.delta_redshift = delta_redshift
        return self

    def set_tau(self, tau, delta_redshift=None):
        """
        Set the optical depth

        :param tau: optical depth
        :param delta_redshift: delta z for reionization
        :return: self
        """
        self.use_optical_depth = True
        self.optical_depth = tau
        if delta_redshift is not None:
            self.delta_redshift = delta_redshift
        return self

    def get_zre(self, params, tau=None):
        """
        Get the midpoint redshift of reionization.

        :param params: :class:`.model.CAMBparams` instance with cosmological parameters
        :param tau: if set, calculate the redshift for optical depth tau, otherwise uses curently set parameters
        :return: reionization mid-point redshift
        """
        if self.use_optical_depth or tau:
            from .camb import CAMBparams
            assert isinstance(params, CAMBparams)
            return self.f_GetZreFromTau(byref(params), c_double(tau or self.optical_depth))
        else:
            return self.redshift
