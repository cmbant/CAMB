from .baseconfig import CAMB_Structure
from ctypes import c_bool, c_double


class ReionizationParams(CAMB_Structure):
    """
    Holds parameters for the reionization model. The default (unphysical) model tanh parameterization is described in
    Appendix B of `arXiv:0804.3865 <http://arxiv.org/abs/0804.3865>`_
    This should become a changeable class at some point.
    """
    _fields_ = [
        ("Reionization", c_bool, "Is there reionization? (can be off for matter power which is independent of it)"),
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
         "Accuracy boosting parameter for the minimum number of time sampling steps through reionization")]

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
