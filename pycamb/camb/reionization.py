from baseconfig import camblib, CAMB_Structure
from ctypes import c_bool, c_int, c_double

# ---Variables in reionization.f90
# To set the value please just put 
# variablename.value = newvalue

# logical
include_helium_fullreion = c_bool.in_dll(camblib, "__reionization_MOD_include_helium_fullreion")
# include_helium_fullreion.value = True

# logical
Reionization_AccuracyBoost = c_bool.in_dll(camblib, "__reionization_MOD_reionization_accuracyboost")
# Reionization_AccuracyBoost.value = 1.

Rionization_zexp = c_bool.in_dll(camblib, "__reionization_MOD_rionization_zexp")

# ---Derived Types in reionization.f90

class ReionizationParams(CAMB_Structure):
    """
    Hold sparameters for the reionization model.
    """
    _fields_ = [
        ("Reionization", c_int),  # logical
        ("use_optical_depth", c_int),  # logical
        ("redshift", c_double),
        ("delta_redshift", c_double),
        ("fraction", c_double),
        ("optical_depth", c_double),
        ("helium_redshift", c_double),  # helium_redshift  = 3.5_dl
        ("helium_delta_redshift", c_double),  # helium_delta_redshift  = 0.5
        ("helium_redshiftstart", c_double)  # helium_redshiftstart  = 5._dl
    ]

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


class ReionizationHistory(CAMB_Structure):
    """
    Internally calculated parameters.
    """
    _fields_ = [
        ("tau_start", c_double),
        ("tau_complete", c_double),
        ("akthom", c_double),
        ("fHe", c_double),
        ("WindowVarMid", c_double),
        ("WindowVarDelta", c_double)
    ]
