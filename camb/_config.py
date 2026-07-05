import os
from ctypes import c_bool, c_char, c_double, c_int

from .baseconfig import CAMBError, import_property

lensing_method_curv_corr = 1
lensing_method_flat_corr = 2
lensing_method_harmonic = 3
lensing_method_curv_corr_full = 4
lensing_method_optimized = 5


class _config:
    # print feedback if > 0 (note in Jupyter notebook this will appear in the terminal, not the notebook)
    FeedbackLevel = import_property(c_int, "config", "FeedbackLevel")

    # if True, the Fortran code prints warnings (e.g. low-resolution, narrow window, parameter sanity checks)
    print_fortran_warnings = import_property(c_bool, "config", "print_fortran_warnings")

    # enable targeted accuracy improvements if > 0, AccuracyTarget being SO-like.
    AccuracyTarget = import_property(c_int, "config", "AccuracyTarget")

    enable_do_near_flat_integration = import_property(c_bool, "config", "enable_do_near_flat_integration")

    enable_near_flat_smallchi_integration = import_property(c_bool, "config", "enable_near_flat_smallchi_integration")

    enable_shifted_nu_scalar_approx = import_property(c_bool, "config", "enable_shifted_nu_scalar_approx")

    enable_olver_source_integration = import_property(c_bool, "config", "enable_olver_source_integration")

    # print additional timing and progress (when FeedbackLevel>0)
    DebugMsgs = import_property(c_bool, "config", "DebugMsgs")

    global_error_flag = import_property(c_int, "config", "global_error_flag")

    ThreadNum = import_property(c_int, "config", "threadnum")

    DoTensorNeutrinos = import_property(c_bool, "gaugeinterface", "dotensorneutrinos")

    DebugParam = import_property(c_double, "config", "debugparam")

    lensing_method = import_property(c_int, "lensing", "lensing_method")

    lensing_sanity_check_amplitude = import_property(c_double, "lensing", "lensing_sanity_check_amplitude")
    # lensing_sanity_check_amplitude.value = 1e-7 by default, will error if  (2*L+1)L(L+1)/4pi C_phi_phi > lensing_
    # sanity_check_amplitude at L=10
    # increase to large number to prevent sanity check (but lensing requires realistic amplitude as non-linear)

    lensing_includes_tensors = import_property(c_bool, "lensing", "lensing_includes_tensors")

    transfer_power_var = import_property(c_int, "transfer", "transfer_power_var")

    _global_error_message = import_property(c_char * 1024, "config", "global_error_message")

    def global_error_message(self):
        return bytearray(self._global_error_message).decode("ascii").strip()

    def check_global_error(self, reference=""):
        if code := self.global_error_flag:
            self.global_error_flag = 0
            if reference:
                reference = f"Error in Fortran called from {reference}:\n"
            else:
                reference = ""
            if err := config.global_error_message():
                raise CAMBError(reference + f"{err}")
            else:
                raise CAMBError(reference + f"Error code: {code}")

    def __repr__(self):
        s = ""
        for x in dir(self):
            if x[0] != "_":
                value = getattr(self, x)
                if not callable(value):
                    s += f"{x} = {value}\n"
        return s


config = _config()

if os.environ.get("BINDER_LAUNCH_HOST"):
    config.ThreadNum = 1  # binder is very slow with more than 1 CPU, force 1 by default
