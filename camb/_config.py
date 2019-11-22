from .baseconfig import import_property, CAMBError
from ctypes import c_char, c_int, c_bool, c_double

lensing_method_curv_corr = 1
lensing_method_flat_corr = 2
lensing_method_harmonic = 3


class _config(object):
    # print feedback if > 0 (note in Jupyter notebook this will appear in the terminal, not the notebook)
    FeedbackLevel = import_property(c_int, "config", "FeedbackLevel")

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
        return bytearray(self._global_error_message).decode('ascii').strip()

    def check_global_error(self, reference=''):
        code = self.global_error_flag
        if code:
            err = config.global_error_message()
            self.global_error_flag = 0
            if reference:
                reference = 'Error in Fortran called from %s:\n' % reference
            else:
                reference = ''
            if err:
                raise CAMBError(reference + '%s' % err)
            else:
                raise CAMBError(reference + 'Error code: %s' % code)

    def __repr__(self):
        s = ''
        for x in dir(self):
            if x[0] != '_':
                s += '%s = %s\n' % (x, getattr(self, x))
        return s


config = _config()
