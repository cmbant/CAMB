from .baseconfig import import_property, c_int, c_bool, c_double

lensing_method_curv_corr = 1
lensing_method_flat_corr = 2
lensing_method_harmonic = 3


class _config(object):
    # print feedback if > 0 (note in Jupyter notebook this will appear in the terminal, not the notebook)
    FeedbackLevel = import_property(c_int, "config", "FeedbackLevel")

    # print additioanl timing a progress (when FeedbackLevel>0)
    DebugMsgs = import_property(c_bool, "config", "DebugMsgs")

    global_error_flag = import_property(c_int, "config", "global_error_flag")

    ThreadNum = import_property(c_int, "config", "threadnum")

    DoTensorNeutrinos = import_property(c_bool, "gaugeinterface", "dotensorneutrinos")

    DebugParam = import_property(c_double, "config", "debugparam")

    lensing_method = import_property(c_int, "lensing", "lensing_method")

    lensing_sanity_check_amplitude = import_property(c_double, "lensing", "lensing_sanity_check_amplitude")
    # lensing_sanity_check_amplitude.value = 1e-7 by default, will error if  (2*l+1)l(l+1)/4pi C_phi_phi > lensing_sanity_check_amplitude at L=10
    # increase to large number to prevent sanity check (but lensing requires realistic amplitude as non-linear)

    lensing_includes_tensors = import_property(c_bool, "lensing", "lensing_includes_tensors")

    transfer_power_var = import_property(c_int, "transfer", "transfer_power_var")

    def __repr__(self):
        s = ''
        for x in dir(self):
            if x[0] != '_':
                s += '%s = %s\n' % (x, getattr(self, x))
        return s

config = _config()
