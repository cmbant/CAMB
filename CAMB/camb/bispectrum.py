from ctypes import c_int, c_double, c_char

# Note currently there is no python interface for getting bispectrum results.

Ini_max_string_len = 1024
max_bispectrum_deltas = 5


class TBispectrumParams:
    _fields_ = [
        ("do_lensing_bispectrum", c_int),  # logical
        ("do_primordial_bispectrum", c_int),  # logical
        ("nfields", c_int),
        ("Slice_Base_L", c_int),
        ("deltas", c_int * max_bispectrum_deltas),
        ("do_parity_odd", c_int),  # logical
        ("DoFisher", c_int),  # logical
        ("export_alpha_beta", c_int),  # logical
        ("FisherNoise", c_double),
        ("FisherNoisePol", c_double),
        ("FisherNoiseFwhmArcmin", c_double),
        ("FullOutputFile", c_char * Ini_max_string_len),
        ("SparseFullOutput", c_int),  # logical
    ]
