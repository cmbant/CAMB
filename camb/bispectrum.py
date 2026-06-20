import ctypes
from ctypes import POINTER, byref, c_bool, c_char, c_double, c_int
from pathlib import Path

from ._config import config
from .baseconfig import (
    CAMB_Structure,
    CAMBParamRangeError,
    CAMBValueError,
    camblib,
    filepath_to_fortran,
    np,
)
from .model import CAMBparams
from .results import CAMBdata

ini_max_string_len = 1024
max_bispectrum_deltas = 5
max_bispectra = 2


class BispectrumParams(CAMB_Structure):
    """
    Settings for calculating local primordial and CMB lensing bispectra.

    The default settings calculate the CMB lensing bispectrum using temperature and polarization fields. To calculate
    the local primordial bispectrum instead, set ``do_lensing_bispectrum=False`` and
    ``do_primordial_bispectrum=True``. The primordial bispectrum is normalized to ``f_NL=1``.

    Full and sparse bispectrum outputs can be too large to return directly to Python, so set
    :attr:`full_output_file` and/or :attr:`Slice_Base_L` to write them from the Fortran calculation. Slice outputs
    write rows for fixed ``L1=Slice_Base_L`` and ``L3-L2`` values from :attr:`deltas`.

    The Fisher matrix calculation is optional and requires building CAMB from source with ``FISHER=Y`` and LAPACK
    linked. On a normal binary or source build without this flag, requesting ``DoFisher=True`` raises
    ``CAMBValueError`` rather than aborting the Python process.
    """

    _fields_ = (
        ("do_lensing_bispectrum", c_bool, "Calculate the CMB lensing bispectrum"),
        ("do_primordial_bispectrum", c_bool, "Calculate the local primordial bispectrum, normalized to f_NL=1"),
        ("nfields", c_int, "1 for T only, 2 for T and E"),
        ("Slice_Base_L", c_int, "Base L for slice output; zero disables slice output"),
        ("ndelta", c_int, "Number of L3-L2 offsets for slice output"),
        ("deltas", c_int * max_bispectrum_deltas, {"size": "ndelta"}, "L3-L2 offsets for slice output"),
        ("do_parity_odd", c_bool, "Calculate parity-odd lensing slices"),
        ("DoFisher", c_bool, "Calculate Fisher matrix outputs; requires a FISHER build"),
        ("export_alpha_beta", c_bool, "Export primordial radial alpha/beta functions"),
        ("FisherNoise", c_double, "Temperature white noise level in micro-Kelvin^2 radians"),
        ("FisherNoisePol", c_double, "Polarization white noise level in micro-Kelvin^2 radians"),
        ("FisherNoiseFwhmArcmin", c_double, "Gaussian beam FWHM in arcmin for Fisher noise"),
        ("_FullOutputFile", c_char * ini_max_string_len, "Base filename for full reduced-bispectrum output"),
        ("SparseFullOutput", c_bool, "Write only sampled/sparse full-bispectrum rows"),
    )

    def __init__(self, **kwargs):
        super().__init__()
        self.do_lensing_bispectrum = True
        self.do_primordial_bispectrum = False
        self.nfields = 2
        self.Slice_Base_L = 0
        self.deltas = []
        self.do_parity_odd = False
        self.DoFisher = False
        self.export_alpha_beta = False
        self.FisherNoise = 0
        self.FisherNoisePol = 0
        self.FisherNoiseFwhmArcmin = 7
        self.full_output_file = ""
        self.SparseFullOutput = False
        for name, value in kwargs.items():
            if name in {"full_output_file", "FullOutputFile"}:
                self.full_output_file = value
            elif name == "deltas":
                self.deltas = value
            elif name in self.get_valid_field_names():
                setattr(self, name, value)
            else:
                raise CAMBValueError(f"Unknown BispectrumParams argument: {name}")

    @property
    def full_output_file(self) -> str:
        return bytes(self._FullOutputFile).split(b"\0", 1)[0].decode("utf-8").strip()

    @full_output_file.setter
    def full_output_file(self, value) -> None:
        encoded = str(value or "").encode("utf-8")
        if len(encoded) > ini_max_string_len:
            raise CAMBParamRangeError(f"full_output_file must be at most {ini_max_string_len} bytes")
        self._FullOutputFile = encoded + b" " * (ini_max_string_len - len(encoded))

    FullOutputFile = full_output_file

    def validate(self) -> None:
        if self.nfields not in (1, 2):
            raise CAMBValueError("BispectrumParams.nfields must be 1 for T only or 2 for T and E")
        if not (self.do_lensing_bispectrum or self.do_primordial_bispectrum):
            raise CAMBValueError("At least one bispectrum type must be enabled")
        if self.ndelta > max_bispectrum_deltas:
            raise CAMBParamRangeError(f"At most {max_bispectrum_deltas} bispectrum deltas are supported")
        if self.Slice_Base_L > 0 and self.ndelta == 0:
            raise CAMBValueError("Slice_Base_L requires at least one delta")
        if (
            self.Slice_Base_L > 0
            and not self.do_parity_odd
            and any((self.Slice_Base_L + delta) % 2 for delta in self.deltas)
        ):
            raise CAMBValueError("Even-parity slice output requires Slice_Base_L + delta to be even")
        if self.do_parity_odd and (not self.do_lensing_bispectrum or self.nfields == 1):
            raise CAMBValueError("do_parity_odd requires lensing bispectrum with polarization")

    def bispectrum_names(self) -> list[str]:
        names = []
        if self.do_primordial_bispectrum:
            names.append("fnl")
        if self.do_lensing_bispectrum:
            names.append("lensing")
        return names

    def expected_output_files(self, output_root="") -> list[Path]:
        """
        Return the normal slice/full-output filenames requested by these settings.

        This helper mirrors the public large-output options and does not attempt to predict internal Fortran
        diagnostic file-tag variants.
        """
        root = str(output_root or "")
        files = []
        if self.Slice_Base_L > 0:
            for name in self.bispectrum_names():
                for delta in self.deltas:
                    if name != "lensing" and (self.Slice_Base_L + delta) % 2:
                        continue
                    files.append(Path(f"{root}bispectrum_{name}_base_{self.Slice_Base_L}_delta_{delta}.dat"))
        if self.full_output_file:
            for name in self.bispectrum_names():
                files.append(Path(f"{root}{self.full_output_file}_{name}.dat"))
        return files


class BispectrumResult(CAMB_Structure):
    """
    Small in-memory summary returned by :func:`get_bispectrum`.

    Large bispectrum tables are written directly to files requested by :class:`BispectrumParams`. This object contains
    Fisher matrices and corresponding one-sigma errors when the Fisher calculation was requested and the library was
    built with ``FISHER=Y``. Use ``has_fisher`` and ``has_optimal_fisher`` to check which arrays were filled.
    """

    _fields_ = (
        ("nbispectra", c_int),
        ("nfields", c_int),
        ("has_fisher", c_bool),
        ("has_optimal_fisher", c_bool),
        ("has_lensing_variance", c_bool),
        ("Fisher", c_double * (max_bispectra * max_bispectra)),
        ("OptimalFisher", c_double * (max_bispectra * max_bispectra)),
        ("Sigma", c_double * max_bispectra),
        ("OptimalSigma", c_double * max_bispectra),
        ("LensingFisherWithVariance", c_double),
    )

    def _matrix(self, values) -> np.ndarray:
        return (
            np.ctypeslib.as_array(values)
            .reshape((max_bispectra, max_bispectra), order="F")[: self.nbispectra, : self.nbispectra]
            .copy()
        )

    @property
    def fisher(self) -> np.ndarray:
        return self._matrix(self.Fisher)

    @property
    def optimal_fisher(self) -> np.ndarray:
        return self._matrix(self.OptimalFisher)

    @property
    def sigma(self) -> np.ndarray:
        return np.ctypeslib.as_array(self.Sigma)[: self.nbispectra].copy()

    @property
    def optimal_sigma(self) -> np.ndarray:
        return np.ctypeslib.as_array(self.OptimalSigma)[: self.nbispectra].copy()


CAMBdata_getbispectrum = camblib.__handles_MOD_cambdata_getbispectrum
CAMBdata_getbispectrum.argtypes = [
    POINTER(CAMBdata),
    POINTER(CAMBparams),
    POINTER(BispectrumParams),
    POINTER(BispectrumResult),
    ctypes.c_char_p,
    ctypes.c_long,
]
CAMBdata_getbispectrum.restype = c_int
CAMB_bispectrumfishercompiled = camblib.__handles_MOD_camb_bispectrumfishercompiled
CAMB_bispectrumfishercompiled.restype = c_int


def get_bispectrum(params, bispectrum_params=None, output_root="", _debug_params=False) -> BispectrumResult:
    """
    Calculate local primordial and/or CMB lensing bispectrum outputs.

    The wrapper runs the Fortran bispectrum calculation using a normal :class:`.model.CAMBparams` instance, without
    going through ``.ini`` files. The default :class:`BispectrumParams` calculates the CMB lensing bispectrum, so
    ``params.DoLensing`` must be true. For requested large outputs, pass an ``output_root`` prefix and set
    ``bispectrum_params.full_output_file`` or ``bispectrum_params.Slice_Base_L``.

    For example::

        from camb import bispectrum

        pars = camb.set_params(lmax=600, lens_potential_accuracy=1)
        bpars = bispectrum.BispectrumParams(Slice_Base_L=10, deltas=[0, 2])
        result = bispectrum.get_bispectrum(pars, bpars, output_root="run1_")
        print(bpars.expected_output_files("run1_"))

    Fisher outputs require a source build with the make variable ``FISHER=Y`` and LAPACK available, for example
    ``cd fortran && make python FISHER=Y`` for gfortran builds using the default ``-lblas -llapack`` link flags.

    :param params: :class:`.model.CAMBparams` instance, or a dict passed to :func:`camb.set_params`
    :param bispectrum_params: optional :class:`.BispectrumParams` instance or dict
    :param output_root: filename prefix for requested large file outputs
    :return: :class:`.BispectrumResult`
    """
    if isinstance(params, dict):
        from .camb import set_params

        params = set_params(**params)
    if bispectrum_params is None:
        bispectrum_params = BispectrumParams()
    elif isinstance(bispectrum_params, dict):
        bispectrum_params = BispectrumParams(**bispectrum_params)
    elif not isinstance(bispectrum_params, BispectrumParams):
        raise CAMBValueError("bispectrum_params must be a BispectrumParams instance or dict")
    if not isinstance(params, CAMBparams):
        raise CAMBValueError("Must pass a CAMBparams instance")
    if not params.ombh2:
        raise CAMBValueError("Parameter values not set")
    bispectrum_params.validate()
    if not (params.WantCls and params.WantScalars):
        raise CAMBValueError("get_bispectrum requires params.WantCls and params.WantScalars")
    if bispectrum_params.do_lensing_bispectrum and not params.DoLensing:
        raise CAMBValueError("do_lensing_bispectrum requires params.DoLensing = True")
    if bispectrum_params.DoFisher and not CAMB_bispectrumfishercompiled():
        raise CAMBValueError("Bispectrum Fisher output requires CAMB to be built with FISHER=Y")
    if _debug_params:
        print(params)
    res = CAMBdata()
    result = BispectrumResult()
    output_root_buffer, output_root_len = filepath_to_fortran(output_root or "")
    error = CAMBdata_getbispectrum(
        byref(res),
        byref(params),
        byref(bispectrum_params),
        byref(result),
        output_root_buffer,
        output_root_len,
    )
    if error or config.global_error_flag:
        config.check_global_error("get_bispectrum")
    return result


__all__ = ["BispectrumParams", "BispectrumResult", "get_bispectrum", "max_bispectra", "max_bispectrum_deltas"]
