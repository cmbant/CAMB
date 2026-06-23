"""Check numerical stability of CAMB parameters against a higher-accuracy reference run.

The command-line interface compares spectra and derived parameters from an input
``.ini`` file with a reference calculation using boosted accuracy settings.
Options such as ``--accuracy-boost``, ``--l-sample-boost``,
``--l-accuracy-boost``, ``--int-tol-boost``, ``--do-late-rad-truncation``, and
``--strict-reference`` set the boosted reference run; they do not change the
standard run being compared to that reference. The standard run is the input
``.ini`` calculation, except for shared setup overrides such as ``--lmax``,
``--set-for-lmax``, ``--lens-output-margin``, and
``--lens-potential-accuracy`` that are applied to both runs.

The tool can also plot fractional differences, estimate a fiducial CMB delta
chi-squared, search for minimal top-level boosts that pass against the same
high-accuracy reference, and refine which underlying component accuracy settings
are most relevant.

The main options generally probe numerical convergence rather than physical
convergence: the physical ``k_max`` is not normally changed between the input and
reference.  Exceptions include optional lensing parameters that increase the
convolution margin or ``k_max``, and nonlinear lensing runs where
``AccuracyBoost * NonlinSourceBoost`` can increase the nonlinear source
``Transfer.kmax`` through ``P%kmax = max(P%kmax, 5 * NL_Boost)`` in
``results.f90``.  Nonlinear power-spectrum modelling uncertainties may of course
be much larger than either numerical or ``k_max`` effects.

Examples::

    camb check_accuracy inifiles/params.ini
    camb check_accuracy inifiles/params.ini --plot-dir accuracy_plots
    camb check_accuracy inifiles/params.ini --strict-reference --reference-lens-output-margin 1500
    camb check_accuracy inifiles/params.ini --find-minimal-boosts --refine-accuracy-components
"""

from __future__ import annotations

import argparse
import itertools
import math
import time
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DEFAULT_ACCURACY_SETTINGS = {
    "AccuracyBoost": 2.0,
    "lSampleBoost": 2.0,
    "lAccuracyBoost": 2.0,
    "IntTolBoost": 2.0,
    "DoLateRadTruncation": True,
    # Force linear l-sampling all the way to max_l in the reference; the default min_l_logl_sampling
    # of 5000 leaves a wide gap in the log block right above 5000 that makes the reference itself
    # under-sampled at high l.
    "min_l_logl_sampling": 100000,
}

REFERENCE_MIN_TRANSFER_K_PER_LOGINT = 100

STRICT_REFERENCE_SETTINGS = {
    "AccuracyBoost": 3.0,
    "lSampleBoost": 3.0,
    "lAccuracyBoost": 3.0,
    "min_l_logl_sampling": 100000,
}

DEFAULT_NOISE_CONFIGS = {
    "planck": {
        "noise_muK_arcmin_T": 45.0,
        "noise_muK_arcmin_P": 70.0,
        "fwhm_arcmin": 7.0,
        "fsky": 0.6,
    },
    "so": {  # optimistic extended 2034 numbers
        "noise_muK_arcmin_T": 2.6,
        "noise_muK_arcmin_P": 3.67,
        "fwhm_arcmin": 1.4,
        "fsky": 0.6,
    },
}

TT_EE_TOLERANCES = [(0, 3e-3), (600, 1e-3), (3500, 3e-3), (6000, 2e-2)]
UNLENSED_TT_EE_TOLERANCES = [(0, 3e-3), (600, 1e-3), (3000, 3e-3), (3500, 5e-3), (6000, 5e-2)]
TENSOR_ONLY_TT_EE_TOLERANCES = [(0, 1e-1), (600, 3e-2), (2500, 1e-1)]
LOW_LMAX_LENSED_EDGE_LMAX = 2500
LOW_LMAX_LENSED_EDGE_TOLERANCE = 1e-2
LOW_LMAX_LENSED_EDGE_MIN_WIDTH = 300
LOW_LMAX_LENSED_EDGE_FRACTION = 0.15
LOW_LMAX_LENSED_EDGE_LABEL = "low-lmax edge"
BB_TOLERANCES = [(0, 5e-3), (1000, 1e-2), (6000, 2e-2), (8000, 1e-1)]
LENSING_TOLERANCES = [(0, 9e-3), (5, 5e-3), (2500, 1e-2), (4000, 2e-2)]
TENSOR_ONLY_RELATIVE_LMAX = 2000
TENSOR_ONLY_ABSOLUTE_LMIN = 600
TENSOR_ONLY_TAIL_RATIO_TOLERANCE = 1e-2
MPK_TOLERANCE = 1e-3
MPK_TOLERANCE_RANGES = [
    (None, 5e-3, 3e-3),
    (5e-3, 2.0, MPK_TOLERANCE),
    (2.0, None, 3e-3),
]
LOW_PRECISION_MPK_TOLERANCE_RANGES = [
    (None, 1.0, 3e-3),
    (1.0, None, 5e-3),
]
NONLINEAR_MPK_TOLERANCE_RANGES = [
    (None, 1.0, 3e-3),
    (1.0, 10.0, 5e-3),
    (10.0, None, 2e-2),
]
HIGH_PRECISION_NONLINEAR_MPK_TOLERANCE_RANGES = [
    (None, 5e-3, 3e-3),
    (5e-3, 2.0, MPK_TOLERANCE),
    (2.0, 10.0, 3e-3),
    (10.0, None, 1e-2),
]
DERIVED_TOLERANCE = 1e-3

CL_COLUMNS = {"TT": 0, "EE": 1, "BB": 2, "TE": 3}
OMEGA_K_FLAT = 5e-7
BASE_ACCURACY_KEYS = ("AccuracyBoost", "lSampleBoost", "lAccuracyBoost", "IntTolBoost")
GLOBAL_BOOST_COMPONENT_KEYS = (
    "TimeStepBoost",
    "BackgroundTimeStepBoost",
    "TimeSwitchBoost",
    "IntTolBoost",
    "SourcekAccuracyBoost",
    "IntkAccuracyBoost",
    "TransferkBoost",
    "NonFlatIntAccuracyBoost",
    "BessIntBoost",
    "LensingBoost",
    "NonlinSourceBoost",
    "BesselBoost",
    "LimberBoost",
    "SourceLimberBoost",
    "KmaxBoost",
    "neutrino_q_boost",
)
COMPONENT_REFINEMENT_KEYS = GLOBAL_BOOST_COMPONENT_KEYS
DISCRETE_ACCURACY_VALUES = {"neutrino_q_boost": (1.0, 1.5, 2.5)}

MatterPowerToleranceSpec = list[tuple[float | None, float | None, float]]


@dataclass(frozen=True)
class RangeTolerance:
    start: int
    stop: int | None
    tolerance: float
    note: str = ""

    @property
    def label(self) -> str:
        suffix = f" ({self.note})" if self.note else ""
        if self.stop is None:
            return f"L >= {self.start}{suffix}"
        return f"{self.start} <= L < {self.stop}{suffix}"


@dataclass
class RunOutput:
    label: str
    params: object
    results: object
    cpu_time: float
    wall_time: float
    derived: dict[str, float]
    lensed_cls: np.ndarray | None
    lens_potential_cls: np.ndarray | None
    matter_power: MatterPowerData | None


@dataclass(frozen=True)
class MatterPowerData:
    k: np.ndarray
    z: np.ndarray
    pk: np.ndarray
    requested_kmax: float
    npoints: int
    compare_at_input_nodes: bool = False


@dataclass
class StatRow:
    quantity: str
    range_label: str
    max_abs: float
    rms: float
    tolerance: float
    passed: bool
    location: str = ""


@dataclass
class ComparisonResult:
    derived_rows: list[StatRow]
    cl_rows: list[StatRow]
    lensing_rows: list[StatRow]
    mpk_rows: list[StatRow]

    @property
    def passed(self) -> bool:
        rows = self.derived_rows + self.cl_rows + self.lensing_rows + self.mpk_rows
        return all(row.passed for row in rows)

    @property
    def worst_failure(self) -> StatRow | None:
        failures = [
            row for row in self.derived_rows + self.cl_rows + self.lensing_rows + self.mpk_rows if not row.passed
        ]
        if not failures:
            return None
        return max(failures, key=stat_score)


def comparison_rows(comparison: ComparisonResult) -> list[StatRow]:
    return comparison.derived_rows + comparison.cl_rows + comparison.lensing_rows + comparison.mpk_rows


def stat_score(row: StatRow) -> float:
    if not math.isfinite(row.max_abs):
        return math.inf
    if row.tolerance == 0:
        return 0.0 if row.max_abs == 0 else math.inf
    score = abs(row.max_abs) / abs(row.tolerance)
    return score if math.isfinite(score) else math.inf


def comparison_score(comparison: ComparisonResult) -> float:
    rows = comparison_rows(comparison)
    if not rows:
        return 0.0
    return max(stat_score(row) for row in rows)


@dataclass(frozen=True)
class NoiseConfig:
    name: str
    noise_muK_arcmin_T: float
    noise_muK_arcmin_P: float
    fwhm_arcmin: float
    fsky: float
    lmin: int
    lmax: int | None
    fields: tuple[str, ...]


@dataclass(frozen=True)
class ChiSquaredResult:
    delta_chi2: float
    lmin: int
    lmax: int
    fields: tuple[str, ...]
    fsky: float
    noise_muK_arcmin_T: float
    noise_muK_arcmin_P: float
    fwhm_arcmin: float
    worst_ell: int
    worst_delta_chi2: float


@dataclass
class AccuracyCheckResult:
    standard: RunOutput
    reference: RunOutput
    comparison: ComparisonResult
    chi2: ChiSquaredResult | None = None


@dataclass(frozen=True)
class PassingCandidate:
    settings: dict[str, float | bool]
    comparison: ComparisonResult
    run: RunOutput | None


@dataclass(frozen=True)
class SearchResult:
    settings: dict[str, float | bool]
    comparison: ComparisonResult
    run: RunOutput | None = None
    fastest: PassingCandidate | None = None

    def __iter__(self):
        yield self.settings
        yield self.comparison


def tolerance_ranges(spec: list[tuple[int, float]]) -> list[RangeTolerance]:
    return [
        RangeTolerance(start, spec[index + 1][0] if index + 1 < len(spec) else None, tolerance)
        for index, (start, tolerance) in enumerate(spec)
    ]


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def comma_separated_floats(value: str) -> list[float]:
    values = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("must contain at least one value")
    return sorted(set(values))


def cmb_fields(value: str) -> tuple[str, ...]:
    fields = tuple(field.strip().upper() for field in value.replace(",", " ").split() if field.strip())
    allowed = {"T", "E", "B"}
    if not fields or any(field not in allowed for field in fields):
        raise argparse.ArgumentTypeError("fields must be a comma- or space-separated subset of T,E,B")
    if len(set(fields)) != len(fields):
        raise argparse.ArgumentTypeError("fields must not contain duplicates")
    return fields


def build_parser(prog: str | None = None) -> argparse.ArgumentParser:
    """Build the command-line parser for the accuracy checker."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Check stability of one CAMB ini file against a boosted-accuracy calculation.",
    )
    parser.add_argument("ini_file", help="CAMB ini file to load with camb.read_ini")
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="pass no_validate=True to camb.read_ini",
    )
    parser.add_argument(
        "--lmax",
        type=positive_int,
        help=(
            "maximum L to compare; default is the requested lensed-CMB/lensing-potential lmax "
            "(l_max_scalar-lens_output_margin) and common calculated lmax for other spectra"
        ),
    )
    parser.add_argument(
        "--set-for-lmax",
        type=positive_int,
        help="call params.set_for_lmax before running; if --lmax is unset, this is also used for comparison",
    )
    parser.add_argument(
        "--lens-output-margin",
        type=int,
        help="lens output margin to use for both runs; can be set with or without --set-for-lmax",
    )
    parser.add_argument(
        "--lens-potential-accuracy",
        type=float,
        help="lens_potential_accuracy to use for both runs; can be set with or without --set-for-lmax",
    )
    parser.add_argument(
        "--reference-lens-output-margin",
        type=int,
        help="lens output margin override for the boosted reference only; use this with sufficiently converged reference accuracy settings to probe lensed-spectrum support sensitivity near the output cutoff",
    )
    parser.add_argument(
        "--reference-lens-potential-accuracy",
        type=float,
        help="lens_potential_accuracy override for the boosted reference only",
    )
    parser.add_argument(
        "--mpk-npoints",
        type=positive_int,
        default=500,
        help="number of k samples for matter power comparison grid",
    )
    parser.add_argument(
        "--mpk-kmin",
        type=float,
        default=1e-4,
        help="minimum k/h for matter power comparison",
    )
    parser.add_argument(
        "--mpk-tolerance",
        type=float,
        default=None,
        help="override matter power tolerance for all k/h ranges",
    )
    parser.add_argument("--derived-tolerance", type=float, default=DERIVED_TOLERANCE)
    parser.add_argument(
        "--plot-dir",
        type=Path,
        help="write plots of fractional errors to this directory",
    )
    parser.add_argument(
        "--chi2",
        action="store_true",
        help="calculate a fiducial Gaussian CMB delta chi-squared for standard vs boosted C_L",
    )
    parser.add_argument(
        "--chi2-config",
        choices=tuple(DEFAULT_NOISE_CONFIGS) + ("custom",),
        default="so",
        help="noise configuration for --chi2",
    )
    parser.add_argument("--chi2-lmin", type=positive_int, default=2)
    parser.add_argument(
        "--chi2-lmax",
        type=positive_int,
        help="maximum L for --chi2; defaults to comparison lmax",
    )
    parser.add_argument("--chi2-fields", type=cmb_fields, default=("T", "E", "B"))
    parser.add_argument("--noise-muk-arcmin-t", type=float, help="temperature white noise for --chi2")
    parser.add_argument("--noise-muk-arcmin-p", type=float, help="polarization white noise for --chi2")
    parser.add_argument("--beam-fwhm-arcmin", type=float, help="Gaussian beam FWHM for --chi2")
    parser.add_argument("--fsky", type=float, help="sky fraction for --chi2")
    parser.add_argument(
        "--find-minimal-boosts",
        action="store_true",
        help="search for the lowest-cost boosted settings that pass against the high-accuracy run",
    )
    parser.add_argument(
        "--exhaustive-boost-search",
        action="store_true",
        help="continue boost search after a passing single-parameter path is found",
    )
    parser.add_argument(
        "--refine-accuracy-components",
        action="store_true",
        help="after --find-minimal-boosts, try AccuracyBoost=1 with underlying accuracy fields tweaked instead",
    )
    parser.add_argument(
        "--component-search-size",
        type=positive_int,
        default=4,
        help="maximum number of component accuracy fields to try together",
    )
    parser.add_argument(
        "--search-grid",
        type=comma_separated_floats,
        help="optional comma-separated seed values for the four boost parameters; the default search does not need this",
    )
    parser.add_argument(
        "--search-tolerance",
        type=float,
        default=0.02,
        help="target precision for refining numeric boost values",
    )
    parser.add_argument(
        "--max-search-runs",
        type=positive_int,
        help="maximum candidate runs for boost search",
    )

    parser.add_argument(
        "--accuracy-boost",
        type=float,
        default=DEFAULT_ACCURACY_SETTINGS["AccuracyBoost"],
    )
    parser.add_argument(
        "--strict-reference",
        action="store_true",
        help="use a stricter boosted reference preset with AccuracyBoost=lSampleBoost=lAccuracyBoost=3 and min_l_logl_sampling=10000",
    )
    parser.add_argument(
        "--l-sample-boost",
        type=float,
        default=DEFAULT_ACCURACY_SETTINGS["lSampleBoost"],
    )
    parser.add_argument(
        "--l-accuracy-boost",
        type=float,
        default=DEFAULT_ACCURACY_SETTINGS["lAccuracyBoost"],
    )
    parser.add_argument("--int-tol-boost", type=float, default=DEFAULT_ACCURACY_SETTINGS["IntTolBoost"])
    parser.add_argument(
        "--do-late-rad-truncation",
        action=argparse.BooleanOptionalAction,
        default=DEFAULT_ACCURACY_SETTINGS["DoLateRadTruncation"],
    )
    return parser


def requested_accuracy_settings(args: argparse.Namespace) -> dict[str, float | bool]:
    settings: dict[str, float | bool] = {
        "AccuracyBoost": args.accuracy_boost,
        "lSampleBoost": args.l_sample_boost,
        "lAccuracyBoost": args.l_accuracy_boost,
        "IntTolBoost": args.int_tol_boost,
        "DoLateRadTruncation": args.do_late_rad_truncation,
    }
    if args.strict_reference:
        for key, value in STRICT_REFERENCE_SETTINGS.items():
            if key in settings:
                settings[key] = max(float(settings[key]), float(value))
            else:
                settings[key] = value
    return settings


def apply_accuracy_settings(params, settings: dict[str, float | bool], *, boost_from_raw: bool = False) -> None:
    for key, setting in settings.items():
        if key == "DoLateRadTruncation":
            params.DoLateRadTruncation = bool(setting)
            continue
        if hasattr(params.Accuracy, key):
            current_value = getattr(params.Accuracy, key)
            if isinstance(current_value, bool):
                value = bool(setting)
                if boost_from_raw:
                    value = current_value or value
            else:
                value = float(setting)
                if boost_from_raw:
                    value = max(float(current_value), value)
            setattr(params.Accuracy, key, value)
            continue
        if key == "min_l_logl_sampling":
            value = int(setting)
            if boost_from_raw:
                value = max(int(params.min_l_logl_sampling), value)
            params.min_l_logl_sampling = value
            continue
        raise ValueError(f"Unknown accuracy setting {key}")


def copy_params(params):
    return params.copy()


def apply_lensing_settings(
    params,
    *,
    set_for_lmax: int | None = None,
    lens_output_margin: int | None = None,
    lens_potential_accuracy: float | None = None,
) -> None:
    if set_for_lmax is None and lens_output_margin is None and lens_potential_accuracy is None:
        return

    lens_accuracy = 0.0 if lens_potential_accuracy is None else lens_potential_accuracy
    if set_for_lmax is not None:
        params.set_for_lmax(
            set_for_lmax,
            lens_output_margin=200 if lens_output_margin is None else lens_output_margin,
            lens_potential_accuracy=lens_accuracy,
            nonlinear=current_uses_nonlinear_lensing(params),
        )
        return

    margin = 0 if lens_output_margin is None else lens_output_margin
    target_lmax = params.max_l - margin if params.DoLensing else params.max_l
    max_eta_k = params.max_eta_k if lens_potential_accuracy is None else None
    params.set_for_lmax(
        max(1, target_lmax),
        max_eta_k=max_eta_k,
        lens_output_margin=margin,
        lens_potential_accuracy=lens_accuracy,
        nonlinear=current_uses_nonlinear_lensing(params),
    )


def current_uses_nonlinear_lensing(params) -> bool:
    return str(params.NonLinear) in {"NonLinear_lens", "NonLinear_both"}


def default_mpk_tolerance(params) -> float | MatterPowerToleranceSpec:
    high_precision = bool(params.Transfer.high_precision)
    nonlinear = str(params.NonLinear) in {"NonLinear_pk", "NonLinear_both"}
    if nonlinear:
        return HIGH_PRECISION_NONLINEAR_MPK_TOLERANCE_RANGES if high_precision else NONLINEAR_MPK_TOLERANCE_RANGES
    return MPK_TOLERANCE_RANGES if high_precision else LOW_PRECISION_MPK_TOLERANCE_RANGES


def load_params(
    ini_file: Path,
    *,
    no_validate: bool,
    settings: dict[str, float | bool] | None = None,
    set_for_lmax: int | None = None,
    lens_output_margin: int | None = None,
    lens_potential_accuracy: float | None = None,
):
    import camb

    params = camb.read_ini(str(ini_file), no_validate=no_validate)
    apply_lensing_settings(
        params,
        set_for_lmax=set_for_lmax,
        lens_output_margin=lens_output_margin,
        lens_potential_accuracy=lens_potential_accuracy,
    )
    if settings:
        apply_accuracy_settings(params, settings)
    return params


def run_params_case(
    params,
    label: str,
    *,
    lmax: int | None,
    mpk_kmin: float,
    mpk_npoints: int,
) -> RunOutput:
    import camb

    requested_matter_kmax = params.Transfer.kmax / (params.H0 / 100.0) if params.WantTransfer else None
    cpu_start = time.process_time()
    wall_start = time.perf_counter()
    results = camb.get_results(params)
    cpu_time = time.process_time() - cpu_start
    wall_time = time.perf_counter() - wall_start

    derived = results.get_derived_params() if params.WantDerivedParameters else {}
    lensed_cls = get_lensed_cls(results, lmax)
    lens_potential_cls = get_lens_potential_cls(results, lmax)
    matter_power = get_matter_power(results, mpk_kmin, mpk_npoints, requested_kmax=requested_matter_kmax)

    return RunOutput(
        label=label,
        params=params,
        results=results,
        cpu_time=cpu_time,
        wall_time=wall_time,
        derived=derived,
        lensed_cls=lensed_cls,
        lens_potential_cls=lens_potential_cls,
        matter_power=matter_power,
    )


def run_case(
    ini_file: Path,
    label: str,
    *,
    no_validate: bool,
    accuracy_settings: dict[str, float | bool] | None,
    lmax: int | None,
    set_for_lmax: int | None,
    lens_output_margin: int | None,
    lens_potential_accuracy: float | None,
    mpk_kmin: float,
    mpk_npoints: int,
) -> RunOutput:
    params = load_params(
        ini_file,
        no_validate=no_validate,
        settings=accuracy_settings,
        set_for_lmax=set_for_lmax,
        lens_output_margin=lens_output_margin,
        lens_potential_accuracy=lens_potential_accuracy,
    )
    return run_params_case(
        params,
        label,
        lmax=lmax,
        mpk_kmin=mpk_kmin,
        mpk_npoints=mpk_npoints,
    )


def run_timing_summary(run: RunOutput) -> str:
    return f"cpu={run.cpu_time:.2f}s wall={run.wall_time:.2f}s"


def run_timing_key(run: RunOutput) -> tuple[float, float]:
    return run.cpu_time, run.wall_time


def get_lensed_cls(results, lmax: int | None) -> np.ndarray | None:
    params = results.Params
    if not (params.WantCls and params.Want_CMB):
        return None
    if bool(params.WantTensors) and not bool(params.WantScalars):
        tensor_lmax = int(params.max_l_tensor)
        lmax = tensor_lmax if lmax is None else min(lmax, tensor_lmax)
    elif bool(params.DoLensing):
        requested_lmax = max(0, int(params.max_l) - int(params.lens_output_margin))
        lmax = requested_lmax if lmax is None else min(lmax, requested_lmax)
    return results.get_total_cls(lmax)


def get_lens_potential_cls(results, lmax: int | None) -> np.ndarray | None:
    params = results.Params
    if not (params.WantCls and params.Want_CMB_lensing):
        return None
    if bool(params.DoLensing):
        requested_lmax = max(0, int(params.max_l) - int(params.lens_output_margin))
        lmax = requested_lmax if lmax is None else min(lmax, requested_lmax)
    return results.get_lens_potential_cls(lmax)


def get_matter_power(
    results,
    minkh: float,
    npoints: int,
    requested_kmax: float | None = None,
) -> MatterPowerData | None:
    params = results.Params
    if not params.WantTransfer:
        return None
    requested_kmax = requested_kmax if requested_kmax is not None else params.Transfer.kmax / (params.H0 / 100.0)
    if requested_kmax <= minkh:
        return None
    nonlinear = str(params.NonLinear) in {"NonLinear_pk", "NonLinear_both"}
    if params.Transfer.k_per_logint != 0:
        if nonlinear:
            k, z, pk = results.get_nonlinear_matter_power_spectrum()
        else:
            k, z, pk = results.get_linear_matter_power_spectrum()
        mask = (k >= minkh) & (k <= requested_kmax)
        if not np.any(mask):
            return None
        k = k[mask]
        pk = pk[:, mask]
        return MatterPowerData(
            k=k,
            z=z,
            pk=pk,
            requested_kmax=float(k[-1]),
            npoints=len(k),
            compare_at_input_nodes=True,
        )
    pk_interpolator, z, _ = results.get_matter_power_interpolator(
        nonlinear=nonlinear,
        return_z_k=True,
        silent=True,
    )
    kmin = max(minkh, pk_interpolator.kmin)
    kmax = min(requested_kmax, pk_interpolator.kmax)
    if kmax <= kmin:
        return None
    k = np.exp(np.linspace(np.log(kmin), np.log(kmax), npoints))
    pk = pk_interpolator.P(z, k)
    return MatterPowerData(k=k, z=z, pk=pk, requested_kmax=kmax, npoints=npoints)


def fractional_delta(values: np.ndarray, reference: np.ndarray, *, floor: float = 0.0) -> np.ndarray:
    denominator = np.abs(reference)
    if floor:
        denominator = np.maximum(denominator, floor)
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = (values - reference) / denominator
    delta[(denominator == 0) & (values == reference)] = 0.0
    return delta


def normalized_te_delta(values: np.ndarray, reference: np.ndarray) -> np.ndarray:
    denominator = np.sqrt(np.maximum(reference[:, CL_COLUMNS["TT"]] * reference[:, CL_COLUMNS["EE"]], 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = (values[:, CL_COLUMNS["TE"]] - reference[:, CL_COLUMNS["TE"]]) / denominator
    delta[(denominator == 0) & (values[:, CL_COLUMNS["TE"]] == reference[:, CL_COLUMNS["TE"]])] = 0.0
    return delta


def finite_stats(
    quantity: str,
    range_label: str,
    errors: np.ndarray,
    tolerance: float,
    locations: np.ndarray | None = None,
) -> StatRow:
    finite = np.isfinite(errors)
    if not np.all(finite):
        bad_index = int(np.flatnonzero(~finite)[0])
        location = ""
        if locations is not None:
            location = str(locations[bad_index])
        return StatRow(
            quantity,
            range_label,
            math.inf,
            math.inf,
            tolerance,
            False,
            location or "non-finite sample",
        )
    if not np.any(finite):
        return StatRow(
            quantity,
            range_label,
            math.inf,
            math.inf,
            tolerance,
            False,
            "no finite samples",
        )
    valid_errors = errors[finite]
    max_index = int(np.argmax(np.abs(valid_errors)))
    max_abs = float(np.abs(valid_errors[max_index]))
    rms = float(np.sqrt(np.mean(valid_errors**2)))
    location = ""
    if locations is not None:
        location = str(locations[finite][max_index])
    return StatRow(quantity, range_label, max_abs, rms, tolerance, max_abs <= tolerance, location)


def compare_derived(standard: RunOutput, reference: RunOutput, tolerance: float) -> list[StatRow]:
    rows = []
    for name in reference.derived:
        if name not in standard.derived:
            continue
        ref = float(reference.derived[name])
        value = float(standard.derived[name])
        denominator = abs(ref) if ref else 1.0
        error = abs((value - ref) / denominator)
        rows.append(
            StatRow(
                f"derived:{name}",
                "fractional",
                error,
                error,
                tolerance,
                error <= tolerance,
            )
        )
    return rows


def compare_cls(standard: RunOutput, reference: RunOutput) -> list[StatRow]:
    if standard.lensed_cls is None or reference.lensed_cls is None:
        return []
    lmax = min(standard.lensed_cls.shape[0], reference.lensed_cls.shape[0]) - 1
    standard_cls = standard.lensed_cls[: lmax + 1]
    reference_cls = reference.lensed_cls[: lmax + 1]
    ell = np.arange(lmax + 1)
    if bool(standard.params.WantTensors) and not bool(standard.params.WantScalars):
        return compare_tensor_only_cls(standard_cls, reference_cls, ell)

    rows = []
    tt_ee_tolerances = tt_ee_tolerance_ranges_for_cls(standard.params, lmax)

    for quantity in ("TT", "EE"):
        errors = fractional_delta(
            standard_cls[:, CL_COLUMNS[quantity]],
            reference_cls[:, CL_COLUMNS[quantity]],
        )
        rows.extend(compare_l_ranges(quantity, errors, ell, tt_ee_tolerances, min_ell=2))

    te_errors = normalized_te_delta(standard_cls, reference_cls)
    rows.extend(compare_l_ranges("TE", te_errors, ell, tt_ee_tolerances, min_ell=2))

    bb_errors = fractional_delta(standard_cls[:, CL_COLUMNS["BB"]], reference_cls[:, CL_COLUMNS["BB"]])
    rows.extend(
        compare_l_ranges(
            "BB",
            bb_errors,
            ell,
            tolerance_ranges(bb_tolerances(standard.params)),
            min_ell=2,
        )
    )
    return rows


def tt_ee_tolerance_ranges_for_cls(params, lmax: int) -> list[RangeTolerance]:
    tolerances = tolerance_ranges(tt_ee_tolerances_for_params(params))
    if not use_low_lmax_lensed_edge_tolerance(params, lmax):
        return tolerances

    edge_width = max(
        LOW_LMAX_LENSED_EDGE_MIN_WIDTH,
        math.ceil(LOW_LMAX_LENSED_EDGE_FRACTION * lmax),
    )
    edge_start = max(0, lmax - edge_width + 1)
    ranges = [range_tolerance for range_tolerance in tolerances if range_tolerance.start < edge_start]
    if ranges:
        last_range = ranges[-1]
        if last_range.stop is None or last_range.stop > edge_start:
            ranges[-1] = RangeTolerance(last_range.start, edge_start, last_range.tolerance)
    ranges = [
        range_tolerance
        for range_tolerance in ranges
        if range_tolerance.stop is None or range_tolerance.start < range_tolerance.stop
    ]
    ranges.append(
        RangeTolerance(
            edge_start,
            None,
            LOW_LMAX_LENSED_EDGE_TOLERANCE,
            LOW_LMAX_LENSED_EDGE_LABEL,
        )
    )
    return ranges


def use_low_lmax_lensed_edge_tolerance(params, lmax: int) -> bool:
    return (
        lmax < LOW_LMAX_LENSED_EDGE_LMAX
        and bool(params.DoLensing)
        and bool(params.WantScalars)
        and bool(params.WantCls)
        and bool(params.Want_CMB)
        and not (bool(params.WantTensors) and not bool(params.WantScalars))
    )


def compare_tensor_only_cls(standard_cls: np.ndarray, reference_cls: np.ndarray, ell: np.ndarray) -> list[StatRow]:
    rows = []
    relative_lmax = min(int(ell[-1]), TENSOR_ONLY_RELATIVE_LMAX)
    low_mask = (ell >= 2) & (ell <= relative_lmax)
    relative_ranges = [
        RangeTolerance(0, 600, 1e-1),
        RangeTolerance(600, relative_lmax + 1, 3e-2),
    ]
    relative_ranges = [range_tolerance for range_tolerance in relative_ranges if range_tolerance.start <= relative_lmax]

    for quantity in ("TT", "EE"):
        errors = tensor_only_relative_or_small_delta(
            standard_cls[:, CL_COLUMNS[quantity]],
            reference_cls[:, CL_COLUMNS[quantity]],
            ell,
            low_mask,
        )
        rows.extend(compare_l_ranges(quantity, errors, ell, relative_ranges, min_ell=2))

    te_errors = tensor_only_relative_or_small_delta(
        standard_cls[:, CL_COLUMNS["TE"]],
        reference_cls[:, CL_COLUMNS["TE"]],
        ell,
        low_mask,
        relative_errors=normalized_te_delta(standard_cls, reference_cls),
    )
    rows.extend(compare_l_ranges("TE", te_errors, ell, relative_ranges, min_ell=2))

    bb_errors = tensor_only_relative_or_small_delta(
        standard_cls[:, CL_COLUMNS["BB"]],
        reference_cls[:, CL_COLUMNS["BB"]],
        ell,
        low_mask,
    )
    rows.extend(
        compare_l_ranges(
            "BB",
            bb_errors,
            ell,
            [RangeTolerance(0, relative_lmax + 1, 2e-2)],
            min_ell=2,
        )
    )

    rows.extend(compare_tensor_only_tail(standard_cls, reference_cls, ell, relative_lmax))
    return rows


def tensor_only_relative_or_small_delta(
    values: np.ndarray,
    reference: np.ndarray,
    ell: np.ndarray,
    scale_mask: np.ndarray,
    *,
    relative_errors: np.ndarray | None = None,
) -> np.ndarray:
    errors = fractional_delta(values, reference) if relative_errors is None else relative_errors.copy()
    if not np.any(scale_mask):
        return errors
    low_scale = max(float(np.max(np.abs(values[scale_mask]))), float(np.max(np.abs(reference[scale_mask]))))
    if low_scale == 0:
        return errors
    small_signal_errors = np.maximum(np.abs(values), np.abs(reference)) / low_scale
    use_absolute_scale = (ell >= TENSOR_ONLY_ABSOLUTE_LMIN) & (small_signal_errors < np.abs(errors))
    errors[use_absolute_scale] = np.copysign(small_signal_errors[use_absolute_scale], errors[use_absolute_scale])
    return errors


def compare_tensor_only_tail(
    standard_cls: np.ndarray,
    reference_cls: np.ndarray,
    ell: np.ndarray,
    relative_lmax: int,
) -> list[StatRow]:
    tail_mask = ell > relative_lmax
    if not np.any(tail_mask):
        return []

    low_mask = (ell >= 2) & (ell <= relative_lmax)
    rows = []
    for quantity, column in CL_COLUMNS.items():
        low_scale = max(
            float(np.max(np.abs(standard_cls[low_mask, column]))),
            float(np.max(np.abs(reference_cls[low_mask, column]))),
        )
        tail_values = np.maximum(np.abs(standard_cls[tail_mask, column]), np.abs(reference_cls[tail_mask, column]))
        if low_scale == 0:
            errors = tail_values
            tolerance = 0.0
        else:
            errors = tail_values / low_scale
            tolerance = TENSOR_ONLY_TAIL_RATIO_TOLERANCE
        rows.append(
            finite_stats(
                f"{quantity} tensor tail",
                f"L > {relative_lmax}",
                errors,
                tolerance,
                locations=ell[tail_mask],
            )
        )
    return rows


def tt_ee_tolerances_for_params(params) -> list[tuple[int, float]]:
    if bool(params.WantTensors) and not bool(params.WantScalars):
        return TENSOR_ONLY_TT_EE_TOLERANCES
    if not bool(params.DoLensing):
        return UNLENSED_TT_EE_TOLERANCES
    return TT_EE_TOLERANCES


def bb_tolerances(params) -> list[tuple[int, float]]:
    if not bool(params.Accuracy.AccurateBB):
        return [(0, 2e-2), (8000, 1e-1)]
    return BB_TOLERANCES


def compare_lensing(standard: RunOutput, reference: RunOutput) -> list[StatRow]:
    if standard.lens_potential_cls is None or reference.lens_potential_cls is None:
        return []
    lmax = min(standard.lens_potential_cls.shape[0], reference.lens_potential_cls.shape[0]) - 1
    standard_pp = standard.lens_potential_cls[: lmax + 1, 0]
    reference_pp = reference.lens_potential_cls[: lmax + 1, 0]
    ell = np.arange(lmax + 1)
    errors = fractional_delta(standard_pp, reference_pp)
    return compare_l_ranges("lens PP", errors, ell, tolerance_ranges(LENSING_TOLERANCES), min_ell=2)


def compare_l_ranges(
    quantity: str,
    errors: np.ndarray,
    ell: np.ndarray,
    ranges: Iterable[RangeTolerance],
    *,
    min_ell: int,
) -> list[StatRow]:
    rows = []
    for range_tolerance in ranges:
        mask = ell >= max(min_ell, range_tolerance.start)
        if range_tolerance.stop is not None:
            mask &= ell < range_tolerance.stop
        if not np.any(mask):
            continue
        selected_ell = ell[mask]
        rows.append(
            finite_stats(
                quantity,
                ell_range_label(selected_ell, range_tolerance),
                errors[mask],
                range_tolerance.tolerance,
                locations=ell[mask],
            )
        )
    return rows


def ell_range_label(selected_ell: np.ndarray, range_tolerance: RangeTolerance) -> str:
    start = int(selected_ell[0])
    stop = int(selected_ell[-1])
    suffix = f" ({range_tolerance.note})" if range_tolerance.note else ""
    if range_tolerance.stop is None:
        return f"L >= {start}{suffix}"
    if stop < range_tolerance.stop - 1:
        return f"{start} <= L <= {stop}{suffix}"
    return range_tolerance.label


def matter_power_range_label(start: float | None, stop: float | None) -> str:
    if start is None:
        return f"all z, k/h < {stop:.4g}"
    if stop is None:
        return f"all z, {start:.4g} <= k/h"
    return f"all z, {start:.4g} <= k/h < {stop:.4g}"


def compare_matter_power(
    standard: RunOutput,
    reference: RunOutput,
    tolerance: float | MatterPowerToleranceSpec,
) -> list[StatRow]:
    if standard.matter_power is None or reference.matter_power is None:
        return []
    common_k, z, standard_pk, reference_pk = common_matter_power_grid(standard.matter_power, reference.matter_power)
    errors = fractional_delta(standard_pk, reference_pk)
    if isinstance(tolerance, float):
        z_grid, k_grid = np.meshgrid(z, common_k, indexing="ij")
        locations = np.array([f"z={z_value:.6g}, k/h={k:.6g}" for z_value, k in zip(z_grid.ravel(), k_grid.ravel())])
        range_label = f"all z, {common_k[0]:.4g} <= k/h <= {common_k[-1]:.4g}"
        return [
            finite_stats(
                "matter P(k)",
                range_label,
                errors.ravel(),
                tolerance,
                locations=locations,
            )
        ]

    rows = []
    for start, stop, range_tolerance in tolerance:
        mask = np.ones_like(common_k, dtype=bool)
        if start is not None:
            mask &= common_k >= start
        if stop is not None:
            mask &= common_k < stop
        if not np.any(mask):
            continue
        selected_k = common_k[mask]
        selected_errors = errors[:, mask]
        z_grid, k_grid = np.meshgrid(z, selected_k, indexing="ij")
        locations = np.array([f"z={z_value:.6g}, k/h={k:.6g}" for z_value, k in zip(z_grid.ravel(), k_grid.ravel())])
        rows.append(
            finite_stats(
                "matter P(k)",
                matter_power_range_label(start, stop),
                selected_errors.ravel(),
                range_tolerance,
                locations=locations,
            )
        )
    return rows


def common_matter_power_grid(
    standard_matter_power: MatterPowerData,
    reference_matter_power: MatterPowerData,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    std_k, std_z, std_pk = (
        standard_matter_power.k,
        standard_matter_power.z,
        standard_matter_power.pk,
    )
    ref_k, ref_z, ref_pk = (
        reference_matter_power.k,
        reference_matter_power.z,
        reference_matter_power.pk,
    )
    if std_pk.shape[0] != ref_pk.shape[0] or not np.allclose(std_z, ref_z):
        raise ValueError("Matter power redshift grids differ")

    kmin = max(float(std_k[0]), float(ref_k[0]))
    requested_kmax = min(standard_matter_power.requested_kmax, reference_matter_power.requested_kmax)
    kmax = min(float(std_k[-1]), float(ref_k[-1]), requested_kmax)
    if kmax <= kmin:
        raise ValueError("Matter power k grids do not overlap")

    if standard_matter_power.compare_at_input_nodes and reference_matter_power.compare_at_input_nodes:
        std_mask = (std_k >= kmin) & (std_k <= kmax)
        common_k = std_k[std_mask]
        if common_k.size == 0:
            raise ValueError("Matter power input-node k grids do not overlap")
        standard_pk = std_pk[:, std_mask]
        reference_pk = interpolate_pk_to_grid(ref_k, ref_pk, common_k)
        return common_k, std_z, standard_pk, reference_pk

    npoints = min(standard_matter_power.npoints, reference_matter_power.npoints)
    common_k = np.exp(np.linspace(np.log(kmin), np.log(kmax), npoints))
    standard_pk = interpolate_pk_to_grid(std_k, std_pk, common_k)
    reference_pk = interpolate_pk_to_grid(ref_k, ref_pk, common_k)
    return common_k, std_z, standard_pk, reference_pk


def interpolate_pk_to_grid(k: np.ndarray, pk: np.ndarray, common_k: np.ndarray) -> np.ndarray:
    log_k = np.log(k)
    log_common_k = np.log(common_k)
    interpolated = np.empty((pk.shape[0], len(common_k)))
    for index, values in enumerate(pk):
        if np.all(values > 0):
            interpolated[index] = np.exp(np.interp(log_common_k, log_k, np.log(values)))
        else:
            interpolated[index] = np.interp(log_common_k, log_k, values)
    return interpolated


def compare_runs(
    standard: RunOutput,
    reference: RunOutput,
    *,
    derived_tolerance: float,
    mpk_tolerance: float | MatterPowerToleranceSpec,
) -> ComparisonResult:
    return ComparisonResult(
        derived_rows=compare_derived(standard, reference, derived_tolerance),
        cl_rows=compare_cls(standard, reference),
        lensing_rows=compare_lensing(standard, reference),
        mpk_rows=compare_matter_power(standard, reference, mpk_tolerance),
    )


def compare_params_accuracy(
    params,
    *,
    reference_accuracy_settings: dict[str, float | bool] | None = None,
    comparator_accuracy_settings: dict[str, float | bool] | None = None,
    lmax: int | None = None,
    set_for_lmax: int | None = None,
    lens_output_margin: int | None = None,
    lens_potential_accuracy: float | None = None,
    reference_lens_output_margin: int | None = None,
    reference_lens_potential_accuracy: float | None = None,
    mpk_kmin: float = 1e-4,
    mpk_npoints: int = 500,
    derived_tolerance: float = DERIVED_TOLERANCE,
    mpk_tolerance: float | MatterPowerToleranceSpec | None = None,
    chi2_config: NoiseConfig | None = None,
    plot_dir: Path | None = None,
) -> AccuracyCheckResult:
    """Compare a :class:`~camb.model.CAMBparams` object against a higher-accuracy copy.

    The input parameters are copied before modification. The returned result
    includes both CAMB runs, row-wise comparison statistics, and optionally a
    fiducial CMB delta chi-squared estimate.
    """
    if mpk_tolerance is None:
        mpk_tolerance = default_mpk_tolerance(params)
    reference_accuracy_settings = reference_accuracy_settings or DEFAULT_ACCURACY_SETTINGS

    standard_params = copy_params(params)
    apply_lensing_settings(
        standard_params,
        set_for_lmax=set_for_lmax,
        lens_output_margin=lens_output_margin,
        lens_potential_accuracy=lens_potential_accuracy,
    )
    if comparator_accuracy_settings:
        apply_accuracy_settings(standard_params, comparator_accuracy_settings)

    reference_params = copy_params(params)
    apply_lensing_settings(
        reference_params,
        set_for_lmax=set_for_lmax,
        lens_output_margin=reference_lens_output_margin
        if reference_lens_output_margin is not None
        else lens_output_margin,
        lens_potential_accuracy=(
            reference_lens_potential_accuracy
            if reference_lens_potential_accuracy is not None
            else lens_potential_accuracy
        ),
    )
    apply_accuracy_settings(reference_params, reference_accuracy_settings, boost_from_raw=True)
    if int(reference_params.Transfer.k_per_logint) > 0:
        # Fixed transfer-output checks compare on the standard requested nodes.
        # Use a denser boosted native grid so reference interpolation is not the
        # limiting error for sparse requested output grids.
        reference_params.Transfer.k_per_logint = max(
            int(reference_params.Transfer.k_per_logint),
            REFERENCE_MIN_TRANSFER_K_PER_LOGINT,
        )

    lmax = lmax or set_for_lmax
    standard = run_params_case(
        standard_params,
        "standard",
        lmax=lmax,
        mpk_kmin=mpk_kmin,
        mpk_npoints=mpk_npoints,
    )
    reference = run_params_case(
        reference_params,
        "boosted",
        lmax=lmax,
        mpk_kmin=mpk_kmin,
        mpk_npoints=mpk_npoints,
    )
    comparison = compare_runs(
        standard,
        reference,
        derived_tolerance=derived_tolerance,
        mpk_tolerance=mpk_tolerance,
    )
    chi2 = calculate_cmb_delta_chi2(standard, reference, chi2_config) if chi2_config else None
    if plot_dir:
        plot_errors(standard, reference, plot_dir)
    return AccuracyCheckResult(standard=standard, reference=reference, comparison=comparison, chi2=chi2)


def chi2_noise_config(args: argparse.Namespace) -> NoiseConfig:
    if args.chi2_config == "custom":
        defaults = {}
    else:
        defaults = DEFAULT_NOISE_CONFIGS[args.chi2_config]

    noise_t = args.noise_muk_arcmin_t
    if noise_t is None:
        noise_t = defaults.get("noise_muK_arcmin_T")
    noise_p = args.noise_muk_arcmin_p
    if noise_p is None:
        noise_p = defaults.get("noise_muK_arcmin_P", noise_t)
    fwhm = args.beam_fwhm_arcmin
    if fwhm is None:
        fwhm = defaults.get("fwhm_arcmin")
    fsky = args.fsky
    if fsky is None:
        fsky = defaults.get("fsky", 1.0)

    missing = [
        name
        for name, value in (
            ("noise_muK_arcmin_T", noise_t),
            ("noise_muK_arcmin_P", noise_p),
            ("fwhm_arcmin", fwhm),
        )
        if value is None
    ]
    if missing:
        raise ValueError(f"--chi2-config custom requires {', '.join(missing)}")
    if fsky <= 0:
        raise ValueError("--fsky must be positive")

    return NoiseConfig(
        name=args.chi2_config,
        noise_muK_arcmin_T=float(noise_t),
        noise_muK_arcmin_P=float(noise_p),
        fwhm_arcmin=float(fwhm),
        fsky=float(fsky),
        lmin=args.chi2_lmin,
        lmax=args.chi2_lmax or effective_lmax(args),
        fields=args.chi2_fields,
    )


def white_noise_from_muK_arcmin(noise_muK_arcmin: float) -> float:
    return (noise_muK_arcmin * np.pi / (180.0 * 60.0)) ** 2


def make_cmb_noise_spectra(ell: np.ndarray, config: NoiseConfig) -> dict[str, np.ndarray]:
    noise_var_t = white_noise_from_muK_arcmin(config.noise_muK_arcmin_T)
    noise_var_p = white_noise_from_muK_arcmin(config.noise_muK_arcmin_P)
    fwhm_degrees = config.fwhm_arcmin / 60.0
    xlc = 180.0 * np.sqrt(8.0 * np.log(2.0)) / np.pi
    sigma2 = (fwhm_degrees / xlc) ** 2
    beam_fac = np.exp(ell * (ell + 1.0) * sigma2)
    dl_fac = ell * (ell + 1.0) / (2.0 * np.pi)
    return {
        "T": dl_fac * noise_var_t * beam_fac,
        "E": dl_fac * noise_var_p * beam_fac,
        "B": dl_fac * noise_var_p * beam_fac,
    }


def calculate_cmb_delta_chi2(standard: RunOutput, reference: RunOutput, config: NoiseConfig) -> ChiSquaredResult | None:
    if not (
        standard.params.WantCls and reference.params.WantCls and standard.params.Want_CMB and reference.params.Want_CMB
    ):
        return None

    lmax_available = min(standard.results.Params.max_l, reference.results.Params.max_l)
    if config.lmax is not None:
        lmax_available = min(lmax_available, config.lmax)
    lmin = max(2, config.lmin)
    if lmax_available < lmin:
        return None

    standard_cls = standard.results.get_total_cls(lmax_available, CMB_unit="muK")
    reference_cls = reference.results.get_total_cls(lmax_available, CMB_unit="muK")
    ell = np.arange(lmax_available + 1)
    noise = make_cmb_noise_spectra(ell, config)

    delta_chi2 = 0.0
    worst_ell = lmin
    worst_delta_chi2 = -math.inf
    for multipole in range(lmin, lmax_available + 1):
        model_cov = cmb_covariance_at_l(standard_cls[multipole], noise, multipole, config.fields)
        fiducial_cov = cmb_covariance_at_l(reference_cls[multipole], noise, multipole, config.fields)
        ell_chi2 = gaussian_delta_chi2_at_l(model_cov, fiducial_cov, multipole, config.fsky)
        delta_chi2 += ell_chi2
        if ell_chi2 > worst_delta_chi2:
            worst_delta_chi2 = ell_chi2
            worst_ell = multipole

    return ChiSquaredResult(
        delta_chi2=delta_chi2,
        lmin=lmin,
        lmax=lmax_available,
        fields=config.fields,
        fsky=config.fsky,
        noise_muK_arcmin_T=config.noise_muK_arcmin_T,
        noise_muK_arcmin_P=config.noise_muK_arcmin_P,
        fwhm_arcmin=config.fwhm_arcmin,
        worst_ell=worst_ell,
        worst_delta_chi2=worst_delta_chi2,
    )


def cmb_covariance_at_l(cls: np.ndarray, noise: dict[str, np.ndarray], ell: int, fields: tuple[str, ...]) -> np.ndarray:
    covariance = np.zeros((len(fields), len(fields)))
    for i, field_i in enumerate(fields):
        for j, field_j in enumerate(fields):
            if i == j:
                covariance[i, j] = cls[CL_COLUMNS[field_i + field_i]] + noise[field_i][ell]
            elif {field_i, field_j} == {"T", "E"}:
                covariance[i, j] = cls[CL_COLUMNS["TE"]]
    return covariance


def gaussian_delta_chi2_at_l(
    model_cov: np.ndarray,
    fiducial_cov: np.ndarray,
    ell: int,
    fsky: float,
) -> float:
    sign_model, logdet_model = np.linalg.slogdet(model_cov)
    sign_fiducial, logdet_fiducial = np.linalg.slogdet(fiducial_cov)
    if sign_model <= 0 or sign_fiducial <= 0:
        return math.inf
    trace = np.trace(np.linalg.solve(model_cov, fiducial_cov))
    return float(fsky * (2 * ell + 1) * (trace + logdet_model - logdet_fiducial - model_cov.shape[0]))


def print_chi2_result(result: ChiSquaredResult | None, config_name: str) -> None:
    if result is None:
        print("\nFiducial CMB delta chi-squared: not calculated")
        return
    fields = ",".join(result.fields)
    print("\nFiducial CMB delta chi-squared:")
    print(
        f"  config={config_name}, fields={fields}, L={result.lmin}..{result.lmax}, "
        f"fsky={result.fsky:g}, beam={result.fwhm_arcmin:g} arcmin"
    )
    print(f"  noise T={result.noise_muK_arcmin_T:g} muK-arcmin, P={result.noise_muK_arcmin_P:g} muK-arcmin")
    print(
        f"  delta chi2={result.delta_chi2:.6g}; largest per-L contribution={result.worst_delta_chi2:.6g} at L={result.worst_ell}"
    )


def print_run_summary(standard: RunOutput, reference: RunOutput) -> None:
    print("Runs:")
    for run in (standard, reference):
        print(f"  {run.label}: {run_timing_summary(run)}")
    if reference.cpu_time:
        print(f"  cpu ratio {standard.label}/{reference.label}: {standard.cpu_time / reference.cpu_time:.3g}")
    if reference.wall_time:
        print(f"  wall ratio {standard.label}/{reference.label}: {standard.wall_time / reference.wall_time:.3g}")


def print_table(title: str, rows: list[StatRow]) -> None:
    if not rows:
        print(f"\n{title}: not calculated")
        return
    print(f"\n{title}:")
    print(f"  {'quantity':<20} {'range':<20} {'max abs':>12} {'rms':>12} {'tol':>12} status")
    for row in rows:
        max_abs = "nan" if math.isnan(row.max_abs) else f"{row.max_abs:.4g}"
        rms = "nan" if math.isnan(row.rms) else f"{row.rms:.4g}"
        status = "OK" if row.passed else "FAIL"
        location = f" at {row.location}" if row.location and row.location != "no finite samples" else ""
        print(
            f"  {row.quantity:<20} {row.range_label:<20} {max_abs:>12} {rms:>12} "
            f"{row.tolerance:>12.4g} {status}{location}"
        )


def print_derived_table(rows: list[StatRow]) -> None:
    title = "Derived parameter fractional changes"
    if not rows:
        print(f"\n{title}: not calculated")
        return
    print(f"\n{title}:")
    print(f"  {'quantity':<20} {'frac change':>12} {'tol':>12} status")
    for row in rows:
        value = "nan" if math.isnan(row.max_abs) else f"{row.max_abs:.4g}"
        status = "OK" if row.passed else "FAIL"
        print(f"  {row.quantity:<20} {value:>12} {row.tolerance:>12.4g} {status}")


def print_comparison(comparison: ComparisonResult) -> None:
    print_derived_table(comparison.derived_rows)
    print_table("Lensed CMB spectra errors", comparison.cl_rows)
    print_comparison_notes(comparison)
    print_table("Lensing potential errors", comparison.lensing_rows)
    print_table("Matter power errors", comparison.mpk_rows)


def print_comparison_notes(comparison: ComparisonResult) -> None:
    low_lmax_edge_rows = [row for row in comparison.cl_rows if LOW_LMAX_LENSED_EDGE_LABEL in row.range_label]
    if low_lmax_edge_rows:
        quantities = ", ".join(sorted({row.quantity for row in low_lmax_edge_rows}))
        print(
            "\nNotes:"
            f"\n  Low-lmax lensed CMB edge tolerance {LOW_LMAX_LENSED_EDGE_TOLERANCE:.4g} "
            f"was used for {quantities}; this applies only to the final "
            f"max({LOW_LMAX_LENSED_EDGE_MIN_WIDTH}, ceil({LOW_LMAX_LENSED_EDGE_FRACTION:g}*lmax)) "
            f"multipoles when lmax < {LOW_LMAX_LENSED_EDGE_LMAX}."
        )


def failure_kmax_note(
    standard: RunOutput,
    reference: RunOutput,
    comparison: ComparisonResult,
    *,
    lens_output_margin: int | None = None,
    lens_potential_accuracy: float | None = None,
    reference_lens_output_margin: int | None = None,
    reference_lens_potential_accuracy: float | None = None,
) -> str | None:
    if comparison.passed:
        return None
    if any(
        value is not None
        for value in (
            lens_output_margin,
            lens_potential_accuracy,
            reference_lens_output_margin,
            reference_lens_potential_accuracy,
        )
    ):
        return None

    standard_params = getattr(standard.results, "Params", standard.params)
    reference_params = getattr(reference.results, "Params", reference.params)
    if not current_uses_nonlinear_lensing(standard_params):
        return None

    standard_kmax = float(getattr(standard_params.Transfer, "kmax", 0.0))
    reference_kmax = float(getattr(reference_params.Transfer, "kmax", 0.0))
    if not (reference_kmax > standard_kmax * (1.0 + 1e-12)):
        return None

    return (
        "The boosted nonlinear-lensing reference used a larger nonlinear-source "
        f"Transfer.kmax ({reference_kmax:.6g}) than the standard run ({standard_kmax:.6g}). "
        "Some failing high-L lensing/CMB rows may therefore reflect the changed nonlinear "
        "matter-power k range, not just numerical convergence at fixed kmax."
    )


def failure_summary(comparison: ComparisonResult) -> str:
    worst = comparison.worst_failure
    if worst is None:
        return "FAIL: at least one comparison exceeded its tolerance."
    location = f" at {worst.location}" if worst.location else ""
    return (
        "FAIL: worst failure is "
        f"{worst.quantity} ({worst.range_label}) max abs {worst.max_abs:.4g} > tol {worst.tolerance:.4g}"
        f"{location}."
    )


def comparison_status(comparison: ComparisonResult) -> str:
    score = comparison_score(comparison)
    if comparison.passed:
        return f"PASS score={score:.4g}"
    worst = comparison.worst_failure
    if worst is None:
        return f"FAIL score={score:.4g}"
    location = f" at {worst.location}" if worst.location else ""
    return (
        f"FAIL score={score:.4g}; worst={worst.quantity} {worst.range_label}{location}, "
        f"{worst.max_abs:.4g}/{worst.tolerance:.4g}"
    )


def search_timing_summary(result: SearchResult) -> list[str]:
    lines = []
    if result.run is not None:
        lines.append(f"  selected timing: {run_timing_summary(result.run)}")
    fastest = result.fastest
    if (
        fastest is not None
        and fastest.run is not None
        and settings_key(fastest.settings) != settings_key(result.settings)
    ):
        lines.append(
            "  fastest passing candidate seen: "
            f"{run_timing_summary(fastest.run)} with {format_settings(changed_settings(fastest.settings, result.settings))}"
        )
    return lines


def print_search_timing_summary(result: SearchResult) -> None:
    for line in search_timing_summary(result):
        print(line)


def plot_errors(standard: RunOutput, reference: RunOutput, plot_dir: Path) -> None:
    import matplotlib.pyplot as plt

    plot_dir.mkdir(parents=True, exist_ok=True)
    if standard.lensed_cls is not None and reference.lensed_cls is not None:
        lmax = min(standard.lensed_cls.shape[0], reference.lensed_cls.shape[0]) - 1
        ell = np.arange(lmax + 1)
        std_cls = standard.lensed_cls[: lmax + 1]
        ref_cls = reference.lensed_cls[: lmax + 1]
        plt.figure(figsize=(9, 5))
        for quantity in ("TT", "EE", "BB"):
            errors = fractional_delta(std_cls[:, CL_COLUMNS[quantity]], ref_cls[:, CL_COLUMNS[quantity]])
            plot_log_l_errors(plt, ell, errors, quantity)
        plot_log_l_errors(plt, ell, normalized_te_delta(std_cls, ref_cls), "TE/sqrt(TT*EE)")
        plt.axhline(0, color="black", linewidth=0.8)
        plt.xlabel("L")
        plt.ylabel("fractional error")
        plt.legend()
        plt.tight_layout()
        path = plot_dir / "lensed_cls_errors.png"
        plt.savefig(path, dpi=150)
        plt.close()
        print(f"Wrote {path}")

    if standard.lens_potential_cls is not None and reference.lens_potential_cls is not None:
        lmax = (
            min(
                standard.lens_potential_cls.shape[0],
                reference.lens_potential_cls.shape[0],
            )
            - 1
        )
        ell = np.arange(lmax + 1)
        errors = fractional_delta(
            standard.lens_potential_cls[: lmax + 1, 0],
            reference.lens_potential_cls[: lmax + 1, 0],
        )
        plt.figure(figsize=(9, 5))
        plot_log_l_errors(plt, ell, errors, "PP")
        plt.axhline(0, color="black", linewidth=0.8)
        plt.xlabel("L")
        plt.ylabel("fractional error")
        plt.legend()
        plt.tight_layout()
        path = plot_dir / "lens_potential_errors.png"
        plt.savefig(path, dpi=150)
        plt.close()
        print(f"Wrote {path}")

    if standard.matter_power is not None and reference.matter_power is not None:
        std_k, std_z, std_pk, ref_pk = common_matter_power_grid(standard.matter_power, reference.matter_power)
        errors = fractional_delta(std_pk, ref_pk)
        plt.figure(figsize=(9, 5))
        for index, z in enumerate(std_z):
            plt.semilogx(std_k, errors[index], label=f"z={z:.4g}")
        plt.axhline(0, color="black", linewidth=0.8)
        plt.xlabel("k/h")
        plt.ylabel("fractional error")
        plt.legend()
        plt.tight_layout()
        path = plot_dir / "matter_power_errors.png"
        plt.savefig(path, dpi=150)
        plt.close()
        print(f"Wrote {path}")


def plot_log_l_errors(plt, ell: np.ndarray, errors: np.ndarray, label: str) -> None:
    mask = (ell >= 2) & np.isfinite(errors)
    if np.any(mask):
        plt.semilogx(ell[mask], errors[mask], label=label, alpha=0.85)


def settings_cost(settings: dict[str, float | bool], raw_settings: dict[str, float | bool]) -> tuple[float, int, tuple]:
    changed = 0
    product = 1.0
    values = []
    numeric_keys = sorted(key for key in set(settings) | set(raw_settings) if key != "DoLateRadTruncation")
    for key in numeric_keys:
        value = float(settings[key])
        raw_value = float(raw_settings[key])
        changed += int(not math.isclose(value, raw_value))
        product *= value
        values.append(value)
    if "DoLateRadTruncation" in settings or "DoLateRadTruncation" in raw_settings:
        changed += int(settings.get("DoLateRadTruncation") != raw_settings.get("DoLateRadTruncation"))
        values.append(int(bool(settings.get("DoLateRadTruncation"))))
    return product, changed, tuple(values)


def raw_accuracy_settings(params) -> dict[str, float | bool]:
    return {
        "AccuracyBoost": float(params.Accuracy.AccuracyBoost),
        "lSampleBoost": float(params.Accuracy.lSampleBoost),
        "lAccuracyBoost": float(params.Accuracy.lAccuracyBoost),
        "IntTolBoost": float(params.Accuracy.IntTolBoost),
        "DoLateRadTruncation": bool(params.DoLateRadTruncation),
    }


def component_accuracy_settings(params) -> dict[str, float | bool]:
    settings = {"AccuracyBoost": 1.0}
    for key in component_refinement_keys(params):
        settings[key] = float(getattr(params.Accuracy, key))
    settings["DoLateRadTruncation"] = bool(params.DoLateRadTruncation)
    return settings


def component_refinement_keys(params) -> tuple[str, ...]:
    keys = list(COMPONENT_REFINEMENT_KEYS)
    if abs(float(getattr(params, "omk", 0.0))) <= OMEGA_K_FLAT:
        keys.remove("NonFlatIntAccuracyBoost")
    if not uses_nonlinear_sources(params):
        keys.remove("NonlinSourceBoost")
    if not has_redshift_windows(params):
        keys.remove("KmaxBoost")
        keys.remove("SourceLimberBoost")
    return tuple(keys)


def uses_nonlinear_sources(params) -> bool:
    return str(getattr(params, "NonLinear", "NonLinear_none")) in {
        "NonLinear_pk",
        "NonLinear_both",
    }


def has_redshift_windows(params) -> bool:
    source_windows = getattr(params, "SourceWindows", None)
    if source_windows is not None:
        try:
            return len(source_windows) >= 1
        except TypeError:
            pass
    return int(getattr(params, "num_redshiftwindows", 0)) >= 1


def component_target_settings(
    raw_settings: dict[str, float | bool],
    passing_settings: dict[str, float | bool],
) -> dict[str, float | bool]:
    global_boost = float(passing_settings.get("AccuracyBoost", 1.0))
    target = dict(raw_settings)
    target["AccuracyBoost"] = 1.0
    target["DoLateRadTruncation"] = passing_settings.get(
        "DoLateRadTruncation", raw_settings.get("DoLateRadTruncation", True)
    )
    for key in GLOBAL_BOOST_COMPONENT_KEYS:
        if key not in raw_settings:
            continue
        direct_value = float(passing_settings.get(key, raw_settings[key]))
        boosted_value = float(raw_settings[key]) * global_boost
        target_value = max(
            float(raw_settings[key]),
            direct_value,
            boosted_value,
        )
        target[key] = next_meaningful_accuracy_value(key, float(raw_settings[key]), target_value)
    for key in ("lSampleBoost", "lAccuracyBoost"):
        if key in raw_settings:
            target[key] = max(
                float(raw_settings[key]),
                float(passing_settings.get(key, raw_settings[key])),
            )
    return target


def accuracy_search_target_settings(
    raw_settings: dict[str, float | bool],
    requested_settings: dict[str, float | bool],
) -> dict[str, float | bool]:
    target = dict(requested_settings)
    for key in BASE_ACCURACY_KEYS:
        target[key] = max(float(raw_settings[key]), float(requested_settings[key]))
    target["DoLateRadTruncation"] = bool(requested_settings["DoLateRadTruncation"])
    return target


def next_meaningful_accuracy_value(key: str, raw_value: float, target_value: float) -> float:
    if key not in DISCRETE_ACCURACY_VALUES:
        return target_value
    for value in DISCRETE_ACCURACY_VALUES[key]:
        if value >= raw_value and value >= target_value:
            return value
    return max(raw_value, DISCRETE_ACCURACY_VALUES[key][-1])


def search_candidates(
    raw_settings: dict[str, float | bool],
    high_accuracy_settings: dict[str, float | bool],
    grid: Iterable[float] | None,
) -> list[dict[str, float | bool]]:
    if grid is None:
        return []
    candidate_values = {}
    for key in BASE_ACCURACY_KEYS:
        raw_value = float(raw_settings[key])
        high_value = max(raw_value, float(high_accuracy_settings[key]))
        candidate_values[key] = sorted({max(raw_value, value) for value in grid if value <= high_value} | {high_value})
    late_rad_values = sorted(
        {
            bool(raw_settings["DoLateRadTruncation"]),
            bool(high_accuracy_settings["DoLateRadTruncation"]),
        },
        key=int,
    )
    candidates = [
        {
            "AccuracyBoost": accuracy_boost,
            "lSampleBoost": l_sample_boost,
            "lAccuracyBoost": l_accuracy_boost,
            "IntTolBoost": int_tol_boost,
            "DoLateRadTruncation": do_late_rad_truncation,
        }
        for accuracy_boost, l_sample_boost, l_accuracy_boost, int_tol_boost, do_late_rad_truncation in itertools.product(
            candidate_values["AccuracyBoost"],
            candidate_values["lSampleBoost"],
            candidate_values["lAccuracyBoost"],
            candidate_values["IntTolBoost"],
            late_rad_values,
        )
    ]
    return sorted(candidates, key=lambda settings: settings_cost(settings, raw_settings))


def bracket_then_refine_boost_subset(
    subset: tuple[str, ...],
    raw_settings: dict[str, float | bool],
    target_settings: dict[str, float | bool],
    late_rad_value: bool,
    evaluator,
    tolerance: float,
) -> tuple[dict[str, float | bool], ComparisonResult | None] | None:
    if not any(float(target_settings[key]) > float(raw_settings[key]) for key in subset):
        return None

    settings = subset_target_settings(raw_settings, target_settings, subset, late_rad_value)
    comparison = evaluator(settings)
    if comparison is None:
        return settings, None
    if not comparison.passed:
        return None
    print(f"Found passing subset {subset}; refining within bracket.")
    return refine_numeric_boosts(
        settings,
        raw_settings,
        evaluator,
        tolerance,
        cost_function=settings_cost,
    )


def subset_target_settings(
    raw_settings: dict[str, float | bool],
    target_settings: dict[str, float | bool],
    subset: tuple[str, ...],
    late_rad_value: bool,
) -> dict[str, float | bool]:
    settings = dict(raw_settings)
    settings["DoLateRadTruncation"] = late_rad_value
    for key in subset:
        settings[key] = target_settings[key]
    return settings


def find_minimal_boosts(
    ini_file: Path,
    args: argparse.Namespace,
    reference: RunOutput,
    high_accuracy_settings: dict[str, float | bool],
    raw_comparison: ComparisonResult | None = None,
    raw_run: RunOutput | None = None,
) -> SearchResult | tuple[None, None]:
    raw_params = load_params(
        ini_file,
        no_validate=args.no_validate,
        set_for_lmax=args.set_for_lmax,
        lens_output_margin=args.lens_output_margin,
        lens_potential_accuracy=args.lens_potential_accuracy,
    )
    raw_settings = raw_accuracy_settings(raw_params)
    target_settings = accuracy_search_target_settings(raw_settings, high_accuracy_settings)
    candidates = search_candidates(raw_settings, target_settings, args.search_grid)
    comparisons: dict[tuple, ComparisonResult] = {}
    runs: dict[tuple, RunOutput] = {}
    comparison_runs: dict[int, RunOutput] = {}
    if raw_comparison is not None:
        raw_key = settings_key(raw_settings)
        comparisons[raw_key] = raw_comparison
        if raw_run is not None:
            runs[raw_key] = raw_run
            comparison_runs[id(raw_comparison)] = raw_run
    fastest: PassingCandidate | None = None
    runs_done = 0

    def make_result(settings: dict[str, float | bool], comparison: ComparisonResult) -> SearchResult:
        key = settings_key(refine_discrete_accuracy_settings(settings, raw_settings))
        return SearchResult(
            settings,
            comparison,
            runs.get(key) or comparison_runs.get(id(comparison)),
            fastest,
        )

    def run_candidate(settings: dict[str, float | bool]) -> ComparisonResult | None:
        nonlocal fastest, runs_done
        settings = refine_discrete_accuracy_settings(settings, raw_settings)
        key = settings_key(settings)
        if key in comparisons:
            return comparisons[key]
        if args.max_search_runs is not None and runs_done >= args.max_search_runs:
            return None
        runs_done += 1
        print(f"  [{runs_done}] {format_settings(settings)}")
        candidate = run_case(
            ini_file,
            "candidate",
            no_validate=args.no_validate,
            accuracy_settings=settings,
            lmax=effective_lmax(args),
            set_for_lmax=args.set_for_lmax,
            lens_output_margin=args.lens_output_margin,
            lens_potential_accuracy=args.lens_potential_accuracy,
            mpk_kmin=args.mpk_kmin,
            mpk_npoints=args.mpk_npoints,
        )
        runs[key] = candidate
        print(f"      {run_timing_summary(candidate)}")
        comparison = compare_runs(
            candidate,
            reference,
            derived_tolerance=args.derived_tolerance,
            mpk_tolerance=args.mpk_tolerance,
        )
        comparisons[key] = comparison
        comparison_runs[id(comparison)] = candidate
        if comparison.passed:
            passing = PassingCandidate(settings, comparison, candidate)
            if fastest is None or run_timing_key(candidate) < run_timing_key(fastest.run):
                fastest = passing
        print(f"      {comparison_status(comparison)}")
        return comparison

    raw_comparison = run_candidate(raw_settings)
    if raw_comparison is None:
        return None, None
    if raw_comparison.passed:
        return make_result(raw_settings, raw_comparison)

    numeric_keys = BASE_ACCURACY_KEYS
    late_rad_values = [bool(raw_settings["DoLateRadTruncation"])]
    if target_settings["DoLateRadTruncation"] != raw_settings["DoLateRadTruncation"]:
        late_rad_values.append(bool(target_settings["DoLateRadTruncation"]))
    exhaustive = getattr(args, "exhaustive_boost_search", False) is True
    if target_settings["DoLateRadTruncation"] != raw_settings["DoLateRadTruncation"]:
        settings = dict(raw_settings)
        settings["DoLateRadTruncation"] = target_settings["DoLateRadTruncation"]
        comparison = run_candidate(settings)
        if comparison is None:
            return None, None
        if comparison.passed:
            return make_result(settings, comparison)

    print("\nSearching single boost parameters from values close to raw...")
    best_settings = None
    best_comparison = None
    best_cost = None
    for subset_size in range(1, len(numeric_keys) + 1):
        if subset_size == 2 and best_settings is not None and not exhaustive:
            return make_result(best_settings, best_comparison)
        if subset_size == 2:
            if best_settings is None:
                print("\nNo single-parameter path passed; searching boost parameter combinations...")
            else:
                print("\nExhaustive search requested; searching boost parameter combinations...")
        elif subset_size > 2:
            print(f"\nSearching {subset_size}-parameter boost combinations...")
        for late_rad_value in late_rad_values:
            for subset in itertools.combinations(numeric_keys, subset_size):
                result = bracket_then_refine_boost_subset(
                    subset,
                    raw_settings,
                    target_settings,
                    late_rad_value,
                    run_candidate,
                    args.search_tolerance,
                )
                if result is None:
                    continue
                refined_settings, refined_comparison = result
                if refined_comparison is None:
                    if best_settings is not None:
                        print("Search run budget exhausted; returning the best passing subset found so far.")
                        return make_result(best_settings, best_comparison)
                    print("Search run budget exhausted before finding a passing subset.")
                    return None, None
                cost = settings_cost(refined_settings, raw_settings)
                if best_cost is None or cost < best_cost:
                    best_settings = refined_settings
                    best_comparison = refined_comparison
                    best_cost = cost
        if best_settings is not None and not exhaustive:
            return make_result(best_settings, best_comparison)
    if best_settings is not None:
        return make_result(best_settings, best_comparison)

    if not candidates:
        return None, None

    print(f"\nAdaptive combination search did not pass; trying {len(candidates)} optional seed candidates...")
    for settings in candidates:
        comparison = run_candidate(settings)
        if comparison is None:
            if best_settings is not None:
                print("Search run budget exhausted; returning the best passing settings found so far.")
                return make_result(best_settings, best_comparison)
            print("Search run budget exhausted before finding a passing coarse candidate.")
            return None, None
        if comparison.passed:
            print("Found passing coarse settings; refining numeric boosts.")
            refined_settings, refined_comparison = refine_numeric_boosts(
                settings,
                raw_settings,
                run_candidate,
                args.search_tolerance,
                cost_function=settings_cost,
            )
            return make_result(refined_settings, refined_comparison)
    return None, None


def refine_accuracy_components(
    ini_file: Path,
    args: argparse.Namespace,
    reference: RunOutput,
    reference_accuracy_settings: dict[str, float | bool],
) -> SearchResult | tuple[None, ComparisonResult | None]:
    raw_params = load_params(
        ini_file,
        no_validate=args.no_validate,
        set_for_lmax=args.set_for_lmax,
        lens_output_margin=args.lens_output_margin,
        lens_potential_accuracy=args.lens_potential_accuracy,
    )
    target_accuracy_settings = accuracy_search_target_settings(
        raw_accuracy_settings(raw_params), reference_accuracy_settings
    )
    raw_settings = component_accuracy_settings(raw_params)
    target_settings = component_target_settings(raw_settings, target_accuracy_settings)
    comparisons: dict[tuple, ComparisonResult] = {}
    runs: dict[tuple, RunOutput] = {}
    comparison_runs: dict[int, RunOutput] = {}
    fastest: PassingCandidate | None = None
    runs_done = 0

    def make_result(settings: dict[str, float | bool], comparison: ComparisonResult) -> SearchResult:
        key = settings_key(refine_discrete_accuracy_settings(settings, raw_settings))
        return SearchResult(
            settings,
            comparison,
            runs.get(key) or comparison_runs.get(id(comparison)),
            fastest,
        )

    def run_candidate(settings: dict[str, float | bool]) -> ComparisonResult | None:
        nonlocal fastest, runs_done
        settings = refine_discrete_accuracy_settings(settings, raw_settings)
        key = settings_key(settings)
        if key in comparisons:
            return comparisons[key]
        if args.max_search_runs is not None and runs_done >= args.max_search_runs:
            return None
        runs_done += 1
        print(f"  [component {runs_done}] {format_settings(changed_settings(settings, raw_settings))}")
        candidate = run_case(
            ini_file,
            "component-candidate",
            no_validate=args.no_validate,
            accuracy_settings=settings,
            lmax=effective_lmax(args),
            set_for_lmax=args.set_for_lmax,
            lens_output_margin=args.lens_output_margin,
            lens_potential_accuracy=args.lens_potential_accuracy,
            mpk_kmin=args.mpk_kmin,
            mpk_npoints=args.mpk_npoints,
        )
        runs[key] = candidate
        print(f"      {run_timing_summary(candidate)}")
        comparison = compare_runs(
            candidate,
            reference,
            derived_tolerance=args.derived_tolerance,
            mpk_tolerance=args.mpk_tolerance,
        )
        comparisons[key] = comparison
        comparison_runs[id(comparison)] = candidate
        if comparison.passed:
            passing = PassingCandidate(settings, comparison, candidate)
            if fastest is None or run_timing_key(candidate) < run_timing_key(fastest.run):
                fastest = passing
        print(f"      {comparison_status(comparison)}")
        return comparison

    keys = tuple(
        key for key in raw_settings if key != "DoLateRadTruncation" and target_settings[key] > raw_settings[key]
    )
    max_size = min(args.component_search_size, len(keys))
    print("\nRefining with AccuracyBoost=1 and component accuracy parameters...")
    if not keys:
        return None, None

    all_target = dict(raw_settings)
    all_target["DoLateRadTruncation"] = target_settings["DoLateRadTruncation"]
    for key in keys:
        all_target[key] = target_settings[key]
    print("Checking all component targets first...")
    all_target_comparison = run_candidate(all_target)
    if all_target_comparison is None:
        print("Component search run budget exhausted.")
        return None, None
    if not all_target_comparison.passed:
        print(
            "All component targets with AccuracyBoost=1 still fail; "
            "the passing run likely depends on direct AccuracyBoost effects."
        )
        print(failure_summary(all_target_comparison))
        return None, all_target_comparison

    print("Pruning unnecessary components from the all-target pass...")
    pruned_settings, pruned_comparison = greedy_prune_component_settings(
        all_target,
        raw_settings,
        keys,
        run_candidate,
    )
    if pruned_comparison is None:
        print("Component search run budget exhausted.")
        return None, None
    pruned_keys = component_changed_numeric_keys(pruned_settings, raw_settings)
    if component_settings_cost(pruned_settings, raw_settings) < component_settings_cost(all_target, raw_settings):
        print(f"Greedy pruning reduced component count to {len(pruned_keys)}.")
    if not pruned_keys:
        return make_result(pruned_settings, pruned_comparison)

    best_settings = None
    best_cost = None
    max_size = min(max_size, len(pruned_keys))
    if max_size >= 1:
        print("Checking target-valued subsets within the pruned component set...")
    for subset_size in range(1, max_size + 1):
        for subset in itertools.combinations(pruned_keys, subset_size):
            settings = subset_target_settings(
                raw_settings,
                target_settings,
                subset,
                bool(target_settings["DoLateRadTruncation"]),
            )
            comparison = run_candidate(settings)
            if comparison is None:
                print("Component search run budget exhausted.")
                return None, None
            if not comparison.passed:
                continue
            cost = component_settings_cost(settings, raw_settings)
            if best_cost is None or cost < best_cost:
                best_settings = settings
                best_cost = cost
        if best_settings is not None:
            print("Found passing component subset; refining values.")
            refined_settings, refined_comparison = refine_numeric_boosts(
                best_settings,
                raw_settings,
                run_candidate,
                args.search_tolerance,
                cost_function=component_settings_cost,
            )
            return make_result(refined_settings, refined_comparison)

    print("No smaller target-valued subset passed; refining the pruned component set.")
    refined_settings, refined_comparison = refine_numeric_boosts(
        pruned_settings,
        raw_settings,
        run_candidate,
        args.search_tolerance,
        cost_function=component_settings_cost,
    )
    return make_result(refined_settings, refined_comparison)


def greedy_prune_component_settings(
    settings: dict[str, float | bool],
    raw_settings: dict[str, float | bool],
    keys: tuple[str, ...],
    evaluator,
) -> tuple[dict[str, float | bool], ComparisonResult | None]:
    current = dict(settings)
    current_comparison = evaluator(current)
    if current_comparison is None or not current_comparison.passed:
        return current, current_comparison

    removable_keys = [
        key for key in keys if key in current and not math.isclose(float(current[key]), float(raw_settings[key]))
    ]
    improved = True
    while improved:
        improved = False
        best_trial = None
        best_comparison = None
        best_cost = None
        for key in removable_keys:
            if math.isclose(float(current[key]), float(raw_settings[key])):
                continue
            trial = dict(current)
            trial[key] = raw_settings[key]
            comparison = evaluator(trial)
            if comparison is None:
                return current, None
            if not comparison.passed:
                continue
            cost = component_settings_cost(trial, raw_settings)
            if best_cost is None or cost < best_cost:
                best_trial = trial
                best_comparison = comparison
                best_cost = cost
        if best_trial is not None:
            current = best_trial
            current_comparison = best_comparison
            improved = True
    return current, current_comparison


def component_changed_numeric_keys(
    settings: dict[str, float | bool], raw_settings: dict[str, float | bool]
) -> tuple[str, ...]:
    return tuple(
        key
        for key, value in settings.items()
        if key not in {"AccuracyBoost", "DoLateRadTruncation"}
        and not math.isclose(float(value), float(raw_settings[key]))
    )


def changed_settings(
    settings: dict[str, float | bool], raw_settings: dict[str, float | bool]
) -> dict[str, float | bool]:
    return {
        key: value
        for key, value in settings.items()
        if key == "AccuracyBoost"
        or key == "DoLateRadTruncation"
        or not math.isclose(float(value), float(raw_settings[key]))
    }


def component_settings_cost(settings: dict[str, float | bool], raw_settings: dict[str, float | bool]) -> tuple:
    changed = changed_settings(settings, raw_settings)
    numeric = [
        float(value) / max(float(raw_settings[key]), 1e-30)
        for key, value in changed.items()
        if key != "DoLateRadTruncation"
    ]
    return len(numeric), math.prod(numeric) if numeric else 1.0, tuple(sorted(changed))


def refine_numeric_boosts(
    settings: dict[str, float | bool],
    raw_settings: dict[str, float | bool],
    evaluator,
    tolerance: float,
    cost_function=None,
) -> tuple[dict[str, float | bool], ComparisonResult]:
    cost_function = settings_cost if cost_function is None else cost_function
    initial_settings = refine_discrete_accuracy_settings(settings, raw_settings)
    initial_comparison = evaluator(initial_settings)
    if initial_comparison is None or not initial_comparison.passed:
        raise RuntimeError("refine_numeric_boosts requires an initially passing settings point")

    numeric_keys = tuple(
        key
        for key in initial_settings
        if key != "DoLateRadTruncation" and float(initial_settings[key]) - float(raw_settings[key]) > tolerance
    )
    best_settings, best_comparison = maybe_restore_late_rad_truncation(
        initial_settings, raw_settings, evaluator, initial_comparison
    )

    if len(numeric_keys) == 1:
        key = numeric_keys[0]
        best_settings, best_comparison = refine_single_numeric_boost(
            best_settings,
            raw_settings,
            evaluator,
            tolerance,
            key,
            best_comparison,
        )
        best_settings, best_comparison = maybe_restore_late_rad_truncation(
            best_settings,
            raw_settings,
            evaluator,
            best_comparison,
        )
        return rounded_numeric_settings(best_settings), best_comparison

    starts: list[tuple[dict[str, float | bool], ComparisonResult]] = [(best_settings, best_comparison)]
    balanced_settings, balanced_comparison = refine_balanced_path(
        initial_settings, raw_settings, evaluator, tolerance, numeric_keys
    )
    if balanced_comparison is not None:
        starts.append(
            maybe_restore_late_rad_truncation(balanced_settings, raw_settings, evaluator, balanced_comparison)
        )

    key_orders = refinement_key_orders(numeric_keys)
    for start_settings, start_comparison in starts:
        for key_order in key_orders:
            refined_settings, refined_comparison = coordinate_refine_numeric_boosts(
                start_settings,
                raw_settings,
                evaluator,
                tolerance,
                key_order,
                start_comparison,
            )
            if cost_function(refined_settings, raw_settings) < cost_function(best_settings, raw_settings):
                best_settings = refined_settings
                best_comparison = refined_comparison

    return rounded_numeric_settings(best_settings), best_comparison


def maybe_restore_late_rad_truncation(
    settings: dict[str, float | bool],
    raw_settings: dict[str, float | bool],
    evaluator,
    comparison: ComparisonResult,
) -> tuple[dict[str, float | bool], ComparisonResult]:
    current = dict(settings)
    current_comparison = comparison
    if current.get("DoLateRadTruncation") != raw_settings.get("DoLateRadTruncation"):
        trial = dict(current)
        trial["DoLateRadTruncation"] = raw_settings["DoLateRadTruncation"]
        comparison = evaluator(trial)
        if comparison is not None and comparison.passed:
            current = trial
            current_comparison = comparison
    return current, current_comparison


def refine_balanced_path(
    settings: dict[str, float | bool],
    raw_settings: dict[str, float | bool],
    evaluator,
    tolerance: float,
    numeric_keys: tuple[str, ...],
) -> tuple[dict[str, float | bool], ComparisonResult | None]:
    continuous_keys = tuple(key for key in numeric_keys if key not in DISCRETE_ACCURACY_VALUES)
    if len(continuous_keys) != len(numeric_keys):
        settings = refine_discrete_accuracy_settings(settings, raw_settings)
    if not numeric_keys:
        comparison = evaluator(settings)
        return dict(settings), comparison
    if not continuous_keys:
        comparison = evaluator(settings)
        return dict(settings), comparison

    span = max(float(settings[key]) - float(raw_settings[key]) for key in continuous_keys)
    if span <= tolerance:
        comparison = evaluator(settings)
        return dict(settings), comparison

    low = 0.0
    high = 1.0
    current = dict(settings)
    current_comparison = evaluator(current)
    while (high - low) * span > tolerance:
        mid = (low + high) / 2
        trial = dict(settings)
        for key in continuous_keys:
            trial[key] = float(raw_settings[key]) + mid * (float(settings[key]) - float(raw_settings[key]))
        comparison = evaluator(trial)
        if comparison is None:
            break
        if comparison.passed:
            current = trial
            current_comparison = comparison
            high = mid
        else:
            low = mid
    return current, current_comparison


def refinement_key_orders(numeric_keys: tuple[str, ...]) -> list[tuple[str, ...]]:
    if not numeric_keys:
        return [()]
    orders = [numeric_keys, tuple(reversed(numeric_keys))]
    orders.extend(numeric_keys[index:] + numeric_keys[:index] for index in range(1, len(numeric_keys)))
    return list(dict.fromkeys(orders))


def coordinate_refine_numeric_boosts(
    settings: dict[str, float | bool],
    raw_settings: dict[str, float | bool],
    evaluator,
    tolerance: float,
    key_order: tuple[str, ...],
    comparison: ComparisonResult,
) -> tuple[dict[str, float | bool], ComparisonResult]:
    current = dict(settings)
    current_comparison = comparison
    improved = True
    while improved:
        improved = False
        for key in key_order:
            before = float(current[key])
            current, current_comparison = refine_single_numeric_boost(
                current, raw_settings, evaluator, tolerance, key, current_comparison
            )
            improved = improved or float(current[key]) < before - tolerance / 2
    return current, current_comparison


def refine_single_numeric_boost(
    settings: dict[str, float | bool],
    raw_settings: dict[str, float | bool],
    evaluator,
    tolerance: float,
    key: str,
    comparison: ComparisonResult,
) -> tuple[dict[str, float | bool], ComparisonResult]:
    current = dict(settings)
    current_comparison = comparison
    if key in DISCRETE_ACCURACY_VALUES:
        return refine_discrete_numeric_boost(current, raw_settings, evaluator, key, current_comparison)
    low = float(raw_settings[key])
    high = float(current[key])
    if high - low <= tolerance:
        return current, current_comparison
    while high - low > tolerance:
        mid = (low + high) / 2
        trial = dict(current)
        trial[key] = mid
        trial_comparison = evaluator(trial)
        if trial_comparison is None:
            break
        if trial_comparison.passed:
            current = trial
            current_comparison = trial_comparison
            high = mid
        else:
            low = mid
    return current, current_comparison


def refine_discrete_accuracy_settings(settings: dict[str, float | bool], raw_settings: dict[str, float | bool]):
    refined = dict(settings)
    for key in DISCRETE_ACCURACY_VALUES:
        if key in refined and key in raw_settings:
            refined[key] = next_meaningful_accuracy_value(key, float(raw_settings[key]), float(refined[key]))
    return refined


def refine_discrete_numeric_boost(
    settings: dict[str, float | bool],
    raw_settings: dict[str, float | bool],
    evaluator,
    key: str,
    comparison: ComparisonResult,
) -> tuple[dict[str, float | bool], ComparisonResult]:
    current = dict(settings)
    current_comparison = comparison
    candidates = discrete_accuracy_candidates(key, float(raw_settings[key]), float(current[key]))
    for value in candidates:
        trial = dict(current)
        trial[key] = value
        trial_comparison = evaluator(trial)
        if trial_comparison is not None and trial_comparison.passed:
            current = trial
            current_comparison = trial_comparison
            break
    return current, current_comparison


def discrete_accuracy_candidates(key: str, raw_value: float, high_value: float) -> tuple[float, ...]:
    values = [raw_value]
    values.extend(value for value in DISCRETE_ACCURACY_VALUES[key] if raw_value < value <= high_value)
    if high_value not in values:
        values.append(high_value)
    return tuple(dict.fromkeys(values))


def rounded_numeric_settings(
    settings: dict[str, float | bool],
) -> dict[str, float | bool]:
    rounded = dict(settings)
    for key, value in rounded.items():
        if key != "DoLateRadTruncation":
            rounded[key] = round(float(value), 6)
    return rounded


def settings_key(settings: dict[str, float | bool]) -> tuple:
    return tuple(
        (key, bool(value) if key == "DoLateRadTruncation" else round(float(value), 8))
        for key, value in sorted(settings.items())
    )


def format_settings(settings: dict[str, float | bool]) -> str:
    preferred_order = (
        "AccuracyBoost",
        "lSampleBoost",
        "lAccuracyBoost",
        "IntTolBoost",
        "DoLateRadTruncation",
    )
    ordered_keys = [key for key in preferred_order if key in settings] + sorted(
        key for key in settings if key not in preferred_order
    )
    return ", ".join(f"{key}={settings[key]}" for key in ordered_keys)


def effective_lmax(args: argparse.Namespace) -> int | None:
    return args.lmax or args.set_for_lmax


def main(argv: list[str] | None = None, *, prog: str | None = None) -> int:
    """Run the accuracy checker command-line interface."""
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)
    ini_file = Path(args.ini_file)
    high_accuracy_settings = requested_accuracy_settings(args)

    print(f"Input ini: {ini_file}")
    print(f"High-accuracy settings: {format_settings(high_accuracy_settings)}")

    base_params = load_params(ini_file, no_validate=args.no_validate)
    if args.mpk_tolerance is None:
        args.mpk_tolerance = default_mpk_tolerance(base_params)
    if base_params.DoLensing and base_params.WantCls and base_params.Want_CMB:
        if args.reference_lens_output_margin is not None:
            print(
                "Note: --reference-lens-output-margin probes support sensitivity near the output cutoff, but it is not by itself "
                "a truth test unless the high-l reference spectra are also well converged."
            )
    noise_config = None
    if args.chi2:
        try:
            noise_config = chi2_noise_config(args)
        except ValueError as exc:
            parser.error(str(exc))
    result = compare_params_accuracy(
        base_params,
        reference_accuracy_settings=high_accuracy_settings,
        lmax=effective_lmax(args),
        set_for_lmax=args.set_for_lmax,
        lens_output_margin=args.lens_output_margin,
        lens_potential_accuracy=args.lens_potential_accuracy,
        reference_lens_output_margin=args.reference_lens_output_margin,
        reference_lens_potential_accuracy=args.reference_lens_potential_accuracy,
        mpk_kmin=args.mpk_kmin,
        mpk_npoints=args.mpk_npoints,
        derived_tolerance=args.derived_tolerance,
        mpk_tolerance=args.mpk_tolerance,
        chi2_config=noise_config,
        plot_dir=args.plot_dir,
    )
    standard = result.standard
    reference = result.reference
    comparison = result.comparison

    print_run_summary(standard, reference)
    print_comparison(comparison)

    if args.chi2:
        print_chi2_result(result.chi2, noise_config.name)

    if args.find_minimal_boosts:
        search_result = find_minimal_boosts(
            ini_file,
            args,
            reference,
            high_accuracy_settings,
            raw_comparison=comparison,
            raw_run=standard,
        )
        settings, search_comparison = search_result
        if settings is None:
            print("\nNo passing candidate settings found.")
            note = failure_kmax_note(
                standard,
                reference,
                comparison,
                lens_output_margin=args.lens_output_margin,
                lens_potential_accuracy=args.lens_potential_accuracy,
                reference_lens_output_margin=args.reference_lens_output_margin,
                reference_lens_potential_accuracy=args.reference_lens_potential_accuracy,
            )
            if note:
                print(f"\nNote: {note}")
            print(f"\n{failure_summary(comparison)}")
            return 1
        else:
            print(f"\nMinimal passing settings from search: {format_settings(settings)}")
            if isinstance(search_result, SearchResult):
                print_search_timing_summary(search_result)
            if search_comparison:
                print_comparison(search_comparison)
            if args.refine_accuracy_components:
                component_result = refine_accuracy_components(ini_file, args, reference, high_accuracy_settings)
                component_settings, component_comparison = component_result
                if component_settings is None:
                    print("\nNo passing AccuracyBoost=1 component settings found.")
                else:
                    print(
                        "\nPassing AccuracyBoost=1 component settings: "
                        f"{format_settings(changed_settings(component_settings, component_accuracy_settings(base_params)))}"
                    )
                    if isinstance(component_result, SearchResult):
                        print_search_timing_summary(component_result)
                    if component_comparison:
                        print_comparison(component_comparison)
            return 0

    if comparison.passed:
        print("\nPASS: standard results are stable against boosted accuracy settings.")
        return 0
    note = failure_kmax_note(
        standard,
        reference,
        comparison,
        lens_output_margin=args.lens_output_margin,
        lens_potential_accuracy=args.lens_potential_accuracy,
        reference_lens_output_margin=args.reference_lens_output_margin,
        reference_lens_potential_accuracy=args.reference_lens_potential_accuracy,
    )
    if note:
        print(f"\nNote: {note}")
    print(f"\n{failure_summary(comparison)}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
