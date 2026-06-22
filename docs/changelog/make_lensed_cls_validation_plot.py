"""Make a lensed-CMB validation plot for the CAMB2 paper draft.

The plotted residual metric is the pointwise envelope

    max(|Delta TT|/TT, |Delta EE|/EE, |Delta TE|/sqrt(TT EE)).

The full CLASS ``cl_ref.pre`` lensed run is very expensive at L=4000, so the
default reference-precision CLASS curve is generated only to ``--class-ref-lmax``
and masked above that multipole. Use ``--class-ref-lmax 4000`` for an overnight
full-reference version.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

import camb
from camb import model
from camb.check_accuracy import DEFAULT_NOISE_CONFIGS, NoiseConfig, make_cmb_noise_spectra

try:
    from classy import Class
except ImportError as exc:  # pragma: no cover - exercised only when CLASS is absent
    raise SystemExit(
        "Could not import classy. Set PYTHONPATH to the local CLASS wrapper, e.g. "
        "PYTHONPATH=/tmp/classy-comparisons/classy-local:$PYTHONPATH"
    ) from exc


DEFAULT_OUTDIR = Path("docs/changelog/lensed_cls_validation")
DEFAULT_CLASS_WORKDIR = Path(os.environ.get("CLASS_WORKDIR", os.environ.get("TMPDIR", "/tmp") + "/classy-comparisons"))
CLASS_PRECISION_IGNORED_KEYS = {"recfast_Nz0", "l_max_ur_ten"}
SPECTRUM_COLUMNS = {"TT": 0, "EE": 1, "BB": 2, "TE": 3}
CLASS_PRECISION_PROFILES = ("cl_permille", "cl_targeted", "cl_targeted_tight", "cl_ref")
RECOMBINATION_TAG = "planck_recfast"
CAMB_RECFAST_APPROX_MODEL = "planck"
CLASS_PLANCK_RECFAST_PARAMS = {
    "recombination": "RECFAST",
    "recfast_photoion_dependence": "Tmat",
    "recfast_Heswitch": 6,
    "recfast_fudge_He": 0.86,
    "recfast_Hswitch": 1,
    "recfast_fudge_H": 1.14,
    "recfast_delta_fudge_H": -0.015,
    "recfast_AGauss1": -0.14,
    "recfast_AGauss2": 0.079,
    "recfast_zGauss1": 7.28,
    "recfast_zGauss2": 6.73,
    "recfast_wGauss1": 0.18,
    "recfast_wGauss2": 0.33,
}
COSMOLOGY = {
    "H0": 67.5,
    "ombh2": 0.02237,
    "omch2": 0.12,
    "As": 2.1e-9,
    "ns": 0.965,
    "tau": 0.0544,
    "YHe": 0.245,
    "TCMB": 2.7255,
    "nnu": 3.044,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmax", type=int, default=4000)
    parser.add_argument("--class-ref-lmax", type=int, default=500)
    parser.add_argument(
        "--class-ref-profile",
        choices=CLASS_PRECISION_PROFILES,
        default="cl_targeted",
        help=(
            "CLASS precision profile for the higher-accuracy comparison curve. "
            "cl_targeted and cl_targeted_tight are moderate CMB/lensing profiles between "
            "cl_permille.pre and cl_ref.pre."
        ),
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--class-workdir", type=Path, default=DEFAULT_CLASS_WORKDIR)
    parser.add_argument("--reference-lens-potential-accuracy", type=float, default=8.0)
    parser.add_argument("--reference-lens-output-margin", type=int, default=1000)
    parser.add_argument(
        "--class-lensing-margin",
        type=int,
        default=2000,
        help=(
            "Fallback CLASS lensed-spectrum margin above the plotted/reference maximum. "
            "Use the more specific CLASS margin options below for the plotted curves."
        ),
    )
    parser.add_argument(
        "--class-permille-lensing-margin",
        type=int,
        default=4000,
        help=(
            "Compute the plotted CLASS cl_permille curve to this many multipoles above the plotted maximum "
            "before truncating."
        ),
    )
    parser.add_argument(
        "--class-ref-lensing-margin",
        type=int,
        help=(
            "Compute the higher-accuracy CLASS reference curve to this many multipoles above --class-ref-lmax. "
            "Defaults to --class-lensing-margin."
        ),
    )
    parser.add_argument(
        "--class-kmax-mode",
        choices=("boosted-camb", "precision-file"),
        default="boosted-camb",
        help=(
            "For CLASS runs, either override k_max_tau0_over_l_max to match the boosted CAMB "
            "max_eta_k/lmax, or leave the precision file value unchanged"
        ),
    )
    parser.add_argument(
        "--class-k-max-tau0-over-l-max",
        type=float,
        help="Explicit CLASS k_max_tau0_over_l_max override; takes precedence over --class-kmax-mode",
    )
    parser.add_argument("--feedback-logt-agn", type=float, default=7.8)
    parser.add_argument("--force", action="store_true", help="rerun cases even if cached data exist")
    return parser


def find_class_root(class_workdir: Path) -> Path:
    roots = [path for path in class_workdir.glob("class_public-[0-9]*") if path.is_dir()]
    if not roots:
        raise FileNotFoundError(f"No CLASS source tree found under {class_workdir}")
    return sorted(roots)[-1]


def load_class_precision_file(path: Path) -> dict[str, float | int | str]:
    overrides: dict[str, float | int | str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if not value or key in CLASS_PRECISION_IGNORED_KEYS:
            continue
        try:
            parsed: float | int | str = float(value)
            if parsed.is_integer() and "." not in value and "e" not in value.lower():
                parsed = int(parsed)
        except ValueError:
            parsed = value
        overrides[key] = parsed
    return overrides


def class_precision_profiles(class_root: Path) -> dict[str, dict[str, float | int | str]]:
    cl_permille = load_class_precision_file(class_root / "cl_permille.pre")
    cl_ref = load_class_precision_file(class_root / "cl_ref.pre")
    cl_targeted = {
        **cl_permille,
        # Tighter k and perturbation sampling, but much lighter than cl_ref.pre.
        "k_step_sub": 0.025,
        "k_step_super": 0.0005,
        "tol_perturbations_integration": 3.0e-6,
        "perturbations_sampling_stepsize": 0.03,
        # Denser transfer/Bessel sampling at high l.
        "l_logstep": 1.05,
        "l_linstep": 25,
        # Use the more accurate lensed-Cl angular quadrature without the full cl_ref source settings.
        "l_switch_limber": 40.0,
        "accurate_lensing": 1,
        "num_mu_minus_lmax": 400,
        "delta_l_max": 1000,
        # Keep scalar sources to higher k before truncating tails.
        "transfer_neglect_delta_k_S_t0": 5.0,
        "transfer_neglect_delta_k_S_t1": 5.0,
        "transfer_neglect_delta_k_S_t2": 5.0,
        "transfer_neglect_delta_k_S_e": 5.0,
    }
    cl_targeted_tight = {
        **cl_targeted,
        # Ref-like scalar sampling/integration, while keeping the all-purpose cl_ref source tail lighter.
        "k_step_sub": 0.015,
        "k_step_super": 0.0002,
        "tol_perturbations_integration": 1.0e-6,
        "perturbations_sampling_stepsize": 0.015,
        "l_logstep": 1.035,
        "l_linstep": 15,
        "l_max_g": 35,
        "l_max_pol_g": 22,
        "l_max_ur": 40,
        "radiation_streaming_approximation": 2,
        "radiation_streaming_trigger_tau_over_tau_k": 160.0,
        "radiation_streaming_trigger_tau_c_over_tau": 60.0,
        "ur_fluid_approximation": 2,
        "ur_fluid_trigger_tau_over_tau_k": 45.0,
        "hyper_sampling_flat": 10.0,
        "hyper_flat_approximation_nu": 1.0e6,
        "num_mu_minus_lmax": 800,
        "transfer_neglect_delta_k_S_t0": 20.0,
        "transfer_neglect_delta_k_S_t1": 20.0,
        "transfer_neglect_delta_k_S_t2": 20.0,
        "transfer_neglect_delta_k_S_e": 20.0,
        "neglect_CMB_sources_below_visibility": 1.0e-20,
        "transfer_neglect_late_source": 1500.0,
    }
    return {
        "cl_permille": cl_permille,
        "cl_targeted": cl_targeted,
        "cl_targeted_tight": cl_targeted_tight,
        "cl_ref": cl_ref,
    }


def class_profile_label(profile: str) -> str:
    labels = {
        "cl_permille": "CLASS cl_permille",
        "cl_targeted": "CLASS cl_permille tuned",
        "cl_targeted_tight": "CLASS cl_permille boosted",
        "cl_ref": "CLASS cl_ref",
    }
    return labels[profile]


def make_camb_params(
    *,
    lmax: int,
    lens_potential_accuracy: float | None,
    lens_output_margin: int,
    boosted: bool,
    feedback_logt_agn: float | None = None,
) -> camb.CAMBparams:
    params = camb.CAMBparams()
    params.set_classes(recombination_model="Recfast")
    params.Recomb.set_params(recfast_approx_model=CAMB_RECFAST_APPROX_MODEL)
    params.set_cosmology(
        H0=COSMOLOGY["H0"],
        ombh2=COSMOLOGY["ombh2"],
        omch2=COSMOLOGY["omch2"],
        mnu=0.0,
        num_massive_neutrinos=0,
        nnu=COSMOLOGY["nnu"],
        YHe=COSMOLOGY["YHe"],
        tau=COSMOLOGY["tau"],
        TCMB=COSMOLOGY["TCMB"],
    )
    params.InitPower.set_params(As=COSMOLOGY["As"], ns=COSMOLOGY["ns"])
    params.NonLinear = model.NonLinear_both
    if feedback_logt_agn is None:
        params.NonLinearModel.set_params(halofit_version="mead2020")
    else:
        params.NonLinearModel.set_params(halofit_version="mead2020_feedback", HMCode_logT_AGN=feedback_logt_agn)
    params.set_for_lmax(
        lmax,
        lens_potential_accuracy=lens_potential_accuracy,
        lens_output_margin=lens_output_margin,
        nonlinear=True,
    )
    if boosted:
        params.set_accuracy(
            AccuracyBoost=2.0,
            lSampleBoost=2.0,
            lAccuracyBoost=2.0,
            min_l_logl_sampling=100000,
        )
        params.Accuracy.IntTolBoost = 2.0
    params.WantCls = True
    params.Want_CMB = True
    params.WantTransfer = False
    params.WantTensors = False
    params.WantDerivedParameters = False
    return params


def run_camb_case(params: camb.CAMBparams, lmax: int) -> dict[str, Any]:
    start_cpu = time.process_time()
    start_wall = time.perf_counter()
    results = camb.get_results(params)
    cpu_s = time.process_time() - start_cpu
    wall_s = time.perf_counter() - start_wall
    return {
        "dls": results.get_lensed_scalar_cls(lmax=lmax, CMB_unit="muK"),
        "raw_cls": results.get_lensed_scalar_cls(lmax=lmax, CMB_unit="muK", raw_cl=True),
        "cpu_s": cpu_s,
        "wall_s": wall_s,
        "max_l": int(params.max_l),
        "max_eta_k": float(params.max_eta_k),
    }


def class_lensed_dls(params: dict[str, float | int | str], lmax: int) -> tuple[np.ndarray, np.ndarray, float, float]:
    cosmo = Class()
    cosmo.set(params)
    start_cpu = time.process_time()
    start_wall = time.perf_counter()
    try:
        cosmo.compute()
        cls = cosmo.lensed_cl(lmax)
    finally:
        cpu_s = time.process_time() - start_cpu
        wall_s = time.perf_counter() - start_wall
        cosmo.struct_cleanup()
        cosmo.empty()

    ell = np.arange(lmax + 1, dtype=np.float64)
    factor = ell * (ell + 1.0) / (2.0 * math.pi) * (COSMOLOGY["TCMB"] * 1.0e6) ** 2
    dls = np.zeros((lmax + 1, 4), dtype=np.float64)
    raw_cls = np.zeros((lmax + 1, 4), dtype=np.float64)
    for name, column in SPECTRUM_COLUMNS.items():
        key = name.lower()
        if key in cls:
            raw_cls[:, column] = cls[key][: lmax + 1] * (COSMOLOGY["TCMB"] * 1.0e6) ** 2
            dls[:, column] = cls[key][: lmax + 1] * factor
    return dls, raw_cls, cpu_s, wall_s


def make_class_params(
    *,
    lmax: int,
    precision_overrides: dict[str, float | int | str],
) -> dict[str, float | int | str]:
    params: dict[str, float | int | str] = {
        "output": "tCl,pCl,lCl",
        "lensing": "yes",
        "non linear": "hmcode",
        "hmcode_version": "2020",
        "l_max_scalars": lmax,
        "H0": COSMOLOGY["H0"],
        "omega_b": COSMOLOGY["ombh2"],
        "omega_cdm": COSMOLOGY["omch2"],
        "Omega_k": 0.0,
        "A_s": COSMOLOGY["As"],
        "n_s": COSMOLOGY["ns"],
        "tau_reio": COSMOLOGY["tau"],
        "YHe": COSMOLOGY["YHe"],
        "T_cmb": COSMOLOGY["TCMB"],
        "N_ur": COSMOLOGY["nnu"],
        "N_ncdm": 0,
    }
    params.update(CLASS_PLANCK_RECFAST_PARAMS)
    params.update(precision_overrides)
    params.update(CLASS_PLANCK_RECFAST_PARAMS)
    params["l_max_scalars"] = lmax
    return params


def run_class_case(
    *,
    lmax: int,
    precision_overrides: dict[str, float | int | str],
    k_max_tau0_over_l_max: float | None,
) -> dict[str, Any]:
    params = make_class_params(lmax=lmax, precision_overrides=precision_overrides)
    if k_max_tau0_over_l_max is not None:
        params["k_max_tau0_over_l_max"] = k_max_tau0_over_l_max
    dls, raw_cls, cpu_s, wall_s = class_lensed_dls(params, lmax)
    return {
        "dls": dls,
        "raw_cls": raw_cls,
        "cpu_s": cpu_s,
        "wall_s": wall_s,
        "class_params": params,
    }


def save_case(path: Path, payload: dict[str, Any]) -> None:
    np.savez_compressed(
        path,
        dls=payload["dls"],
        raw_cls=payload["raw_cls"],
        metadata=np.array(json.dumps({k: v for k, v in payload.items() if k not in {"dls", "raw_cls"}})),
    )


def load_case(path: Path) -> dict[str, Any]:
    data = np.load(path, allow_pickle=False)
    metadata = json.loads(str(data["metadata"]))
    return {"dls": data["dls"], "raw_cls": data["raw_cls"], **metadata}


def cached_case(path: Path, force: bool, runner: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    if path.exists() and not force:
        print(f"Using cached {path}")
        return load_case(path)
    print(f"Running {path.stem}...")
    payload = runner()
    save_case(path, payload)
    return payload


def kmax_tag(k_max_tau0_over_l_max: float | None) -> str:
    return "precision" if k_max_tau0_over_l_max is None else f"kmaxtau{k_max_tau0_over_l_max:.6g}".replace(".", "p")


def float_tag(value: float) -> str:
    return f"{value:.6g}".replace("-", "m").replace(".", "p")


def residual_components(candidate: np.ndarray, reference: np.ndarray) -> dict[str, np.ndarray]:
    with np.errstate(divide="ignore", invalid="ignore"):
        tt = (candidate[:, 0] - reference[:, 0]) / reference[:, 0]
        ee = (candidate[:, 1] - reference[:, 1]) / reference[:, 1]
        te = (candidate[:, 3] - reference[:, 3]) / np.sqrt(np.maximum(reference[:, 0] * reference[:, 1], 0.0))
    return {"TT": tt, "EE": ee, "TE": te}


def residual_envelope(candidate: np.ndarray, reference: np.ndarray, lmax: int) -> np.ndarray:
    envelope = np.full(lmax + 1, np.nan)
    max_l = min(lmax, candidate.shape[0] - 1, reference.shape[0] - 1)
    components = residual_components(candidate[: max_l + 1], reference[: max_l + 1])
    stacked = np.vstack([np.abs(values) for values in components.values()])
    finite_columns = np.any(np.isfinite(stacked), axis=0)
    envelope[: max_l + 1][finite_columns] = np.nanmax(stacked[:, finite_columns], axis=0)
    envelope[:2] = np.nan
    return envelope


def mask_above(values: np.ndarray, lmax: int) -> np.ndarray:
    masked = values.copy()
    if masked.shape[0] > lmax + 1:
        masked[lmax + 1 :] = np.nan
    return masked


def amplitude_fisher_per_ell(dls: np.ndarray, config: NoiseConfig, lmax: int) -> np.ndarray:
    ell = np.arange(lmax + 1)
    noise = make_cmb_noise_spectra(ell, config)
    fisher = np.zeros(lmax + 1)
    for multipole in range(max(2, config.lmin), lmax + 1):
        signal = np.array(
            [
                [dls[multipole, 0], dls[multipole, 3]],
                [dls[multipole, 3], dls[multipole, 1]],
            ]
        )
        total = signal + np.diag([noise["T"][multipole], noise["E"][multipole]])
        try:
            inv_total = np.linalg.inv(total)
        except np.linalg.LinAlgError:
            continue
        fisher[multipole] = config.fsky * (2 * multipole + 1) * np.trace(inv_total @ signal @ inv_total @ signal) / 2.0
    return fisher


def binned_amplitude_sigma(
    dls: np.ndarray, config: NoiseConfig, lmax: int, fractional_bin_width: float = 0.5
) -> np.ndarray:
    fisher_per_ell = amplitude_fisher_per_ell(dls, config, lmax)
    cumulative = np.concatenate(([0.0], np.cumsum(fisher_per_ell)))
    fisher = np.zeros(lmax + 1)
    for multipole in range(max(2, config.lmin), lmax + 1):
        half_width = fractional_bin_width * multipole / 2.0
        lo = max(config.lmin, int(math.ceil(multipole - half_width)))
        hi = min(lmax, int(math.floor(multipole + half_width)))
        if hi >= lo:
            fisher[multipole] = cumulative[hi + 1] - cumulative[lo]
    sigma = np.full(lmax + 1, np.nan)
    mask = fisher > 0
    sigma[mask] = 1.0 / np.sqrt(fisher[mask])
    return sigma


def summarize_curve(values: np.ndarray, ranges: dict[str, tuple[int, int]]) -> dict[str, dict[str, float | int | None]]:
    summary = {}
    for name, (lo, hi) in ranges.items():
        hi = min(hi, values.shape[0] - 1)
        segment = values[lo : hi + 1]
        if not np.any(np.isfinite(segment)):
            summary[name] = {"max": None, "ell": -1}
            continue
        offset = int(np.nanargmax(segment))
        summary[name] = {"max": float(segment[offset]), "ell": lo + offset}
    return summary


def make_plot(
    path: Path,
    ell: np.ndarray,
    curves: dict[str, np.ndarray],
    sigma_a: np.ndarray,
    summaries: dict[str, dict[str, dict[str, float | int]]],
) -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    def color_for(label: str) -> str:
        if label.startswith("CAMB default"):
            return "#111111"
        if label.startswith("CLASS cl_permille"):
            return "#2f8f46" if "boosted /" in label or "tuned /" in label else "#1b6ca8"
        if label.startswith("CLASS cl_ref"):
            return "#5b2a86"
        if label.startswith("HMCode feedback"):
            return "#b23a48"
        if "boosted" in label or "tuned" in label:
            return "#2f8f46"
        return "0.2"

    def linestyle_for(label: str):
        if label.startswith("CLASS cl_permille") and "boosted /" not in label and "tuned /" not in label:
            return "-"
        if "boosted" in label or "tuned" in label or label.startswith("CLASS cl_ref"):
            return (0, (4, 2))
        if label.startswith("HMCode feedback"):
            return "-."
        return "-"

    fig, ax = plt.subplots(figsize=(7.15, 4.15))
    ymin = 1.0e-4
    ax.fill_between(
        ell,
        ymin,
        sigma_a,
        where=np.isfinite(sigma_a) & (sigma_a > ymin),
        color="0.86",
        alpha=0.85,
        linewidth=0,
        zorder=0,
    )
    ax.axhline(1.0e-3, color="0.35", linewidth=0.9, linestyle=":")
    for label, values in curves.items():
        ax.plot(
            ell,
            values,
            color=color_for(label),
            linestyle=linestyle_for(label),
            linewidth=1.25,
            label=label,
        )
    ax.set_xscale("function", functions=(np.sqrt, np.square))
    ax.set_yscale("log")
    ax.set_xlim(2, ell[-1])
    ax.set_ylim(ymin, 1.0e-2)
    xticks = np.array([2, 10, 30, 100, 300, 1000, 2000, 3000, 4000, 6000])
    xticks = xticks[(xticks >= 2) & (xticks <= ell[-1])]
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(int(tick)) for tick in xticks])
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"lensed $C_\ell$ residual envelope")
    ax.grid(True, which="major", color="0.88", linewidth=0.6)
    ax.grid(True, which="minor", color="0.93", linewidth=0.35)
    ax.legend(loc="upper left", frameon=True, framealpha=0.96, borderpad=0.6)
    ax.text(
        0.018,
        0.745,
        r"$\max(|\Delta TT|/TT,\ |\Delta EE|/EE,\ |\Delta TE|/\sqrt{TT\,EE})$",
        transform=ax.transAxes,
        color="0.25",
    )
    fig.tight_layout(pad=0.6)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)


def main() -> None:
    args = build_parser().parse_args()
    if args.lmax <= 2:
        raise ValueError("--lmax must be greater than 2")
    if args.class_ref_lmax > args.lmax:
        raise ValueError("--class-ref-lmax cannot exceed --lmax")

    args.outdir.mkdir(parents=True, exist_ok=True)
    class_root = find_class_root(args.class_workdir)
    class_version = class_root.name.removeprefix("class_public-")
    precision = class_precision_profiles(class_root)
    class_ref_lensing_margin = (
        args.class_lensing_margin if args.class_ref_lensing_margin is None else args.class_ref_lensing_margin
    )
    class_permille_compute_lmax = args.lmax + args.class_permille_lensing_margin
    class_ref_compute_lmax = args.class_ref_lmax + class_ref_lensing_margin

    print(f"CAMB {camb.__version__}: {Path(camb.__file__).resolve()}")
    print(f"CLASS {class_version}: {class_root}")

    default = cached_case(
        args.outdir / f"camb_default_{RECOMBINATION_TAG}.npz",
        args.force,
        lambda: run_camb_case(
            make_camb_params(
                lmax=args.lmax,
                lens_potential_accuracy=None,
                lens_output_margin=200,
                boosted=False,
            ),
            args.lmax,
        ),
    )
    boosted = cached_case(
        args.outdir / f"camb_boosted_reference_{RECOMBINATION_TAG}.npz",
        args.force,
        lambda: run_camb_case(
            make_camb_params(
                lmax=args.lmax,
                lens_potential_accuracy=args.reference_lens_potential_accuracy,
                lens_output_margin=args.reference_lens_output_margin,
                boosted=True,
            ),
            args.lmax,
        ),
    )
    if args.class_k_max_tau0_over_l_max is not None:
        class_permille_kmax = args.class_k_max_tau0_over_l_max
        class_ref_kmax = args.class_k_max_tau0_over_l_max
    elif args.class_kmax_mode == "boosted-camb":
        class_permille_kmax = float(boosted["max_eta_k"]) / float(class_permille_compute_lmax)
        class_ref_kmax = float(boosted["max_eta_k"]) / float(class_ref_compute_lmax)
    else:
        class_permille_kmax = None
        class_ref_kmax = None
    class_permille_kmax_tag = kmax_tag(class_permille_kmax)
    class_ref_kmax_tag = kmax_tag(class_ref_kmax)
    feedback = cached_case(
        args.outdir / f"camb_hmcode_feedback_{RECOMBINATION_TAG}_logt{float_tag(args.feedback_logt_agn)}.npz",
        args.force,
        lambda: run_camb_case(
            make_camb_params(
                lmax=args.lmax,
                lens_potential_accuracy=args.reference_lens_potential_accuracy,
                lens_output_margin=args.reference_lens_output_margin,
                boosted=True,
                feedback_logt_agn=args.feedback_logt_agn,
            ),
            args.lmax,
        ),
    )
    class_permille = cached_case(
        args.outdir
        / f"class_cl_permille_{RECOMBINATION_TAG}_lmax{class_permille_compute_lmax}_{class_permille_kmax_tag}.npz",
        args.force,
        lambda: run_class_case(
            lmax=class_permille_compute_lmax,
            precision_overrides=precision["cl_permille"],
            k_max_tau0_over_l_max=class_permille_kmax,
        ),
    )
    class_ref = cached_case(
        args.outdir
        / (
            f"class_{args.class_ref_profile}_{RECOMBINATION_TAG}_compute_lmax{class_ref_compute_lmax}"
            f"_plot_lmax{args.class_ref_lmax}_{class_ref_kmax_tag}.npz"
        ),
        args.force,
        lambda: run_class_case(
            lmax=class_ref_compute_lmax,
            precision_overrides=precision[args.class_ref_profile],
            k_max_tau0_over_l_max=class_ref_kmax,
        ),
    )

    ell = np.arange(args.lmax + 1)
    curves = {
        "CAMB default / boosted": residual_envelope(default["dls"], boosted["dls"], args.lmax),
        "CLASS cl_permille / boosted CAMB": residual_envelope(class_permille["dls"], boosted["dls"], args.lmax),
        (
            f"{class_profile_label(args.class_ref_profile)} / boosted CAMB"
            if args.class_ref_lmax == args.lmax
            else f"{class_profile_label(args.class_ref_profile)} (L<={args.class_ref_lmax}) / boosted CAMB"
        ): mask_above(residual_envelope(class_ref["dls"], boosted["dls"], args.lmax), args.class_ref_lmax),
        "HMCode feedback / DM-only": residual_envelope(feedback["dls"], boosted["dls"], args.lmax),
    }
    ranges = {
        "low": (2, 599),
        "interior": (600, 3499),
        "tail": (3500, args.lmax),
        "all": (2, args.lmax),
    }
    summaries = {label: summarize_curve(values, ranges) for label, values in curves.items()}
    so_defaults = DEFAULT_NOISE_CONFIGS["so"]
    noise_config = NoiseConfig(
        name="so",
        noise_muK_arcmin_T=so_defaults["noise_muK_arcmin_T"],
        noise_muK_arcmin_P=so_defaults["noise_muK_arcmin_P"],
        fwhm_arcmin=so_defaults["fwhm_arcmin"],
        fsky=so_defaults["fsky"],
        lmin=2,
        lmax=args.lmax,
        fields=("T", "E"),
    )
    amplitude_error_bin_width = 0.5
    sigma_a = binned_amplitude_sigma(boosted["dls"], noise_config, args.lmax, amplitude_error_bin_width)

    plot_path = args.outdir / "lensed_cls_validation_residuals.pdf"
    make_plot(plot_path, ell, curves, sigma_a, summaries)

    metadata = {
        "plot": str(plot_path),
        "png": str(plot_path.with_suffix(".png")),
        "lmax": args.lmax,
        "class_ref_lmax": args.class_ref_lmax,
        "camb_version": camb.__version__,
        "camb_path": str(Path(camb.__file__).resolve()),
        "class_version": class_version,
        "class_root": str(class_root),
        "class_ref_profile": args.class_ref_profile,
        "class_ref_profile_label": class_profile_label(args.class_ref_profile),
        "class_precision_overrides": precision[args.class_ref_profile],
        "class_lensing_margin": args.class_lensing_margin,
        "class_permille_lensing_margin": args.class_permille_lensing_margin,
        "class_ref_lensing_margin": class_ref_lensing_margin,
        "class_permille_compute_lmax": class_permille_compute_lmax,
        "class_ref_compute_lmax": class_ref_compute_lmax,
        "recombination": {
            "tag": RECOMBINATION_TAG,
            "camb_model": "Recfast",
            "camb_recfast_approx_model": CAMB_RECFAST_APPROX_MODEL,
            "class_params": CLASS_PLANCK_RECFAST_PARAMS,
        },
        "cosmology": COSMOLOGY,
        "baseline_nonlinear": "CAMB/CLASS HMCode 2020 dark-matter-only",
        "feedback_curve": {
            "halofit_version": "mead2020_feedback",
            "HMCode_logT_AGN": args.feedback_logt_agn,
        },
        "boosted_reference": {
            "AccuracyBoost": 2.0,
            "lSampleBoost": 2.0,
            "lAccuracyBoost": 2.0,
            "IntTolBoost": 2.0,
            "lens_potential_accuracy": args.reference_lens_potential_accuracy,
            "lens_output_margin": args.reference_lens_output_margin,
            "max_eta_k": boosted["max_eta_k"],
            "max_l": boosted["max_l"],
        },
        "class_permille_k_max_tau0_over_l_max": class_permille_kmax,
        "class_ref_k_max_tau0_over_l_max": class_ref_kmax,
        "class_matched_max_eta_k": boosted["max_eta_k"] if args.class_kmax_mode == "boosted-camb" else None,
        "class_kmax_mode": args.class_kmax_mode,
        "standard_reference_noise": noise_config.__dict__,
        "amplitude_error_band": {
            "kind": "single-bin spectral-amplitude 1 sigma",
            "fractional_bin_width": amplitude_error_bin_width,
            "bin_definition": "for each plotted L, use integer multipoles in [ceil(0.75 L), floor(1.25 L)]",
        },
        "summaries": summaries,
        "run_times": {
            "camb_default": {"cpu_s": default["cpu_s"], "wall_s": default["wall_s"]},
            "camb_boosted_reference": {"cpu_s": boosted["cpu_s"], "wall_s": boosted["wall_s"]},
            "camb_hmcode_feedback": {"cpu_s": feedback["cpu_s"], "wall_s": feedback["wall_s"]},
            "class_cl_permille": {"cpu_s": class_permille["cpu_s"], "wall_s": class_permille["wall_s"]},
            f"class_{args.class_ref_profile}_compute_lmax{class_ref_compute_lmax}_{class_ref_kmax_tag}": {
                "cpu_s": class_ref["cpu_s"],
                "wall_s": class_ref["wall_s"],
            },
        },
        "note": (
            "CLASS lensed spectra are computed above the plotted maximum and then truncated, to reduce high-ell "
            "lensing-convolution edge effects. The CLASS reference curve is masked above class_ref_lmax."
        ),
    }
    (args.outdir / "lensed_cls_validation_residuals.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {plot_path}")
    print(f"Wrote {plot_path.with_suffix('.png')}")
    print(f"Wrote {args.outdir / 'lensed_cls_validation_residuals.json'}")


if __name__ == "__main__":
    main()
