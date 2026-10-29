"""Regenerate the lag_corrections accuracy table and the checks on its ground truth.

The table gives RECFAST's peak x_e error, ln theta* error, lensed TT error and run time with the
lag_corrections option off and on, each against a ground truth that the script builds from the
option at its tightest settings. Self-contained: needs only `camb` (the build under test) and
`numpy`. The changelog entries docs/changelog/2026-10-23-recfast-lag-corrections.md,
2026-10-23-recfast-heiii-saha-continued.md and 2026-10-23-recfast-snap-saha-z.md quote its output.

    python scripts/recfast_lag_corrections.py            # the four-row table
    python scripts/recfast_lag_corrections.py --all      # and the supporting comparisons
    python scripts/recfast_lag_corrections.py --latex    # and the LaTeX tabular body
    python scripts/recfast_lag_corrections.py --ede      # and the early dark energy scan, about 70 s more

The supporting comparisons are the gains at CAMB's own settings and at equal cost, the
solver-tolerance and grid scans, CAMB's best TT over its own settings with the HeIII step,
self-consistency, the T_m offset, the implicit route, the redshift at which the delta x_e
perturbation starts, and the reference that docs/changelog/recfast_rosenbrock.md tunes against.
The early dark energy scan measures both paths in 16 AxionEffectiveFluid backgrounds, each against its
own ground truth.

Everything is at CAMB's default accuracy -- lmax 2500, lens_potential_accuracy 4 (which
set_for_lmax picks for this lmax), TimeStepBoost 1, Nz 2046 -- and the settings are printed with
the numbers.

RECFAST settings it sets: use_rosenbrock, rosenbrock_handoff_xH, rosenbrock_tol, Nz,
heIII_saha_continued, lag_corrections, Tm_handover_frac, hydrogen_saha_xH, Saha_He_end_xHe and
rk45_tol. It leaves snap_saha_z alone, since lag_corrections turns it on.

Every candidate is measured against the reference on the *same* Nz; see `Reference`.

Where a redshift is needed it is derived from x_e; see `switch_redshifts`.
"""

from __future__ import annotations

import argparse
import time

import numpy as np

import camb
from camb import constants
from camb.dark_energy import AxionEffectiveFluid

# ---------------------------------------------------------------- configuration

COSMOLOGY = {"H0": 67.4, "ombh2": 0.0224, "omch2": 0.12, "mnu": 0.06}
PRIMORDIAL = {"tau": 0.0544, "As": 2.1e-9, "ns": 0.9649}
LMAX = 2500
L_MIN = 30  # the TT comparison window is L_MIN <= l <= LMAX
NZ_DEFAULT = 2046  # RECFAST's own grid; z nodes are 1e4 - i*(1e4/Nz)
TT_RECHECK = 1.0e-6  # a table row's TT error below this is recomputed at TimeStepBoost 4 and the change printed
X_WINDOW = (400.0, 3400.0)  # the x_e comparison window
X_WINDOW_HEIII = (400.0, 8000.0)  # the same, including the HeIII step at z = 5000

# The option. Two routes: the explicit solver (Rosenbrock off) and the implicit one.
OPTION_EXPLICIT = {"lag_corrections": True, "use_rosenbrock": False}
OPTION_IMPLICIT = {"lag_corrections": True}
LEGACY: dict = {}  # the new settings off, CAMB's default path

# Ground truth: the option at its tightest setting that runs, on the explicit route. Three to
# four decades tighter in every threshold than the tightest table row. Three things make it
# defensible, and `cross_check` prints all three:
#   (1) it converges -- one step looser in each threshold moves it by 6.1e-9 in x_e, 4.8e-10 in TT;
#   (2) the two solvers agree -- Rosenbrock on and off differ by 1.3e-8 in x_e, 1.3e-9 in TT;
#   (3) the start is forgotten -- at Tm_handover_frac <= 1e-7 the handover has moved back to the
#       first integrated row, so the full stiff T_m equation runs throughout and the history is
#       bitwise independent of the handover (identical at 1e-8, 1e-10, 1e-12 and 0).
# hydrogen_saha_xH cannot go past 0.99999: at 0.999999 the explicit solver fails at z ~ 2350,
# where x_H crosses the Kholupenko switch-on at 1 - 1e-7 (the x_H < 0.9999999d0 tests in recfast.f90).
REFERENCE = dict(
    OPTION_EXPLICIT,
    Saha_He_end_xHe=0.999999,
    hydrogen_saha_xH=0.99999,
    Tm_handover_frac=1.0e-8,
    rk45_tol=1.0e-10,
)
REFERENCE_LOOSER = dict(  # one step looser in each threshold, for the self-convergence bound
    OPTION_EXPLICIT,
    Saha_He_end_xHe=0.99999,
    hydrogen_saha_xH=0.9999,
    Tm_handover_frac=1.0e-7,
    rk45_tol=1.0e-10,
)
REFERENCE_IMPLICIT = dict(  # the same point on the implicit route; 1.2e-8 is its floor
    OPTION_IMPLICIT,
    Saha_He_end_xHe=0.999999,
    hydrogen_saha_xH=0.99999,
    Tm_handover_frac=1.0e-8,
    rk45_tol=1.0e-10,
    rosenbrock_tol=1.2e-8,
)
GRID_REFERENCE_NZ = 16 * NZ_DEFAULT  # the grid-converged reference, for the grid comparison only
# The reference docs/changelog/recfast_rosenbrock.md tunes against, for section 10
CHANGELOG_REFERENCE = {"use_rosenbrock": True, "rosenbrock_handoff_xH": 0.0, "rosenbrock_tol": 1e-7, "Nz": 160000}

ZNODES = 1.0e4 - np.arange(1, NZ_DEFAULT + 1) * (1.0e4 / NZ_DEFAULT)

# ---------------------------------------------------------------- machinery


def params(recomb: dict, *, cls: bool = True, boost: float = 1.0, fluid: tuple | None = None):
    """CAMB parameters for one RECFAST configuration; fluid = (w_n, fde_zc, zc) adds an AxionEffectiveFluid."""
    p = camb.set_params(WantCls=cls, lmax=LMAX if cls else 0, **COSMOLOGY, **PRIMORDIAL)
    if fluid is not None:
        de = AxionEffectiveFluid()
        de.set_params(w_n=fluid[0], fde_zc=fluid[1], zc=fluid[2])
        p.DarkEnergy = de
    if boost != 1.0:
        p.Accuracy.TimeStepBoost = boost
    for key, value in recomb.items():
        if not hasattr(p.Recomb, key):
            raise SystemExit(f"this CAMB build has no Recfast setting {key!r}")
        setattr(p.Recomb, key, value)
    return p


def solve(recomb: dict, *, boost: float = 1.0, nodes: np.ndarray | None = None, fluid: tuple | None = None) -> dict:
    """x_e and T_m at the grid nodes, the lensed spectra, and theta*."""
    p = params(recomb, boost=boost, fluid=fluid)
    results = camb.get_results(p)
    z = ZNODES if nodes is None else nodes
    hist = results.get_background_redshift_evolution(z, ["x_e", "T_b"], format="array")
    derived = results.get_derived_params()
    return {
        "z": z,
        "x_e": hist[:, 0],
        "T_m": hist[:, 1],
        "tt": results.get_cmb_power_spectra(p, CMB_unit="muK")["total"][:, 0],
        "thetastar": derived["thetastar"],
        "zstar": derived["zstar"],
    }


def peak_absolute(a: np.ndarray, b: np.ndarray, z: np.ndarray, window) -> tuple[float, float]:
    inside = (z > window[0]) & (z < window[1])
    diff = np.abs(a[inside] - b[inside])
    i = int(np.argmax(diff))
    return float(diff[i]), float(z[inside][i])


def peak_fractional(a: np.ndarray, b: np.ndarray, z: np.ndarray, window) -> float:
    inside = (z > window[0]) & (z < window[1])
    return float(np.max(np.abs(a[inside] / b[inside] - 1.0)))


def tt_error(tt: np.ndarray, reference: np.ndarray) -> tuple[float, int, float]:
    ell = np.arange(reference.size)
    inside = (ell >= L_MIN) & (ell <= LMAX)
    frac = np.zeros_like(reference)
    nonzero = reference != 0
    frac[nonzero] = (tt[nonzero] - reference[nonzero]) / reference[nonzero]
    i = int(np.abs(frac[inside]).argmax())
    return (float(np.abs(frac[inside]).max()), int(ell[inside][i]), float(np.sqrt((frac[inside] ** 2).mean())))


def timings(configs: list[dict], repeats: int = 30) -> list[dict]:
    """get_background, interleaved so machine drift hits every entry alike. configs[0] is the
    baseline. Returns best, median and the best-relative-to-baseline for each."""
    built = [params(c, cls=False) for c in configs]
    samples: list[list[float]] = [[] for _ in built]
    for _ in range(repeats):
        for i, p in enumerate(built):
            t0 = time.perf_counter()
            camb.get_background(p)
            samples[i].append(time.perf_counter() - t0)
    best = [min(s) for s in samples]
    return [
        {"best_ms": b * 1e3, "median_ms": float(np.median(s)) * 1e3, "relative": b / best[0]}
        for b, s in zip(best, samples)
    ]


class Reference:
    """The ground truth. Solved once per (Nz, TimeStepBoost) actually asked for.

    Matched-grid rule: a candidate is always measured against the reference on the *same* Nz.
    Comparing a refined candidate with a reference on the default grid measures the difference
    between the two grids, which is far larger than either one's error -- at Nz 2046 against
    Nz x16 it is 7e-7 in x_e, a ringing packet over 2339 < z < 2375 where x_H crosses the
    Kholupenko switch-on (the x_H < 0.9999999d0 tests in recfast.f90), present identically in
    every solver. Only the deliberate grid comparison overrides this, via ref_nz, and it says so
    when it does.
    """

    def __init__(self) -> None:
        self.cache: dict[tuple[int, float], dict] = {}
        self.r = self.at()

    def at(self, nz: int = NZ_DEFAULT, boost: float = 1.0) -> dict:
        key = (int(nz), float(boost))
        if key not in self.cache:
            cfg = REFERENCE if nz == NZ_DEFAULT else dict(REFERENCE, Nz=int(nz))
            self.cache[key] = solve(cfg, boost=boost)
        return self.cache[key]

    def measure(self, recomb: dict, *, ref_nz: int | None = None) -> dict:
        nz = int(recomb.get("Nz", NZ_DEFAULT)) if ref_nz is None else int(ref_nz)
        self.r = base = self.at(nz)
        s = solve(recomb)
        dxe, zxe = peak_absolute(s["x_e"], base["x_e"], ZNODES, X_WINDOW)
        dxe_heiii, zheiii = peak_absolute(s["x_e"], base["x_e"], ZNODES, X_WINDOW_HEIII)
        tt, at_l, rms = tt_error(s["tt"], base["tt"])
        row = {
            "d_xe": dxe,
            "d_xe_at_z": zxe,
            "d_xe_heIII": dxe_heiii,
            "d_xe_heIII_at_z": zheiii,
            "d_Tm": peak_fractional(s["T_m"], base["T_m"], ZNODES, X_WINDOW),
            "tt": tt,
            "tt_at_l": at_l,
            "tt_rms": rms,
            "dln_thetastar": s["thetastar"] / base["thetastar"] - 1.0,
            "ref_nz": nz,
            "boost4": None,
        }
        if tt < TT_RECHECK:  # at this size CAMB's own time-step sampling can set the TT error
            tt4, at_l4, _ = tt_error(solve(recomb, boost=4.0)["tt"], self.at(nz, 4.0)["tt"])
            moved = abs(tt4 / tt - 1.0) if tt else 0.0
            row["boost4"] = {
                "tt_boost4": tt4,
                "at_l": at_l4,
                "moved": moved,
                "suspect": bool(moved > 0.05 or at_l < 100 or at_l4 < 100),
            }
        return row


def switch_redshifts(recomb: dict) -> float:
    """The delta-x_e perturbation switch, derived from x_e because Recombination_saha_z is not
    reachable from Python: get_Saha_z is a Fortran type-bound procedure that recombination.py does
    not expose, and recombination_saha_tau lives in TThermoData, which CAMBdata does not mirror.

    With snap_saha_z the recorder takes the last grid node with x_H > 0.985, whichever branch
    produced it; helium is neutral there, so the last node with x_e > 0.985 locates it to better
    than one node. The value recorded with snap_saha_z off cannot be derived this way, as the note
    printed in section 9 explains.
    """
    s = solve(dict(recomb), nodes=ZNODES)
    window = (ZNODES > 1200.0) & (ZNODES < 3000.0)  # excludes reionization, where x_e rises again
    above_985 = ZNODES[window & (s["x_e"] > 0.985)]
    return float(above_985.min()) if above_985.size else float("nan")


def heIII_saha_step() -> tuple[float, float, float]:
    """The size of the HeIII discontinuity itself: the He++ <-> He+ Saha root at z = 5000 minus
    1 + fHe, which is what the row below the node is set to. This is the number the continuation
    removes, and unlike a reading of x_e it does not depend on how the output is sampled. The
    formula is the He++ <-> He+ Saha root in TRecfast_init, verbatim, with the run's own Y_p, T_CMB
    and n_H and
    CAMB's own constants (camb/constants.py, the values fortran/constants.f90 uses)."""
    p = params(LEGACY, cls=False)
    background = camb.get_background(p)
    fHe, Tnow = background.fHe, float(p.TCMB)
    h_P, c_light, k_B, m_e = constants.h_p, constants.c, constants.k_B, constants.m_e
    G_newton, m_H, L_He2_ion = constants.G, constants.m_H, 4.389088863e7  # L_He2_ion as in recfast.f90
    CB1_He2 = h_P * c_light * L_He2_ion / k_B
    CR = 2.0 * np.pi * (m_e / h_P) * (k_B / h_P)
    bigH = 100.0e3 / constants.Mpc
    mu_H = 1.0 / (1.0 - float(p.YHe))
    n_H = 3.0 * bigH**2 * p.ombh2 / (8.0 * np.pi * G_newton * mu_H * m_H)
    ainv = 1.0 + 5000.0
    rhs = np.exp(1.5 * np.log(CR * Tnow / ainv) - CB1_He2 / (Tnow * ainv)) / n_H
    x0 = 0.5 * (np.sqrt((rhs - 1.0 - fHe) ** 2 + 4.0 * (1.0 + 2.0 * fHe) * rhs) - (rhs - 1.0 - fHe))
    return float(x0 - 1.0 - fHe), float(x0), float(1.0 + fHe)


# ---------------------------------------------------------------- the table

TABLE = [
    ("CAMB's own settings, option off", LEGACY),
    (
        "option on, Rosenbrock off",
        dict(OPTION_EXPLICIT, Saha_He_end_xHe=0.99, hydrogen_saha_xH=0.99, Tm_handover_frac=3.7e-4),
    ),
    (
        "option on, equal cost",
        dict(OPTION_EXPLICIT, Saha_He_end_xHe=0.995, hydrogen_saha_xH=0.995, Tm_handover_frac=1.25e-4),
    ),
    (
        "option on, tighter still",
        dict(OPTION_EXPLICIT, Saha_He_end_xHe=0.999, hydrogen_saha_xH=0.997, Tm_handover_frac=4.5e-5, rk45_tol=1e-8),
    ),
]


def banner() -> None:
    p = params(LEGACY, cls=False)
    print("RECFAST lag_corrections: accuracy table and ground-truth checks")
    print(f"  camb {camb.__version__} from {camb.__file__}")
    print(
        f"  cosmology  H0 {COSMOLOGY['H0']}, ombh2 {COSMOLOGY['ombh2']}, "
        f"omch2 {COSMOLOGY['omch2']}, mnu {COSMOLOGY['mnu']}, "
        f"tau {PRIMORDIAL['tau']}, As {PRIMORDIAL['As']:g}, ns {PRIMORDIAL['ns']}"
    )
    print(f"  derived    YHe {p.YHe:.10f}, N_eff {p.N_eff:.4f}, TCMB {p.TCMB} K")
    print(
        f"  accuracy   lmax {LMAX}, lens_potential_accuracy 4 (set_for_lmax default at this "
        f"lmax), TimeStepBoost 1, Nz {NZ_DEFAULT}"
    )
    print(
        f"  errors     against the option at {REFERENCE['Saha_He_end_xHe']} / "
        f"{REFERENCE['hydrogen_saha_xH']} / {REFERENCE['Tm_handover_frac']:g} / "
        f"rk45_tol {REFERENCE['rk45_tol']:g}, Rosenbrock off,"
    )
    print("             solved on the same Nz as whatever it is compared with")
    print(
        f"  windows    x_e over {X_WINDOW[0]:.0f} < z < {X_WINDOW[1]:.0f} "
        f"(and < {X_WINDOW_HEIII[1]:.0f} where stated), TT over {L_MIN} <= l <= {LMAX}"
    )
    print()


def _against(label: str, a: dict, b: dict) -> dict:
    dxe, at_z = peak_absolute(a["x_e"], b["x_e"], ZNODES, X_WINDOW)
    tt, at_l, _ = tt_error(a["tt"], b["tt"])
    dln = a["thetastar"] / b["thetastar"] - 1.0
    print(
        f"  {label:38s} x_e {dxe:9.3e} at z {at_z:6.0f}   "
        f"T_m {peak_fractional(a['T_m'], b['T_m'], ZNODES, X_WINDOW):9.3e}   "
        f"ln theta* {dln:+10.3e}   "
        f"TT {tt:9.3e} at l {at_l}"
    )
    return {"d_xe": dxe, "tt": tt, "dln_thetastar": dln}


def cross_check(ref: Reference) -> None:
    """Why the reference can be trusted: (1) its own convergence, (2) the agreement of the two
    solvers and (3) its independence of the T_m start. A check that cannot run is skipped, with a note."""
    print(f"the reference and its error budget   (all at Nz {NZ_DEFAULT})")
    checks = [_against("(1) one step looser in each threshold", ref.r, solve(REFERENCE_LOOSER))]
    try:
        implicit = solve(REFERENCE_IMPLICIT)
    except SystemExit as exc:
        print(f"  (2) Rosenbrock on against off           skipped: {exc}")
    else:
        checks.append(_against("(2) Rosenbrock on against off", implicit, ref.r))
    checks.append(_against("(3) Tm_handover_frac -> 0", solve(dict(REFERENCE, Tm_handover_frac=0.0)), ref.r))
    print(
        f"  {'bound, worst of (1)-(3)':38s} "
        f"x_e {max(b['d_xe'] for b in checks):9.3e}             "
        f"ln theta* {max(abs(b['dln_thetastar']) for b in checks):10.3e}   "
        f"TT {max(b['tt'] for b in checks):9.3e}"
    )
    print()


def print_table(rows: list[tuple[str, dict]], times: list[dict]) -> None:
    print(f"{'row':34s} {'peak |dx_e|':>12s} {'at z':>6s} {'ln theta*':>11s} {'TT max':>10s} {'at l':>5s} {'time':>6s}")
    for (label, _), row, t in zip(TABLE, rows, times):
        print(
            f"{label:34s} {row['d_xe']:12.3e} {row['d_xe_at_z']:6.0f} "
            f"{row['dln_thetastar']:+11.3e} {row['tt']:10.3e} {row['tt_at_l']:5d} "
            f"{t['relative']:6.3f}"
        )
    print()
    for (label, _), row in zip(TABLE, rows):
        g = row["boost4"]
        if g:
            mark = "  <-- CHECK" if g["suspect"] else ""
            print(
                f"  boost 4  {label:30s} TT {row['tt']:.3e} at l {row['tt_at_l']} -> "
                f"{g['tt_boost4']:.3e} at l {g['at_l']} at TimeStepBoost 4, "
                f"moved {g['moved'] * 100:.1f}%{mark}"
            )
    print()
    for (label, _), t in zip(TABLE, times):
        print(
            f"  timing {label:32s} best {t['best_ms']:.3f} ms, median {t['median_ms']:.3f} ms, "
            f"relative {t['relative']:.3f}"
        )
    print()


def sci(value: float, digits: int = 1) -> str:
    """LaTeX scientific notation, e.g. 5.8\\times 10^{-5}."""
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = abs(value) / 10.0**exponent
    return f"{sign}{mantissa:.{digits}f}\\times 10^{{{exponent}}}"


def print_latex(rows: list[dict], times: list[dict]) -> None:
    print("% tabular body: setting & peak |dx_e| & ln theta* & TT & time")
    for (label, _), row, t in zip(TABLE, rows, times):
        print(
            f"{label} & ${sci(row['d_xe'])}$ & ${sci(row['dln_thetastar'])}$ & "
            f"${sci(row['tt'])}$ & ${t['relative']:.2f}$ \\\\"
        )
    print()


# ---------------------------------------------------------------- the supporting comparisons


def text_numbers(ref: Reference, rows: list[dict], times: list[dict]) -> None:
    legacy, option, equal, tight = rows
    print("=" * 78)
    print("supporting comparisons")
    print("=" * 78)

    print("\n1. the option on against off, at CAMB's own settings")
    print(f"   x_e          {legacy['d_xe']:.3e} -> {option['d_xe']:.3e}   {legacy['d_xe'] / option['d_xe']:.1f}x")
    step, saha, one_plus_fHe = heIII_saha_step()
    print(
        f"   x_e + HeIII  {step:.3e} -> {option['d_xe']:.3e}   {step / option['d_xe']:.1f}x"
        f"   (the Saha discontinuity at z = 5000: {saha:.9f} - {one_plus_fHe:.9f})"
    )
    print(
        f"     the same from a reading of x_e, which depends on the output sampling: "
        f"{legacy['d_xe_heIII']:.3e} at z = {legacy['d_xe_heIII_at_z']:.0f} -> "
        f"{option['d_xe_heIII']:.3e}   {legacy['d_xe_heIII'] / option['d_xe_heIII']:.1f}x"
    )
    print(f"   TT           {legacy['tt']:.3e} -> {option['tt']:.3e}   {legacy['tt'] / option['tt']:.1f}x")
    print(
        f"   ln theta*    {legacy['dln_thetastar']:+.3e} -> {option['dln_thetastar']:+.3e}   "
        f"{abs(legacy['dln_thetastar'] / option['dln_thetastar']):.1f}x"
    )
    print(
        f"   time         {times[1]['relative']:.3f}  "
        f"({(1 - times[1]['relative']) * 100:.1f}% faster; medians "
        f"{times[0]['median_ms']:.3f} and {times[1]['median_ms']:.3f} ms)"
    )

    print("\n2. at equal cost")
    print(
        f"   x_e  {legacy['d_xe'] / equal['d_xe']:.1f}x      TT  "
        f"{legacy['tt'] / equal['tt']:.1f}x      time {times[2]['relative']:.3f}"
    )

    print("\n3. the same solver tolerance on both paths (Rosenbrock on, handoff 0.976)")
    for name, base in (
        ("option off", LEGACY),
        (
            "option on",
            dict(OPTION_IMPLICIT, Saha_He_end_xHe=0.99, rosenbrock_handoff_xH=0.976, Tm_handover_frac=3.7e-4),
        ),
    ):
        loose = ref.measure(dict(base, rosenbrock_tol=3e-4))
        tightened = ref.measure(dict(base, rosenbrock_tol=3e-6))
        ratio = tightened["tt"] / loose["tt"]
        verb = "worse" if ratio > 1 else "better"
        print(f"   {name:14s} TT {loose['tt']:.3e} -> {tightened['tt']:.3e}   {max(ratio, 1 / ratio):.1f}x {verb}")

    print("\n4. an 8x finer RECFAST grid")
    print(
        f"   (the one place the matched-grid rule is overridden: both sides are measured "
        f"against\n    the grid-converged reference at Nz {GRID_REFERENCE_NZ}, because the "
        f"question is what\n    refinement does to each path, not what it does to the grid)"
    )
    for name, base in (("option off", LEGACY), ("option on", TABLE[1][1])):
        coarse = ref.measure(dict(base), ref_nz=GRID_REFERENCE_NZ)
        fine = ref.measure(dict(base, Nz=8 * NZ_DEFAULT), ref_nz=GRID_REFERENCE_NZ)
        ratio = fine["tt"] / coarse["tt"]
        print(
            f"   {name:14s} TT {coarse['tt']:.3e} -> {fine['tt']:.3e} "
            f"({max(ratio, 1 / ratio):.2f}x {'worse' if ratio > 1 else 'better'}), "
            f"x_e {coarse['d_xe']:.3e} -> {fine['d_xe']:.3e}, "
            f"ln theta* {coarse['dln_thetastar']:+.2e} -> {fine['dln_thetastar']:+.2e}"
        )

    print("\n5. the option-off path's best TT over its own knobs, and the HeIII step")
    best = (float("inf"), "")
    for name, cfg in (
        ("default", LEGACY),
        ("helium threshold 0.999", {"Saha_He_end_xHe": 0.999}),
        ("rosenbrock_tol 3e-6", {"rosenbrock_tol": 3e-6}),
        ("handoff 0.99", {"rosenbrock_handoff_xH": 0.99}),
        ("Nz 8x", {"Nz": 8 * NZ_DEFAULT}),
    ):
        m = ref.measure(dict(LEGACY, **cfg))
        best = min(best, (m["tt"], name))
        print(
            f"   {name:24s} TT {m['tt']:.3e}   x_e {m['d_xe']:.3e}   "
            f"HeIII-window x_e {m['d_xe_heIII']:.3e} at z = {m['d_xe_heIII_at_z']:.0f}"
        )
    print(f"   best: {best[0]:.3e} at {best[1]}")
    with_heiii = ref.measure(dict(LEGACY, heIII_saha_continued=True))
    step, _, _ = heIII_saha_step()
    print(f"   the step itself, from the Saha root at z = 5000: {step:.4e}")
    print(
        f"   read off x_e instead, the continuation removes it: {legacy['d_xe_heIII']:.3e} -> "
        f"{with_heiii['d_xe_heIII']:.3e}   (output-sampling dependent; see heIII_saha_step)"
    )

    print("\n6. self-consistency")
    print(
        f"   the tightest table row against the reference: x_e {tight['d_xe']:.3e}, "
        f"T_m {tight['d_Tm']:.3e}, ln theta* {tight['dln_thetastar']:+.3e}"
    )
    # rosenbrock_tol 1.2e-8 is the implicit solver's floor; below it the integration fails. Left
    # at its 3e-4 default the implicit route is nowhere near converged and the comparison measures
    # that instead.
    on = dict(
        OPTION_IMPLICIT,
        Saha_He_end_xHe=0.9999,
        rosenbrock_handoff_xH=0.999,
        Tm_handover_frac=1e-5,
        rk45_tol=1e-9,
        rosenbrock_tol=1.2e-8,
    )
    off = dict(OPTION_EXPLICIT, Saha_He_end_xHe=0.9999, hydrogen_saha_xH=0.999, Tm_handover_frac=1e-5, rk45_tol=1e-9)
    a, b = solve(on), solve(off)
    dxe, at_z = peak_absolute(a["x_e"], b["x_e"], ZNODES, X_WINDOW)
    print(
        f"   Rosenbrock on against off at that setting: x_e {dxe:.3e} at z = {at_z:.0f}, "
        f"ln theta* {a['thetastar'] / b['thetastar'] - 1:+.3e}"
    )

    print("\n7. the matter-temperature offset and CAMB's largest x_e error")
    sl = solve(LEGACY)
    i1100 = int(np.argmin(np.abs(ZNODES - 1100.0)))
    print(
        f"   |T_m - reference| at z = {ZNODES[i1100]:.1f}: "
        f"{abs(sl['T_m'][i1100] - ref.r['T_m'][i1100]):.3e} K "
        f"(T_m there is {sl['T_m'][i1100]:.4f} K)"
    )
    print(f"   largest x_e error with the option off: {legacy['d_xe']:.3e} at z = {legacy['d_xe_at_z']:.0f}")

    print("\n8. the implicit route")
    for handoff in (0.976, 0.995):
        cfg = dict(OPTION_IMPLICIT, Saha_He_end_xHe=0.99, rosenbrock_handoff_xH=handoff, Tm_handover_frac=3.7e-4)
        m = ref.measure(cfg)
        t = timings([LEGACY, cfg], repeats=10)[1]
        print(f"   handoff {handoff}: TT {m['tt']:.3e}, x_e {m['d_xe']:.3e}, time {t['relative']:.3f}")

    print("\n9. the delta-x_e perturbation switch, derived from x_e (see switch_redshifts)")
    print(
        f"   node with snap_saha_z, handoff 0.996: z = "
        f"{switch_redshifts(dict(LEGACY, rosenbrock_handoff_xH=0.996)):.2f}"
    )
    print(f"   node with snap_saha_z, default handoff 0.976: z = {switch_redshifts(LEGACY):.2f}")
    print("   the value recorded with snap_saha_z off at handoff 0.996 is NOT printed here, and")
    print("   cannot be: it is written inside ION at an integrator stage point, which")
    print("   overwrites the node written on the Rosenbrock rows. It lands near z = 1555, part")
    print("   way between nodes, so no function of the node values can reproduce it. At the")
    print("   default handoff the ION write never fires and the two recorders agree on the node.")

    print("\n10. the reference that docs/changelog/recfast_rosenbrock.md tunes against")
    print("   (full-range Rosenbrock, handoff 0, rosenbrock_tol 1e-7; it and the ground truth both on Nz 160000)")
    truth = solve(dict(REFERENCE, Nz=CHANGELOG_REFERENCE["Nz"]))["x_e"]
    i2815 = int(np.argmin(np.abs(ZNODES - 2815.2)))
    inside = (ZNODES > X_WINDOW[0]) & (ZNODES < X_WINDOW[1])
    for label, cfg in (
        ("that reference", CHANGELOG_REFERENCE),
        ("with its helium Saha end at 0.99999", dict(CHANGELOG_REFERENCE, Saha_He_end_xHe=0.99999)),
    ):
        d = solve(cfg)["x_e"] - truth
        j = int(np.argmax(np.abs(d[inside])))
        print(
            f"   {label:36s} x_e {d[i2815]:+.3e} at z = {ZNODES[i2815]:.1f}, "
            f"largest {d[inside][j]:+.3e} at z = {ZNODES[inside][j]:.1f} over 400 < z < 3400"
        )
    print()


# ---------------------------------------------------------------- the early dark energy scan

# AxionEffectiveFluid backgrounds (w_n, f_de at z_c, z_c) for --ede
EDE_FLUIDS = [
    (w_n, fde, zc) for w_n in (1.0 / 3.0, 1.0) for fde in (0.05, 0.1) for zc in (1000.0, 1700.0, 3500.0, 5000.0)
]


def ede_scan() -> None:
    """Both paths against the ground truth in LCDM and in the 16 EDE_FLUIDS backgrounds, each with its own ground
    truth. With the option off, RECFAST's tight-coupling temperature equation takes dH/dz from a fit in matter,
    radiation and curvature that leaves the fluid out, while the option needs no dH/dz. The scan runs at
    TimeStepBoost 1 and does not recompute its TT errors at TimeStepBoost 4."""
    print("=" * 78)
    print("early dark energy: AxionEffectiveFluid backgrounds, each against its own ground truth")
    print("=" * 78)
    print(f"   {'background':33s} {'path':11s} {'peak |dx_e|':>11s} {'ln theta*':>11s} {'TT max':>10s} {'at l':>5s}")
    result, truths = {}, {}
    for fluid in [None, *EDE_FLUIDS]:
        name = "LCDM" if fluid is None else f"w_n {fluid[0]:.3g}, f_de {fluid[1]:g} at z_c {fluid[2]:g}"
        truths[fluid] = truth = solve(REFERENCE, fluid=fluid)
        for path, cfg in (("option off", LEGACY), ("option on", OPTION_EXPLICIT)):
            s = solve(cfg, fluid=fluid)
            dxe, _ = peak_absolute(s["x_e"], truth["x_e"], ZNODES, X_WINDOW)
            tt, at_l, _ = tt_error(s["tt"], truth["tt"])
            dln = s["thetastar"] / truth["thetastar"] - 1
            result[(fluid, path)] = (dxe, dln, tt)
            print(f"   {name:33s} {path:11s} {dxe:11.3e} {dln:+11.3e} {tt:10.3e} {at_l:5d}")
    lcdm_off, lcdm_on = result[(None, "option off")], result[(None, "option on")]
    off = [result[(f, "option off")] for f in EDE_FLUIDS]
    on = [result[(f, "option on")] for f in EDE_FLUIDS]
    worst = max(EDE_FLUIDS, key=lambda f: result[(f, "option off")][2])
    above = sum(r[2] > lcdm_off[2] for r in off)
    print(
        f"\n   option off: TT {min(r[2] for r in off):.3e} to {max(r[2] for r in off):.3e}, above LCDM's "
        f"{lcdm_off[2]:.3e} in {above} of {len(off)}, at most {result[(worst, 'option off')][2] / lcdm_off[2]:.2f} "
        f"times LCDM (w_n {worst[0]:.3g}, f_de {worst[1]:g} at z_c {worst[2]:g})"
    )
    print(f"   option on:  TT {min(r[2] for r in on):.3e} to {max(r[2] for r in on):.3e} (LCDM {lcdm_on[2]:.3e})")
    print(
        f"   option off: |ln theta*| {min(abs(r[1]) for r in off):.3e} to {max(abs(r[1]) for r in off):.3e} "
        f"(LCDM {abs(lcdm_off[1]):.3e})"
    )
    implicit = solve(REFERENCE_IMPLICIT, fluid=worst)
    dxe, _ = peak_absolute(implicit["x_e"], truths[worst]["x_e"], ZNODES, X_WINDOW)
    tt, _, _ = tt_error(implicit["tt"], truths[worst]["tt"])
    print(f"   ground truth in that background, Rosenbrock on against off: x_e {dxe:.3e}, TT {tt:.3e}")
    print()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--all", action="store_true", help="also print the supporting comparisons")
    ap.add_argument("--latex", action="store_true", help="print the tabular body as well")
    ap.add_argument("--ede", action="store_true", help="also scan 16 early dark energy fluids (about 70 s more)")
    args = ap.parse_args()

    started = time.time()
    banner()
    ref = Reference()
    cross_check(ref)
    rows = [ref.measure(cfg) for _, cfg in TABLE]
    times = timings([cfg for _, cfg in TABLE])
    print_table(rows, times)
    if args.latex:
        print_latex(rows, times)
    if args.all:
        text_numbers(ref, rows, times)
    if args.ede:
        ede_scan()
    print(f"total wall clock: {time.time() - started:.1f} s")


if __name__ == "__main__":
    main()
