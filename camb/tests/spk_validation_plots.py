import argparse
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from _spk_helpers import compute_suppression_data

_RELATION_NAMES = {1: "power_law", 2: "cosmo_power_law", 3: "double_power_law"}
_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]


def _plot_validation(data, relation_kinds, redshifts, so, out_dir):
    """3 rows (param sets) x 3 cols (relations), each cell has suppression + residual."""
    n_params = 3
    fig, axes = plt.subplots(
        n_params * 2,
        3,
        figsize=(18, n_params * 5),
        squeeze=False,
        gridspec_kw={"height_ratios": [3, 1] * n_params, "hspace": 0.08},
    )

    for icol, relation_kind in enumerate(relation_kinds):
        for ip in range(n_params):
            ax_top = axes[ip * 2, icol]
            ax_bot = axes[ip * 2 + 1, icol]

            for ic, z in enumerate(redshifts):
                k, sup_camb, sup_ref, rel = data[(relation_kind, so, ip, z)]
                color = _COLORS[ic % len(_COLORS)]
                ax_top.plot(k, sup_camb, lw=1.6, color=color, label=f"z={z}" if icol == 0 and ip == 0 else None)
                ax_top.plot(
                    k, sup_ref, lw=1.0, ls="--", color="k", label="pySPK" if ic == 0 and icol == 0 and ip == 0 else None
                )
                ax_bot.plot(k, rel, lw=1.2, color=color)

            ax_top.axvline(8.0, color="grey", ls=":", lw=0.9, alpha=0.6)
            ax_top.set_xscale("log")
            ax_top.grid(True, alpha=0.25)
            ax_top.tick_params(labelbottom=False)

            if ip == 0:
                ax_top.set_title(f"{_RELATION_NAMES[relation_kind]} (kind={relation_kind})")

            ax_bot.axhline(0.0, color="k", lw=0.6, alpha=0.5)
            ax_bot.axvline(8.0, color="grey", ls=":", lw=0.9, alpha=0.6)
            ax_bot.set_xscale("log")
            ax_bot.grid(True, alpha=0.25)

            if ip == n_params - 1:
                ax_bot.set_xlabel("k [h/Mpc]")
            else:
                ax_bot.tick_params(labelbottom=False)

        for ip in range(n_params):
            if icol == 0:
                axes[ip * 2, 0].set_ylabel(f"param_set {ip + 1}\n" + R"$P_{\mathrm{hydro}}/P_{\mathrm{DM}}$")
                axes[ip * 2 + 1, 0].set_ylabel(R"$S_{\mathrm{CAMB}}/S_{\mathrm{pySPK}} - 1$")

    axes[0, 0].legend(fontsize=8, loc="lower left")

    fig.suptitle(
        f"SP(k) validation: CAMB Fortran (colour) vs pySPK (black dashed), SO={so}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_dir / f"spk_validation_SO{so}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _write_summary(data, out_dir):
    """Write JSON summary of max/mean/p95 relative errors across all combinations."""
    rels = []
    for (_rk, _so, _ip, _z), (k, _sup_camb, _sup_ref, rel) in data.items():
        finite = rel[np.isfinite(rel)]
        if len(finite):
            rels.append(finite)
    all_rel = np.concatenate(rels) if rels else np.array([])
    summary = {
        "evaluations": len(data),
        "max_rel_err": float(np.max(np.abs(all_rel))) if len(all_rel) else float("nan"),
        "mean_rel_err": float(np.mean(np.abs(all_rel))) if len(all_rel) else float("nan"),
        "p95_rel_err": float(np.quantile(np.abs(all_rel), 0.95)) if len(all_rel) else float("nan"),
    }
    json_path = out_dir / "spk_summary.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return json_path, summary


def main():
    parser = argparse.ArgumentParser(description="Generate SP(k) validation plots and summary.")
    parser.add_argument("--out-dir", default="spk_validation")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    relation_kinds = [1, 2, 3]
    sos = [200, 500]
    redshifts = [0.125, 0.5, 1.0, 2.0, 3.0]

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        data = compute_suppression_data(relation_kinds, sos, redshifts)

    for so in sos:
        _plot_validation(data, relation_kinds, redshifts, so=so, out_dir=out_dir)

    json_path, summary = _write_summary(data, out_dir)

    print(f"Wrote plots to {out_dir}")
    print(f"Wrote {json_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
