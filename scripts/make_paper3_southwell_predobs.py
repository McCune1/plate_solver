# -*- coding: utf-8 -*-
"""Southwell Table 2.31 predicted vs observed relative error (Paper 3 SM).

Eighteen points, log-log, identity line, MATCH/WEAK/MISS coded.
Observed errors for the three window-endpoint rows are the wide re-search
values in Table S.5 (daggers), not the stored boundary lambda^2.
Numbers are locked to PAPER3_RING_DISK_SUPPLEMENTARY.tex Table S.5 /
p3_phaseD_inverse_ba_close.json P3 with Sec. 18.120 corrections.

No plate_solver import.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "paper" / "paper_figures"
OUT_PDF = FIG / "southwell_pred_vs_obs.pdf"
OUT_PNG = FIG / "southwell_pred_vs_obs.png"
OUT_TXT = FIG / "southwell_pred_vs_obs.txt"

# n, b/a, printed lambda^2, pred, obs, gate, dagger (wide re-search obs)
ROWS = [
    (0, 0.276, 6.25, 0.001206, 0.00119, "MATCH", False),
    (0, 0.642, 25.0, 0.02959, 0.02981, "WEAK", False),
    (0, 0.840, 81.0, 0.4472, 0.6419, "MISS", True),
    (1, 0.060, 2.82, 0.001706, 0.001838, "MATCH", False),
    (1, 0.397, 9.00, 0.002322, 0.002299, "MATCH", False),
    (1, 0.603, 21.2, 0.002782, 0.002746, "MATCH", False),
    (1, 0.634, 25.0, 0.002768, 0.00274, "MATCH", False),
    (1, 0.771, 64.0, 0.01198, 0.01191, "WEAK", False),
    (1, 0.827, 121.0, 0.05979, 0.05642, "MISS", False),
    (2, 0.186, 6.25, 0.007324, 0.007153, "MATCH", False),
    (2, 0.349, 9.00, 0.00331, 0.003267, "MATCH", False),
    (2, 0.522, 16.0, 0.000732, 0.0007167, "MATCH", False),
    (2, 0.769, 64.0, 0.02465, 0.02473, "WEAK", False),
    (2, 0.81, 100.0, 0.03623, 0.03474, "MISS", False),
    (3, 0.43, 16.0, 0.01425, 0.01342, "WEAK", False),
    (3, 0.59, 25.0, 0.00149, 0.001459, "MATCH", False),
    (3, 0.71, 49.0, 0.08081, 0.07457, "MISS", True),
    (3, 0.82, 100.0, 0.1033, 0.1104, "MISS", True),
]

STYLE = {
    "MATCH": dict(marker="o", color="#2166ac", facecolor="#2166ac",
                  label="MATCH", zorder=3, s=38),
    "WEAK": dict(marker="^", color="#f4a582", facecolor="#f4a582",
                 label="WEAK", zorder=3, s=42),
    "MISS": dict(marker="s", color="#b2182b", facecolor="#b2182b",
                 label="MISS", zorder=3, s=36),
}


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Southwell Table 2.31 predicted vs observed relative error",
        "# pred = (d lambda^2 / d(b/a)) * |delta(b/a)| / lambda^2",
        "# obs  = |lambda^2_solver - printed| / printed  (wide re-search if dagger=1)",
        "# n  b/a  printed  pred  obs  gate  dagger",
    ]
    for n, ba, printed, pred, obs, gate, dag in ROWS:
        lines.append(
            "%d  %.3f  %g  %.6g  %.6g  %s  %d"
            % (n, ba, printed, pred, obs, gate, int(dag))
        )
    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
        "mathtext.fontset": "dejavuserif",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })

    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    lo, hi = 4.5e-4, 1.05
    xx = np.array([lo, hi])
    ax.fill_between(xx, xx / 2.0, xx * 2.0, color="0.88", lw=0, zorder=0)
    ax.plot(xx, xx, color="0.15", lw=1.0, zorder=1)
    ax.plot(xx, xx * 2.0, color="0.45", lw=0.7, ls="--", zorder=1)
    ax.plot(xx, xx / 2.0, color="0.45", lw=0.7, ls="--", zorder=1)

    seen = set()
    for n, ba, printed, pred, obs, gate, dag in ROWS:
        st = STYLE[gate]
        lab = st["label"] if (gate not in seen and not dag) else None
        if lab:
            seen.add(gate)
        ax.scatter(
            [pred], [obs], marker=st["marker"], s=st["s"],
            facecolors=st["facecolor"] if not dag else "white",
            edgecolors=st["color"], linewidths=1.15,
            label=lab, zorder=st["zorder"],
        )
    for gate in ("MATCH", "WEAK", "MISS"):
        if gate not in seen:
            st = STYLE[gate]
            ax.scatter([], [], marker=st["marker"], s=st["s"],
                       facecolors=st["facecolor"], edgecolors=st["color"],
                       linewidths=1.15, label=st["label"])

    # Load-bearing cell: n=0, b/a=0.840
    ax.annotate(
        r"$n=0,\ b/a=0.840$",
        xy=(0.4472, 0.6419), xytext=(0.018, 0.55),
        fontsize=8, color="#b2182b",
        arrowprops=dict(arrowstyle="-", color="#b2182b", lw=0.6),
        ha="left", va="center",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"predicted relative error")
    ax.set_ylabel(r"observed relative error")
    ax.legend(frameon=False, loc="upper left", fontsize=8, handletextpad=0.4)
    ax.tick_params(which="both", top=True, right=True)
    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    ratios = [obs / pred for (_, _, _, pred, obs, _, _) in ROWS]
    print("wrote", OUT_PDF)
    print("wrote", OUT_PNG)
    print("wrote", OUT_TXT)
    print("n=%d  ratio min=%.2f max=%.2f  within x2: %d/%d"
          % (len(ROWS), min(ratios), max(ratios),
             sum(0.5 <= r <= 2.0 for r in ratios), len(ROWS)))


if __name__ == "__main__":
    main()
