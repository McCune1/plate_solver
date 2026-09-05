# -*- coding: utf-8 -*-
"""Build the JSV graphical abstract (standalone PDF + PNG).

Reads the published MAC-calibration sidecar in paper/paper_figures/ so the
strip plot is the same 24 candidates as Figure 3(a). Residual-screen ranges
are the flexural values stated in the paper §5. No plate_solver import.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Wedge
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
SIDECAR = ROOT / "paper" / "paper_figures" / "mac_separation_summary.txt"
OUT_PDF = ROOT / "paper" / "GRAPHICAL_ABSTRACT.pdf"
OUT_PNG = ROOT / "paper" / "GRAPHICAL_ABSTRACT.png"


def _load_calibration(path: Path):
    reals, arts, geoms = [], [], []
    in_a = False
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("Panel (a):"):
                in_a = True
                continue
            if line.startswith("Panel (b):"):
                break
            if not in_a or line.startswith("geometry") or not line.strip():
                continue
            geom, _cand, _om, label, mac, _ok = line.rstrip("\n").split("\t")
            mac = float(mac)
            if label == "REAL":
                reals.append((geom, mac))
            else:
                arts.append((geom, mac))
            if geom not in geoms:
                geoms.append(geom)
    return geoms, reals, arts


def _draw_sector(ax):
    ax.set_xlim(-0.35, 1.55)
    ax.set_ylim(-1.05, 1.12)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Completely free annular sector", fontsize=10, pad=4)

    ri, ro = 0.42, 1.05
    th1, th2 = -45.0, 45.0
    ax.add_patch(Wedge((0, 0), ro, th1, th2, width=ro - ri,
                       facecolor="#cfd8dc", edgecolor="none"))
    ax.add_patch(Wedge((0, 0), ro, th1, th2, width=0.0,
                       facecolor="none", edgecolor="0.15", lw=1.4))
    ax.add_patch(Wedge((0, 0), ri, th1, th2, width=0.0,
                       facecolor="none", edgecolor="0.15", lw=1.4))
    t1, t2 = np.deg2rad(th1), np.deg2rad(th2)
    ax.plot([ri * np.cos(t1), ro * np.cos(t1)],
            [ri * np.sin(t1), ro * np.sin(t1)], color="0.15", lw=1.6)
    ax.plot([ri * np.cos(t2), ro * np.cos(t2)],
            [ri * np.sin(t2), ro * np.sin(t2)], color="0.15", lw=1.6)
    ax.plot([0, ro + 0.08], [0, 0], ls="--", color="0.55", lw=0.8)
    ax.annotate(r"$\theta=+\Theta$ free", xy=(ro * np.cos(t2), ro * np.sin(t2)),
                xytext=(0.78, 0.82), fontsize=8,
                arrowprops=dict(arrowstyle="-", color="0.3", lw=0.7))
    ax.annotate(r"$\theta=-\Theta$ free", xy=(ro * np.cos(t1), ro * np.sin(t1)),
                xytext=(0.78, -0.92), fontsize=8,
                arrowprops=dict(arrowstyle="-", color="0.3", lw=0.7))
    ax.text(0.22, 0.04, r"$r=R_i,R_o$" + "\nexact", fontsize=8, color="0.2",
            ha="center", va="bottom")
    ax.text(0.95, 0.04, r"$2\Theta$", fontsize=8, ha="left", color="0.2")


def _draw_residual(ax):
    ax.set_title("Flexural residual screen (portable)", fontsize=10, pad=4)
    ax.set_xscale("log")
    ax.set_xlim(8e-4, 2e3)
    ax.set_ylim(-0.55, 1.55)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["spurious zero", "physical mode"], fontsize=8)
    ax.set_xlabel(r"edge residual ratio $r_M$", fontsize=8)
    ax.tick_params(axis="x", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # All-basis envelope from §5 (n=16–44): no overlap, gap ~300×.
    # Production n=20 ranges sit inside that envelope and give >1000×.
    phys_all, spur_all = (0.002, 0.112), (37.0, 505.0)
    phys_n20, spur_n20 = (0.0025, 0.037), (74.9, 348.0)
    ax.add_patch(FancyBboxPatch((phys_all[0], 0.78), phys_all[1] - phys_all[0], 0.44,
                                boxstyle="round,pad=0.0,rounding_size=0.0",
                                facecolor="#2c7bb6", edgecolor="none",
                                alpha=0.35, transform=ax.transData))
    ax.add_patch(FancyBboxPatch((phys_n20[0], 0.86), phys_n20[1] - phys_n20[0], 0.28,
                                boxstyle="round,pad=0.0,rounding_size=0.0",
                                facecolor="#2c7bb6", edgecolor="none",
                                alpha=0.95, transform=ax.transData))
    ax.add_patch(FancyBboxPatch((spur_all[0], -0.22), spur_all[1] - spur_all[0], 0.44,
                                boxstyle="round,pad=0.0,rounding_size=0.0",
                                facecolor="#d7191c", edgecolor="none",
                                alpha=0.35, transform=ax.transData))
    ax.add_patch(FancyBboxPatch((spur_n20[0], -0.14), spur_n20[1] - spur_n20[0], 0.28,
                                boxstyle="round,pad=0.0,rounding_size=0.0",
                                facecolor="#d7191c", edgecolor="none",
                                alpha=0.95, transform=ax.transData))
    ax.text(np.sqrt(phys_all[0] * phys_all[1]), 1.28, r"$0.002$–$0.11$ ($n=16$–$44$)",
            ha="center", va="bottom", fontsize=7.5, color="#2c7bb6")
    ax.text(np.sqrt(spur_all[0] * spur_all[1]), 0.32, r"$37$–$505$ ($n=16$–$44$)",
            ha="center", va="bottom", fontsize=7.5, color="#d7191c")
    ax.annotate("no overlap across $n=16$–$44$\n$>1000\\times$ at production $n=20$",
                xy=(1.5, 0.5), xytext=(1.5, 0.5),
                ha="center", va="center", fontsize=7.5, color="0.25")


def _draw_mac(ax, geoms, reals, arts):
    ax.set_title("Extensional MAC bar (threshold 0.744)", fontsize=10, pad=4)
    xpos = {g: i for i, g in enumerate(geoms)}
    rng = np.random.default_rng(0)
    for geom, mac in reals:
        x = xpos[geom] + rng.uniform(-0.12, 0.12)
        ax.scatter(x, mac, s=28, c="#2c7bb6", marker="o",
                   edgecolors="0.15", linewidths=0.4, zorder=3)
    for geom, mac in arts:
        x = xpos[geom] + rng.uniform(-0.12, 0.12)
        ax.scatter(x, mac, s=32, c="#d7191c", marker="x",
                   linewidths=1.2, zorder=3)
    ax.axhline(0.744, color="0.2", ls="--", lw=1.0, zorder=2)
    ax.text(2.48, 0.76, "0.744", fontsize=7.5, va="bottom", ha="left")
    ax.set_ylim(-0.05, 1.08)
    ax.set_xlim(-0.45, 2.55)
    ax.set_xticks(range(len(geoms)))
    ax.set_xticklabels([g.replace(" ", "") for g in geoms], fontsize=7.5)
    ax.set_ylabel("MAC", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#2c7bb6",
               markeredgecolor="0.15", markersize=6, label="REAL-like"),
        Line2D([0], [0], marker="x", color="#d7191c", markersize=6,
               linestyle="none", label="ARTIFACT-like"),
    ], loc="center right", fontsize=7.5, frameon=False, bbox_to_anchor=(1.02, 0.55))


def main():
    geoms, reals, arts = _load_calibration(SIDECAR)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    })
    fig = plt.figure(figsize=(13.28, 5.31), facecolor="white")
    fig.suptitle(
        "Exact wave-function formulation for completely free annular sector plates",
        fontsize=13, fontweight="bold", y=0.98,
    )
    fig.text(
        0.5, 0.915,
        "Seok–Tiersten boundary determinant extended to FFFF  ·  "
        "spurious zeros discriminated, not ignored",
        ha="center", va="top", fontsize=9, color="0.25",
    )

    gs = fig.add_gridspec(
        1, 3, left=0.03, right=0.98, top=0.84, bottom=0.18,
        wspace=0.28,
    )
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    _draw_sector(ax0)
    _draw_residual(ax1)
    _draw_mac(ax2, geoms, reals, arts)

    fig.text(
        0.5, 0.075,
        "Flexure: two-sided residual screen, FE 0.64% (isotropic) / 1.08% (orthotropic), "
        "6 of 6 literature targets at a third geometry.\n"
        "Extension: 26 solver roots covering 27 FE modes (worst unique-partner 0.052%); "
        "MAC bar classifies 24 calibration candidates and all 24 sweep keys; "
        "two independent semi-analytical tables.",
        ha="center", va="center", fontsize=8.2, color="0.15",
        linespacing=1.45,
    )
    fig.text(
        0.5, 0.018,
        "McCune & Stutts  ·  graphical abstract (not embedded in the manuscript)",
        ha="center", va="bottom", fontsize=7, color="0.45",
    )
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_PDF}")
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
