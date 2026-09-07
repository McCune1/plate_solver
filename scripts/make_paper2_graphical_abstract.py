# -*- coding: utf-8 -*-
"""Standalone JSV graphical abstract for Paper 2 from existing figures.

Left: four-corner checkerboard (Fig. 1). Right: MAC overlay at ell/b=1.5
first SYM (Fig. 2). No plate_solver import; no new compute.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.image import imread

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "paper" / "paper_figures"
OUT_PDF = ROOT / "paper" / "PAPER2_GRAPHICAL_ABSTRACT.pdf"
OUT_PNG = ROOT / "paper" / "PAPER2_GRAPHICAL_ABSTRACT.png"


def main():
    left = imread(FIG / "rect_ffff_corner_checkerboard.png")
    right = imread(FIG / "rect_ff_oop_mac_lb15_sym0964.png")

    fig = plt.figure(figsize=(7.2, 3.15), dpi=200)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.35], wspace=0.08)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax0.imshow(left)
    ax1.imshow(right)
    for ax in (ax0, ax1):
        ax.set_axis_off()
    ax0.set_title("FFFF rectangle: Kirchhoff checkerboard", fontsize=8, pad=3)
    ax1.set_title(r"First SYM, $\ell/b=1.5$: FE $U_Z$ vs $w$, MAC $=1.000$",
                  fontsize=8, pad=3)
    fig.suptitle(
        "Exact wave-function FFFF rectangular plates  |  Screen B + unre-tuned MAC 0.744",
        fontsize=9, y=0.98,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.82, bottom=0.02)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("wrote", OUT_PDF)
    print("wrote", OUT_PNG)


if __name__ == "__main__":
    main()
