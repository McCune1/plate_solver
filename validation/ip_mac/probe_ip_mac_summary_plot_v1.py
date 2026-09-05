# -*- coding: utf-8 -*-
"""
probe_ip_mac_summary_plot_v1.py -- PLOTTING PASS ONLY. Reads existing
validation/ip_mac/*.out cluster logs and renders the two-sided extensional
MAC-separation summary figure for main-text Sec.~6.8 (sec:ipmacbar). Makes
no solver call, writes no checkpoint, computes no new Omega, and has no
SOLVER_VERSION implication -- pure post-processing of numbers already
printed in the repository's own run logs.

CONTEXT: punch-list item 3 (paper review 2026-08-20) -- the subsection's
strongest quantitative result (the (0.5801, 0.9079) window that
"classifies all 24 candidates correctly") is currently a table of three
(min REAL, max ARTIFACT, gap) triples plus a second table pooling the
24-key nu=0.35 sweep by sector angle. This script turns both tables into
one two-panel strip-plot figure: every calibration candidate's own MAC
value (not just the extrema) and every tight-matched sweep candidate's
own MAC value, both against the production threshold 0.744.

PANEL (a) -- calibration: the three fully-adjudicated geometries backing
the main-text three-row table, one point per candidate (8 candidates x 3
geometries = 24 points), parsed from the SUMMARY block of each
probe_ip_block_selection_rule_v1_*.out log (RULE_LOCALMIN_RESID d=0.001
column -- the same column the paper's floor/ceiling numbers are quoted
from).

PANEL (b) -- sweep: the 24-key nu=0.35 geometry sweep backing the
main-text six-row table (pooled by 2*Theta/pi, matching that table's own
pooling convention), one point per TIGHT (FE match < 0.5%) candidate,
parsed from the SUMMARY block of each probe_ip_rank14_mac_bar_v1_*.out
log (the "max-sc" column -- the same same-class-restricted maximum the
production CALL is scored from, not the unrestricted "MAC@n" column).
Non-tight candidates are excluded from the plot exactly as they are
excluded from the paper's own pooled REAL/ART counts.

INPUT: LOGDIR env var (default: the directory this script lives in),
expected to contain the 3 calibration logs and the 24 sweep logs named
per validation/ip_mac/README.md. Missing files are reported and skipped,
not fatal -- a partial run (e.g. only the 3 calibration logs present)
still produces panel (a) alone.

OUTPUT: OUTDIR env var (default: LOGDIR/../../paper/paper_figures),
writing mac_separation_summary.pdf (for \\includegraphics), a matching
.png preview, and a mac_separation_summary.txt sidecar dumping every
parsed (geometry, candidate/Omega, label, MAC) tuple, mirroring the
existing sigmin_trace_2p55_2p90.{pdf,png,txt} convention in that
directory.

Standalone: no plate_solver import, no cluster submission. Needs only
matplotlib (Agg backend) and the standard library. Runs in seconds.
"""
import os
import re
import sys

MAC_THRESHOLD = 0.744
CAL_WINDOW = (0.5801, 0.9079)  # combined REAL floor / ARTIFACT ceiling, Sec 6.8

HERE = os.path.dirname(os.path.abspath(__file__))
LOGDIR = os.environ.get("LOGDIR", HERE)
OUTDIR = os.environ.get(
    "OUTDIR", os.path.normpath(os.path.join(LOGDIR, "..", "..", "paper", "paper_figures"))
)

# --- calibration logs (Sec 6.8's three-row table) --------------------------
CAL_FILES = [
    ("ffp1", "probe_ip_block_selection_rule_v1_ffp1_2431872.out"),
    ("r200a100", "probe_ip_block_selection_rule_v1_r200a100_2431873.out"),
    ("r200a50", "probe_ip_block_selection_rule_v1_r200a50_2432028.out"),
]

# --- sweep logs (Sec 6.8's six-row, angle-pooled table) ---------------------
SWEEP_FILES = [
    "probe_ip_rank14_mac_bar_v1_2432508.out",
    "probe_ip_rank14_mac_bar_v1_2433081.out",
    "probe_ip_rank14_mac_bar_v1_2433083.out",
    "probe_ip_rank14_mac_bar_v1_2433085.out",
] + [f"probe_ip_rank14_mac_bar_v1_{job}.out" for job in range(2433497, 2433517)]


def _slice_block(lines, start_pat, end_pat):
    """Returns the lines strictly between the first line matching
    start_pat and the first subsequent line matching end_pat (both
    matched with `in`, not full regex -- the logs' section banners are
    plain substrings). Returns [] if either boundary is not found, so a
    malformed/truncated log degrades to "nothing parsed", not a crash."""
    started = False
    out = []
    for ln in lines:
        if not started:
            if start_pat in ln:
                started = True
            continue
        if end_pat in ln:
            break
        out.append(ln)
    return out


def parse_calibration_log(path, geom_key):
    """Parses one probe_ip_block_selection_rule_v1_*.out SUMMARY block.
    Each row is whitespace-delimited:
      cand  label  sigma_blk  resid_blk  q#:MAC:C|W  [q#:MAC:C|W ...]
    The RULE_LOCALMIN_RESID d=0.001 column (first "rule d=..." column) is
    used -- the same one Sec 6.8's floor/ceiling numbers are quoted from.
    Returns a list of dicts: geom, cand, Omega, label, mac, correct."""
    with open(path, "r") as f:
        lines = f.readlines()
    header = "".join(lines[:20])
    m = re.search(r"r0/2b=([\d.]+)\s+2Theta/pi=([\d.]+)", header)
    geom_label = f"({m.group(1)}, {m.group(2)})" if m else geom_key

    block = _slice_block(
        lines,
        "SUMMARY -- rule scoring against published labels",
        "READING (per pre-registered criteria",
    )
    rows = []
    for ln in block:
        toks = ln.split()
        if len(toks) < 6:
            continue
        cand, label = toks[0], toks[1]
        if label not in ("REAL", "ARTIFACT"):
            continue
        rule001 = toks[4]  # "q#:MAC:C" or "q#:MAC:W"
        parts = rule001.split(":")
        if len(parts) != 3:
            continue
        try:
            mac = float(parts[1])
        except ValueError:
            continue
        correct = parts[2] == "C"
        om_m = re.search(r"_(\d+p\d+)$", cand)
        omega = float(om_m.group(1).replace("p", ".")) if om_m else None
        rows.append(dict(geom=geom_label, cand=cand, Omega=omega, label=label,
                          mac=mac, correct=correct))
    return geom_label, rows


def parse_sweep_log(path):
    """Parses one probe_ip_rank14_mac_bar_v1_*.out SUMMARY block. Each row:
      Om  Hz  blk  near  rel%  MAC@n  max-sc  at  tight  call
    Only tight==1 rows are returned (matching the paper's own pooled
    REAL/ART counts, which are counted over tight matches only), using
    the max-sc column (the same-class-restricted maximum the production
    CALL/threshold-0.744 decision is actually made from). Returns
    (two_theta_pi, [rows]) or (None, []) if the header can't be read."""
    with open(path, "r") as f:
        lines = f.readlines()
    header = "".join(lines[:10])
    m = re.search(r"2Theta/pi=([\d.]+)", header)
    two_theta_pi = float(m.group(1)) if m else None
    ratio_m = re.search(r"r0/2b=([\d.]+)", header)
    r0_2b = float(ratio_m.group(1)) if ratio_m else None

    block = _slice_block(lines, "SUMMARY", "READING (this key only")
    rows = []
    for ln in block:
        toks = ln.split()
        if len(toks) != 10:
            continue
        try:
            omega = float(toks[0])
            maxsc = float(toks[6])
            tight = toks[8]
            call = toks[9]
        except ValueError:
            continue
        if call not in ("REAL-like", "ARTIFACT-like"):
            continue
        if tight != "1":
            continue
        label = "REAL" if call == "REAL-like" else "ARTIFACT"
        rows.append(dict(r0_2b=r0_2b, two_theta_pi=two_theta_pi, Omega=omega,
                          label=label, mac=maxsc))
    return two_theta_pi, rows


def load_all():
    cal_groups = []  # [(geom_label, [rows])]
    for geom_key, fname in CAL_FILES:
        path = os.path.join(LOGDIR, fname)
        if not os.path.isfile(path):
            print(f"  [calibration] MISSING: {fname}", flush=True)
            continue
        geom_label, rows = parse_calibration_log(path, geom_key)
        print(f"  [calibration] {fname}: {len(rows)} candidates ({geom_label})",
              flush=True)
        cal_groups.append((geom_label, rows))

    sweep_by_angle = {}  # two_theta_pi -> [rows]
    n_sweep_files = 0
    n_sweep_missing = 0
    for fname in SWEEP_FILES:
        path = os.path.join(LOGDIR, fname)
        if not os.path.isfile(path):
            n_sweep_missing += 1
            continue
        n_sweep_files += 1
        two_theta_pi, rows = parse_sweep_log(path)
        if two_theta_pi is None:
            print(f"  [sweep] {fname}: header not recognized, skipped", flush=True)
            continue
        sweep_by_angle.setdefault(two_theta_pi, []).extend(rows)
    print(f"  [sweep] {n_sweep_files} of {len(SWEEP_FILES)} log files found "
          f"({n_sweep_missing} missing)", flush=True)
    for angle in sorted(sweep_by_angle):
        rows = sweep_by_angle[angle]
        n_real = sum(1 for r in rows if r["label"] == "REAL")
        n_art = sum(1 for r in rows if r["label"] == "ARTIFACT")
        print(f"    2Theta/pi={angle:.2f}: {len(rows)} tight candidates "
              f"({n_real} REAL / {n_art} ARTIFACT)", flush=True)
    return cal_groups, sweep_by_angle


def make_figure(cal_groups, sweep_by_angle, out_base):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["pdf.fonttype"] = 42  # embed real (not Type-3) fonts
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.size"] = 9

    real_style = dict(marker="o", s=32, facecolors="none",
                       edgecolors="#1a7a1a", linewidths=1.2, alpha=0.65,
                       label="REAL-like")
    art_style = dict(marker="x", s=34, color="#b3261e", linewidths=1.4,
                      alpha=0.55, label="ARTIFACT-like")

    have_cal = any(rows for _, rows in cal_groups)
    have_sweep = any(sweep_by_angle.values())
    n_panels = int(have_cal) + int(have_sweep)
    if n_panels == 0:
        print("Nothing parsed from any log -- no figure written.", flush=True)
        return None

    fig, axes = plt.subplots(1, n_panels, figsize=(5.2 * n_panels, 3.4),
                              squeeze=False)
    axes = axes[0]
    ax_i = 0

    def _spread(n, center, half_width=0.11):
        """n deterministic x-offsets from center, evenly fanned out over
        [center-half_width, center+half_width] (n=1 -> just center)."""
        if n <= 1:
            return [center]
        step = (2 * half_width) / (n - 1)
        return [center - half_width + k * step for k in range(n)]

    if have_cal:
        ax = axes[ax_i]; ax_i += 1
        xt_labels = [g for g, _ in cal_groups]
        for gi, (geom_label, rows) in enumerate(cal_groups):
            real_rows = [r for r in rows if r["label"] == "REAL"]
            art_rows = [r for r in rows if r["label"] == "ARTIFACT"]
            if real_rows:
                xs = _spread(len(real_rows), gi - 0.16)
                kw = dict(real_style)
                if gi != 0:
                    kw.pop("label", None)
                ax.scatter(xs, [r["mac"] for r in real_rows], **kw)
            if art_rows:
                xs = _spread(len(art_rows), gi + 0.16)
                kw = dict(art_style)
                if gi != 0:
                    kw.pop("label", None)
                ax.scatter(xs, [r["mac"] for r in art_rows], **kw)
        ax.axhspan(CAL_WINDOW[0], CAL_WINDOW[1], color="#f0c419", alpha=0.15,
                   zorder=0)
        ax.axhline(MAC_THRESHOLD, color="black", linestyle="--", linewidth=1.0,
                   zorder=1)
        ax.text(len(cal_groups) - 0.52, MAC_THRESHOLD + 0.02,
                f"threshold {MAC_THRESHOLD}", fontsize=7.5, va="bottom")
        ax.set_xticks(range(len(cal_groups)))
        ax.set_xticklabels([f"$r_0/2b,\\ 2\\Theta/\\pi$\n{g}" for g in xt_labels],
                            fontsize=7.5)
        ax.set_ylabel("MAC")
        ax.set_ylim(-0.03, 1.03)
        ax.set_title("(a) three-geometry calibration, "
                      f"{sum(len(r) for _, r in cal_groups)} candidates",
                      fontsize=9)
        ax.legend(loc="center left", fontsize=7.5, framealpha=0.9)

    if have_sweep:
        ax = axes[ax_i]; ax_i += 1
        angles = sorted(sweep_by_angle)
        for gi, angle in enumerate(angles):
            rows = sweep_by_angle[angle]
            real_rows = [r for r in rows if r["label"] == "REAL"]
            art_rows = [r for r in rows if r["label"] == "ARTIFACT"]
            if real_rows:
                xs = _spread(len(real_rows), gi - 0.18, half_width=0.13)
                kw = dict(real_style)
                if gi != 0:
                    kw.pop("label", None)
                ax.scatter(xs, [r["mac"] for r in real_rows], **kw)
            if art_rows:
                xs = _spread(len(art_rows), gi + 0.18, half_width=0.13)
                kw = dict(art_style)
                if gi != 0:
                    kw.pop("label", None)
                ax.scatter(xs, [r["mac"] for r in art_rows], **kw)
        ax.axhline(MAC_THRESHOLD, color="black", linestyle="--", linewidth=1.0,
                   zorder=1)
        ax.text(len(angles) - 0.55, MAC_THRESHOLD + 0.02,
                f"threshold {MAC_THRESHOLD}", fontsize=7.5, va="bottom")
        ax.set_xticks(range(len(angles)))
        ax.set_xticklabels([f"{a:.2f}" for a in angles], fontsize=8)
        ax.set_xlabel("$2\\Theta/\\pi$ (pooled over 4 radius ratios)")
        ax.set_ylabel("MAC")
        ax.set_ylim(-0.03, 1.03)
        n_tight = sum(len(r) for r in sweep_by_angle.values())
        ax.set_title(f"(b) 24-key $\\nu{{=}}0.35$ sweep, {n_tight} tight "
                      "candidates", fontsize=9)
        # upper-right is empty here: the 1.25/1.50 keys have no REAL-like
        # tight candidate at all, so nothing sits above the threshold there
        ax.legend(loc="upper right", fontsize=7.5, framealpha=0.9)

    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    pdf_path = out_base + ".pdf"
    png_path = out_base + ".png"
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    return pdf_path, png_path


def write_sidecar(cal_groups, sweep_by_angle, out_base):
    with open(out_base + ".txt", "w") as f:
        f.write("mac_separation_summary -- parsed source data\n")
        f.write("=" * 70 + "\n")
        f.write("Panel (a): calibration (RULE_LOCALMIN_RESID d=0.001 MAC)\n")
        f.write("geometry\tcandidate\tOmega\tlabel\tmac\tcorrect\n")
        for geom_label, rows in cal_groups:
            for r in rows:
                f.write(f"{geom_label}\t{r['cand']}\t{r['Omega']}\t"
                        f"{r['label']}\t{r['mac']:.4f}\t{r['correct']}\n")
        f.write("\nPanel (b): sweep (tight matches only, same-class max MAC)\n")
        f.write("2Theta/pi\tr0/2b\tOmega\tlabel\tmac\n")
        for angle in sorted(sweep_by_angle):
            for r in sweep_by_angle[angle]:
                f.write(f"{angle}\t{r['r0_2b']}\t{r['Omega']}\t"
                        f"{r['label']}\t{r['mac']:.4f}\n")


def main():
    print("=" * 78)
    print("  probe_ip_mac_summary_plot_v1 -- start", flush=True)
    print(f"  LOGDIR={LOGDIR}")
    print(f"  OUTDIR={OUTDIR}")
    print("=" * 78, flush=True)

    cal_groups, sweep_by_angle = load_all()

    out_base = os.path.join(OUTDIR, "mac_separation_summary")
    result = make_figure(cal_groups, sweep_by_angle, out_base)
    if result is None:
        raise SystemExit(1)
    pdf_path, png_path = result
    write_sidecar(cal_groups, sweep_by_angle, out_base)

    print(f"\nWrote {pdf_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {out_base}.txt (sidecar, parsed source data)")
    print("\nPlotting pass only -- no solver call, no checkpoint, no "
          "SOLVER_VERSION action, no package changes.", flush=True)


if __name__ == "__main__":
    main()
