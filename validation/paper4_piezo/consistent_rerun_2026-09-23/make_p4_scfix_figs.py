"""Regenerate the Paper 4 figures affected by the 2026-09-23 consistent-
projection / reciprocal-charge fix (LESSONS Sec 18.231):
piezo_p4_admittance_zoom.png (driven_ff_qv, h1/2h=1/12) and
piezo_p4_lever16_noise_floor.png (from the rerun JSON)."""
import json, math, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.environ.get("PKG_PATH", "."))
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver as P

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
H = 0.01
s = P(r_i=0.1, r_o=0.6, h=H, E=200e9, nu=0.3, rho=7800.0, h1=2 * H / 12,
      C11E=132e9, C12E=71e9, C13E=73e9, C33E=115e9, e31=4.1, e33=14.1,
      X11=7.124e-9, X33=5.841e-9, rho_pzt=7500.0, dps=30)
fr = s.elastic_bisect(720.0, 780.0, 0) / (2 * math.pi)
fa = s.oc_ff_bisect(720.0, 800.0, 0) / (2 * math.pi)
f = np.linspace(118.6, 120.0, 281)
c = np.array([s.driven_ff_qv(2 * math.pi * x, 0)["q_over_v"] for x in f]) / s.clamped_C0()
c[np.abs(c) > 50] = np.nan
fig, ax = plt.subplots(figsize=(7.2, 4.0))
ax.plot(f, c, color="k", lw=1.4)
ax.axhline(1.0, color="0.6", ls="--", lw=0.9)
ax.axhline(0.0, color="0.5", lw=0.8)
ax.axvline(fr, color="#1f77b4", ls=":", lw=1.3, label=r"$f_r$ (SC)")
ax.axvline(fa, color="#d62728", ls=":", lw=1.3, label=r"$f_a$ (OC)")
ax.set_xlim(118.6, 120.0); ax.set_ylim(-6, 6)
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel(r"$C_{\mathrm{eff}}/C_0$")
ax.legend(frameon=False, fontsize=9, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "piezo_p4_admittance_zoom.png"), dpi=140)
plt.close()
print("zoom: f_r=%.6f Hz f_a=%.6f Hz" % (fr, fa))

d = json.load(open(os.environ.get("LEVER16_JSON", "piezo_p4_lever16_freq_noise_floor_results.json")))
fig, ax = plt.subplots(figsize=(7.2, 4.5))
for key, lab, col, mk in (("ff", r"F–F open-circuit $e_{31}$", "#b2182b", "o"),
                          ("cf", r"C–F open-circuit $e_{31}$", "#2166ac", "^")):
    cur = [r for r in d[key]["curve"] if not r["is_lab_floor"]]
    ax.loglog([r["rel_prec_pct"] for r in cur], [r["rel_unc_pct"] for r in cur],
              marker=mk, color=col, lw=1.6, ms=7, label=lab)
lab = 100 * d["rel_prec_lab"]
ax.axvline(lab, color="#1a8a3a", ls="--", lw=1.6)
ax.axvline(0.01, color="0.7", ls=":", lw=0.9)
ax.text(lab * 1.07, 0.55, "Boeringa–McCune\nimpedance/SLDV floor\n(%.3f%%)" % lab,
        color="#1a8a3a", fontsize=9)
ax.set_xlabel("assumed relative frequency error (%)")
ax.set_ylabel(r"propagated relative uncertainty in $e_{31}$ (%)")
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.7)
ax.legend(frameon=False, fontsize=10, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "piezo_p4_lever16_noise_floor.png"), dpi=150)
plt.close()
print("noise floor written")
