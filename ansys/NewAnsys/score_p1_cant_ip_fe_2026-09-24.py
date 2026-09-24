#!/usr/bin/env python3
# Scorer for ansys_p1_cant_ip_r125_{a025,a050,a100}_2026-09-24.inp (Paper 1
# cantilever in-plane rows, r0/(2b)=1.25, 2Theta/pi = 0.25/0.50/1.00).
# Python 3.6-safe. Reads raw Hz only, converts here (never a deck-computed column).
#
# PRE-REGISTERED 2026-09-24, before any ANSYS output existed. (A first draft,
# a025 only, used the SM Table S.9 mode-1 value 0.339343 as a control. An
# in-house Q9 plane-stress FE run the same day -- probe_p1_cant_ip_q9fe_2026-09-24.py,
# Richardson over 16/32/64, validated against PLANE183 cc_ip_mesh96 to 0.005% --
# put that mode at 0.348738, so the control would have failed for a reason that
# has nothing to do with FE correctness. The ANSYS run's job is to CONFIRM the
# in-house FE reference for every in-plane row.)
#  F0 mesh:   |f96-f192|/f192 <= 0.05% for modes 1-4 of each deck.
#  F1 limit:  Richardson from 48/96/192 if monotone, else f192.
#  F2 decisive: |Omega_ANSYS - Omega_Q9| <= 0.05% for modes 1-4 of each deck.
#     PASS -> the Q9 values are externally confirmed and become the rows' FE
#     reference. FAIL -> one of the two FE models is wrong; use neither.
#  Report only (no gate): printed (production), Seok-Tiersten and SM Table S.9
#     "converged" values against Omega_ANSYS.
import math, os, re
WBAR = 2480.106584155221   # (pi/2b) sqrt(c66/rho), b=2, E=210e9, nu=.35, rho=7800
ROWS = {  # tag: (printed, Seok-Tiersten, S.9 converged (None = n.r.), Q9 modes 1-4)
 "a025": ([0.341243, 0.808929, 0.940032], [0.34589, 0.81921, 0.94012], [0.339343, None, None],
          [0.348738, 0.823991, 0.940601, 1.551622]),
 "a050": ([0.120677, 0.347033, 0.529239], [0.12088, 0.35191, 0.52930], [0.120472, 0.350357, 0.529361],
          [0.121449, 0.353054, 0.529883, 0.943960]),
 "a100": ([0.039466, 0.106529, 0.263049], [0.03948, 0.10556, 0.26129], [0.039520, 0.106136, 0.260998],
          [0.039555, 0.105826, 0.261755, 0.458798]),
}
def read(tag, n):
    fn = "p1_cant_ip_r125_%s_mesh%d.txt" % (tag, n)
    if not os.path.exists(fn): return None
    f = []
    for line in open(fn):
        m = re.match(r"\s*(\d+)\.?\s+([-+0-9.Ee]+)", line)
        if m: f.append(float(m.group(2)))
    return f if len(f) >= 4 else None
allok = True; anyrun = False
for tag in ("a025", "a050", "a100"):
    pr, st, s9, q9 = ROWS[tag]
    M = {n: read(tag, n) for n in (48, 96, 192)}
    print("==== %s" % tag)
    if any(v is None for v in M.values()):
        print("  NOT_RUN: missing", [n for n, v in M.items() if v is None]); allok = False; continue
    anyrun = True
    for n in (48, 96, 192): print("  mesh%-3d" % n, " ".join("%.6f" % x for x in M[n][:5]))
    om = []
    for k in range(4):
        a, b, c = M[48][k], M[96][k], M[192][k]
        d = abs(b - c) / c * 100; f0 = d <= 0.05
        if (a - b) * (b - c) > 0 and b != c:
            p = math.log(abs(a - b) / abs(b - c)) / math.log(2); finf = c + (c - b) / (2 ** p - 1); how = "p=%.2f" % p
        else:
            finf = c; how = "f192"
        O = 2 * math.pi * finf / WBAR; om.append(O)
        e = (O - q9[k]) / q9[k] * 100; f2 = abs(e) <= 0.05
        allok &= f0 and f2
        print("  mode %d: F0 %.4f%% %s | Omega_ANSYS %.6f (%s) vs Q9 %.6f %+.4f%% F2 %s"
              % (k + 1, d, "PASS" if f0 else "FAIL", O, how, q9[k], e, "PASS" if f2 else "FAIL"))
    for k in range(3):
        O = om[k]
        s9s = ("%+.3f%%" % ((s9[k] - O) / O * 100)) if s9[k] else "n.r."
        print("  row m%d vs ANSYS: printed %+.3f%%  Seok-Tiersten %+.3f%%  S.9 %s"
              % (k + 1, (pr[k] - O) / O * 100, (st[k] - O) / O * 100, s9s))
print("VERDICT:", ("CONFIRMED: Q9 reference stands for all in-plane rows" if (allok and anyrun)
                   else ("NOT_RUN" if not anyrun else "UNRESOLVED -- read F0/F2 above")))
