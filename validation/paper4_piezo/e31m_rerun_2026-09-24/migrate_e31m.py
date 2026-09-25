# Writes *_e31m_2026-09-24.py copies of the Paper 4 probes at e31 = -4.1.
import re, sys, os
SCALE = (13.0504347826 / 4.8504347826) ** 2   # ebar31 ratio squared, -4.1 vs +4.1
BANNER = ("# e31 = -4.1 COPY (LESSONS Sec 18.246) of {parent}.\n"
          "# Written by migrate_e31m.py. Changes, all fixed before running:\n"
          "#   * E31_TRUE = -4.1 (Liu/Duan Table 1 as printed, Sec 18.244);\n"
          "#   * e31 inverse bracket [1, 10] -> [-10, -1] (e31bar = e31 - 8.95 is\n"
          "#     monotone there; the reflection twin sits at 22.0);\n"
          "#   * relative uncertainties divide by |E31_TRUE|;\n"
          "#   * coupling-magnitude sanity windows (quantities ~ e31bar^2) scaled\n"
          "#     by (13.05/4.85)^2 = {scale:.3f};\n"
          "#   * output files carry the suffix _e31m.\n"
          "# Anchor gates that pin +4.1 manuscript numbers are reported against\n"
          "# those numbers and are NOT expected to pass at -4.1.\n")
COMMON = [
    ("E31_TRUE, E33 = 4.1, 14.1", "E31_TRUE, E33 = -4.1, 14.1"),
    ("E31_LO, E31_HI = 1.0, 10.0", "E31_LO, E31_HI = -10.0, -1.0"),
    ('"e31": (1.0, 10.0)', '"e31": (-10.0, -1.0)'),
    ("/ E31_TRUE * 100", "/ abs(E31_TRUE) * 100"),
    ("(delta_e31_neg_001 / E31_TRUE)", "(delta_e31_neg_001 / abs(E31_TRUE))"),
]
def sc(x):
    return "%.4g" % (x * SCALE)
SPECIFIC = {
    "probe_piezo_p4_ff_open_circuit_2026-09-16.py": [
        ("(0.002 < drel < 0.005)", "(%s < drel < %s)" % (sc(0.002), sc(0.005))),
        ("g3_pass = 0.002 < k2_el < 0.02", "g3_pass = %s < k2_el < %s" % (sc(0.002), sc(0.02))),
    ],
    "probe_piezo_p4_cf_mixed_edge_2026-09-16.py": [
        ("(0.002 < drel_cf < 0.005)", "(%s < drel_cf < %s)" % (sc(0.002), sc(0.005))),
        ("(1e-6 < drel_fc < 0.01)", "(%s < drel_fc < %s)" % (sc(1e-6), sc(0.01))),
        ("abs(drel_cf_sc) < 0.0005", "abs(drel_cf_sc) < %s" % sc(0.0005)),
        ("abs(drel_fc_sc) < 0.0005", "abs(drel_fc_sc) < %s" % sc(0.0005)),
    ],
}
for f in sorted(os.listdir(".")):
    if not (f.startswith("probe_piezo_p4_") and f.endswith("2026-09-16.py")):
        continue
    s = open(f).read()
    if "E31_TRUE" not in s:
        continue
    n = 0
    for a, b in COMMON + SPECIFIC.get(f, []):
        if a in s:
            s = s.replace(a, b); n += 1
    # outputs
    s, k1 = re.subn(r'(_results)(\.json)', r'\1_e31m\2', s)
    s, k2 = re.subn(r'("[A-Za-z0-9_]+)(\.(?:csv|png))"', r'\1_e31m\2"', s)
    s = BANNER.format(parent=f, scale=SCALE) + s
    out = f.replace("2026-09-16.py", "e31m_2026-09-24.py")
    open(out, "w", newline="\n").write(s)
    left = [l for l in s.splitlines() if re.search(r"\b4\.1\b", l) and "E31_TRUE" in l]
    print("%-58s edits=%d json=%d csvpng=%d" % (out, n, k1, k2))
