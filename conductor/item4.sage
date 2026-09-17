# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# item4.sage -- CM type when e=2: type = d - W1 + kappa - 1 ?  on every e=2 case of the archived sweep gor_sweep_OUT.txt (gor.sage list 100 20)
#               plus A_22(69); type by fresh colon lattice (lat.sage), by F_p colon (capas), and archived gor.sage value.
# Also item 5: HS 12.2.3  l(O/A) - l(A/f) >= type - 1, and almost-Gorenstein equality, on all archived rows.
# Run (in this directory): sage item4.sage
import re, sys, time
load("lat.sage")
rows = []
for ln in open("gor_sweep_OUT.txt"):
    mm = re.match(r"m=(\d+) q=(\d+) p=(\d+) E=(\d+) d=(\d+) g=(\d+) e=(\d+) W=\[([^\]]*)\] lOA=(\d+) lAC=(\d+) type=(\d+)", ln)
    if mm:
        g = mm.groups()
        rows.append(dict(m=int(g[0]), q=int(g[1]), p=int(g[2]), E=int(g[3]), d=int(g[4]), g=int(g[5]), e=int(g[6]),
                         W=[int(x) for x in g[7].split(",")], lOA=int(g[8]), lAC=int(g[9]), type=int(g[10])))
print("archived rows parsed:", len(rows))
# ---- item 5 on archived data (all e, and e=2)
for label, sel in (("all e", rows), ("e=2", [r for r in rows if r["e"] == 2])):
    viol = [r for r in sel if r["lOA"] - r["lAC"] < r["type"] - 1]
    eq = [r for r in sel if r["lOA"] - r["lAC"] == r["type"] - 1]
    gor = [r for r in sel if r["type"] == 1]
    print("[item5 archived %s] rows=%d  HS12.2.3 violations=%d  equality(almost Gorenstein, incl. Gorenstein)=%d  of which Gorenstein(type 1)=%d  strict=%d" % (
        label, len(sel), len(viol), len(eq), len(gor), len(sel) - len(eq)))
e2 = [r for r in rows if r["e"] == 2]
print("[item5] e=2 almost-Gorenstein non-Gorenstein cases:")
for r in e2:
    if r["lOA"] - r["lAC"] == r["type"] - 1 and r["type"] > 1:
        print("   m=%d q=%d p=%d d=%d W=%s lOA=%d lAf=%d type=%d" % (r["m"], r["q"], r["p"], r["d"], r["W"], r["lOA"], r["lAC"], r["type"]))
print("[item5] e=2 strict (not almost Gorenstein):")
for r in e2:
    if r["lOA"] - r["lAC"] > r["type"] - 1:
        print("   m=%d q=%d p=%d d=%d W=%s lOA=%d lAf=%d type=%d" % (r["m"], r["q"], r["p"], r["d"], r["W"], r["lOA"], r["lAC"], r["type"]))
sys.stdout.flush()
# ---- item 4
cases = [(r["m"], r["q"], r["p"], r["type"]) for r in e2] + [(22, 69, 3, None)]
fails = 0; n = 0
table = []
for (m, q, p, tarch) in cases:
    res = analyse(m, q, p, colon=True)
    if res["e"] != 2:
        print("  !! e != 2 in fresh computation:", m, q, p, res["e"]); fails += 1; continue
    pred = res["d"] - res["W1"] + res["kappa"] - 1
    ok = (res["typeC"] == pred) and (res["typeF"] == res["typeC"]) and (tarch is None or tarch == res["typeC"]) and res["W"][1] == res["W1"]
    n += 1; fails += (not ok)
    table.append((m, q, p, res["d"], res["W1"], res["kappa"], res["typeC"], res["typeF"], tarch, pred, ok))
print("\nITEM4 TABLE  m q p | d W1 kappa | type_colon type_Fp type_archived | pred=d-W1+kappa-1 | ok")
for t in table:
    print("  %3d %4d %2d | %2d %2d %2d | %2d %2d %s | %2d | %s" % (t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], "OK" if t[10] else "FAIL"))
print("ITEM4 TOTAL e=2 cases=%d fails=%d" % (n, fails))
