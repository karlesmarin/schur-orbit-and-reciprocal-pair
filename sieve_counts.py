# -*- coding: utf-8 -*-
"""The endpoint sieve, quantified: the numbers the paper quotes when it says that testing
Psi_r(lambda)|_{z_i=1} = 0 is a filter and, past one reciprocal pair, nothing more.

Reads fig_data_new.json (the same exact data the Section-8 figures are drawn from) and prints, per r:
the number of shapes in range, how many vanish identically, how many vanish only at the endpoint, and
therefore how many of the sieve's verdicts are spurious.

It also runs the consistency control the counts depend on: a shape that vanishes identically must
also vanish at the endpoint, so the two classes are nested and the sieve never misses a genuine zero.

Authors: Carles Marin, Claude (AI assistant)."""
import io, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(io.open(os.path.join(HERE, "fig_data_new.json"), encoding="utf-8"))["loci"]

print("control: shapes with z_exact but not z_end (must be 0 -- the classes are nested)")
bad_total = 0
for k in ["r1", "r2", "r3"]:
    bad = sum(1 for x in D[k]["rows"] if x["z_exact"] and not x["z_end"])
    bad_total += bad
    print("   %s: %d" % (k, bad))
print("   total: %d" % bad_total)

print("\nthe endpoint sieve, per r")
print("   %-4s %-4s %8s %8s %10s %8s %8s" %
      ("r", "N", "|lam|<=", "shapes", "identical", "endpt", "flagged"))
tot = []
for k in ["r1", "r2", "r3"]:
    b = D[k]
    rows = b["rows"]
    r = (b["N"] - 2) // 2
    mx = max(x["size"] for x in rows)
    ident = sum(1 for x in rows if x["z_exact"])
    endp = sum(1 for x in rows if (not x["z_exact"]) and x["z_end"])
    flag = ident + endp
    tot.append((r, len(rows), ident, endp, flag))
    print("   %-4d %-4d %8d %8d %10d %8d %8d" % (r, b["N"], mx, len(rows), ident, endp, flag))

print("\n   the sieve flags `flagged` shapes and `endpt` of them do not vanish identically:")
for r, n, ident, endp, flag in tot:
    pct = 100.0 * endp / flag if flag else 0.0
    print("      r=%d : flags %3d, spurious %d  ->  %.0f%%" % (r, flag, endp, pct))

print("\n   endpoint-only shapes, named:")
for k in ["r1", "r2", "r3"]:
    b = D[k]
    r = (b["N"] - 2) // 2
    names = [tuple(x["lam"]) for x in b["rows"] if (not x["z_exact"]) and x["z_end"]]
    print("      r=%d : %s" % (r, ", ".join(str(t) for t in names) if names else "none"))
