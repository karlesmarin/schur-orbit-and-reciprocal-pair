# -*- coding: utf-8 -*-
"""Checks the claim  |nu| >= N+1 = 2r+3  for every non-standard label nu, against the labels actually
printed in the archived outputs (typeD_rule.txt, typeD_residue.txt, prove_W.txt).  Also checks the
consequence: every residue shape lambda has |lambda| >= 2r+3.
Authors: Carles Marin, Claude (AI assistant)."""
import io, os, re

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
TUP = re.compile(r"\((\s*\d+\s*(?:,\s*\d+\s*)*,?)\)")


def tuples_in(line):
    out = []
    for m in TUP.finditer(line):
        parts = [p for p in m.group(1).replace(" ", "").split(",") if p != ""]
        if parts:
            out.append(tuple(int(p) for p in parts))
    return out


def scan(fname, r_of_header):
    """walk a saved output, tracking the current r from its headers, and collect labels nu."""
    path = os.path.join(OUT, fname)
    if not os.path.exists(path):
        print("  MISSING %s" % fname); return {}
    per_r, r = {}, None
    for line in io.open(path, encoding="utf-8", errors="replace"):
        m = r_of_header.search(line)
        if m:
            r = int(m.group(1)); continue
        if "nu=" not in line:
            continue
        for t in tuples_in(line.split("nu=", 1)[1]):
            if r is not None and len(t) > 1:          # skip (c1,c2) pairs by taking the first tuple
                per_r.setdefault(r, []).append(t)
                break
    return per_r


print("labels nu printed in the archived outputs, smallest |nu| per r:")
hdr = re.compile(r"^\s*r\s*=\s*(\d+)")
bad = 0
for fname in ["typeD_rule.txt", "typeD_residue.txt"]:
    per_r = scan(fname, hdr)
    for r in sorted(per_r):
        sizes = [sum(t) for t in per_r[r]]
        thr = 2 * r + 3
        mn = min(sizes)
        who = per_r[r][sizes.index(mn)]
        ok = "ok" if mn >= thr else "VIOLATION"
        if mn < thr:
            bad += 1
        print("  %-22s r=%d : %3d labels, min |nu| = %2d (threshold 2r+3 = %2d)  %-9s  smallest: %s"
              % (fname, r, len(sizes), mn, thr, ok, str(who)))

# the consequence, on the residue shapes prove_W prints
print("\nresidue shapes lambda (no proved isolating witness), smallest |lambda| per r:")
path = os.path.join(OUT, "prove_W.txt")
r = None
for line in io.open(path, encoding="utf-8", errors="replace"):
    m = hdr.search(line)
    if m:
        r = int(m.group(1))
    if "residue examples" not in line:
        continue
    shapes = [t for t in tuples_in(line.split("residue examples", 1)[1]) if len(t) > 1]
    # entries are (lambda, ell) pairs flattened by the printer; keep tuples of length >= 2
    sizes = [sum(t) for t in shapes]
    thr = 2 * r + 3
    if sizes:
        mn = min(sizes)
        ok = "ok" if mn >= thr else "VIOLATION"
        if mn < thr:
            bad += 1
        print("  r=%d : %d shapes shown, min |lambda| = %d (threshold %d)  %s"
              % (r, len(shapes), mn, thr, ok))

print("\nVIOLATIONS: %d" % bad)
