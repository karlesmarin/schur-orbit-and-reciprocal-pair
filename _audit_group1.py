# -*- coding: utf-8 -*-
"""Re-runs every ancillary script that backs the FIRST group of the paper's verification table and
saves the full output of each, so each row can be checked against what the script actually prints
rather than against a total aggregated by hand.
Carles Marin + Claude (AI assistant)."""
import io, os, subprocess, sys

SCRIPTS = ["theorem_full.py", "law_control.py", "falsify.py", "d_from_quotient.py",
           "sign_ayyer_idiom.py", "single_char.py", "extra_locus.py", "extra_structure.py",
           "rect_degeneracy.py", "enumeration.py", "ak53_consistency.py"]

os.makedirs("_out", exist_ok=True)
for s in SCRIPTS:
    if not os.path.exists(s):
        print("MISSING: %s" % s); continue
    print("running %s" % s); sys.stdout.flush()
    try:
        r = subprocess.run([sys.executable, s], capture_output=True, text=True, timeout=1800)
        out = r.stdout + ("\n--- STDERR ---\n" + r.stderr if r.stderr.strip() else "")
        io.open("_out/%s.txt" % s[:-3], "w", encoding="utf-8").write(out)
        print("   exit=%d lines=%d" % (r.returncode, out.count("\n")))
    except subprocess.TimeoutExpired:
        print("   TIMEOUT")
