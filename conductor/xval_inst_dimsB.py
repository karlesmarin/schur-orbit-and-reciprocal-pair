# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# xval_inst_dimsB.py -- recompute a subset of the layer profiles of descent/dos/l2r outputs with the project's dims_B (perfil_W.py).
# What: independent instrument (numpy, Dickson recursion on y=zeta_q) vs inst.sage (FLINT product tree).
# Run: python xval_inst_dimsB.py   (imports perfil_W.py from ../analisis, or from this directory if it travels here)
import contextlib
import glob
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.join(HERE, "..", "analisis")
if not os.path.exists(os.path.join(AN, "perfil_W.py")):
    AN = HERE
sys.path.insert(0, AN)
os.chdir(AN)
sys.argv = ["x"]
with contextlib.redirect_stdout(io.StringIO()):
    import perfil_W as PW
from math import gcd


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def W_of(q, p, m, n, K):
    dB = PW.dims_B(q, p, m, n, K)
    return [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, len(dB))]


MMAX, DMAX = 70, 130
tot = bad = 0
for f in sorted(glob.glob(os.path.join(HERE, "descent_OUT_all.txt"))):
    for ln in open(f):
        mo = re.match(r"p=(\d+) a=(\d+) n=\s*(\d+) M=\s*(\d+) v=(\d+) q=(\d+) m=(\d+) d=\d+ \| W=(\[[^\]]*\]).*W1\(A_(\d+)\((\d+)\)\)=(\d+)", ln)
        if not mo:
            continue
        p, a, n, M, v, q, m = map(int, mo.groups()[:7])
        W = eval(mo.group(8)); m0, q0, w10 = int(mo.group(9)), int(mo.group(10)), int(mo.group(11))
        K = p ** a + 1
        if m > MMAX or K * phi(n) > DMAX:
            continue
        Wd = W_of(q, p, m, n, K); W0 = W_of(q0, p, m0, n, 2)
        ok = (Wd == W) and (W0[1] == w10)
        tot += 1; bad += (not ok)
        print("descent q=%d p=%d m=%d K=%d dims_B W=%s inst W=%s | A0 W1 %d vs %d %s" % (q, p, m, K, Wd, W, W0[1], w10, "OK" if ok else "MISMATCH"))
n = None
for ln in open(os.path.join(HERE, "dos_OUT.txt")):
    mo = re.match(r"n=(\d+) ", ln)
    if mo:
        n = int(mo.group(1)); continue
    mo = re.match(r"\s+p=\s*(\d+) N=(\d+) v=(\d+) q=(\d+) m=(\d+) W=(\[[^\]]*\])", ln)
    if not mo:
        continue
    p, N, v, q, m = map(int, mo.groups()[:5]); W = eval(mo.group(6))
    if m > MMAX or 3 * phi(n) > DMAX:
        continue
    Wd = W_of(q, p, m, n, 3)
    ok = (Wd == W); tot += 1; bad += (not ok)
    print("dos n=%d q=%d p=%d m=%d dims_B W=%s inst W=%s %s" % (n, q, p, m, Wd, W, "OK" if ok else "MISMATCH"))
for ln in open(os.path.join(HERE, "l2r_OUT.txt")):
    mo = re.match(r"r=\s*(\d+) q'=(\d+) p=(\d+) \S+\s+N/B=\d+ m=(\d+) v=(\d+) q=(\d+) \| W=(\[[^\]]*\])", ln)
    if not mo:
        continue
    r, qp, p, m, v, q = map(int, mo.groups()[:6]); W = eval(mo.group(7))
    if m > MMAX or 2 * phi(qp) > DMAX:
        continue
    Wd = W_of(q, p, m, qp, 2)
    ok = (Wd == W); tot += 1; bad += (not ok)
    print("l2r q'=%d q=%d p=%d m=%d dims_B W=%s inst W=%s %s" % (qp, q, p, m, Wd, W, "OK" if ok else "MISMATCH"))
print("XVAL dims_B TOTAL=%d mismatches=%d" % (tot, bad))
