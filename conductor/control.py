# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Controles cercanos: (1) segmento no centrado 0..2m (anillo no real): cabe pi^P en B?  (2) V0 con pesos r en vez de r-gorro.
import sys
import numpy as np
from check_descenso import *

def gens_shift(R, m, q):
    p = R.p
    coef = [R.one()]
    for x in range(0, 2*m+1):
        Mx = R.mat(R.ypow(x, q))
        new = [np.zeros(R.D, dtype=np.int64) for _ in range(len(coef)+1)]
        for i, c in enumerate(coef):
            new[i+1] = (new[i+1] + c) % p
            new[i] = (new[i] - Mx @ c) % p
        coef = new
    return coef[:-1]

def control(p, a, n, M):
    P = p**a; m = (P*M*n-1)//2; v = a+1
    while p**v*n <= 2*m:
        v += 1
    q = p**v*n; R = Ring(p, n, P+1)
    erho = crt(0, n, 1, p**v); eom = crt(1, n, 0, p**v)
    one = R.one(); Mpi = R.mat((R.ypow(erho, q) - one) % p)
    piP = one
    for _ in range(P):
        piP = (Mpi @ piP) % p
    Bs = algebra(R, gens_shift(R, m, q))
    Ws = layers(Bs, p, n, P+1)
    B = algebra(R, gens_A(R, m, q))
    om = R.ypow(eom, q); Mom = R.mat(om); ompow = [one]
    for _ in range(n-1):
        ompow.append((Mom @ ompow[-1]) % p)
    MpiP = R.mat(piP)
    Gbad = []
    for j in range(n):
        g = np.zeros(R.D, dtype=np.int64)
        for r in range(n):
            g = (g + (M*r) * ompow[(P*r*j) % n]) % p
        Gbad.append((MpiP @ g) % p)
    print("p=%d a=%d n=%d M=%d | no centrado: W=%s piP in B_shift: %s | V0 con pesos r (mal): in B: %s" % (p, a, n, M, Ws, rank(Bs + [piP], p) == len(Bs), rank(B + Gbad, p) == len(B)), flush=True)

for s in sys.argv[1:]:
    control(*[int(x) for x in s.split(",")])
