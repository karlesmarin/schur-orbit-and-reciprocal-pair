# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# item2_log.py -- finite logarithm: prod_{j=-m0}^{m0} (T - xi^j) mod (p, pi^{s+1}) in F_p[omega][pi]/(pi^{s+1}),
#                 xi = omega*rho, pi = rho-1, rho of order p^v (v > a), omega of order n; nP = 2 m0 + 1, P = p^a, s = phi(P).
# Claim: = T^{nP} - 1 + pi^s * L_p(T^{n p^{a-1}}),  L_p(X) = sum_{r=1}^{p-1} X^r / r.
# Run: python item2_log.py
import numpy as np
from math import comb

def cyclo_mod(n, p):
    # Phi_n over Z by division, then mod p ; coefficients low->high
    def polydiv(a, b):
        a = a[:]; out = [0] * (len(a) - len(b) + 1)
        for i in range(len(a) - len(b), -1, -1):
            c = a[i + len(b) - 1] // b[-1]
            out[i] = c
            for k in range(len(b)):
                a[i + k] -= c * b[k]
        assert all(x == 0 for x in a)
        return out
    phis = {}
    def phi_poly(k):
        if k in phis: return phis[k]
        num = [-1] + [0] * (k - 1) + [1]
        for dd in range(1, k):
            if k % dd == 0:
                num = polydiv(num, phi_poly(dd))
        phis[k] = num
        return num
    return [c % p for c in phi_poly(n)]

def reduce_w(vec, phin, p):
    v = [int(x) % p for x in vec]
    dn = len(phin) - 1
    for i in range(len(v) - 1, dn - 1, -1):
        c = v[i]
        if c:
            for k in range(dn + 1):
                v[i - dn + k] = (v[i - dn + k] - c * phin[k]) % p
    return v[:dn]

def run(p, a, n, v, ctrl=None):
    P = p**a; s = (p - 1) * p**(a - 1); m0 = (n * P - 1) // 2
    E = (p - 1) * p**(v - 1)
    assert s + 1 <= E
    S1 = s + 1; deg = n * P
    X = np.zeros((deg + 1, S1, n), dtype=np.int64)
    X[0, 0, 0] = 1
    cur = 0
    order = p**v
    rng = range(0, n * P) if ctrl == "shift" else range(-m0, m0 + 1)
    for j in rng:
        jj = j % order
        b = np.array([comb(jj, u) % p for u in range(S1)], dtype=np.int64)
        # c * X[:cur+1]
        Y = np.roll(X[:cur + 1], j % n, axis=2)
        CY = np.zeros_like(Y)
        for u in range(S1):
            if b[u]:
                CY[:, u:, :] = (CY[:, u:, :] + b[u] * Y[:, :S1 - u, :]) % p
        new = np.zeros_like(X)
        new[1:cur + 2] = X[:cur + 1]
        new[:cur + 1] = (new[:cur + 1] - CY) % p
        X = new; cur += 1
    phin = cyclo_mod(n, p)
    dn = len(phin) - 1
    # expected
    exp = np.zeros((deg + 1, S1, dn), dtype=np.int64)
    exp[deg, 0, 0] = 1
    exp[0, 0, 0] = (-1) % p
    for r in range(1, p):
        exp[r * n * p**(a - 1), s, 0] = (exp[r * n * p**(a - 1), s, 0] + (-1 if ctrl == "neg" else 1) * pow(r, -1, p)) % p
    bad = 0; omega_dep = 0
    got_s = {}
    for k in range(deg + 1):
        for t in range(S1):
            red = reduce_w(X[k, t], phin, p)
            if any(red[1:]):
                omega_dep += 1
            if t == s and any(red):
                got_s[k] = red
            if list(red) != list(exp[k, t]):
                bad += 1
    print("[ctrl=%s] " % ctrl + "p=%d a=%d n=%d v=%d  s=%d nP=%d : mismatched coefficients=%d, omega-dependent coefficients=%d ; pi^s part (T-degree: [const, w, ...]) = %s" % (
        p, a, n, v, s, deg, bad, omega_dep, {k: w[0] if not any(w[1:]) else w for k, w in sorted(got_s.items())}))

for (p, a, n) in [(3, 1, 5), (3, 2, 5), (5, 1, 7), (3, 1, 7), (5, 2, 7), (7, 1, 5), (3, 3, 5), (3, 1, 11)]:
    for v in (a + 1, a + 2):
        run(p, a, n, v)
print("controls (must show mismatches):")
for (p, a, n) in [(3, 1, 5), (5, 1, 7), (3, 2, 5)]:
    run(p, a, n, a + 1, "neg")
    run(p, a, n, a + 1, "shift")
