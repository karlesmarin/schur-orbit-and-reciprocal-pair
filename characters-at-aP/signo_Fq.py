# -*- coding: utf-8 -*-
# VERIFICACION DEL INFORME: LA FORMULA CON SIGNO.   8 de septiembre de 2026.
#
# El informe propone, para tipo C y TODA potencia:
#     chi = 0                                     si  N_q(lambda) > r_q
#     chi = (-1)^{E_{p,q}(lambda)} Z_q(ell)/Z_q(rho)   si  N_q(lambda) = r_q
# con  Z_q(x) = prod_{alpha>0, q | (x,alpha)} (x,alpha)   y
#      E_{p,q}(lambda) = sum_{alpha>0} ( floor(p(ell,a)/q) - floor(p(rho,a)/q) ),   p = k/d.
#
# ESTO CHOCA APARENTEMENTE con sign_k2.py, que declaro que el signo NO es una paridad de quince
# estadisticos naturales.  No es contradiccion: E no estaba entre los quince --- lleva PISOS, no
# conteos.  Pero por eso mismo hay que medirlo, no creerlo.
#
# ADEMAS se comprueban aqui tres objeciones concretas del informe:
#   O1  el contraejemplo al enunciado literal del Teorema del valor absoluto en k=2:
#       C_2, q=3, lambda=(4,2) -> chi = 0 pero la formula por casos da 1.
#   O2  "la clase que queda vacia SI esta determinada por lambda": la unica cuya ocupacion esta una
#       unidad por debajo de su capacidad.
#   O3  el recuento cerrado en regimen regular:  S_{C_m}(q) = prod (q-(2i-1)) si q impar,
#       prod (q-2i) si q par --- y por tanto S(q) = S(q-1) para q par, en TODO rango.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python signo_Fq.py > signo_Fq_OUT.txt 2>&1

import sys
from itertools import product as iproduct
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 100


def dominants(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for rest in rec(k - 1, v):
                yield (v,) + rest
    return rec(m, top)


def h_closed(q, d, N):
    base = [0] * (N + 1)
    k = 0
    while q * k <= N:
        c = 1
        for i in range(1, d):
            c = c * (k + i) // i
        base[q * k] = c
        k += 1
    h = [0] * (N + 1)
    for n in range(N + 1):
        if d % 2 == 0:
            v = base[n] - (2 * base[n - 1] if n >= 1 else 0) + (base[n - 2] if n >= 2 else 0)
        else:
            v = base[n] - (base[n - 2] if n >= 2 else 0)
        h[n] = v
    return h


def det_int(M):
    n = len(M)
    A = [row[:] for row in M]
    sign, prev = 1, 1
    for c in range(n):
        if A[c][c] == 0:
            piv = next((r for r in range(c + 1, n) if A[r][c] != 0), None)
            if piv is None:
                return 0
            A[c], A[piv] = A[piv], A[c]
            sign = -sign
        for r in range(c + 1, n):
            for kk in range(c + 1, n):
                A[r][kk] = (A[r][kk] * A[c][c] - A[r][c] * A[c][kk]) // prev
            A[r][c] = 0
        prev = A[c][c]
    return sign * A[n - 1][n - 1]


def chi(lam, m, h):
    H = lambda n: (h[n] if 0 <= n < len(h) else 0)
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    return det_int(M) // 2


def pares_C(v):
    m = len(v)
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            out.append(v[i] - v[j])
            out.append(v[i] + v[j])
    for i in range(m):
        out.append(2 * v[i])
    return out


def cls(x, q):
    return min(x % q, (-x) % q)


def prod(v):
    r = 1
    for x in v:
        r *= x
    return r


# --------------------------------------------------------------------------- la formula
print(SEP)
print("LA FORMULA DEL INFORME:  chi = (-1)^E * Z_q(ell)/Z_q(rho)   sobre los supervivientes")
print(SEP)
print(f"{'m':>3} {'k':>3} {'d':>3} {'q':>4} {'p':>3} {'vivos':>7} {'formula':>12} "
       f"{'|.| solo':>10} {'signo':>10}")

T = {"ok": 0, "bad": 0, "mod": 0, "modn": 0, "sig": 0, "sign": 0}
fallos = []
for m in range(2, 9):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    Prho = pares_C(rho)
    top = 9 if m <= 4 else (7 if m <= 6 else 5)
    for k in range(1, t):
        d = gcd(k, t)
        q = t // d
        p = k // d
        if q < 2:
            continue
        h = h_closed(q, d, 3 * (top + m) + 6)
        r_q = sum(1 for b in Prho if b % q == 0)
        Zrho = prod([b for b in Prho if b % q == 0])
        Erho = sum((p * b) // q for b in Prho)
        vivos = ok = bad = mod = modn = sig = sign_ = 0
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            P = pares_C(ell)
            if sum(1 for a in P if a % q == 0) != r_q:
                continue
            vivos += 1
            v = chi(lam, m, h)
            Zell = prod([a for a in P if a % q == 0])
            E = sum((p * a) // q for a in P) - Erho
            pred = (-1) ** (E % 2) * (Zell // Zrho)
            ok, bad = (ok + 1, bad) if pred == v else (ok, bad + 1)
            if pred != v and len(fallos) < 6:
                fallos.append((m, k, d, q, lam, v, pred))
            mod, modn = (mod + 1, modn) if abs(pred) == abs(v) else (mod, modn + 1)
            s_pred = 1 if (E % 2 == 0) else -1
            s_real = 1 if v > 0 else -1
            sig, sign_ = (sig + 1, sign_) if s_pred == s_real else (sig, sign_ + 1)
        if vivos == 0:
            continue
        T["ok"] += ok; T["bad"] += bad
        T["mod"] += mod; T["modn"] += modn
        T["sig"] += sig; T["sign"] += sign_
        print(f"{m:>3} {k:>3} {d:>3} {q:>4} {p:>3} {vivos:>7} {ok:>5}/{vivos:<6} "
              f"{mod:>4}/{vivos:<5} {sig:>4}/{vivos:<5}"
              + ("" if bad == 0 else "   <-- FALLA"))

print(SEP)
print(f"   FORMULA COMPLETA (valor con signo): {T['ok']} aciertos, {T['bad']} fallos")
print(f"   solo el modulo ..................: {T['mod']} aciertos, {T['modn']} fallos")
print(f"   solo el signo (-1)^E ............: {T['sig']} aciertos, {T['sign']} fallos")
if fallos:
    print("   primeros fallos (m,k,d,q,lambda,chi,pred):")
    for f in fallos:
        print("    ", f)

# --------------------------------------------------------------------------- O1
print(SEP)
print("O1 -- contraejemplo del informe al enunciado literal del teorema de k=2")
m, k, q, d = 2, 2, 3, 2
rho = [2, 1]
lam = (4, 2)
ell = [lam[i] + rho[i] for i in range(m)]
h = h_closed(q, d, 40)
v = chi(lam, m, h)
P = pares_C(ell)
r_q = sum(1 for b in pares_C(rho) if b % q == 0)
Nq = sum(1 for a in P if a % q == 0)
fijas = {0} | ({q // 2} if q % 2 == 0 else set())
porclase = {}
for e in ell:
    porclase.setdefault(cls(e, q), []).append(e)
D = 1
for c, f in porclase.items():
    if c in fijas and len(f) == 1:
        D *= (2 * f[0]) // q
    elif c not in fijas and len(f) == 2:
        a, b = f
        D *= (abs(a - b) // q) if (a - b) % q == 0 else ((a + b) // q)
print(f"   lambda={lam}  ell={ell}  (ell,alpha) = {sorted(P)}")
print(f"   N_q={Nq}  r_q={r_q}   ->  {'NO superviviente' if Nq != r_q else 'superviviente'}")
print(f"   chi calculado = {v}   ;   el producto por casos del teorema da {D}")
print("   " + ("OK  el informe TIENE RAZON: el enunciado literal necesita la hipotesis"
               if (v == 0 and D != 0) else "el contraejemplo NO se reproduce"))

# --------------------------------------------------------------------------- O2
print(SEP)
print("O2 -- '¿esta determinada por lambda la clase que queda vacia?'")
mal = tot = 0
for m in range(2, 8):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    Prho = pares_C(rho)
    for k in range(1, t):
        d = gcd(k, t)
        q = t // d
        if q < 2 or d % 2:
            continue
        r_q = sum(1 for b in Prho if b % q == 0)
        fijas = {0} | ({q // 2} if q % 2 == 0 else set())
        top = 8 if m <= 4 else 6
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            if sum(1 for a in pares_C(ell) if a % q == 0) != r_q:
                continue
            occ = {}
            for e in ell:
                occ[cls(e, q)] = occ.get(cls(e, q), 0) + 1
            faltan = [c for c in range(q // 2 + 1)
                      if occ.get(c, 0) < ((d // 2) if c in fijas else d)]
            tot += 1
            if len(faltan) != 1:
                mal += 1
print(f"   supervivientes con d par: {tot}; en {tot - mal} hay EXACTAMENTE una clase por debajo de")
print(f"   su capacidad, y en {mal} no.")
print("   " + ("OK  el informe tiene razon: la plaza vacia esta determinada"
               if mal == 0 else "el informe NO acierta aqui: hay casos con varias"))

# --------------------------------------------------------------------------- O3
print(SEP)
print("O3 -- recuento cerrado en regimen regular:  prod(q-(2i-1)) impar / prod(q-2i) par")
print(f"   {'m':>3} {'q':>4} {'S medido':>10} {'formula':>10} {'=?':>4}")
malo3 = 0
for m in range(1, 5):
    for q in range(2 * m + 1, 2 * m + 9):
        S = 0
        for x in iproduct(range(q), repeat=m):
            # x en coordenadas de pesos fundamentales -> ell en Z^m: ell_i = sum_{j>=i} x_j
            ell = [sum(x[j] for j in range(i, m)) for i in range(m)]
            P = pares_C(ell)
            rho = [m - i for i in range(m)]
            r_q = sum(1 for b in pares_C(rho) if b % q == 0)
            if sum(1 for a in P if a % q == 0) == r_q:
                S += 1
        f = (prod([q - (2 * i - 1) for i in range(1, m + 1)]) if q % 2
             else prod([q - 2 * i for i in range(1, m + 1)]))
        marca = "SI" if S == f else "no"
        if S != f:
            malo3 += 1
        print(f"   {m:>3} {q:>4} {S:>10} {f:>10} {marca:>4}")
print("   " + ("OK  el informe acierta en regimen regular" if malo3 == 0
               else f"   {malo3} discrepancias"))

print(SEP)
if T["bad"] == 0:
    print("VEREDICTO: la formula con signo del informe SE CUMPLE en todo lo medido.")
else:
    print("VEREDICTO: la formula con signo NO se cumple; leer los fallos.")
