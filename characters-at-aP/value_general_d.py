# -*- coding: utf-8 -*-
# EL VALOR ABSOLUTO PARA d >= 3.   8 de septiembre de 2026.
#
# POR QUE.  El articulo cierra con un problema: evaluar el limite
#       chi_lambda(a_P^k) = xi^{-(lambda,rho)} lim_{z->xi} prod_{alpha>0} (1-z^{(ell,alpha)})/(1-z^{(rho,alpha)})
# para d >= 3, y dice que ahi "los exponentes que se anulan en el denominador ya no valen todos q
# sino que recorren q, 2q, ..., asi que el recuento ya no es un factor por clase".  Eso es cierto
# como objecion a la FORMA de la formula de d=2, pero no como objecion al limite: si se separan los
# factores que se anulan de los que no, con z = xi(1+eps) y q | a se tiene 1-z^a = -a*eps + O(eps^2),
# luego los eps se cancelan porque hay tantos arriba como abajo (eso ES el criterio) y queda
#
#       CONJETURA:   |chi_lambda(a_P^k)|  =  prod_{alpha>0, q | (ell,alpha)} (ell,alpha)
#                                            ---------------------------------------------
#                                            prod_{alpha>0, q | (rho,alpha)} (rho,alpha)
#
# siempre que la parte que NO se anula tenga modulo 1 --- que en d=2 es la Proposicion de los M_s.
# Esta gate mide si eso vale para TODO d.  En d=2 el denominador vale q^{r_q} y la conjetura SE
# REDUCE al teorema ya publicado, asi que la gate contiene su propio control de coherencia.
#
# EXACTITUD.  El alfabeto de a_P^k es d*mu_q menos dos puntos, luego
#       H(z) = (1-z)^2/(1-z^q)^d   si d es par,     H(z) = (1-z^2)/(1-z^q)^d   si d es impar,
# con coeficientes ENTEROS: el caracter es un determinante entero y no hay ningun flotante en toda
# la gate.  Y el alfabeto es estable por Galois, asi que chi es un entero racional.
#
# CONTROLES Y SENUELOS
#   C0  la serie cerrada H(z) coincide con la calculada a pelo desde los autovalores
#   C1  el criterio: chi = 0 exactamente cuando N_q(lambda) > r_q   (si falla, no se mide nada)
#   D1  SENUELO: la misma formula sin el denominador
#   D2  SENUELO: la misma formula con el estadistico de CORRAICES en vez del de raices
#   D3  SENUELO: el producto sobre TODAS las raices y no solo sobre las que se anulan
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python value_general_d.py > value_general_d_OUT.txt 2>&1

import sys
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
    """H(z) = (1-z)^2/(1-z^q)^d  (d par);  (1-z^2)/(1-z^q)^d  (d impar).  Enteros."""
    # 1/(1-z^q)^d  =  sum_k C(k+d-1, d-1) z^{qk}
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


def h_from_eigs_mp(m, k, t, N):
    """control: la misma serie desde los autovalores, en coma flotante de 40 digitos."""
    from mpmath import mp, mpc, exp, pi
    mp.dps = 40
    z = exp(2j * pi / mpf_t(t))
    xs = [z ** ((j * k) % t) for j in range(1, m + 1)]
    h = [mpc(0)] * (N + 1)
    h[0] = mpc(1)
    for x in xs:
        for y in (x, 1 / x):
            acc = mpc(0)
            new = [mpc(0)] * (N + 1)
            for n in range(N + 1):
                acc = acc * y + h[n]
                new[n] = acc
            h = new
    return h


def mpf_t(x):
    from mpmath import mpf
    return mpf(x)


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
    """los (v,alpha) sobre las raices positivas de C_m, en el estadistico de RAICES."""
    m = len(v)
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            out.append(v[i] - v[j])
            out.append(v[i] + v[j])
    for i in range(m):
        out.append(2 * v[i])
    return out


def coraices_C(v):
    """el estadistico de CORRAICES: <v, alpha^vee>; la raiz larga 2e_i tiene corraiz e_i."""
    m = len(v)
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            out.append(v[i] - v[j])
            out.append(v[i] + v[j])
    for i in range(m):
        out.append(v[i])
    return out


def prod(v):
    r = 1
    for x in v:
        r *= x
    return r


# --------------------------------------------------------------------------- C0
print(SEP)
print("C0 -- la serie cerrada H(z) contra la calculada desde los autovalores (40 digitos)")
peor = 0.0
for m in (3, 4, 5):
    t = 2 * m + 2
    for k in range(1, t):
        d = gcd(k, t)
        q = t // d
        if q < 2:
            continue
        N = 14
        hc = h_closed(q, d, N)
        hv = h_from_eigs_mp(m, k, t, N)
        for n in range(N + 1):
            peor = max(peor, abs(complex(hv[n]) - hc[n]))
print("   peor discrepancia: %.3e" % peor + ("   OK" if peor < 1e-20 else "   FALLA"))
if peor >= 1e-20:
    sys.exit(1)

# --------------------------------------------------------------------------- medida
print(SEP)
print("LA CONJETURA:  |chi| = prod_{q | (ell,a)} (ell,a)  /  prod_{q | (rho,a)} (rho,a)")
print(SEP)
print(f"{'m':>3} {'k':>3} {'d':>3} {'q':>4} {'r_q':>4} {'vivos':>7} {'conjetura':>10} "
      f"{'D1 sin den':>11} {'D2 corr.':>9} {'D3 todas':>9} {'crit':>6}")

tot = {"ok": 0, "bad": 0, "d1": 0, "d1n": 0, "d2": 0, "d2n": 0, "d3": 0, "d3n": 0,
       "crit_ok": 0, "crit_bad": 0, "d3grande": 0}
fallos = []
for m in range(2, 9):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    top = 9 if m <= 4 else (7 if m <= 6 else 5)
    for k in range(1, t):
        d = gcd(k, t)
        q = t // d
        if q < 2 or d < 3:
            continue
        h = h_closed(q, d, 3 * (top + m) + 6)
        den = prod([b for b in pares_C(rho) if b % q == 0])
        den_cor = prod([b for b in coraices_C(rho) if b % q == 0])
        r_q = sum(1 for b in pares_C(rho) if b % q == 0)
        ok = bad = d1 = d1n = d2 = d2n = d3 = d3n = vivos = 0
        cok = cbad = 0
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            P = pares_C(ell)
            Nq = sum(1 for a in P if a % q == 0)
            v = chi(lam, m, h)
            vivo = (Nq == r_q)
            cok, cbad = (cok + 1, cbad) if vivo == (v != 0) else (cok, cbad + 1)
            if not vivo:
                continue
            vivos += 1
            num = prod([a for a in P if a % q == 0])
            pred = num // den if den and num % den == 0 else None
            ok, bad = (ok + 1, bad) if pred == abs(v) else (ok, bad + 1)
            if pred != abs(v) and len(fallos) < 6:
                fallos.append((m, k, d, q, lam, abs(v), pred))
            d1, d1n = (d1 + 1, d1n) if num == abs(v) else (d1, d1n + 1)
            Pc = coraices_C(ell)
            numc = prod([a for a in Pc if a % q == 0])
            pc = numc // den_cor if den_cor and numc % den_cor == 0 else None
            d2, d2n = (d2 + 1, d2n) if pc == abs(v) else (d2, d2n + 1)
            p3 = prod(P) // prod(pares_C(rho))
            d3, d3n = (d3 + 1, d3n) if p3 == abs(v) else (d3, d3n + 1)
        if vivos == 0:
            continue
        tot["ok"] += ok; tot["bad"] += bad
        tot["d1"] += d1; tot["d1n"] += d1n
        tot["d2"] += d2; tot["d2n"] += d2n
        tot["d3"] += d3; tot["d3n"] += d3n
        tot["crit_ok"] += cok; tot["crit_bad"] += cbad
        print(f"{m:>3} {k:>3} {d:>3} {q:>4} {r_q:>4} {vivos:>7} {ok:>5}/{vivos:<4} "
              f"{d1:>5}/{vivos:<5} {d2:>4}/{vivos:<4} {d3:>4}/{vivos:<4} "
              f"{'OK' if cbad == 0 else 'MAL':>6}"
              + ("" if bad == 0 else "   <-- FALLA"))

print(SEP)
print(f"   CONJETURA ........ {tot['ok']} aciertos, {tot['bad']} fallos")
print(f"   criterio (control) {tot['crit_ok']} aciertos, {tot['crit_bad']} fallos")
print(f"   D1 sin denominador {tot['d1']} / {tot['d1'] + tot['d1n']}"
      f"   ({'discrimina' if tot['d1n'] else 'NO DISCRIMINA'})")
print(f"   D2 corraices ..... {tot['d2']} / {tot['d2'] + tot['d2n']}"
      f"   ({'discrimina' if tot['d2n'] else 'NO DISCRIMINA'})")
print(f"   D3 todas las raices {tot['d3']} / {tot['d3'] + tot['d3n']}"
      f"   ({'discrimina' if tot['d3n'] else 'NO DISCRIMINA'})")
if fallos:
    print()
    print("   primeros fallos (m,k,d,q,lambda,|chi|,conjetura):")
    for f in fallos:
        print("    ", f)
print()
if tot["bad"] == 0 and tot["crit_bad"] == 0 and tot["d1n"] and tot["d2n"] and tot["d3n"]:
    print("RESULTADO: la forma cerrada del valor absoluto vale para TODO d medido, no solo d=2,")
    print("           y los tres senuelos fallan, luego la medida discrimina.")
else:
    print("RESULTADO: NO cuadra.  Leer los fallos antes de escribir nada.")
