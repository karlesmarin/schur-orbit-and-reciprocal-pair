# -*- coding: utf-8 -*-
# LAS POTENCIAS DEL SEGUNDO ELEMENTO DE KOSTANT.   25 de agosto de 2026.
#
# POR QUE.  Kostant (1976, §3) distingue dos elementos regulares, a_Q (Coxeter, orden h) y a_P
# (orden r*h^v; en C_m vale h+2 = 2m+2).  Para las POTENCIAS de a_Q existe teoria: NPP
# (arXiv:2504.14684) introducen el sistema de raices R(m) = {raices de ALTURA divisible por m} y
# una constante d_m que es la dimension de una representacion del grupo dual; Polo
# (arXiv:2504.09204) completa su clasificacion en todos los casos.  Para las potencias de a_P no
# hay nada escrito que hayamos encontrado.  Esta gate mide si la hay.
#
# LOS OBJETOS.  En Sp(2m):  t = 2m+2,  a_P = rho(zeta_t) con rho = (m, m-1, ..., 1) --- salvo Weyl,
# el vector escalera a = (1, 2, ..., m).  Para la potencia k-esima,  d = gcd(k,t),  q = t/d es su
# orden, y
#         alpha(a_P^k) = 1   <=>   q | <alpha, rho>.
# Llamamos R_P(q) = { alpha en Phi(C_m) : q | <alpha, rho> } al sistema del centralizador.
# Notese que NPP usan la ALTURA ht(beta) = <beta, rho^v>; nosotros <alpha, rho>.  Es el mismo eje
# raiz/corraiz que separa a_P de a_Q, y por eso los dos sistemas no tienen por que coincidir.
#
# QUE SE MIDE, y en este orden --- los controles primero.
#   C0  QUE VARIANTE de la formula de Jacobi-Trudi simplectica es la correcta, contra el
#       bialternante en un punto REGULAR (donde el bialternante si vale).  La variante mala se
#       reporta al lado: si las dos pasaran, el instrumento no mediria nada.
#   C1  El alfabeto de a_P^2 colapsa:  {x_j^{+-1}} = cada raiz (m+1)-esima NO trivial, DOS veces,
#       luego  H(z) = ((1-z)/(1-z^{m+1}))^2  y los h_n son ENTEROS.  Se comprueba contra la serie
#       calculada a pelo desde los autovalores.
#   V1  REGLA DE CAPACIDAD.  chi_lambda(a_P^k) != 0  <=>  cada clase plegada mod q recibe a lo sumo
#       d filas si no es fija (c != 0, q/2), y a lo sumo floor(d/2) si lo es.  Para k=1 esto es
#       (1,0), que es exactamente el Teorema 4.1 de la nota (clase fija -> 0, colision -> 0).
#   V2  EL VALOR, para k=2:  |chi| = producto de una dimension por bloque, que es la formula de
#       dimensiones de Weyl del centralizador  prod |<w.ell, alpha>| / <rho_P, alpha>.
#       CONTROL: sin el factor de las clases fijas la formula debe EMPEORAR en todos los rangos.
#   V3  LA DUALIDAD.  ¿Es R_P(q) en C_m, leido por corraices, el R(q) de NPP en el dual B_m?
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python aP_potencias.py > aP_potencias_OUT.txt 2>&1

import cmath
import itertools
import random
import sys
from collections import Counter, defaultdict
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SEP = "=" * 92


# --------------------------------------------------------------------------- utilidades

def dominants(m, top):
    for lam in itertools.product(range(top + 1), repeat=m):
        if all(lam[i] >= lam[i + 1] for i in range(m - 1)):
            yield list(lam)


def det_generic(M, zero, is_small):
    """Eliminacion gaussiana valida para float y para complex."""
    n = len(M)
    M = [row[:] for row in M]
    d = type(M[0][0])(1)
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if is_small(M[p][c]):
            return zero
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        inv = 1 / M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] * inv
            if f:
                for kk in range(c, n):
                    M[r][kk] -= f * M[c][kk]
    return d


def det_f(M):
    return det_generic(M, 0.0, lambda x: abs(x) < 1e-12)


def det_c(M):
    return det_generic(M, 0j, lambda x: abs(x) < 1e-11)


def det_int(M):
    """Bareiss: determinante exacto de matriz entera."""
    n = len(M)
    if n == 0:
        return 1
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for c in range(n - 1):
        if M[c][c] == 0:
            piv = next((r for r in range(c + 1, n) if M[r][c] != 0), None)
            if piv is None:
                return 0
            M[c], M[piv] = M[piv], M[c]
            sign = -sign
        for r in range(c + 1, n):
            for kk in range(c + 1, n):
                M[r][kk] = (M[r][kk] * M[c][c] - M[r][c] * M[c][kk]) // prev
            M[r][c] = 0
        prev = M[c][c]
    return sign * M[n - 1][n - 1]


def h_from_eigs(xs, N, one):
    """h_n, n=0..N, del alfabeto {x, 1/x : x in xs}, por serie de potencias."""
    h = [one * 0] * (N + 1)
    h[0] = one
    for x in xs:
        for y in (x, 1 / x):
            acc = one * 0
            new = [one * 0] * (N + 1)
            for n in range(N + 1):
                acc = acc * y + h[n]
                new[n] = acc
            h = new
    return h


def jt_matrix(lam, H, m, variant):
    """Weyl / Koike-Terada para Sp(2m).  Devuelve (matriz, divisor)."""
    if variant == "A":
        return [[H(lam[i] - (i + 1) + (j + 1)) - H(lam[i] - (i + 1) - (j + 1))
                 for j in range(m)] for i in range(m)], 1
    return [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
             for j in range(m)] for i in range(m)], 2


def sp_bialt(lam, xs):
    m = len(lam)
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    num = [[x ** e - x ** (-e) for x in xs] for e in ell]
    den = [[x ** e - x ** (-e) for x in xs] for e in rho]
    return det_f(num) / det_f(den)


def cls(x, q):
    return min(x % q, (-x) % q)


# --------------------------------------------------------------------------- C0

print(SEP)
print("C0 -- QUE VARIANTE DE JACOBI-TRUDI SIMPLECTICA ES LA CORRECTA (punto regular, bialternante)")
print("     A: det(h_{l_i-i+j} - h_{l_i-i-j})        B: (1/2) det(h_{l_i-i+j} + h_{l_i-i-j+2})")
random.seed(7)
for m in (2, 3, 4):
    xs = [random.uniform(1.3, 2.7) for _ in range(m)]
    h = h_from_eigs(xs, 60, 1.0)
    Hf = lambda n: (h[n] if 0 <= n < len(h) else 0.0)
    okA = okB = tot = 0
    for lam in dominants(m, 4):
        tot += 1
        ref = sp_bialt(lam, xs)
        for v in ("A", "B"):
            M, div = jt_matrix(lam, Hf, m, v)
            if abs(det_f(M) / div - ref) <= 1e-6 * max(1.0, abs(ref)):
                okA, okB = (okA + 1, okB) if v == "A" else (okA, okB + 1)
    print(f"     m={m}: {tot} pesos   variante A: {okA}/{tot}   variante B: {okB}/{tot}"
          + ("   <- B es la buena, y A falla: el test puede fallar" if okB == tot and okA < tot else "   <- REVISAR"))

# --------------------------------------------------------------------------- C1

print(SEP)
print("C1 -- EL ALFABETO DE a_P^2 COLAPSA:  H(z) = ((1-z)/(1-z^{m+1}))^2,  h_n ENTEROS")


def h_aP2_closed(m):
    """De la forma cerrada:  (1-z)^2 * sum_j (j+1) z^{(m+1)j}."""
    q = m + 1

    def c(n):
        return 0 if (n < 0 or n % q) else n // q + 1

    return lambda n: 0 if n < 0 else c(n) - 2 * c(n - 1) + c(n - 2)


for m in (2, 3, 4, 5, 6):
    t = 2 * m + 2
    z = cmath.exp(2j * cmath.pi / t)
    xs = [z ** (2 * j) for j in range(1, m + 1)]           # autovalores de a_P^2
    N = 30
    hnum = h_from_eigs(xs, N, 1 + 0j)
    Hcl = h_aP2_closed(m)
    peor = max(abs(hnum[n] - Hcl(n)) for n in range(N + 1))
    enteros = all(abs(hnum[n].imag) < 1e-8 for n in range(N + 1))
    print(f"     m={m}: max|h_serie - h_forma_cerrada| = {peor:.2e}   h_n reales: {enteros}"
          f"   h_0..h_8 = {[Hcl(n) for n in range(9)]}")

# --------------------------------------------------------------------------- V1

print(SEP)
print("V1 -- REGLA DE CAPACIDAD.   d = gcd(k,t),  q = t/d = orden de a_P^k")
print("     chi != 0  <=>  cada clase plegada mod q recibe <= d filas (<= floor(d/2) si es fija)")
print("     k=1 da (1,0), que es el Teorema 4.1 de la nota: clase fija -> 0, colision -> 0")
print()


def h_power_num(m, k, N):
    t = 2 * m + 2
    z = cmath.exp(2j * cmath.pi / t)
    return h_from_eigs([z ** ((j * k) % t) for j in range(1, m + 1)], N, 1 + 0j)


def chi_num(lam, m, h):
    H = lambda n: (h[n] if 0 <= n < len(h) else 0j)
    M, div = jt_matrix(lam, H, m, "B")
    v = det_c(M) / div
    if abs(v.imag) > 1e-5 or abs(v.real - round(v.real)) > 1e-5:
        return None
    return int(round(v.real))


print(f"{'m':>3} {'k':>3} {'q':>4} {'cap':>5} {'capF':>5} {'pesos':>7} {'aciertos':>9} {'fallos':>7}")
tot_ok = tot_bad = 0
for m in range(2, 8):
    t = 2 * m + 2
    for k in range(1, 7):
        d = gcd(k, t)
        q = t // d
        if q < 3:
            continue
        capN, capF = d, d // 2
        fixed = {0} | ({q // 2} if q % 2 == 0 else set())
        top = 8 if m <= 4 else 6
        h = h_power_num(m, k, max(20, 2 * (top + m) + 4))
        ok = bad = tot = 0
        for lam in dominants(m, top):
            ell = [lam[i] + (m - i) for i in range(m)]
            v = chi_num(lam, m, h)
            if v is None:
                continue
            cnt = Counter(cls(e, q) for e in ell)
            pred = all(n <= (capF if c in fixed else capN) for c, n in cnt.items())
            tot += 1
            ok, bad = (ok + 1, bad) if pred == (v != 0) else (ok, bad + 1)
        tot_ok += ok
        tot_bad += bad
        print(f"{m:>3} {k:>3} {q:>4} {capN:>5} {capF:>5} {tot:>7} {ok:>9} {bad:>7}"
              + ("" if bad == 0 else "   <-- FALLA"))
print(f"\n     TOTAL: {tot_ok} aciertos, {tot_bad} fallos")

# --------------------------------------------------------------------------- V2

print(SEP)
print("V2 -- EL VALOR EN a_P^2.   |chi| = prod_bloques (dimension de Weyl)")
print("     clase no fija con dos filas -> |ell_i -+ ell_j| / q   (el signo que hace entero)")
print("     clase fija con una fila     -> 2*ell_i / q")
print("     CONTROL: quitando el factor de las clases fijas la formula debe EMPEORAR")
print()


def pred_valor(ell, q, fixed, con_fijas=True):
    por = defaultdict(list)
    for e in ell:
        por[cls(e, q)].append(e)
    for c, v in por.items():
        if len(v) > (1 if c in fixed else 2):
            return 0, 0
    p, nf = 1, 0
    for c, es in por.items():
        if c in fixed:
            if con_fijas:
                p *= 2 * es[0] // q
                nf += 1
        elif len(es) == 2:
            a, b = es
            dd = (a - b) if (a - b) % q == 0 else (a + b)
            p *= abs(dd) // q
            nf += 1
    return p, nf


print(f"{'m':>3} {'q':>4} {'pesos':>7} {'|chi| = pred':>14} {'CONTROL sin fijas':>19} {'#fact = |R_P^+|':>17}")
for m in range(2, 8):
    q = m + 1
    Hcl = h_aP2_closed(m)
    fixed = {0} | ({q // 2} if q % 2 == 0 else set())
    nphi = len([1 for i in range(1, m + 1) for j in range(i + 1, m + 1) if (i + j) % q == 0])
    nphi += 1 if q % 2 == 0 else 0
    ok = ctrl = tot = nfok = 0
    for lam in dominants(m, 8 if m <= 5 else 6):
        ell = [lam[i] + (m - i) for i in range(m)]
        M, div = jt_matrix(lam, lambda n: Hcl(n), m, "B")
        dd = det_int(M)
        assert dd % div == 0, (m, lam)
        val = abs(dd // div)
        tot += 1
        p, nf = pred_valor(ell, q, fixed, True)
        if p == val:
            ok += 1
            if val and nf == nphi:
                nfok += 1
        p2, _ = pred_valor(ell, q, fixed, False)
        ctrl += (p2 == val)
    print(f"{m:>3} {q:>4} {tot:>7} {str(ok)+'/'+str(tot):>14} {str(ctrl)+'/'+str(tot):>19} {nfok:>17}")

# --------------------------------------------------------------------------- V3

print(SEP)
print("V3 -- LA DUALIDAD.   R_P(q) en C_m, leido por corraices,  frente a  R(q) = {ht = 0 mod q} en B_m")


def roots_C(m):
    R = []
    for i in range(m):
        R.append(("long", tuple(2 if k == i else 0 for k in range(m))))
        for j in range(i + 1, m):
            for s in (-1, 1):
                R.append(("short", tuple((1 if k == i else 0) + s * (1 if k == j else 0) for k in range(m))))
    return R


def roots_B(m):
    R = []
    for i in range(m):
        R.append(("short", tuple(1 if k == i else 0 for k in range(m))))
        for j in range(i + 1, m):
            for s in (-1, 1):
                R.append(("long", tuple((1 if k == i else 0) + s * (1 if k == j else 0) for k in range(m))))
    return R


lin = lambda v, m: sum(c * (m - i) for i, c in enumerate(v))   # <., rho_C> = ht_B, misma forma lineal

print("     control: alturas de las raices simples de B_m (deben ser todo 1)")
for m in (3, 5):
    s = [tuple((1 if k == i else 0) - (1 if k == i + 1 else 0) for k in range(m)) for i in range(m - 1)]
    s.append(tuple(1 if k == m - 1 else 0 for k in range(m)))
    print(f"        m={m}: {[lin(x, m) for x in s]}")
print()
print(f"{'m':>3} {'k':>3} {'q':>4} {'|R_P(q)|':>9} {'|R(q) en B|':>12} {'iguales':>8}   diferencia")
for m in range(2, 9):
    t = 2 * m + 2
    for k in (1, 2, 3, 4):
        if t % k:
            continue
        q = t // k
        if q < 2:
            continue
        RP = {(tuple(c // 2 for c in a) if kind == "long" else a)
              for kind, a in roots_C(m) if lin(a, m) % q == 0}
        RB = {b for kind, b in roots_B(m) if lin(b, m) % q == 0}
        extra = sorted(RP - RB)
        falta = sorted(RB - RP)
        tipo = lambda S: Counter("corta" if sum(abs(c) for c in x) == 1 else "larga" for x in S)
        dif = "" if not (extra or falta) else f"solo a_P: {dict(tipo(extra))}   solo Coxeter: {dict(tipo(falta))}"
        print(f"{m:>3} {k:>3} {q:>4} {len(RP):>9} {len(RB):>12} {'SI' if RP == RB else 'NO':>8}   {dif}")
print()
print("     LECTURA: coinciden exactamente cuando q es IMPAR.  Cuando q es par, R_P(q) es")
print("     ESTRICTAMENTE MAYOR, y lo que sobra son siempre corraices e_i --- las raices LARGAS")
print("     2e_i de C_m con q | 2<e_i,rho> pero q no divide <e_i,rho>, o sea <e_i,rho> = q/2.")
print("     Es la clase plegada FIJA q/2: la misma unica raiz larga que en la nota bloquea t=h y")
print("     crea el segundo mecanismo de anulacion.  Un solo objeto, ahora en tres sitios.")
print(SEP)
