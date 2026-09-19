# -*- coding: utf-8 -*-
# PUENTE REAL / NO REAL -- el INDICE por reticulos, no por capas.
#
# El puente predice   v_p([O:B_n(q)]) = 3d - 2 + 2*delta,   d = phi(r)/2, delta = d - rango M_r.
# El plano del papel (thm:sombra con e=2) daria  2D - 1 - W_1 = 3d - 2 + delta  con D=2d,
# W_1 = d+1-delta.  O sea: la sombra se separa del plano EXACTAMENTE delta.
#
# Aqui el indice se calcula por CIERRE DE RETICULOS ENTEROS en O = Z[zeta_q], sin usar la teoria
# de capas, para que no sea circular.  Se contrasta tambien el orden real R = A_{r-1}(q).
#
# Authors: Carles Marin, Claude (AI assistant).
import sys

SOLO = None
for a in sys.argv:
    if a.startswith("solo="):
        SOLO = a.split("=", 1)[1]


def indice_reticulo(exps, M, ngen, vueltas=60):
    """indice [O_M : Z[e_1..e_ngen]] con los autovalores zeta_M^e."""
    K = CyclotomicField(M)
    z = K.gen()
    S = PolynomialRing(K, 'T')
    T = S.gen()
    f = prod([T - z ** (e % M) for e in exps])
    co = f.list()
    N = len(exps)
    g = [(-1) ** k * co[N - k] for k in range(1, ngen + 1)]
    g = [y for y in g if y != 0]
    L = matrix(ZZ, [K(1).list()]).hermite_form(include_zero_rows=False)
    for _ in range(vueltas):
        filas = [list(fi) for fi in L.rows()]
        for fi in L.rows():
            el = K(list(fi))
            for gg in g:
                filas.append((el * gg).list())
        L2 = matrix(ZZ, filas).hermite_form(include_zero_rows=False)
        if L2 == L:
            break
        L = L2
    return prod([x for x in L.elementary_divisors() if x != 0])


def indice_real(m, M, vueltas=60):
    """indice [O_+ : A_m(M)] dentro del subcuerpo real, por el mismo cierre."""
    K = CyclotomicField(M)
    z = K.gen()
    Kp, mp = K.subfield(z + z ** -1)
    S = PolynomialRing(K, 'T')
    T = S.gen()
    f = prod([T - z ** j for j in range(1, m + 1)] + [T - z ** (-j) for j in range(1, m + 1)])
    co = f.list()
    N = 2 * m
    g = [(-1) ** k * co[N - k] for k in range(1, m + 1)]
    g = [Kp(y) for y in g if y != 0]
    O = Kp.maximal_order()
    base = O.basis()
    Bmat = matrix(QQ, [Kp(b).list() for b in base]).inverse()
    def coord(x):
        return vector(QQ, Kp(x).list()) * Bmat
    L = matrix(ZZ, [coord(Kp(1))]).hermite_form(include_zero_rows=False)
    dim = Kp.degree()
    for _ in range(vueltas):
        filas = [list(fi) for fi in L.rows()]
        for fi in L.rows():
            el = sum([ZZ(fi[i]) * Kp(base[i]) for i in range(dim)])
            for gg in g:
                filas.append(list(coord(el * gg)))
        L2 = matrix(ZZ, filas).hermite_form(include_zero_rows=False)
        if L2 == L:
            break
        L = L2
    return prod([x for x in L.elementary_divisors() if x != 0])


def rango_Mr(r, p):
    s = [0] * r
    for c in range(1, r):
        s[c] = 2 * c - r
    units = [a for a in range(1, r) if gcd(a, r) == 1]
    return matrix(GF(p), [[s[(inverse_mod(a, r) * c) % r] for c in range(r)] for a in units]).rank()


# (p, v, r, n)  --  n = N r o N r + 1
CELDAS = [(3, 2, 5, 5), (3, 2, 5, 6), (3, 2, 5, 10), (3, 2, 7, 7), (3, 2, 7, 8),
          (5, 1, 7, 7), (5, 1, 7, 8), (3, 2, 11, 11), (3, 2, 11, 12), (7, 1, 9, 9),
          (3, 2, 13, 13), (5, 1, 11, 11), (3, 2, 23, 23)]
if SOLO == "grande":
    CELDAS = [(3, 2, 23, 23)]
elif SOLO == "chico":
    CELDAS = [c for c in CELDAS if c[2] < 23]

print("%-5s %-3s %-4s %-5s %-3s %-4s %-7s %-9s %-9s %-9s %s" % (
    "q", "p", "r", "n", "d", "delt", "v_p(B)", "3d-2+2del", "plano", "v_p(R)", "veredicto"))
ok = mal = 0
for (p, v, r, n) in CELDAS:
    q = p ** v * r
    if not (2 <= n < q):
        continue
    d = euler_phi(r) // 2
    delta = d - rango_Mr(r, p)
    iB = indice_reticulo(list(range(n)), q, n)
    vB = valuation(iB, p)
    pred = 3 * d - 2 + 2 * delta
    plano = 3 * d - 2 + delta
    try:
        iR = indice_real(r - 1, q)
        vR = valuation(iR, p)
    except Exception as exc:
        vR = -1
    bien = (vB == pred) and (vR < 0 or 2 * vR + d == vB)
    ok += 1 if bien else 0
    mal += 0 if bien else 1
    print("%-5d %-3d %-4d %-5d %-3d %-4d %-7d %-9d %-9d %-9s %s" % (
        q, p, r, n, d, delta, vB, pred, plano,
        str(vR) if vR >= 0 else "?", "ok" if bien else "*** DISCREPA ***"))
    sys.stdout.flush()
print()
print("v_p([O:B]) == 3d-2+2delta  y  v_p = 2 v_p(R) + d :  %d ok, %d discrepan" % (ok, mal))
