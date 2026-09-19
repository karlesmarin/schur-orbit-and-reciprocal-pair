# -*- coding: utf-8 -*-
# ¿SON EL MISMO ANILLO, O SOLO TIENEN EL MISMO CONDUCTOR?
#
# Anadir un autovalor igual a 1 al segmento cambia los generadores por
#     e_k(x u {1}) = e_k(x) + e_{k-1}(x),
# que es triangular sobre Z con unos en la diagonal.  Luego  Z[e_1^new, ..., e_n^new] = Z[e_1,...,e_r],
# y el anillo de SU(2m+1) en la recta de rho seria LITERALMENTE el de Sp(2m).  Si es asi, no hay que
# trasladar ningun teorema: todo el paper de tipo C se aplica al mismo objeto con otro nombre.
#
# Se comprueba comparando los RETICULOS, no los indices: dos anillos distintos pueden compartir indice.
#
# Authors: Carles Marin, Claude (AI assistant).
load("inst_tipoA.sage")


def reticulo(exps, M, r):
    K = CyclotomicField(M)
    z = K.gen()
    xs = [z ** (e % M) for e in exps]
    S = PolynomialRing(K, 'T')
    T = S.gen()
    f = prod([T - x for x in xs])
    nn = len(xs)
    g = [(-1) ** k * f.list()[nn - k] for k in range(1, r + 1)]
    g = [y for y in g if y != 0]
    L = matrix(ZZ, [K(1).list()]).hermite_form(include_zero_rows=False)
    for _ in range(40):
        filas = [list(fi) for fi in L.rows()]
        for fi in L.rows():
            el = K(list(fi))
            for gg in g:
                filas.append((el * gg).list())
        L2 = matrix(ZZ, filas).hermite_form(include_zero_rows=False)
        if L2 == L:
            return L
        L = L2
    return L


print("A_SU(2m+1)(q)  ==  A_Sp(2m)(q)  como RETICULOS dentro de Z[zeta_q] ?")
print("%-3s %-4s %-7s %-7s %-9s %s" % ("m", "q", "rango A", "rango C", "iguales", "indice comun"))
ok = malo = 0
for m in range(1, 8):
    n = 2 * m + 1
    for q in range(3, 50):
        if euler_phi(q) > 28:
            continue
        LA = reticulo(exps_su(n), q, n - 1)
        LC = reticulo(exps_C(m), q, m)
        igual = (LA == LC)
        if igual:
            ok += 1
        else:
            malo += 1
        if not igual or q <= 9:
            idx = prod([x for x in LA.elementary_divisors() if x != 0])
            print("%-3d %-4d %-7d %-7d %-9s %s" % (
                m, q, LA.nrows(), LC.nrows(), "si" if igual else "*** NO ***", idx))
print()
print("casillas: %d | iguales: %d | distintas: %d" % (ok + malo, ok, malo))
if malo == 0:
    print("*** SON EL MISMO ANILLO: no hay teorema que trasladar, hay objeto que reconocer ***")

print()
print("Y n PAR: ¿el de SU(2m) se relaciona igual con algun segmento de tipo C?")
print("   (los exponentes son los impares +-1..+-(2m-1) sobre M = 2q; se compara con el segmento")
print("    de Sp(2m-1)?? no existe; se compara con el DESPLAZADO de tipo C que el paper ya trata)")
for n in (4, 6):
    for q in (7, 9, 11, 15):
        M = 2 * q
        if euler_phi(M) > 28:
            continue
        L = reticulo(exps_su(n), M, n - 1)
        idx = prod([x for x in L.elementary_divisors() if x != 0])
        print("   n=%d q=%-3d  M=%-3d rango=%-3d indice=%s   exps mod q = %s" % (
            n, q, M, L.nrows(), idx, sorted([e % q for e in exps_su(n)])))
