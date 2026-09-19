# -*- coding: utf-8 -*-
# DATOS PARA LA FIGURA DEL PLANO, EN RANGO ANCHO.
#
# La figura  v_p = 2 dim - 1 - W_1  es la tesis hecha imagen, pero con el barrido corto las 140
# casillas ocupan solo 11 posiciones distintas y un 3D con once puntos no se gana las tres
# dimensiones.  Aqui se ensancha el rango para que el plano tenga puntos de verdad.
#
# Imprime, para los dos anillos:  n q p q' dim W_1 v_p   con dim = phi(q')/2 en el real y phi(q')
# en la sombra.  El v_p sale de los divisores elementales del reticulo, no de la formula.
#
# Authors: Carles Marin, Claude (AI assistant).
load("inst_tipoA.sage")


def indice(exps, M, r):
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
            break
        L = L2
    return prod([x for x in L.elementary_divisors() if x != 0])


print("# tipo n q p qprima dim W1 vp")
for n in range(3, 15):
    for q in range(3, 81):
        # --- real
        M = q if n % 2 == 1 else 2 * q
        if euler_phi(M) <= 48:
            try:
                idx = indice(exps_su(n), M, n - 1)
            except Exception:
                idx = 1
            if idx > 1:
                for p, _ in factor(idx):
                    v = valuation(M, p)
                    qp = M // p ** v
                    if qp < 3 or euler_phi(p ** v) < 2:
                        continue
                    try:
                        W = ImgSeg(M, p, exps_su(n), n - 1, 2).W()
                    except Exception:
                        continue
                    if len(W) > 1:
                        print("real %d %d %d %d %s %d %d" % (
                            n, q, p, qp, euler_phi(qp) / 2, W[1], valuation(idx, p)))
        # --- sombra
        if euler_phi(q) <= 48:
            try:
                idg = indice([j for j in range(n)], q, n)
            except Exception:
                idg = 1
            if idg > 1:
                for p, _ in factor(idg):
                    v = valuation(q, p)
                    qp = q // p ** v
                    if qp < 3 or euler_phi(p ** v) < 2:
                        continue
                    try:
                        W = ImgSeg(q, p, [j for j in range(n)], n, 2).W()
                    except Exception:
                        continue
                    if len(W) > 1:
                        print("sombra %d %d %d %d %d %d %d" % (
                            n, q, p, qp, euler_phi(qp), W[1], valuation(idg, p)))
