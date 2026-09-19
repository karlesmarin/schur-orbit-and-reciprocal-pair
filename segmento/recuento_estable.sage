# -*- coding: utf-8 -*-
# EL RECUENTO BUENO: la descripcion correcta, restringida a la hipotesis del teorema.
#
# El 285/317 de esta manana mezclaba dos cosas que no debia: usaba el PROXY (rango de la matriz de
# momentos en el espacio de funciones) y contaba casillas SIN exigir estabilidad.  Aqui se rehace:
#
#   objeto correcto : dim_{F_p} del modulo F_p[G].g dentro de R = F_p[X]/Phi_{q'},  g = sum S(c) w^c
#   proxy           : rango de [S(u^{-1}c)] en el espacio de funciones
#   hipotesis       : segmento ESTABLE en q'  (la del Teorema)
#
# Se dan las cuatro cifras -- las dos descripciones, dentro y fuera de la hipotesis -- para que el
# numero que se publique sea el que corresponde y se vea lo que aporta la hipotesis.
#
# Authors: Carles Marin, Claude (AI assistant).
load("inst_tipoA.sage")


def S_de(exps, qp):
    S = [0] * qp
    for e in exps:
        S[e % qp] += e
    return S


def estable(exps, qp):
    base = sorted([e % qp for e in exps])
    for u in range(1, qp):
        if gcd(u, qp) == 1 and sorted([(u * e) % qp for e in exps]) != base:
            return False
    return True


def dim_en_R(S, qp, p):
    F = GF(p)
    R = PolynomialRing(F, 'y')
    y = R.gen()
    phi = R(cyclotomic_polynomial(qp))
    D = phi.degree()
    units = [u for u in range(1, qp) if gcd(u, qp) == 1]
    filas = []
    for u in units:
        el = sum(F(S[c]) * y ** ((u * c) % qp) for c in range(qp)) % phi
        v = el.list()
        filas.append(v + [F(0)] * (D - len(v)))
    return matrix(F, filas).rank()


def rango_proxy(S, qp, p):
    units = [u for u in range(1, qp) if gcd(u, qp) == 1 and u < qp - u]
    if not units:
        return 0
    filas = [[S[(inverse_mod(u, qp) * c) % qp] for c in range(qp)] for u in units]
    return matrix(GF(p), filas).rank()


cnt = {("est", "R"): 0, ("est", "px"): 0, ("ine", "R"): 0, ("ine", "px"): 0}
tot = {"est": 0, "ine": 0}
fallosR = []
for n in range(3, 12):
    for q in range(3, 50):
        M = q if n % 2 == 1 else 2 * q
        if euler_phi(M) > 28:
            continue
        for p, _ in factor(M):
            v = valuation(M, p)
            qp = M // p ** v
            if qp < 3 or euler_phi(p ** v) < 2:
                continue
            ex = exps_su(n)
            try:
                I = ImgSeg(M, p, ex, n - 1, 2)
                W = I.W()
                if len(W) < 2:
                    continue
                w = W[1]
            except Exception:
                continue
            S = S_de(ex, qp)
            est = estable(ex, qp)
            k = "est" if est else "ine"
            tot[k] += 1
            if dim_en_R(S, qp, p) == w:
                cnt[(k, "R")] += 1
            elif est:
                fallosR.append((n, q, p, qp, w, dim_en_R(S, qp, p), rango_proxy(S, qp, p)))
            if rango_proxy(S, qp, p) == w:
                cnt[(k, "px")] += 1

print("RECUENTO, cuatro cifras:")
print()
print("%-34s %-16s %s" % ("", "BAJO estabilidad", "fuera de ella"))
print("%-34s %-16s %s" % ("modulo dentro de R (lo correcto)",
                          "%d / %d" % (cnt[("est", "R")], tot["est"]),
                          "%d / %d" % (cnt[("ine", "R")], tot["ine"])))
print("%-34s %-16s %s" % ("proxy en el espacio de funciones",
                          "%d / %d" % (cnt[("est", "px")], tot["est"]),
                          "%d / %d" % (cnt[("ine", "px")], tot["ine"])))
print()
if fallosR:
    print("fallos de la descripcion correcta DENTRO de la hipotesis (%d):" % len(fallosR))
    for f in fallosR[:12]:
        print("   n=%d q=%d p=%d q'=%d  W_1=%d  dim_R=%d  proxy=%d" % f)
else:
    print("*** la descripcion correcta no falla ni una vez bajo estabilidad ***")
