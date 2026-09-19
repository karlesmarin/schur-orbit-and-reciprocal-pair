# -*- coding: utf-8 -*-
# W_1 EN TIPO A, Y SI LA LEY LOCAL DE TIPO C SE TRASLADA.
#
# En tipo C el paper prueba  v_p([O:A]) = 2d - 1 - W_1  con d = phi(q')/2, bajo W_1 > d/2 (o
# W_1 > 0 y W_2 = d).  Aqui se calcula W_1 en tipo A con `inst_tipoA.sage` y se contrasta contra
# el v_p MEDIDO por reticulos, que es independiente y ya esta validado contra 48 casillas del paper.
#
# CONTROL C0, primero y sin el cual lo demas no vale: el instrumento nuevo, alimentado con el
# segmento de TIPO C, tiene que devolver los mismos W que `inst.sage` y cumplir la ley publicada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: docker ... sage capa1_tipoA.sage
load("inst_tipoA.sage")


def cierre_indice(exps, M, r):
    """v_p del indice por reticulos: el instrumento independiente."""
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


print("=" * 96)
print("C0  CONTROL: el instrumento nuevo con segmento de TIPO C reproduce la ley 2d-1-W_1")
print("=" * 96)
print("%-3s %-4s %-3s %-4s %-4s %-6s %-8s %-8s %-8s %s" % (
    "m", "q", "p", "q'", "d", "W_1", "2d-1-W1", "v_p real", "W_1>d/2", "veredicto"))
ok = malo = fuera = 0
for m in range(2, 5):
    for q in range(5, 46):
        if euler_phi(q) > 24:
            continue
        idx = cierre_indice(exps_C(m), q, m)
        if idx == 1:
            continue
        for p, _ in factor(idx):
            v = valuation(q, p)
            qp = q // p ** v
            if qp < 3 or euler_phi(p ** v) < 2:
                continue
            d = euler_phi(qp) / 2
            try:
                I = ImgSeg(q, p, exps_C(m), m, 2)
            except Exception as e:
                print("   m=%d q=%d p=%d ERROR %s" % (m, q, p, type(e).__name__))
                continue
            W = I.W()
            W1 = W[1] if len(W) > 1 else None
            vp = valuation(idx, p)
            pred = 2 * d - 1 - W1 if W1 is not None else None
            hip = (W1 is not None and W1 > d / 2)
            ver = "-"
            if pred is not None:
                if hip:
                    ver = "ok" if pred == vp else "*** DISCREPA ***"
                    if pred == vp:
                        ok += 1
                    else:
                        malo += 1
                else:
                    ver = "fuera de hipotesis"
                    fuera += 1
            print("%-3d %-4d %-3d %-4d %-4s %-6s %-8s %-8d %-8s %s" % (
                m, q, p, qp, d, W1, pred, vp, "si" if hip else "no", ver))
print("   bajo hipotesis: %d ok, %d discrepan | fuera de hipotesis: %d" % (ok, malo, fuera))
if malo == 0 and ok > 0:
    print("   el instrumento nuevo reproduce la ley local de tipo C.")

print()
print("=" * 96)
print("TIPO A (SU).  W_1 medido, y la ley de tipo C puesta a prueba fuera de su tipo")
print("=" * 96)
print("%-3s %-4s %-3s %-4s %-4s %-6s %-8s %-8s %-8s %s" % (
    "n", "q", "p", "q'", "d", "W_1", "2d-1-W1", "v_p real", "W_1>d/2", "familia"))
okA = maloA = fueraA = 0
for n in range(3, 10):
    for q in range(3, 61):
        M = q if n % 2 == 1 else 2 * q
        if euler_phi(M) > 32:
            continue
        idx = cierre_indice(exps_su(n), M, n - 1)
        if idx == 1:
            continue
        for p, _ in factor(idx):
            v = valuation(M, p)
            qp = M // p ** v
            if qp < 3 or euler_phi(p ** v) < 2:
                continue
            d = euler_phi(qp) / 2
            try:
                I = ImgSeg(M, p, exps_su(n), n - 1, 2)
            except Exception as e:
                print("   n=%d q=%d p=%d ERROR %s" % (n, q, p, type(e).__name__))
                continue
            W = I.W()
            W1 = W[1] if len(W) > 1 else None
            vp = valuation(idx, p)
            pred = 2 * d - 1 - W1 if W1 is not None else None
            hip = (W1 is not None and W1 > d / 2)
            qpq = q // p ** valuation(q, p)
            fam = "/".join([nom for nom, x in (("n-1", n - 1), ("n", n), ("n+1", n + 1))
                            if qpq >= 3 and x % qpq == 0]) or "-"
            if pred is not None and hip:
                if pred == vp:
                    okA += 1
                else:
                    maloA += 1
            elif pred is not None:
                fueraA += 1
            print("%-3d %-4d %-3d %-4d %-4s %-6s %-8s %-8d %-8s %s" % (
                n, q, p, qp, d, W1, pred, vp, "si" if hip else "no", fam))
print()
print("   TIPO A bajo la hipotesis W_1 > d/2: %d aciertos, %d discrepancias | fuera: %d" % (
    okA, maloA, fueraA))
if maloA == 0 and okA > 0:
    print("   *** la ley local de tipo C se traslada al tipo A en toda la muestra ***")
elif maloA:
    print("   la ley NO se traslada tal cual: hay %d discrepancias, que es el dato." % maloA)
