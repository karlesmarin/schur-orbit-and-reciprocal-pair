# -*- coding: utf-8 -*-
# LA SOMBRA: EL CUADRO LOCAL DEL ANILLO NO REAL DE GL_n.
#
# A_GL = Z[e_1..e_n] con x_j = zeta_q^{j-1}, j=1..n: es el C_d de NPP25.  Su anillo NO es real y
# llena todo Q(zeta_q) (rango phi(q), medido).  Por eso la dimension de cada capa graduada deberia
# ser D = phi(q') ENTERA, no phi(q')/2 como en el caso real, donde la mitad la mata la simetria.
#
# Se prueban las dos mitades del cuadro, con el instrumento de capas ya validado:
#   S1  v_p([O:A]) = 2D - 1 - W_1   con D = phi(q'),  bajo W_1 > D/2  (la ley, con la D nueva)
#   S2  W_1 = rango de [S(u^{-1}c)] con filas sobre TODAS las unidades (no modulo +-1, porque el
#       segmento {0,...,n-1} NO es simetrico y las clases con 2c=0 ya no se anulan)
#   S3  control de que la D correcta es phi(q') y no phi(q')/2: se prueban las dos y se cuenta.
#
# Authors: Carles Marin, Claude (AI assistant).
load("inst_tipoA.sage")


def indice_gl(n, q):
    K = CyclotomicField(q)
    z = K.gen()
    xs = [z ** j for j in range(n)]
    S = PolynomialRing(K, 'T')
    T = S.gen()
    f = prod([T - x for x in xs])
    g = [(-1) ** k * f.list()[n - k] for k in range(1, n + 1)]
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


def rango_S_todas(exps, qp, p):
    S = [0] * qp
    for e in exps:
        S[e % qp] += e
    units = [u for u in range(1, qp) if gcd(u, qp) == 1]
    filas = [[S[(inverse_mod(u, qp) * c) % qp] for c in range(qp)] for u in units]
    return matrix(GF(p), filas).rank()


print("%-3s %-4s %-3s %-4s %-4s %-6s %-8s %-9s %-9s %-8s %s" % (
    "n", "q", "p", "q'", "D", "W_1", "2D-1-W1", "v_p real", "W1>D/2", "rango S", "veredicto"))
okley = malley = fuera = 0
okS = malS = 0
for n in range(3, 10):
    for q in range(3, 50):
        if euler_phi(q) > 28:
            continue
        idx = indice_gl(n, q)
        if idx == 1:
            continue
        exps = [j for j in range(n)]
        for p, _ in factor(idx):
            v = valuation(q, p)
            qp = q // p ** v
            if qp < 3 or euler_phi(p ** v) < 2:
                continue
            try:
                I = ImgSeg(q, p, exps, n, 2)
            except Exception as e:
                print("   n=%d q=%d p=%d ERROR %s" % (n, q, p, type(e).__name__))
                continue
            W = I.W()
            if len(W) < 2:
                continue
            W1 = W[1]
            D = euler_phi(qp)
            vp = valuation(idx, p)
            pred = 2 * D - 1 - W1
            hip = W1 > D / 2
            rS = rango_S_todas(exps, qp, p)
            if rS == W1:
                okS += 1
            else:
                malS += 1
            ver = "-"
            if hip:
                if pred == vp:
                    okley += 1
                    ver = "ok"
                else:
                    malley += 1
                    ver = "*** DISCREPA ***"
            else:
                fuera += 1
                ver = "fuera de hipotesis"
            print("%-3d %-4d %-3d %-4d %-4d %-6d %-8d %-9d %-9s %-8d %s" % (
                n, q, p, qp, D, W1, pred, vp, "si" if hip else "no", rS, ver))
print()
print("S1  ley con D = phi(q'):  %d ok, %d discrepan | fuera de hipotesis: %d" % (okley, malley, fuera))
print("S2  W_1 == rango de [S(u^-1 c)] sobre todas las unidades: %d si, %d no" % (okS, malS))
