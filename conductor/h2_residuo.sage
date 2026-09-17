# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# H2 Y H2' (prediccion escrita antes de los datos): LA IMAGEN DE A EN LOS CUERPOS RESIDUALES DEL CONDUCTOR.
#
# Para cada casilla (m,q) con conductor no trivial y cada p | N(c):
#   a = v_p([O+:A]) ,  b = v_p([A:c])  -- por reticulos (og_orden + conductor, como factorizar_conductor)
#   los P | c sobre p con (e,f)       -- maximizando SOLO en p
#   dim = dim_Fp de la imagen de A en prod_{P|c, P|p} O/P  -- cierre de F_p-subalgebra, en cuerpos
#         finitos: no usa ni el reticulo del conductor ni HNF.
#
# Imprime una linea maquina por (casilla,p):
#   H2|m|q|p|a|b|dim|forma|claseH2|claseH2p
# forma = "f.e+f.e..." ; claseH2 in {na (algun e>=2), forzado, libre} ; claseH2p in {forzado, libre}.
# La clase se decide SOLO con (e,f) -- antes de mirar b y dim.
#
# CONTROLES (resumen al final de cada casilla; el lector los suma)
#   T    teorema: si todo e = 1, b == dim.  Fallo = instrumento roto.
#   S1   senuelo: theta con reduccion generadora de F_{p^k}, 1<k<f (f compuesto) -> dim k.
#   S2   senuelo: theta con reduccion generadora de F_{p^f} -> dim f.
#        Se hacen en el primer P | c sobre p con f >= 2.  Si el instrumento devolviera 1 siempre, H2'
#        "pasaria" sin medir nada.
#
# REPL-SAFE.  Nombres sin pisar globales de sage (nada de prod/sum/list/vector/matrix = ...).
#
# Authors: Carles Marin, Claude (AI assistant).
import os
import sys
load("orden_generado.sage")
x = polygen(QQ, 'x')
def psi(q):
    K = CyclotomicField(q)
    z = K.gen()
    return (z + z**(-1)).minpoly(x)
def simetricas(m, L):
    a = L.gen()
    et = [a]
    if m >= 2:
        et.append(a*a - 2)
    for j in range(2, m):
        et.append(a*et[j-1] - et[j-2])
    et = et[:m]
    R = PolynomialRing(L, 'T')
    T = R.gen()
    P = prod([T - w for w in et])
    return [(-1)**i * P[m-i] for i in range(1, m+1)]
def reticulo_conductor(L, A):
    d = L.degree()
    a = L.gen()
    BA = matrix(QQ, [list(L(bb)) for bb in A.basis()])
    Tm = matrix(QQ, [list(a*(a**k)) for k in range(d)])
    blocks = []
    Tk = identity_matrix(QQ, d)
    iBA = BA.inverse()
    for k in range(d):
        blocks.append(Tk * iBA)
        Tk = Tk * Tm
    Mall = block_matrix(1, d, blocks)
    N = lcm([QQ(t).denominator() for t in Mall.list()])
    if N == 1:
        return None
    B = matrix(ZZ, N*Mall)
    n = B.ncols()
    St = block_matrix(ZZ, 2, 1, [B, -N*identity_matrix(ZZ, n)])
    K = St.left_kernel().basis_matrix()
    return K[:, 0:d].hermite_form(include_zero_rows=False)
def vp(n, p):
    n = ZZ(n)
    return n.valuation(p) if n != 0 else 0
def dim_imagen(gens, primos, p):
    # gens: elementos de L2 enteros en cada P.  primos: lista de (P, kP, V, to_V).
    # Espacio ambiente: prod kP como F_p-espacio (suma de dimensiones).  Cierre por multiplicacion.
    F = GF(p)
    def red(z):
        return [kP(z) for (P, kP, V, to_V) in primos]
    def vec(r):
        cs = []
        for (t, (P, kP, V, to_V)) in zip(r, primos):
            cs.extend(list(to_V(t)))
        return vector(F, cs)
    def mul(r1, r2):
        return [u*v for (u, v) in zip(r1, r2)]
    uno = [kP(1) for (P, kP, V, to_V) in primos]
    gr = [red(g) for g in gens]
    elems = [uno]
    W = matrix(F, [vec(uno)]).row_space()
    while True:
        nuevos = list(elems)
        for r in elems:
            for g in gr:
                nuevos.append(mul(r, g))
        # base de elementos (no de vectores) que genera el span
        filas = []
        base = []
        for r in nuevos:
            v = vec(r)
            if len(filas) == 0 or matrix(F, filas + [v]).rank() > len(filas):
                filas.append(v)
                base.append(r)
        if len(base) == W.dimension():
            break
        elems = base
        W = matrix(F, filas).row_space()
    return W.dimension()
def preparar(Ps):
    out = []
    for P in Ps:
        kP = P.residue_field()
        V, from_V, to_V = kP.vector_space(map=True)
        out.append((P, kP, V, to_V))
    return out
CAS = []
for t in os.environ.get("CASILLAS", "6:39,7:56,4:27").split(","):
    ab = t.strip().split(":")
    CAS.append((ZZ(ab[0]), ZZ(ab[1])))
for (m, q) in CAS:
    f = psi(q)
    L = NumberField(f, 'a')
    d = L.degree()
    a = L.gen()
    A = og_orden(L, simetricas(m, L))
    OL = L.order(a)
    idx = ZZ(A.index_in(OL))
    Lat = reticulo_conductor(L, A)
    if Lat is None:
        print("H2NADA|%d|%d|conductor trivial" % (m, q))
        sys.stdout.flush()
        continue
    nc = ZZ(abs(Lat.determinant()))
    tmal = 0
    tmed = 0
    for (p, vv) in nc.factor():
        L2 = NumberField(f, 'b', maximize_at_primes=[p])
        c = L2.ideal([L2(list(v)) for v in Lat.rows()])
        fac = [(P, ZZ(ee)) for (P, ee) in c.factor() if ZZ(P.residue_field().characteristic()) == p]
        ef = [(ZZ(P.residue_class_degree()), ee) for (P, ee) in fac]
        aa = vp(idx, p)
        bb = vp(nc, p) - aa
        primos = preparar([P for (P, ee) in fac])
        gens2 = [L2(list(s)) for s in simetricas(m, L)]
        dm = dim_imagen(gens2, primos, p)
        forma = "+".join(["%d.%d" % (ff, ee) for (ff, ee) in ef])
        r = len(ef)
        todos_e1 = all(ee == 1 for (ff, ee) in ef)
        fs = [ff for (ff, ee) in ef]
        if not todos_e1:
            cH2 = "na"
        elif (r == 1 and fs[0].is_prime()) or (r == 2 and fs == [1, 1]):
            cH2 = "forzado"
        else:
            cH2 = "libre"
        if (r == 1 and fs[0] == 1) or (todos_e1 and r == 1 and fs[0].is_prime()) or (todos_e1 and r == 2 and fs == [1, 1]):
            cH2p = "forzado"
        else:
            cH2p = "libre"
        print("H2|%d|%d|%d|%d|%d|%d|%s|%s|%s" % (m, q, p, aa, bb, dm, forma, cH2, cH2p))
        if todos_e1:
            tmed += 1
            if bb != dm:
                tmal += 1
                print("H2T_FALLO|%d|%d|%d|b=%d|dim=%d" % (m, q, p, bb, dm))
        # senuelos en el primer P con f >= 2
        for (P, kP, V, to_V) in primos:
            ff = ZZ(kP.degree())
            if ff < 2:
                continue
            g = kP.multiplicative_generator()
            th_full = kP.lift(g)
            d2 = dim_imagen([L2(th_full)], [(P, kP, V, to_V)], p)
            print("H2S2|%d|%d|%d|f=%d|dim=%d|%s" % (m, q, p, ff, d2, "ok" if d2 == ff else "MAL"))
            ks = [k for k in ff.divisors() if 1 < k < ff]
            if ks:
                k = ks[0]
                th = kP.lift(g**((p**ff - 1) // (p**k - 1)))
                d1 = dim_imagen([L2(th)], [(P, kP, V, to_V)], p)
                print("H2S1|%d|%d|%d|f=%d|k=%d|dim=%d|%s" % (m, q, p, ff, k, d1, "ok" if d1 == k else "MAL"))
            break
        sys.stdout.flush()
    print("H2FIN|%d|%d|idx=%d|Nc=%d|T_medidas=%d|T_fallos=%d" % (m, q, idx, nc, tmed, tmal))
    sys.stdout.flush()
