# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# FACTORIZA EL CONDUCTOR EN LAS CASILLAS CARAS, Y PRUEBA P1 DONDE DUELE.
#
# POR QUE EXISTE.  `ley_conductor.sage` muere con OOM (137) en grado alto, y muere en el mismo sitio
# las dos veces: `L.maximal_order()` y `ideal.factor()`, que fuerzan la maximalizacion global por
# round4.  Pero el conductor se calcula SIN orden maximal (potencias de a, y O+ = Z[a] es teorema),
# y su soporte son unos pocos primos.  Luego se puede maximizar SOLO en esos primos.
#
# LO QUE SE PRUEBA.  P1, pre-registrada ANTES de mirar estos datos:
#
#     [A : c]  =  prod_{P | c}  p_P        (un factor por IDEAL primo, exponente 1)
#     [O+ : A] =  prod_P p^(f_P e_P - 1)   (consecuencia, P2)
#
# y la ley VIEJA -- prod_p p^ceil(v_p/2) -- coincide con P2 exactamente cuando f.e <= 3 para todo P.
# Por eso aguantaba 26 de 26 en un barrido que solo tenia f.e <= 3.  No era una ley.
#
# EL CASO QUE DECIDE: (4,63), grado 18, N(c) = 343 = 7^3 con [A:c] medido = 7.  Si c es UN ideal con
# f.e = 3, P1 predice 7 y acierta.  Si fueran TRES ideales de f.e = 1, P1 predice 343 y MUERE.  La
# medida distingue las dos cosas; la fila del barrido ancho no, porque no pudo factorizar.
#
# CONTROLES
#   F1  COHERENCIA prod_P N(P)^e == N(c), con N(c) leido del determinante de la HNF.  Dos rutas.
#   F2  COHERENCIA [O+:A].[A:c] == N(c).
#   F3  P1 :  [A:c] == prod_{P|c} p.  Se cuenta SOLO sobre las casillas NO Gorenstein: en las
#       Gorenstein [A:c] = [O+:A] = raiz(N(c)) es teorema y P1 acertaria sin medir nada.
#   F4  P2 :  [O+:A] == prod_P p^(f e - 1), que es F3 leida del otro lado via F2 -- NO es un control
#       independiente y se dice aqui para que no se cuente dos veces.
#   F5  SENUELO, la ley vieja prod_p p^ceil(v_p/2).  Tiene que FALLAR en las casillas con f.e >= 4.
#       Si no hay ninguna casilla con f.e >= 4 en la lista, F5 no distingue nada y lo dice.
#
# REPL-SAFE: sin lineas en blanco dentro de bloques indentados.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run (detached):
#   MSYS_NO_PATHCONV=1 docker run -d --name facc --memory=10g -e CASILLAS="4:63,4:72,5:33,5:44" \
#     -v "$PWD:/work" -w /work sage-normaliz:local \
#     bash -c "sage < factorizar_conductor.sage > factorizar_conductor_OUT.txt 2>&1"
import os
import sys
load("orden_generado.sage")
x = polygen(QQ, 'x')
def psi(q):
    K = CyclotomicField(q)
    z = K.gen()
    return (z + z**(-1)).minpoly(x)
def etas_de(m, L):
    a = L.gen()
    out = [a]
    if m >= 2:
        out.append(a*a - 2)
    for j in range(2, m):
        out.append(a*out[j-1] - out[j-2])
    return out[:m]
def simetricas(m, L):
    R = PolynomialRing(L, 'T')
    T = R.gen()
    P = prod([T - w for w in etas_de(m, L)])
    return [(-1)**i * P[m-i] for i in range(1, m+1)]
def reticulo_conductor(L, A):
    # c = { x en Z[a] : x.Z[a] contenido en A }.  Z[a] = O+ es TEOREMA (Washington Prop. 2.16) y es
    # la hipotesis que hace valida esta cuenta sin tocar el orden maximal.
    d = L.degree()
    a = L.gen()
    BA = matrix(QQ, [list(L(b)) for b in A.basis()])
    T = matrix(QQ, [list(a*(a**k)) for k in range(d)])
    blocks = []
    Tk = identity_matrix(QQ, d)
    for k in range(d):
        blocks.append(Tk * BA.inverse())
        Tk = Tk * T
    Mall = block_matrix(1, d, blocks)
    N = lcm([QQ(t).denominator() for t in Mall.list()])
    if N == 1:
        return None
    B = matrix(ZZ, N*Mall)
    n = B.ncols()
    St = block_matrix(ZZ, 2, 1, [B, -N*identity_matrix(ZZ, n)])
    K = St.left_kernel().basis_matrix()
    return K[:, 0:d].hermite_form(include_zero_rows=False)
CAS = []
for t in os.environ.get("CASILLAS", "4:63,4:72,5:33,5:44").split(","):
    ab = t.split(":")
    CAS.append((ZZ(ab[0]), ZZ(ab[1])))
print("=" * 118)
print("F1-F5   el conductor factorizado en las casillas caras, y P1 puesta a prueba donde puede morir")
print("=" * 118)
print("   casillas pedidas: %s" % os.environ.get("CASILLAS", "4:63,4:72,5:33,5:44"))
print("   m    q   d   [O+:A]     N(c)   [A:c]  Goren  primos de c (p,f,e)                P1    P2   leyVieja")
print("   " + "-" * 112)
f1mal = 0
f2mal = 0
f3med = 0
f3mal = 0
f5med = 0
f5mal = 0
nogor = 0
gor = 0
fe_altos = 0
saltadas = {}
for (m, q) in CAS:
    try:
        f = psi(q)
        L = NumberField(f, 'a')
        d = L.degree()
        a = L.gen()
        # og_orden y no L.order: L.order construye d^m monomios y ESO eran los OOM de cel6/7/8
        A = og_orden(L, simetricas(m, L))
        OL = L.order(a)
        idx = ZZ(A.index_in(OL))
        Lat = reticulo_conductor(L, A)
    except Exception as ex:
        saltadas["excepcion: " + str(ex)[:50]] = saltadas.get("excepcion: " + str(ex)[:50], 0) + 1
        continue
    if Lat is None:
        saltadas["conductor trivial"] = saltadas.get("conductor trivial", 0) + 1
        continue
    nc = ZZ(abs(Lat.determinant()))
    acx = nc // idx
    # maximalizar SOLO en los primos del soporte de c: es lo unico que hace falta para factorizarlo,
    # y es lo que evita el round4 global que reventaba la memoria.
    ps = [pp for (pp, vv) in nc.factor()]
    try:
        L2 = NumberField(f, 'b', maximize_at_primes=ps)
        c = L2.ideal([L2(list(v)) for v in Lat.rows()])
        fac = [(ZZ(P.residue_field().characteristic()), ZZ(P.residue_class_degree()), ZZ(e)) for (P, e) in c.factor()]
    except Exception as ex:
        saltadas["factorizacion fallo: " + str(ex)[:40]] = saltadas.get("factorizacion fallo: " + str(ex)[:40], 0) + 1
        continue
    normaf = ZZ(1)
    p1 = ZZ(1)
    p2 = ZZ(1)
    maxfe = 0
    for (pp, ff, ee) in fac:
        normaf = normaf * pp**(ff*ee)
        p1 = p1 * pp
        p2 = p2 * pp**(ff*ee - 1)
        if ff*ee > maxfe:
            maxfe = ff*ee
    if normaf != nc:
        f1mal += 1
    if idx*acx != nc:
        f2mal += 1
    vieja = ZZ(1)
    for (pp, vv) in nc.factor():
        vieja = vieja * pp**ceil(vv/2)
    es_gor = (idx == acx)
    if es_gor:
        gor += 1
    else:
        nogor += 1
        f3med += 1
        if p1 != acx:
            f3mal += 1
    if maxfe >= 4:
        fe_altos += 1
        if not es_gor:
            f5med += 1
            if vieja != idx:
                f5mal += 1
    fs = "  ".join(["(%d,%d,%d)" % t for t in fac])
    print("  %3d %4d %3d %8s %8s %7s  %-5s  %-40s %-5s %-5s %-5s" % (m, q, d, str(idx), str(nc), str(acx), "SI" if es_gor else "no", fs, "ok" if p1 == acx else "NO", "ok" if p2 == idx else "NO", "ok" if vieja == idx else "NO"))
    sys.stdout.flush()
print("")
print("   REGLA: cada control cuenta sus MEDIDAS.  Un veredicto con 0 medidas no es un veredicto.")
for (k, v) in sorted(saltadas.items()):
    print("     descartadas por %-52s %d" % (k, v))
print("   F1  COHERENCIA prod N(P)^e == N(c) :  fallos %d" % f1mal)
print("   F2  COHERENCIA [O+:A].[A:c] == N(c) :  fallos %d" % f2mal)
if f3med == 0:
    vf3 = "*** 0 MEDIDAS: todas las casillas de la lista son Gorenstein y P1 no se ha probado ***"
elif f3mal == 0:
    vf3 = "AGUANTA sobre %d casillas NO Gorenstein" % f3med
else:
    vf3 = "*** MUERE ***"
print("   F3  P1  [A:c] == prod_{P|c} p :  medidas %d (solo NO Gorenstein)  fallos %d   %s" % (f3med, f3mal, vf3))
print("   F4  P2  [O+:A] == prod p^(fe-1) :  NO es control aparte -- es F3 dividida por F2.  No se cuenta.")
if f5med == 0:
    vf5 = "*** sin casilla no-Gorenstein con f.e >= 4, la ley vieja y P2 son la MISMA formula aqui: F5 no distingue ***"
elif f5mal == f5med:
    vf5 = "PASA -- la ley vieja falla en las %d casillas donde las dos formulas difieren" % f5med
else:
    vf5 = "*** la ley vieja ACIERTA en %d casillas con f.e >= 4: P2 no la sustituye ***" % (f5med - f5mal)
print("   F5  SENUELO ley vieja prod p^ceil(v/2) :  casillas con algun f.e >= 4: %d ; medidas (y no Gorenstein) %d ; falla en %d   %s" % (fe_altos, f5med, f5mal, vf5))
print("   Gorenstein %d / no Gorenstein %d en la lista." % (gor, nogor))
