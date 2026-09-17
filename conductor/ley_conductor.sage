# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# LA LEY DEL INDICE, PARTIDA EN TEOREMA + RESTO:  donde A es Gorenstein no hay ley que ajustar.
#
# DE DONDE VIENE.  `ley_ideal.sage` midio  [O+:A] = prod_P p^ceil(f_P e_P / 2)  sobre los IDEALES
# PRIMOS P del conductor, y MURIO en una sola casilla de 16: (3,56), que es la unica de la muestra
# donde un primo racional aporta DOS ideales -- c = P_2^2 . P_7 . P_7'.  Ahi la version por ideal
# predice 7.7 = 49 y el indice solo se lleva 7.  La correccion es agregar POR PRIMO RACIONAL:
#
#     [O+ : A] = prod_p p^ceil( v_p([O+:c]) / 2 )   =   el menor n tal que [O+:c] divide n^2
#
# que no menciona f ni e ni ideales: solo la valoracion de un ENTERO, que es donde vive un ceil.
#
# EL TEOREMA QUE HAY DEBAJO, Y POR QUE ESTO NO ES TODO AJUSTE.  Como c es un ideal de O+ contenido
# en A, la cadena c <= A <= O+ da  [O+:c] = [O+:A].[A:c].  La ley dice que esos dos factores son
# iguales salvo los primos de valoracion IMPAR.  Y "[O:A] = [A:c]" NO es una observacion: es el
# criterio clasico de que el orden A sea GORENSTEIN.  Luego:
#
#   - en las casillas Gorenstein la ley NO es una ley ajustada, es ese teorema, y  [O+:A]^2 = [O+:c]
#     EXACTAMENTE, sin techo;
#   - el techo solo hace trabajo en las casillas NO Gorenstein.
#
# Esto es lo que el gate separa.  Un ajuste que resulta ser un teorema en unas casillas y un ajuste
# en otras se escribe asi, no como "una ley de 16 de 16".
#
# CONTROLES  (cada uno cuenta lo que MIDIO, no lo que recorrio)
#   G1  COHERENCIA: prod_P p^(f.e) == [O+:c] leido de la HNF del reticulo.  Dos rutas al mismo numero.
#   G2  GORENSTEIN: reparte las casillas en [O+:A] == [A:c] (teorema, sin techo) y el resto.
#       Las dos clases tienen que estar POBLADAS: si todas fueran Gorenstein el techo no se mide,
#       y si ninguna lo fuera el teorema no estaria en juego.
#   G3  LA LEY AGREGADA: [O+:A] == prod_p p^ceil(v_p/2).  Se imprime cada fallo.
#   G4  SENUELO, la version POR IDEAL (la que murio):  prod_P p^ceil(f.e/2).  Tiene que fallar.
#       Y tiene que fallar EN LAS CASILLAS PARTIDAS: se cuenta aparte cuantas casillas tienen algun
#       primo racional con >= 2 ideales en c.  Si ese contador es 0, G3 y G4 son la MISMA formula
#       sobre esta muestra y G3 no ha sido probada contra nada.  Ese es el modo de fallo de
#       `ley_ideal.sage`, que tenia UNA sola casilla partida.
#   G5  SENUELO floor.  Tiene que fallar en las no-Gorenstein.
#   G6  SENUELO sin la mitad, prod p^v_p == [O+:c].  Tiene que fallar (es el conductor, no el indice).
#
# REPL-SAFE: sin lineas en blanco dentro de bloques indentados.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run (detached):
#   docker run -d --name leyc -v "$PWD:/work" \
#     -w /work sage-normaliz:local bash -c "sage < ley_conductor.sage > ley_conductor_OUT.txt 2>&1"
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
def analiza(m, q):
    L = NumberField(psi(q), 'a')
    d = L.degree()
    a = L.gen()
    ss = simetricas(m, L)
    try:
        A = L.order(ss)
    except Exception:
        return dict(causa="orden degenerado (rango < d)")
    # O+ = Z[zeta+zeta^-1] es TEOREMA CLASICO (Washington, *Cyclotomic Fields*, Prop. 2.16), y el
    # guion lo NECESITA: el conductor de abajo se calcula con las potencias de a en los dos lados,
    # o sea dentro de Z[a].  Si Z[a] != O+ eso mide el conductor en OTRO anillo.  Hasta el 16-sep
    # esa hipotesis estaba implicita.  Ahora: se usa Z[a], y se VERIFICA donde es barato.
    # (maximal_order() en grado >= 14 reventaba la memoria -- los dos 137 del barrido ancho.)
    OL = L.order(a)
    verif = "asumida"
    if d <= GVERIF:
        if ZZ(OL.index_in(L.maximal_order())) != 1:
            return dict(causa="Z[a] NO es el orden maximal: la hipotesis del guion FALLA aqui")
        verif = "SI"
    idx = ZZ(A.index_in(OL))
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
        return dict(d=d, idx=idx, nc=ZZ(1), fac=[], verif=verif)
    # { u in Z^d : u.Mall in Z^n }: sobre Z, no sobre Zmod(N), que no es local cuando N tiene dos
    # primos distintos y devuelve un modulo sin `.basis()` utilizable justo en las casillas partidas.
    B = matrix(ZZ, N*Mall)
    n = B.ncols()
    St = block_matrix(ZZ, 2, 1, [B, -N*identity_matrix(ZZ, n)])
    K = St.left_kernel().basis_matrix()
    Lat = K[:, 0:d].hermite_form(include_zero_rows=False)
    nc = abs(Lat.determinant())
    # La factorizacion en ideales primos FUERZA el orden maximal por dentro (es lo que reventaba la
    # memoria).  G3 -- la ley -- solo necesita nc e idx, asi que en grado alto se devuelve fac=None
    # y la fila mide G3 pero NO G1/G4.  Un control que no se pudo correr se dice, no se silencia.
    # Se maximiza SOLO en el soporte de c.  El round4 GLOBAL era el que reventaba la memoria (137
    # dos veces el 16-sep); el conductor vive sobre cuatro primos y no hace falta mas.  Que el resto
    # sea maximal es el mismo teorema de Washington que ya usa G0, no una suposicion nueva.
    ps = [pp for (pp, vv) in ZZ(nc).factor()]
    try:
        L2 = NumberField(L.defining_polynomial(), 'b', maximize_at_primes=ps)
        c = L2.ideal([L2(list(v)) for v in Lat.rows()])
        fac = []
        for (P, e) in c.factor():
            fac.append((ZZ(P.residue_field().characteristic()), ZZ(P.residue_class_degree()), ZZ(e)))
    except Exception:
        return dict(d=d, idx=idx, nc=ZZ(nc), fac=None, verif=verif)
    return dict(d=d, idx=idx, nc=ZZ(nc), fac=fac, verif=verif)
import os
import sys
GRADO = ZZ(os.environ.get("GRADO_MAX", "12"))
QMAX = ZZ(os.environ.get("Q_MAX", "70"))
# Por encima de este grado NO se calcula el orden maximal ni se factoriza el conductor en ideales:
# las dos cosas fuerzan la maximalizacion y son las que mataron el barrido ancho con OOM (137).
# Esas filas miden G3 (la ley) y G5, y NO miden G1/G4/G6.  Se cuentan aparte.
GVERIF = ZZ(os.environ.get("GRADO_VERIF", "12"))
CASOS = []
MMAX = ZZ(os.environ.get("M_MAX", "5"))
# M_MIN existe porque la memoria se ACUMULA entre casillas (pila de PARI), no porque una casilla
# sea cara: el barrido murio con OOM en (5,44), que ya se habia calculado suelto sin problema.
# Trocear por m acota la acumulacion y de paso paraleliza.  Un OOM no siempre es una casilla grande.
MMIN = ZZ(os.environ.get("M_MIN", "2"))
for m in range(MMIN, MMAX + 1):
    for q in range(7, QMAX + 1):
        if q % 4 == 2:
            continue
        if euler_phi(q) // 2 > GRADO:
            continue
        if q <= 2*m:
            continue
        CASOS.append((m, q))
for t in [(3, 56), (4, 45), (4, 40)]:
    if t not in CASOS:
        CASOS.append(t)
print("=" * 118)
print("G1-G6   el indice como raiz del conductor:  [O+:A] = menor n con [O+:c] | n^2 ; y el teorema Gorenstein debajo")
print("=" * 118)
print("   barrido: m en 2..5, q <= %s, grado phi(q)/2 <= %s  -- el barrido es CAPADO y lo dice." % (QMAX, GRADO))
print("   m    q   d   [O+:A]   [O+:c]   [A:c]  Goren  primos de c (p,f,e)                 G1  ley  p-ideal  floor")
print("   " + "-" * 112)
g1med = 0
g1mal = 0
nmed = 0
gor = 0
nogor = 0
g3mal = 0
g3fallos = []
g4med = 0
g4mal = 0
g7med = 0
g7mal = 0
g7alto = 0
g7alto_vieja = 0
g7fallos = []
g5mal = 0
sin_fac = 0
g6mal = 0
partidas = 0
g4mal_part = 0
saltadas = {}
for (m, q) in CASOS:
    try:
        r = analiza(m, q)
    except Exception as ex:
        saltadas["excepcion: " + str(ex)[:40]] = saltadas.get("excepcion: " + str(ex)[:40], 0) + 1
        continue
    if "causa" in r:
        saltadas[r["causa"]] = saltadas.get(r["causa"], 0) + 1
        continue
    if r["nc"] == 1:
        saltadas["conductor trivial (indice 1)"] = saltadas.get("conductor trivial (indice 1)", 0) + 1
        continue
    idx = r["idx"]
    nc = r["nc"]
    prod_fe = ZZ(0)
    porideal = ZZ(0)
    p1pred = ZZ(0)
    fe_alto = False
    cuenta_p = {}
    if r["fac"] is None:
        sin_fac += 1
    else:
        prod_fe = ZZ(1)
        porideal = ZZ(1)
        p1pred = ZZ(1)
        for (pp, ff, ee) in r["fac"]:
            prod_fe = prod_fe * pp**(ff*ee)
            porideal = porideal * pp**ceil(ff*ee/2)
            p1pred = p1pred * pp
            if ff*ee >= 4:
                fe_alto = True
            cuenta_p[pp] = cuenta_p.get(pp, 0) + 1
        g1med += 1
        if prod_fe != nc:
            g1mal += 1
    ley = ZZ(1)
    flo = ZZ(1)
    for (pp, vv) in ZZ(nc).factor():
        ley = ley * pp**ceil(vv/2)
        flo = flo * pp**floor(vv/2)
    acx = ZZ(nc) // idx
    es_gor = (idx == acx)
    nmed += 1
    if es_gor:
        gor += 1
    else:
        nogor += 1
    if ley != idx:
        g3mal += 1
        g3fallos.append((m, q, idx, ley))
    partida = (len(cuenta_p) > 0) and (max(cuenta_p.values()) >= 2)
    if partida:
        partidas += 1
    if len(cuenta_p) > 0:
        g4med += 1
        if porideal != idx:
            g4mal += 1
            if partida:
                g4mal_part += 1
        if prod_fe != idx:
            g6mal += 1
    if flo != idx:
        g5mal += 1
    # G7: P1, pre-registrada antes de los datos.  Se mide SOLO en las no-Gorenstein:
    # en las Gorenstein [A:c] = raiz(N(c)) es teorema y cualquier formula que lo reproduzca no mide.
    if (len(cuenta_p) > 0) and (not es_gor):
        g7med += 1
        if p1pred != acx:
            g7mal += 1
            g7fallos.append((m, q, acx, p1pred))
        if fe_alto:
            g7alto += 1
            if ley != idx:
                g7alto_vieja += 1
    if r["fac"] is None:
        fs = "(sin factorizar: grado > GRADO_VERIF)"
    else:
        fs = "  ".join(["(%d,%d,%d)" % t for t in r["fac"]])
    if partida:
        fs = fs + "  <<PARTIDA"
    print("  %3d %4d %3d %8s %8s %7s  %-5s  %-40s %3s %4s %7s %6s" % (m, q, r["d"], str(idx), str(nc), str(acx), "SI" if es_gor else "no", fs, "ok" if prod_fe == nc else ("--" if r["fac"] is None else "NO"), "ok" if ley == idx else "NO", "ok" if porideal == idx else ("--" if r["fac"] is None else "no"), "ok" if flo == idx else "no"))
    sys.stdout.flush()
print("")
print("   REGLA: cada control cuenta sus MEDIDAS.  Un veredicto con 0 medidas no es un veredicto.")
print("   casillas con conductor no trivial: %d   (de %d recorridas)" % (nmed, len(CASOS)))
for (k, v) in sorted(saltadas.items()):
    print("     descartadas por %-44s %d" % (k, v))
print("   G1  prod p^(f.e) == [O+:c] :  medidas %d   fallos %d   %s" % (g1med, g1mal, ("PASA" if g1mal == 0 else "NO PASA") if g1med > 0 else "*** 0 MEDIDAS ***"))
if gor > 0 and nogor > 0:
    vg2 = "PASA -- las dos clases pobladas: el teorema cubre %d casillas y el techo hace trabajo en %d" % (gor, nogor)
else:
    vg2 = "*** una clase vacia: o el techo no se mide, o el teorema no esta en juego ***"
print("   G2  GORENSTEIN [O+:A] == [A:c] :  SI %d   no %d   %s" % (gor, nogor, vg2))
print("   G3  LEY  [O+:A] == prod p^ceil(v_p/2) :  medidas %d   fallos %d   %s" % (nmed, g3mal, ("AGUANTA" if g3mal == 0 else "*** MUERE ***") if nmed > 0 else "*** 0 MEDIDAS ***"))
for t in g3fallos:
    print("       m=%d q=%d  indice %s  la ley predice %s" % t)
print("       DONDE TIENE CONTENIDO: en las %d casillas Gorenstein la ley y G2 son la MISMA afirmacion" % gor)
print("       (N(c) cuadrado + Gorenstein => [O+:A] = raiz(N(c)) = prod p^ceil(v/2), sin eleccion).")
print("       El techo solo se mide en las %d NO Gorenstein.  Ese es el tamano real de la muestra." % nogor)
if partidas > 0 and g4mal_part > 0:
    vg4 = "PASA -- G3 esta probada donde las dos formulas DIFIEREN"
else:
    vg4 = "*** sin casilla partida que falle, G3 y G4 coinciden sobre esta muestra: G3 no se ha probado contra nada ***"
print("   G0  HIPOTESIS  Z[a] == O+ :  verificada en los %d casos de grado <= %d ; ASUMIDA (teorema de" % (g4med, GVERIF))
print("       Washington Prop. 2.16, no recomprobada) en las %d filas de grado mayor." % sin_fac)
print("   G4  SENUELO por IDEAL :  fallos %d de %d MEDIDAS ; casillas PARTIDAS %d, y falla en %d de ellas   %s" % (g4mal, g4med, partidas, g4mal_part, vg4))
if sin_fac > 0:
    print("       (%d filas SIN factorizacion de ideales por grado > %d: miden G3/G5 y no G1/G4/G6.)" % (sin_fac, GVERIF))
print("   G5  SENUELO floor :  fallos %d de %d   %s" % (g5mal, nmed, ("PASA -- el techo esta identificado" if g5mal > 0 else "FALLA -- floor y ceil no se distinguen")))
if g7med == 0:
    vg7 = "*** 0 MEDIDAS: sin casilla no-Gorenstein factorizada, P1 no se ha probado ***"
elif g7mal == 0:
    vg7 = "AGUANTA sobre %d casillas NO Gorenstein" % g7med
else:
    vg7 = "*** MUERE ***"
print("   G7  P1  [A:c] == prod_{P|c} p  (pre-registrada) :  medidas %d   fallos %d   %s" % (g7med, g7mal, vg7))
for t in g7fallos:
    print("       m=%d q=%d  [A:c]=%s  P1 predice %s" % t)
print("       Y donde P1 y la ley VIEJA difieren, que es f.e >= 4:  %d casillas ; la vieja falla en %d." % (g7alto, g7alto_vieja))
if g7alto == 0:
    print("       (0 casillas con f.e >= 4: sobre esta muestra P1 y la ley vieja son la MISMA formula")
    print("        y G7 no la ha sustituido -- es exactamente el modo de fallo del barrido de anoche.)")
print("   G6  SENUELO sin la mitad :  fallos %d de %d MEDIDAS   %s" % (g6mal, g4med,("PASA" if g6mal > 0 else "FALLA -- la mitad no hace nada")))
