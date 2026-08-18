# -*- coding: utf-8 -*-
# Newt(N_delta) POR SUMAS DE MINKOWSKI, ENTERO.   16 de agosto de 2026.
#
# prop:newtden prueba UN vertice de Newt(N_delta), el dominante.  La vuelta 26 pedia el politopo
# entero por sumas de Minkowski, y marcaba el teorema como el mas barato de la hoja de ruta.
#
# LA CUENTA, antes de programar.  N_delta = +- prod_{i<j} (x_i - x_j) sobre el alfabeto
# (1, zeta, ..., zeta^{t-1}, z_1^{+-1}, ..., z_r^{+-1}).  Cada factor es un binomio, luego su Newton
# en z es un SEGMENTO, y por Ostrowski el total es la suma de Minkowski.  Agrupando:
#
#   factores (zeta^a - zeta^b)        : punto {0}                       -- no aportan
#   por cada k:  t.[0,e_k] + t.[-e_k,0] + [-e_k,e_k]  =  [-(t+1)e_k, (t+1)e_k]
#   por cada k<l: los cuatro emparejamientos (z_k^{+-1} - z_l^{+-1}) suman el zonotopo generado
#                 por  e_k - e_l  y  e_k + e_l,  cada uno con radio 1
#
# O sea:  Newt(N_delta) es el ZONOTOPO generado por las raices positivas de C_r, con las cortas
# e_k +- e_l de radio 1 y las largas de radio t+1.  Y entonces sus vertices son la orbita de W(C_r)
# del vertice dominante de prop:newtden, o sea  2^r . r!  vertices.
#
# LO QUE SE MIDE
#   V1  FATAL: el zonotopo == el politopo de orbita  conv( W(C_r) . (N-1, N-3, ..., N-2r+1) ).
#   V2  numero de vertices == 2^r r!.
#   V3  FATAL donde se pueda expandir: el casco convexo del SOPORTE de N_delta expandido == V1.
#       (Ostrowski dice que no hay encogimiento por cancelacion; aqui se comprueba sobre un cuerpo
#       ciclotomico, que es donde uno se lo creeria menos.)
#   V4  el vertice dominante coincide con prop:newtden.
#
# SENUELOS
#   D1  radio t en vez de t+1 en las largas.  Debe fallar.
#   D2  las raices de D_r (sin las largas) en vez de las de C_r.  Debe fallar.
#   D3  las cortas con radio 2 en vez de 1.  Debe fallar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage newt_zonotope.sage

import json
import sys
import itertools

CASOS = [(t, r) for t in range(2, 13) for r in range(1, 5)]


def seg(v):
    """el segmento [-v, v] como politopo."""
    return Polyhedron(vertices=[[-x for x in v], list(v)])


def zonotopo(t, r, radio_larga=None, cortas=True, radio_corta=1, tipo="C"):
    if radio_larga is None:
        radio_larga = t + 1
    P = Polyhedron(vertices=[[0] * r])
    for k in range(r):
        v = [0] * r
        v[k] = radio_larga
        if tipo == "C":
            P = P + seg(v)
    if cortas:
        for k in range(r):
            for l in range(k + 1, r):
                for s in (1, -1):
                    v = [0] * r
                    v[k] = radio_corta
                    v[l] = s * radio_corta
                    P = P + seg(v)
    return P


def orbita_C(v):
    r = len(v)
    pts = set()
    for perm in itertools.permutations(range(r)):
        for eps in itertools.product((1, -1), repeat=r):
            pts.add(tuple(eps[i] * v[perm[i]] for i in range(r)))
    return Polyhedron(vertices=[list(p) for p in pts])


def newt_expandido(t, r):
    """el casco del soporte de N_delta expandido de verdad.  Solo para N pequeno."""
    N = t + 2 * r
    K = CyclotomicField(t) if t > 2 else QQ
    zeta = K.gen() if t > 2 else K(-1)
    L = LaurentPolynomialRing(K, r, 'z')
    zs = L.gens()
    alfabeto = [L(K(zeta) ** a) for a in range(t)]
    for g in zs:
        alfabeto += [g, g ** -1]
    P = L(1)
    for i in range(len(alfabeto)):
        for j in range(i + 1, len(alfabeto)):
            P *= (alfabeto[i] - alfabeto[j])
    if P == 0:
        return None
    pts = []
    for e, c in zip(P.exponents(), P.coefficients()):
        if c == 0:
            continue
        pts.append([int(v) for v in (e if hasattr(e, '__iter__') else (e,))])
    return Polyhedron(vertices=pts)


print("=" * 108)
print("Newt(N_delta) COMO ZONOTOPO DE LAS RAICES DE C_r")
print("=" * 108)
print("")
print("   t   r    N   V1 zonotopo==orbita   V2 vertices   2^r r!   V3 expandido   V4 top")
print("   " + "-" * 92)
sys.stdout.flush()

RES = []
d1 = d2 = d3 = nd = 0
for (t, r) in CASOS:
    N = t + 2 * r
    Z = zonotopo(t, r)
    top = [N - 1 - 2 * k for k in range(r)]
    O = orbita_C(top)
    v1 = (Z == O)
    nv = Z.n_vertices()
    esperado = 2 ** r * factorial(r)
    v2 = (nv == esperado)
    v4 = tuple(max(Z.vertices_list(), key=lambda p: sum((r - i) * p[i] for i in range(r)))) == tuple(top)
    v3 = "---"
    if N <= 8:
        E = newt_expandido(t, r)
        v3 = "si" if (E is not None and E == Z) else "NO"
    # senuelos
    nd += 1
    d1 += 1 if zonotopo(t, r, radio_larga=t) == O else 0
    d2 += 1 if zonotopo(t, r, tipo="D") == O else 0
    d3 += 1 if zonotopo(t, r, radio_corta=2) == O else 0
    print("  %2d  %2d  %3d        %-5s          %5d      %5d      %-5s        %-5s"
          % (t, r, N, str(v1), nv, esperado, v3, str(v4)))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "N": int(N), "V1": bool(v1), "vertices": int(nv),
                "esperado": int(esperado), "V2": bool(v2), "V3": v3, "V4": bool(v4)})

print("")
print("  V1  zonotopo == politopo de orbita : %d de %d" % (sum(1 for x in RES if x["V1"]), len(RES)))
print("  V2  vertices == 2^r r!             : %d de %d" % (sum(1 for x in RES if x["V2"]), len(RES)))
print("  V3  contra el expandido            : %d de %d  (los casos con N <= 8)"
      % (sum(1 for x in RES if x["V3"] == "si"), sum(1 for x in RES if x["V3"] != "---")))
print("  V4  el vertice dominante           : %d de %d" % (sum(1 for x in RES if x["V4"]), len(RES)))
print("  D1  SENUELO radio t en las largas  : %d de %d  (debe ser 0)" % (d1, nd))
print("  D2  SENUELO sin las largas (D_r)   : %d de %d  (debe ser 0)" % (d2, nd))
print("  D3  SENUELO cortas de radio 2      : %d de %d  (debe ser 0)" % (d3, nd))

json.dump({"casos": RES, "senuelos": {"D1": int(d1), "D2": int(d2), "D3": int(d3), "n": int(nd)}},
          open("newt_zonotope_DUMP.json", "w"), indent=1)
print("")
print("=" * 108)
print("DONE")
