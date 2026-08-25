# -*- coding: utf-8 -*-
# EL EMPAREJAMIENTO AL NIVEL (mu,nu), DENTRO DE CADA EXPONENTE DE z --- su plan de manana, punto 1
# 19 de agosto de 2026, noche.  Continuacion de pI_involucion_refinada.py.
#
# LO QUE ESTA GATE MIDE, Y LO QUE NO.  No construye la involucion: mide su OBSTRUCCION MINIMA.  Dentro
# de un exponente e fijo hay n+(e) terminos con sgn_t = +1 y n-(e) con -1; CUALQUIER involucion que
# cumpla (2) y (3) deja al menos |n+(e) - n-(e)| puntos fijos ahi, y existe una que deja exactamente
# eso.  Luego
#
#     F(lambda) = sum_e |n+(e) - n-(e)|
#
# es el numero de puntos fijos MINIMO POSIBLE, alcanzable, y no depende de que toggle se elija.  Si
# F es grande o crece con lambda, la condicion (4) es inalcanzable y NINGUN toggle --- ni el de primer
# cruce ni otro --- puede cumplir su lista.  Si F es pequeno y estable, la condicion (4) tiene sitio y
# el trabajo siguiente (construir el toggle canonico) esta justificado.
#
# LA PREDICCION, ESCRITA ANTES DE MEDIR.  Su condicion (4) pide que los supervivientes sean las
# transversales de Laplace, que son 3 o 4.  El producto de tres f = z^a - z^{-a} tiene 8 monomios.
# Luego:
#     (P1)  F(lambda) esta ACOTADO --- no crece con |lambda| --- y su cota esta en el entorno de 8.
#     (P2)  F(lambda) << el residuo al nivel de mu solo (hasta 14 en t=2 y 24 en t=3).
# FALSADOR DECLARADO: si F crece con lambda, (P1) muere y con ella la condicion (4) tal como esta
# escrita; habria que reescribirla o admitir puntos fijos en numero variable.
#
# ############################################################################################
# CORRECCION, ESCRITA AL LEER LA PRIMERA SALIDA Y ANTES DE CONTARLA COMO RESULTADO.
#
# G1 SE CAYO EN TAUTOLOGIA.  El coeficiente de z^e en Phi ES n+(e) - n-(e).  Luego
#
#         F(lambda) = sum_e |n+(e) - n-(e)| = sum_e |coef_e(Phi)| = || Phi ||_1
#
# identicamente, para toda lambda y todo t.  O sea F NO ES UNA MEDIDA SOBRE LA INVOLUCION: es la
# norma L1 del propio resultado, calculable sin mirar J_lambda.  La primera salida lo enseñaba a la
# cara y yo no lo vi: G3 daba «F = numero de monomios» en 109 de 329, que es exactamente 329 menos
# los 220 con algun coeficiente distinto de +-1.
#
# ESO NO DEJA LA GATE SIN CONTENIDO --- LE CAMBIA EL ENUNCIADO, y a mejor:
#
#     El NUMERO de puntos fijos no es una eleccion de diseno.  Cualquier involucion que cumpla sus
#     condiciones (2) y (3) deja AL MENOS ||Phi||_1 puntos fijos, sea cual sea el toggle.  Luego su
#     condicion (4) no puede leerse como «deja pocos puntos fijos»:
#     es enteramente una afirmacion sobre CUALES sobreviven, y es falsable con un solo numero ---
#     el lado de Laplace, antes de L5, tiene que tener la MISMA norma L1.
#
# Por eso: (a) G1 se reetiqueta como IDENTIDAD y se VERIFICA como tal (F == ||Phi||_1 en todas), que
# es un control que si puede fallar --- fallaria si mi algebra estuviera mal; (b) el senuelo D1 que
# escribi («sin signos F seria |J|») NO PUEDE FALLAR y se retira, no se maquilla; (c) la pregunta
# viva pasa a ser G4, la norma que el lado de Laplace debe reproducir, y esa es la gate siguiente.
# (P1) y (P2) quedan como estaban escritas, para que se vea que (P1) era una prediccion mal planteada.
#
# SEGUNDA CORRECCION, 43_IN, Y ES SUYA: LO DE ARRIBA ES UNA DESIGUALDAD, NO UNA IGUALDAD.
# F(lambda) = ||Phi||_1 es el MINIMO sobre involuciones, y es alcanzable --- eso sigue siendo cierto
# y es lo que esta gate calcula.  Pero UNA involucion arbitraria que cumpla (2) y (3) puede dejar
# MAS: dentro de un exponente e, si F_e+ y F_e- son los fijos de cada signo, los 2-ciclos dan
# n+(e) - n-(e) = F_e+ - F_e-, luego F_e+ + F_e- >= |n+(e) - n-(e)|, con igualdad SOLO si todos los
# fijos de ese peso tienen el mismo signo (cancelacion maximal en cada weight space).  Nada impide
# que queden un + y un - los dos fijos.  Lo correcto es
#
#         |Fix| >= ||Phi||_1
#
# y la conclusion de 42_OUT §3 sobrevive igual, con el argumento mas limpio: con ||Phi||_1 = 30,
# cualquier involucion monomial-preserving sobre J_lambda necesita AL MENOS 30 puntos fijos, y
# Laplace tiene tres o cuatro transversales.
# ############################################################################################
#
# QUE SE MIDE
#   G0  CALIBRACION: la identidad refinada contra el valor directo, 329/329 y 791/791 (pI_involucion
#       _refinada.py).  Si no reproduce, no estamos midiendo lo mismo.
#   G1  IDENTIDAD, verificada y no supuesta: F(lambda) == ||Phi||_1 en todas las lambda.  Si fallara,
#       mi lectura del coeficiente estaria mal y todo lo demas caeria.
#   G2  la compresion: |J_lambda| (terminos totales) contra F.  Cuanto cancela el refinamiento.
#   G3  .son todos los coeficientes de Phi iguales a +-1?  Si lo son, F = numero de monomios de Phi.
#   G4  EL NUMERO QUE HAY QUE DEVOLVERLE: el histograma y el MAXIMO de ||Phi||_1, que es lo que el
#       lado de Laplace tiene que reproducir para que su condicion (4) sea posible.
#   D2  SENUELO DE ALCANCE: emparejar SIN respetar el exponente (o sea tirando su condicion (2)).
#       Tiene que dar F menor --- si diera lo mismo, la condicion (2) no estaria haciendo nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_involucion_toggle.py > pI_involucion_toggle_OUT.txt 2>&1

import itertools
import json
import sys
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CASOS = [2, 3]
LMAX = 7


def padd(a, b):
    o = dict(a)
    for k in b:
        o[k] = o.get(k, 0) + b[k]
        if o[k] == 0:
            del o[k]
    return o


def pneg(a):
    return {k: -v for k, v in a.items()}


def pmul(a, b):
    o = {}
    for k1 in a:
        for k2 in b:
            k = k1 + k2
            o[k] = o.get(k, 0) + a[k1] * b[k2]
    return {k: v for k, v in o.items() if v != 0}


def chi(r):
    if r < 0:
        return {}
    return {r - 2 * i: 1 for i in range(r + 1)}


def det(M):
    n = len(M)
    if n == 0:
        return {0: 1}
    if n == 1:
        return dict(M[0][0])
    acc = {}
    for j in range(n):
        if not M[0][j]:
            continue
        men = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        tm = pmul(M[0][j], det(men))
        acc = padd(acc, tm if j % 2 == 0 else pneg(tm))
    return acc


def beta_set(lam, n):
    l = list(lam) + [0] * (n - len(lam))
    return [l[i] + n - 1 - i for i in range(n)]


def part_from_beta(B):
    S = sorted(B, reverse=True)
    n = len(S)
    return tuple(x for x in [S[i] - (n - 1 - i) for i in range(n)] if x > 0)


def core_quot_sgn(lam, n, t):
    B = set(beta_set(lam, n))
    sg = 1
    seguir = True
    while seguir:
        seguir = False
        for x in sorted(B, reverse=True):
            if x - t >= 0 and (x - t) not in B:
                sg *= (-1) ** len([y for y in B if x - t < y < x])
                B.discard(x)
                B.add(x - t)
                seguir = True
                break
    core = part_from_beta(B)
    orig = sorted(beta_set(lam, n), reverse=True)
    quot = []
    for r in range(t):
        xs = sorted([x for x in orig if x % t == r], reverse=True)
        m = len(xs)
        quot.append(tuple(x for x in [(xs[k] - r) // t - (m - 1 - k) for k in range(m)] if x > 0))
    return core, quot, sg


def sub(mu, lam):
    m = list(mu) + [0] * 10
    l = list(lam) + [0] * 10
    return all(m[i] <= l[i] for i in range(10))


def hstrip(sup, inf):
    a = list(sup) + [0] * 10
    b = list(inf) + [0] * 10
    if any(b[i] > a[i] for i in range(10)):
        return False
    return all(a[i + 1] <= b[i] for i in range(9))


def particiones(maxpart, maxlen):
    out = [()]

    def rec(pref, resto, tope):
        for p in range(min(tope, maxpart), 0, -1):
            nu = pref + [p]
            out.append(tuple(nu))
            if resto > 1:
                rec(nu, resto - 1, p)

    rec([], maxlen, maxpart)
    return sorted(set(out))


print("=" * 100)
print("EL EMPAREJAMIENTO AL NIVEL (mu,nu) DENTRO DE CADA EXPONENTE --- .cuantos puntos fijos MINIMOS?")
print("=" * 100)
print("")

resumen = {}

for T in CASOS:
    N = T + 2

    def hA(k):
        if k < 0:
            return {}
        acc, m = {}, 0
        while m * T <= k:
            acc = padd(acc, chi(k - m * T))
            m += 1
        return acc

    def schur_directo(lam):
        n = max(1, len(lam))
        l = list(lam) + [0] * n
        return det([[hA(l[i] - (i + 1) + (j + 1)) for j in range(n)] for i in range(n)])

    lams = [l for l in particiones(LMAX, N) if l]
    g0_ok, g0_bad = 0, []
    Fhist = Counter()
    F_por_tam = defaultdict(list)
    Jtot, Ftot = 0, 0
    coef_no_unidad = 0
    F_igual_monomios = 0
    d1_igual, d1_no = 0, 0
    d2_menor, d2_igual = 0, 0
    peor = []

    for lam in lams:
        mus = []
        for mu in particiones(lam[0], len(lam) + 1):
            if not sub(mu, lam):
                continue
            core, quot, sg = core_quot_sgn(mu, N, T)
            if len(core) != 0 or any(len(q) > 1 for q in quot):
                continue
            mus.append((mu, sg))
        porexp = defaultdict(lambda: [0, 0])
        npos, nneg = 0, 0
        for (mu, sg) in mus:
            for nu in particiones(lam[0], len(lam) + 1):
                if not (sub(mu, nu) and sub(nu, lam)):
                    continue
                if not (hstrip(nu, mu) and hstrip(lam, nu)):
                    continue
                e = 2 * sum(nu) - sum(mu) - sum(lam)
                porexp[e][0 if sg > 0 else 1] += 1
                if sg > 0:
                    npos += 1
                else:
                    nneg += 1
        refi = {}
        for e, (a, b) in porexp.items():
            if a - b:
                refi[e] = a - b
        directo = schur_directo(lam)
        if refi == directo:
            g0_ok += 1
        else:
            g0_bad.append(lam)
        F = sum(abs(a - b) for (a, b) in porexp.values())
        J = npos + nneg
        Fhist[F] += 1
        F_por_tam[sum(lam)].append(F)
        Jtot += J
        Ftot += F
        if any(abs(v) != 1 for v in directo.values()):
            coef_no_unidad += 1
        if F == len(directo):
            F_igual_monomios += 1
        # G1: la identidad F == ||Phi||_1, VERIFICADA (no supuesta)
        if F == sum(abs(v) for v in directo.values()):
            d1_igual += 1
        else:
            d1_no += 1
        # D2: emparejar ignorando el exponente
        F2 = abs(npos - nneg)
        if F2 < F:
            d2_menor += 1
        elif F2 == F:
            d2_igual += 1
        if len(peor) < 4 and F >= 10:
            peor.append((tuple(lam), J, F))

    n = len(lams)
    tams = sorted(F_por_tam)
    print("  t=%d  (%d lambda, lambda_1 <= %d, hasta %d filas)" % (T, n, LMAX, N))
    print("    G0  CALIBRACION identidad refinada contra el valor directo : %d de %d   %s"
          % (g0_ok, n, "PASA" if not g0_bad else "*** FALLA en %s ***" % g0_bad[:4]))
    print("    G1  F(lambda) == ||Phi||_1, el MINIMO alcanzable (verificado) : %d de %d   %s"
          % (d1_igual, n, "PASA" if d1_no == 0 else "*** FALLA en %d ***" % d1_no))
    print("        --> TODA involucion que cumpla (2) y (3) deja AL MENOS ||Phi||_1 puntos fijos.")
    print("            (43_IN: es una DESIGUALDAD; la igualdad exige cancelacion maximal en cada")
    print("             weight space.  La conclusion es la misma y el argumento mas limpio.)")
    print("            Su condicion (4) es por tanto SOLO sobre CUALES sobreviven.")
    print("    G4  la norma que el lado de Laplace tiene que reproducir")
    print("        histograma de ||Phi||_1 : %s" % dict(sorted(Fhist.items())))
    print("        maximo : %d     media : %.2f" % (max(Fhist), Ftot / n))
    print("        ||Phi||_1 por |lambda| (min..max) :")
    for s in tams:
        v = F_por_tam[s]
        print("           |lambda| = %2d : n=%4d   F de %d a %d" % (s, len(v), min(v), max(v)))
    print("    G2  compresion : %d terminos (mu,nu) en total --> %d puntos fijos  (%.1f%%)"
          % (Jtot, Ftot, 100.0 * Ftot / Jtot))
    print("    G3  lambda con algun coeficiente de Phi distinto de +-1 : %d de %d" % (coef_no_unidad, n))
    print("        F = numero de monomios de Phi : %d de %d" % (F_igual_monomios, n))
    print("    (el senuelo D1 que habia escrito --- «sin signos F seria |J|» --- NO PUEDE FALLAR")
    print("     y se ha retirado, no maquillado.  Ver la nota de correccion en la cabecera.)")
    print("    D2  SENUELO DE ALCANCE (emparejar sin respetar el exponente):")
    print("        da MENOS puntos fijos en %d de %d lambda ; igual en %d" % (d2_menor, n, d2_igual))
    print("        --> %s" % ("la condicion (2) HACE TRABAJO" if d2_menor > 0 else
                              "*** la condicion (2) no aporta nada ***"))
    if peor:
        print("    peores casos (lambda, |J|, F) : %s" % peor)
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = {"n": int(n), "G0": [int(g0_ok), len(g0_bad)],
                           "Fhist": {str(k): int(v) for k, v in Fhist.items()},
                           "Fmax": int(max(Fhist)), "Jtot": int(Jtot), "Ftot": int(Ftot),
                           "F_por_tam": {str(s): [min(v), max(v)] for s, v in F_por_tam.items()},
                           "coef_no_unidad": int(coef_no_unidad),
                           "F_igual_monomios": int(F_igual_monomios),
                           "D2_menor": int(d2_menor)}

json.dump(resumen, open("pI_involucion_toggle_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
