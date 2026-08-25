# -*- coding: utf-8 -*-
# EL LADO DE LAPLACE, ANTES DE L5 --- y la condicion (4) de 41_IN, bien posada
# 19 de agosto de 2026, noche.  Continuacion de pI_involucion_toggle.py.
#
# ############################################################################################
# LO PRIMERO, Y ES UNA CORRECCION DE METODO: ESTO YA ESTA EN EL PAPER I.
#
# La identidad que yo estaba re-derivando es el LEMA L3 del propio Paper I (paper/orbit_pair.tex,
# \label{lem:L3}), literal:
#
#     Phi_t(lambda;z) = (-1)^{t + C(N+1,2)} / [(z^t-1)(z^{-t}-1)(z-z^{-1})]
#                       . sum_{U in Ucal} (-1)^{j1+j2+inv(b_{S_U})} f(beta_{j1} - beta_{j2})
#
# con las columnas indexadas DESDE 1.  Y el parrafo inmediatamente posterior al lema ya demuestra
# --- no conjetura --- que |Ucal| es 0, 3 o 4, que es exactamente lo que yo habia escrito como
# «prediccion escrita antes de medir».  No era una prediccion: era un teorema nuestro, publicado.
#
# La primera version de esta gate se lo salto y calculo los signos a ojo.  Le faltaba el factor
# global (-1)^{t+C(N+1,2)} y le sobraba un +1 en el exponente, o sea estaba globalmente cambiada de
# signo: G1 daba 40/329 y 168/791, que son EXACTAMENTE las lambda con los dos lados nulos.  El fallo
# fue del tipo que la casa ya tiene fichado: re-derivar en vez de abrir el .tex.
#
# QUE QUEDA ENTONCES DE ESTA GATE.  Deja de ser un descubrimiento y pasa a ser dos cosas utiles:
#   (i)  CALIBRACION de un lema publicado con un instrumento nuevo --- que es lo que la Regla 2 de
#        ERRATA_PENDIENTE.md obliga a hacer el dia que el instrumento nace;
#   (ii) el numero que la condicion (4) de 41_IN puede exigir de verdad, que es el objeto del que
#        se habla mas abajo y NO es ||Phi||_1.
# ############################################################################################
#
# LA CORRECCION QUE ORIGINA ESTA GATE.  Al cerrar pI_involucion_toggle.py escribi que el paso
# siguiente era «comparar ||Phi||_1 contra la norma L1 del lado de Laplace».  ESO ESTA MAL POSADO:
# el lado de Laplace es un NUMERADOR, y Phi es ese numerador dividido por el Vandermonde.  Comparar
# sus normas L1 es comparar dos objetos distintos.  La coincidencia que señale ---maximo 30 contra
# 4 transversales x 8 monomios = 32--- era coincidencia de talla, y avise de que lo era.
#
# LA CUENTA QUE LO ARREGLA, y es lo que hace la gate posible sin tocar una sola raiz de la unidad.
# Con x = (1, w, ..., w^{t-1}, z, z^{-1}), N = t+2, beta = beta(lambda; N) y Phi = a_beta / a_rho:
#
#     a_rho = prod_{i<j<t}(w^i - w^j) . prod_{i<t}(w^i - z) . prod_{i<t}(w^i - z^{-1}) . (z - z^{-1})
#           = Delta_orbita . (-1)^t (z^t - 1) . (-1)^t (z^{-t} - 1) . (z - z^{-1})
#           = Delta_orbita . (z^t - 1)(z^{-t} - 1)(z - z^{-1})
#
# y en la expansion de Laplace de a_beta por las DOS filas libres, el menor complementario t x t es
# el Vandermonde de {w^{beta_j} : j en S}.  Cuando los beta_j de S tocan las t clases de residuo,
# ese conjunto ES el de todas las raices t-esimas, luego el menor vale +- Delta_orbita y AL DIVIDIR
# SE VA.  Queda una identidad con coeficientes enteros:
#
#     Phi_t(lambda;z) . (z^t - 1)(z^{-t} - 1)(z - z^{-1})  =  sum_{ {j1<j2} admisibles } eps . f(beta_j1 - beta_j2)
#
# ---y el eps es el del Lema L3, no uno mio: eps = (-1)^{t+C(N+1,2)} . (-1)^{j1+j2+inv(b_S)}, con
#    las columnas desde 1.  (La primera version de esta gate escribio (-1)^{1+j1+j2}.sgn(pi), que es
#    el mismo objeto globalmente cambiado de signo.  Ver la caja de arriba.)
#
# ESO es «(A) antes de L5»: una suma de 3 o 4 terminos con coeficientes +-1.  Los signos NO se
# ajustan a posteriori --- se toman del Lema L3 y se COMPRUEBAN.
#
# LO QUE ESTO LE HACE A SU CONDICION (4).  Su (4) pide que los puntos fijos de la involucion sean
# las transversales de Laplace.  Con la cuenta de arriba, eso NO puede ser un enunciado sobre Phi:
#     - una involucion sobre J_lambda que cumpla (2) y (3) deja ||Phi||_1 puntos fijos (hasta 30);
#     - el lado de Laplace tiene a lo sumo 3 o 4 terminos, o sea a lo sumo 8 monomios.
# Los dos numeros viven en objetos distintos, separados por el factor (z^t-1)(z^{-t}-1)(z-z^{-1}).
# Luego (4) solo puede leerse de una de estas dos maneras, y la gate mide cual:
#     (L-a) la involucion vive sobre J_lambda y la proyeccion a las transversales es MUCHOS-A-UNO;
#     (L-b) la involucion vive sobre el indice del NUMERADOR, no sobre J_lambda.
#
# QUE SE MIDE
#   G1  FATAL: el Lema L3 del Paper I, lado a lado, con SUS signos --- no con los mios.  Es una
#       calibracion de un lema publicado, no un hallazgo.
#   G2  el numero de transversales admisibles por lambda.  El Paper I ya DEMUESTRA que es 0, 3 o 4
#       (parrafo posterior a lem:L3); aqui solo se comprueba que el instrumento lo reproduce.
#   G3  ||lado de Laplace||_1, que es el numero que su condicion (4) puede exigir de verdad.
#   G4  el contraste con ||Phi||_1 --- los dos objetos, uno al lado del otro.
#   D1  SENUELO: incluir TODOS los pares {j1<j2}, sin la condicion de residuos.  Tiene que romper.
#   D2  SENUELO: todos los signos +1.  Tiene que romper.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_laplace_antes_de_L5.py > pI_laplace_antes_de_L5_OUT.txt 2>&1

import itertools
import json
import sys
from collections import Counter

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


def f(d):
    if d == 0:
        return {}
    return {d: 1, -d: -1}


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


def sgn_perm(seq):
    """Signo de la permutacion que ordena seq (valores distintos) crecientemente."""
    s, L = 1, list(seq)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            if L[i] > L[j]:
                s = -s
    return s


def norma1(p):
    return sum(abs(v) for v in p.values())


print("=" * 100)
print("EL LADO DE LAPLACE, ANTES DE L5 --- Phi . (z^t-1)(z^-t-1)(z-z^-1) = suma sobre transversales")
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

    # D = (z^t - 1)(z^-t - 1)(z - z^-1)
    D = pmul(pmul({T: 1, 0: -1}, {-T: 1, 0: -1}), {1: 1, -1: -1})

    lams = [l for l in particiones(LMAX, N) if l]
    g1_ok, g1_bad = 0, []
    g2 = Counter()
    g3 = Counter()
    g4_pares = []
    d1_rompe, d2_rompe = 0, 0
    d1_n, d2_n = 0, 0

    for lam in lams:
        l = list(lam) + [0] * N
        beta = [l[i] + N - 1 - i for i in range(N)]
        # transversales admisibles: quitar {j1,j2} deja los t restantes con residuos todos distintos
        trans = []
        todos = []
        for (j1, j2) in itertools.combinations(range(N), 2):
            S = [j for j in range(N) if j not in (j1, j2)]
            res = [beta[j] % T for j in S]
            todos.append((j1, j2, S, res))
            if len(set(res)) == T:
                trans.append((j1, j2, S, res))
        # Lema L3 del Paper I, con SUS signos.  Columnas indexadas desde 1 (aqui j es 0-indexado,
        # luego (j1+1)+(j2+1) = j1+j2 en paridad).  Factor global (-1)^{t + C(N+1,2)}.
        glob = (-1) ** (T + (N + 1) * N // 2)
        lado = {}
        for (j1, j2, S, res) in trans:
            eps = glob * (-1) ** (j1 + j2) * sgn_perm(res)
            lado = padd(lado, {k: eps * v for k, v in f(beta[j1] - beta[j2]).items()})
        phi = schur_directo(lam)
        izq = pmul(phi, D)
        if izq == lado:
            g1_ok += 1
        else:
            g1_bad.append(tuple(lam))
        g2[len(trans)] += 1
        g3[norma1(lado)] += 1
        g4_pares.append((norma1(phi), norma1(lado)))
        # D1: sin la condicion de residuos
        d1 = {}
        for (j1, j2, S, res) in todos:
            eps = glob * (-1) ** (j1 + j2) * (sgn_perm(res) if len(set(res)) == T else 1)
            d1 = padd(d1, {k: eps * v for k, v in f(beta[j1] - beta[j2]).items()})
        d1_n += 1
        if d1 != izq:
            d1_rompe += 1
        # D2: todos los signos +1
        d2 = {}
        for (j1, j2, S, res) in trans:
            d2 = padd(d2, f(beta[j1] - beta[j2]))
        d2_n += 1
        if d2 != izq:
            d2_rompe += 1

    n = len(lams)
    print("  t=%d  (%d lambda, lambda_1 <= %d, %d columnas de beta)" % (T, n, LMAX, N))
    print("    G1  FATAL  Phi.(z^t-1)(z^-t-1)(z-z^-1) == suma de transversales con signo Laplace")
    print("        %d de %d   %s" % (g1_ok, n, "PASA" if not g1_bad else
                                     "*** FALLA en %s ***" % g1_bad[:5]))
    print("    G2  numero de transversales admisibles : %s" % dict(sorted(g2.items())))
    print("        el Paper I ya lo demuestra (0, 3 o 4); el instrumento lo reproduce --> %s"
          % ("SE CUMPLE" if set(g2) <= {0, 3, 4} else "*** FALSA, aparece %s ***"
             % sorted(set(g2) - {0, 3, 4})))
    print("    G3  ||lado de Laplace||_1 : %s     maximo %d" % (dict(sorted(g3.items())), max(g3)))
    print("    G4  el contraste, ||Phi||_1 contra ||Laplace||_1 :")
    print("        ||Phi||_1   de %d a %d" % (min(p for p, _ in g4_pares), max(p for p, _ in g4_pares)))
    print("        ||Lapl||_1  de %d a %d" % (min(q for _, q in g4_pares), max(q for _, q in g4_pares)))
    print("        lambda con ||Phi||_1 > ||Lapl||_1 : %d de %d"
          % (sum(1 for p, q in g4_pares if p > q), n))
    print("    D1  SENUELO sin la condicion de residuos : rompe en %d de %d  --> %s"
          % (d1_rompe, d1_n, "muerde" if d1_rompe > 0 else "*** NO MUERDE ***"))
    print("    D2  SENUELO con todos los signos +1     : rompe en %d de %d  --> %s"
          % (d2_rompe, d2_n, "muerde" if d2_rompe > 0 else "*** NO MUERDE ***"))
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = {"n": int(n), "G1": [int(g1_ok), len(g1_bad)],
                           "G1_fallos": g1_bad[:20],
                           "G2": {str(k): int(v) for k, v in g2.items()},
                           "G3": {str(k): int(v) for k, v in g3.items()},
                           "D1_rompe": int(d1_rompe), "D2_rompe": int(d2_rompe)}

json.dump(resumen, open("pI_laplace_antes_de_L5_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
