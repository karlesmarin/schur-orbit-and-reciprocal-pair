# -*- coding: utf-8 -*-
# BARRIDO ANCHO DE LA UNIMODULARIDAD.   16 de agosto de 2026.
#
# _probe_determinante.py establecio  c(Lambda,X) = eps_t . det M(Lambda,X)  en 87 652 de 87 652, y
# conjeturo que M es TOTALMENTE UNIMODULAR -- lo que implicaria (L1).  Pero esa medida vive en cuatro
# pares (t,r) con t <= 7 y r <= 3, y la unimodularidad total estaba MUESTREADA a 1000 matrices por
# poblacion, no agotada.  Es la evidencia mas fina del paper sosteniendo el enunciado mas pesado.
#
# Aqui se ensancha: mas t, mas r, cajas mayores, y la unimodularidad AGOTADA sobre cada matriz viva.
#
# EL PREFILTRO, y por que es licito.  Recorrer toda la caja de X y calcular un determinante R'xR' en
# cada punto es lo que hacia lento el barrido anterior.  Pero si una columna de M es ENTERAMENTE
# NULA, det M = 0 sin mas cuenta -- y eso pasa en la inmensa mayoria de los X.  El prefiltro
# descarta esos, y se comprueba (P0) que en una muestra de descartados el determinante es
# efectivamente 0, para que el atajo no se crea sin control.
#
# CONTROLES
#   D0 (FATAL)  c(Lambda,X) == eps_t . det M(Lambda,X) sobre los X VIVOS (los que sobreviven al
#               prefiltro).  Puntuar sobre toda la caja infla el marcador: la inmensa mayoria son
#               0 == 0 y cualquier matriz rota acierta ahi.
#   P0          en una muestra de X descartados por el prefiltro, det M == 0 y c == 0.
#   TU (FATAL)  unimodularidad total AGOTADA sobre cada matriz viva -- todos los menores.
#   D2          reparto de |det M|;  D3  reparto de las entradas de M.
#   D1          el signo global es eps_t en cada rama.
#   S1 senuelo  progresion par;  S2 senuelo  sin los signos del plegado.  Puntuados solo donde hay
#               algo que puntuar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python unimodularidad_barrido.py

import itertools
import json
import random
import sys
from collections import Counter

from divided_differences import plegar, nu_de, nu_extendida, eps_t, enderezar_D

# (t, r, cota).  Los cuatro primeros son la poblacion vieja; el resto es lo nuevo.
CASOS = [(3, 2, 6), (5, 2, 4), (7, 2, 3), (3, 3, 4),
         (3, 2, 12), (5, 2, 7), (7, 2, 5), (9, 2, 3), (11, 2, 2),
         (3, 3, 6), (5, 3, 3), (3, 4, 3)]

# La unimodularidad ya NO se muestrea: se agota sobre cada matriz DISTINTA, que es un conjunto
# pequeno.  El tope antiguo desaparecio con la deduplicacion.
TOPE_CLASES = 5         # rango R' maximo para el que se calcula la forma canonica (coste R'!2^R')


def det_entero(M):
    n = len(M)
    s = 0
    for perm in itertools.permutations(range(n)):
        term = 1
        for i in range(n):
            term *= M[i][perm[i]]
            if term == 0:
                break
        if term == 0:
            continue
        p, visto = 1, [False] * n
        for i in range(n):
            if visto[i]:
                continue
            j, L = i, 0
            while not visto[j]:
                visto[j] = True
                j = perm[j]
                L += 1
            if L % 2 == 0:
                p = -p
        s += p * term
    return s


def es_TU(M):
    n = len(M)
    for k in range(1, n + 1):
        for filas in itertools.combinations(range(n), k):
            for cols in itertools.combinations(range(n), k):
                if abs(det_entero([[M[i][j] for j in cols] for i in filas])) > 1:
                    return False
    return True


def canon_signada(M):
    """Forma canonica de M bajo  M -> P M Q  con P, Q permutaciones CON SIGNO.

    Es la equivalencia exacta que la unimodularidad total no distingue: todo menor de P M Q es
    +- un menor de M, luego TU y |det| son constantes en cada clase.  Por eso basta testar un
    representante por clase, que es el ahorro que sugiere el algebra de Galois: escalar los residuos
    por una unidad es precisamente una permutacion con signo de filas y columnas.

    Se canonizan las filas por separado (min(v,-v)) porque las negaciones de fila son independientes,
    y se minimiza sobre las permutaciones con signo de COLUMNA.
    """
    n = len(M)
    mejor = None
    for perm in itertools.permutations(range(n)):
        for signos in itertools.product((1, -1), repeat=n):
            filas = []
            for i in range(n):
                v = tuple(signos[j] * M[i][perm[j]] for j in range(n))
                w = tuple(-x for x in v)
                filas.append(min(v, w))
            k = tuple(sorted(filas))
            if mejor is None or k < mejor:
                mejor = k
    return mejor


def matriz(V, X, t, mp, r, paridad="impar", con_signo_plegado=True):
    Rp = mp + r
    M = [[0] * Rp for _ in range(Rp)]
    for i in range(Rp):
        cl, ep = plegar(V[i], t)
        for p in range(mp):
            if cl == mp - p:
                M[i][p] = ep if con_signo_plegado else 1
        for j in range(r):
            tot = 0
            for e in (1, -1):
                num = e * V[i] - X[j]
                if num > 0 and num % t == 0:
                    k = num // t
                    if (k % 2 == 1) if paridad == "impar" else (k % 2 == 0):
                        tot += e
            M[i][mp + j] = tot
    return M


def columna_nula(M):
    n = len(M)
    for j in range(n):
        if all(M[i][j] == 0 for i in range(n)):
            return True
    return False


def c_cerrada(nu, X, t, r, M):
    ks = []
    for j in range(r):
        L, k = [], 1
        while X[j] + t * k <= M:
            if X[j] + t * k >= -M:
                L.append(k)
            k += 2
        if not L:
            return 0
        ks.append(L)
    s = 0
    for k in itertools.product(*ks):
        s += nu_extendida(nu, tuple(X[j] + t * k[j] for j in range(r)))
    return s


if __name__ == "__main__":
    print("=" * 108)
    print("BARRIDO ANCHO DE LA UNIMODULARIDAD")
    print("=" * 108)
    print("")
    sys.stdout.flush()

    D0 = D0n = TU = TUn = P0 = P0n = 0
    S1 = S1n = S2 = S2n = 0
    G1 = G1n = 0
    distintas = {}
    repeticiones = Counter()
    galois_fallos = []
    signos = Counter()
    absdet = Counter()
    entradas = Counter()
    malos = []
    no_tu = []
    por_caso = []
    rnd = random.Random(20260816)

    for (t, r, cota) in CASOS:
        mp = (t - 1) // 2
        Rp = mp + r
        d0 = d0n = tu = tun = vivos = pesos = 0
        sg_caso = Counter()
        for Lam in itertools.product(range(cota + 1), repeat=Rp):
            if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
                continue
            nu = nu_de(list(Lam), t, r)
            if not nu:
                continue
            pesos += 1
            V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
            Mx = max(max(abs(v) for v in k) for k in nu)
            E = eps_t(t, mp)
            # Los X con TODAS las columnas libres no nulas: se generan coordenada a coordenada en vez de
            # barrer la caja entera.  Es el mismo conjunto que sobrevive al prefiltro, y para r=4 la caja
            # completa son decenas de millones de matrices.
            permitidos = set()
            for i in range(Rp):
                for e in (1, -1):
                    k = 1
                    while True:
                        x = e * V[i] - t * k
                        if x < -Mx:
                            break
                        if x <= Mx:
                            permitidos.add(x)
                        k += 2
            permitidos = sorted(permitidos)
            descartados = []
            for X in itertools.product(permitidos, repeat=r):
                if any(X[j] <= X[j + 1] for j in range(r - 2)):
                    continue
                if not (X[r - 2] > abs(X[r - 1])):
                    continue
                Mm = matriz(V, X, t, mp, r)
                if columna_nula(Mm):
                    if len(descartados) < 40:
                        descartados.append((X, Mm))
                    continue
                vivos += 1
                for fila in Mm:
                    for e in fila:
                        entradas[e] += 1
                dd = det_entero(Mm)
                cc = c_cerrada(nu, X, t, r, Mx)
                absdet[abs(dd)] += 1
                d0n += 1
                D0n += 1
                if cc == E * dd:
                    d0 += 1
                    D0 += 1
                    if dd:
                        sg_caso[E] += 1
                        signos[(t, r, E)] += 1
                elif len(malos) < 3:
                    malos.append({"t": t, "r": r, "Lambda": list(Lam), "X": list(X),
                                  "det": dd, "c": cc, "eps": E})
                # DEDUPLICACION.  Las entradas de M solo dependen de los residuos y del orden, asi que
                # la misma matriz reaparece miles de veces al barrer (Lambda, X).  Guardando una por
                # forma, la unimodularidad se AGOTA en vez de muestrearse.
                clave = tuple(tuple(fila) for fila in Mm)
                repeticiones[(t, r)] += 1
                if clave not in distintas.setdefault((t, r), {}):
                    distintas[(t, r)][clave] = (list(Lam), list(X))
                # senuelos, solo donde hay algo que puntuar
                d1 = det_entero(matriz(V, X, t, mp, r, paridad="par"))
                if cc != 0 or d1 != 0:
                    S1n += 1
                    if abs(d1) == abs(cc):
                        S1 += 1
                d2 = det_entero(matriz(V, X, t, mp, r, con_signo_plegado=False))
                if cc != 0 or d2 != 0:
                    S2n += 1
                    if abs(d2) == abs(cc):
                        S2 += 1
            # P0  el prefiltro, comprobado sobre una muestra de descartados
            for (X, Mm) in rnd.sample(descartados, min(6, len(descartados))):
                P0n += 1
                if det_entero(Mm) == 0 and c_cerrada(nu, X, t, r, Mx) == 0:
                    P0 += 1
        # La unimodularidad, AGOTADA sobre las matrices distintas de este caso.
        dd = distintas.get((t, r), {})
        for clave, origen in dd.items():
            TUn += 1
            tun += 1
            if es_TU([list(f) for f in clave]):
                TU += 1
                tu += 1
            elif len(no_tu) < 2:
                no_tu.append({"t": t, "r": r, "Lambda": origen[0], "X": origen[1],
                              "M": [list(f) for f in clave]})
        # G1  y cuantas quedan SALVO permutaciones con signo -- el numero irreducible de tests.
        clases = 0
        if Rp <= TOPE_CLASES and dd:
            clases = len({canon_signada([list(f) for f in clave]) for clave in dd})
        etiq = ("%4d" % clases) if clases else ("  --" if Rp > TOPE_CLASES else "   0")
        print("  t=%2d r=%d Lambda_i<=%2d : %4d pesos | %7d X vivos | %5d distintas | %s clases | "
              "D0 %7d/%7d | TU %5d/%5d"
              % (t, r, cota, pesos, vivos, len(dd), etiq, d0, d0n, tu, tun))
        sys.stdout.flush()
        por_caso.append({"t": t, "r": r, "cota": cota, "pesos": pesos, "vivos": vivos,
                         "distintas": len(dd),
                         "clases_salvo_signo": (clases if Rp <= TOPE_CLASES else None),
                         "repeticiones": repeticiones[(t, r)],
                         "D0": [d0, d0n], "TU": [tu, tun], "signo": dict(sg_caso)})

    print("")
    print("-" * 108)
    print("  D0  FATAL  c == eps_t . det M sobre los X vivos : %d de %d" % (D0, D0n))
    print("  TU  FATAL  unimodularidad total AGOTADA sobre cada matriz distinta : %d de %d"
          % (TU, TUn))
    print("  P0  el prefiltro: descartado => det 0 y c 0     : %d de %d" % (P0, P0n))
    print("  D1  signo global por (t,r) : %s" % {str(k): v for k, v in sorted(signos.items())})
    print("  D2  reparto de |det M|     : %s" % dict(sorted(absdet.items())))
    print("  D3  entradas de M          : %s" % dict(sorted(entradas.items())))
    print("")
    print("  S1  SENUELO progresion par    : %d de %d" % (S1, S1n))
    print("  S2  SENUELO sin signo plegado : %d de %d" % (S2, S2n))
    if malos:
        print("")
        print("  !! desacuerdos en D0:")
        for m in malos:
            print("    " + json.dumps(m)[:300])
    if no_tu:
        print("")
        print("  !! matrices NO totalmente unimodulares:")
        for m in no_tu:
            print("    t=%d r=%d Lambda=%s X=%s" % (m["t"], m["r"], m["Lambda"], m["X"]))
            for fila in m["M"]:
                print("       %s" % fila)

    json.dump({"por_caso": por_caso, "D0": [D0, D0n], "TU": [TU, TUn], "tope_clases": TOPE_CLASES,
               "P0": [P0, P0n], "signos": {str(k): v for k, v in sorted(signos.items())},
               "abs_det": {str(k): v for k, v in sorted(absdet.items())},
               "entradas": {str(k): v for k, v in sorted(entradas.items())},
               "senuelos": {"par": [S1, S1n], "sin_signo": [S2, S2n]},
               "malos": malos, "no_TU": no_tu},
              open("unimodularidad_barrido_DUMP.json", "w"), indent=1)
    print("")
    print("=" * 108)
    print("DONE")
