# -*- coding: utf-8 -*-
# LA DESCOMPOSICION QUE FALTABA: c(X) COMO UN DETERMINANTE.
#
# Otra descomposicion de W^1.  Las probadas hasta ahora parten el SIGNO en
# factores (enderezado, barajado, delta) y ninguna se concentra.  El fallo es de fondo: todas suman
# sobre SUBCONJUNTOS S y luego intentan explicar el signo aparte.
#
# Pero un elemento de W^1 no es un subconjunto: es una permutacion con signo de V, que reparte los
# R' valores entre m' ranuras congeladas y r ranuras libres.  Y la formula cerrada
#
#       c(X) = sum_{k impar} nu~(X + t k)
#
# es, escrita asi, una suma CON SIGNO sobre emparejamientos perfectos entre los V_i y las ranuras.
# Una suma con signo sobre emparejamientos perfectos es un DETERMINANTE.
#
# La matriz.  Filas = los V_i.  Columnas = las ranuras:
#   - m' ranuras congeladas, la p-esima pide la clase plegada  c = m' - p  (el orden decreciente de
#     lem:T).  Entrada = el signo del plegado si  plegar(V_i) = (c, eps),  y 0 si no.
#   - r ranuras libres, la j-esima pide el valor  X_j.  Entrada = suma de los eps en {+1,-1} tales
#     que  eps*V_i = X_j + t*k  con k impar positivo;  0 si ninguno.
#
# Si  c(X) = +- det(M),  entonces (L1) deja de ser una pregunta sobre cancelaciones y pasa a ser:
#
#       ESTA MATRIZ 0/+-1 TIENE DETERMINANTE EN {0,+-1}
#
# que es un enunciado de unimodularidad, con literatura entera detras (Ghouila-Houri, matrices de
# intervalos, matrices de red).  Ese seria el cambio de categoria que buscamos.
#
# CONTROLES
#   D0 (FATAL)  det(M) == c(X) salvo un signo global constante por (t, r).
#   D1          el signo global, .es constante?  Si cambia caso a caso, la identidad es falsa y lo
#               que hay es una coincidencia de valores absolutos.
#   D2          reparto de |det(M)| -- tiene que ser {0,1} si D0 sale.
#   D3          las entradas de M estan en {0,+-1} (si aparece un +-2 la lectura de arriba falla).
#   S1 senuelo  la misma matriz con las ranuras libres pidiendo k PAR: debe fallar.
#   S2 senuelo  la misma matriz sin los signos del plegado (todo +1): debe fallar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_determinante.py

import itertools
import json
from collections import Counter

from divided_differences import CASOS, plegar, nu_de, nu_extendida, eps_t


def det_entero(M):
    n = len(M)
    s = 0
    for perm in itertools.permutations(range(n)):
        p = 1
        visto = [False] * n
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
        term = p
        for i in range(n):
            term *= M[i][perm[i]]
            if term == 0:
                break
        s += term
    return s


def es_TU(M):
    """Totalmente unimodular: TODO menor cuadrado tiene determinante en {0,+-1}."""
    n = len(M)
    for k in range(1, n + 1):
        for filas in itertools.combinations(range(n), k):
            for cols in itertools.combinations(range(n), k):
                sub = [[M[i][j] for j in cols] for i in filas]
                if abs(det_entero(sub)) > 1:
                    return False
    return True


def matriz(V, X, t, mp, r, paridad="impar", con_signo_plegado=True):
    Rp = mp + r
    M = [[0] * Rp for _ in range(Rp)]
    for i in range(Rp):
        for p in range(mp):
            c = mp - p
            cl, ep = plegar(V[i], t)
            if cl == c:
                M[i][p] = ep if con_signo_plegado else 1
        for j in range(r):
            tot = 0
            for ep in (1, -1):
                num = ep * V[i] - X[j]
                if num > 0 and num % t == 0:
                    k = num // t
                    if (k % 2 == 1) if paridad == "impar" else (k % 2 == 0):
                        tot += ep
            M[i][mp + j] = tot
    return M


def c_referencia(nu, t, r):
    """La formula cerrada ya validada (256/256 contra el pelado)."""
    if not nu:
        return {}
    M = max(max(abs(v) for v in k) for k in nu)
    out = {}
    for X in itertools.product(range(-M, M + 1), repeat=r):
        if any(X[j] <= X[j + 1] for j in range(r - 2)):
            continue
        if not (X[r - 2] > abs(X[r - 1])):
            continue
        ks = []
        for j in range(r):
            L, k = [], 1
            while X[j] + t * k <= M:
                if X[j] + t * k >= -M:
                    L.append(k)
                k += 2
            ks.append(L)
        if any(not L for L in ks):
            continue
        s = 0
        for k in itertools.product(*ks):
            s += nu_extendida(nu, tuple(X[j] + t * k[j] for j in range(r)))
        if s:
            out[X] = s
    return out


print("=" * 104)
print("c(X) COMO DETERMINANTE:  .es (L1) un enunciado de unimodularidad?")
print("=" * 104)
print("")

D0 = D0n = 0
S1 = S1n = S2 = S2n = 0
signos = Counter()
absdet = Counter()
entradas = Counter()
testigos = []
no_tu = []
TU = [0, 0]
TUn = [0, 0]
por_caso = []

for (t, r, cota) in CASOS:
    mp = (t - 1) // 2
    Rp = mp + r
    d0 = d0n = s1 = s1n = s2 = s2n = 0
    sg_caso = Counter()
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu = nu_de(list(Lam), t, r)
        if not nu:
            continue
        V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
        ref = c_referencia(nu, t, r)
        Mx = max(max(abs(v) for v in k) for k in nu)
        for X in itertools.product(range(-Mx, Mx + 1), repeat=r):
            if any(X[j] <= X[j + 1] for j in range(r - 2)):
                continue
            if not (X[r - 2] > abs(X[r - 1])):
                continue
            M = matriz(V, X, t, mp, r)
            for fila in M:
                for e in fila:
                    entradas[e] += 1
            dd = det_entero(M)
            cc = ref.get(X, 0)
            absdet[abs(dd)] += 1
            d0n += 1
            if dd == 0 and cc == 0:
                d0 += 1
            elif dd != 0 and cc != 0 and abs(dd) == abs(cc):
                d0 += 1
                sg_caso[cc // dd] += 1
                signos[(t, r, cc // dd)] += 1
            elif len(testigos) < 4:
                testigos.append({"t": t, "r": r, "Lambda": list(Lam), "X": list(X), "V": V,
                                 "det": dd, "c": cc, "M": M})
            # SENUELOS.  Puntuarlos sobre TODO X los infla: 86636 de 87652 tienen c = 0 y
            # cualquier matriz rota acierta ahi por no hacer nada.  Se puntuan solo donde hay algo
            # que acertar: c != 0, o el senuelo inventa un det != 0 donde no lo hay.
            d1 = det_entero(matriz(V, X, t, mp, r, paridad="par"))
            if cc != 0 or d1 != 0:
                s1n += 1
                if abs(d1) == abs(cc):
                    s1 += 1
            d2 = det_entero(matriz(V, X, t, mp, r, con_signo_plegado=False))
            if cc != 0 or d2 != 0:
                s2n += 1
                if abs(d2) == abs(cc):
                    s2 += 1
            # C4  .es M TOTALMENTE unimodular?  Si lo es, hay caracterizacion clasica que invocar.
            # Se examinan las dos poblaciones por separado: las 1016 con c != 0 y una muestra
            # igual de grande con c == 0, que son la mayoria y donde un fallo tambien contaria.
            cual = 0 if cc != 0 else 1
            if TUn[cual] < 1000:
                TUn[cual] += 1
                if not es_TU(M):
                    TU[cual] += 1
                    if len(no_tu) < 2:
                        no_tu.append({"t": t, "r": r, "X": list(X), "M": M})
    D0 += d0; D0n += d0n; S1 += s1; S1n += s1n; S2 += s2; S2n += s2n
    print("  t=%d r=%d :  det(M) reproduce c en %5d de %5d   | signos globales vistos: %s"
          % (t, r, d0, d0n, dict(sg_caso)))
    por_caso.append({"t": t, "r": r, "D0": d0, "n": d0n, "signos": {str(k): v for k, v in sg_caso.items()}})

print("")
print("-" * 104)
print("  D0  FATAL  det(M) == c(X) salvo signo global : %d de %d" % (D0, D0n))
print("  D1  signos globales por (t,r) : %s" % {str(k): v for k, v in sorted(signos.items())})
print("      contra eps_t (cor:oddsign), calculado y no leido a ojo:")
ok_eps = True
for (tt, rr, sg), n in sorted(signos.items()):
    e = eps_t(tt, (tt - 1) // 2)
    marca = "coincide" if sg == e else "NO COINCIDE"
    ok_eps = ok_eps and sg == e
    print("        t=%d r=%d : signo global %+d, eps_t = %+d  -> %s" % (tt, rr, sg, e, marca))
print("      el signo global ES eps_t : %s" % ("SI" if ok_eps else "NO"))
print("  C4  M totalmente unimodular:  con c != 0 : %d de %d   |   con c == 0 : %d de %d"
      % (TUn[0] - TU[0], TUn[0], TUn[1] - TU[1], TUn[1]))
print("      (tope declarado: 1000 matrices por poblacion, no es exhaustivo)")
print("  D2  reparto de |det(M)|       : %s" % dict(sorted(absdet.items())))
print("  D3  entradas de M             : %s" % dict(sorted(entradas.items())))
print("")
print("  S1  SENUELO progresion par    : %d de %d" % (S1, S1n))
print("  S2  SENUELO sin signo plegado : %d de %d" % (S2, S2n))
if testigos:
    print("")
    print("  !! primeros desacuerdos:")
    for tt in testigos[:2]:
        print("    t=%d r=%d Lambda=%s X=%s V=%s  det=%d  c=%d"
              % (tt["t"], tt["r"], tt["Lambda"], tt["X"], tt["V"], tt["det"], tt["c"]))
        for fila in tt["M"]:
            print("       %s" % fila)
print("")
print("  LECTURA: si D0 y D1 salen, (L1) es 'esta matriz 0/+-1 tiene determinante en {0,+-1}', o sea")
print("  unimodularidad, y eso tiene literatura entera detras.  Si D2 ensena un 2, (L1) es FALSA y")
print("  tengo un contraejemplo, que seria aun mas informativo.")

json.dump({"por_caso": por_caso, "D0": [D0, D0n],
           "signos": {str(k): v for k, v in sorted(signos.items())},
           "abs_det": {str(k): v for k, v in sorted(absdet.items())},
           "entradas": {str(k): v for k, v in sorted(entradas.items())},
           "senuelos": {"par": [S1, S1n], "sin_signo": [S2, S2n]},
           "testigos": testigos[:2], "TU_fallos": {"c_no_nulo": [TU[0], TUn[0]], "c_nulo": [TU[1], TUn[1]]}, "no_TU": no_tu},
          open("_probe_determinante_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
