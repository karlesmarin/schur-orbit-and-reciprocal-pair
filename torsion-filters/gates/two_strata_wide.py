# -*- coding: utf-8 -*-
# EL BARRIDO ANCHO, y las dos lecturas mias que mata.  12 de agosto de 2026 (noche).
#
# En una sola noche escribi DOS veces la misma clase de error: leer una firma de una poblacion que
# estaba truncada por el rango del barrido.
#
#   "todos t = 6"   -> FALSO: subiendo M salen en t = 4 (r = 2 y r = 3) y en t = 8 (r = 3).
#   "e = t en todos" -> FALSO: t = 8, r = 3, M = 21 los da con e = 6 < t = 8.
#
# Los dos los mato subir M.  Aqui se mide la poblacion en serio y, sobre todo, se somete a su PEOR
# poblacion el unico enunciado que sigue en pie:
#
#     "Phi_t != 0  =>  alguno de los TRES PRIMEROS estratos es no nulo"    (noche del 12)
#
# COLUMNAS
#   W1  poblacion por configuracion y por M: ¿crece?  ¿donde aparece la primera de cada t?
#   W2  la distribucion de e frente a t.  Si e = t fuera la firma, la columna e < t seria 0.
#   W3  PROFUNDIDAD por arriba de CADA contraejemplo, con la expansion ENTERA.  Prediccion <= 4.
#   W4  cuantos son ESENCIALMENTE distintos (canonicos bajo traslacion + complemento).
#   C0  CONTROL: sobre formas con [Phi]_top != 0 la profundidad debe ser 0.
#
# Guion AUTOCONTENIDO a proposito: two_strata_depth.py importaba de two_strata_audit.py y eso
# ejecuta el guion entero (sus tablas salen pegadas al principio de two_strata_depth_OUT.txt).
# Defecto mio, sin consecuencia en los numeros, pero no se repite.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python two_strata_wide.py

import itertools
import sys
from collections import defaultdict


def perm_sign(seq):
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def topdeg_dict(T, r):
    D = defaultdict(int)
    n = 2 * r
    for a in itertools.permutations(range(r)):
        for b in itertools.permutations(range(r)):
            q = [0] * n
            e = [0] * r
            for i in range(r):
                q[i] = 2 * a[i]
                e[a[i]] += T[i]
            for i in range(r):
                q[r + i] = 2 * b[i] + 1
                e[b[i]] -= T[r + i]
            D[tuple(e)] += perm_sign(q)
    return dict((k, v) for k, v in D.items() if v != 0)


def full_expansion(tr, r):
    n = 2 * r
    D = defaultdict(int)
    for x in tr:
        T, w = x[0], x[1]
        for q in itertools.permutations(range(n)):
            e = [0] * r
            for a in range(n):
                e[q[a] // 2] += (T[a] if q[a] % 2 == 0 else -T[a])
            D[tuple(e)] += w * perm_sign(list(q))
    return dict((k, v) for k, v in D.items() if v != 0)


def analyse(beta, t, r):
    N = len(beta)
    cl = defaultdict(list)
    for i, b in enumerate(beta):
        cl[b % t].append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    if not E:
        return None
    S = sorted((beta[i] for k in E for i in cl[k]), reverse=True)
    C = S[0] + S[-1]
    cond_i = set(C - v for v in S) == set(S)
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    cond_ii = (len(fixed) == 2 and all(k in E for k in fixed))
    tr = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        dif = tuple(T[2 * i] - T[2 * i + 1] for i in range(r))
        tr.append((T, w, sum(T[:r]) - sum(T[r:]), sum(dif), tuple(sorted(dif, reverse=True))))
    Dmax = max(x[2] for x in tr)
    Dmin = min(x[3] for x in tr)
    G = [x for x in tr if x[2] == Dmax]
    Gb = [x for x in tr if x[3] == Dmin]
    if len(G) == 1:
        top_zero = False
    else:
        top = defaultdict(int)
        for x in G:
            for k, v in topdeg_dict(list(x[0]), r).items():
                top[k] += x[1] * v
        top_zero = not any(v for v in top.values())
    byD = defaultdict(int)
    for x in Gb:
        byD[x[4]] += x[1]
    return dict(t=t, r=r, e=len(E), Dmax=Dmax, Dmin=Dmin, tr=tr,
                crit=(cond_i and cond_ii), cond_i=cond_i,
                top_zero=top_zero, bot_zero=all(v == 0 for v in byD.values()))


def is_fail(a):
    return a is not None and (not a['crit']) and a['top_zero'] and a['bot_zero']


def canon(beta):
    b = tuple(sorted((v - min(beta) for v in beta), reverse=True))
    c = tuple(sorted((max(b) - v for v in b), reverse=True))
    return min(b, c)


CONFIGS = [(4, 2, 17), (4, 2, 19), (4, 3, 17), (6, 3, 17), (6, 3, 19),
           (8, 3, 19), (8, 3, 20), (8, 3, 21), (8, 4, 19), (10, 5, 21)]

print("=" * 108)
print("W1/W2  la poblacion por configuracion, y e frente a t")
print("=" * 108)
print("")
print("     t   r    M       N | formas | contraejemplos | e = t | e < t | canonicos")
print("  " + "-" * 100)
FOUND = defaultdict(list)
for (t, r, M) in CONFIGS:
    N = t + 2 * r
    if M < N - 1:
        print("  %4d %3d %4d %7d | SALTADA (M < N-1)" % (t, r, M, N))
        continue
    n = nexc = net = nlt = 0
    cset = set()
    for comb in itertools.combinations(range(M + 1), N):
        beta = sorted(comb, reverse=True)
        a = analyse(beta, t, r)
        if a is None:
            continue
        n += 1
        if is_fail(a):
            nexc += 1
            if a['e'] == t:
                net += 1
            else:
                nlt += 1
            cset.add(canon(beta))
            FOUND[(t, r, M)].append((tuple(beta), a))
    print("  %4d %3d %4d %7d | %6d | %14d | %5d | %5d | %d"
          % (t, r, M, N, n, nexc, net, nlt, len(cset)))
    sys.stdout.flush()

print("")
print("=" * 108)
print("W3  PROFUNDIDAD de cada contraejemplo, con la expansion ENTERA.  Prediccion de la noche: <= 4")
print("=" * 108)
print("")
prof = defaultdict(int)
profb = defaultdict(int)
deep = []
byconf = defaultdict(lambda: defaultdict(int))
nzero = 0
for key in sorted(FOUND):
    for beta, a in FOUND[key]:
        FE = full_expansion(a['tr'], a['r'])
        if not FE:
            nzero += 1
            continue
        d = a['Dmax'] - max(sum(k) for k in FE)
        db = min(sum(abs(x) for x in k) for k in FE) - a['Dmin']
        prof[d] += 1
        profb[db] += 1
        byconf[key][d] += 1
        if d > 4:
            deep.append((key, beta, d))
    sys.stdout.flush()

print("     configuracion | profundidad por arriba")
print("  " + "-" * 100)
for key in sorted(byconf):
    print("     t=%d r=%d M=%-3d | %s" % (key[0], key[1], key[2],
                                          "  ".join("%d:%d" % (k, byconf[key][k])
                                                    for k in sorted(byconf[key]))))
print("")
print("  TOTAL por arriba: %s" % "  ".join("%d:%d" % (k, prof[k]) for k in sorted(prof)))
print("  TOTAL por abajo : %s" % "  ".join("%d:%d" % (k, profb[k]) for k in sorted(profb)))
print("  Phi_t == 0 entre los contraejemplos (imposible por construccion): %d" % nzero)
print("")
if deep:
    print("  *** %d CONTRAEJEMPLOS POR DEBAJO DEL TERCER ESTRATO ***" % len(deep))
    for key, beta, d in deep[:10]:
        print("      t=%d r=%d  beta = %s  profundidad %d" % (key[0], key[1], list(beta), d))
else:
    print("  NINGUNO baja del tercer estrato.")

# ---------------------------------------------------------------- C0 ----------------------------
print("")
print("=" * 108)
print("C0  CONTROL: con [Phi]_top != 0 la profundidad debe ser 0")
print("=" * 108)
print("")
n = bad = 0
for comb in itertools.combinations(range(18), 10):
    beta = sorted(comb, reverse=True)
    a = analyse(beta, 4, 3)
    if a is None or a['crit'] or a['top_zero']:
        continue
    n += 1
    if n > 200:
        break
    FE = full_expansion(a['tr'], 3)
    if not FE or a['Dmax'] - max(sum(k) for k in FE) != 0:
        bad += 1
print("  %d formas : %d con profundidad != 0  (debe ser 0)" % (min(n, 200), bad))
print("")
print("DONE")
