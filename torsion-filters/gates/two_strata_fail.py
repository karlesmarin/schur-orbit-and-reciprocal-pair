# -*- coding: utf-8 -*-
# LAS OCHO QUE ROMPEN EL ENUNCIADO DE DOS ESTRATOS.  Verificacion dura, 12 de agosto de 2026 (noche, tras el commit 48c231e).
#
# QUE SE VERIFICA.  bottom_sees_gcom.py encontro 8 formas con  [Phi]_top = 0,  [Phi]_bot = 0  y
# sigma_C(S) != S.  Si son reales, el enunciado
#       [Phi]_top = 0  y  [Phi]_bot = 0   =>   sigma_C(S) = S
# es FALSO, y con el la ruta de dos estratos.  Antes de creerselo hay que descartar que sea un
# defecto de mis rebanadas: aqui se expande Phi_t ENTERO, monomio a monomio, sin usar ninguna de las
# dos formulas de estrato.
#
# LAS DOS GRADUACIONES, que NO son la misma y el pipeline usa las dos:
#   grado total   sum e_j        -- el estrato de ARRIBA es su maximo (Dmax)
#   grado absoluto sum |e_j|     -- el estrato de ABAJO es su minimo (Dmin)
# Ambas dan condiciones necesarias para Phi_t = 0.  Se comprueban las dos contra la expansion.
#
# COLUMNAS, cada una capaz de fallar
#   V1  Phi_t != 0 de verdad (el diccionario de monomios no es vacio).  Si fuera 0, no habria
#       contraejemplo: seria una forma anulante con S no concentrico, que contradiria la necesidad
#       de (i) -- medida en 9071 formas.  Las dos salidas son informativas y hay que decir cual es.
#   V2  max sum e_j de la expansion == Dmax del pipeline, y su rebanada es VACIA (top = 0).
#   V3  min sum |e_j| de la expansion == Dmin del pipeline, y su rebanada es VACIA (bot = 0).
#   V4  donde esta el PRIMER estrato no nulo, en las dos graduaciones.  La noche del 12 midio
#       "profundidad <= 4" (D1-0, D1-2, D1-4) en 100288 formas: la prediccion es D1-2 o D1-4.
#   V5  CONTROL QUE PUEDE FALLAR: la misma expansion sobre una forma de branch (b) construida a
#       mano (autocomplementaria) debe dar Phi_t == 0.  Si no, el expansor esta mal y V1 no vale.
#   V6  CONTROL: sobre una forma generica cualquiera, la rebanada Dmax de la expansion debe
#       coincidir con la formula del estrato de arriba (sector todo-positivo).
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python two_strata_fail.py

import itertools
from collections import defaultdict

FAILS = [
    (6, 3, [16, 15, 14, 12, 11, 7, 5, 4, 3, 2, 1, 0]),
    (6, 3, [16, 15, 14, 13, 12, 11, 9, 5, 4, 2, 1, 0]),
    (6, 3, [16, 14, 13, 11, 9, 7, 6, 5, 4, 3, 2, 0]),
    (6, 3, [16, 14, 13, 12, 11, 10, 9, 7, 5, 3, 2, 0]),
    (6, 3, [17, 16, 15, 13, 12, 8, 6, 5, 4, 3, 2, 1]),
    (6, 3, [17, 16, 15, 14, 13, 12, 10, 6, 5, 3, 2, 1]),
    (6, 3, [17, 15, 14, 12, 10, 8, 7, 6, 5, 4, 3, 1]),
    (6, 3, [17, 15, 14, 13, 12, 11, 10, 8, 6, 4, 3, 1]),
]


def perm_sign(seq):
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def transversals(beta, t, r):
    N = len(beta)
    cl = defaultdict(list)
    for i, b in enumerate(beta):
        cl[b % t].append(i)
    if len(cl) < t:
        return None, None
    out = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        out.append((w, T))
    return out, dict(cl)


def full_expansion(tm, r):
    """Phi_t entero: sum_T w(T) * A(T), A(T) = det sobre las columnas (z_j, 1/z_j)."""
    n = 2 * r
    D = defaultdict(int)
    for w, T in tm:
        for q in itertools.permutations(range(n)):
            e = [0] * r
            for a in range(n):
                e[q[a] // 2] += (T[a] if q[a] % 2 == 0 else -T[a])
            D[tuple(e)] += w * perm_sign(list(q))
    return dict((k, v) for k, v in D.items() if v != 0)


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


def report(t, r, beta, label):
    beta = sorted(beta, reverse=True)
    tm, cl = transversals(beta, t, r)
    if tm is None:
        print("  %s : CLASE VACIA (branch (a)), Phi_t = 0 por palomar" % label)
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    S = sorted((beta[i] for k in E for i in cl[k]), reverse=True)
    C = S[0] + S[-1]
    cond_i = set(C - v for v in S) == set(S)
    Dmax = max(sum(T[:r]) - sum(T[r:]) for _, T in tm)
    Dmin = min(sum(T[2 * i] - T[2 * i + 1] for i in range(r)) for _, T in tm)

    FE = full_expansion(tm, r)
    if not FE:
        print("  %-22s Phi_t == 0   (i) concentrico: %s" % (label, cond_i))
        return dict(zero=True, cond_i=cond_i)

    tot = defaultdict(int)
    ab = defaultdict(int)
    for k in FE:
        tot[sum(k)] += 1
        ab[sum(abs(x) for x in k)] += 1
    mx, mn = max(tot), min(ab)
    top_slice = sum(1 for k in FE if sum(k) == Dmax)
    bot_slice = sum(1 for k in FE if sum(abs(x) for x in k) == Dmin)
    d_top = Dmax - mx
    d_bot = mn - Dmin
    print("  %-22s Phi_t != 0 (%4d monomios) | Dmax %3d, primer no nulo %3d (baja %d) | "
          "Dmin %2d, primer no nulo %2d (sube %d) | (i) %s"
          % (label, len(FE), Dmax, mx, d_top, Dmin, mn, d_bot, "SI" if cond_i else "no"))
    return dict(zero=False, cond_i=cond_i, Dmax=Dmax, mx=mx, Dmin=Dmin, mn=mn,
                top_slice=top_slice, bot_slice=bot_slice, d_top=d_top, d_bot=d_bot,
                tm=tm, FE=FE)


print("=" * 108)
print("V1-V4  las ocho formas, expandidas ENTERAS")
print("=" * 108)
print("")
res = []
for i, (t, r, beta) in enumerate(FAILS):
    res.append(report(t, r, beta, "t=%d r=%d  #%d" % (t, r, i + 1)))

bad = 0
for x in res:
    if x is None or x['zero']:
        bad += 1
        continue
    if x['top_slice'] or x['bot_slice']:
        bad += 1

print("")
print("  V1  Phi_t != 0 en las ocho              : %d de 8" % sum(1 for x in res if x and not x['zero']))
print("  V2  la rebanada Dmax es VACIA           : %d de 8"
      % sum(1 for x in res if x and not x['zero'] and x['top_slice'] == 0))
print("  V3  la rebanada Dmin es VACIA           : %d de 8"
      % sum(1 for x in res if x and not x['zero'] and x['bot_slice'] == 0))
print("  V4  profundidad por arriba (Dmax - primer no nulo) : %s"
      % sorted(set(x['d_top'] for x in res if x and not x['zero'])))
print("      profundidad por abajo  (primer no nulo - Dmin) : %s"
      % sorted(set(x['d_bot'] for x in res if x and not x['zero'])))

# ---------------------------------------------------------------- V5, V6 -------------------------
print("")
print("=" * 108)
print("V5  CONTROL: branch (b) construida a mano -- la expansion DEBE dar 0")
print("=" * 108)
print("")
# t=2, r=2, N=6: lambda autocomplementaria con w impar.  beta = lambda_i + N - i.
# lambda = (3,3,3,0,0,0) -> lambda_i + lambda_{N+1-i} = 3 (impar), C par.
for lam in [(3, 3, 3, 0, 0, 0), (5, 5, 5, 0, 0, 0), (2, 2, 2, 2, 0, 0)]:
    N = 6
    beta = [lam[i] + N - 1 - i for i in range(N)]
    w = lam[0] + lam[-1]
    lab = "lambda=%s (w=%d, %s)" % (str(lam), w, "impar -> branch (b)" if w % 2 else "par")
    report(2, 2, beta, lab)

print("")
print("=" * 108)
print("V6  CONTROL: la rebanada Dmax de la expansion == la formula del estrato de arriba")
print("=" * 108)
print("")
nchk = nbad = 0
for (t, r, M) in [(4, 2, 12), (6, 3, 17)]:
    N = t + 2 * r
    seen = 0
    for comb in itertools.combinations(range(M + 1), N):
        beta = sorted(comb, reverse=True)
        tm, cl = transversals(beta, t, r)
        if tm is None:
            continue
        seen += 1
        if seen > 12:
            break
        Dmax = max(sum(T[:r]) - sum(T[r:]) for _, T in tm)
        FE = full_expansion(tm, r)
        sl = dict((k, v) for k, v in FE.items() if sum(k) == Dmax)
        top = defaultdict(int)
        for w, T in tm:
            if sum(T[:r]) - sum(T[r:]) == Dmax:
                for k, v in topdeg_dict(list(T), r).items():
                    top[k] += w * v
        top = dict((k, v) for k, v in top.items() if v != 0)
        nchk += 1
        if top != sl:
            nbad += 1
    print("  t=%d r=%d : %d formas comprobadas" % (t, r, min(seen, 12)))
print("")
print("  V6 : %d desacuerdos de %d" % (nbad, nchk))

print("")
if bad or nbad:
    print("VEREDICTO SUSPENDIDO -- alguna columna de aceptacion fallo.")
else:
    print("VEREDICTO: las ocho son REALES.  '[Phi]_top = 0 y [Phi]_bot = 0 => sigma_C(S) = S' es FALSO.")
print("DONE")
