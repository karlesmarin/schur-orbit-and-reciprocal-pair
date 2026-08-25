# -*- coding: utf-8 -*-
# DONDE PUEDE VIVIR LA INVOLUCION: LOS ESPACIOS DE PESO DE K_lambda.   20 de agosto de 2026.
#
# POR QUE ESTA GATE EXISTE.  pI_numerator_lift.py construyo K_lambda = J_lambda x Omega_t y midio
# que el minimo de puntos fijos es F_K = ||D_t.Phi||_1 = 2|U|, el numero de monomios de Laplace.
# Eso dice que la condicion (4) de 41_IN es POSIBLE.  No dice DONDE estan los puntos fijos, que es
# lo unico que hace falta para construir el toggle.  Esta gate mira dentro de los espacios de peso.
#
# LO QUE YA ES TAUTOLOGICO, Y NO SE MIDE COMO SI FUERA UN HALLAZGO.  Como sum_{K} sgn z^wt = D_t.Phi
# (verificado 329/329 y 791/791), el desequilibrio n+(e) - n-(e) de cada espacio de peso ES el
# coeficiente de z^e en D_t.Phi.  Luego «los pesos desequilibrados son los de Laplace» y «los signos
# coinciden» NO son medidas: son la identidad otra vez.  Se anotan como controles, no como resultado.
#
# LO QUE SI ES CONTENIDO
#   W1  ⚠️ PREDICCION MIA, Y LA GATE LA MATO.  Escribi que de F_K = 2|U| se seguia |coef| <= 1 en
#       todo espacio de peso, o sea que los 2|U| monomios de Laplace ocupan 2|U| pesos DISTINTOS.
#       FALSO: ||suma||_1 = #terminos no exige pesos distintos --- dos monomios que coinciden en peso
#       CON EL MISMO SIGNO se refuerzan y la norma L1 no lo nota.  Medido: |coef| llega a 3 en t=2
#       (lambda = (1,1,1,1), peso -1) y a 2 en t=3.  Se conserva el numero como medida.
#   W1b LA FORMA LOCAL BUENA, que es lo que W1 queria decir y no decia:  para CADA peso e,
#            |coef_e(D_t.Phi)|  ==  #{monomios de Laplace de peso e} ,  y todos con el mismo signo.
#       Sale de que sum_e |coef_e| = 2|U| = total de monomios y de la desigualdad triangular por
#       peso: la igualdad en la suma la fuerza en cada sumando.  Eso convierte (4) en una condicion
#       LOCAL --- «en el espacio de peso e hay que dejar fijos exactamente los monomios de Laplace de
#       peso e» --- que es lo unico que un toggle puede cumplir espacio a espacio.
#   W2  EL TAMANO de los espacios de peso --- eso es de K_lambda, no del polinomio, y es la holgura
#       que tiene el toggle.  Reparto general y reparto RESTRINGIDO a los desequilibrados.
#   W3  .hay espacios desequilibrados de tamano 1?  Ahi el punto fijo esta FORZADO y no hay nada que
#       elegir: seria el esqueleto del toggle canonico.  Cuantos, y que fraccion de los 2|U|.
#   W4  .son los pesos desequilibrados los EXTREMOS del soporte de K_lambda?  Si lo fueran, «los
#       puntos fijos son los extremos» seria una regla candidata, y barata de comprobar.
#   W5  la ANATOMIA de un espacio desequilibrado: .de que parte de Omega_t vienen sus elementos?
#       Si los desequilibrados solo tocan ciertas elecciones de Omega, el toggle se puede definir
#       ahi y propagar.
#   V1  SENUELO: los mismos repartos sobre J_lambda (sin levantar).  Alli W1 tiene que FALLAR ---
#       ||Phi||_1 llega a 30 con coeficientes grandes --- y es el contraste que justifica el lift.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_toggle_espacios_de_peso.py > pI_toggle_espacios_de_peso_OUT.txt 2>&1

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


def sgn_perm(seq):
    s, L = 1, list(seq)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            if L[i] > L[j]:
                s = -s
    return s


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


def omega(t):
    """Los 8 terminos de (z^t-1)(z^-t-1)(z-z^-1) sin expandir; etiqueta = la eleccion por factor."""
    out = []
    for c in itertools.product([(t, +1, "z^t"), (0, -1, "-1")],
                               [(-t, +1, "z^-t"), (0, -1, "-1")],
                               [(1, +1, "z"), (-1, -1, "-z^-1")]):
        out.append((sum(x[0] for x in c), c[0][1] * c[1][1] * c[2][1],
                    "|".join(x[2] for x in c)))
    return out


print("=" * 100)
print("LOS ESPACIOS DE PESO DE K_lambda --- donde puede vivir el toggle")
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

    D = pmul(pmul({T: 1, 0: -1}, {-T: 1, 0: -1}), {1: 1, -1: -1})
    OM = omega(T)
    lams = [l for l in particiones(LMAX, N) if l]

    w1_max = 0
    w1_mal = []
    w1b_ok, w1b_tot, w1b_mal = 0, 0, []
    lap_sin_deseq = 0
    tam_todos = Counter()
    tam_deseq = Counter()
    w3_forzados = 0
    w3_deseq_tot = 0
    w4_extremos = 0
    w4_deseq = 0
    w4_todos_extremos = 0
    w5 = Counter()
    v1_max = 0
    v1_coef = Counter()
    calib = 0

    for lam in lams:
        J = []
        for mu in particiones(lam[0], len(lam) + 1):
            if not sub(mu, lam):
                continue
            core, quot, sg = core_quot_sgn(mu, N, T)
            if len(core) != 0 or any(len(q) > 1 for q in quot):
                continue
            for nu in particiones(lam[0], len(lam) + 1):
                if not (sub(mu, nu) and sub(nu, lam)):
                    continue
                if not (hstrip(nu, mu) and hstrip(lam, nu)):
                    continue
                J.append((2 * sum(nu) - sum(mu) - sum(lam), sg))

        phi = schur_directo(lam)
        izq = pmul(phi, D)

        # los 2|U| monomios de Laplace, uno a uno (peso, signo), SIN sumarlos
        l_ = list(lam) + [0] * N
        beta = [l_[i] + N - 1 - i for i in range(N)]
        glob = (-1) ** (T + (N + 1) * N // 2)
        monos = []
        for (j1, j2) in itertools.combinations(range(N), 2):
            S = [j for j in range(N) if j not in (j1, j2)]
            res = [beta[j] % T for j in S]
            if len(set(res)) != T:
                continue
            eps = glob * (-1) ** (j1 + j2) * sgn_perm(res)
            for k, v in f(beta[j1] - beta[j2]).items():
                monos.append((k, eps * v))
        porpeso_lap = defaultdict(list)
        for (k, s) in monos:
            porpeso_lap[k].append(s)

        espacios = defaultdict(list)     # e -> [(sgn, etiqueta de Omega)]
        for (e, s) in J:
            for (w, sw, et) in OM:
                espacios[e + w].append((s * sw, et))

        # control de calibracion: el desequilibrio ES el coeficiente
        ok = True
        for e, lst in espacios.items():
            c = sum(x[0] for x in lst)
            if c != izq.get(e, 0):
                ok = False
        if ok and all(e in espacios for e in izq):
            calib += 1

        if espacios:
            emin, emax = min(espacios), max(espacios)
        else:
            emin = emax = None

        deseq = []
        for e, lst in sorted(espacios.items()):
            c = sum(x[0] for x in lst)
            tam_todos[len(lst)] += 1
            if c:
                deseq.append(e)
                tam_deseq[len(lst)] += 1
                if abs(c) > w1_max:
                    w1_max = abs(c)
                if abs(c) > 1 and len(w1_mal) < 5:
                    w1_mal.append((tuple(lam), e, c))
                # W1b: la forma local
                lap = porpeso_lap.get(e, [])
                w1b_tot += 1
                if abs(c) == len(lap) and lap and all(x == lap[0] for x in lap) and \
                        (1 if c > 0 else -1) == lap[0]:
                    w1b_ok += 1
                elif len(w1b_mal) < 5:
                    w1b_mal.append((tuple(lam), e, c, lap))
                if len(lst) == 1:
                    w3_forzados += 1
                w3_deseq_tot += 1
                if e in (emin, emax):
                    w4_extremos += 1
                w4_deseq += 1
                # W5: de que elecciones de Omega vienen los elementos del espacio desequilibrado
                w5[tuple(sorted(set(x[1] for x in lst)))] += 1
        if deseq and set(deseq) == {emin, emax}:
            w4_todos_extremos += 1
        # .hay pesos de Laplace que NO son desequilibrados?  (cancelacion entre transversales)
        for e in porpeso_lap:
            if izq.get(e, 0) == 0:
                lap_sin_deseq += 1

        # V1: lo mismo sobre J_lambda, sin levantar
        v1_coef.update(abs(v) for v in phi.values())
        if phi:
            v1_max = max(v1_max, max(abs(v) for v in phi.values()))

    n = len(lams)
    print("-" * 100)
    print("t=%d   %d lambda,  lambda_1 <= %d" % (T, n, LMAX))
    print("-" * 100)
    print("  ctrl  el desequilibrio de cada espacio ES el coeficiente de D_t.Phi : %d de %d   %s"
          % (calib, n, "CALIBRA" if calib == n else "*** FALLA ***"))
    print("  W1  ⚠️ prediccion MUERTA: |coef| <= 1 en todo espacio.  Maximo medido %d   %s"
          % (w1_max, "(habria pasado)" if w1_max <= 1 else "FALSA -- p.ej. %s" % (w1_mal[:3],)))
    print("      dos monomios de Laplace pueden caer en el MISMO peso con el MISMO signo:")
    print("      la norma L1 no lo nota, luego F_K = 2|U| no implicaba pesos distintos")
    print("  W1b FATAL, la forma local buena:  |coef_e| == #monomios de Laplace de peso e, mismo signo")
    print("        %d de %d espacios desequilibrados   %s"
          % (w1b_ok, w1b_tot, "EXACTO -- (4) es una condicion LOCAL"
             if w1b_ok == w1b_tot else "*** FALLA en %d: %s ***" % (w1b_tot - w1b_ok, w1b_mal[:3])))
    print("        pesos de Laplace SIN desequilibrio (transversales que se cancelan) : %d"
          % lap_sin_deseq)
    print("  W2  tamano de los espacios de peso de K_lambda")
    print("        todos         : %s" % dict(sorted(tam_todos.items())))
    print("        desequilibrados: %s" % dict(sorted(tam_deseq.items())))
    print("  W3  espacios desequilibrados de tamano 1 (punto fijo FORZADO) : %d de %d   (%.1f%%)"
          % (w3_forzados, w3_deseq_tot, 100.0 * w3_forzados / max(1, w3_deseq_tot)))
    print("  W4  .son los pesos desequilibrados los EXTREMOS del soporte?")
    print("        desequilibrados que son extremo : %d de %d ; lambda donde SON exactamente los dos extremos : %d de %d"
          % (w4_extremos, w4_deseq, w4_todos_extremos, n))
    print("  W5  de que elecciones de Omega_t vienen los espacios desequilibrados (top 6)")
    for k, v in sorted(w5.items(), key=lambda x: -x[1])[:6]:
        print("        %5d x  %s" % (v, " + ".join(k)))
    print("  V1  SENUELO sobre J_lambda sin levantar: |coef de Phi| maximo %d ; reparto %s"
          % (v1_max, dict(sorted(v1_coef.items()))))
    print("        --> %s" % ("*** W1 no discriminaria ***" if v1_max <= 1
                              else "alli (4) es imposible; el lift es lo que la hace local"))
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = dict(n=int(n), calib=int(calib), W1_max=int(w1_max),
                               tam_todos={str(k): int(v) for k, v in tam_todos.items()},
                               tam_deseq={str(k): int(v) for k, v in tam_deseq.items()},
                               W3=[int(w3_forzados), int(w3_deseq_tot)],
                               W4=[int(w4_extremos), int(w4_deseq), int(w4_todos_extremos)],
                               V1_max=int(v1_max))

json.dump(resumen, open("pI_toggle_espacios_de_peso_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
