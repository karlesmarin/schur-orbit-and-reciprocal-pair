# -*- coding: utf-8 -*-
# EL ESQUELETO FORZADO: QUE ELEMENTO DE K_lambda ES CADA MONOMIO DE LAPLACE.  20 de agosto de 2026.
#
# POR QUE ESTA GATE EXISTE.  pI_toggle_espacios_de_peso.py dejo (4) en forma LOCAL --- en el espacio
# de peso e hay que fijar exactamente los |coef_e| monomios de Laplace de ese peso ---, pero (4)
# habla de OBJETOS y lo que tenemos es un conteo.  Para construir el toggle hace falta la aplicacion
#      monomio de Laplace  -->  elemento de K_lambda
# y no la tenemos.  Aqui esta el unico sitio donde esa aplicacion NO es una eleccion nuestra: los
# espacios de peso desequilibrados de TAMANO 1 (23.5%% en t=2, 28.9%% en t=3).  Alli el punto fijo
# esta forzado, se puede LEER, y si hay regla tiene que verse en ellos.
#
# LAS HIPOTESIS, ESCRITAS ANTES DE MIRAR
#   H1  el omega de un punto fijo forzado es siempre uno de los DOS elementos de peso extremo de
#       Omega_t --- `z^t|-1|z` (peso t+1) y `-1|z^-t|-z^-1` (peso -t-1).  Salio de W5 y aqui se
#       comprueba como implicacion, no como frecuencia.
#   H2  el mu de un punto fijo forzado tiene beta-set relacionado con el de lambda QUITANDO las dos
#       cuentas j1, j2 de la transversal de Laplace de ese peso.  Es lo que dice L3: la transversal
#       borra dos columnas y deja t con residuos distintos.
#   H3  en un punto fijo forzado, nu = lambda (la tira de arriba es vacia) o nu = mu (la de abajo).
#       Es la degeneracion que uno espera de un extremo.
#   H4  sgn_t(mu) del punto fijo forzado == el signo del monomio de Laplace, salvo el global.
#       Si falla, la aplicacion no puede ser «el punto fijo ES la transversal» sin un twist.
#
# QUE SE MIDE
#   S0  CALIBRACION: los espacios de tamano 1 y desequilibrio 1 son exactamente los que W3 conto.
#   S1  H1, como implicacion sobre TODOS los forzados.
#   S2  H2: reparto de beta(mu) contra beta(lambda) menos {j1,j2}, en las cuatro formas que puede
#       tomar la relacion.  Se imprime el reparto, no un si/no --- no se de antemano cual es.
#   S3  H3 y H4.
#   S4  y lo que hace falta para el toggle: .cuantos monomios de Laplace quedan SIN forzar, o sea en
#       espacios de tamano >= 2?  Esa es la parte que una regla tendra que decidir.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_toggle_esqueleto.py > pI_toggle_esqueleto_OUT.txt 2>&1

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
    s, L = 1, list(seq)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            if L[i] > L[j]:
                s = -s
    return s


def omega(t):
    out = []
    for c in itertools.product([(t, +1, "z^t"), (0, -1, "-1")],
                               [(-t, +1, "z^-t"), (0, -1, "-1")],
                               [(1, +1, "z"), (-1, -1, "-z^-1")]):
        out.append((sum(x[0] for x in c), c[0][1] * c[1][1] * c[2][1],
                    "|".join(x[2] for x in c)))
    return out


print("=" * 100)
print("EL ESQUELETO FORZADO --- .que elemento de K_lambda es cada monomio de Laplace?")
print("=" * 100)
print("")

resumen = {}

for T in CASOS:
    N = T + 2
    OM = omega(T)
    EXTREMOS = {et for (w, s, et) in OM if abs(w) == T + 1}
    D = pmul(pmul({T: 1, 0: -1}, {-T: 1, 0: -1}), {1: 1, -1: -1})

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

    forzados = 0
    s1_ok = 0
    s1_et = Counter()
    s1b = Counter()
    s2 = Counter()
    s3_nu = Counter()
    s4_ok, s4_tot = 0, 0
    libres = 0
    lap_tot = 0
    ej = []

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
                J.append((2 * sum(nu) - sum(mu) - sum(lam), sg, mu, nu))

        phi = schur_directo(lam)
        izq = pmul(phi, D)

        l_ = list(lam) + [0] * N
        beta = [l_[i] + N - 1 - i for i in range(N)]
        glob = (-1) ** (T + (N + 1) * N // 2)
        # peso -> lista de (j1, j2, signo del monomio)
        lap = defaultdict(list)
        for (j1, j2) in itertools.combinations(range(N), 2):
            S = [j for j in range(N) if j not in (j1, j2)]
            res = [beta[j] % T for j in S]
            if len(set(res)) != T:
                continue
            eps = glob * (-1) ** (j1 + j2) * sgn_perm(res)
            for k, v in f(beta[j1] - beta[j2]).items():
                lap[k].append((j1, j2, eps * v))
        lap_tot += sum(len(v) for v in lap.values())

        espacios = defaultdict(list)
        for (e, s, mu, nu) in J:
            for (w, sw, et) in OM:
                espacios[e + w].append((s * sw, et, mu, nu, s))

        for e, lst in espacios.items():
            c = sum(x[0] for x in lst)
            if not c:
                continue
            if len(lst) == 1:
                forzados += 1
                sgnk, et, mu, nu, sgt = lst[0]
                # S1
                if et in EXTREMOS:
                    s1_ok += 1
                s1_et[et] += 1
                # H1b: omega toma EXACTAMENTE UNO de los dos factores de torsion
                partes = et.split("|")
                cuantos = (1 if partes[0] != "-1" else 0) + (1 if partes[1] != "-1" else 0)
                s1b[cuantos] += 1
                # S2: beta(mu) contra beta(lambda) sin {j1,j2}
                bmu = sorted(beta_set(mu, N), reverse=True)
                cand = lap[e]
                if len(cand) == 1:
                    j1, j2, sl = cand[0]
                    resto = sorted([beta[j] for j in range(N) if j not in (j1, j2)], reverse=True)
                    if bmu == resto:
                        s2["beta(mu) == beta(lambda) sin {j1,j2}"] += 1
                    elif set(resto).issubset(set(bmu)):
                        s2["beta(lambda) sin {j1,j2} contenido en beta(mu)"] += 1
                    elif set(bmu).issubset(set(beta)):
                        s2["beta(mu) contenido en beta(lambda), otra cosa"] += 1
                    else:
                        s2["ninguna de las anteriores"] += 1
                    # S4: signo
                    s4_tot += 1
                    if (1 if sgnk > 0 else -1) == (1 if sl > 0 else -1):
                        s4_ok += 1
                    if len(ej) < 6:
                        ej.append((tuple(lam), e, et, tuple(mu), tuple(nu), sgt, (j1, j2), sl))
                else:
                    s2["el peso tiene %d monomios de Laplace" % len(cand)] += 1
                # S3
                if tuple(nu) == tuple(lam):
                    s3_nu["nu == lambda"] += 1
                elif tuple(nu) == tuple(mu):
                    s3_nu["nu == mu"] += 1
                else:
                    s3_nu["nu intermedio"] += 1
            else:
                libres += abs(c)

    n = len(lams)
    print("-" * 100)
    print("t=%d   %d lambda,  lambda_1 <= %d" % (T, n, LMAX))
    print("-" * 100)
    print("  S0  espacios desequilibrados de tamano 1 (punto fijo FORZADO) : %d" % forzados)
    print("  S1  H1  su omega es uno de los DOS extremos de Omega_t : %d de %d   %s"
          % (s1_ok, forzados, "SIN EXCEPCION" if s1_ok == forzados else "*** FALLA ***"))
    print("        etiquetas que aparecen : %s" % dict(s1_et))
    print("  S1b H1b .cuantos factores de torsion (z^t, z^-t) toma su omega? : %s"
          % dict(sorted(s1b.items())))
    print("        --> %s" % ("EXACTAMENTE UNO, sin excepcion" if set(s1b) == {1}
                              else "*** no siempre uno ***"))
    print("  S2  H2  beta(mu) contra beta(lambda) menos las dos cuentas de la transversal")
    for k, v in sorted(s2.items(), key=lambda x: -x[1]):
        print("        %6d  %s" % (v, k))
    print("  S3  H3  que es nu en un punto fijo forzado : %s" % dict(s3_nu))
    print("  S4  H4  signo del elemento == signo del monomio de Laplace : %d de %d   %s"
          % (s4_ok, s4_tot, "SIEMPRE" if s4_ok == s4_tot else "coincide en %.1f%%"
             % (100.0 * s4_ok / max(1, s4_tot))))
    print("  S5  monomios de Laplace TOTALES %d ; forzados %d ; por decidir (espacios >= 2) %d"
          % (lap_tot, forzados, libres))
    print("  ejemplos (lambda, peso, omega, mu, nu, sgn_t(mu), (j1,j2), signo Laplace):")
    for x in ej[:6]:
        print("        %s" % (x,))
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = dict(n=int(n), forzados=int(forzados), S1=int(s1_ok),
                               S2={k: int(v) for k, v in s2.items()},
                               S3={k: int(v) for k, v in s3_nu.items()},
                               S4=[int(s4_ok), int(s4_tot)],
                               lap_tot=int(lap_tot), libres=int(libres))

json.dump(resumen, open("pI_toggle_esqueleto_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
