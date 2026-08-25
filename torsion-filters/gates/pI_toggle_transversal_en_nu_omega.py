# -*- coding: utf-8 -*-
# .ESTA LA TRANSVERSAL EN (nu, omega)?   20 de agosto de 2026.
#
# ⚠️ LO QUE NO CUENTA COMO RESPUESTA.  El peso e es funcion de (mu,nu,omega), y entre los puntos
# fijos forzados de peso con UN solo monomio de Laplace, la transversal es funcion del peso.  Luego
# «(mu,nu,omega) determina la transversal» es TAUTOLOGICO y no se mide.  Lo que NO es tautologico es
# que (nu,omega) baste --- ahi falta |mu| --- ni, mucho menos, que haya una REGLA legible.
#
# QUE SE MIDE
#   N1  FATAL, y es la pregunta: .determina (nu,omega) la transversal?  O sea, .hay dos puntos fijos
#       forzados de la MISMA lambda con el mismo (nu,omega) y transversales DISTINTAS?  Y lo mismo
#       para nu solo y para omega solo, que son los dos casos degenerados.
#   N2  si la determina, .con que REGLA?  Se prueban once lecturas candidatas de {j1,j2} como
#       conjunto de filas o de cuentas, y se imprime el acierto de cada una.  No se cual es: por eso
#       se prueban todas y se publica la tabla, no la ganadora.
#   N3  el cruce omega x «que transversal», por si omega selecciona un tipo y nu el resto.
#   D1  SENUELO: la misma pregunta sobre puntos fijos ELEGIDOS AL AZAR en espacios de tamano >= 2.
#       Si (nu,omega) «determinara» la transversal tambien ahi, no estaria midiendo nada del
#       esqueleto sino una coincidencia de rango.
#
# ####################################################################################################
# N4  LA COMPOSICION, y es de esta misma manana.  El criterio de A_r = B_r salio de una sola
# observacion: un determinante de INDICADORES [t | a_i - b_j] solo es no nulo si el multiconjunto de
# residuos de las filas coincide con el de las columnas.  Y la condicion de transversal de L3 es
# LA MISMA ESPECIE de enunciado: quitar dos cuentas j1, j2 de beta(lambda) de modo que las t
# restantes tengan residuos TODOS DISTINTOS.  Compuestas, se predice el conteo del Lema L3:
#
#   beta(lambda) tiene N = t+2 cuentas sobre t residuos, luego un EXCESO de 2.  Para que al borrar
#   dos quede un sistema completo de residuos hacen falta multiplicidades m_r >= 1 para todo r, y
#   entonces el exceso se reparte de dos maneras y solo dos:
#       un residuo con m = 3   -->  C(3,2) = 3 transversales
#       dos residuos con m = 2 -->  2 x 2  = 4 transversales
#       algun residuo con m = 0 -->  0 transversales
#
#   O sea |U| in {0, 3, 4}, que es EXACTAMENTE lo que el Paper I demuestra --- y aqui sale de contar
#   multiplicidades en el abaco, en tres lineas.  Se mide como implicacion, casilla a casilla.
#   Y de paso da una lectura candidata de {j1,j2}: las dos cuentas borradas VIVEN forzosamente en
#   los residuos sobre-representados.  Eso acota donde puede estar la transversal; cual de las 3 o 4
#   es, sigue siendo lo que (nu,omega) tendria que decidir.
# ####################################################################################################
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_toggle_transversal_en_nu_omega.py > pI_toggle_transversal_en_nu_omega_OUT.txt 2>&1

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


def pad(p, n):
    return list(p) + [0] * (n - len(p))


def lecturas(lam, mu, nu, N, beta, et):
    """Once lecturas candidatas de {j1,j2} a partir de (lambda, mu, nu, omega)."""
    L, M, V = pad(lam, N), pad(mu, N), pad(nu, N)
    A = frozenset(i for i in range(N) if L[i] > V[i])          # filas de la tira lambda/nu
    B = frozenset(i for i in range(N) if V[i] > M[i])          # filas de la tira nu/mu
    C = frozenset(i for i in range(N) if L[i] > M[i])
    todo = frozenset(range(N))
    bl = set(beta)
    bnu = set(beta_set(nu, N))
    bmu = set(beta_set(mu, N))
    return {
        "A  filas de lambda/nu": A,
        "B  filas de nu/mu": B,
        "C  filas de lambda/mu": C,
        "~A complemento de A": todo - A,
        "~B complemento de B": todo - B,
        "~C complemento de C": todo - C,
        "A u B": A | B,
        "A ^ B  (dif. simetrica)": A ^ B,
        "A n B": A & B,
        "filas j con beta_j fuera de beta(nu)": frozenset(j for j in range(N) if beta[j] not in bnu),
        "filas j con beta_j fuera de beta(mu)": frozenset(j for j in range(N) if beta[j] not in bmu),
    }


print("=" * 100)
print(".ESTA LA TRANSVERSAL EN (nu, omega)?")
print("=" * 100)
print("")

resumen = {}

for T in CASOS:
    N = T + 2
    OM = omega(T)
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

    n1_nuom_col, n1_nu_col, n1_om_col = 0, 0, 0
    n1_tot = 0
    n2 = Counter()
    n2_tot = 0
    n3 = Counter()
    d1_nuom_col, d1_tot = 0, 0
    ej_col = []
    n4_ok, n4_tot, n4_mal, n4_dentro = 0, 0, [], 0
    n4_forma = Counter()

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
        l_ = pad(lam, N)
        beta = [l_[i] + N - 1 - i for i in range(N)]
        glob = (-1) ** (T + (N + 1) * N // 2)
        lap = defaultdict(list)
        for (j1, j2) in itertools.combinations(range(N), 2):
            S = [j for j in range(N) if j not in (j1, j2)]
            res = [beta[j] % T for j in S]
            if len(set(res)) != T:
                continue
            eps = glob * (-1) ** (j1 + j2) * sgn_perm(res)
            for k, v in f(beta[j1] - beta[j2]).items():
                lap[k].append((j1, j2, eps * v))

        # --- N4: la composicion.  |U| predicho por las multiplicidades de residuo del abaco
        mult = Counter(b % T for b in beta)
        nU = sum(1 for (j1, j2) in itertools.combinations(range(N), 2)
                 if len(set(beta[j] % T for j in range(N) if j not in (j1, j2))) == T)
        if any(mult.get(r, 0) == 0 for r in range(T)):
            pred, forma = 0, "algun residuo con m=0"
        elif 3 in mult.values():
            pred, forma = 3, "un residuo con m=3"
        elif list(mult.values()).count(2) == 2:
            pred, forma = 4, "dos residuos con m=2"
        else:
            pred, forma = -1, "perfil no previsto %s" % sorted(mult.values())
        n4_tot += 1
        if pred == nU:
            n4_ok += 1
        else:
            if len(n4_mal) < 5:
                n4_mal.append((tuple(lam), sorted(mult.values()), pred, nU))
        n4_forma[(forma, nU)] += 1
        sobre = frozenset(j for j in range(N) if mult[beta[j] % T] >= 2)

        espacios = defaultdict(list)
        for (e, s, mu, nu) in J:
            for (w, sw, et) in OM:
                espacios[e + w].append((s * sw, et, mu, nu))

        # --- los forzados
        por_nuom = defaultdict(set)
        por_nu = defaultdict(set)
        por_om = defaultdict(set)
        # --- senuelo: uno cualquiera de cada espacio grande, para contrastar
        d_por_nuom = defaultdict(set)

        for e, lst in espacios.items():
            c = sum(x[0] for x in lst)
            if not c:
                continue
            cand = lap.get(e, [])
            if len(cand) != 1:
                continue
            j1, j2, sl = cand[0]
            tr = (j1, j2)
            if len(lst) == 1:
                sgnk, et, mu, nu = lst[0]
                n1_tot += 1
                por_nuom[(tuple(nu), et)].add(tr)
                por_nu[tuple(nu)].add(tr)
                por_om[et].add(tr)
                n3[(et, tr)] += 1
                # N2: las once lecturas
                n2_tot += 1
                for nombre, S in lecturas(lam, mu, nu, N, beta, et).items():
                    if S == frozenset(tr):
                        n2[nombre] += 1
                if frozenset(tr) <= sobre:
                    n4_dentro += 1
            else:
                # senuelo: TODOS los elementos del espacio grande, no uno elegido
                for (sgnk, et, mu, nu) in lst:
                    d1_tot += 1
                    d_por_nuom[(tuple(nu), et)].add(tr)

        for k, v in por_nuom.items():
            if len(v) > 1:
                n1_nuom_col += 1
                if len(ej_col) < 4:
                    ej_col.append((tuple(lam), k, sorted(v)))
        for k, v in por_nu.items():
            if len(v) > 1:
                n1_nu_col += 1
        for k, v in por_om.items():
            if len(v) > 1:
                n1_om_col += 1
        for k, v in d_por_nuom.items():
            if len(v) > 1:
                d1_nuom_col += 1

    n = len(lams)
    print("-" * 100)
    print("t=%d   %d lambda,  lambda_1 <= %d   ---  %d puntos fijos FORZADOS" % (T, n, LMAX, n1_tot))
    print("-" * 100)
    print("  N1  FATAL  .determina (nu,omega) la transversal, dentro de cada lambda?")
    print("        colisiones (mismo (nu,omega), transversales distintas) : %d   %s"
          % (n1_nuom_col, "NINGUNA -- (nu,omega) LA DETERMINA" if n1_nuom_col == 0
             else "SI las hay, p.ej. %s" % (ej_col[:2],)))
    print("        y los dos casos degenerados:  solo nu : %d colisiones    solo omega : %d colisiones"
          % (n1_nu_col, n1_om_col))
    print("  D1  SENUELO la misma pregunta sobre los espacios de tamano >= 2 (%d elementos)" % d1_tot)
    print("        colisiones : %d   %s"
          % (d1_nuom_col, "*** tampoco colisiona: N1 no mide el esqueleto ***" if d1_nuom_col == 0
             else "colisiona --> N1 dice algo del esqueleto"))
    print("  N2  .con que REGLA?  aciertos de once lecturas de {j1,j2}, sobre %d forzados" % n2_tot)
    for nombre in sorted(n2, key=lambda x: -n2[x]):
        print("        %5d / %-5d  (%5.1f%%)   %s" % (n2[nombre], n2_tot,
                                                      100.0 * n2[nombre] / max(1, n2_tot), nombre))
    if not n2:
        print("        NINGUNA de las once acierta ni una vez")
    print("  N4  COMPOSICION  |U| predicho por las multiplicidades de residuo del abaco : %d de %d   %s"
          % (n4_ok, n4_tot, "SIN EXCEPCION" if n4_ok == n4_tot
             else "*** FALLA en %d: %s ***" % (n4_tot - n4_ok, n4_mal[:3])))
    print("        perfil -> |U| : %s" % dict(sorted(n4_forma.items())))
    print("        y {j1,j2} contenido en las cuentas de residuo sobre-representado : %d de %d   %s"
          % (n4_dentro, n2_tot, "SIEMPRE" if n4_dentro == n2_tot else "*** no siempre ***"))
    print("  N3  cruce omega x transversal (top 8)")
    for k, v in sorted(n3.items(), key=lambda x: -x[1])[:8]:
        print("        %5d x  omega=%-18s  (j1,j2)=%s" % (v, k[0], k[1]))
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = dict(n=int(n), forzados=int(n1_tot),
                               N1=[int(n1_nuom_col), int(n1_nu_col), int(n1_om_col)],
                               D1=[int(d1_nuom_col), int(d1_tot)],
                               N2={k: int(v) for k, v in n2.items()}, N2_tot=int(n2_tot))

json.dump(resumen, open("pI_toggle_transversal_en_nu_omega_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
