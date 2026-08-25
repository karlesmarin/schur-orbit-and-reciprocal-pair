# -*- coding: utf-8 -*-
# EL NUMERATOR LIFT  K_lambda = J_lambda x Omega_t  --- CONSTRUIDO, no descrito.  20 de agosto de 2026.
#
# POR QUE ESTA GATE EXISTE.  F8 lleva desde el 19 como «el lado Paper I: el toggle LGV de primer
# cruce, sobre el numerator lift K_lambda = J_lambda x Omega_t (no sobre J_lambda).  Sigue sin gate».
# Y lo que hay medido son las dos puntas y no el objeto de en medio:
#   - sobre J_lambda,  |Fix| >= ||Phi||_1,  y ||Phi||_1 llega a 30 (t=2) y 21 (t=3);
#   - el lado de Laplace, antes de L5, tiene ||.||_1 en {0,6,8}, maximo 8.
# 43_IN §2 dice cual es el objeto correcto: multiplicar COMBINATORIAMENTE por el denominador,
#     K_lambda := J_lambda x Omega_t ,  wt((mu,nu),w) = e(mu,nu) + d(w) ,  sgn = sgn_t(mu).sgn(w)
# con Omega_t un conjunto signado que represente los terminos de D_t = (z^t-1)(z^-t-1)(z-z^-1), y
# «sin expandir D_t prematuramente»: cada factor aporta DOS elecciones, luego |Omega_t| = 8 y no 6.
#
# LO QUE DECIDE ESTA GATE.  Sobre J_lambda la condicion (4) de 41_IN es ARITMETICAMENTE IMPOSIBLE:
# ninguna involucion monomial-preserving puede dejar 3 o 4 objetos de Laplace cuando el minimo de
# puntos fijos es 30.  La pregunta que F8 tenia abierta es si el levantamiento arregla eso o solo lo
# mueve de sitio.  El minimo sobre K_lambda es, por el mismo argumento de weight spaces,
#     F_K(lambda) = sum_e |n+(e) - n-(e)|   medido SOBRE K_lambda   =   ||D_t . Phi||_1,
# y ese numero hay que compararlo con 2.|U|, que son los monomios de las |U| transversales de Laplace
# (cada una lleva un menor 2x2 f(beta_j1-beta_j2), que son DOS monomios --- la sutileza de 43_IN §2).
#
# QUE SE MIDE
#   L0  CALIBRACION 1: la suma signada sobre J_lambda reproduce Phi_t(lambda;z).  Es la calibracion
#       de pI_involucion_toggle.py, repetida aqui para que esta gate no dependa de otra.
#   L1  FATAL, y es la construccion: sum_{K_lambda} sgn(x) z^{wt(x)} == D_t(z) . Phi_t(lambda;z).
#       Si esto no pasa, el objeto que 43_IN describe no es el que hemos construido.
#   L2  |K_lambda| = 8 . |J_lambda|, y el reparto de signos.
#   L3  IDENTIDAD verificada: F_K(lambda) == ||D_t . Phi||_1.
#   L4  FATAL, y es la pregunta de F8:  F_K(lambda) == 2 . |U(lambda)|, con U las transversales
#       admisibles de L3 (Paper I).  Si pasa, el minimo de puntos fijos sobre K_lambda es EXACTAMENTE
#       el numero de monomios de Laplace, y la condicion (4) deja de ser imposible por conteo.
#   L5  el contraste que es el resultado: histograma de F_J = ||Phi||_1 (sobre J_lambda) contra F_K
#       (sobre K_lambda), y cuantas lambda tenian F_J > 8.
#   E1  SENUELO: Omega_t sin el factor (z - z^-1), o sea 4 elementos.  L1 tiene que romperse.
#   E2  SENUELO: Omega_t con todos los signos +1.  L1 tiene que romperse.
#   E3  SENUELO: Omega_t con D_t EXPANDIDO a sus 6 monomios (que da el mismo polinomio pero otro
#       conjunto).  L1 pasa --- tiene que pasar --- pero |K| deja de ser 8|J| y F_K NO cambia:
#       sirve para ver que lo que decide es el POLINOMIO y no la particion en 8 o en 6 piezas.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_numerator_lift.py > pI_numerator_lift_OUT.txt 2>&1

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


def norma1(p):
    return sum(abs(v) for v in p.values())


def omega(t, sin_ultimo=False, signos_positivos=False):
    """Omega_t: los terminos de (z^t-1)(z^-t-1)(z-z^-1) SIN expandir --- una eleccion por factor.
    Devuelve lista de (peso, signo).  |Omega_t| = 8."""
    factores = [[(t, +1), (0, -1)], [(-t, +1), (0, -1)]]
    if not sin_ultimo:
        factores.append([(1, +1), (-1, -1)])
    out = []
    for combo in itertools.product(*factores):
        w = sum(c[0] for c in combo)
        s = 1
        for c in combo:
            s *= c[1]
        out.append((w, 1 if signos_positivos else s))
    return out


print("=" * 100)
print("EL NUMERATOR LIFT  K_lambda = J_lambda x Omega_t  --- CONSTRUIDO")
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
    OM_e1 = omega(T, sin_ultimo=True)
    OM_e2 = omega(T, signos_positivos=True)
    OM_e3 = [(k, v) for k, v in sorted(D.items())]      # D_t EXPANDIDO: 6 piezas, no 8

    lams = [l for l in particiones(LMAX, N) if l]

    l0_ok, l0_bad = 0, []
    l1_ok, l1_bad = 0, []
    l2_ok = 0
    l3_ok, l3_bad = 0, []
    l4_ok, l4_bad = 0, []
    l4b_ok, l4c_tot, l4c_phi0 = 0, 0, 0
    FJ = Counter()
    FK = Counter()
    FJ_mayor_8 = 0
    e1_rompe, e2_rompe, e3_ok = 0, 0, 0
    e3_tam_igual = 0
    sgn_rep = Counter()
    Ktot, Jtot = 0, 0

    for lam in lams:
        # ---- J_lambda
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

        # L0: la suma signada sobre J reproduce Phi
        sJ = {}
        for (e, s) in J:
            sJ = padd(sJ, {e: s})
        if sJ == phi:
            l0_ok += 1
        else:
            l0_bad.append(tuple(lam))

        # ---- K_lambda = J x Omega
        def suma_sobre(OMx):
            acc = {}
            for (e, s) in J:
                for (w, sw) in OMx:
                    acc = padd(acc, {e + w: s * sw})
            return acc

        izq = pmul(phi, D)
        sK = suma_sobre(OM)
        if sK == izq:
            l1_ok += 1
        else:
            l1_bad.append(tuple(lam))

        # L2
        Jtot += len(J)
        nK = len(J) * len(OM)
        Ktot += nK
        if nK == 8 * len(J):
            l2_ok += 1
        for (e, s) in J:
            for (w, sw) in OM:
                sgn_rep[s * sw] += 1

        # L3: F_K por weight spaces sobre K
        porexp = defaultdict(lambda: [0, 0])
        for (e, s) in J:
            for (w, sw) in OM:
                porexp[e + w][0 if s * sw > 0 else 1] += 1
        F_K = sum(abs(a - b) for (a, b) in porexp.values())
        if F_K == norma1(izq):
            l3_ok += 1
        else:
            l3_bad.append((tuple(lam), F_K, norma1(izq)))

        # L4: F_K == 2 |U|
        l = list(lam) + [0] * N
        beta = [l[i] + N - 1 - i for i in range(N)]
        nU = 0
        for (j1, j2) in itertools.combinations(range(N), 2):
            S = [j for j in range(N) if j not in (j1, j2)]
            if len(set(beta[j] % T for j in S)) == T:
                nU += 1
        if F_K == 2 * nU:
            l4_ok += 1
        else:
            l4_bad.append((tuple(lam), F_K, nU, norma1(phi)))
        if F_K <= 2 * nU:
            l4b_ok += 1
        if F_K != 2 * nU:
            l4c_tot += 1
            if phi == {}:
                l4c_phi0 += 1

        # L5
        F_J = norma1(phi)
        FJ[F_J] += 1
        FK[F_K] += 1
        if F_J > 8:
            FJ_mayor_8 += 1

        # senuelos
        if suma_sobre(OM_e1) != izq:
            e1_rompe += 1
        if suma_sobre(OM_e2) != izq:
            e2_rompe += 1
        if suma_sobre(OM_e3) == izq:
            e3_ok += 1
        if len(J) * len(OM_e3) == 8 * len(J):
            e3_tam_igual += 1

    n = len(lams)
    print("-" * 100)
    print("t=%d   %d lambda,  lambda_1 <= %d,  %d columnas de beta   |Omega_t| = %d"
          % (T, n, LMAX, N, len(OM)))
    print("-" * 100)
    print("  L0  CALIBRACION  suma signada sobre J_lambda == Phi_t : %d de %d   %s"
          % (l0_ok, n, "CALIBRA" if not l0_bad else "*** FALLA %s ***" % l0_bad[:4]))
    print("  L1  FATAL  sum_{K} sgn . z^wt == D_t . Phi_t : %d de %d   %s"
          % (l1_ok, n, "EL OBJETO ES ESE" if not l1_bad else "*** FALLA %s ***" % l1_bad[:4]))
    print("  L2  |K_lambda| = 8 |J_lambda| : %d de %d ;  total %d = 8 x %d ;  signos %s"
          % (l2_ok, n, Ktot, Jtot, dict(sgn_rep)))
    print("  L3  IDENTIDAD  F_K(lambda) == ||D_t . Phi||_1 : %d de %d   %s"
          % (l3_ok, n, "PASA" if not l3_bad else "*** FALLA %s ***" % l3_bad[:4]))
    print("  L4  FATAL  F_K(lambda) == 2 . |U(lambda)|  (monomios de Laplace) : %d de %d   %s"
          % (l4_ok, n, "EXACTO" if not l4_bad else "no exacto en %d, ej. (lambda, F_K, |U|, ||Phi||_1) %s"
             % (len(l4_bad), l4_bad[:4])))
    print("      L4b  la COTA  F_K <= 2|U| : %d de %d   %s"
          % (l4b_ok, n, "SIEMPRE" if l4b_ok == n else "*** ROTA ***"))
    print("      L4c  de las %d sin igualdad, cuantas tienen Phi == 0 : %d   %s"
          % (l4c_tot, l4c_phi0,
             "TODAS -- la igualdad vale siempre que Phi != 0" if l4c_tot == l4c_phi0
             else "*** hay %d con Phi != 0 ***" % (l4c_tot - l4c_phi0)))
    print("  L5  EL CONTRASTE")
    print("        F_J = ||Phi||_1  sobre J_lambda : %s    maximo %d" % (dict(sorted(FJ.items())), max(FJ)))
    print("        F_K              sobre K_lambda : %s    maximo %d" % (dict(sorted(FK.items())), max(FK)))
    print("        lambda con F_J > 8 (donde (4) es IMPOSIBLE sobre J) : %d de %d" % (FJ_mayor_8, n))
    print("  E1  SENUELO Omega sin el factor (z - z^-1) : rompe en %d de %d   %s"
          % (e1_rompe, n, "DISCRIMINA" if e1_rompe else "*** NO ***"))
    print("  E2  SENUELO Omega con todos los signos +1 : rompe en %d de %d   %s"
          % (e2_rompe, n, "DISCRIMINA" if e2_rompe else "*** NO ***"))
    print("  E3  CONTROL Omega = D_t EXPANDIDO (6 piezas) : la identidad sigue pasando en %d de %d,"
          % (e3_ok, n))
    print("        y |K| = 8|J| deja de valer en %d de %d  --> lo que decide es el POLINOMIO, no la"
          % (n - e3_tam_igual, n))
    print("        particion; las 8 piezas son la eleccion de 43_IN para conservar la procedencia")
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = dict(n=int(n), L0=int(l0_ok), L1=int(l1_ok), L2=int(l2_ok),
                               L3=int(l3_ok), L4=int(l4_ok),
                               FJ={str(k): int(v) for k, v in FJ.items()},
                               FK={str(k): int(v) for k, v in FK.items()},
                               FJ_mayor_8=int(FJ_mayor_8), Ktot=int(Ktot), Jtot=int(Jtot),
                               E1=int(e1_rompe), E2=int(e2_rompe), E3=int(e3_ok))

json.dump(resumen, open("pI_numerator_lift_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
