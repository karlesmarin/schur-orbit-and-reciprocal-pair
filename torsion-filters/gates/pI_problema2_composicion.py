# -*- coding: utf-8 -*-
# PROBLEMA 2 DEL PAPER I: .COMPONEN LOS DOS FORMALISMOS?   19 de agosto de 2026, noche.
#
# EL PROBLEMA, LITERAL (paper/orbit_pair.tex, prob:versch):
#   «The zeta-orbits are the province of phi_t, in Albion's formulation [Alb23].  Alphabets closed
#    under inversion are the province of the universal symplectic and orthogonal functions ...  Our
#    alphabet is one of each ...  WE KNOW OF NOTHING THAT ACTS ON BOTH HALVES AT ONCE, and the
#    question is whether the two formalisms compose.»
#
# LO QUE LA LECTURA DEL 19 DE AGOSTO APORTA, Y EL PROBLEMA NO SABIA:
#   (1) [Alb23] §6.1 SI define un operador que ve una letra libre: phi^q_t, con
#       phi^q_t h_{at+b} := q^b sum_k q^{kt} h_{a-k},  y su Prop. 6.2 lo calcula sobre s_lambda.
#       O sea que para UNA letra el objeto existe y esta publicado --- en el mismo articulo que
#       nuestro Problema cita como el formalismo de la mitad torcida.
#   (2) Y la composicion para DOS letras es inmediata en cuanto se escribe:  phi_t actua solo sobre
#       X, y la regla de ramificacion separa el alfabeto, luego
#
#         s_lambda(mu_t U {z,1/z}) = sum_{mu subset lambda} (phi_t s_mu)(1) * s_{lambda/mu}(z,1/z)
#
#       donde (phi_t s_mu)(1) es el Teorema 3.1 de [Alb23] evaluado en 1: CERO salvo que core_t(mu)
#       sea vacio, y entonces sgn_t(mu) * prod_r s_{mu^{(r)}}(1), que a su vez vale 1 exactamente
#       cuando CADA componente del t-cociente tiene a lo sumo UNA FILA.
#
# LA PREGUNTA QUE ESTA GATE CONTESTA.  .Es esa suma igual a Phi_t(lambda;z)?  Si lo es, la respuesta
# al Problema 2 es SI --- los dos formalismos componen, y la composicion es explicita ---, y lo que
# queda abierto deja de ser «existe algo?» para ser «por que el resultado FACTORIZA en tres senos
# cuando la suma no lo sugiere».
#
# QUE SE MIDE
#   C0  CALIBRACION: s_nu(1) = 1 si l(nu) <= 1 y 0 si no; y h_k del alfabeto contra el producto.
#   C1  FATAL: la composicion contra el valor DIRECTO, igualdad exacta de polinomios de Laurent.
#   C2  cuantos mu sobreviven por lambda --- si fueran pocos, la suma es un objeto manejable y ese
#       es justamente el camino que el Problema pide.
#   C3  la lectura del sumando: .es cierto que solo sobreviven los mu cuyo t-cociente tiene todas
#       las componentes de a lo sumo una fila?
#   D1  SENUELO: sin la condicion de core vacio (sumar TODOS los mu).  Tiene que fallar.
#   D2  SENUELO: con s_{lambda/mu}(z,z) en vez de (z,1/z).  Tiene que fallar.
#   D3  SENUELO: sin el signo sgn_t(mu).  Tiene que fallar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_problema2_composicion.py > pI_problema2_composicion_OUT.txt 2>&1

import json
import sys
from collections import Counter

CASOS = [2, 3, 4, 5]
LMAX = 9


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
    """t-core, t-cociente y sgn_t, por el descenso del beta-set (la ruta de [Alb23] §2.2)."""
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
    # el cociente, de las clases de residuo del beta-set ORIGINAL
    orig = sorted(beta_set(lam, n), reverse=True)
    quot = []
    for r in range(t):
        xs = sorted([x for x in orig if x % t == r], reverse=True)
        m = len(xs)
        q = [(xs[k] - r) // t - (m - 1 - k) for k in range(m)]
        quot.append(tuple(x for x in q if x > 0))
    return core, quot, sg


def sub(mu, lam):
    m = list(mu) + [0] * 8
    l = list(lam) + [0] * 8
    return all(m[i] <= l[i] for i in range(8))


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
print("PROBLEMA 2 DEL PAPER I --- .componen phi_t y el par reciproco?")
print("=" * 100)
print("")

resumen = {}
for T in CASOS:
    N = T + 2                      # el alfabeto del paper: t raices + un par
    def hA(k):
        """h_k de mu_t U {z,1/z}: sum sobre los multiplos de t."""
        if k < 0:
            return {}
        acc = {}
        m = 0
        while m * T <= k:
            acc = padd(acc, chi(k - m * T))
            m += 1
        return acc
    def hP(k):
        return chi(k) if k >= 0 else {}
    def schur_directo(lam):
        n = max(1, len(lam))
        l = list(lam) + [0] * n
        return det([[hA(l[i] - (i + 1) + (j + 1)) for j in range(n)] for i in range(n)])
    def skew_par(lam, mu):
        n = max(1, len(lam))
        l = list(lam) + [0] * n
        m = list(mu) + [0] * n
        return det([[hP(l[i] - m[j] - (i + 1) + (j + 1)) for j in range(n)] for i in range(n)])
    def skew_par_zz(lam, mu):
        """senuelo: (z,z) en vez de (z,1/z): h_k = (k+1) z^k."""
        def h2(k):
            return {k: k + 1} if k >= 0 else {}
        n = max(1, len(lam))
        l = list(lam) + [0] * n
        m = list(mu) + [0] * n
        return det([[h2(l[i] - m[j] - (i + 1) + (j + 1)) for j in range(n)] for i in range(n)])
    # ---- C0 -----------------------------------------------------------------------------------
    c0 = 0
    for k in range(0, LMAX + 3):
        c0 += (hA(k) == padd(chi(k), hA(k - T)) if k >= T else hA(k) == chi(k))
    print("  t=%d,  alfabeto de %d letras" % (T, N))
    print("    C0  la recursion h_k(A) = chi_k + h_{k-t}(A) : %d de %d" % (c0, LMAX + 3))
    lams = [l for l in particiones(LMAX, N) if l]
    c1_ok, c1_bad = 0, []
    d1_ok = d2_ok = d3_ok = 0
    sup = Counter()
    c3_ok = 0
    for lam in lams:
        directo = schur_directo(lam)
        comp = {}
        comp_sin_core = {}
        comp_sin_sgn = {}
        comp_zz = {}
        nsup = 0
        for mu in particiones(lam[0], len(lam) + 1):
            if not sub(mu, lam):
                continue
            core, quot, sg = core_quot_sgn(mu, N, T)
            vale = (len(core) == 0) and all(len(q) <= 1 for q in quot)
            sk = skew_par(lam, mu)
            if vale:
                nsup += 1
                comp = padd(comp, sk if sg > 0 else pneg(sk))
                comp_sin_sgn = padd(comp_sin_sgn, sk)
                comp_zz = padd(comp_zz, skew_par_zz(lam, mu) if sg > 0 else pneg(skew_par_zz(lam, mu)))
            comp_sin_core = padd(comp_sin_core, sk if sg > 0 else pneg(sk))
        if comp == directo:
            c1_ok += 1
        else:
            c1_bad.append(lam)
        sup[nsup] += 1
        d1_ok += (comp_sin_core == directo)
        d2_ok += (comp_zz == directo)
        d3_ok += (comp_sin_sgn == directo)
    n = len(lams)
    print("    C1  FATAL la composicion contra el valor directo : %d de %d   %s"
          % (c1_ok, n, "PASA" if (n > 0 and not c1_bad) else "FALLA en %s" % c1_bad[:5]))
    print("    C2  mu supervivientes por lambda : %s" % dict(sorted(sup.items())))
    print("    D1  senuelo sin la condicion de core vacio : %d de %d   %s"
          % (d1_ok, n, "NO DISCRIMINA" if d1_ok == n else "falla en %d (bien)" % (n - d1_ok)))
    print("    D2  senuelo con (z,z) en vez de (z,1/z)    : %d de %d   %s"
          % (d2_ok, n, "NO DISCRIMINA" if d2_ok == n else "falla en %d (bien)" % (n - d2_ok)))
    print("    D3  senuelo sin el signo sgn_t(mu)         : %d de %d   %s"
          % (d3_ok, n, "NO DISCRIMINA" if d3_ok == n else "falla en %d (bien)" % (n - d3_ok)))
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = dict(n=int(n), C1=[int(c1_ok), len(c1_bad)],
                               sup={str(k): int(v) for k, v in sup.items()},
                               D1=int(d1_ok), D2=int(d2_ok), D3=int(d3_ok))

json.dump(resumen, open("pI_problema2_composicion_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
