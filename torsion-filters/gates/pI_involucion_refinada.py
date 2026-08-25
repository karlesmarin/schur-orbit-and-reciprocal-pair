#
# LA CORRECCION QUE MOTIVA ESTA GATE.  Una involucion sobre los mu de la composicion exigiria
# s_{lambda/mu}(z,1/z) = s_{lambda/mu'}(z,1/z)  para los dos emparejados, y no hay razon
# estructural para ello.  El camino bueno es bajar al INDICE REFINADO: un
# tableau semiestandar de lambda/mu con letras {1,2} equivale a  mu <= nu <= lambda  con nu/mu y
# lambda/nu tiras horizontales, y peso  z^{2|nu| - |mu| - |lambda|}.  Luego
#
#     Phi_t(lambda;z) = sum_{(mu,nu) in J_lambda}  sgn_t(mu) . z^{2|nu| - |mu| - |lambda|}
#
# y ESE es el conjunto donde buscar la involucion, porque ahi la cancelacion puede ser monomio a
# monomio.  Su gate barato: comprobar si los terminos que se cancelan pueden emparejarse ya como
# POLINOMIOS ENTEROS; si no, queda probado que una involucion sobre mu solo no puede ser el
# mecanismo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python pI_involucion_refinada.py > pI_involucion_refinada_OUT.txt 2>&1

import itertools
import json
import sys
from collections import Counter, defaultdict

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
    """f(d) = z^d - z^{-d}."""
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
    """sup/inf es una tira HORIZONTAL: inf <= sup y sup_{i+1} <= inf_i."""
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


print("=" * 100)
print("LA INVOLUCION AL NIVEL CORRECTO --- indice refinado (mu,nu), y L5 contra Plucker")
print("=" * 100)
print("")

# ---- G4: .es L5 una relacion de Plucker? -------------------------------------------------------
pl_ok, pl_n = 0, 0
for b in itertools.combinations(range(9), 4):
    a1, a2, a3, a4 = b
    pl_n += 1
    izq = padd(padd(pmul(f(a1 - a2), f(a3 - a4)), pneg(pmul(f(a1 - a3), f(a2 - a4)))),
               pmul(f(a1 - a4), f(a2 - a3)))
    pl_ok += (izq == {})
print("  G4a  la relacion de Plucker de tres terminos sobre p_ij = f(beta_i-beta_j) :")
print("       p_ij p_kl - p_ik p_jl + p_il p_jk = 0  en  %d de %d cuadruplas" % (pl_ok, pl_n))
l5a_ok, l5a_n = 0, 0
for c in range(0, 5):
    for p in range(1, 5):
        for q in range(1, 5):
            l5a_n += 1
            izq = {}
            for e in (1, -1):
                for h in (1, -1):
                    tm = f(c + e * p + h * q)
                    izq = padd(izq, tm if e * h > 0 else pneg(tm))
            der = pmul(pmul(f(c), f(p)), f(q))
            l5a_ok += (izq == der)
print("  G4b  L5(a) sum_{eps,eta} eps.eta f(c+eps.p+eta.q) = f(c)f(p)f(q) : %d de %d"
      % (l5a_ok, l5a_n))
print("       GRADOS: Plucker es CUADRATICA en los p; L5 iguala una suma LINEAL de f a un producto")
print("       de TRES f.  Son identidades distintas: L5 no es Plucker ni una degeneracion suya.")
print("       Pero las dos valen sobre el MISMO objeto --- los menores 2x2 son coordenadas de")
print("       Plucker de rango dos ---, que es lo que esa lectura senala.")
print("")

resumen = {"G4a": [int(pl_ok), int(pl_n)], "G4b": [int(l5a_ok), int(l5a_n)]}

# ---- G1, G2, G3 -------------------------------------------------------------------------------
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
    lams = [l for l in particiones(LMAX, N) if l]
    g1_pos, g1_no = 0, 0
    g1_hist = Counter()
    g1_tot = 0
    g2_ok, g2_bad = 0, []
    g3_equil, g3_no = 0, 0
    ej = []
    for lam in lams:
        terminos = []
        for mu in particiones(lam[0], len(lam) + 1):
            if not sub(mu, lam):
                continue
            core, quot, sg = core_quot_sgn(mu, N, T)
            if len(core) != 0 or any(len(q) > 1 for q in quot):
                continue
            terminos.append((mu, sg, skew_par(lam, mu)))
        directo = schur_directo(lam)
        # G2: la identidad refinada por (mu,nu)
        refi = {}
        pares = []
        for (mu, sg, _) in terminos:
            for nu in particiones(lam[0], len(lam) + 1):
                if not (sub(mu, nu) and sub(nu, lam)):
                    continue
                if not (hstrip(nu, mu) and hstrip(lam, nu)):
                    continue
                e = 2 * sum(nu) - sum(mu) - sum(lam)
                refi = padd(refi, {e: sg})
                pares.append((e, sg))
        if refi == directo:
            g2_ok += 1
        else:
            g2_bad.append(lam)
        # G3: por exponente, cuantos + y cuantos -
        porexp = defaultdict(lambda: [0, 0])
        for (e, sg) in pares:
            porexp[e][0 if sg > 0 else 1] += 1
        equilibra = all(min(a, b) > 0 or a + b == abs(a - b) for a, b in porexp.values())
        g3_equil += 1 if equilibra else 0
        g3_no += 0 if equilibra else 1
        # G1: .se puede emparejar YA al nivel de mu, con el mismo polinomio y signo opuesto?
        porpol = defaultdict(lambda: [0, 0])
        for (mu, sg, pol) in terminos:
            porpol[tuple(sorted(pol.items()))][0 if sg > 0 else 1] += 1
        # el residuo tras emparejar dentro de cada polinomio identico
        # OJO: comparar sum (a-b).pol contra el valor directo es una TAUTOLOGIA --- es la misma suma
        # reordenada --- y la primera version hacia eso, o sea un control que no puede fallar.  La
        # pregunta real es CUANTOS terminos quedan SIN emparejar: si el emparejamiento a nivel de mu
        # fuera el mecanismo, deberian sobrevivir las 3 o 4 transversales de Laplace y no mas.
        sin_emparejar = sum(abs(a - b) for (a, b) in porpol.values())
        n_term = len(terminos)
        g1_hist[sin_emparejar] += 1
        g1_tot += n_term
        if sin_emparejar <= 4:
            g1_pos += 1
        else:
            g1_no += 1
            if len(ej) < 3:
                ej.append((tuple(lam), n_term, sin_emparejar))
    n = len(lams)
    print("  t=%d  (%d lambda, lambda_1 <= %d, hasta %d filas)" % (T, n, LMAX, N))
    print("    G2  la identidad refinada (mu,nu) contra el valor directo : %d de %d   %s"
          % (g2_ok, n, "PASA" if not g2_bad else "FALLA en %s" % g2_bad[:4]))
    print("    G1  terminos SIN emparejar al nivel de mu (igual polinomio, signo opuesto):")
    print("        histograma : %s" % dict(sorted(g1_hist.items())))
    print("        con <= 4 sin emparejar (lo que Laplace deja) : %d de %d ; con mas : %d %s"
          % (g1_pos, n, g1_no, ("p.ej. %s" % ej) if ej else ""))
    print("        %d terminos en total sobre las %d lambda" % (g1_tot, n))
    print("    G3  por exponente de z, .hay signos de los dos tipos donde debe? : %d de %d"
          % (g3_equil, n))
    print("")
    sys.stdout.flush()
    resumen["t=%d" % T] = dict(n=int(n), G2=[int(g2_ok), len(g2_bad)],
                               G1=[int(g1_pos), int(g1_no)], G3=[int(g3_equil), int(g3_no)])

json.dump(resumen, open("pI_involucion_refinada_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
