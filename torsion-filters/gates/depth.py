# -*- coding: utf-8 -*-
# EL ENUNCIADO DE PROFUNDIDAD.  Baja algun estrato de D1 - 4?
#
# LO QUE SE QUIERE.  Todo lo medido hasta ahora dice que el primer grado NO nulo de Phi_t esta en
# D1, D1-2 o D1-4 y nunca mas abajo.  Si eso es un teorema, cierra t >= 4 de un golpe -- en vez de
# un argumento por piso, que es la forma que tendria una regresion infinita.  Aqui se busca el
# contraejemplo.
#
# COMO SE MIDE BARATO.  Laplace de A(T) por las columnas PARES (las de z_j, que son 0,2,..,2r-2):
#
#     A(T) = sum_{A subset T, |A| = r}  eps(A) * a_A(z) * a_{T\A}(1/z),
#     eps(A) = (-1)^{suma de indices de FILA de A}      [el r(r-1) del bloque de columnas es par]
#
# y cada sumando es HOMOGENEO de grado  2*sum(A) - sum(T).  Luego TODO estrato -- no solo el de
# arriba y el segundo -- es una suma de productos de dos alternantes, y el certificado de dimension
# se le aplica igual:
#
#     Delta_D = sum_{g, A : 2 sum A - sum T_g = D}  w(g) eps(A) dim(atil) dim(astar)
#     Delta_D != 0  =>  [Phi_t]_D != 0.
#
# Sano por construccion: es la evaluacion del caracter en z = 1, y (prod z)^d vale 1 ahi, asi que
# el d que cambia de sumando a sumando no estorba (torcer por det^d no cambia la dimension).
# Ojo: es COTA SUPERIOR de la profundidad.  Delta_D = 0 no prueba [Phi]_D = 0, asi que las formas
# donde Delta se anula en los tres primeros pisos hay que recalcularlas EXACTAS.
#
# CONTROLES
#   Y0  VALIDACION: a profundidad 0 el certificado tiene que coincidir con el criterio del estrato
#       de arriba (|G| = 2, mismo INV, signos opuestos).  Si no coincide, el signo eps esta mal y
#       todo lo demas sobra.  Se cuentan los desacuerdos en las DOS direcciones.
#   Y1  distribucion de la profundidad sobre todas las formas.
#   Y2  formas con Delta = 0 en TODOS los pisos hasta D1-12: las candidatas a romper el enunciado,
#       y el UNICO sitio donde puede esconderse un contraejemplo.  Se comprueban TODAS.
#       *** DEFECTO MIO, cazado en la primera pasada: recalculaba solo deep[:400] de 858, y el tope
#       no estaba escrito en ninguna parte.  Un tope silencioso justo sobre las formas que deciden
#       el enunciado.  Levantado.  Y el recalculo exacto (permutaciones) no escalaba, asi que ahora
#       va en dos fases: primero det M sobre GF(p) en puntos aleatorios -- control que PUEDE fallar,
#       porque encontraria un valor no nulo si Phi_t no fuera identicamente cero -- y solo se
#       expande exacto lo que salga no nulo. ***
#   Y3  no vacuidad: n de cada profundidad, y n de las formas con Phi_t == 0 (branch (a)/(b)), que
#       NO cuentan -- ahi no hay "primer grado no nulo".
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python depth.py

import itertools
from collections import Counter, defaultdict
from fractions import Fraction

from second_stratum import setup, all_transversals, inv_of, perm_sign

CONFIGS = [(2, 2, 13), (2, 3, 15), (4, 2, 15), (4, 3, 16), (6, 2, 16),
           (6, 3, 17), (8, 2, 17), (6, 4, 18), (8, 3, 18)]
MAXDEPTH = 12


_DIM = {}


def dim_gl(lam):
    v = _DIM.get(lam)
    if v is not None:
        return v
    r = len(lam)
    num = den = 1
    for i in range(r):
        for j in range(i + 1, r):
            num *= lam[i] - lam[j] + j - i
            den *= j - i
    v = num // den
    _DIM[lam] = v
    return v


def dims(A, B, r):
    """dim(atil) * dim(astar) para el par (A arriba, B abajo), ambos decrecientes."""
    alpha = [A[i] - (r - 1 - i) for i in range(r)]
    atil = tuple(a - alpha[-1] for a in alpha)
    Ls = [B[0] - B[r - 1 - i] for i in range(r)]
    astar = tuple(Ls[i] - (r - 1 - i) for i in range(r))
    return dim_gl(atil) * dim_gl(astar)


def deltas(beta, t, r):
    """{grado: Delta_D} y D1."""
    st = setup(beta, t)
    if st is None:
        return None, None
    cl, E, Cd = st
    if not E:
        return None, None
    tr = all_transversals(beta, cl, r, t)
    D1 = max(x[3] for x in tr)
    acc = defaultdict(int)
    idx = range(2 * r)
    for (_, T, w, _) in tr:
        sT = sum(T)
        for R in itertools.combinations(idx, r):
            A = [T[a] for a in R]
            B = [T[a] for a in idx if a not in R]
            D = 2 * sum(A) - sT
            eps = -1 if (sum(R) % 2) else 1
            acc[D] += w * eps * dims(A, B, r)
    return acc, D1


P = 1000033          # primo con P = 1 (mod 24), asi que hay raices t-esimas para t = 2,4,6,8
assert P % 24 == 1   # (el assert ya cazo un P mal elegido: 1000003 = 19 mod 24)


def _root(t):
    """una raiz t-esima primitiva de la unidad en GF(P)."""
    for g in range(2, 200):
        z = pow(g, (P - 1) // t, P)
        if all(pow(z, d, P) != 1 for d in range(1, t)):
            return z
    raise RuntimeError


def det_modp(rows):
    n = len(rows)
    m = [r[:] for r in rows]
    det = 1
    for c in range(n):
        p = next((i for i in range(c, n) if m[i][c]), None)
        if p is None:
            return 0
        if p != c:
            m[c], m[p] = m[p], m[c]
            det = -det
        det = det * m[c][c] % P
        inv = pow(m[c][c], P - 2, P)
        for i in range(c + 1, n):
            f = m[i][c] * inv % P
            if f:
                for j in range(c, n):
                    m[i][j] = (m[i][j] - f * m[c][j]) % P
    return det % P


def nonzero_at_random(beta, t, r, tries=6, seed=20260812):
    """True si det M != 0 en algun punto aleatorio => Phi_t no es identicamente cero.
    Es un control que PUEDE fallar: si Phi_t no fuera cero, aqui saldria un valor no nulo."""
    import random
    rnd = random.Random(seed + sum(beta))
    zeta = _root(t)
    N = len(beta)
    for _ in range(tries):
        zs = [rnd.randrange(2, P - 1) for _ in range(r)]
        alpha = [pow(zeta, k, P) for k in range(t)]
        for z in zs:
            alpha.append(z)
            alpha.append(pow(z, P - 2, P))
        rows = [[pow(a, b, P) for a in alpha] for b in beta]
        if det_modp(rows):
            return True
    return False


def exact_top_degree(beta, t, r):
    """primer grado NO nulo de Phi_t, exacto (expansion completa).  None si Phi_t == 0."""
    st = setup(beta, t)
    cl, E, Cd = st
    tr = all_transversals(beta, cl, r, t)
    acc = defaultdict(int)
    n = 2 * r
    for (_, T, w, _) in tr:
        for s in itertools.permutations(range(n)):
            e = [0] * r
            for a in range(n):
                b = s[a]
                e[b // 2] += T[a] if b % 2 == 0 else -T[a]
            acc[tuple(e)] += w * perm_sign(s)
    nz = [sum(k) for k, v in acc.items() if v]
    return max(nz) if nz else None


def main():
    y0_bad = y0_n = 0
    depth = Counter()
    deep = []
    phizero = 0
    tot = 0

    for (t, r, M) in CONFIGS:
        N = t + 2 * r
        for comb in itertools.combinations(range(M + 1), N):
            beta = tuple(sorted(comb, reverse=True))
            acc, D1 = deltas(beta, t, r)
            if acc is None:
                continue
            tot += 1

            # Y0: validacion del signo contra el criterio del estrato de arriba
            st = setup(beta, t)
            cl, E, Cd = st
            trv = all_transversals(beta, cl, r, t)
            G = [x for x in trv if x[3] == D1]
            crit_top_zero = (len(G) == 2 and inv_of(G[0][1], r) == inv_of(G[1][1], r)
                             and G[0][2] == -G[1][2])
            y0_n += 1
            if crit_top_zero != (acc.get(D1, 0) == 0):
                y0_bad += 1

            d = None
            for k in range(0, MAXDEPTH + 1, 2):
                if acc.get(D1 - k, 0) != 0:
                    d = k
                    break
            if d is None:
                deep.append((t, r, beta, D1))
                phizero += 1
            else:
                depth[d] += 1
        print("   hecho t=%d r=%d M=%d  (acumulado %d formas)" % (t, r, M, tot), flush=True)

    print("")
    print("Y0 VALIDACION del signo: certificado a profundidad 0  vs  criterio del estrato de arriba")
    print("     desacuerdos: %d de %d" % (y0_bad, y0_n))
    print("")
    print("Y1 PROFUNDIDAD (D1 menos el primer grado con Delta != 0), sobre %d formas:" % tot)
    for k in sorted(depth):
        print("      D1 - %-3d : %d" % (k, depth[k]))
    print("")
    print("Y2 formas donde Delta se anula en TODOS los pisos hasta D1-%d : %d" % (MAXDEPTH, phizero))
    if deep:
        print("   se comprueban TODAS (el certificado solo da cota superior).  Fase 1: det M sobre")
        print("   GF(%d) en 6 puntos aleatorios; fase 2: expansion exacta de las que salgan != 0." % P)
        real = Counter()
        alive = []
        for (t, r, beta, D1) in deep:
            if nonzero_at_random(beta, t, r):
                alive.append((t, r, beta, D1))
            else:
                real["Phi == 0"] += 1
        print("      Phi_t identicamente cero (branch (a)/(b)) : %d de %d" % (real["Phi == 0"],
                                                                             len(deep)))
        print("      con algun valor NO nulo, o sea Phi_t != 0 : %d" % len(alive))
        for (t, r, beta, D1) in alive:
            dt = exact_top_degree(beta, t, r)
            d = "Phi == 0 (contradice la fase 1!)" if dt is None else D1 - dt
            print("         t=%d r=%d beta=%s  PROFUNDIDAD REAL %s" % (t, r, list(beta), d))
        if not alive:
            print("      *** ninguna forma con Delta = 0 hasta D1-%d tiene Phi_t != 0 ***"
                  % MAXDEPTH)


if __name__ == "__main__":
    main()
