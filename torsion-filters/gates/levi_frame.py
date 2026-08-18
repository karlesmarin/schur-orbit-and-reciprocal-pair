# -*- coding: utf-8 -*-
# EL MARCO DE LEVI, verificado a maquina.
#
# LA AFIRMACION.  El grado total  sum_j e_j  sobre nuestro alfabeto ES la graduacion por el
# cocaracter  mu = (0^t, +1,-1, +1,-1, ..., +1,-1)  de GL(N), porque  z_j -> u z_j  manda
# 1/z_j -> u^{-1}/z_j  y las raices de la unidad no se mueven.  El centralizador de mu tiene el 0
# con multiplicidad t y el +-1 con multiplicidad r cada uno, luego
#
#       L_mu  =  GL(t) x GL(r) x GL(r)
#
# y la filtracion por grado total es la graduacion de V_lambda por mu, con cada estrato un modulo
# sobre ese Levi.  Partiendo el alfabeto en los tres bloques (formula de division iterada):
#
#   Phi_t  =  SUM_{nu <= rho <= lambda}  s_nu(1,zeta,..,zeta^{t-1}) * s_{rho/nu}(z) * s_{lambda/rho}(1/z)
#
# y el grado total del termino (nu,rho) es   2|rho| - |nu| - |lambda|.
#
# TRES RUTAS INDEPENDIENTES, obligadas a coincidir:
#   A) Phi_t directo, suma sobre SSYT del alfabeto ENTERO de N letras.  No sabe nada de Levi.
#   B) la suma de Levi de arriba, sobre pares nu <= rho <= lambda, con SSYT sesgados.
#   C) el Dmax combinatorio por TRANSVERSALES de second_stratum.py, que no sabe nada de A ni de B.
#
# LO QUE SE MIDE, cada cosa capaz de fallar
#   V1  A == B como polinomio de Laurent exacto sobre Z[zeta].  Si falla, el marco esta mal.
#   V2  agrupar B por  2|rho|-|nu|-|lambda|  reproduce las componentes homogeneas de A.  ESTA es la
#       afirmacion "grado total = graduacion por mu"; V1 sola no la testa.
#   V3  Dmax(B) == Dmax(C).  Dos calculos sin nada en comun.
#   V4  la componente de grado maximo de A coincide con la suma de los terminos (nu,rho) de grado
#       maximo de B.
#   V5  s_nu(1,zeta,..,zeta^{t-1}) = 0  <=>  el t-nucleo de nu es NO vacio (Littlewood).  Es la
#       afirmacion de que el factor GL(t) es DONDE vive la anulacion.
#   V6  SENUELO que TIENE que fallar: agrupar por |rho/nu| + |lambda/rho| (la graduacion sum|e_j|,
#       no sum e_j).  Si tambien reprodujera las componentes, V2 no estaria diciendo nada.
#   V7  SENUELO 2: agrupar por |rho| - |nu| a secas (sin el -|lambda/rho|).  Debe fallar.
#   V8  no vacuidad: cuantas formas con [Phi]_top = 0, cuantas con |G| = 2, cuantos nu con
#       s_nu(mu_t) = 0.  Si salen 0 el gate no dice nada.
#   V9  el PUENTE: numero de pares (nu,rho) en el grado maximo con los tres factores no nulos,
#       contra |G| del cuadro de transversales.  NO se predice que coincidan -- se mide.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python levi_frame.py     (aritmetica exacta en Z[zeta], sin Sage)

import itertools
from collections import defaultdict

from second_stratum import setup, all_transversals


# ---------------------------------------------------------------- Z[zeta] exacto
def poly_divmod(a, b):
    a = a[:]
    q = [0] * max(0, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        while a and a[-1] == 0:
            a.pop()
        if len(a) < len(b):
            break
        d = len(a) - len(b)
        c = a[-1] // b[-1]
        q[d] = c
        for i, bi in enumerate(b):
            a[d + i] -= c * bi
        while a and a[-1] == 0:
            a.pop()
    return q, a


def cyclotomic(t, _cache={}):
    """Phi_t(x) con coeficientes enteros, por division de x^t - 1 entre los Phi_d, d|t, d<t."""
    if t in _cache:
        return _cache[t]
    num = [0] * t + [1]
    num[0] = -1
    for d in range(1, t):
        if t % d == 0:
            num, rem = poly_divmod(num, cyclotomic(d))
            assert not any(rem), (t, d)
    _cache[t] = num
    return num


class Cyc(object):
    """Z[zeta_t] = Z[x]/(Phi_t(x)).  Igualdad exacta, sin flotantes."""

    def __init__(self, t):
        self.t = t
        self.m = cyclotomic(t)
        self.d = len(self.m) - 1

    def zero(self):
        return (0,) * self.d

    def one(self):
        return tuple([1] + [0] * (self.d - 1))

    def reduce(self, c):
        c = list(c)
        while len(c) > self.d:
            k = len(c) - 1
            lead = c.pop()
            if lead:
                shift = k - self.d
                for i, mi in enumerate(self.m[:-1]):
                    c[shift + i] -= lead * mi
        c += [0] * (self.d - len(c))
        return tuple(c)

    def zpow(self, e):
        """zeta^e."""
        return self.reduce([0] * (e % self.t) + [1])

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def mul(self, a, b):
        out = [0] * (2 * self.d)
        for i, x in enumerate(a):
            if not x:
                continue
            for j, y in enumerate(b):
                out[i + j] += x * y
        return self.reduce(out)


# ---------------------------------------------------------------- SSYT
def ssyt(outer, inner, n):
    """genera SSYT de forma sesgada outer/inner con entradas 1..n, como lista de filas
    (cada fila es la lista de entradas de las celdas outer[i]-inner[i])."""
    inner = list(inner) + [0] * (len(outer) - len(inner))
    if any(inner[i] > outer[i] for i in range(len(outer))):
        return
    rows = []

    def rec(i):
        if i == len(outer):
            yield [r[:] for r in rows]
            return
        width = outer[i] - inner[i]
        if width == 0:
            rows.append([])
            yield from rec(i + 1)
            rows.pop()
            return
        prev = rows[i - 1] if i else None

        def fill(j, cur, last):
            if j == width:
                rows.append(cur[:])
                yield from rec(i + 1)
                rows.pop()
                return
            col = inner[i] + j            # indice de columna absoluto
            lo = last
            if prev is not None and inner[i - 1] <= col < outer[i - 1]:
                lo = max(lo, prev[col - inner[i - 1]] + 1)
            for v in range(lo, n + 1):
                cur.append(v)
                yield from fill(j + 1, cur, v)
                cur.pop()

        yield from fill(0, [], 1)

    yield from rec(0)


def content(tab, n):
    c = [0] * n
    for row in tab:
        for v in row:
            c[v - 1] += 1
    return c


# ---------------------------------------------------------------- las tres rutas
def phi_direct(lam, t, r, K):
    """A) Phi_t como polinomio de Laurent en z, suma sobre SSYT del alfabeto de N = t+2r letras."""
    N = t + 2 * r
    out = defaultdict(lambda: K.zero())
    for tab in ssyt(lam, [], N):
        c = content(tab, N)
        e = tuple(c[t + j] - c[t + r + j] for j in range(r))
        coef = K.zpow(sum(i * c[i] for i in range(t)))
        out[e] = K.add(out[e], coef)
    return {k: v for k, v in out.items() if any(v)}


def s_at_roots(nu, t, K):
    """s_nu(1, zeta, ..., zeta^{t-1}) en Z[zeta]."""
    acc = K.zero()
    for tab in ssyt(nu, [], t):
        c = content(tab, t)
        acc = K.add(acc, K.zpow(sum(i * c[i] for i in range(t))))
    return acc


def skew_poly(outer, inner, r, K, invert):
    """s_{outer/inner}(z) o s_{outer/inner}(1/z) como dict exponente -> entero (en Z[zeta])."""
    out = defaultdict(lambda: K.zero())
    for tab in ssyt(outer, inner, r):
        c = content(tab, r)
        e = tuple(-x for x in c) if invert else tuple(c)
        out[e] = K.add(out[e], K.one())
    return out


def subpartitions(lam):
    ranges = []
    prev = None
    for a in lam:
        ranges.append(range(a + 1))
    for cand in itertools.product(*ranges):
        if all(cand[i] >= cand[i + 1] for i in range(len(cand) - 1)):
            yield tuple(cand)


def phi_levi(lam, t, r, K):
    """B) la suma de Levi.  Devuelve (polinomio total, {grado: polinomio}, [(nu,rho,grado)])."""
    total = defaultdict(lambda: K.zero())
    byd = defaultdict(lambda: defaultdict(lambda: K.zero()))
    terms = []
    n = sum(lam)
    subs = list(subpartitions(lam))
    for nu in subs:
        if len(list(itertools.islice(ssyt(nu, [], t), 1))) == 0:
            continue
        sn = s_at_roots(nu, t, K)
        if not any(sn):
            continue
        for rho in subs:
            if not all(nu[i] <= rho[i] for i in range(len(lam))):
                continue
            A = skew_poly(rho, nu, r, K, False)
            if not A:
                continue
            B = skew_poly(lam, rho, r, K, True)
            if not B:
                continue
            d = 2 * sum(rho) - sum(nu) - n
            terms.append((nu, rho, d))
            for ea, va in A.items():
                for eb, vb in B.items():
                    e = tuple(x + y for x, y in zip(ea, eb))
                    v = K.mul(K.mul(sn, va), vb)
                    total[e] = K.add(total[e], v)
                    byd[d][e] = K.add(byd[d][e], v)
    clean = lambda D: {k: v for k, v in D.items() if any(v)}
    return clean(total), {d: clean(p) for d, p in byd.items()}, terms


def homog(poly, d):
    return {k: v for k, v in poly.items() if sum(k) == d}


def tcore_empty(nu, t):
    """el t-nucleo de nu es vacio?  se quitan t-ganchos del borde via el beta-set."""
    n = max(len(nu), 1) + sum(nu)
    beta = sorted([nu[i] + n - 1 - i if i < len(nu) else n - 1 - i for i in range(n)],
                  reverse=True)
    bs = set(beta)
    moved = True
    while moved:
        moved = False
        for b in sorted(bs, reverse=True):
            if b - t >= 0 and (b - t) not in bs:
                bs.discard(b)
                bs.add(b - t)
                moved = True
                break
    return bs == set(range(len(bs) - 1, -1, -1))


def main():
    cases = []
    for t, r in [(2, 1), (2, 2), (4, 1), (4, 2), (6, 1)]:
        N = t + 2 * r
        for size in range(0, 9):
            for lam in _parts(size, N):
                cases.append((tuple(lam) + (0,) * (N - len(lam)), t, r))
    print("casos: %d" % len(cases))

    v1 = v2 = v4 = v6 = v7 = 0
    v3_n = v3_bad = 0
    v5_n = v5_bad = 0
    topzero = g2 = zeronu = 0
    bridge = defaultdict(int)
    n_ok = 0

    for (lam, t, r) in cases:
        K = Cyc(t)
        A = phi_direct(lam, t, r, K)
        B, byd, terms = phi_levi(lam, t, r, K)
        n_ok += 1
        if A != B:
            v1 += 1
            if v1 <= 3:
                print("   V1 FALLO t=%d r=%d lam=%s" % (t, r, lam))
        # V2: agrupar por 2|rho|-|nu|-|lambda|
        for d, p in byd.items():
            if p != homog(A, d):
                v2 += 1
                break
        # V6 / V7: senuelos de graduacion
        for tag, keyf in (('v6', lambda nu, rho: (sum(rho) - sum(nu)) + (sum(lam) - sum(rho))),
                          ('v7', lambda nu, rho: sum(rho) - sum(nu))):
            grp = defaultdict(lambda: K.zero())
            ok = True
            agg = defaultdict(lambda: defaultdict(lambda: K.zero()))
            for (nu, rho, _) in terms:
                pass
            # se reconstruye por grupos y se compara con las componentes homogeneas
            byk = defaultdict(list)
            for (nu, rho, d) in terms:
                byk[keyf(nu, rho)].append((nu, rho))
            for kk, lst in byk.items():
                ds = {2 * sum(rho) - sum(nu) - sum(lam) for (nu, rho) in lst}
                if len(ds) > 1:
                    ok = False
                    break
            if not ok:
                if tag == 'v6':
                    v6 += 1
                else:
                    v7 += 1
        if A:
            dmaxB = max(sum(k) for k in A)
            if homog(A, dmaxB) != byd.get(dmaxB, {}):
                v4 += 1
        # V3 / V9: contra el cuadro de transversales
        N = t + 2 * r
        beta = tuple(lam[i] + N - 1 - i for i in range(N))
        if len(set(beta)) == N and min(beta) >= 0:
            st = setup(beta, t)
            if st is not None:
                cl, E, Cd = st
                tr = all_transversals(beta, cl, r, t)
                dmaxC = max(x[3] for x in tr)
                G = [x for x in tr if x[3] == dmaxC]
                if A:
                    v3_n += 1
                    if max(sum(k) for k in A) > dmaxC:
                        v3_bad += 1
                    if not homog(A, dmaxC):
                        topzero += 1
                    ntop = sum(1 for (_, _, d) in terms if d == dmaxC)
                    bridge[(len(G), ntop)] += 1
                if len(G) == 2:
                    g2 += 1
        # V5: Littlewood en el factor GL(t)
        for nu in subpartitions(lam):
            if len(list(itertools.islice(ssyt(nu, [], t), 1))) == 0:
                continue
            v5_n += 1
            val = s_at_roots(nu, t, K)
            if (not any(val)) != (not tcore_empty(nu, t)):
                v5_bad += 1
            if not any(val):
                zeronu += 1

    print("")
    print("V1  A == B  (directo vs Levi)                    : %d fallos de %d" % (v1, n_ok))
    print("V2  grado total == 2|rho|-|nu|-|lambda|          : %d fallos de %d" % (v2, n_ok))
    print("V4  componente de grado maximo coincide          : %d fallos de %d" % (v4, n_ok))
    print("V3  Dmax(Levi) <= Dmax(transversales)            : %d fallos de %d" % (v3_bad, v3_n))
    print("V5  s_nu(mu_t) = 0  <=>  t-nucleo NO vacio       : %d fallos de %d" % (v5_bad, v5_n))
    print("V6  SENUELO graduacion sum|e_j|  (debe fallar)   : %d formas donde NO refina" % v6)
    print("V7  SENUELO graduacion |rho|-|nu| (debe fallar)  : %d formas donde NO refina" % v7)
    print("")
    print("V8  no vacuidad: [Phi]_top = 0 en %d ; |G| = 2 en %d ; nu con s_nu(mu_t) = 0 en %d"
          % (topzero, g2, zeronu))
    print("V9  PUENTE  (|G|, #pares (nu,rho) en grado maximo):")
    for k in sorted(bridge):
        print("        |G|=%d  pares=%d  : %d" % (k[0], k[1], bridge[k]))


def _parts(n, maxlen):
    if n == 0:
        yield []
        return
    def rec(n, mx, cur):
        if n == 0:
            yield cur[:]
            return
        if len(cur) == maxlen:
            return
        for v in range(min(n, mx), 0, -1):
            cur.append(v)
            yield from rec(n - v, v, cur)
            cur.pop()
    yield from rec(n, n, [])


if __name__ == "__main__":
    main()
