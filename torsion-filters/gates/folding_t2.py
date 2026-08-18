# -*- coding: utf-8 -*-
# ============================================================================================
#  EL PLEGADO EN t = 2, VERSION 2.  14 de agosto de 2026.
#
#  HISTORIA (la v1 esta archivada en folding_t2_OUT.txt, y su cabecera en el commit anterior):
#    - v1 predijo  s_lambda(1,-1,z^{+-1}) = +- monomio * UN caracter simplectico.  FALSA: 33 de 62.
#    - la correccion "combinacion entera NO NEGATIVA" (folding_t2_refino.py) tambien FALLA:
#      145 de 290, y las 145 que fallan son MENOS una combinacion no negativa.  Reparto 145/145.
#    - consulta externa del 13: confirma que la factorizacion en un solo caracter NO esta publicada
#      y NO PUEDE EXISTIR, con el mismo contraejemplo que ya habia salido aqui, y con la familia
#      cerrada  s_(m)(1,-1,z,1/z) = sum_{j<=m/2} sp_(m-2j)  via  1/((1-u^2)(1-zu)(1-u/z)).
#      Ruta correcta: branching GL_N -> O_N (Littlewood; Frohmader Thm 5.6 fuera del rango estable)
#      y DESPUES plegado D_{r+1} -> C_r (Jantzen; Kumar-Lusztig-Prasad).
#
#  LO QUE ESTE GUION AÑADE, y que no estaba en la consulta: DOS TEOREMAS GRATIS que fijan el signo.
#
#    (T1)  prod(alfabeto) = 1 * (-1) * prod(z_i/z_i) = -1,  luego  s_{lambda + 1^N} = - s_lambda.
#          El signo global TIENE que invertirse al sumar 1^N.  Como N = 2r+2 es PAR, (-1)^{|lambda|}
#          no se invierte: queda REFUTADO por algebra, sin medir.  (Explica el 145/290 = azar.)
#    (T2)  el alfabeto es cerrado por inversion, luego s_lambda(u^{-1}) = s_lambda(u), y por
#          complementacion en el rectangulo k^N:   s_{lambda*} = (-1)^k s_lambda.
#
#  Y LA REGLA DE SIGNO, que resulta NO ser desconocida.  Laplace por las dos columnas de paridad
#  (C1' = (C1+C2)/2 = [beta_i par], C2' = (C1-C2)/2 = [beta_i impar], det(M) = -2 det(M')) da
#
#      A(beta) = -2 * sum_{p<q, paridad opuesta} (-1)^{p+q+1} * (+1 si beta_p par else -1)
#                       * det_pares(las 2r filas restantes)
#
#  El termino de grado maximo lo aporta el par {p,q} = {1,2} (0-indexado), y de ahi:
#
#      r = 1, si beta_1 != beta_2 (mod 2):   signo = +  <=>  beta_1 PAR.
#
#  Que es el signo del ESTRATO DE ARRIBA -- el mismo epsilon del maximizador de deg(T) = sumH - sumL
#  que ya esta en el paper con |G| <= 2 PROBADO.  El plegado no abre una regla nueva: devuelve al
#  cuadro de transversales.  Esta es la prediccion que N2 contrasta, con su caso no generico como
#  denominador.
#
#  COLUMNAS
#    C0  ACEPTACION, fatal.  a) valores conocidos.  b) sympy independiente sobre una muestra.
#        c) LA FAMILIA CERRADA de la consulta externa (prediccion externa, no nuestra).
#        d) dim sp_mu por la formula de Weyl de C_r.  e) invariancia de Weyl de la restriccion.
#        f) T1 y T2 medidos (son teoremas: si fallan, el guion esta mal, no el teorema).
#    N1  LA DICOTOMIA: s_lambda(1,-1,z^{+-1}) = +- (combinacion NO NEGATIVA de caracteres de Sp(2r)).
#        Se cuenta todo-no-negativo / todo-no-positivo / MEZCLADO.  Un solo MEZCLADO la refuta.
#    N2  LA REGLA DE SIGNO.  Prediccion combinatoria desde beta contra el signo MEDIDO, con
#        Contingencia (que se niega a certificar si la casilla contrafactual esta vacia).
#    N3  COTAS: la consulta dice que no hay ninguna.  Se mide el crecimiento de los coeficientes y
#        se verifica el caso cerrado [z^0] s_(2M)(1,-1,z,1/z) = M+1.
#    N4  EL CRITERIO t=2 (cero <=> (i) y (ii)) por esta via completamente independiente -- sin
#        Laplace por estratos, sin probe(), sin scan() -- y ahora tambien para r = 2 y r = 3.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python folding_t2.py        (desde gates/)
# ============================================================================================

import json
import os
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, permutations, product

from _control import Contingencia

OUT_JSON = "folding_t2_RESULT.json"

# ---------------------------------------------------------------- Laurent en r variables --------
# Un polinomio es un dict {(e_1,...,e_r): coef entero}.  Sin sympy, sin division de fracciones.


def padd(a, b):
    o = dict(a)
    for e, c in b.items():
        v = o.get(e, 0) + c
        if v:
            o[e] = v
        elif e in o:
            del o[e]
    return o


def pscale(a, k):
    if k == 0:
        return {}
    return {e: c * k for e, c in a.items()}


def psub(a, b):
    return padd(a, pscale(b, -1))


def lex_lead(a):
    """monomio maximo en orden lex (los exponentes se comparan como tuplas)."""
    return max(a)


def ldiv(num, den):
    """division exacta de Laurent por Laurent.  Falla ruidosamente si no divide."""
    if not num:
        return {}
    q = {}
    rem = dict(num)
    le = lex_lead(den)
    lc = den[le]
    guard = 0
    while rem:
        guard += 1
        if guard > 200000:
            raise RuntimeError("ldiv no termina: no es division exacta")
        me = lex_lead(rem)
        mc = rem[me]
        if mc % lc:
            raise RuntimeError("ldiv: coeficiente no divisible (%r / %r)" % (mc, lc))
        f = mc // lc
        ke = tuple(x - y for x, y in zip(me, le))
        q[ke] = q.get(ke, 0) + f
        for e, c in den.items():
            t = tuple(x + y for x, y in zip(e, ke))
            v = rem.get(t, 0) - f * c
            if v:
                rem[t] = v
            elif t in rem:
                del rem[t]
    return {e: c for e, c in q.items() if c}


def _sgn(p):
    s = 1
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                s = -s
    return s


_PERMS = {}


def perms_signed(n):
    if n not in _PERMS:
        _PERMS[n] = [(p, _sgn(p)) for p in permutations(range(n))]
    return _PERMS[n]


def det_pares(filas, r):
    """det de la matriz 2r x 2r con filas indexadas por 'filas' (exponentes beta) y columnas
    (z_1^b, z_1^{-b}, ..., z_r^b, z_r^{-b}).  Cada entrada es un MONOMIO: el determinante es una
    suma sobre permutaciones de sumas de exponentes, sin una sola multiplicacion de polinomios."""
    n = 2 * r
    out = {}
    for p, s in perms_signed(n):
        e = [0] * r
        for i in range(n):
            c = p[i]
            a, lo = c >> 1, c & 1
            e[a] += -filas[i] if lo else filas[i]
        k = tuple(e)
        v = out.get(k, 0) + s
        if v:
            out[k] = v
        elif k in out:
            del out[k]
    return out


def alt(beta, r):
    """A(beta) = det(u_j^{beta_i}) con u = (1, -1, z_1, 1/z_1, ..., z_r, 1/z_r).

    Por Laplace sobre las dos primeras columnas tras C1'=(C1+C2)/2, C2'=(C1-C2)/2 (indicadores de
    paridad), que hace det(M) = -2 det(M').  Solo sobreviven los pares de filas de paridad OPUESTA:
    esas dos clases mod 2 son exactamente las del criterio."""
    N = 2 * r + 2
    tot = {}
    for p in range(N):
        for q in range(p + 1, N):
            pp, pq = beta[p] & 1, beta[q] & 1
            if pp == pq:
                continue
            m2 = 1 if pp == 0 else -1            # E_p O_q - O_p E_q
            s = -1 if (p + q + 1) & 1 else 1
            resto = [beta[i] for i in range(N) if i != p and i != q]
            sub = det_pares(resto, r)
            if sub:
                tot = padd(tot, pscale(sub, s * m2))
    return pscale(tot, -2)


_DEN = {}


def denominador(r):
    if r not in _DEN:
        N = 2 * r + 2
        _DEN[r] = alt([N - 1 - i for i in range(N)], r)
    return _DEN[r]


def restriccion(lam, r):
    """s_lambda(1, -1, z_1^{+-1}, ..., z_r^{+-1}) exacta, como Laurent entero."""
    N = 2 * r + 2
    beta = [lam[i] + (N - 1 - i) for i in range(N)]
    A = alt(beta, r)
    if not A:
        return {}, beta
    return ldiv(A, denominador(r)), beta


# ---------------------------------------------------------------- caracteres de Sp(2r) ----------
_SPC = {}


def sp_char(mu, r):
    """caracter irreducible de Sp(2r) de peso dominante mu, por la formula de Weyl:
    det(z_i^{a_j} - z_i^{-a_j}) / det(z_i^{r-j+1} - z_i^{-(r-j+1)}),  a_j = mu_j + r - j + 1."""
    key = (tuple(mu), r)
    if key in _SPC:
        return _SPC[key]
    a = [mu[j] + r - j for j in range(r)]
    num = _det_binom(a, r)
    den = _det_binom([r - j for j in range(r)], r)
    out = ldiv(num, den)
    _SPC[key] = out
    return out


def _det_binom(a, r):
    """det de (i,j) |-> z_i^{a_j} - z_i^{-a_j}."""
    out = {}
    for p, s in perms_signed(r):
        for ch in product((1, -1), repeat=r):
            e = [0] * r
            c = s
            for i in range(r):
                e[i] = ch[i] * a[p[i]]
                if ch[i] == -1:
                    c = -c
            k = tuple(e)
            v = out.get(k, 0) + c
            if v:
                out[k] = v
            elif k in out:
                del out[k]
    return out


def dom(e):
    """representante dominante del peso e para C_r: valores absolutos, ordenados decreciente."""
    return tuple(sorted((abs(x) for x in e), reverse=True))


def to_sp(P, r):
    """expansion en la base {sp_mu}.  Devuelve {mu: coef entero}.  Greedy por el peso dominante
    maximo (orden (|mu|, lex)), que es correcto porque sp_mu = z^mu + pesos estrictamente menores."""
    coef = {}
    rem = dict(P)
    guard = 0
    while rem:
        guard += 1
        if guard > 20000:
            raise RuntimeError("to_sp no termina")
        best = max((dom(e) for e in rem), key=lambda d: (sum(d), d))
        c = rem.get(best)
        if not c:
            raise RuntimeError("to_sp: el polinomio no es invariante de Weyl (falta z^%r)" % (best,))
        coef[best] = coef.get(best, 0) + c
        rem = psub(rem, pscale(sp_char(list(best), r), c))
    return {m: c for m, c in coef.items() if c}


def weyl_invariante(P, r):
    """comprueba invariancia bajo z_i <-> 1/z_i y bajo permutar las z_i."""
    for a in range(r):
        for e, c in P.items():
            f = list(e)
            f[a] = -f[a]
            if P.get(tuple(f), 0) != c:
                return False
    for e, c in P.items():
        if P.get(tuple(sorted(e, reverse=True)), 0) == 0 and c:
            pass
    for p, _ in perms_signed(r):
        for e, c in P.items():
            f = tuple(e[p[i]] for i in range(r))
            if P.get(f, 0) != c:
                return False
    return True


# ---------------------------------------------------------------- el criterio t = 2 -------------
def condiciones(lam, t=2):
    """None = rama (a) (alguna clase vacia).  Si no, (i) and (ii)."""
    N = len(lam)
    beta = [lam[i] + (N - 1 - i) for i in range(N)]
    cl = {}
    for b in beta:
        cl.setdefault(b % t, []).append(b)
    if len(cl) < t:
        return None
    E = [k for k in cl if len(cl[k]) >= 2]
    if not E:
        return None
    S = sorted({v for k in E for v in cl[k]})
    C = S[0] + S[-1]
    i_ok = set(C - v for v in S) == set(S)
    sols = [k for k in range(t) if (2 * k - C) % t == 0]
    return i_ok and len(sols) == 2 and all(k in E for k in sols)


def signo_predicho(beta, r):
    """LA REGLA, deducida del Laplace y escrita ANTES de correr nada:  el termino de grado maximo
    de A(beta) lo aporta el par de filas de paridad opuesta {p,q} que deja el resto mas "abierto";
    el generico es {p,q} = {1,2} y entonces el signo es (+1 si beta[1] es par else -1).
    Devuelve (signo, generico) -- generico False marca el caso que la regla NO cubre."""
    if (beta[1] & 1) != (beta[2] & 1):
        return (1 if (beta[1] & 1) == 0 else -1), True
    return None, False


def signo_de(coefs):
    """signo global de una expansion: +1 si todo >= 0, -1 si todo <= 0, 0 si MEZCLADA."""
    pos = any(c > 0 for c in coefs.values())
    neg = any(c < 0 for c in coefs.values())
    if pos and neg:
        return 0
    return 1 if pos else (-1 if neg else 1)


def particiones(N, lmax):
    for c in combinations_with_replacement(range(lmax, -1, -1), N):
        yield list(c)


# ===================================================================== C0 ========================
print("=" * 112)
print("C0  ACEPTACION -- fatal.  Si algo de aqui falla, NADA de lo de abajo vale.")
print("=" * 112)
print("")
FALLOS = 0

# --- C0a  valores conocidos ---------------------------------------------------------------------
print("  C0a  valores conocidos del bialternante")
casos = [
    (1, [0, 0, 0, 0], {(0,): 1}),                       # s_vacia = 1
    (1, [1, 0, 0, 0], {(1,): 1, (-1,): 1}),             # e_1 = 1 - 1 + z + 1/z
    (1, [1, 1, 1, 1], {(0,): -1}),                      # e_4 = producto del alfabeto = -1
    (1, [1, 1, 0, 0], {}),                              # e_2 = 0
    (2, [0] * 6, {(0, 0): 1}),
    (2, [1] + [0] * 5, {(1, 0): 1, (-1, 0): 1, (0, 1): 1, (0, -1): 1}),
    (2, [1] * 6, {(0, 0): -1}),
]
for r, lam, esp in casos:
    got, _ = restriccion(lam, r)
    ok = got == esp
    FALLOS += not ok
    print("       r=%d  s_%-14s -> %-46s %s" % (r, str(tuple(lam)), str(got)[:46],
                                                "ok" if ok else "*** FALLA (esperado %s) ***" % esp))
print("")

# --- C0b  sympy independiente -------------------------------------------------------------------
print("  C0b  contraste contra una implementacion INDEPENDIENTE (sympy, cociente de alternantes)")
try:
    import sympy as _sp

    def _sympy_restr(lam, r):
        N = 2 * r + 2
        zs = list(_sp.symbols('w1:%d' % (r + 1))) if r > 1 else [_sp.Symbol('w1')]
        u = [_sp.Integer(1), _sp.Integer(-1)]
        for z in zs:
            u += [z, 1 / z]
        beta = [lam[i] + (N - 1 - i) for i in range(N)]
        num = _sp.Matrix(N, N, lambda i, j: u[j] ** beta[i])
        den = _sp.Matrix(N, N, lambda i, j: u[j] ** (N - 1 - i))
        return _sp.simplify(_sp.cancel(_sp.together(num.det() / den.det()))), zs

    def _a_dict(expr, zs, r):
        expr = _sp.expand(_sp.simplify(expr))
        if expr == 0:
            return {}
        d = {}
        desp = 1
        for z in zs:
            desp *= z ** 30
        p = _sp.Poly(_sp.expand(_sp.cancel(expr * desp)), *zs)
        for mon, c in zip(p.monoms(), p.coeffs()):
            d[tuple(int(m) - 30 for m in mon)] = int(c)
        return {k: v for k, v in d.items() if v}

    muestra = [(1, l) for l in [[3, 1, 0, 0], [4, 2, 1, 0], [5, 5, 2, 0], [2, 2, 2, 1], [6, 3, 3, 1]]]
    malos = 0
    for r, lam in muestra:
        mio, _ = restriccion(lam, r)
        e, zs = _sympy_restr(lam, r)
        suyo = _a_dict(e, zs, r)
        ok = mio == suyo
        malos += not ok
        print("       r=%d  lambda=%-18s %s" % (r, str(tuple(lam)), "ok" if ok else
                                                "*** DISCREPA ***  mio=%s  sympy=%s" % (mio, suyo)))
    FALLOS += malos
    print("       %d discrepancias sobre %d" % (malos, len(muestra)))
except ImportError:
    print("       sympy no disponible -- C0b SALTADO (y eso DEBILITA la aceptacion)")
print("")

# --- C0b2  LA DEFINICION: suma sobre tableaux semiestandar ---------------------------------------
# sympy con N = 6 y 8 tarda decenas de minutos (fue lo que colgo la v1), asi que para r >= 2 el
# contraste independiente es la DEFINICION combinatoria de s_lambda, que es lo mas independiente
# que hay: no usa alternantes, ni Laplace, ni division.
print("  C0b2 contraste contra LA DEFINICION: s_lambda = suma sobre SSYT con entradas en el alfabeto")


def ssyt_restr(lam, r):
    N = 2 * r + 2
    filas = [x for x in lam if x > 0]
    if not filas:
        return {tuple([0] * r): 1}
    coords = [(i, j) for i in range(len(filas)) for j in range(filas[i])]
    fill = {}
    out = {}

    def rec(k):
        if k == len(coords):
            cnt = [0] * N
            for v in fill.values():
                cnt[v] += 1
            e = tuple(cnt[2 + 2 * a] - cnt[3 + 2 * a] for a in range(r))
            c = -1 if cnt[1] & 1 else 1
            v = out.get(e, 0) + c
            if v:
                out[e] = v
            elif e in out:
                del out[e]
            return
        i, j = coords[k]
        lo = 0
        if j > 0:
            lo = max(lo, fill[(i, j - 1)])
        if i > 0:
            lo = max(lo, fill[(i - 1, j)] + 1)
        for v in range(lo, N):
            fill[(i, j)] = v
            rec(k + 1)

    rec(0)
    return out


malos = 0
muestra2 = [(1, [3, 1, 0, 0]), (1, [4, 4, 2, 1]), (1, [2, 2, 0, 0]),
            (2, [2, 1, 0, 0, 0, 0]), (2, [3, 2, 1, 1, 0, 0]), (2, [2, 2, 2, 1, 1, 0]),
            (2, [1, 1, 1, 1, 0, 0]), (3, [2, 1, 0, 0, 0, 0, 0, 0]), (3, [1, 1, 1, 1, 1, 1, 0, 0]),
            (3, [2, 2, 1, 1, 0, 0, 0, 0])]
for r, lam in muestra2:
    mio, _ = restriccion(lam, r)
    suyo = ssyt_restr(lam, r)
    ok = mio == suyo
    malos += not ok
    print("       r=%d  lambda=%-26s %s" % (r, str(tuple(lam)), "ok" if ok else
                                            "*** DISCREPA ***  mio=%s  ssyt=%s" % (mio, suyo)))
FALLOS += malos
print("       %d discrepancias sobre %d" % (malos, len(muestra2)))
print("")

# --- C0b3  evaluacion numerica exacta ------------------------------------------------------------
print("  C0b3 evaluacion NUMERICA exacta (Fraction): determinante por eliminacion gaussiana en")
print("       puntos z racionales, contra mi polinomio de Laurent evaluado ahi.  Cubre r = 4.")


def det_frac(M):
    n = len(M)
    M = [row[:] for row in M]
    d = Fraction(1)
    for c in range(n):
        p = None
        for i in range(c, n):
            if M[i][c]:
                p = i
                break
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        inv = Fraction(1) / M[c][c]
        for i in range(c + 1, n):
            if M[i][c]:
                f = M[i][c] * inv
                for j in range(c, n):
                    M[i][j] -= f * M[c][j]
    return d


def num_restr(lam, r, zv):
    N = 2 * r + 2
    u = [Fraction(1), Fraction(-1)]
    for z in zv:
        u += [z, Fraction(1) / z]
    beta = [lam[i] + (N - 1 - i) for i in range(N)]
    num = det_frac([[u[j] ** beta[i] for j in range(N)] for i in range(N)])
    den = det_frac([[u[j] ** (N - 1 - i) for j in range(N)] for i in range(N)])
    return num / den


def eval_laurent(P, zv):
    tot = Fraction(0)
    for e, c in P.items():
        v = Fraction(c)
        for a, ea in enumerate(e):
            v *= zv[a] ** ea
        tot += v
    return tot


PUNTOS = [[Fraction(2), Fraction(3), Fraction(5), Fraction(7)],
          [Fraction(3, 2), Fraction(5, 2), Fraction(7, 3), Fraction(11, 4)]]
malos = 0
n3 = 0
muestra3 = [(1, [5, 3, 1, 0]), (2, [3, 2, 2, 1, 0, 0]), (2, [4, 3, 1, 1, 1, 0]),
            (3, [2, 2, 1, 1, 1, 0, 0, 0]), (3, [3, 1, 1, 0, 0, 0, 0, 0]),
            (4, [2, 1, 1, 0, 0, 0, 0, 0, 0, 0]), (4, [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]),
            (4, [2, 2, 1, 1, 1, 1, 0, 0, 0, 0])]
for r, lam in muestra3:
    mio, _ = restriccion(lam, r)
    for pt in PUNTOS:
        zv = pt[:r]
        a = eval_laurent(mio, zv)
        b = num_restr(lam, r, zv)
        n3 += 1
        if a != b:
            malos += 1
            print("       r=%d  lambda=%-26s *** DISCREPA ***  mio=%s  num=%s"
                  % (r, str(tuple(lam)), a, b))
    print("       r=%d  lambda=%-26s ok en %d puntos" % (r, str(tuple(lam)), len(PUNTOS)))
FALLOS += malos
print("       %d discrepancias sobre %d evaluaciones" % (malos, n3))
print("")

# --- C0c  la familia cerrada de la consulta externa ----------------------------------------------
print("  C0c  PREDICCION EXTERNA (no nuestra):  s_(m)(1,-1,z,1/z) = sum_{j<=m/2} sp_(m-2j)")
print("       viene de  sum_m h_m u^m = 1/((1-u^2)(1-zu)(1-u/z))  y  1/((1-zu)(1-u/z)) = sum sp_(n) u^n")
malos = 0
for m in range(0, 15):
    P, _ = restriccion([m, 0, 0, 0], 1)
    got = to_sp(P, 1)
    esp = {(m - 2 * j,): 1 for j in range(m // 2 + 1)}
    ok = got == esp
    malos += not ok
    if m <= 6 or not ok:
        print("       m=%-3d -> %-44s %s" % (m, str(sorted(got.items(), reverse=True))[:44],
                                             "ok" if ok else "*** FALLA (esperado %s) ***" % esp))
print("       %d fallos sobre 15 valores de m" % malos)
FALLOS += malos
print("")

# --- C0d  dimensiones de Sp(2r) por la formula de Weyl -------------------------------------------
print("  C0d  sp_mu(1,...,1) contra la formula de dimension de Weyl para C_r")


def dim_weyl_C(mu, r):
    a = [Fraction(mu[j] + r - j) for j in range(r)]
    b = [Fraction(r - j) for j in range(r)]
    num, den = Fraction(1), Fraction(1)
    for i in range(r):
        num *= a[i]
        den *= b[i]
        for j in range(i + 1, r):
            num *= a[i] ** 2 - a[j] ** 2
            den *= b[i] ** 2 - b[j] ** 2
    return num / den


malos = 0
pruebas = [((3,), 1), ((0,), 1), ((7,), 1), ((2, 1), 2), ((3, 3), 2), ((1, 0), 2),
           ((2, 1, 1), 3), ((1, 1, 1), 3)]
for mu, r in pruebas:
    val = sum(sp_char(list(mu), r).values())
    esp = dim_weyl_C(list(mu), r)
    ok = Fraction(val) == esp
    malos += not ok
    print("       sp_%-10s r=%d :  suma de coeficientes = %-8d  Weyl = %-8s %s"
          % (str(mu), r, val, esp, "ok" if ok else "*** FALLA ***"))
FALLOS += malos
print("")

# --- C0e / C0f  invariancia de Weyl, y los dos teoremas gratis -----------------------------------
print("  C0e  invariancia de Weyl de la restriccion (z_i <-> 1/z_i y permutar las z_i)")
print("  C0f  T1:  s_{lambda + 1^N} = - s_lambda     (porque prod(alfabeto) = -1)")
print("       T2:  s_{lambda*}      = (-1)^k s_lambda,  lambda* = complemento en el rectangulo k^N")
malos_e = malos_1 = malos_2 = 0
n_e = n_1 = n_2 = 0
for r, LM in [(1, 8), (2, 5), (3, 3)]:
    N = 2 * r + 2
    for lam in particiones(N, LM):
        P, _ = restriccion(lam, r)
        n_e += 1
        if P and not weyl_invariante(P, r):
            malos_e += 1
        # T1
        Q, _ = restriccion([x + 1 for x in lam], r)
        n_1 += 1
        if Q != pscale(P, -1):
            malos_1 += 1
        # T2, con k = lambda_1
        k = lam[0]
        lstar = [k - lam[N - 1 - i] for i in range(N)]
        R, _ = restriccion(lstar, r)
        n_2 += 1
        if R != pscale(P, (-1) ** k):
            malos_2 += 1
print("       C0e  invariancia de Weyl : %d fallos sobre %d" % (malos_e, n_e))
print("       C0f  T1                  : %d fallos sobre %d" % (malos_1, n_1))
print("       C0f  T2                  : %d fallos sobre %d" % (malos_2, n_2))
FALLOS += malos_e + malos_1 + malos_2
print("")
print("       COROLARIO DE T1, SIN MEDIR NADA: N = 2r+2 es PAR, luego (-1)^{|lambda|} NO cambia al")
print("       sumar 1^N mientras que el signo SI.  (-1)^{|lambda|} queda REFUTADA como regla de")
print("       signo por algebra -- y eso explica su 145/290 de ayer, que era exactamente el azar.")
print("")
if FALLOS:
    print("  C0 FALLA (%d fallos) -- veredicto SUSPENDIDO." % FALLOS)
    print("DONE (suspendido)")
    raise SystemExit(1)
print("  C0 PASA")
sys.stdout.flush()

RES = {}
CFG = [(1, 14), (2, 8), (3, 5), (4, 2)]

# ===================================================================== N1, N2, N4 ================
print("")
print("=" * 112)
print("N1  LA DICOTOMIA  |  N2  LA REGLA DE SIGNO  |  N4  el criterio, por via independiente")
print("=" * 112)

GLOBAL = Counter()
for r, LM in CFG:
    N = 2 * r + 2
    print("")
    print("  r=%d  (N=%d, plegado D_%d -> C_%d = Sp(%d), lambda_1 <= %d)" % (r, N, r + 1, r, 2 * r, LM))
    cnt = Counter()
    cont_sig = Contingencia("signo MEDIDO = +", "signo PREDICHO de beta = +")
    cont_cri = Contingencia("s_lambda == 0", "(i) y (ii)")
    mezcladas = []
    no_generico = 0
    maxcoef = 0
    maxsuma = 0
    ejemplos = []
    n = 0
    for lam in particiones(N, LM):
        P, beta = restriccion(lam, r)
        n += 1
        cero = not P
        crit = condiciones(lam)
        if crit is not None:
            cont_cri.add(cero, crit, tuple(lam))
        if cero:
            cnt["cero"] += 1
            continue
        coefs = to_sp(P, r)
        s = signo_de(coefs)
        if s == 0:
            cnt["MEZCLADA"] += 1
            if len(mezcladas) < 6:
                mezcladas.append((tuple(lam), sorted(coefs.items(), reverse=True)))
        else:
            cnt["no negativa" if s > 0 else "no positiva"] += 1
            pred, gen = signo_predicho(beta, r)
            if gen:
                cont_sig.add(s > 0, pred > 0, tuple(lam))
            else:
                no_generico += 1
        maxcoef = max(maxcoef, max(abs(c) for c in coefs.values()))
        maxsuma = max(maxsuma, sum(abs(c) for c in coefs.values()))
        if len(ejemplos) < 8 and len(coefs) > 1:
            ejemplos.append((tuple(lam), sorted(coefs.items(), reverse=True)[:5], s))
    print("     formas: %d" % n)
    print("     ejemplos (lambda -> multiplicidades en la base sp_mu, hasta 5 terminos):")
    for lam, c, s in ejemplos:
        print("        %-22s %-56s signo %+d" % (str(lam), str(c)[:56], s))
    print("")
    print("     N1  DICOTOMIA")
    for k in ("cero", "no negativa", "no positiva", "MEZCLADA"):
        print("        %-14s : %d" % (k, cnt[k]))
    print("        -> %s" % ("*** REFUTADA: %d expansiones con signos MEZCLADOS ***" % cnt["MEZCLADA"]
                             if cnt["MEZCLADA"] else
                             "AGUANTA: toda expansion no nula es +- una combinacion no negativa"))
    for lam, c in mezcladas:
        print("           contraejemplo %-20s %s" % (str(lam), str(c)[:70]))
    print("")
    print("     N2  REGLA DE SIGNO  (generico = beta_1 y beta_2 de paridad distinta)")
    print("        fuera del caso generico (la regla NO opina): %d" % no_generico)
    cont_sig.informe(indent="        ")
    print("")
    print("     N4  CRITERIO  cero <=> (i) y (ii)   [via bialternante, sin estratos ni probe()]")
    cont_cri.informe(indent="        ")
    print("")
    print("     N3  coeficiente maximo |a_mu| = %d ;  suma de |a_mu| maxima = %d" % (maxcoef, maxsuma))
    RES["r%d" % r] = {"n": n, "cnt": dict(cnt), "no_generico": no_generico,
                      "desac_signo": cont_sig.desacuerdos(), "signo_valido": cont_sig.valido(),
                      "desac_criterio": cont_cri.desacuerdos(), "criterio_valido": cont_cri.valido(),
                      "maxcoef": maxcoef, "maxsuma": maxsuma}
    GLOBAL["mezcladas"] += cnt["MEZCLADA"]
    GLOBAL["desac_signo"] += cont_sig.desacuerdos()
    GLOBAL["desac_criterio"] += cont_cri.desacuerdos()
    sys.stdout.flush()

# ===================================================================== N3 ========================
print("")
print("=" * 112)
print("N3  COTAS: la consulta externa dice que NO hay ninguna.  Se verifica el caso cerrado.")
print("=" * 112)
print("")
print("     [z^0] s_(2M)(1,-1,z,1/z) tiene que valer M+1, y s_(2M)(1,-1,1,1) = (M+1)^2")
malos = 0
for M in range(0, 9):
    P, _ = restriccion([2 * M, 0, 0, 0], 1)
    c0 = P.get((0,), 0)
    val = sum(P.values())
    ok = (c0 == M + 1) and (val == (M + 1) ** 2)
    malos += not ok
    print("        M=%-3d  [z^0] = %-4d (esperado %-4d)   s(1,-1,1,1) = %-6d (esperado %-6d)  %s"
          % (M, c0, M + 1, val, (M + 1) ** 2, "ok" if ok else "*** FALLA ***"))
print("        %d fallos sobre 9" % malos)
print("")
print("     -> no hay analogo de {0,+-1,+-2}: los coeficientes crecen sin cota ya en r = 1.")
RES["cotas_fallos"] = malos

# ===================================================================== VEREDICTO =================
print("")
print("=" * 112)
print("VEREDICTO")
print("=" * 112)
print("")
print("  1. MUERTO Y ENTERRADO: 's_lambda(1,-1,z^{+-1}) = +- monomio * UN caracter simplectico'.")
print("     Contraejemplo minimo propio y externo coincidentes: s_(2)(1,-1,z,1/z) = sp_(2) + sp_(0).")
print("")
print("  2. LO QUE QUEDA EN PIE: +- combinacion NO NEGATIVA de caracteres de Sp(2r).")
print("     expansiones con signos MEZCLADOS en total: %d  -> %s"
      % (GLOBAL["mezcladas"], "REFUTADA" if GLOBAL["mezcladas"] else "MEDIDA, sin excepcion"))
print("")
print("  3. LA REGLA DE SIGNO NO ERA DESCONOCIDA: es el signo del ESTRATO DE ARRIBA.")
print("     desacuerdos totales prediccion-vs-medida: %d" % GLOBAL["desac_signo"])
print("     Y (-1)^{|lambda|} esta REFUTADA por algebra (T1), no por conteo.")
print("")
print("  4. EL CRITERIO t=2 por via independiente: %d desacuerdos en total." % GLOBAL["desac_criterio"])
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1)
print("DONE")
