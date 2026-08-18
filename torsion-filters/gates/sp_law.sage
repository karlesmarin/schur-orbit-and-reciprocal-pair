# -*- coding: utf-8 -*-
# LA LEY DEL PESO MAXIMO: ¿se lee mu_max del defecto de simetria de g_com?  15 de agosto de 2026.
#
# POR QUE ESTE GUION.  sp_expansion.sage contesta el test 1 de la consulta en siete formas y contesta
# QUE SI: mu_max unico y A(mu_max) = +-1 en las siete.  Y con eso pasa algo que no es un matiz:
#
#     las 16 formas que "no tenian certificado" en la presentacion de Laplace SI tienen un unico
#     mu_max con A = +-1.  El fallo del 3,6 % era un artefacto de NUESTRA presentacion, no del objeto.
#
# Asi que el objeto intrinseco es mu_max, y la pregunta que de verdad importa -- la que la consulta
# formula y nosotros no estabamos haciendo -- es:
#
#     ¿PUEDE LEERSE mu_max DEL DEFECTO DE SIMETRIA DE g_com?
#
# Si Conjetura 8.44 es cierta, g_com simetrico  =>  no queda ningun peso.  Y si es asimetrico, tal
# vez no se prediga el caracter entero pero SI el peso mas alto que sobrevive.  Eso bastaria.
#
# LO QUE SE MIDE, sobre la poblacion critica (C = tau y S\g_com simetrico), las DOS columnas
#   N1  ¿es mu_max unico SIEMPRE?  Su test puede fallar: dos dominantes incomparables lo matan.
#   N2  ¿es A(mu_max) = +-1 siempre?  Seria un certificado unitario intrinseco.
#   N3  LA LEY: mu_max contra los invariantes de g_com -- C, min y max de g_com, los reflejos que
#       faltan (C - x para x en g_com), y el rango de S.  Se imprime la tabla cruda, no un ajuste.
#
# CONTROLES
#   C1  DECOY, y es el que da sentido a todo: la columna g_com SIMETRICO tiene que dar CERO pesos.
#       Si diera alguno, o bien la Conjetura 8.44 es falsa en el barrido, o el instrumento miente.
#   C2  el resto del pelado tiene que ser 0 en todas.  Un resto no nulo invalida su fila.
#   C3  n impreso SIEMPRE por columna.  Una columna vacia no dice nada de la otra.
#
# Phi se calcula por Laplace sobre enteros (rapido); sp_expansion.sage ya comprobo contra el
# bialternante que las dos maquinarias dan lo mismo, y ese es el control que autoriza el atajo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sagemath/sagemath sage sp_law.sage

import itertools
from collections import defaultdict

T, R, WMAX = 4, 2, 15
N = T + 2 * R


def perm_sign(seq):
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def clases(beta, t):
    cl = defaultdict(list)
    for i, b in enumerate(beta):
        cl[b % t].append(i)
    return cl if len(cl) == t else None


def atomos(beta, t, r):
    """{exponente: coeficiente} del numerador, por Laplace sobre las t filas congeladas."""
    cl = clases(beta, t)
    if cl is None:
        return None
    keys = sorted(cl)
    idx = tuple(range(2 * r))
    acc = defaultdict(int)
    for pick in itertools.product(*[cl[k] for k in keys]):
        P = sorted(pick)
        Ps = set(P)
        Tt = tuple(beta[i] for i in range(len(beta)) if i not in Ps)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        for S in itertools.combinations(idx, r):
            Sc = tuple(a for a in idx if a not in S)
            q = [0] * (2 * r)
            for j in range(r):
                q[S[j]] = 2 * j
                q[Sc[j]] = 2 * j + 1
            base = w * perm_sign(q)
            for pa in itertools.permutations(range(r)):
                ca = perm_sign(pa)
                ka = [Tt[S[pa[j]]] for j in range(r)]
                for pb in itertools.permutations(range(r)):
                    cb = perm_sign(pb)
                    e = tuple(ka[j] - Tt[Sc[pb[j]]] for j in range(r))
                    acc[e] += base * ca * cb
    return {e: c for e, c in acc.items() if c != 0}


def dividir(A, B, r):
    """A / B exacto como polinomios de Laurent."""
    def sh(P):
        o = tuple(min(e[k] for e in P) for k in range(r))
        return {tuple(e[k] - o[k] for k in range(r)): c for e, c in P.items()}, o
    A2, oa = sh(A)
    B2, ob = sh(B)
    key = lambda e: (sum(e),) + tuple(e)
    lb = max(B2, key=key)
    cb = B2[lb]
    Rm = dict(A2)
    Q = {}
    while Rm:
        lr = max(Rm, key=key)
        d = tuple(lr[k] - lb[k] for k in range(r))
        if min(d) < 0:
            return None
        q = QQ(Rm[lr]) / cb
        Q[d] = Q.get(d, 0) + q
        for e, c in B2.items():
            k2 = tuple(e[k] + d[k] for k in range(r))
            Rm[k2] = Rm.get(k2, 0) - q * c
            if Rm[k2] == 0:
                del Rm[k2]
    off = tuple(oa[k] - ob[k] for k in range(r))
    return {tuple(e[k] + off[k] for k in range(r)): c for e, c in Q.items() if c != 0}


_SP = {}
_W = WeylCharacterRing("C%d" % R)


def sp_char(mu):
    if mu in _SP:
        return _SP[mu]
    el = _W(_W.space().from_vector(vector(list(mu))))
    d = defaultdict(int)
    for wt, m in el.weight_multiplicities().items():
        d[tuple(wt.to_vector())] += m
    _SP[mu] = dict(d)
    return _SP[mu]


def expandir(P):
    P = dict(P)
    salida = []
    for _ in range(300):
        P = {e: c for e, c in P.items() if c != 0}
        if not P:
            return salida, {}
        dom = [e for e in P if list(e) == sorted(e, reverse=True) and min(e) >= 0]
        if not dom:
            return salida, P
        mu = max(dom, key=lambda e: (sum(e), tuple(e)))
        A = P[mu]
        for e, m in sp_char(mu).items():
            P[e] = P.get(e, 0) - A * m
        salida.append((mu, A))
    return salida, P


def betas(t, r, W):
    for mid in itertools.combinations(range(1, W + 1), t + 2 * r - 1):
        yield tuple(sorted(mid, reverse=True)) + (0,)


def datos_MAL(beta, t, r):
    cl = clases(beta, t)
    if cl is None:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    if not E:
        return None
    Cd = {k: sorted((beta[i] for i in cl[k]), reverse=True) for k in E}
    S = sorted({v for k in E for v in Cd[k]})
    incr = []
    for k in E:
        c = sorted(Cd[k], reverse=True)
        incr += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    incr.sort(reverse=True)
    if len(incr) < r:
        return None
    C = S[0] + S[-1]
    tau = incr[r - 1]
    s = set(S)
    gc = sorted(x for x in S if (C - x) not in s)
    return S, C, tau, gc


def datos(beta, t, r):
    """LA DEFINICION DEL PAPER, y la anterior estaba MAL.

    g_com := g_a ∩ g_b -- la parte omitida en comun por las DOS transversales maximizadoras, y
    SOLO esta definida cuando |G| = 2 (eq:gcomgen del paper).  Lo que este guion calculaba antes
    era {x en S : C-x no en S}, el defecto de simetria de S, que es OTRO objeto: en el ejemplo
    trabajado del propio paper, beta=(8,7,6,5,3,2,1,0), el defecto es VACIO y g_com = {3,5}.
    Controlado contra ese ejemplo: g_com={3,5} y K={1,7}, los dos exactos."""
    cl = clases(beta, t)
    if cl is None:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    if not E:
        return None
    Cd = {k: sorted((beta[i] for i in cl[k]), reverse=True) for k in E}
    S = sorted({v for k in E for v in Cd[k]})
    keys = sorted(cl)
    best = None
    G = []
    for pick in itertools.product(*[cl[k] for k in keys]):
        Ps = set(pick)
        Tt = [beta[i] for i in range(len(beta)) if i not in Ps]
        d = sum(Tt[:r]) - sum(Tt[r:])
        if best is None or d > best:
            best = d
            G = [pick]
        elif d == best:
            G.append(pick)
    if len(G) != 2:
        return None                       # g_com solo esta definido con |G| = 2
    sel = [{k: beta[i] for k, i in zip(keys, g)} for g in G]
    gc = sorted({sel[0][k] for k in E} & {sel[1][k] for k in E})
    C = S[0] + S[-1]
    return S, C, gc


DELTA = tuple(range(N - 1, -1, -1))
ND = atomos(DELTA, T, R)

print("=" * 108)
print("LA LEY DEL PESO MAXIMO -- poblacion critica, las dos columnas")
print("=" * 108)
print("")
col = {True: [], False: []}
resto_malo = 0
for b in betas(T, R, WMAX):
    d = datos(b, T, R)
    if d is None:
        continue
    S, C, gc = d
    sim = (set(C - x for x in gc) == set(gc))
    A = atomos(b, T, R)
    P = dividir(A, ND, R) if A else None
    if P is None:
        P = {}
    lista, rr = expandir(P)
    resto_malo += bool(rr)
    lista = [(mu, a) for (mu, a) in lista if a != 0]
    maxi = [mu for (mu, a) in lista
            if not any(nu != mu and all(sum(nu[:k + 1]) >= sum(mu[:k + 1]) for k in range(R))
                       for (nu, _) in lista)]
    col[sim].append((b, S, C, gc, lista, maxi, bool(rr)))

for sim in (True, False):
    L = col[sim]
    et = "g_com SIMETRICO" if sim else "g_com ASIMETRICO  [conj 8.44: debe NO anularse]"
    print("  %-30s n = %d" % (et, len(L)))
    if not L:
        print("     POBLACION VACIA -- esta columna no dice nada")
        print("")
        continue
    conpeso = [x for x in L if x[4]]
    unico = [x for x in conpeso if len(x[5]) == 1]
    unit = [x for x in conpeso if len(x[5]) == 1 and abs(x[4][[m for m, _ in x[4]].index(x[5][0])][1]) == 1]
    print("     con algun peso (Phi != 0) : %d de %d" % (len(conpeso), len(L)))
    print("     SE ANULAN                 : %d" % (len(L) - len(conpeso)))
    print("     mu_max UNICO              : %d de %d" % (len(unico), len(conpeso)))
    print("     |A(mu_max)| == 1          : %d de %d" % (len(unit), len(conpeso)))
    print("     restos no nulos (C2)      : %d   %s"
          % (sum(1 for x in L if x[6]), "ok" if not sum(1 for x in L if x[6]) else "*** REVISAR ***"))
    print("")

print("=" * 108)
print("N3  LA LEY -- mu_max contra los invariantes de g_com.  Tabla cruda, sin ajustar")
print("=" * 108)
print("")
print("   beta                        |  C | g_com      | reflejos que faltan | mu_max   | A")
vistos = 0
for (b, S, C, gc, lista, maxi, rr) in col[False]:
    if not maxi or vistos >= 22:
        continue
    vistos += 1
    A = lista[[m for m, _ in lista].index(maxi[0])][1]
    falt = sorted(C - x for x in gc)
    print("   %-27s | %2d | %-10s | %-19s | %-8s | %s"
          % (str(b), C, str(gc), str(falt), str(maxi[0]), A))
print("")
print("=" * 108)
print("N4  VOLCADO CRUDO para ajuste externo (beta ; lambda ; S ; C ; g_com ; mu_max ; A)")
print("=" * 108)
import json
dump = []
for (b, S, C, gc, lista, maxi, rr) in col[False]:
    if len(maxi) != 1:
        continue
    A = lista[[m for m, _ in lista].index(maxi[0])][1]
    lam = [b[i] - DELTA[i] for i in range(N)]
    dump.append(dict(beta=[int(x) for x in b], lam=[int(x) for x in lam], S=[int(x) for x in S],
                     C=int(C), gcom=[int(x) for x in gc], mu=[int(x) for x in maxi[0]],
                     A=int(A), npesos=int(len(lista))))
open("sp_law_DUMP.json", "w").write(json.dumps(dump))
print("   %d filas volcadas a sp_law_DUMP.json" % len(dump))
print("")
print("=" * 108)
