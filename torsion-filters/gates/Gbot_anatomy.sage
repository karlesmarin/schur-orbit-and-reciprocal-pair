# -*- coding: utf-8 -*-
# LA ANATOMIA DE Gbot.  Arriba el maximizador era casi unico; abajo hay hasta 16.  Que los ordena?
#
# DONDE ESTAMOS.  Phi_t = 0  <=>  [Phi]_top = 0 Y [Phi]_bot = 0 (0 excepciones, 10289 formas, t=2
# incluido).  Arriba: |G| <= 2, probado, y ya da la condicion (ii).  Abajo: [Phi]_bot = 0 es un
# RECUENTO CON SIGNO, sum_{T in Gbot, sort d(T) = D} w(T) = 0 para cada multiconjunto D, probado a
# partir de que la permutacion de bloques de tamano 2 tiene signo +1.  Falta la anatomia de Gbot.
#
# LOS DOS EXTREMOS SON EL MISMO TIPO DE OBJETO, con distinto patron de signos por POSICION:
#       deg(T)    = +..+ - ..-   (bloques)      = sum_{i<=r} u_i - sum_{i>r} u_i
#       degmin(T) = + - + - ...  (alternante)   = sum_i (-1)^{i-1} u_i
# Arriba el patron de bloques se dejaba escribir como max_A (2 sum A - sum T) y de ahi salia la
# separabilidad.  Abajo el signo de cada x depende de la PARIDAD del numero de elementos quitados por
# encima de el: con j(x) = #{g por encima de x},   degmin = sum_{x no en g} (-1)^{p(x)-1+j(x)} x.
# Eso ya NO es separable de la misma manera, y el tamano de Gbot lo confirma.  Pero 16 = 2^4 sugiere
# otra estructura, y es la que se prueba aqui:
#
#   H7   Gbot es un PRODUCTO:  Gbot = prod_k A_k  con A_k subset de la clase k.  Es decir, la
#        eleccion optima en cada clase es independiente de la de las demas.
#   H8   y cada factor tiene |A_k| <= 2, de donde |Gbot| = 2^{#clases empatadas}.
#   H9   cuando |A_k| = 2, los dos elementos son ADYACENTES en la clase, como arriba.
#
# Si H7+H8 valen, la minimizacion es separable en el sentido util —una decision por clase— y el
# recuento con signo se vuelve un producto, que es exactamente lo que hace falta para atacarlo.
#
# MEDIDO, cada cosa capaz de fallar:
#   G1  H7: Gbot es un producto.  Se compara Gbot con prod_k (proyeccion de Gbot a la clase k).
#   G2  H8: max |A_k|.   G3  H9: los dos de un factor son adyacentes en la clase.
#   G4  |Gbot| = prod_k |A_k| = 2^{#empatadas}, y su distribucion.
#   G5  la reflexion: cuando S es concentrico, sigma_C actua sobre Gbot y w(C-T) = -w(T)?  Es el
#       mecanismo de la direccion probada y tiene que verse.
#   G6  contingencia: Gbot cerrado bajo sigma_C  vs  [Phi]_bot = 0.  Arriba la clausura era
#       NECESARIA; aqui hay que mirar, no suponer.
#   G7  SENUELO: la misma pregunta para el estrato de ARRIBA -- alli Gtop tambien deberia ser un
#       producto con factores <= 2, y si sale que si, H7 no es una propiedad del estrato de abajo
#       sino de la construccion, y hay que decirlo.
#   G8  no vacuidad: tiene que haber |A_k| = 2 de verdad, y clases con |A_k| = 1.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

CONF = [(2, 2, 22), (2, 3, 16), (4, 2, 22), (4, 3, 16), (6, 2, 18), (6, 3, 14), (8, 2, 16)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def perm_sign(q):
    n = len(q)
    seen = [False] * n
    s = 1
    for i in range(n):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = q[j]
            ln += 1
        if ln % 2 == 0:
            s = -s
    return s


def analyse(beta, t, r):
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    Es = set(E)
    CC = dict((k, sorted((beta[i] for i in cl[k]), reverse=True)) for k in E)
    S = sorted((b for k in E for b in CC[k]), reverse=True)
    C = S[0] + S[-1]
    conc = set(C - b for b in S) == set(S)
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    cond_ii = (len(fixed) == 2 and all(k in E for k in fixed))
    rows = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        gE = tuple(beta[i] for i in P if beta[i] % t in Es)
        db = sum(T[2 * i] - T[2 * i + 1] for i in range(r))
        dt = sum(T[:r]) - sum(T[r:])
        rows.append((w, T, gE, db, dt))
    Dmin = min(x[3] for x in rows)
    Dmax = max(x[4] for x in rows)
    Gb = [x for x in rows if x[3] == Dmin]
    Gt = [x for x in rows if x[4] == Dmax]
    return dict(conc=conc, cond_ii=cond_ii, crit=conc and cond_ii, E=E, CC=CC, C=C,
                Gb=Gb, Gt=Gt)


def product_test(G, E, CC):
    """is the set of kept excess-values a product over the classes?  returns (is_product, factors)."""
    gs = [dict((v % len(CC) if False else v, None) for v in x[2]) for x in G]
    keep = [set() for _ in E]
    fac = []
    for idx, k in enumerate(E):
        s = set()
        for x in G:
            for v in x[2]:
                if v in CC[k]:
                    s.add(v)
        fac.append(sorted(s, reverse=True))
    prodsize = prod([len(f) for f in fac])
    actual = set(frozenset(x[2]) for x in G)
    isprod = (prodsize == len(actual))
    if isprod:
        for combo in itertools.product(*fac):
            if frozenset(combo) not in actual:
                isprod = False
                break
    return isprod, fac


print("=" * 106)
print("G1-G4  is Gbot a product over the classes, with factors of size <= 2?")
print("=" * 106)
print("")
print("     t   r  shapes | G1 not a product  G2 max|A_k|  G3 non-adjacent  G4 |Gbot| != prod |A_k|"
      " | max|Gbot|")
print("  " + "-" * 102)

TOT = dict(sh=0, g1=0, g3=0, g4=0, two=0, one=0)
MX = 0
MXA = 0
for t, r, MAX in CONF:
    N = t + 2 * r
    nsh = g1 = g3 = g4 = 0
    mxa = mxg = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            nsh += 1
            isprod, fac = product_test(a['Gb'], a['E'], a['CC'])
            if not isprod:
                g1 += 1
            mxa = max(mxa, max(len(f) for f in fac))
            mxg = max(mxg, len(set(frozenset(x[2]) for x in a['Gb'])))
            for idx, k in enumerate(a['E']):
                f = fac[idx]
                if len(f) == 2:
                    TOT['two'] += 1
                    p1, p2 = CC_idx = (a['CC'][k].index(f[0]), a['CC'][k].index(f[1]))
                    if abs(p1 - p2) != 1:
                        g3 += 1
                elif len(f) == 1:
                    TOT['one'] += 1
            if prod([len(f) for f in fac]) != len(set(frozenset(x[2]) for x in a['Gb'])):
                g4 += 1
    MX = max(MX, mxg)
    MXA = max(MXA, mxa)
    TOT['sh'] += nsh
    TOT['g1'] += g1
    TOT['g3'] += g3
    TOT['g4'] += g4
    print("  %4d %3d %7d | %16d %12d %16d %25d | %9d"
          % (t, r, nsh, g1, mxa, g3, g4, mxg))

print("")
print("  G1 shapes where Gbot is NOT a product: %d of %d" % (TOT['g1'], TOT['sh']))
print("  G2 largest factor seen: %d      G4 shapes where |Gbot| != prod|A_k|: %d" % (MXA, TOT['g4']))
print("  G3 factors of size 2 whose members are NOT adjacent in their class: %d" % TOT['g3'])
print("  G8 non-vacuity: factors of size 2 seen %d times, of size 1 %d times."
      % (TOT['two'], TOT['one']))
print("  largest |Gbot| in range: %d" % MX)

# ---------------------------------------------------------------- G5, G6, G7 ---------------------
print("")
print("=" * 106)
print("G5/G6  the reflection on Gbot, and G7 the same product question for the TOP stratum")
print("=" * 106)
print("")
print("     t   r | crit shapes  sigma_C acts on Gbot  w(C-T) = -w(T) | Gbot closed & bot=0"
      "  closed only  bot=0 only | G7 top not a product")
print("  " + "-" * 102)

CT = {}
g7 = 0
for t, r, MAX in CONF:
    N = t + 2 * r
    nc = act = sgn = a11 = a10 = a01 = 0
    g7l = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            C = a['C']
            Ts = set(x[1] for x in a['Gb'])
            W = dict((x[1], x[0]) for x in a['Gb'])
            closed = all(tuple(sorted((C - y for y in T), reverse=True)) in Ts for T in Ts)
            byD = {}
            for w, T, gE, db, dt in a['Gb']:
                d = tuple(sorted((T[2 * i] - T[2 * i + 1] for i in range(r)), reverse=True))
                byD[d] = byD.get(d, 0) + w
            botzero = all(v == 0 for v in byD.values())
            if a['crit']:
                nc += 1
                if closed:
                    act += 1
                    if all(W[tuple(sorted((C - y for y in T), reverse=True))] == -W[T] for T in Ts):
                        sgn += 1
            if closed and botzero:
                a11 += 1
            elif closed:
                a10 += 1
            elif botzero:
                a01 += 1
            isprod, _ = product_test(a['Gt'], a['E'], a['CC'])
            if not isprod:
                g7l += 1
    g7 += g7l
    CT[(t, r)] = (a11, a10, a01)
    print("  %4d %3d | %12d %21d %16d | %19d %13d %12d | %20d"
          % (t, r, nc, act, sgn, a11, a10, a01, g7l))

print("")
print("  G7: shapes where the TOP maximiser set is not a product: %d." % g7)
print("     If this is 0 too, being a product is a property of the CONSTRUCTION (one choice per")
print("     class), not something special about the bottom -- and it must be said that way.")
print("")
print("DONE")
