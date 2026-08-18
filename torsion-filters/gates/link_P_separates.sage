# -*- coding: utf-8 -*-
# EL ESLABON QUE QUEDA: [Phi]_top = 0 => T_B = C - T_A.  Se reduce a "P separa", y eso se testea.
#
# ESTADO.  proof_G_le_2 probo |G| <= 2 y que los dos maximizadores difieren en UN intercambio entre
# las clases k y k' = k + t/2, con
#       T_A = K u {a1, a2},   T_B = K u {b1, b2},   a1 = c_{k,j},   b1 = c_{k,j+1},
#                                                   b2 = c_{k',j'}, a2 = c_{k',j'+1},
# y el empate Delta_k(j) = Delta_{k'}(j') = V dice exactamente
#       a1 + b1 = a2 + b2 = V.
# O sea sigma_V YA intercambia los dos elementos que difieren, sin hipotesis ninguna.  Luego
#       T_B = V - T_A   <=>   K = V - K,
# y el eslabon entero se reduce a la simetria de la parte COMUN.  Y lemma_V_eq_C prueba que si la
# reflexion cruza los dos entonces V = C, asi que basta trabajar con V.
#
# [Phi]_top = w_A P(T_A) + w_B P(T_B) con w = +-1 y P(T) = a_H(z) a_L(1/z).  Luego
#       [Phi]_top = 0   <=>   P(T_A) = +- P(T_B).
# Y P tiene DOS invariancias exactas, las dos de una linea:
#       P(T + m) = P(T)      porque a_{H+m}(z) = (prod z_j^m) a_H(z) y a_{L+m}(1/z) lo cancela;
#       P(c - T) = P(T)      el lema de reflexion, con (-1)^{r(r-1)} = +1.
# Asi que P es constante en la orbita de T bajo el grupo afin generado por traslaciones y la
# reflexion.  LA PREGUNTA, y es la unica que queda:
#
#       A1   P(T) = +- P(T')  =>  T' esta en la orbita de T ?
#
# Si la respuesta es SI, el eslabon cae: P(T_A) = +-P(T_B) obliga a T_B = T_A + m o T_B = c - T_A, y
# el caso traslacion hay que descartarlo aparte (A3).  Si la respuesta es NO, hay coincidencias de
# Schur y el eslabon necesita otra idea -- y el contraejemplo mas pequeno es el objeto a mirar.
#
# LO QUE SE MIDE, cada cosa capaz de fallar:
#   A1  BUSQUEDA EXHAUSTIVA de colisiones: todos los T de tamano 2r en [0,M], agrupados por P.  Se
#       reporta cuantos grupos NO son una sola orbita.  Es la pregunta de arriba, sin rodeos.
#   A2  no vacuidad: las orbitas de tamano 2 (T no simetrico) deben existir y deben detectarse como
#       colision, o el test no esta mirando nada.
#   A3  sobre las formas reales con |G| = 2: contingencia entre [Phi]_top = 0, T_B = V - T_A,
#       T_B = T_A + m, y K = V - K.  Si la traslacion nunca ocurre, el eslabon cae del todo.
#   A4  control forzado: a1 + b1 = a2 + b2 = V en TODA forma con |G| = 2.  Sale del empate; si
#       fallara una sola vez, la descripcion de los dos maximizadores estaria mal.
#   A5  el signo: cuando T_B = V - T_A, w_B debe ser -w_A.  Es la otra mitad del eslabon y esta
#       medida, no probada; se cuenta aparte.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

# ------------------------------------------------------------------ P and its orbit --------------


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


def P_of(T, r):
    """a_H(z) a_L(1/z) as {exponent tuple: coeff}; T given in DECREASING order."""
    D = {}
    n = 2 * r
    for a in itertools.permutations(range(r)):
        for b in itertools.permutations(range(r)):
            q = [0] * n
            e = [0] * r
            for i in range(r):
                q[i] = 2 * a[i]
                e[a[i]] += T[i]
            for i in range(r):
                q[r + i] = 2 * b[i] + 1
                e[b[i]] -= T[r + i]
            k = tuple(e)
            D[k] = D.get(k, 0) + perm_sign(q)
    return dict((k, v) for k, v in D.items() if v != 0)


def key_of(T, r):
    """canonical key of P(T), up to global sign."""
    D = P_of(T, r)
    a = tuple(sorted(D.items()))
    b = tuple(sorted((k, -v) for k, v in D.items()))
    return min(a, b)


def orbit(T):
    """T normalised to min 0, and its reflection, both as decreasing tuples."""
    u = sorted(T, reverse=True)
    m = u[-1]
    n1 = tuple(x - m for x in u)
    M = u[0]
    n2 = tuple(sorted((M - x for x in u), reverse=True))
    m2 = n2[-1]
    n2 = tuple(x - m2 for x in n2)
    return frozenset([n1, n2])


# ------------------------------------------------------------------ A1, A2 -----------------------
print("=" * 104)
print("A1  exhaustive collision search: does P(T) = +- P(T') force T' into the orbit of T?")
print("=" * 104)
print("")
print("     r   2r   M   sets tested   distinct P   groups   groups that are ONE orbit   BAD groups")
print("  " + "-" * 100)

BADEX = []
A2ok = 0
for r, M in ((2, 15), (3, 12), (4, 10)):
    n = 2 * r
    if M + 1 < n:
        continue
    buckets = {}
    cnt = 0
    for T in itertools.combinations(range(M + 1), n):
        Td = tuple(sorted(T, reverse=True))
        if Td[-1] != 0:
            continue                       # normalise translations away: keep min = 0
        cnt += 1
        buckets.setdefault(key_of(Td, r), []).append(Td)
    good = bad = 0
    twoorbits = 0
    for k, v in buckets.items():
        orbs = set()
        for T in v:
            orbs.add(orbit(T))
        if len(orbs) == 1:
            good += 1
            if len(list(orbs)[0]) == 2:
                twoorbits += 1
        else:
            bad += 1
            if len(BADEX) < 8:
                BADEX.append((r, sorted(v)[:4]))
    A2ok += twoorbits
    print("  %4d %4d %3d %13d %12d %8d %27d %12d"
          % (r, n, M, cnt, len(buckets), len(buckets), good, bad))

print("")
print("  A2 non-vacuity: groups whose single orbit has TWO members (T not self-reflective) and were")
print("     therefore a real collision that the test had to absorb: %d.  Must be nonzero." % A2ok)
if BADEX:
    print("")
    print("  BAD groups -- P collides across different orbits.  The link needs another idea:")
    for r, v in BADEX:
        print("     r=%d : %s" % (r, ", ".join(str(list(x)) for x in v)))
else:
    print("")
    print("  NO bad group in the whole search: P separates 2r-sets exactly up to translation and")
    print("  reflection.  Then P(T_A) = +-P(T_B) forces T_B = T_A + m or T_B = c - T_A.")

# ------------------------------------------------------------------ A3, A4, A5 -------------------
CONF = [(4, 1, 30), (4, 2, 24), (4, 3, 18), (6, 2, 18), (6, 3, 14), (8, 2, 16), (10, 2, 18)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


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
    best = None
    G = {}
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        Pp = sorted(pick)
        Ps = set(Pp)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        d = sum(T[:r]) - sum(T[r:])
        if best is None or d > best:
            best, G = d, {}
        if d == best:
            w = perm_sign([beta[i] % t for i in Pp])
            if sum(Pp) % 2:
                w = -w
            G[frozenset(beta[i] for i in Pp if beta[i] % t in Es)] = (T, w)
    inc = []
    for k in E:
        ck = CC[k]
        for j in range(1, len(ck)):
            inc.append((ck[j - 1] + ck[j], k))
    inc.sort(key=lambda z: -z[0])
    V = inc[r - 1][0]
    S = sorted((b for k in E for b in CC[k]), reverse=True)
    return G, V, S[0] + S[-1]


print("")
print("=" * 104)
print("A3/A4/A5  the real shapes with |G| = 2")
print("=" * 104)
print("")
print("     t   r  |G|=2 | A4 bad | top=0 | T_B = V-T_A | K = V-K | translation T_B = T_A+m"
      " | A5 sign bad")
print("  " + "-" * 100)

TT = dict(g=0, a4=0, z=0, refl=0, ksym=0, tr=0, a5=0, zr=0, rz=0)
for t, r, MAX in CONF:
    N = t + 2 * r
    g2 = a4 = nz = nrefl = nksym = ntr = a5 = zr = rz = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            an = analyse(beta, t, r)
            if an is None:
                continue
            G, V, Cv = an
            if len(G) != 2:
                continue
            g2 += 1
            (TA, wA), (TB, wB) = [G[g] for g in G]
            K = set(TA) & set(TB)
            DA = sorted(set(TA) - K)
            DB = sorted(set(TB) - K)
            # A4: the two differing pairs must both sum to V
            if len(DA) != 2 or len(DB) != 2 or \
                    sorted(x + y for x, y in zip(DA, sorted(DB, reverse=True))) != [V, V]:
                a4 += 1
            top = {}
            for TT_, ww in ((TA, wA), (TB, wB)):
                for kk, vv in P_of(list(TT_), r).items():
                    top[kk] = top.get(kk, 0) + ww * vv
            z = not any(v != 0 for v in top.values())
            refl = tuple(sorted((V - x for x in TA), reverse=True)) == tuple(sorted(TB, reverse=True))
            ksym = (set(V - x for x in K) == K)
            tr = any(tuple(sorted((x + m for x in TA), reverse=True)) == tuple(sorted(TB, reverse=True))
                     for m in set(b - a for a in TA for b in TB) if m != 0)
            nz += 1 if z else 0
            nrefl += 1 if refl else 0
            nksym += 1 if ksym else 0
            ntr += 1 if tr else 0
            if z != refl:
                zr += 1
            if refl and wB != -wA:
                a5 += 1
    print("  %4d %3d %6d | %6d | %5d | %11d | %7d | %23d | %11d"
          % (t, r, g2, a4, nz, nrefl, nksym, ntr, a5))
    for a, b in (('g', g2), ('a4', a4), ('z', nz), ('refl', nrefl), ('ksym', nksym),
                 ('tr', ntr), ('a5', a5), ('zr', zr)):
        TT[a] += b

print("")
print("  totals over %d shapes with |G| = 2:" % TT['g'])
print("     A4 (a1+b1 = a2+b2 = V) violations: %d  -- must be 0, it is forced by the tie." % TT['a4'])
print("     [Phi]_top = 0 on %d, T_B = V - T_A on %d, K = V - K on %d, a translation on %d."
      % (TT['z'], TT['refl'], TT['ksym'], TT['tr']))
print("     shapes where [Phi]_top = 0 and T_B = V - T_A DISAGREE: %d" % TT['zr'])
print("     A5: reflection holds but the signs are NOT opposite: %d" % TT['a5'])
print("")
if TT['a4'] == 0 and TT['zr'] == 0 and TT['a5'] == 0:
    print("  READING.  On every shape with two maximisers, [Phi]_top = 0 is EXACTLY T_B = V - T_A,")
    print("  which by A4 is exactly K = V - K, and the sign is then automatically opposite.")
    print("  With lemma_V_eq_C (reflection => V = C) the chain closes to condition (ii) as soon as")
    print("  A1 above holds, plus ruling out the translation branch.")
else:
    print("  SOMETHING DISAGREES -- read the columns.")
print("")
print("DONE")
