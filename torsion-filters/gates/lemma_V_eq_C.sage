# -*- coding: utf-8 -*-
# LEMA: si la reflexion intercambia los dos maximizadores, el valor de empate es EXACTAMENTE C.
# Con eso, Phi_t = 0 => (ii) queda probado salvo un unico eslabon medido.
#
# proof_G_le_2.sage probo |G| <= 2 y que, cuando |G| = 2, el empate esta entre las clases k y k+t/2 y
# el valor de empate cumple V = 2k (mod t).  V7 midio que **cuando [Phi]_top se anula** ese V cumple
# ademas V = C (mod t), 209 de 209.  Esto lo demuestra, y en la forma fuerte V = C sin el modulo.
#
# ================================================================================================
# LEMA.  Sea |G| = 2 con maximizadores T_A, T_B, y supongase T_B = C - T_A.  Entonces V = C.
#
# Prueba.  Los dos difieren en un solo intercambio de incremento entre las clases k y k' = k + t/2:
#   g_A contiene c_{k,j+1} y c_{k',j'},      luego  T_A contiene c_{k,j}   y c_{k',j'+1};
#   g_B contiene c_{k,j}   y c_{k',j'+1},    luego  T_B contiene c_{k,j+1} y c_{k',j'},
# y en las demas clases g_A = g_B, asi que T_A y T_B comparten todo lo demas.  Sea K = T_A cap T_B.
#
#   sigma_C es una INVOLUCION y sigma_C(T_A) = T_B, luego sigma_C(T_B) = T_A y
#         sigma_C(K) = sigma_C(T_A cap T_B) = sigma_C(T_A) cap sigma_C(T_B) = T_B cap T_A = K.
#   Por tanto sigma_C manda T_A \ K = {c_{k,j}, c_{k',j'+1}} SOBRE T_B \ K = {c_{k,j+1}, c_{k',j'}}.
#
# Solo hay dos maneras:
#   (1)  C - c_{k,j} = c_{k,j+1}   y   C - c_{k',j'+1} = c_{k',j'}.
#        La primera da  V = Delta_k(j) = c_{k,j} + c_{k,j+1} = C.
#   (2)  C - c_{k,j} = c_{k',j'}   y   C - c_{k',j'+1} = c_{k,j+1}.
#        Sumando:  (c_{k,j} + c_{k,j+1}) + (c_{k',j'} + c_{k',j'+1}) = 2C,  o sea  V + V = 2C,  V = C.
# En ambos casos V = C.  QED
#
# Y EL CASO (2) NO PUEDE OCURRIR, que es lo que dice su columna con un 0.  En el caso (2)
# C = c_{k,j} + c_{k',j'} = k + k' = 2k + t/2 (mod t).  Pero acabamos de obtener V = C y el paso 5 de
# proof_G_le_2 da V = 2k (mod t), luego C = 2k (mod t).  Las dos cosas exigen t/2 = 0 (mod t), falso
# para todo t >= 2.  Asi que la reflexion nunca cruza las dos clases empatadas: intercambia los dos
# elementos DENTRO de cada una.  La columna case(2) es un control forzado a 0, no una rama sin usar.
#
# COROLARIO (condicion (ii)).  V = 2k (mod t) por el paso 5 de proof_G_le_2, asi que V = C da
# 2k = C (mod t): las dos clases empatadas k y k+t/2 SON las dos clases fijas de sigma_C, y son de
# exceso porque tienen incrementos.  Eso es exactamente la condicion (ii).
#
# ESTADO DE LA CADENA, dicho sin adorno:
#   Phi_t = 0  =>  [Phi]_top = 0                                  trivial
#   [Phi]_top = 0  =>  |G| = 2 y T_B = C - T_A                    *** MEDIDO ***, 0 excepciones
#   |G| = 2 y T_B = C - T_A  =>  V = C                            probado aqui
#   V = C  =>  condicion (ii)                                     probado (paso 5 + este lema)
# El unico eslabon sin probar es el segundo, y es una afirmacion sobre productos de alternantes.
# ================================================================================================
#
# VERIFICACIONES, cada una capaz de fallar:
#   L1  sigma_C(K) = K en cada forma con |G| = 2 cerrada bajo la reflexion.
#   L2  sigma_C manda T_A \ K sobre T_B \ K, y el caso que ocurre es (1) o (2).  El caso (2) esta
#       PROHIBIDO por el argumento de residuos de arriba, asi que su columna debe salir 0: es un
#       control forzado, y un solo caso (2) mataria el argumento.
#   L3  V = C EXACTAMENTE, no solo modulo t, en toda forma cerrada.
#   L4  el senuelo: entre las formas con |G| = 2 NO cerradas, V = C debe fallar a menudo -- si V = C
#       saliera siempre, el lema no estaria usando la reflexion para nada.
#   L5  no vacuidad: debe haber formas con |G| = 2 cerradas y formas con |G| = 2 no cerradas.
#   L6  control forzado: en toda forma cerrada las dos clases empatadas deben ser las dos clases
#       fijas de sigma_C, y ambas de exceso -- o sea la condicion (ii) debe valer.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

# F1, de la auditoria socratica: este fichero EXCLUIA t = 2, que es el unico caso del que va la
# Conjetura 9.4.  Los pasos se comprobaban solo en un regimen que no es el objetivo.  Anadido.
CONF = [(2, 1, 34), (2, 2, 26), (2, 3, 20), (4, 1, 30), (4, 2, 24), (4, 3, 18),
        (6, 2, 18), (6, 3, 14), (8, 2, 16), (10, 2, 18)]


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


def setup(beta, t, r):
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    C = dict((k, sorted((beta[i] for i in cl[k]), reverse=True)) for k in E)
    return cl, E, C


def maximisers(beta, t, r, cl, E):
    N = len(beta)
    Es = set(E)
    best = None
    out = {}
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        d = sum(T[:r]) - sum(T[r:])
        if best is None or d > best:
            best, out = d, {}
        if d == best:
            out[frozenset(beta[i] for i in P if beta[i] % t in Es)] = T
    return best, out


def tie_value(C, E, r):
    inc = []
    for k in E:
        ck = C[k]
        for j in range(1, len(ck)):
            inc.append((ck[j - 1] + ck[j], k))
    inc.sort(key=lambda z: -z[0])
    V = inc[r - 1][0]
    return V, sorted(set(k for v, k in inc if v == V))


print("=" * 104)
print("L1-L6  the lemma, its case split, its decoy, and the forced consequence")
print("=" * 104)
print("")
print("     t   r  |G|=2 | closed | L1 bad  L2 bad  case(1)  case(2) | L3 bad | L6 bad || NOT closed:"
      "  V=C   V!=C")
print("  " + "-" * 100)

TOT = dict(g2=0, cl=0, l1=0, l2=0, c1=0, c2=0, l3=0, l6=0, nc=0, ncy=0, ncn=0)
for t, r, MAX in CONF:
    N = t + 2 * r
    g2 = cls = l1 = l2 = c1 = c2 = l3 = l6 = nc = ncy = ncn = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            st = setup(beta, t, r)
            if st is None:
                continue
            cl, E, Cc = st
            _, G = maximisers(beta, t, r, cl, E)
            if len(G) != 2:
                continue
            g2 += 1
            S = sorted((b for k in E for b in Cc[k]), reverse=True)
            Cv = S[0] + S[-1]
            V, tied = tie_value(Cc, E, r)
            TA, TB = [G[g] for g in G]
            refl = tuple(sorted((Cv - x for x in TA), reverse=True))
            closed = (refl == tuple(sorted(TB, reverse=True)))
            if not closed:
                nc += 1
                if V == Cv:
                    ncy += 1
                else:
                    ncn += 1
                continue
            cls += 1
            # L1
            K = set(TA) & set(TB)
            if set(Cv - x for x in K) != K:
                l1 += 1
            # L2
            DA = sorted(set(TA) - K)
            DB = sorted(set(TB) - K)
            if len(DA) != 2 or len(DB) != 2 or sorted(Cv - x for x in DA) != DB:
                l2 += 1
            else:
                # which case: (1) the reflection stays inside each class, (2) it swaps them
                if (Cv - DA[0]) % t == DA[0] % t:
                    c1 += 1
                else:
                    c2 += 1
            # L3
            if V != Cv:
                l3 += 1
            # L6
            fixed = [k for k in range(t) if (2 * k - Cv) % t == 0]
            if not (len(fixed) == 2 and all(k in E for k in fixed) and set(fixed) == set(tied)):
                l6 += 1
    print("  %4d %3d %6d | %6d | %6d %7d %8d %8d | %6d | %6d || %11d %5d %6d"
          % (t, r, g2, cls, l1, l2, c1, c2, l3, l6, nc, ncy, ncn))
    for a, b in (('g2', g2), ('cl', cls), ('l1', l1), ('l2', l2), ('c1', c1), ('c2', c2),
                 ('l3', l3), ('l6', l6), ('nc', nc), ('ncy', ncy), ('ncn', ncn)):
        TOT[a] += b

print("")
print("  totals: |G|=2 on %d shapes, of which %d are closed under the reflection." % (TOT['g2'], TOT['cl']))
print("  L1 bad %d, L2 bad %d, L3 bad (V != C exactly) %d, L6 bad %d -- all must be 0."
      % (TOT['l1'], TOT['l2'], TOT['l3'], TOT['l6']))
print("  the case split: case (1) %d, case (2) %d." % (TOT['c1'], TOT['c2']))
print("  case (2) is FORBIDDEN by residues -- it would need C = 2k and C = 2k + t/2 at once -- so")
print("  its column is a forced 0, and a single hit would kill the argument.")
print("  L4 decoy: among the %d NOT-closed shapes with |G|=2, V = C on %d and V != C on %d."
      % (TOT['nc'], TOT['ncy'], TOT['ncn']))
print("  L5 non-vacuity: closed %s, not closed %s."
      % ("yes" if TOT['cl'] else "NO", "yes" if TOT['nc'] else "NO"))
print("")
ok = (TOT['l1'] == TOT['l2'] == TOT['l3'] == TOT['l6'] == 0 and TOT['c2'] == 0
      and TOT['c1'] > 0 and TOT['cl'] > 0 and TOT['ncn'] > 0)
if ok:
    print("  PROVED AND NON-VACUOUS: the reflection fixes the common part, the two differing")
    print("  elements map onto each other, V equals C exactly, and the tied classes are the two")
    print("  fixed classes of sigma_C -- which is condition (ii).  And V = C is NOT automatic:")
    print("  it fails on %d of the not-closed shapes, so the reflection is carrying the argument."
          % TOT['ncn'])
else:
    print("  SOMETHING FAILED -- read the columns.")
print("")
print("DONE")
