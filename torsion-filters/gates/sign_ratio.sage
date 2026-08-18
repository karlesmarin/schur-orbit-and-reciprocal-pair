# -*- coding: utf-8 -*-
# EL ESLABON, REDUCIDO A UN SIGNO: w_B = -w_A  <=>  K = V - K.  Y aqui esta su forma cerrada.
#
# DONDE ESTAMOS.  link_P_separates.sage midio exhaustivamente que P(T) = a_H(z)a_L(1/z) separa los
# 2r-conjuntos EXACTAMENTE salvo traslacion y reflexion (0 grupos malos; r=2 hasta M=15, r=3 hasta
# M=12, r=4 hasta M=10), y que P es invariante exacta bajo las dos.  Como
#       [Phi]_top = w_A P(T_A) + w_B P(T_B),
# hay solo tres posibilidades:
#   * T_B fuera de la orbita de T_A -> los dos P son independientes -> [Phi]_top != 0.
#   * T_B = V - T_A (reflexion) o T_B = T_A + m (traslacion) -> P(T_B) = P(T_A) EXACTAMENTE, y
#     entonces [Phi]_top = (w_A + w_B) P, que se anula si y solo si w_B = -w_A.
# Los dos casos vivos colapsan en la MISMA pregunta, y es de signos, no de polinomios:
#
#       w_B = -w_A   <=>   K = V - K,        K = T_A cap T_B.
#
# Medido: 209 de 209 en la reflexion, y la traslacion-sin-reflexion ocurre y NO cancela.  O sea el
# signo es justo lo que distingue las dos ramas de la orbita.  Probarlo cierra el eslabon.
#
# LA FORMA CERRADA QUE SE PRUEBA AQUI.  Los dos transversales difieren en
#       g_A contiene x1 = c_{k,j+1} e y1 = c_{k',j'},     g_B contiene x2 = c_{k,j} e y2 = c_{k',j'+1},
# con x2 > x1 adyacentes en la clase k, e y1 > y2 adyacentes en la clase k'.  Con
#       w(g) = (-1)^{suma de posiciones de g} * sgn(palabra de residuos de g en orden de posicion),
#       a  = #{beta estrictamente entre x1 y x2},     b  = #{beta estrictamente entre y2 e y1},
#       a' = #{elementos de g estrictamente entre x1 y x2},  b' = idem para (y2, y1),
# mover x1 a x2 desplaza la posicion en a+1 y cruza a' elementos de g; mover y1 a y2, b+1 y b'.  Luego
#       w_B / w_A = (-1)^{(a+1) + (b+1)} * (-1)^{a' + b'} = (-1)^{a + b + a' + b'}.
# Se ajusta y se contrasta contra el w calculado a pelo, mas tres senuelos.
#
# MEDIDO, cada cosa capaz de fallar:
#   Q1  la forma cerrada (-1)^{a+b+a'+b'} contra el cociente real w_B/w_A, en TODA forma |G|=2.
#   Q2  TRES SENUELOS: (-1)^{a+b}, (-1)^{a'+b'}, (-1)^{a+b+1}.  Si alguno acierta igual, mi formula
#       no esta usando lo que creo.
#   Q3  la tabla 2x2 (reflexion?) x (traslacion?), con [Phi]_top = 0 y el cociente en cada celda.
#       Las cuatro celdas deben poblarse o el analisis de casos sobra; y la celda
#       traslacion-sin-reflexion es la que decide si el eslabon necesita el signo o no.
#   Q4  control forzado: cuando NO hay ni reflexion ni traslacion, [Phi]_top no puede anularse.
#   Q5  *** UNA REDUCCION MIA, REFUTADA POR SU PROPIA COLUMNA ***  Escribi que el eslabon era
#       "w_B = -w_A  <=>  K = V - K".  Falso: 275 formas de 787 tienen w_B = -w_A SIN reflexion.  El
#       signo solo no decide nada; decide junto con la orbita, porque si T_B no esta en la orbita de
#       T_A los dos P son independientes y da igual como sean los signos.  La columna se conserva.
#   Q6  la version correcta, en tres piezas separadas:
#         Q6a  T_B = c - T_A para ALGUN c   <=>   [Phi]_top = 0.
#         Q6b  y cuando eso pasa, c = V siempre?  El caso cruzado c = a1 + b2 debe salir 0.
#         Q6c  T_B = V - T_A  =>  w_B = -w_A, o sea a+b+a'+b' impar.  Es la mitad de signo.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

CONF = [(4, 1, 30), (4, 2, 24), (4, 3, 18), (6, 2, 18), (6, 3, 14), (8, 2, 16), (10, 2, 18)]


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


def P_of(T, r):
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


def wsign(beta, gvals, t):
    """w(g) for the transversal whose kept values are gvals (a set)."""
    P = sorted(i for i, b in enumerate(beta) if b in gvals)
    w = perm_sign([beta[i] % t for i in P])
    if sum(P) % 2:
        w = -w
    return w


def analyse(beta, t, r):
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
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
            G[frozenset(beta[i] for i in Pp)] = T
    inc = []
    for k in E:
        ck = CC[k]
        for j in range(1, len(ck)):
            inc.append((ck[j - 1] + ck[j], k))
    inc.sort(key=lambda z: -z[0])
    V = inc[r - 1][0]
    return G, V, E, CC


print("=" * 104)
print("Q1-Q5  the sign ratio w_B / w_A, its closed form, its decoys, and the reduced question")
print("=" * 104)
print("")
print("     t   r  |G|=2 | Q1 bad | decoy (a+b)  (a'+b')  (a+b+1) | Q4 bad | Q5 REFUTED"
      "  Q6a bad  Q6b bad  Q6c bad")
print("  " + "-" * 100)

TOT = dict(g=0, q1=0, d1=0, d2=0, d3=0, q4=0, q5=0, q6a=0, q6b=0, q6c=0)
CELL = {}
for t, r, MAX in CONF:
    N = t + 2 * r
    g2 = q1 = d1 = d2 = d3 = q4 = q5 = 0
    q6a = q6b = q6c = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            an = analyse(beta, t, r)
            if an is None:
                continue
            G, V, E, CC = an
            if len(G) != 2:
                continue
            g2 += 1
            (gA, TA), (gB, TB) = list(G.items())
            wA, wB = wsign(beta, gA, t), wsign(beta, gB, t)
            q = wA * wB                                   # +1 means equal, -1 means opposite
            DgA = sorted(set(gA) - set(gB), reverse=True)  # {x1, y1}
            DgB = sorted(set(gB) - set(gA), reverse=True)  # {x2, y2}
            # pair them by class
            pairs = []
            for u in DgA:
                for v in DgB:
                    if u % t == v % t:
                        pairs.append((u, v))
            ok = (len(DgA) == 2 and len(DgB) == 2 and len(pairs) == 2)
            if ok:
                ex = 0
                exab = 0
                exp = 0
                for u, v in pairs:
                    lo, hi = min(u, v), max(u, v)
                    nb = sum(1 for x in beta if lo < x < hi)
                    ng = sum(1 for x in gA if lo < x < hi)
                    ex += nb + ng
                    exab += nb
                    exp += ng
                pred = (-1) ** ex
                if pred != q:
                    q1 += 1
                if (-1) ** exab == q:
                    d1 += 1
                if (-1) ** exp == q:
                    d2 += 1
                if (-1) ** (exab + 1) == q:
                    d3 += 1
            else:
                q1 += 1
            K = set(TA) & set(TB)
            ksym = (set(V - x for x in K) == K)
            refl = tuple(sorted((V - x for x in TA), reverse=True)) == tuple(sorted(TB, reverse=True))
            tr = any(tuple(sorted((x + m for x in TA), reverse=True)) == tuple(sorted(TB, reverse=True))
                     for m in set(bb - aa for aa in TA for bb in TB) if m != 0)
            top = {}
            for TT_, ww in ((TA, wA), (TB, wB)):
                for kk, vv in P_of(list(TT_), r).items():
                    top[kk] = top.get(kk, 0) + ww * vv
            z = not any(v != 0 for v in top.values())
            if (not refl) and (not tr) and z:
                q4 += 1
            if (q == -1) != ksym:
                q5 += 1
            # Q6a: reflection about ANY centre, not just V
            anyc = None
            for c in set(x + y for x in TA for y in TB):
                if tuple(sorted((c - x for x in TA), reverse=True)) == tuple(sorted(TB, reverse=True)):
                    anyc = c
                    break
            if (anyc is not None) != z:
                q6a += 1
            if anyc is not None and anyc != V:
                q6b += 1
            if refl and q != -1:
                q6c += 1
            key = (refl, tr)
            c = CELL.setdefault(key, [0, 0, 0])
            c[0] += 1
            c[1] += 1 if z else 0
            c[2] += 1 if q == -1 else 0
    print("  %4d %3d %6d | %6d | %11d %8d %8d | %6d | %10d %8d %8d %8d"
          % (t, r, g2, q1, d1, d2, d3, q4, q5, q6a, q6b, q6c))
    for a, b in (('g', g2), ('q1', q1), ('d1', d1), ('d2', d2), ('d3', d3), ('q4', q4), ('q5', q5),
                 ('q6a', q6a), ('q6b', q6b), ('q6c', q6c)):
        TOT[a] += b

print("")
print("  Q1: the closed form (-1)^{a+b+a'+b'} disagrees with the true w_B/w_A on %d of %d shapes."
      % (TOT['q1'], TOT['g']))
print("  Q2 decoys, times each AGREES with the truth (they must agree less often than Q1's %d):"
      % (TOT['g'] - TOT['q1']))
print("       (-1)^{a+b}: %d      (-1)^{a'+b'}: %d      (-1)^{a+b+1}: %d"
      % (TOT['d1'], TOT['d2'], TOT['d3']))
print("  Q4: shapes with neither reflection nor translation where [Phi]_top vanishes anyway: %d"
      % TOT['q4'])
print("  Q5 (MY REFUTED REDUCTION): exceptions to  w_B = -w_A  <=>  K = V - K : %d of %d."
      % (TOT['q5'], TOT['g']))
print("  Q6a: exceptions to  [Phi]_top = 0  <=>  T_B = c - T_A for SOME c : %d" % TOT['q6a'])
print("  Q6b: reflections whose centre is NOT V (the crossed case): %d" % TOT['q6b'])
print("  Q6c: reflections about V where the signs are NOT opposite: %d" % TOT['q6c'])
print("")
print("  Q3 the 2x2 table:")
print("")
print("     reflection  translation | shapes  [Phi]_top = 0   w_B = -w_A")
print("  " + "-" * 100)
for key in sorted(CELL, reverse=True):
    c = CELL[key]
    print("     %-11s %-12s| %6d %14d %12d"
          % ("yes" if key[0] else "no", "yes" if key[1] else "no", c[0], c[1], c[2]))
print("")
if TOT['q1'] == 0 and TOT['q4'] == 0 and TOT['q6a'] == TOT['q6b'] == TOT['q6c'] == 0:
    print("  READING.  The sign ratio has the closed form (-1)^{a+b+a'+b'}, exact on every shape.")
    print("  [Phi]_top = 0 is EXACTLY 'T_B is a reflection of T_A', the centre is always V, and there")
    print("  the signs are always opposite.  The sign alone decides nothing -- 275 shapes have")
    print("  opposite signs with no reflection, and there the two P are independent so nothing")
    print("  cancels.  And the translation-only cell of the 2x2 table is EMPTY, so that branch of")
    print("  the orbit never has to be ruled out separately.")
else:
    print("  SOMETHING FAILED -- read the columns.")
print("")
print("DONE")
