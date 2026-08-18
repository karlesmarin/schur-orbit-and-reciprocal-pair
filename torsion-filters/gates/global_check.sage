# -*- coding: utf-8 -*-
# DOS COMPROBACIONES GLOBALES que la mirada de conjunto exige antes de opinar de nada.
#
# X1  NO VACUIDAD EN EL RANGO QUE FALTABA.  El paper prueba el reciproco en ell(lambda) <= N/2
#     (Thm stable) y dice que lo que queda abierto es justamente ell(lambda) > N/2.  Nuestra prueba
#     NO menciona ell(lambda) en ningun sitio.  Eso es o una ruta distinta, o una vacuidad: si las
#     216 formas verificadas fueran todas rectangulos de altura N/2, no habriamos probado nada nuevo.
#     Hay que contar cuantas caen ESTRICTAMENTE por encima de N/2 y ensenar una.
#
# X2  EL SESGO DEL MUESTREO.  Todo lo de hoy se barrio acotando |lambda|, que muestrea betas
#     concentrados.  La auditoria socratica ya seniala que acotar la PARTE MAYOR (o sea, tomar todos
#     los N-subconjuntos de {0..M}) alcanza betas mucho mas dispersos, y con esa rebanada rehizo
#     |G| <= 2.  Los resultados NUEVOS -- el estrato de abajo, H6 y la cadena de t = 2 -- solo se han
#     visto con el muestreo por |lambda|.  Aqui se rehace la conclusion de t = 2 con el OTRO
#     muestreo, que es una rebanada independiente y no un rango mas grande del mismo.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools


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


def topdeg_dict(T, r):
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


def from_beta(beta, t, r):
    """beta strictly decreasing; returns criterion, [Phi]_top and branch (b)."""
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
    best = None
    G = {}
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        d = sum(T[:r]) - sum(T[r:])
        if best is None or d > best:
            best, G = d, {}
        if d == best:
            w = perm_sign([beta[i] % t for i in P])
            if sum(P) % 2:
                w = -w
            G[frozenset(beta[i] for i in P if beta[i] % t in Es)] = (T, w)
    top = {}
    for g, (T, w) in G.items():
        for k, v in topdeg_dict(list(T), r).items():
            top[k] = top.get(k, 0) + w * v
    lam = [beta[i] - (N - 1 - i) for i in range(N)]
    w_ = lam[0] + lam[N - 1]
    brB = all(lam[i] + lam[N - 1 - i] == w_ for i in range(N)) and (w_ % 2 == 1)
    ell = sum(1 for x in lam if x > 0)
    return dict(crit=conc and cond_ii, top=dict((k, v) for k, v in top.items() if v != 0),
                brB=brB, ell=ell, lam=lam, N=N)


print("=" * 106)
print("X1  the 216 verified vanishing shapes at t = 2: how many are ABOVE Littlewood's range?")
print("=" * 106)
print("")
print("     t   r    N | vanishing | ell < N/2  ell = N/2  ell > N/2 | of the ell > N/2, rectangles?")
print("  " + "-" * 102)

EX = []
TOT = [0, 0, 0, 0]
for t, r, MAX in ((2, 1, 40), (2, 2, 32), (2, 3, 24), (2, 4, 18), (2, 5, 14)):
    N = t + 2 * r
    nv = lo = eq = hi = rect = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            lam = list(l) + [0] * (N - len(l))
            beta = [lam[i] + N - 1 - i for i in range(N)]
            a = from_beta(beta, t, r)
            if a is None or not a['crit']:
                continue
            nv += 1
            if 2 * a['ell'] < N:
                lo += 1
            elif 2 * a['ell'] == N:
                eq += 1
            else:
                hi += 1
                k = lam[0]
                if all(x == k for x in lam[:N // 2]) and all(x == 0 for x in lam[N // 2:]):
                    rect += 1
                if len(EX) < 6:
                    EX.append((t, r, list(l), a['ell'], N))
    TOT[0] += nv
    TOT[1] += lo
    TOT[2] += eq
    TOT[3] += hi
    print("  %4d %3d %4d | %9d | %9d %10d %10d | %28d" % (t, r, N, nv, lo, eq, hi, rect))

print("")
print("  totals: %d vanishing shapes; %d below N/2, %d at N/2, %d STRICTLY ABOVE N/2."
      % (TOT[0], TOT[1], TOT[2], TOT[3]))
print("  Theorem (stable) of the paper already covers ell <= N/2.  The %d above it are the ones the"
      % TOT[3])
print("  paper leaves open, and they are what the new argument has to be earning.")
print("")
print("  examples above the range:")
for t, r, lam, ell, N in EX:
    print("     t=%d r=%d N=%d  ell=%d > %d   lambda = %s" % (t, r, N, ell, N // 2, str(lam)))

# ---------------------------------------------------------------- X2 -----------------------------
print("")
print("=" * 106)
print("X2  the t = 2 conclusion re-run on the OTHER slice: all beta = N-subsets of {0..M}")
print("=" * 106)
print("")
print("     t   r    N    M   shapes | [Phi]top=0  criterion  branch (b) | disagreements"
      "  | max |lambda| reached")
print("  " + "-" * 102)

bad2 = 0
for t, r, M in ((2, 1, 26), (2, 2, 20), (2, 3, 16), (2, 4, 13), (2, 5, 12)):
    N = t + 2 * r
    nsh = nz = ncr = nb = dis = 0
    mxs = 0
    for T in itertools.combinations(range(M + 1), N):
        beta = sorted(T, reverse=True)
        a = from_beta(beta, t, r)
        if a is None:
            continue
        nsh += 1
        z = not a['top']
        nz += 1 if z else 0
        ncr += 1 if a['crit'] else 0
        nb += 1 if a['brB'] else 0
        if not (z == a['crit'] == a['brB']):
            dis += 1
        mxs = max(mxs, sum(a['lam']))
    bad2 += dis
    print("  %4d %3d %4d %4d %8d | %10d %10d %11d | %14d | %20d"
          % (t, r, N, M, nsh, nz, ncr, nb, dis, mxs))

print("")
print("  X2 disagreements between [Phi]_top = 0, the criterion, and branch (b): %d" % bad2)
print("  This slice bounds the LARGEST PART, not |lambda|, so it reaches beta sets that are far more")
print("  spread out than the |lambda|-bounded sweep; it is an independent slice, not a bigger range.")
print("")
if TOT[3] > 0 and bad2 == 0:
    print("  BOTH CHECKS PASS: the new argument is NOT vacuous -- it covers %d shapes strictly above"
          % TOT[3])
    print("  Littlewood's range, exactly where the paper says the converse was missing -- and the")
    print("  conclusion survives a completely different sampling of beta.")
else:
    print("  SOMETHING FAILED -- read the columns.")
print("")
print("DONE")
