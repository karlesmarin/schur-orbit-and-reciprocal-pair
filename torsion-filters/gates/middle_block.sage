# -*- coding: utf-8 -*-
# THE MIDDLE BLOCK: a PROVED sufficient condition for non-vanishing, and where it stops.
#
# WHAT topdeg_gate.sage FOUND.  The degree-Dmax part of Phi_t survives in 12404 of the 12424
# criterion-failing shapes, and all 20 exceptions have condition (ii) HOLDING and (i) failing.  So
# the top-degree obstruction is not a heuristic: it nearly IS the converse.  This gate extracts the
# part of it that is a theorem, and measures exactly how much of the converse that part buys.
#
# ------------------------------------------------------------------------------------------------
# THEOREM (new, and the proof is three lines).  Let S = s_1 > .. > s_n be the excess values,
# n = 2r + e, e = |E| the number of excess classes, and let
#       M  =  { s_{r+1}, .., s_{r+e} }                      -- the MIDDLE BLOCK.
# If M contains exactly one element of each excess class, then  Phi_t(lambda; z) is NOT identically
# zero.
#
# Proof.  Each Laplace term is indexed by the transversal g of kept rows, one per excess class, and
# carries T = S \ g.  Writing H, L for the top and bottom halves of T,
#       deg(T) = sum H - sum L  <=  (s_1 + .. + s_r) - (s_{n-r+1} + .. + s_n),
# because the top r elements of T are dominated termwise by s_1..s_r and its bottom r dominate
# s_{n-r+1}..s_n.  Equality forces T to contain both extreme blocks, i.e. g is contained in M, and
# |g| = |M| = e forces g = M.  So if M is a transversal the bound is attained by EXACTLY ONE term,
# and the degree-Dmax part of Phi_t is that single term's.  A polynomial with a nonzero homogeneous
# part is nonzero.  QED
#
# *** SOBRE-CORRECCION MIA, y la debo a la SEGUNDA auditoria socratica.  Aqui puse antes que decir
# "la parte de grado Dmax ES +- a_H(z)a_L(1/z)" era FALSO por los 2^r sectores de orientacion.  Eso
# solo es cierto bajo la graduacion por sum_j |e_j|; bajo la graduacion por el GRADO TOTAL sum_j e_j
# -- que es la que usa todo el pipeline, porque topdeg_dict() calcula el sector todo-positivo -- la
# frase original era EXACTA.  O sea retire una frase verdadera.  Comprobado por el auditor sobre 89
# formas: la componente de grado total maximo del determinante N x N completo es
# 2 * sum_{g en G} w_g a_{H_g}(z) a_{L_g}(1/z), 0 discrepancias.
# CONVENIO, y de aqui en adelante uno solo: [Phi]_top es la componente de GRADO TOTAL maximo, o sea
# el sector todo-positivo.  Bajo sum|e_j| hay 2^r sectores, son disjuntos porque ningun e_j es 0 en
# grado maximo, y z_j -> 1/z_j los permuta con signo uniforme, asi que se anulan todos a la vez: por
# eso da igual cual se mire, pero hay que decir cual se mira. ***
#
# The same argument proves the strictly more general
#       COROLLARY.  If the maximiser of deg over transversals is unique, Phi_t is not zero,
# and the middle-block case is the one where the maximiser is forced by inspection.
# ------------------------------------------------------------------------------------------------
#
# MEASURED HERE:
#   M1  the coverage: of the shapes where the criterion fails, how many does the theorem settle,
#       how many does the corollary settle, and how many need the full top-degree part.
#   M2  the stratification of the failures: (i) alone, (ii) alone, both.  The 20 top-degree
#       exceptions all had (ii) HOLDING, so the sharp question is whether
#             (ii) fails  =>  [Phi_t]_top != 0
#       has any exception at all.  If it has none, then Phi_t = 0 => (ii) is the provable half.
#   M3  for the exceptions, how far below Dmax the first surviving degree sits.
#
# CONTROLS, each able to fail:
#   C1  forced: if the criterion HOLDS the middle block can never be a transversal, the maximiser
#       can never be unique, and [Phi]_top must vanish.  Any hit kills the theorem.
#   C2  acceptance: on the shapes the theorem settles, the top-degree part must have exactly one
#       contributing term and its coefficients must all be +-1.
#   C3  non-vacuity: the theorem must settle a nonzero number of shapes, and must FAIL to settle a
#       nonzero number, or it is either empty or trivially everything.
#   C4  the stratum "(ii) fails" must be non-empty in every configuration, or M2 says nothing.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

# F1, de la auditoria socratica: este fichero EXCLUIA t = 2, y el corolario del maximizador unico
# que el paso 2 de conj_crit_t2 usa se apoyaba en el.  Anadido.
CONF = [(2, 1, 34), (2, 2, 26), (2, 3, 20), (4, 1, 30), (4, 2, 24), (4, 3, 18),
        (6, 2, 18), (6, 3, 14), (8, 2, 16)]


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


def terms(beta, t, r):
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    out = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        Ps = set(P)
        out.append((w, [beta[i] for i in range(N) if i not in Ps]))
    return out


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
    return {k: v for k, v in D.items() if v != 0}


def full_expansion(tm, r):
    n = 2 * r
    D = {}
    for w, T in tm:
        for q in itertools.permutations(range(n)):
            e = [0] * r
            for a in range(n):
                e[q[a] // 2] += (T[a] if q[a] % 2 == 0 else -T[a])
            k = tuple(e)
            D[k] = D.get(k, 0) + w * perm_sign(list(q))
    return {k: v for k, v in D.items() if v != 0}


def analyse(beta, t, r):
    cl = {}
    for b in beta:
        cl.setdefault(b % t, []).append(b)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    S = sorted((b for k in E for b in cl[k]), reverse=True)
    n = len(S)
    e = len(E)
    C = S[0] + S[-1]
    conc = set(C - b for b in S) == set(S)
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    cond_ii = (len(fixed) == 2 and all(k in E for k in fixed))
    crit = conc and cond_ii
    mid = S[r:r + e]
    mid_transversal = (len(set(b % t for b in mid)) == e)
    tm = terms(beta, t, r)
    degs = [sum(T[:r]) - sum(T[r:]) for _, T in tm]
    Dmax = max(degs)
    G = [i for i, d in enumerate(degs) if d == Dmax]
    top = {}
    for i in G:
        w, T = tm[i]
        for k, v in topdeg_dict(T, r).items():
            top[k] = top.get(k, 0) + w * v
    top = {k: v for k, v in top.items() if v != 0}
    bound = sum(S[:r]) - sum(S[n - r:])
    return dict(crit=crit, conc=conc, cond_ii=cond_ii, e=e, S=S, tm=tm,
                mid_transversal=mid_transversal, Dmax=Dmax, bound=bound,
                nG=len(G), top=top, G=G)


# ---------------------------------------------------------------- C1, C2, C3 ---------------------
print("=" * 104)
print("C1  forced consequences on the shapes where the criterion HOLDS      C2  acceptance")
print("=" * 104)
print("")
print("     t   r  crit shapes  mid transv.  unique max.  [Phi]top!=0 | C2: Dmax!=bound  #max!=1"
      "  top=0  coef!=+-1")
print("  " + "-" * 100)

c1fail = c2fail = 0
for t, r, MAX in CONF:
    N = t + 2 * r
    ncrit = m1 = m2 = m3 = 0
    d1 = d2 = d3 = d4 = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            if a['crit']:
                ncrit += 1
                m1 += 1 if a['mid_transversal'] else 0
                m2 += 1 if a['nG'] == 1 else 0
                m3 += 1 if a['top'] else 0
            elif a['mid_transversal']:
                d1 += 1 if a['Dmax'] != a['bound'] else 0
                d2 += 1 if a['nG'] != 1 else 0
                d3 += 0 if a['top'] else 1
                d4 += 1 if any(abs(v) != 1 for v in a['top'].values()) else 0
    c1fail += m1 + m2 + m3
    c2fail += d1 + d2 + d3
    print("  %4d %3d %12d %12d %12d %12d | %14d %8d %6d %10d"
          % (t, r, ncrit, m1, m2, m3, d1, d2, d3, d4))

print("")
print("  the last column is a PREDICTION OF MINE, and it is refuted: I wrote that the surviving")
print("  coefficients would be +-1.  They need not be -- a_H(z) a_L(1/z) is a PRODUCT of two")
print("  alternants, and two different (alpha,beta) can give the same exponent vector, e.g.")
print("  H={5,3}, L={2,0} sends both (5-2,3-0) and (3-0,5-2) to (3,3).  The theorem does not use")
print("  the value of the coefficient, only that the polynomial is nonzero, so it survives intact.")
print("")
if c1fail or c2fail:
    print("  C1/C2 FAILED -- the theorem is wrong.  Stop here.")
    raise SystemExit(1)
print("  C1 PASS: on every shape the criterion calls a zero, the middle block is never a")
print("  transversal, the degree maximiser is never unique, and the top-degree part vanishes.")
print("  C2 PASS: whenever the middle block IS a transversal, Dmax equals the a priori bound, the")
print("  maximiser is unique, and the top-degree part is nonzero -- which is the whole theorem.")

# ---------------------------------------------------------------- M1, M2 -------------------------
print("")
print("=" * 104)
print("M1/M2  the converse, stratified by WHICH condition fails and by WHAT settles it")
print("=" * 104)
print("")
print("     t   r  crit fails | (i) only  (ii) only    both | theorem  corollary  [Phi]_top | LEFT")
print("  " + "-" * 100)

LEFT = []
TOT = dict(f=0, i=0, ii=0, b=0, th=0, co=0, tp=0)
STRAT = {}
for t, r, MAX in CONF:
    N = t + 2 * r
    nf = ni = nii = nb = nth = nco = ntp = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None or a['crit']:
                continue
            nf += 1
            key = (a['conc'], a['cond_ii'])
            st = STRAT.setdefault(key, [0, 0])
            st[0] += 1
            if a['top']:
                st[1] += 1
            if a['conc'] and not a['cond_ii']:
                nii += 1
            elif (not a['conc']) and a['cond_ii']:
                ni += 1
            else:
                nb += 1
            nth += 1 if a['mid_transversal'] else 0
            nco += 1 if a['nG'] == 1 else 0
            if a['top']:
                ntp += 1
            else:
                LEFT.append((t, r, list(l), beta, a))
    print("  %4d %3d %11d | %8d %9d %7d | %7d %10d %10d | %4d"
          % (t, r, nf, ni, nii, nb, nth, nco, ntp, nf - ntp))
    TOT['f'] += nf
    TOT['i'] += ni
    TOT['ii'] += nii
    TOT['b'] += nb
    TOT['th'] += nth
    TOT['co'] += nco
    TOT['tp'] += ntp

print("")
print("  totals: criterion fails %d = (i) only %d + (ii) only %d + both %d"
      % (TOT['f'], TOT['i'], TOT['ii'], TOT['b']))
print("          settled by the THEOREM (middle block a transversal): %d  (%.1f%%)"
      % (TOT['th'], 100.0 * TOT['th'] / TOT['f']))
print("          settled by the COROLLARY (unique degree maximiser):  %d  (%.1f%%)"
      % (TOT['co'], 100.0 * TOT['co'] / TOT['f']))
print("          settled by the full top-degree part:                 %d  (%.1f%%)"
      % (TOT['tp'], 100.0 * TOT['tp'] / TOT['f']))
print("          C3 non-vacuity: theorem settles some (%s) and not all (%s)"
      % ("yes" if TOT['th'] > 0 else "NO", "yes" if TOT['th'] < TOT['f'] else "NO"))
print("")
print("  M2  the sharp question -- does the top-degree part ever miss when (ii) FAILS?")
print("")
print("     concentric (i)   condition (ii)   shapes   [Phi]_top survives   misses")
print("  " + "-" * 100)
for key in sorted(STRAT):
    n, s = STRAT[key]
    print("     %-16s %-16s %7d %20d %8d"
          % ("holds" if key[0] else "fails", "holds" if key[1] else "fails", n, s, n - s))
print("")
miss_ii = sum(STRAT[k][0] - STRAT[k][1] for k in STRAT if not k[1])
n_ii = sum(STRAT[k][0] for k in STRAT if not k[1])
print("  (ii) FAILS: %d shapes, top-degree part misses %d of them.  C4 non-vacuity: %s"
      % (n_ii, miss_ii, "PASS" if n_ii > 0 else "VACUOUS"))

# ---------------------------------------------------------------- M3 -----------------------------
print("")
print("=" * 104)
print("M3  the shapes the top-degree part misses: how far down is the first surviving degree?")
print("=" * 104)
print("")
print("     t   r  lambda                     |E|  Dmax  bound  #maximisers  first surviving deg")
print("  " + "-" * 100)
gaps = {}
for t, r, lam, beta, a in LEFT:
    FE = full_expansion(a['tm'], r)
    dd = sorted(set(sum(abs(x) for x in k) for k in FE), reverse=True)
    first = dd[0] if dd else None
    gap = (a['Dmax'] - first) if first is not None else None
    gaps[gap] = gaps.get(gap, 0) + 1
    print("  %4d %3d  %-26s %3d %5d %6d %12d %14s   (gap %s)"
          % (t, r, str(lam), a['e'], a['Dmax'], a['bound'], a['nG'],
             str(first), str(gap)))
print("")
print("  gap below Dmax of the first surviving degree: %s" % (sorted(gaps.items()),))
print("")
print("DONE")
