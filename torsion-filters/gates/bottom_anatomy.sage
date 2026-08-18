# -*- coding: utf-8 -*-
# EL ESTRATO DE ABAJO ES MAS SIMPLE QUE EL DE ARRIBA: no hace falta rigidez de Schur, solo signos.
#
# bottom_stratum.sage midio H6 SIN EXCEPCIONES sobre 10289 formas donde el criterio falla, t = 2,4,6,8:
#       criterio falla   =>   [Phi]_top != 0   O   [Phi]_bot != 0,
# o sea el reciproco entero se reduce a DOS condiciones y no a una escalera de estratos.  Y lo cubre
# tambien en t = 2, que es donde vive conj:crit.  La mitad de arriba ya es teorema.  Esta es la otra.
#
# ------------------------------------------------------------------------------------------------
# LA OBSERVACION QUE LO SIMPLIFICA TODO.  En el estrato de abajo el emparejamiento es el de
# CONSECUTIVOS, con diferencias d_i = u_{2i-1} - u_{2i}, y solo hay que repartir las r parejas entre
# las r variables.  Esa permutacion mueve BLOQUES DE TAMANO 2, y una permutacion de bloques de tamano
# par tiene signo  det(P_sigma (x) I_2) = det(P_sigma)^2 = +1.  Luego
#
#       [A(T)]_bot  =  sum_{sigma in S_r}  prod_j z_{sigma(j)}^{d_j}  =  |Stab(d)| * m_{sort(d)}(z),
#
# TODOS los signos son +.  Un solo A(T) NUNCA se cancela por dentro en el estrato de abajo -- al
# reves que arriba, donde a_H(z)a_L(1/z) es un alternante y se cancela solo en cuanto hay repeticion.
# Y las funciones monomiales m_lambda son linealmente independientes, asi que
#
#       [Phi]_bot == 0   <=>   para cada multiconjunto D:   sum_{T in Gbot, sort(d(T)) = D} w(T) = 0
#
# -- un RECUENTO CON SIGNO, sin polinomios dentro y sin necesitar el teorema de Purbhoo-van
# Willigenburg.  Es una condicion estrictamente mas elemental que la del estrato de arriba.
# ------------------------------------------------------------------------------------------------
#
# MEDIDO, cada cosa capaz de fallar:
#   C1  la observacion del signo: perm_sign de la permutacion de bloques debe ser +1 SIEMPRE, y
#       [A(T)]_bot debe tener todos los coeficientes iguales y positivos, valiendo |Stab(d)|.
#   C2  EL CRITERIO: [Phi]_bot = 0  <=>  todas las sumas con signo por multiconjunto D se anulan.
#       Las dos direcciones, contadas aparte.
#   C3  la anatomia de Gbot: cuantos T minimizan degmin, y como se reparten.  Arriba salio |G| <= 2;
#       aqui hay que MIRAR, no suponer.
#   C4  cuantos multiconjuntos D distintos aparecen dentro de un mismo Gbot -- si es siempre 1, el
#       criterio se reduce a "sum w(T) = 0" a secas.
#   C5  SENUELO: ignorar el reparto por D y pedir solo sum_{T in Gbot} w(T) = 0.  Debe SOBRE-disparar
#       si D reparte de verdad, o el reparto no esta haciendo nada.
#   C6  control forzado: si el criterio (i)+(ii) VALE, el recuento con signo debe anularse.
#   C7  no vacuidad: tiene que haber formas con [Phi]_bot = 0 y formas con [Phi]_bot != 0.
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


def botdeg_dict(T, r):
    D = {}
    n = 2 * r
    for s in itertools.permutations(range(r)):
        q = [0] * n
        e = [0] * r
        for i in range(r):
            q[2 * i] = 2 * s[i]
            q[2 * i + 1] = 2 * s[i] + 1
            e[s[i]] = T[2 * i] - T[2 * i + 1]
        D[tuple(e)] = D.get(tuple(e), 0) + perm_sign(q)
    return dict((k, v) for k, v in D.items() if v != 0)


def block_signs(r):
    out = []
    for s in itertools.permutations(range(r)):
        q = [0] * (2 * r)
        for i in range(r):
            q[2 * i] = 2 * s[i]
            q[2 * i + 1] = 2 * s[i] + 1
        out.append(perm_sign(q))
    return out


def analyse(beta, t, r):
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    S = sorted((beta[i] for k in E for i in cl[k]), reverse=True)
    C = S[0] + S[-1]
    conc = set(C - b for b in S) == set(S)
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    cond_ii = (len(fixed) == 2 and all(k in E for k in fixed))
    tm = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        Ps = set(P)
        tm.append((w, tuple(beta[i] for i in range(N) if i not in Ps)))
    db = [sum(T[2 * i] - T[2 * i + 1] for i in range(r)) for _, T in tm]
    Dmin = min(db)
    Gb = [(w, T) for (w, T), d in zip(tm, db) if d == Dmin]
    bot = {}
    for w, T in Gb:
        for k, v in botdeg_dict(list(T), r).items():
            bot[k] = bot.get(k, 0) + w * v
    bot = dict((k, v) for k, v in bot.items() if v != 0)
    byD = {}
    for w, T in Gb:
        d = tuple(sorted((T[2 * i] - T[2 * i + 1] for i in range(r)), reverse=True))
        byD[d] = byD.get(d, 0) + w
    return dict(crit=conc and cond_ii, conc=conc, cond_ii=cond_ii, Gb=Gb, bot=bot,
                byD=byD, nG=len(Gb), nD=len(byD))


print("=" * 106)
print("C1  the block-permutation sign, and [A(T)]_bot = |Stab(d)| m_sort(d)")
print("=" * 106)
print("")
print("     r   block signs all +1   coefficients all equal & positive   value = |Stab(d)|")
print("  " + "-" * 102)
c1bad = 0
for r in (1, 2, 3, 4):
    sg = block_signs(r)
    ok1 = all(x == 1 for x in sg)
    ok2 = ok3 = True
    for T in itertools.combinations(range(14), 2 * r):
        Td = sorted(T, reverse=True)
        D = botdeg_dict(Td, r)
        d = [Td[2 * i] - Td[2 * i + 1] for i in range(r)]
        st = prod([factorial(list(d).count(x)) for x in set(d)])
        if any(v <= 0 for v in D.values()) or len(set(D.values())) != 1:
            ok2 = False
        if list(D.values())[0] != st:
            ok3 = False
    if not (ok1 and ok2 and ok3):
        c1bad += 1
    print("  %4d %20s %38s %20s"
          % (r, "yes" if ok1 else "NO", "yes" if ok2 else "NO", "yes" if ok3 else "NO"))
print("")
if c1bad:
    print("  C1 FAILED -- stop.")
    raise SystemExit(1)
print("  C1 PASS: permuting blocks of size 2 has sign det(P)^2 = +1, so the bottom part of a single")
print("  A(T) is a POSITIVE sum and can never cancel by itself.  All cancellation is across T.")

print("")
print("=" * 106)
print("C2-C7  the signed-count criterion for [Phi]_bot")
print("=" * 106)
print("")
print("     t   r  shapes | max|Gb|  |Gb| distribution      max #D | C2 bad | C5 decoy bad | C6 bad")
print("  " + "-" * 102)

TOT = dict(sh=0, c2=0, c5=0, c6=0, z=0, nz=0)
for t, r, MAX in CONF:
    N = t + 2 * r
    nsh = c2 = c5 = c6 = nz = nnz = 0
    dist = {}
    mxD = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            nsh += 1
            dist[a['nG']] = dist.get(a['nG'], 0) + 1
            mxD = max(mxD, a['nD'])
            crit_bot = all(v == 0 for v in a['byD'].values())
            zero = not a['bot']
            if crit_bot != zero:
                c2 += 1
            if (sum(a['byD'].values()) == 0) != zero:
                c5 += 1
            if a['crit'] and not zero:
                c6 += 1
            if zero:
                nz += 1
            else:
                nnz += 1
    dd = ", ".join("%d:%d" % kv for kv in sorted(dist.items())[:5])
    print("  %4d %3d %7d | %7d  %-22s %6d | %6d | %12d | %6d"
          % (t, r, nsh, max(dist), dd[:22], mxD, c2, c5, c6))
    TOT['sh'] += nsh
    TOT['c2'] += c2
    TOT['c5'] += c5
    TOT['c6'] += c6
    TOT['z'] += nz
    TOT['nz'] += nnz

print("")
print("  C2 exceptions to  [Phi]_bot = 0  <=>  every signed count by D vanishes : %d of %d"
      % (TOT['c2'], TOT['sh']))
print("  C5 decoy (one global sum of w instead of split by D): wrong %d times -- must be nonzero"
      % TOT['c5'])
print("  C6 criterion-holding shapes where [Phi]_bot survives: %d (must be 0)" % TOT['c6'])
print("  C7 non-vacuity: [Phi]_bot vanishes on %d shapes and survives on %d." % (TOT['z'], TOT['nz']))
print("")
if TOT['c2'] == 0 and TOT['c6'] == 0 and TOT['c5'] > 0 and TOT['z'] > 0 and TOT['nz'] > 0:
    print("  THE BOTTOM STRATUM IS A SIGNED COUNT.  [Phi]_bot = 0 iff, for every multiset D of")
    print("  consecutive-gap differences, the transversals in Gbot realising D carry weights w")
    print("  summing to zero.  No Schur rigidity, no polynomial: just w = +-1 added up.")
else:
    print("  SOMETHING FAILED -- read the columns.")
print("")
print("DONE")
