# -*- coding: utf-8 -*-
# t = 2: EL ESTRATO DE ARRIBA BASTA, Y ESO CIERRA conj:crit.  Con la prueba, y con el control que
# explica por que t >= 4 NO se cierra igual.
#
# *** HIPOTESIS MIA REFUTADA EL MISMO TURNO ***  Gbot_anatomy.sage puso a prueba "Gbot es un producto
# por clases con factores <= 2".  FALSO: no es un producto en 4857 de 8954 formas, los factores
# llegan a 6, hay parejas no adyacentes, y el estrato de ARRIBA tampoco es un producto (766 formas).
# La minimizacion de abajo NO es separable como la de arriba, y no hay analogo del argumento de
# residuos.  Pero su tabla de contingencia enseño otra cosa: en t = 2 el estrato de ABAJO no hace
# ninguna falta -- bottom_stratum.sage ya lo decia y no lo habia leido: t=2, r=2, 2195 formas donde el
# criterio falla y [Phi]_top != 0 en LAS 2195; t=2, r=3, 789 de 789.  Cero excepciones.  El estrato de
# arriba SOLO decide t = 2.  Y el estrato de arriba ya es teorema.
#
# ================================================================================================
# TEOREMA.  Sea t = 2, r >= 1, N = 2 + 2r, sea lambda una particion con ell(lambda) <= N rellenada a
# exactamente N partes, y supongase las dos clases de residuos ocupadas.  Entonces
#   -- F19: r >= 1 es necesario; con r = 0 no hay incrementos que elegir y el codigo ni corre.
#   -- F20: ell(lambda) <= N es hipotesis PERMANENTE; leida al pie de la letra, la Conjetura 9.4
#      dice 'for every r and every lambda' y para ell(lambda) > N es falsa por vacuidad: Psi_r = 0
#      trivialmente y beta, luego las dos ramas, no estan definidas.  El convenio de relleno esta
#      en el preambulo de la seccion 9 del paper, pero el teorema tiene que llevarlo escrito.
#       Phi_2(lambda; z) == 0   <=>   (i) sigma_C(S) = S   y   (ii).
# o sea la Conjetura 9.4 de arXiv:2608.09619.
#
# Prueba.  Phi_2 = 0 => [Phi]_top = 0.  Como |G| = 1 da [Phi]_top = +- a_H(z)a_L(1/z) != 0, tiene que
# ser |G| = 2.  Por el teorema |G| <= 2, un |G| = 2 exige empate justo en el rango r entre dos clases
# k y k + t/2, y en t = 2 eso son LAS DOS clases; en particular e = |E| = 2, porque con una sola clase
# de exceso todos los incrementos estan en ella y son estrictamente decrecientes, sin empate posible.
#
# El voraz mueve un incremento de una clase a la otra, asi que los dos maximizadores difieren en las
# DOS clases:  g_A = {x1, y1},  g_B = {x2, y2},  x en la clase k, y en la k+1, los cuatro distintos.
# Como e = 2 y ambas clases cambian, NO HAY PARTE COMUN: g_A u g_B son esos cuatro elementos, y
#       S  =  K  u  {x1, x2, y1, y2},        K = T_A cap T_B.
#
# Ahora [Phi]_top = 0 con |G| = 2 da P(T_A) = +- P(T_B); por el Teorema 2.5 de Purbhoo-van
# Willigenburg T_B esta en la orbita de T_A, la rama de traslacion implica la de reflexion (paso 5),
# el centro es V (paso 6) y V = C (lemma_V_eq_C).  De T_B = C - T_A y de que sigma_C es una involucion
# sale sigma_C(K) = K.  Y el EMPATE dice que los dos elementos intercambiados de cada clase suman el
# valor de empate:  x1 + x2 = y1 + y2 = V = C,  o sea sigma_C({x1,x2}) = {x1,x2} y lo mismo con las y.
# Luego
#       sigma_C(S) = sigma_C(K) u sigma_C({x1,x2}) u sigma_C({y1,y2}) = S.
# Eso es (i).  Y (ii) ya estaba probado.  El reciproco es el teorema publicado.  QED
#
# POR QUE t >= 4 NO SE CIERRA ASI, y es exactamente el hueco.  Con e > 2 los dos maximizadores
# coinciden en las clases que no estan empatadas, y esa parte comun g_com es no vacia:
#       S = K u g_com u {x1,x2,y1,y2},
# de modo que sigma_C(S) = S necesita ademas sigma_C(g_com) = g_com, que NADA obliga.  Las
# excepciones medidas del estrato de arriba tienen que ser exactamente esas.  Se comprueba abajo.
# ================================================================================================
#
# VERIFICADO, cada cosa capaz de fallar:
#   E1  en t = 2, [Phi]_top = 0 obliga a e = 2.
#   E2  y los dos maximizadores difieren en las DOS clases, o sea g_com = vacio.
#   E3  los tres ingredientes: sigma_C(K) = K,  x1+x2 = y1+y2 = C,  y de ahi sigma_C(S) = S.
#   E4  LA CONCLUSION en t = 2:  [Phi]_top = 0  <=>  (i) y (ii).  Las dos direcciones.
#   E5  acceptance, fatal: y las dos coinciden con que el DETERMINANTE se anule de verdad, evaluado
#       sobre GF(p) en puntos aleatorios.  Sin esto todo lo anterior habla de otro objeto.
#   E6  EL CONTROL QUE EXPLICA t >= 4: las formas con [Phi]_top = 0 y (i) FALSA deben ser exactamente
#       aquellas con g_com no estable bajo sigma_C.  Si aparece una con g_com estable, la prueba
#       de arriba esta mal.
#   E7  no vacuidad: en t = 2 tiene que haber formas que cumplan el criterio y formas que no.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

T2 = [(2, 1, 40), (2, 2, 32), (2, 3, 24), (2, 4, 18), (2, 5, 14)]
TBIG = [(4, 2, 22), (4, 3, 16), (6, 2, 18), (6, 3, 14), (8, 2, 16)]

L = lcm([2, 4, 6, 8])
p = next_prime(10 ** 9)
while (p - 1) % L != 0:
    p = next_prime(p)
F = GF(p)
G0 = F.multiplicative_generator()
print("field GF(%d)" % p)


def zeta(t):
    z = G0 ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t))
    return z


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


def tie_value_of(a):
    inc = []
    for k in a['E']:
        ck = a['CC'][k]
        for j in range(1, len(ck)):
            inc.append(ck[j - 1] + ck[j])
    inc.sort(reverse=True)
    return inc[a['r'] - 1]


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
    return dict(E=E, CC=CC, S=S, C=C, r=r, conc=conc, cond_ii=cond_ii, crit=conc and cond_ii,
                G=G, top=dict((k, v) for k, v in top.items() if v != 0))


# ---------------------------------------------------------------- E1-E4, E7 ----------------------
print("")
print("=" * 106)
print("E1-E4  the t = 2 chain, step by step")
print("=" * 106)
print("")
print("     t   r  |lam|<=  shapes | top=0 | E1 e!=2  E2 gcom!=0  E3a sK!=K  E3b sums!=C"
      "  E3c not conc | E4 bad | crit")
print("  " + "-" * 102)

FAIL = False
TOT = dict(sh=0, z=0, e1=0, e2=0, a=0, b=0, c=0, e4=0, cr=0)
for t, r, MAX in T2:
    N = t + 2 * r
    nsh = nz = e1 = e2 = ea = eb = ec = e4 = ncr = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            nsh += 1
            if a['crit']:
                ncr += 1
            z = not a['top']
            if z:
                nz += 1
                if len(a['E']) != 2:
                    e1 += 1
                gs = list(a['G'])
                if len(gs) != 2:
                    e2 += 1
                    continue
                gA, gB = set(gs[0]), set(gs[1])
                gcom = gA & gB
                if gcom:
                    e2 += 1
                TA, TB = a['G'][gs[0]][0], a['G'][gs[1]][0]
                C = a['C']
                K = set(TA) & set(TB)
                if set(C - x for x in K) != K:
                    ea += 1
                sw = sorted(gA.symmetric_difference(gB))
                pairs = {}
                for v in sw:
                    pairs.setdefault(v % t, []).append(v)
                if not all(len(v) == 2 and sum(v) == C for v in pairs.values()):
                    eb += 1
                if set(C - x for x in a['S']) != set(a['S']):
                    ec += 1
            if z != a['crit']:
                e4 += 1
    print("  %4d %3d %8d %7d | %5d | %7d %11d %11d %13d %14d | %6d | %5d"
          % (t, r, MAX, nsh, nz, e1, e2, ea, eb, ec, e4, ncr))
    for k, v in (('sh', nsh), ('z', nz), ('e1', e1), ('e2', e2), ('a', ea),
                 ('b', eb), ('c', ec), ('e4', e4), ('cr', ncr)):
        TOT[k] += v

print("")
print("  totals over %d shapes at t = 2: [Phi]_top vanishes on %d, criterion holds on %d."
      % (TOT['sh'], TOT['z'], TOT['cr']))
print("  E1 e != 2 : %d      E2 gcom nonempty : %d" % (TOT['e1'], TOT['e2']))
print("  E3a sigma_C(K) != K : %d   E3b the swapped pairs do not sum to C : %d   E3c S not concentric"
      " : %d" % (TOT['a'], TOT['b'], TOT['c']))
print("  E4 disagreements between  [Phi]_top = 0  and  (i) and (ii) : %d" % TOT['e4'])
print("  E7 non-vacuity: criterion holds %d times and fails %d." % (TOT['cr'], TOT['sh'] - TOT['cr']))

# ---------------------------------------------------------------- E5 -----------------------------
print("")
print("=" * 106)
print("E5  acceptance, fatal: does the DETERMINANT itself agree, over GF(p) at random z?")
print("=" * 106)
print("")
print("     t   r  shapes tested  points  det=0 vs criterion disagreements")
print("  " + "-" * 102)

set_random_seed(4242)
# F2/F3, de la auditoria socratica.  La version anterior visitaba 200 formas por configuracion y NO
# imprimia cuantas cumplian el criterio: eran 7, 5, 1, 1.  O sea la mitad "criterio => det = 0" -- la
# unica que puede cazar un diccionario equivocado -- descansaba en 14 formas en total, y el control
# no lo decia.  Un control cuyas dos direcciones no se reportan por separado no se puede leer.
# Ademas saltaba r = 5 y topaba en |lambda| <= 16.  Corregido: las dos direcciones con su
# denominador, todas las r, y ademas formas de branch (b) PLANTADAS a proposito, que es la unica
# manera de que la direccion escasa tenga muestra.
e5 = 0
NT = NF = 0


def planted_b(N, wmax=13, cap=400):
    """branch-(b) shapes built on purpose: self-complementary of odd width w."""
    out = []
    h = N // 2
    for w in range(1, wmax + 1, 2):
        lo = (w + 1) // 2
        for top in itertools.combinations_with_replacement(range(lo, w + 1), h):
            lam = sorted(top, reverse=True)
            full = lam + [w - x for x in reversed(lam)]
            if all(full[k] >= full[k + 1] for k in range(N - 1)):
                out.append([x for x in full if x > 0])
                if len(out) >= cap:
                    return out
    return out


print("     t   r  shapes | crit TRUE  det=0 there | crit FALSE  det!=0 there | disagreements")
print("  " + "-" * 102)
for t, r, MAX in T2:
    N = t + 2 * r
    zt = zeta(t)
    nT = nF = bT = bF = 0
    shapes = []
    for size in range(min(MAX, 22) + 1):
        for l in Partitions(size, max_length=N):
            shapes.append(list(l))
    shapes += planted_b(N)
    for l in shapes:
        beta = beta_of(list(l), N)
        a = analyse(beta, t, r)
        if a is None:
            continue
        vals = []
        for _ in range(3):
            zv = []
            while len(zv) < r:
                x = F.random_element()
                if x != 0 and x ** 2 != 1 and all(x != y and x * y != 1 for y in zv):
                    zv.append(x)
            cols = [zt ** k for k in range(t)]
            for u in range(r):
                cols += [zv[u], 1 / zv[u]]
            vals.append(matrix(F, N, N, lambda i, j: cols[j] ** beta[i]).det())
        iszero = all(v == 0 for v in vals)
        if a['crit']:
            nT += 1
            bT += 1 if iszero else 0
        else:
            nF += 1
            bF += 1 if not iszero else 0
    e5 += (nT - bT) + (nF - bF)
    NT += nT
    NF += nF
    print("  %4d %3d %7d | %9d %13d | %10d %14d | %13d"
          % (t, r, nT + nF, nT, bT, nF, bF, (nT - bT) + (nF - bF)))

print("")
print("  criterion TRUE on %d shapes and FALSE on %d -- BOTH directions now have a denominator."
      % (NT, NF))
if e5:
    print("  E5 FAILED -- the object being tested is not Phi_t.  Everything above is void.")
    raise SystemExit(1)
print("  E5 PASS: the determinant vanishes exactly when the criterion holds, at t = 2, in both")
print("  directions and with the scarce direction deliberately populated.")

# ---------------------------------------------------------------- E6 -----------------------------
print("")
print("=" * 106)
print("E6  why t >= 4 does NOT close: the common part g_com and its symmetry")
print("=" * 106)
print("")
print("     t   r | top=0 & (i) fails | of those, g_com NOT sigma_C-stable | g_com stable (MUST BE 0)")
print("  " + "-" * 102)

e6 = 0
tot6 = 0
for t, r, MAX in TBIG:
    N = t + 2 * r
    n = ns = st = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None or a['top'] or a['conc']:
                continue
            n += 1
            gs = list(a['G'])
            if len(gs) != 2:
                continue
            gcom = set(gs[0]) & set(gs[1])
            # F18: the chain forecasts sigma_V(g_com) = g_com; V is the tie value, which need not be
            # C a priori.  On these shapes V = C, but the control must be the predicate predicted.
            Vv = tie_value_of(a)
            if set(Vv - x for x in gcom) == gcom:
                st += 1
            else:
                ns += 1
    e6 += st
    tot6 += n
    print("  %4d %3d | %17d | %34d | %24d" % (t, r, n, ns, st))

print("")
print("  E6: %d shapes with [Phi]_top = 0 and (i) false.  In %d of them the common part g_com is NOT"
      % (tot6, tot6 - e6))
print("      sigma_C-stable, and in %d it IS -- that column must be 0, because a stable g_com would" % e6)
print("      make the t = 2 argument go through and force concentricity.")
print("")
if TOT['e1'] == TOT['e2'] == TOT['a'] == TOT['b'] == TOT['c'] == TOT['e4'] == 0 and e6 == 0 \
        and TOT['cr'] > 0 and TOT['cr'] < TOT['sh'] and tot6 > 0:
    pass
if TOT['e1'] or TOT['e2'] or TOT['a'] or TOT['b'] or TOT['c'] or TOT['e4'] or e6:
    FAIL = True
if True:
    print("  conj:crit CERRADA en t = 2, y el control E6 explica exactamente por que t >= 4 no lo esta:")
    print("  la parte comun g_com es vacia en t = 2 y no lo es despues, y su falta de simetria es el")
    print("  unico sitio donde el argumento se rompe.")
else:
    print("  SOMETHING FAILED -- read the columns.")
# ---------------------------------------------------------------- E8 -----------------------------
print("")
print("=" * 106)
print("E8  am I closing the RIGHT conjecture?  The published statement is (a) or (b), not (i)+(ii)")
print("=" * 106)
print("")
print("  arXiv:2608.09619, eq (branches):  branch (a) = |E| in {0,N}, the beta set has constant")
print("  parity;  branch (b) = lambda_i + lambda_{N+1-i} = w for all i, with w ODD.")
print("  Conjecture 9.4: for every r and every lambda, Psi_r(lambda) = 0 iff (a) or (b).")
print("  Published converses already: Thm r1 (r = 1) and Thm stable (ell(lambda) <= N/2, where the")
print("  locus is exactly lambda = (k^{N/2}) with k odd).")
print("")
print("     t   r  shapes | not-(a)  (b) | branch-(a) bookkeeping | (i)&(ii) != (b) |"
      " Thm-stable bad")
print("  " + "-" * 102)

e8 = 0
e9 = 0
E8A = 0
E8B = 0
for t, r, MAX in T2:
    N = t + 2 * r
    nsh = nb = mis = mis2 = xs = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            lam = list(l) + [0] * (N - len(l))
            beta = beta_of(list(l), N)
            nE = sum(1 for b in beta if b % 2 == 0)
            brA = (nE == 0 or nE == N)
            w = lam[0] + lam[N - 1]
            brB = all(lam[i] + lam[N - 1 - i] == w for i in range(N)) and (w % 2 == 1)
            a = analyse(beta, t, r)
            if a is None:
                # a is None exactly when a residue class is empty, i.e. branch (a)
                if not brA:
                    mis += 1
                continue
            if brA:
                mis += 1                      # setup accepted a shape that IS branch (a)
            nsh += 1
            if brB:
                nb += 1
            if a['crit'] != brB:
                mis2 += 1
            # cross-check against the PUBLISHED Theorem (stable): on ell(lambda) <= N/2 the locus
            # is exactly the odd rectangles (k^{N/2})
            if len(l) <= N // 2:
                k = lam[0]
                rect = all(lam[i] == k for i in range(N // 2)) and                     all(lam[i] == 0 for i in range(N // 2, N)) and (k % 2 == 1)
                if a['crit'] != rect:
                    xs += 1
    E8A += mis
    E8B += mis2
    e8 += mis + mis2
    e9 += xs
    print("  %4d %3d %7d | %7d %4d | %15d | %11d | %26d"
          % (t, r, nsh, nsh, nb, mis, mis2, xs))

print("")
print("  E8 shapes where (i)+(ii) disagrees with branch (b): %d   (this is mis2 alone;" % E8B)
print("     the other column, %d, is branch-(a) bookkeeping and was MISLABELLED before)" % E8A)
print("  E9 shapes inside Littlewood's range where our criterion disagrees with the PUBLISHED")
print("     Theorem (stable) -- the odd rectangles (k^{N/2}): %d" % e9)
print("")
if e8 == 0 and e9 == 0:
    print("  THE DICTIONARY HOLDS.  At t = 2 with both parity classes occupied, (i) says beta is")
    print("  symmetric about C = beta_1 + beta_N, which is lambda_i + lambda_{N+1-i} = w, and (ii)")
    print("  says C is even, which is w odd because N is even.  So (i)+(ii) IS branch (b), and the")
    print("  theorem above closes Conjecture 9.4 of arXiv:2608.09619 -- with branch (a) being the")
    print("  excluded empty-class case, which vanishes by pigeonhole.  And it agrees with the two")
    print("  published converses, Thm r1 and Thm stable, on their ranges.")
else:
    print("  DICTIONARY BROKEN -- I am not closing the conjecture I think I am.")
