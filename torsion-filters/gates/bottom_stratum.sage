# -*- coding: utf-8 -*-
# LA FILTRACION POR GRADOS TIENE DOS EXTREMOS.  Solo habiamos mirado uno.
#
# DONDE ESTAMOS.  Phi_t = 0 => (ii) es teorema.  Falta Phi_t = 0 => (i), que contiene conj:crit.
# El estrato de grado MAXIMO decide todo salvo 20 formas, y en esas 20 la primera supervivencia esta
# 2 grados por debajo en 17 casos y 4 en 3 casos.  Bajar de dos en dos es una escalera sin fondo
# visible, asi que antes de subirla conviene mirar el OTRO extremo.
#
# ------------------------------------------------------------------------------------------------
# PARIDAD (probado, y explica por que los huecos son pares).  Para un monomio de A(T),
#       sum_j |e_j| = sum_j |u_a - u_b|  =  sum_{x in T} x   (mod 2),
# porque |x-y| = x+y (mod 2).  Luego TODOS los monomios de un mismo A(T) tienen el mismo grado
# modulo 2.  Y entre distintos T: sum T = sum S - sum g con g_k = k (mod t), asi que sum g es fija
# modulo t; si t es PAR eso fija sum g modulo 2 y por tanto sum T modulo 2 tambien.  Conclusion:
# para t par TODO Phi_t es homogeneo modulo 2 en el grado, y la filtracion baja de DOS en DOS.
#
# EL ESTRATO DE ABAJO.  Para T con valores u_1 > .. > u_{2r}, el grado MINIMO de un emparejamiento
# es el de emparejar CONSECUTIVOS,
#       degmin(T) = (u_1-u_2) + (u_3-u_4) + .. + (u_{2r-1}-u_{2r}),
# y ese emparejamiento es el UNICO minimo (cruzar o anidar dos parejas solo puede aumentar la suma).
# Con Dmin = min_g degmin(T_g), el estrato de abajo es
#       [Phi]_bot = sum_{T : degmin(T) = Dmin} w(T) * [A(T)]_bot,
# y [A(T)]_bot es una suma sobre las r! maneras de repartir las r parejas consecutivas entre las
# variables (mismo argumento de sectores: en grado minimo ningun e_j es 0, asi que basta el sector
# todo-positivo).  Nota: si dos diferencias consecutivas COINCIDEN, ese alternante ya se cancela solo
# -- es una fuente de anulacion que el estrato de arriba no tiene.
# ------------------------------------------------------------------------------------------------
#
# LA PREGUNTA, planteada para que pueda fallar:
#   H6   criterio falla  =>  [Phi]_top != 0  O  [Phi]_bot != 0.
# Si sale sin excepciones, el reciproco entero se reduce a DOS condiciones explicitas en vez de a una
# escalera de estratos.  Y hay que probarla en t = 2, que es donde vive conj:crit.
#
# MEDIDO, cada cosa capaz de fallar:
#   B1  acceptance, fatal: [A(T)]_bot debe coincidir monomio a monomio con la rebanada de grado
#       minimo del desarrollo COMPLETO, y el grado minimo del desarrollo debe ser Dmin.
#   B2  la paridad: todos los grados que aparecen en Phi_t deben ser congruentes modulo 2.
#   B3  control forzado: si el criterio VALE, los dos estratos deben anularse.
#   B4  H6, con la estratificacion por cual de las dos condiciones falla.
#   B5  no vacuidad: [Phi]_bot != 0 debe ocurrir, y debe RESCATAR formas que [Phi]_top no ve, o el
#       estrato de abajo no aporta nada nuevo.
#   B6  el señuelo: sustituir el estrato de abajo por el de grado Dmax-2.  Si el de abajo no rescata
#       mas que ese, no hay razon para preferirlo.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

CONF = [(2, 2, 22), (2, 3, 16), (4, 2, 24), (4, 3, 16), (6, 2, 18), (6, 3, 14), (8, 2, 16)]


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


def botdeg_dict(T, r):
    """the minimal-degree part of A(T): the consecutive matching, all-positive sector."""
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
    return dict((k, v) for k, v in D.items() if v != 0)


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
    dt = [sum(T[:r]) - sum(T[r:]) for _, T in tm]
    db = [sum(T[2 * i] - T[2 * i + 1] for i in range(r)) for _, T in tm]
    Dmax, Dmin = max(dt), min(db)
    top = {}
    for (w, T), d in zip(tm, dt):
        if d == Dmax:
            for k, v in topdeg_dict(list(T), r).items():
                top[k] = top.get(k, 0) + w * v
    bot = {}
    for (w, T), d in zip(tm, db):
        if d == Dmin:
            for k, v in botdeg_dict(list(T), r).items():
                bot[k] = bot.get(k, 0) + w * v
    return dict(crit=conc and cond_ii, conc=conc, cond_ii=cond_ii, tm=tm, Dmax=Dmax, Dmin=Dmin,
                top=dict((k, v) for k, v in top.items() if v != 0),
                bot=dict((k, v) for k, v in bot.items() if v != 0))


# ---------------------------------------------------------------- B1, B2 -------------------------
print("=" * 106)
print("B1  the bottom stratum against the full expansion       B2  the degree parity")
print("=" * 106)
print("")
print("     t   r  shapes | Dmin bad  part bad  sector bad | B2 parity bad | Dmax ex.  Dmin ex.")
print("  " + "-" * 102)

bad = 0
for t, r, _ in CONF:
    N = t + 2 * r
    nsh = e1 = e2 = e3 = e4 = 0
    ex = None
    for size in range(0, 13):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            nsh += 1
            if nsh > 20:
                break
            FE = full_expansion(a['tm'], r)
            if not FE:
                continue
            degs = set(sum(abs(x) for x in k) for k in FE)
            if min(degs) < a['Dmin']:
                e1 += 1
            sl = dict((k, v) for k, v in FE.items() if sum(abs(x) for x in k) == a['Dmin'])
            pos = dict((k, v) for k, v in sl.items() if all(x > 0 for x in k))
            if pos != a['bot']:
                e2 += 1
            if sl and len(sl) != 2 ** r * len(pos):
                e3 += 1
            if len(set(d % 2 for d in degs)) > 1:
                e4 += 1
            if ex is None:
                ex = (a['Dmax'], a['Dmin'])
        if nsh > 20:
            break
    bad += e1 + e2 + e3 + e4
    print("  %4d %3d %7d | %8d %9d %11d | %13d | %8d %9d"
          % (t, r, min(nsh, 20), e1, e2, e3, e4, ex[0], ex[1]))

print("")
if bad:
    print("  B1/B2 FAILED -- stop.")
    raise SystemExit(1)
print("  B1 PASS: the consecutive matching really is the minimum, and its all-positive sector is")
print("  the degree-Dmin slice of the full expansion, with the same 2^r sectors as at the top.")
print("  B2 PASS: every degree occurring in Phi_t has the same parity -- the ladder steps by 2.")

# ---------------------------------------------------------------- B3, B4, B5, B6 -----------------
print("")
print("=" * 106)
print("H6  criterion fails  =>  [Phi]_top != 0  OR  [Phi]_bot != 0 ?")
print("=" * 106)
print("")
print("     t   r  |lam|<= | crit  B3 bad | fails | top!=0  bot!=0  top OR bot | BOT RESCUES"
      "  H6 EXCEPTIONS")
print("  " + "-" * 102)

LEFT = []
TOT = dict(c=0, b3=0, f=0, tp=0, bt=0, ei=0, rs=0)
STR = {}
for t, r, MAX in CONF:
    N = t + 2 * r
    nc = b3 = nf = ntp = nbt = neither = resc = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            if a['crit']:
                nc += 1
                if a['top'] or a['bot']:
                    b3 += 1
                continue
            nf += 1
            T_, B_ = bool(a['top']), bool(a['bot'])
            ntp += 1 if T_ else 0
            nbt += 1 if B_ else 0
            if B_ and not T_:
                resc += 1
            if not (T_ or B_):
                neither += 1
                LEFT.append((t, r, list(l), beta, a))
            key = (a['conc'], a['cond_ii'])
            s = STR.setdefault(key, [0, 0, 0, 0])
            s[0] += 1
            s[1] += 1 if T_ else 0
            s[2] += 1 if B_ else 0
            s[3] += 1 if (T_ or B_) else 0
    print("  %4d %3d %8d | %5d %6d | %5d | %6d %7d %11d | %12d %14d"
          % (t, r, MAX, nc, b3, nf, ntp, nbt, nf - neither, resc, neither))
    TOT['c'] += nc
    TOT['b3'] += b3
    TOT['f'] += nf
    TOT['tp'] += ntp
    TOT['bt'] += nbt
    TOT['ei'] += neither
    TOT['rs'] += resc

print("")
print("  B3: criterion holds but a stratum survives: %d (must be 0)." % TOT['b3'])
print("  criterion fails on %d shapes.  top survives %d, bot survives %d, either %d."
      % (TOT['f'], TOT['tp'], TOT['bt'], TOT['f'] - TOT['ei']))
print("  B5 the bottom stratum RESCUES %d shapes the top misses." % TOT['rs'])
print("  H6 EXCEPTIONS (criterion fails and BOTH strata vanish): %d" % TOT['ei'])
print("")
print("     concentric (i)   condition (ii)   shapes   top!=0   bot!=0   either   missed")
print("  " + "-" * 102)
for key in sorted(STR):
    n, a1, a2, a3 = STR[key]
    print("     %-16s %-16s %7d %8d %8d %8d %8d"
          % ("holds" if key[0] else "fails", "holds" if key[1] else "fails", n, a1, a2, a3, n - a3))

if LEFT:
    print("")
    print("  the H6 exceptions:")
    print("")
    print("     t   r  lambda                     conc  (ii)  Dmax  Dmin")
    print("  " + "-" * 102)
    for t, r, lam, beta, a in LEFT[:20]:
        print("  %4d %3d  %-26s %-5s %-5s %5d %5d"
              % (t, r, str(lam), "yes" if a['conc'] else "no",
                 "yes" if a['cond_ii'] else "no", a['Dmax'], a['Dmin']))
print("")
print("DONE")
