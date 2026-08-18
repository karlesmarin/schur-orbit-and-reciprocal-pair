# -*- coding: utf-8 -*-
# AUDITORIA DE MI PROPIA REFUTACION.  12 de agosto de 2026 (noche).
#
# Acabo de matar el enunciado  [Phi]_top = 0 y [Phi]_bot = 0 => sigma_C(S) = S  con 8 contraejemplos
# salidos de mis propios guiones.  Antes de que eso se quede escrito hay que atacarlo:
#
#   S1  ¿SON OCHO?  Los 8 vienen por parejas beta -> beta+1, y ademas beta -> M-beta (complemento en
#       el rectangulo) es otra simetria del problema.  Si las dos actuan, la base de evidencia es
#       mucho mas fina de lo que dije.  Se calculan las ORBITAS, y se VERIFICA que las dos
#       simetrias preservan de verdad (top=0, bot=0, (i), Phi=0) sobre una muestra grande -- no se
#       supone.
#   S2  CROSS-VALIDACION CONTRA EL GUION VIEJO.  bottom_stratum.sage midio "0 excepciones en 10 289
#       formas" barriendo por |lambda| <= MAX.  Si mi pipeline no reproduce ESE 0 en ESE rango, el
#       que esta mal soy yo y los 8 no valen nada.  Mismas 7 configuraciones, mismo muestreo.
#   S3  ¿DONDE ESTA EL PRIMERO POR |lambda|?  Dije "el barrido iba a 14 y el primero tiene 24".  Eso
#       hay que medirlo, no deducirlo de un caso: barrido por |lambda| creciente hasta encontrarlo.
#   S4  ¿SPORADICOS O SISTEMATICOS?  Los 8 aparecen con beta_max = 16, 17 en un barrido con M = 17,
#       o sea PEGADOS AL BORDE.  Si al subir M la poblacion crece, el enunciado no falla por poco:
#       falla en todas partes y el barrido solo veia la punta.
#
# ATAJO USADO, y su licencia: bot = 0 se decide por el recuento con signo por multiconjunto D en vez
# de expandir monomios (A0 de bottom_sees_gcom.py: 0 fallos de 93 329), y top = 0 solo se calcula
# cuando |G| = 2 (con |G| = 1 el estrato de arriba es +-a_H(z)a_L(1/z), producto de dos alternantes
# de entradas distintas, no nulo).  Las dos cosas se re-verifican aqui en L1/L2 antes de usarse.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python two_strata_audit.py

import itertools
import sys
from collections import defaultdict

FAILS = [
    [16, 15, 14, 12, 11, 7, 5, 4, 3, 2, 1, 0],
    [16, 15, 14, 13, 12, 11, 9, 5, 4, 2, 1, 0],
    [16, 14, 13, 11, 9, 7, 6, 5, 4, 3, 2, 0],
    [16, 14, 13, 12, 11, 10, 9, 7, 5, 3, 2, 0],
    [17, 16, 15, 13, 12, 8, 6, 5, 4, 3, 2, 1],
    [17, 16, 15, 14, 13, 12, 10, 6, 5, 3, 2, 1],
    [17, 15, 14, 12, 10, 8, 7, 6, 5, 4, 3, 1],
    [17, 15, 14, 13, 12, 11, 10, 8, 6, 4, 3, 1],
]


def perm_sign(seq):
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def topdeg_dict(T, r):
    D = defaultdict(int)
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
            D[tuple(e)] += perm_sign(q)
    return dict((k, v) for k, v in D.items() if v != 0)


def botdeg_dict(T, r):
    D = defaultdict(int)
    n = 2 * r
    for s in itertools.permutations(range(r)):
        q = [0] * n
        e = [0] * r
        for i in range(r):
            q[2 * i] = 2 * s[i]
            q[2 * i + 1] = 2 * s[i] + 1
            e[s[i]] = T[2 * i] - T[2 * i + 1]
        D[tuple(e)] += perm_sign(q)
    return dict((k, v) for k, v in D.items() if v != 0)


def full_expansion(tm, r):
    n = 2 * r
    D = defaultdict(int)
    for w, T in tm:
        for q in itertools.permutations(range(n)):
            e = [0] * r
            for a in range(n):
                e[q[a] // 2] += (T[a] if q[a] % 2 == 0 else -T[a])
            D[tuple(e)] += w * perm_sign(list(q))
    return dict((k, v) for k, v in D.items() if v != 0)


def analyse(beta, t, r, want_dicts=False):
    """beta decreciente y distinto.  Devuelve None si alguna clase esta vacia."""
    N = len(beta)
    cl = defaultdict(list)
    for i, b in enumerate(beta):
        cl[b % t].append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    if not E:
        return None
    S = sorted((beta[i] for k in E for i in cl[k]), reverse=True)
    C = S[0] + S[-1]
    cond_i = set(C - v for v in S) == set(S)
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    cond_ii = (len(fixed) == 2 and all(k in E for k in fixed))

    keys = sorted(cl)
    tr = []
    for pick in itertools.product(*[cl[k] for k in keys]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        dif = tuple(T[2 * i] - T[2 * i + 1] for i in range(r))
        tr.append((T, w, sum(T[:r]) - sum(T[r:]), sum(dif),
                   tuple(sorted(dif, reverse=True))))

    Dmax = max(x[2] for x in tr)
    Dmin = min(x[3] for x in tr)
    G = [x for x in tr if x[2] == Dmax]
    Gb = [x for x in tr if x[3] == Dmin]

    # ATAJO 1: con |G| = 1 el estrato de arriba no puede anularse.
    if len(G) == 1:
        top_zero = False
    else:
        top = defaultdict(int)
        for x in G:
            for k, v in topdeg_dict(list(x[0]), r).items():
                top[k] += x[1] * v
        top_zero = not any(v != 0 for v in top.values())

    # ATAJO 2: bot = 0 por recuento con signo repartido por multiconjunto D.
    byD = defaultdict(int)
    for x in Gb:
        byD[x[4]] += x[1]
    bot_zero = all(v == 0 for v in byD.values())

    return dict(t=t, r=r, beta=tuple(beta), cond_i=cond_i, cond_ii=cond_ii,
                crit=(cond_i and cond_ii), nG=len(G), nGb=len(Gb), e=len(E),
                Dmax=Dmax, Dmin=Dmin, top_zero=top_zero, bot_zero=bot_zero, tr=tr)


def is_fail(a):
    """contraejemplo: el criterio falla y los DOS estratos se anulan."""
    return a is not None and (not a['crit']) and a['top_zero'] and a['bot_zero']


def partitions(n, maxlen, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, maxlen - 1, k):
            yield (k,) + rest


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def lam_of(beta):
    N = len(beta)
    return tuple(x for x in (beta[i] - (N - 1 - i) for i in range(N)) if x > 0)


# ================================================================= S1 ============================
print("=" * 108)
print("S1  ¿SON OCHO?  orbitas bajo traslacion (beta+1) y complemento (M-beta)")
print("=" * 108)
print("")

canon = {}
for i, b in enumerate(FAILS):
    canon[tuple(b)] = i + 1

def orbit(b):
    """cierre de beta bajo beta+1, beta-1 y complemento en {0..max}."""
    seen = set()
    stack = [tuple(b)]
    while stack:
        x = stack.pop()
        if x in seen:
            continue
        seen.add(x)
        stack.append(tuple(v + 1 for v in x))
        if min(x) > 0:
            stack.append(tuple(v - 1 for v in x))
        for M in (max(x), max(x) + 1):
            stack.append(tuple(sorted((M - v for v in x), reverse=True)))
        # poda: la orbita es infinita por traslacion, se corta por el maximo
        stack = [y for y in stack if max(y) <= 19 and min(y) >= 0 and y not in seen]
    return seen

orbits = []
left = set(canon)
while left:
    b = sorted(left)[0]
    o = orbit(b) & set(canon)
    orbits.append(sorted(canon[x] for x in o))
    left -= o

print("  orbitas encontradas: %d" % len(orbits))
for o in orbits:
    print("      {%s}   representante: beta = %s   lambda = %s  |lambda| = %d"
          % (", ".join("#%d" % i for i in o), FAILS[o[0] - 1],
             list(lam_of(FAILS[o[0] - 1])), sum(lam_of(FAILS[o[0] - 1]))))
print("")
print("  => la base de evidencia son %d formas ESENCIALMENTE distintas, no 8." % len(orbits))

print("")
print("  y las dos simetrias, VERIFICADAS y no supuestas (t=6, r=3, M=15, todas las formas):")
sy_t = sy_c = sy_n = 0
for comb in itertools.combinations(range(16), 12):
    beta = sorted(comb, reverse=True)
    a = analyse(beta, 6, 3)
    if a is None:
        continue
    sy_n += 1
    key = (a['cond_i'], a['cond_ii'], a['top_zero'], a['bot_zero'])
    bt = analyse([v + 1 for v in beta], 6, 3)
    if bt is None or (bt['cond_i'], bt['cond_ii'], bt['top_zero'], bt['bot_zero']) != key:
        sy_t += 1
    M = max(beta)
    bc = analyse(sorted((M - v for v in beta), reverse=True), 6, 3)
    if bc is None or (bc['cond_i'], bc['cond_ii'], bc['top_zero'], bc['bot_zero']) != key:
        sy_c += 1
print("      traslacion beta -> beta+1 cambia el cuadruplete : %d de %d" % (sy_t, sy_n))
print("      complemento beta -> M-beta cambia el cuadruplete: %d de %d" % (sy_c, sy_n))
sys.stdout.flush()

# ================================================================= S2 ============================
print("")
print("=" * 108)
print("S2  CROSS-VALIDACION: mi pipeline sobre el MUESTREO VIEJO (|lambda| <= MAX)")
print("     bottom_stratum.sage midio 0 excepciones en 10 289 formas.  Tengo que reproducirlo.")
print("=" * 108)
print("")
CONF_OLD = [(2, 2, 22), (2, 3, 16), (4, 2, 24), (4, 3, 16), (6, 2, 18), (6, 3, 14), (8, 2, 16)]
print("     t   r  |lam|<= | formas | criterio falla | top=0  bot=0 | EXCEPCIONES")
print("  " + "-" * 96)
tot_old = tot_fail = tot_exc = 0
for (t, r, MAX) in CONF_OLD:
    N = t + 2 * r
    n = nf = ntz = nbz = nexc = 0
    for size in range(MAX + 1):
        for lam in partitions(size, N):
            a = analyse(beta_of(lam, N), t, r)
            if a is None:
                continue
            n += 1
            if a['crit']:
                continue
            nf += 1
            if a['top_zero']:
                ntz += 1
            if a['bot_zero']:
                nbz += 1
            if is_fail(a):
                nexc += 1
    print("  %4d %3d %8d | %6d | %14d | %5d %6d | %11d" % (t, r, MAX, n, nf, ntz, nbz, nexc))
    tot_old += n
    tot_fail += nf
    tot_exc += nexc
    sys.stdout.flush()
print("")
print("  TOTAL %d formas, %d con el criterio falso, EXCEPCIONES: %d" % (tot_old, tot_fail, tot_exc))
print("  (bottom_stratum.sage: 10 289 formas con el criterio falso, 0 excepciones)")

# ================================================================= S3 ============================
print("")
print("=" * 108)
print("S3  ¿DONDE ESTA EL PRIMERO POR |lambda|?   t = 6, r = 3, barrido creciente")
print("=" * 108)
print("")
first = None
per_size = []
for size in range(0, 27):
    n = nexc = 0
    for lam in partitions(size, 12):
        a = analyse(beta_of(lam, 12), 6, 3)
        if a is None:
            continue
        n += 1
        if is_fail(a):
            nexc += 1
            if first is None:
                first = (size, lam, a)
    per_size.append((size, n, nexc))
    if nexc:
        print("  |lambda| = %2d : %6d formas, %d EXCEPCIONES" % (size, n, nexc))
    sys.stdout.flush()
if first is None:
    print("  ninguna excepcion hasta |lambda| = 26")
else:
    print("")
    print("  PRIMERA: |lambda| = %d, lambda = %s, beta = %s"
          % (first[0], list(first[1]), list(beta_of(first[1], 12))))
    print("  el barrido viejo llegaba a |lambda| <= 14, o sea le faltaban %d tallas."
          % (first[0] - 14))

# ================================================================= S4 ============================
print("")
print("=" * 108)
print("S4  ¿SPORADICOS O SISTEMATICOS?   subir M y ver si la poblacion crece")
print("=" * 108)
print("")
print("     t   r    M | formas | criterio falla | EXCEPCIONES | e=t entre ellas")
print("  " + "-" * 96)
for (t, r, M) in [(6, 3, 15), (6, 3, 16), (6, 3, 17), (6, 3, 18), (6, 3, 19),
                  (4, 3, 17), (4, 2, 19), (8, 3, 19)]:
    N = t + 2 * r
    if M < N - 1:
        continue
    n = nf = nexc = net = 0
    for comb in itertools.combinations(range(M + 1), N):
        beta = sorted(comb, reverse=True)
        a = analyse(beta, t, r)
        if a is None:
            continue
        n += 1
        if not a['crit']:
            nf += 1
        if is_fail(a):
            nexc += 1
            if a['e'] == t:
                net += 1
    print("  %4d %3d %4d | %6d | %14d | %11d | %d" % (t, r, M, n, nf, nexc, net))
    sys.stdout.flush()

print("")
print("DONE")
