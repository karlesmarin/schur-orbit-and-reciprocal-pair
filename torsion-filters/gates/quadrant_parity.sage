# -*- coding: utf-8 -*-
# *** AVISO RETIRADO: ESTE GUION ESTABA BIEN. ***
#
# Llevo un aviso, puesto por mi, diciendo que sus numeros eran falsos porque contradecian a
# quadrant_mechanism.sage.  Era al reves.  parity2.sage recalcula lo mismo por DOS rutas
# independientes, pasa un test de aceptacion sobre el caso impreso a mano, y reproduce
# exactamente estos numeros: 93 formas con eps no constante y solo 2 con eps = -1.  El que esta
# mal es quadrant_mechanism.sage.  Dejo la traza del aviso equivocado en lugar de borrarla.
#
# The parity that turns the reflection into a criterion.
#
# Psi ~ sum_T sgn(T) A(T) over the complements of transversals, and A(C-T) = (-1)^r A(T).  So with
# r even the sum collapses exactly when the reflection T -> C-T carries OPPOSITE Laplace signs, and
# doubles when it carries the same.  That is the whole difference between the 30 vanishing shapes
# and the 277 that are reflective and do not vanish.  Write
#
#     eps  =  sgn(C-T) * sgn(T)          (+1 same sign, -1 opposite)
#
# and the question is a closed formula for eps in terms of beta.
#
# This computes eps, checks it is CONSTANT over the terms -- if it is not, the whole framing is
# wrong and nothing else matters -- and then scores candidate closed forms against it.  Each
# candidate must agree on every shape; agreeing on the vanishing ones alone proves nothing, since
# eps = -1 is what vanishing means.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

t, r = 4, 2
N = t + 2 * r
MAX = 22


def shuffle_sign(pick, rest):
    perm = list(pick) + list(rest)
    sg = 1
    for a in range(len(perm)):
        for b in range(a + 1, len(perm)):
            if perm[a] > perm[b]:
                sg = -sg
    return sg


def setup(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for k, b in enumerate(beta):
        cls.setdefault(b % t, []).append(k)
    if len(cls) < t:
        return None
    exc = sorted([i for v in cls.values() if len(v) >= 2 for i in v])
    return beta, cls, exc


def centres(beta, exc):
    S = set(beta[i] for i in exc)
    return [C for C in range(2 * max(S) + 1) if set(C - b for b in S) == S]


def eps_of(beta, cls, C):
    """sgn(C-T)*sgn(T), or None if not constant / not closed."""
    pos = {}
    for i, b in enumerate(beta):
        pos[b] = i
    keys = sorted(cls)
    vals = None
    for pick in itertools.product(*[cls[k] for k in keys]):
        rest = tuple(sorted(set(range(N)) - set(pick)))
        s1 = shuffle_sign(sorted(pick), rest)
        try:
            rest2 = tuple(sorted(pos[C - beta[i]] for i in rest))
        except KeyError:
            return None
        pick2 = tuple(sorted(set(range(N)) - set(rest2)))
        s2 = shuffle_sign(pick2, rest2)
        e = s1 * s2
        if vals is None:
            vals = e
        elif vals != e:
            return 0                      # not constant
    return vals


print("=" * 78)
print("The Laplace parity of the reflection,  t = 4, r = 2,  |lambda| <= %d" % MAX)
print("=" * 78)

data = []
nonconst = 0
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        lam = list(l)
        st = setup(lam)
        if st is None:
            continue
        beta, cls, exc = st
        if not exc:
            continue
        for C in centres(beta, exc):
            e = eps_of(beta, cls, C)
            if e is None:
                continue
            if e == 0:
                nonconst += 1
                continue
            sizes = sorted(len(v) for v in cls.values())
            m = len(exc)
            k = sum(1 for v in cls.values() if len(v) >= 2)
            singles = [i for v in cls.values() if len(v) == 1 for i in v]
            below = sum(1 for i in singles
                        if min(beta[j] for j in exc) < beta[i] < max(beta[j] for j in exc))
            data.append(dict(lam=lam, beta=beta, C=C, eps=e, m=m, k=k,
                             sizes=tuple(sizes), s=t - k, below=below))
            break

print("")
print("  shapes with a reflective excess union : %d" % len(data))
print("  CONTROL, eps not constant over terms  : %d   (must be 0)" % nonconst)
print("  of them eps = -1 (the sum cancels)    : %d" % sum(1 for d in data if d['eps'] == -1))

CAND = {
    "(-1)^(m(m-1)/2)": lambda d: (-1) ** (d['m'] * (d['m'] - 1) // 2),
    "(-1)^(m(m-1)/2 + s)": lambda d: (-1) ** (d['m'] * (d['m'] - 1) // 2 + d['s']),
    "(-1)^(m/2)": lambda d: (-1) ** (d['m'] // 2),
    "(-1)^k": lambda d: (-1) ** d['k'],
    "(-1)^(m(m-1)/2 + below)": lambda d: (-1) ** (d['m'] * (d['m'] - 1) // 2 + d['below']),
    "(-1)^(m(m-1)/2 + k)": lambda d: (-1) ** (d['m'] * (d['m'] - 1) // 2 + d['k']),
    "(-1)^below": lambda d: (-1) ** d['below'],
}
print("")
print("  %-32s %s" % ("candidate closed form", "agrees with eps"))
print("  " + "-" * 54)
for name, f in CAND.items():
    ok = sum(1 for d in data if f(d) == d['eps'])
    print("  %-32s %5d / %-5d%s" % (name, ok, len(data), "   <-- EXACT" if ok == len(data) else ""))

print("")
print("  eps against the coarse invariants, to see what it does depend on:")
tab = {}
for d in data:
    tab.setdefault((d['m'], d['k'], d['s'], d['below']), set()).add(d['eps'])
mixed = [key for key, v in tab.items() if len(v) > 1]
print("    cells (m,k,s,below) with BOTH values of eps: %d of %d" % (len(mixed), len(tab)))
for key in sorted(tab)[:14]:
    print("    m=%d k=%d s=%d below=%d  ->  eps in %s" % (key + (sorted(tab[key]),)))

print("")
print("DONE")
