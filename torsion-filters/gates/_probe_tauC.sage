# -*- coding: utf-8 -*-
# .Tiene el filtro de tipo C una constante de normalizacion, como el de tipo B?
# El paper afirma, en la prueba de prop:galoissign, que en tipo C la constante es +1.  even_transversal
# encontro un testigo en contra: t=6, eta=(3,3) da tau=-1 y delta=+1.  Aqui se mira si es una
# constante por t, o si la discrepancia depende de eta -- que serian dos diagnosticos muy distintos.

def plegar(v, t):
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def sgn_perm(perm):
    n, s, visto = len(perm), 1, [False] * len(perm)
    for i in range(n):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            s = -s
    return int(s)


def delta_C(a, t, m):
    cl, ep = [], []
    for v in a:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, m + 1)):
        return 0
    s = sgn_perm([c - 1 for c in cl])
    for e in ep:
        s *= e
    return int(s)


def tauC(eta, t, m):
    K = CyclotomicField(t)
    z = K.gen()
    W = WeylCharacterRing("C%d" % m)
    el = W(W.space().from_vector(vector([Integer(v) for v in eta])))
    s = K(0)
    for wt, mult in el.weight_multiplicities().items():
        k = tuple(int(v) for v in wt.to_vector())
        s += int(mult) * z ** (sum((i + 1) * k[i] for i in range(m)) % t)
    return int(QQ(s)) if s in QQ else None


import itertools
print("=" * 96)
print("CONSTANTE DE TIPO C:  .es tau^C = delta, o hay un factor?")
print("=" * 96)
for t in [4, 6, 8, 10]:
    m = (t - 2) // 2
    razones = {}
    n = 0
    ejemplos = []
    for eta in itertools.product(range(6), repeat=m):
        if any(eta[j] < eta[j + 1] for j in range(m - 1)):
            continue
        tv = tauC(list(eta), t, m)
        a = [int(eta[j]) + m - j for j in range(m)]
        dv = delta_C(a, t, m)
        if tv is None:
            continue
        n += 1
        if dv == 0 and tv == 0:
            continue
        if dv == 0 or tv == 0:
            razones.setdefault("soporte distinto", []).append((list(eta), tv, dv))
            continue
        razones.setdefault(int(tv) * int(dv), []).append((list(eta), tv, dv))
    print("")
    print("  t=%2d  C_%d :  %d etas" % (t, m, n))
    for k, v in sorted(razones.items(), key=lambda kv: str(kv[0])):
        print("     tau/delta = %-16s : %3d casos, p.ej. %s" % (k, len(v), v[:2]))
print("")
print("=" * 96)
print("DONE")
