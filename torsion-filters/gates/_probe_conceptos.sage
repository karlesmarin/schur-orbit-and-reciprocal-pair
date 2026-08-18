# -*- coding: utf-8 -*-
# DOS CONCEPTOS QUE EL PAPER USA SIN COMPROBAR.   16 de agosto de 2026.
#
# (A) .Es  W^1 = { w : w(u) estrictamente H-dominante }  el W^1 ESTANDAR?
#     El W^1 de Kostant/GKRS se define con rho_G solo y NO depende de Lambda.  El mio se define con
#     u = Lambda + rho_G y SI parece depender.  Si son el mismo conjunto, mi criterio tiene que dar
#     el MISMO conjunto de elementos de Weyl para dos Lambda distintos.  Test falsable:
#     A1  el conjunto de w obtenido con Lambda = 0 y con Lambda arbitrario coincide.
#     A2  su cardinal es |W_G|/|W_H|.
#     Si A1 falla, el indice del multiplete depende de Lambda y la formula que escribimos no es la
#     de GKRS.
#
# (B) .Es  q^B_t  una proyeccion sobre un anillo de RANGO UNO?
#     Solo lo es si TODO superviviente pliega al mismo punto del alcove -- el vacio.  Si alguno
#     pliega a otro punto, el anillo tiene rango > 1 y "tau en {0,+-1}" esta diciendo otra cosa.
#     B1  para todo eta con tau^B_t(eta) != 0, el vector plegado es el de eta = 0.
#     B2  cuantos puntos del alcove se alcanzan en total.
#     Y el contraste con el par, donde el paper ya lo tiene probado por otra via.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage _probe_conceptos.sage

import itertools
from collections import Counter


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


def W1(Lam, mp, r):
    """{ w : w(u) estrictamente H-dominante }, devuelto como conjunto de (perm, signos)."""
    Rp = mp + r
    u2 = [2 * int(Lam[i]) + 2 * (Rp - i) - 1 for i in range(Rp)]   # 2u, entero
    out = set()
    for perm in itertools.permutations(range(Rp)):
        base = [u2[perm[i]] for i in range(Rp)]
        for eps in itertools.product((1, -1), repeat=Rp):
            w = [base[i] * eps[i] for i in range(Rp)]
            if not (all(w[i] > w[i + 1] for i in range(mp - 1)) and w[mp - 1] > 0):
                continue
            f = w[mp:]
            if r >= 2 and not (all(f[i] > f[i + 1] for i in range(r - 1))
                               and f[r - 2] > abs(f[r - 1])):
                continue
            out.add((perm, eps))
    return out


print("=" * 100)
print("(A)  .es W^1 independiente de Lambda?")
print("=" * 100)
a1 = a1n = a2 = 0
for (t, r) in [(3, 2), (5, 2), (7, 2), (3, 3), (5, 3)]:
    mp = (t - 1) // 2
    Rp = mp + r
    base = W1([0] * Rp, mp, r)
    esperado = 2 * binomial(Rp, mp)
    a2 += 1 if len(base) == esperado else 0
    iguales = 0
    pruebas = 0
    for Lam in itertools.product(range(4), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        pruebas += 1
        if W1(list(Lam), mp, r) == base:
            iguales += 1
    a1 += iguales
    a1n += pruebas
    print("  t=%d r=%d :  |W^1| = %3d  (esperado 2C(R',m') = %3d)   mismo conjunto que Lambda=0 : %3d de %3d"
          % (t, r, len(base), esperado, iguales, pruebas))
print("")
print("  A1  W^1 no depende de Lambda : %d de %d" % (a1, a1n))
print("  A2  |W^1| = 2 C(R',m')       : %d de 5" % a2)

print("")
print("=" * 100)
print("(B)  .pliega TODO superviviente al vacio?")
print("=" * 100)


def plegar(v, t):
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def car(typ, rk, mu):
    W = WeylCharacterRing("%s%d" % (typ, rk))
    el = W(W.space().from_vector(vector([Integer(v) for v in mu])))
    d = {}
    for wt, mult in el.weight_multiplicities().items():
        k = tuple(int(v) for v in wt.to_vector())
        d[k] = d.get(k, 0) + int(mult)
    return d


def tau(typ, rk, eta, t):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car(typ, rk, eta).items():
        s += mult * z ** (sum((i + 1) * int(wt[i]) for i in range(rk)) % t)
    return int(QQ(s)) if s in QQ else None


for (t, typ) in [(3, "B"), (5, "B"), (7, "B"), (9, "B"), (4, "C"), (6, "C"), (8, "C"), (10, "C")]:
    rk = (t - 1) // 2 if typ == "B" else (t - 2) // 2
    if rk < 1:
        continue
    puntos = Counter()
    vivos = 0
    for eta in itertools.product(range(5), repeat=rk):
        if any(eta[j] < eta[j + 1] for j in range(rk - 1)):
            continue
        tv = tau(typ, rk, list(eta), t)
        if not tv:
            continue
        vivos += 1
        if typ == "B":
            A = [2 * int(eta[j]) + 2 * (rk - j) - 1 for j in range(rk)]
        else:
            A = [int(eta[j]) + rk - j for j in range(rk)]
        cl = tuple(sorted(plegar(v, t)[0] for v in A))
        # el punto del alcove: las clases plegadas, ordenadas.  El vacio es (1,2,...,rk).
        puntos[cl] += 1
    vacio = tuple(range(1, rk + 1))
    print("  t=%2d tipo %s rango %d :  %3d supervivientes | puntos del alcove alcanzados: %d | "
          "todos al vacio: %s"
          % (t, typ, rk, vivos, len(puntos), str(set(puntos) == {vacio})))
print("")
print("=" * 100)
print("DONE")
