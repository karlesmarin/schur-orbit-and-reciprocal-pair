# -*- coding: utf-8 -*-
# ¿SE CUMPLE LA HIPOTESIS (A) DE LA CONSULTA EXTERNA?  -- 15 de agosto de 2026.
#
# LA AFIRMACION QUE SE VERIFICA.  La consulta propone que la saturacion de CANCELAN no necesita nada
# tropical: bastaria con que la expansion se escriba
#
#       F_s(z) = sum_{omega in Omega} c_omega z^{a_omega + s b_omega}
#
# con Omega FINITO E INDEPENDIENTE DE s y los c_omega INDEPENDIENTES DE s.  De ahi sale un lema
# elemental -- dos funciones afines se cruzan a lo sumo una vez, Omega es finito, luego pasado cierto
# s el patron se congela -- y la saturacion seria INEVITABLE, no un hallazgo.
#
# Es una afirmacion falsable y se verifica antes de creerla.  En nuestra expansion los atomos son
#
#       (transversal T, S subset de posiciones, permutacion de S, permutacion de S^c)
#
# con coeficiente  w(T) * split_sign(S,Sc) * sgn(perm_A) * sgn(perm_B)  y exponente suma de +-beta.
# Los EXPONENTES son afines en j de oficio, porque los beta lo son.  Asi que (A) se cumple si y solo
# si NO cambian con j:
#
#   A1  el numero de transversales
#   A2  el multiconjunto de signos w(T)
#   A3  el numero de atomos y el multiconjunto de coeficientes
#
# Si los tres se quedan quietos, (A) vale y el lema afin explica la saturacion.  Si alguno se mueve,
# su explicacion NO se aplica a esta familia y hay que decirlo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python affine_A_check.py

from collections import Counter
from itertools import combinations

from second_stratum import setup, all_transversals
from depth_histogram import alt, split_sign


def atomos(beta, t, r):
    """(n_transversales, Counter(signos w), n_atomos, Counter(coeficientes)) o None."""
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(beta, cl, r, t)
    ws = Counter()
    coefs = Counter()
    n_at = 0
    idx = tuple(range(2 * r))
    for x in tr:
        w, T = x[2], x[1]
        ws[w] += 1
        for S in combinations(idx, r):
            Sc = tuple(a for a in idx if a not in S)
            A = alt([T[a] for a in S], r)
            B = alt([-T[a] for a in Sc], r)
            base = w * split_sign(S, Sc, r)
            for ka, ca in A.items():
                for kb, cb in B.items():
                    coefs[base * ca * cb] += 1
                    n_at += 1
    return len(tr), ws, n_at, coefs


def extremos_S(beta, t):
    cl, E, Cd = setup(beta, t)
    S = sorted({v for k in E for v in Cd[k]})
    return S[-1], S[0]


# DOS REGIMENES, y la distincion es mia y no de la consulta: si el paso s NO es multiplo de t, al
# mover el extremo se le cambia la CLASE DE RESIDUO en los j impares, y entonces (A) falla por culpa
# del paso elegido y no por la matematica.  Se prueban los dos para no cargarle a su lema un fallo
# nuestro.  s=26 no es multiplo de 4; s=28 si.  s=21 no es multiplo de 6; s=24 si.
CASOS = [(4, 2, 26, (18, 17, 11, 8, 7, 6, 1, 0), 6),
         (4, 2, 28, (18, 17, 11, 8, 7, 6, 1, 0), 6),
         (6, 3, 21, (16, 15, 14, 13, 11, 6, 5, 4, 3, 2, 1, 0), 4),
         (6, 3, 24, (16, 15, 14, 13, 11, 6, 5, 4, 3, 2, 1, 0), 4)]

print("=" * 96)
print("HIPOTESIS (A): ¿son Omega y los c_omega independientes de j?")
print("=" * 96)
for (t, r, s, seed, J) in CASOS:
    print("\n--- t=%d r=%d  paso s=%d  semilla %s" % (t, r, s, seed))
    print("      j |    W | transversales | multiconj. de w | atomos | multiconj. de coeficientes")
    hi, lo = extremos_S(seed, t)
    ref = None
    veredicto = "(A) SE CUMPLE"
    for j in range(J + 1):
        b = tuple(sorted([(x + s * j if x == hi else (x - s * j if x == lo else x)) for x in seed],
                         reverse=True))
        a = atomos(b, t, r)
        if a is None:
            print("    %3d | %4d | (no aplica)" % (j, b[0] - b[-1]))
            veredicto = "INTERRUMPIDO"
            break
        ntr, ws, nat, coefs = a
        firma = (ntr, tuple(sorted(ws.items())), nat, tuple(sorted(coefs.items())))
        if ref is None:
            ref = firma
        elif firma != ref:
            veredicto = "(A) FALLA en j=%d" % j
        print("    %3d | %4d | %13d | %-15s | %6d | %s"
              % (j, b[0] - b[-1], ntr, dict(ws), nat,
                 str(dict(sorted(coefs.items())))[:44]))
    print("    ->", veredicto)

print()
print("=" * 96)
print("  Si (A) se cumple, la saturacion la explica el lema afin y NO es un hallazgo nuestro.")
print("  Si (A) falla, esa explicacion no se aplica a esta familia.")
print("=" * 96)
