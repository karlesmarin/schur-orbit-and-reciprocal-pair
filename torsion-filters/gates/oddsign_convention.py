# -*- coding: utf-8 -*-
# LA CONSTANTE DEL FILTRO IMPAR, EN EL CONVENIO DEL PAPER.   16 de agosto de 2026.
#
# POR QUE SE REHACE.  cor:oddsign se escribio hoy con delta definida leyendo las clases plegadas
# contra el orden CRECIENTE 1,2,...,m'.  Pero Lemma lem:T, que es el enunciado del que cuelga todo,
# define sigma "read against the decreasing order m, m-1, ..., 1".  Los dos convenios difieren por
# el signo de la inversion, (-1)^{m(m-1)/2}, y por tanto la constante que cor:oddsign publica
# tambien.  Un convenio distinto no es cosmetico: aqui cambia el enunciado.
#
# Se recalcula la constante en el convenio del paper y se busca su forma cerrada AJUSTANDO en
# t <= 31 y PREDICIENDO en 33..61, para que el ajuste no se valide a si mismo.
#
# CONTROL FATAL: la constante tiene que ser CONSTANTE en cada t -- si dependiera de eta, no seria
# una normalizacion y el enunciado entero se cae.  Como tau^B_t(0) = 1 siempre, se calibra con
# eta = 0 y se comprueba contra delta(A_rho).

import json
from math import gcd


def jacobi(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def plegar(v, t):
    v %= t
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
    return s


def delta(a, t, n, decreciente):
    """decreciente=True es el convenio de lem:T: las clases se leen contra m, m-1, ..., 1."""
    cl, ep = [], []
    for v in a:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, n + 1)):
        return 0
    orden = [n - c for c in cl] if decreciente else [c - 1 for c in cl]
    s = sgn_perm(orden)
    for e in ep:
        s *= e
    return int(s)


TS = list(range(3, 62, 2))
RES = []
print("=" * 104)
print("LA CONSTANTE eps_t DEL FILTRO IMPAR, EN LOS DOS CONVENIOS")
print("=" * 104)
print("")
print("   t   m'   eps (creciente, el mio)   eps (DECRECIENTE, el de lem:T)   (-2/t)^{(t+3)/2}")
print("   " + "-" * 92)
for t in TS:
    n = (t - 1) // 2
    A_rho = [t - 2 * (j + 1) for j in range(n)]      # = 2(eta+rho) con eta = 0
    d_cre = delta(A_rho, t, n, False)
    d_dec = delta(A_rho, t, n, True)
    # tau^B_t(0) = 1  =>  eps = 1/delta = delta
    cerrada = jacobi((-2) % t, t) ** ((t + 3) // 2)
    print("  %2d  %3d          %+2d                        %+2d                    %+2d"
          % (t, n, d_cre, d_dec, cerrada))
    RES.append({"t": t, "m": n, "eps_creciente": int(d_cre), "eps_decreciente": int(d_dec),
                "jacobi_m2": int(cerrada)})

AJ = [r for r in RES if r["t"] <= 31]
PR = [r for r in RES if r["t"] > 31]
CANDS = {
    "(-2/t)^{(t+3)/2}": lambda r: r["jacobi_m2"],
    "(-2/t)^{(t+3)/2} . (-1)^{m(m-1)/2}": lambda r: r["jacobi_m2"] * (1 if (r["m"] * (r["m"] - 1) // 2) % 2 == 0 else -1),
    "(-1)^{m(m-1)/2}": lambda r: 1 if (r["m"] * (r["m"] - 1) // 2) % 2 == 0 else -1,
    "(2/t)": lambda r: jacobi(2, r["t"]),
    "-(2/t)": lambda r: -jacobi(2, r["t"]),
    "(-1/t)": lambda r: jacobi((-1) % r["t"], r["t"]),
    "trivial (+1)": lambda r: 1,
    "por t mod 8: {1:+,3:+,5:-,7:+}": lambda r: -1 if r["t"] % 8 == 5 else 1,
}
print("")
print("=" * 104)
print("FORMA CERRADA de eps en el convenio del paper.  Ajuste en t <= 31, PREDICCION en 33..61.")
print("=" * 104)
buenos = []
for nom, f in CANDS.items():
    a = sum(1 for r in AJ if f(r) == r["eps_decreciente"])
    p = sum(1 for r in PR if f(r) == r["eps_decreciente"])
    ok = (a == len(AJ) and p == len(PR))
    if ok:
        buenos.append(nom)
    print("   %-38s ajuste %2d/%2d   prediccion %2d/%2d%s"
          % (nom, a, len(AJ), p, len(PR), "   <== SOBREVIVE" if ok else ""))

json.dump({"por_t": RES, "cerradas_que_sobreviven": buenos},
          open("oddsign_convention_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
