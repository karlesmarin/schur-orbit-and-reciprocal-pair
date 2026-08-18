# -*- coding: utf-8 -*-
# (L1) POR SUBCONJUNTOS: la formula que la prueba, puesta a prueba.   16 de agosto de 2026.
#
# EL ARGUMENTO, escrito antes de medir.  Queremos  c(Lambda, mu) = sum_eta B^odd tau^B(eta), que es
# el coeficiente de o_mu(z) en  o_Lambda(p_t)  con  p_t = (mu_t, z^{+-}).  Por Weyl en B_{R'},
#
#     o_Lambda(x) = det( w_i^{A_j} - w_i^{-A_j} ) / det( w_i^{P_j} - w_i^{-P_j} ),
#     A_j = 2 Lambda_j + 2(R'-j) + 1   (impares, DISTINTOS),   P_j = 2(R'-j)+1,
#
# con w_i^2 = x_i, y el punto p_t tiene  x = (zeta^1..zeta^{m'}, z_1..z_r).  Desarrollando por
# LAPLACE a lo largo de las m' primeras filas -- las de las raices de la unidad --:
#
#     N_Lambda = sum_{S, |S|=m'}  +-  det[bloque zeta con exponentes A_S] . det[bloque libre con A_{S^c}]
#
# y el primer factor es el numerador de Weyl de B_{m'} en el elemento de torsion, que por nuestro
# propio lema vale CERO salvo que A_S sea REGULAR mod t, y entonces vale +- el denominador.  Luego
#
#     o_Lambda(p_t) = sum_{S regular} +- (caracter del factor libre con exponentes A_{S^c}).
#
# Y como los A_j son DISTINTOS, subconjuntos distintos dan multiconjuntos complementarios distintos,
# luego pesos distintos: cada mu recibe A LO SUMO UN termino.  De ahi  c en {0,+-1}  --  (L1).
#
# Ademas eso da una FORMULA CERRADA, no solo la cota:
#
#     c(Lambda, mu) = +- 1  si existe S con |S| = m', A_S regular mod t, y el complemento da mu;
#                       0   si no.
#
# LO QUE SE MIDE, y es lo que decide si el argumento vale
#   S1  para cada Lambda: los subconjuntos S regulares, y el peso mu^+ que da cada complemento.
#   S2  ¿es INYECTIVA la aplicacion S -> mu^+?  Si dos S dieran el mismo mu, el argumento falla.
#   S3  el soporte predicho {mu : existe S regular} contra el soporte MEDIDO {mu : c != 0}.
#   S4  el numero de subconjuntos regulares contra el numero de mu con c != 0.
#
# CONTROLES
#   C0  c se calcula por la ruta de siempre (branching + tau), que no sabe nada de subconjuntos.
#   C1  se cuentan los fallos de los dos tipos por separado.
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage L1_subset_formula.sage

import json
import sys
import itertools
from collections import defaultdict

_CH = {}
def car(typ, rk, mu):
    key = (typ, rk, tuple(int(v) for v in mu))
    if key not in _CH:
        W = WeylCharacterRing("%s%d" % (typ, rk))
        el = W(W.space().from_vector(vector([Integer(v) for v in mu])))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + int(mult)
        _CH[key] = d
    return _CH[key]


_BR = {}
def branch(Rp, mp, rr, Lam):
    key = (Rp, mp, rr, tuple(int(v) for v in Lam))
    if key not in _BR:
        W = WeylCharacterRing("B%d" % Rp)
        X = WeylCharacterRing("B%dxD%d" % (mp, rr))
        br = branching_rule("B%d" % Rp, "B%dxD%d" % (mp, rr), "orthogonal_sum")
        el = W(W.space().from_vector(vector([Integer(v) for v in Lam])))
        d = {}
        for wt, c in el.branch(X, rule=br).monomial_coefficients().items():
            v = [int(u) for u in wt.to_vector()]
            d[(tuple(v[:mp]), tuple(v[mp:]))] = int(c)
        _BR[key] = d
    return _BR[key]


_TAU = {}
def tauB(mp, eta, t):
    key = (mp, tuple(int(v) for v in eta), t)
    if key not in _TAU:
        K = CyclotomicField(t)
        z = K.gen()
        s = K(0)
        for wt, mult in car("B", mp, eta).items():
            s += mult * z ** (sum((i + 1) * wt[i] for i in range(mp)) % t)
        _TAU[key] = QQ(s) if s in QQ else None
    return _TAU[key]


def regular(A, t, k):
    """A_S regular mod t: ninguno cero, y las clases +- distintas."""
    c = [v % t for v in A]
    if 0 in c:
        return False
    cl = [min(v, t - v) for v in c]
    return len(set(cl)) == k


print("=" * 122)
print("(L1) POR SUBCONJUNTOS:  ¿reproduce la formula el c medido?")
print("=" * 122)

RES = []
for (t, r) in [(3, 2), (5, 2), (3, 3), (7, 2)]:
    mp = (t - 1) // 2
    Rp = mp + r
    # Lambdas de una caja
    LAMS = []
    for k in range(0, 3 * t + 2):
        for e in Partitions(k, max_length=Rp):
            LAMS.append(tuple(list(e) + [0] * (Rp - len(e))))
    n = iny_fallos = sop_fallos = 0
    npares = 0
    ejemplo = None
    for Lam in LAMS:
        A = [2 * Lam[j] + 2 * (Rp - j - 1) + 1 for j in range(Rp)]
        # medido
        cmed = defaultdict(lambda: 0)
        for (eta, mu), cc in branch(Rp, mp, r, Lam).items():
            tv = tauB(mp, eta, t)
            if tv is None or tv == 0:
                continue
            mplus = tuple(list(mu[:-1]) + [abs(mu[-1])])
            cmed[mplus] += int(cc) * int(tv)
        cmed = {k2: v for k2, v in cmed.items() if v != 0}
        # predicho: subconjuntos S de tamaño mp con A_S regular
        pred = defaultdict(lambda: 0)
        for S in itertools.combinations(range(Rp), mp):
            AS = [A[i] for i in S]
            if not regular(AS, t, mp):
                continue
            comp = sorted([A[i] for i in range(Rp) if i not in S], reverse=True)
            # el peso libre: quitar el rho del factor libre en las mismas coordenadas dobladas.
            # rho de B_r doblado seria (2(r-j)+1); se resta y se divide por 2.
            mu = tuple((comp[j] - (2 * (r - j - 1) + 1)) // 2 for j in range(r))
            if any(v < 0 for v in mu) or list(mu) != sorted(mu, reverse=True):
                continue
            pred[mu] += 1
        n += 1
        npares += len(cmed)
        # S2: inyectividad
        if any(v > 1 for v in pred.values()):
            iny_fallos += 1
            if ejemplo is None:
                ejemplo = ("inyectividad", [int(v) for v in Lam], dict(pred))
        # S3: soportes
        if set(pred) != set(cmed):
            sop_fallos += 1
            if ejemplo is None:
                ejemplo = ("soporte", [int(v) for v in Lam],
                           sorted(set(pred) - set(cmed))[:3], sorted(set(cmed) - set(pred))[:3])
    print("")
    print("  t=%d r=%d  (R'=%d, m'=%d):  %d Lambda, %d pares (Lambda,mu) con c != 0"
          % (t, r, Rp, mp, n, npares))
    print("     S2  inyectividad S -> mu : fallos en %d de %d Lambda" % (iny_fallos, n))
    print("     S3  soporte predicho == soporte medido : fallos en %d de %d" % (sop_fallos, n))
    if ejemplo:
        print("     primer fallo: %s" % str(ejemplo))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "n_Lambda": int(n), "n_pares": int(npares),
                "fallos_inyectividad": int(iny_fallos), "fallos_soporte": int(sop_fallos)})

print("")
print("=" * 122)
print("  LECTURA, escrita ANTES de correr:")
print("   * si la inyectividad y el soporte salen sin fallos, el argumento de Laplace vale y (L1)")
print("     queda PROBADA -- y con formula cerrada, que es mas de lo que pedia.")
print("   * si el soporte falla, el desarrollo tiene mas terminos de los que creo y hay que mirar")
print("     que se me escapa: seguramente el bloque libre no da el caracter que supongo.")
json.dump(RES, open("L1_subset_formula_DUMP.json", "w"), indent=1)
print("=" * 122)
print("DONE")
