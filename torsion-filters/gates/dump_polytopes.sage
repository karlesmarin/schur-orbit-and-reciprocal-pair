# -*- coding: utf-8 -*-
# LOS TRES POLITOPOS DE LA LEY, calculados para dibujarlos.   15 de agosto de 2026.
#
# La ley dice   mu_max = top Newt(N_beta) - 2 rho_{C_r} - (t-1),  y su razon es Ostrowski:
#
#     N_beta = Phi * N_delta     =>     Newt(N_beta) = Newt(Phi) (+) Newt(N_delta)
#
# con (+) la suma de Minkowski.  Para DIBUJARLO hacen falta los tres soportes de verdad, no la ley:
# si se dibujara conv(W(C_r) mu_max) se estaria dibujando el enunciado en vez del dato.  Aqui se
# vuelcan los tres soportes reales, en r=3 para que sean solidos tridimensionales.
#
# CONTROLES
#   C0  FATAL.  Se comprueba la identidad de Minkowski sobre los VERTICES dominantes:
#       top(N_beta) == top(Phi) + top(N_delta).  Si falla, la figura no se dibuja.
#   C1  se vuelca tambien el numero de puntos de cada soporte, para que la figura pueda decir
#       cuantos son y no solo enseñarlos.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage dump_polytopes.sage

import json

t, r = 4, 3
N = t + 2 * r


def soporte(expos_beta):
    """soporte del ALTERNANTE det(x_i^{beta_j}) como polinomio de Laurent en las r variables."""
    K = CyclotomicField(t)
    zeta = K.gen()
    L = LaurentPolynomialRing(K, r, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(t)] + [g ** e for g in zs for e in (1, -1)]
    M = matrix(L, N, N, lambda i, j: x[i] ** expos_beta[j])
    q = M.determinant()
    return sorted({tuple(int(v) for v in e) for e, c in zip(q.exponents(), q.coefficients()) if c != 0})


def soporte_phi(expos_beta):
    K = CyclotomicField(t)
    zeta = K.gen()
    L = LaurentPolynomialRing(K, r, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(t)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(N - 1, -1, -1))
    def alt(ex):
        return matrix(L, N, N, lambda i, j: x[i] ** ex[j]).determinant()
    q = L(alt(list(expos_beta)) / alt(delta))
    return sorted({tuple(int(v) for v in e) for e, c in zip(q.exponents(), q.coefficients()) if c != 0})


def dominante(P):
    """el/los maximos en dominancia entre los representantes dominantes del soporte."""
    dom = {tuple(sorted((abs(v) for v in p), reverse=True)) for p in P}
    def domina(a, b):
        sa = sb = 0
        for u, v in zip(a, b):
            sa += u; sb += v
            if sa < sb:
                return False
        return True
    return sorted(v for v in dom if not any(u != v and domina(u, v) for u in dom))


BETA = (13, 11, 9, 7, 5, 3, 2, 1, 0, -1)     # una forma con Phi != 0 en t=4, r=3
DELTA = tuple(range(N - 1, -1, -1))

print("=" * 96)
print("LOS TRES POLITOPOS  --  t=%d  r=%d  N=%d  beta=%s" % (t, r, N, str(BETA)))
print("=" * 96)

SB = soporte(BETA)
SD = soporte(DELTA)
SP = soporte_phi(BETA)
tb, td, tp = dominante(SB), dominante(SD), dominante(SP)
print("  Newt(N_beta) : %4d puntos, top %s" % (len(SB), tb))
print("  Newt(N_delta): %4d puntos, top %s" % (len(SD), td))
print("  Newt(Phi)    : %4d puntos, top %s" % (len(SP), tp))
print("")
ok = (len(tb) == 1 and len(td) == 1 and len(tp) == 1 and
      tuple(tb[0][i] for i in range(r)) == tuple(tp[0][i] + td[0][i] for i in range(r)))
print("  C0  top(N_beta) == top(Phi) + top(N_delta) : %s" % ("PASA" if ok else "*** FALLA ***"))
print("      %s  ==  %s + %s" % (str(tb[0]) if tb else "-", str(tp[0]) if tp else "-",
                                 str(td[0]) if td else "-"))
sigma = tuple(N - 2 * k + 1 for k in range(1, r + 1))
dosrho = tuple(2 * (r - k + 1) for k in range(1, r + 1))
print("      y el desplazamiento: sigma_r = %s = (t-1) + 2rho_{C_r} = %s"
      % (str(sigma), str(tuple((t - 1) + d for d in dosrho))))

if not ok:
    print("DONE (sin volcado: C0 falla)")
else:
    # int() en CADA entero: los Integer de Sage no son serializables, y es la tercera vez hoy.
    I = lambda v: [int(x) for x in v]
    json.dump({"t": int(t), "r": int(r), "N": int(N), "beta": I(BETA),
               "Nbeta": [I(p) for p in SB], "Ndelta": [I(p) for p in SD],
               "Phi": [I(p) for p in SP],
               "top_Nbeta": I(tb[0]), "top_Ndelta": I(td[0]), "top_Phi": I(tp[0]),
               "sigma": I(sigma)},
              open("polytopes_DUMP.json", "w"))
    print("")
    print("  volcado en polytopes_DUMP.json")
    print("DONE")
