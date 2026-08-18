# -*- coding: utf-8 -*-
# LA LEY DEL PESO SUPERIOR EN EL IMPAR: ¿que rho lleva el desplazamiento?   16 de agosto de 2026.
#
# DE DONDE SALE.  Cruzar dos formulas del paper que nunca se habian mirado juntas.
#
#   (A)  la ley del peso superior, medida SOLO en el caso par:
#            mu_max = top Newt(N_beta)  -  sigma_r,     sigma_r = 2 rho_{C_r} + (t-1) (1,...,1)
#
#   (B)  el pipeline impar, donde el bloque LIBRE no es C_r sino D_r.
#
# Si el desplazamiento es "el rho del grupo libre" -- que es como esta escrito y como se justifica --
# entonces en el impar tiene que ser  2 rho_{D_r} + (t-1)(1,...,1),  y no el de C_r.  Los dos
# difieren en un vector UNIFORME (2,2,...,2), asi que el experimento separa limpio:
#
#        2 rho_{C_r} = (2r, 2r-2, ..., 2)          sigma_C = (N-1, N-3, ..., N-2r+1)
#        2 rho_{D_r} = (2r-2, 2r-4, ..., 0)        sigma_D = (N-3, N-5, ..., t-1)
#
# Y si NINGUNO de los dos acierta, la ley es de las coordenadas de exponente y no del grupo, que
# tambien es una respuesta -- y cambiaria como hay que enunciarla en el paper.
#
# LO QUE SE MIDE, forma a forma del impar
#   L1  top Newt(N_beta): el vertice dominante del NUMERADOR, calculado expandiendo el bialternante
#       (no el politopo abstracto: el polinomio de verdad).
#   L2  mu_max^+ por la ruta de siempre (pelado en D_r, plegado a O(2r)).
#   L3  ¿top - sigma_C = mu_max?  ¿top - sigma_D = mu_max?  Las dos, por separado.
#
# CONTROLES
#   C0  el caso PAR tiene que dar sigma_C acertando: es la ley ya publicada.  Si fallara aqui, el
#       instrumento esta roto y no se mira nada mas.
#   C1  se imprime el top y las dos predicciones, no solo el veredicto.
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage mumax_odd.sage

import json
import sys


def bialt(beta, tt, nvar, num=True):
    """el numerador det(x^beta) (num=True) o el objeto Phi (num=False), como dict de exponentes."""
    Nn = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    Lr = LaurentPolynomialRing(K, nvar, 'z')
    zs = Lr.gens()
    xx = [Lr(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(Nn - 1, -1, -1))

    def alt(expo):
        return matrix(Lr, Nn, Nn, lambda i, j: xx[i] ** expo[j]).determinant()

    q = alt(list(beta)) if num else None
    if not num:
        den = alt(delta)
        if den == 0:
            return None
        q = alt(list(beta)) / den
        try:
            q = Lr(q)
        except Exception:
            return "NO-POL"
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        k = tuple(int(v) for v in e) if hasattr(e, '__iter__') else (int(e),)
        if c != 0:
            out[k] = c
    return out


def top_dominante(P, rr):
    """el vertice dominante de mayor grado del soporte."""
    dom = [e for e in P if list(e) == sorted(e, reverse=True) and min(e) >= 0]
    if not dom:
        return None
    return max(dom, key=lambda e: (sum(e), e))


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


def dominante(e, typ):
    f = list(e)
    if typ in ("B", "C"):
        return f == sorted(f, reverse=True) and min(f) >= 0
    if len(f) < 2:
        return f[0] >= 0
    return all(f[i] >= f[i + 1] for i in range(len(f) - 2)) and f[-2] >= abs(f[-1])


def pelar(P, typ, rk, tope=8000):
    P = {e: QQ(c) for e, c in P.items() if c != 0}
    out = {}
    for _ in range(tope):
        if not P:
            return out, {}
        dom = [e for e in P if dominante(e, typ)]
        if not dom:
            return out, P
        mu = max(dom, key=lambda e: (sum(abs(v) for v in e), e))
        c = P[mu]
        out[mu] = out.get(mu, 0) + c
        for k, v in car(typ, rk, mu).items():
            nv = P.get(k, 0) - c * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P


def mumax_plus(P1, r, typ):
    A1, rest = pelar(P1, typ, r)
    A1 = {k: QQ(v) for k, v in A1.items() if v != 0}
    if rest or not A1:
        return None, None
    Ap = {}
    for mu, c in A1.items():
        Ap[tuple(list(mu[:-1]) + [abs(mu[-1])])] = c
    S = list(Ap)
    maxi = [m for m in S if not any(n != m and all(sum(n[:k + 1]) >= sum(m[:k + 1])
                                                   for k in range(len(m))) for n in S)]
    if len(maxi) != 1:
        return None, None
    return maxi[0], Ap[maxi[0]]


def betas(N, tope):
    return [tuple(sorted(c, reverse=True)) + (0,) for c in Combinations(range(1, tope + 1), N - 1)]


print("=" * 126)
print("LA LEY DEL PESO SUPERIOR EN EL IMPAR:  ¿2 rho_{C_r} o 2 rho_{D_r}?")
print("=" * 126)

RES = []
for (t, r, tope, tipo) in [(2, 2, 9, "C"), (4, 2, 11, "C"), (3, 2, 9, "D"), (5, 2, 10, "D")]:
    N = t + 2 * r
    sig_C = [N - 1 - 2 * k for k in range(r)]              # 2 rho_{C_r} + (t-1)
    sig_D = [N - 3 - 2 * k for k in range(r)]              # 2 rho_{D_r} + (t-1)
    print("")
    print("  t=%d  r=%d  N=%d   (%s)   sigma_C = %s   sigma_D = %s"
          % (t, r, N, "PAR, control" if tipo == "C" else "IMPAR", str(tuple(sig_C)), str(tuple(sig_D))))
    print("  beta                       | top Newt(N_beta) | mu_max^+   | top - sigma_C | top - sigma_D | acierta")
    print("  " + "-" * 118)
    okC = okD = n = 0
    muestra = 0
    for b in betas(N, tope):
        P1 = bialt(b, t, r, num=False)
        if P1 in (None, "NO-POL") or not P1:
            continue
        mm, A = mumax_plus(P1, r, tipo)
        if mm is None:
            continue
        NB = bialt(b, t, r, num=True)
        top = top_dominante(NB, r)
        if top is None:
            continue
        n += 1
        pC = tuple(top[k] - sig_C[k] for k in range(r))
        pD = tuple(top[k] - sig_D[k] for k in range(r))
        if pC == mm:
            okC += 1
        if pD == mm:
            okD += 1
        if muestra < 8:
            muestra += 1
            print("  %-26s | %-16s | %-10s | %-13s | %-13s | %s"
                  % (str(b), str(top), str(mm), str(pC), str(pD),
                     ("C" if pC == mm else "") + ("D" if pD == mm else "") or "NINGUNO"))
            sys.stdout.flush()
    print("  ...")
    print("  FORMAS %d  |  sigma_C acierta %d  |  sigma_D acierta %d" % (n, okC, okD))
    RES.append({"t": int(t), "r": int(r), "tipo": tipo, "n": int(n),
                "sigma_C_acierta": int(okC), "sigma_D_acierta": int(okD),
                "sigma_C": [int(v) for v in sig_C], "sigma_D": [int(v) for v in sig_D]})
    sys.stdout.flush()

print("")
print("=" * 126)
print("  LECTURA, escrita ANTES de correr:")
print("   * si en el par gana C y en el impar gana D -> el desplazamiento es el rho DEL GRUPO LIBRE,")
print("     y la ley se enuncia una sola vez para las dos paridades.")
print("   * si gana C en los dos -> el desplazamiento es de las coordenadas de exponente y NO del")
print("     grupo, y la justificacion que da el paper esta mal aunque la formula acierte.")
print("   * si no gana ninguno en el impar -> la ley no se traslada y hay que medir su forma alli.")
json.dump(RES, open("mumax_odd_DUMP.json", "w"), indent=1)
print("=" * 126)
print("DONE")
