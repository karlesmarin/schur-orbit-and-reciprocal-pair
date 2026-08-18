# -*- coding: utf-8 -*-
# .EN QUE BASE VIVE (L1)?  EL CENTRALIZADOR DEL ELEMENTO DE TORSION.   16 de agosto de 2026.
#
# DE DONDE SALE.  gkrs_L1.sage confirma la identidad GKRS especializada (T1) y deja el hueco entero
# en UN sitio: dividir por  Delta_t = prod_j (z_j^{t/2} - z_j^{-t/2}).  Ese factor es, exactamente, el
# denominador relativo del par  B_r > D_r  dilatado por t.  Y eso apunta a que la base equivocada
# puede ser la nuestra:
#
#   xi_t actua sobre el bloque congelado con autovalores  1, zeta^{+-1}, ..., zeta^{+-m'}.  Su
#   espacio propio de autovalor 1 dentro de  SO(2R'+1)  tiene dimension  1 + 2r.  Luego
#
#         Z_{SO(2R'+1)}(xi_t)  =  T^{m'}  x  SO(2r+1)  =  T^{m'} x B_r ,   NO  D_r.
#
#   Como xi_t es CENTRAL en su centralizador, actua por un escalar en cada constituyente irreducible
#   de  V_Lambda|_{Z(xi_t)},  y por tanto
#
#         ch V_Lambda(xi_t . s)  =  sum_pi  m_pi . lambda_pi(xi_t) . ch pi(s),    s en Z(xi_t).
#
#   Si la teoria del caracter en un elemento de torsion ([NPP25], Kostant) da multiplicidad uno y
#   autovalor racional, los coeficientes son  0, +-1  EN LA BASE DE  B_r.  (L1) estaria enunciada en
#   la base derivada, no en la primitiva.
#
# LO QUE SE MIDE, por Lambda
#   B1  el reparto de  d(Lambda,kappa) := coeficiente de  chi^{B_r}_kappa  en  ch V_Lambda(xi_t, z).
#   B2  el reparto de  c(Lambda,mu)  en la base  D_r  -- la misma serie, otra base, para contrastar.
#   B3  soportes: |supp d|  contra  |supp c|.
#   B4  .es  d  el que se lee del numerador GKRS?  se compara  |supp d|  con  |supp nu| / 2^r.
#
# CONTROLES
#   C0  CONTROL FATAL: el pelado en cada base tiene que cerrar (resto vacio).  Si no cierra, el
#       objeto no es una combinacion de caracteres de ese grupo y el enunciado no existe.
#   C1  SENUELO: se pela tambien en  C_r.  Si el impar viviera en  C_r  no habriamos aprendido nada;
#       tiene que salir peor o no cerrar.
#   C2  n impreso siempre, y los peores casos con su Lambda.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage gkrs_centralizer.sage

import json
import sys
from collections import Counter

CASOS = [(3, 2, 6), (5, 2, 4), (3, 3, 3)]


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


def pelar(P, typ, rk, tope=20000):
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


def serie(Lam, t, mp, r):
    """ch V_Lambda evaluado en (zeta^1,...,zeta^m', z), como dict {exponente en z: coef en QQ}."""
    Rp = mp + r
    K = CyclotomicField(t)
    zeta = K.gen()
    acc = {}
    for wt, mult in car("B", Rp, Lam).items():
        e = K(zeta) ** (sum((i + 1) * int(wt[i]) for i in range(mp)) % t)
        k = tuple(int(wt[mp + j]) for j in range(r))
        acc[k] = acc.get(k, K(0)) + int(mult) * e
    out = {}
    for k, v in acc.items():
        if v not in QQ:
            return None
        if QQ(v) != 0:
            out[k] = QQ(v)
    return out


def dominantes(rk, cota):
    def rec(k, tope):
        if k == 0:
            yield ()
            return
        for a in range(tope, -1, -1):
            for resto in rec(k - 1, a):
                yield (a,) + resto
    return list(rec(rk, cota))


print("=" * 118)
print("EN QUE BASE VIVE (L1):  Z(xi_t) = T^{m'} x B_r,  no D_r")
print("=" * 118)
sys.stdout.flush()

RES = []
for (t, r, cota) in CASOS:
    mp = (t - 1) // 2
    Rp = mp + r
    agg = {"n": 0, "cierraB": 0, "cierraD": 0, "cierraC": 0,
           "dB": Counter(), "dD": Counter(), "dC": Counter(),
           "supB": 0, "supD": 0, "peorB": None, "peorD": None, "no_racional": 0}
    for Lam in dominantes(Rp, cota):
        S = serie(Lam, t, mp, r)
        if S is None:
            agg["no_racional"] += 1
            continue
        agg["n"] += 1
        dB, restB = pelar(dict(S), "B", r)
        dD, restD = pelar(dict(S), "D", r)
        dC, restC = pelar(dict(S), "C", r)
        dB = {k: v for k, v in dB.items() if v != 0}
        dD = {k: v for k, v in dD.items() if v != 0}
        dC = {k: v for k, v in dC.items() if v != 0}
        if not restB:
            agg["cierraB"] += 1
            for v in dB.values():
                agg["dB"][int(v)] += 1
            agg["supB"] += len(dB)
            if agg["peorB"] is None and any(abs(int(v)) > 1 for v in dB.values()):
                agg["peorB"] = {"Lambda": [int(a) for a in Lam],
                                "d": {str(k): int(v) for k, v in dB.items()}}
        if not restD:
            agg["cierraD"] += 1
            for v in dD.values():
                agg["dD"][int(v)] += 1
            agg["supD"] += len(dD)
            if agg["peorD"] is None and any(abs(int(v)) > 1 for v in dD.values()):
                agg["peorD"] = {"Lambda": [int(a) for a in Lam],
                                "c": {str(k): int(v) for k, v in dD.items()}}
        if not restC:
            agg["cierraC"] += 1
            for v in dC.values():
                agg["dC"][int(v)] += 1
    print("")
    print("  t=%d  r=%d  (m'=%d, R'=%d)  Lambda_1 <= %d  ->  %d pesos   (no racionales: %d)"
          % (t, r, mp, Rp, cota, agg["n"], agg["no_racional"]))
    print("     C0  el pelado cierra:   B_r %d/%d   D_r %d/%d   C_r %d/%d  (senuelo)"
          % (agg["cierraB"], agg["n"], agg["cierraD"], agg["n"], agg["cierraC"], agg["n"]))
    print("     B1  reparto de d en base B_r : %s" % dict(sorted(agg["dB"].items())))
    print("     B2  reparto de c en base D_r : %s" % dict(sorted(agg["dD"].items())))
    print("     C1  reparto en base C_r      : %s   (senuelo)" % dict(sorted(agg["dC"].items())))
    print("     B3  soporte total:  B_r %d   D_r %d" % (agg["supB"], agg["supD"]))
    if agg["peorB"]:
        print("     !! primer |d| > 1 en base B_r : %s" % json.dumps(agg["peorB"]))
    else:
        print("     ** ningun |d| > 1 en base B_r")
    if agg["peorD"]:
        print("     !! primer |c| > 1 en base D_r : %s" % json.dumps(agg["peorD"]))
    else:
        print("     ** ningun |c| > 1 en base D_r")
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "cota": int(cota), "n": int(agg["n"]),
                "cierraB": int(agg["cierraB"]), "cierraD": int(agg["cierraD"]),
                "cierraC": int(agg["cierraC"]),
                "dB": {str(k): int(v) for k, v in agg["dB"].items()},
                "dD": {str(k): int(v) for k, v in agg["dD"].items()},
                "dC": {str(k): int(v) for k, v in agg["dC"].items()},
                "supB": int(agg["supB"]), "supD": int(agg["supD"]),
                "peorB": agg["peorB"], "peorD": agg["peorD"]})

print("")
print("=" * 118)
print("  LECTURA, escrita ANTES de correr:")
print("   * si en base B_r todos los coeficientes son 0,+-1 y en D_r no, entonces (L1) esta escrita")
print("     en la base derivada: el enunciado primitivo es sobre el CENTRALIZADOR, y (L1) se sigue")
print("     por la ramificacion B_r -> D_r, que es libre de multiplicidad.")
print("   * si en B_r aparecen coeficientes grandes, el centralizador NO es la explicacion y hay que")
print("     volver a la division por Delta_t.")
print("   * si el pelado en B_r no cierra, la serie no es combinacion de caracteres de B_r y la")
print("     hipotesis muere ahi mismo.")
json.dump(RES, open("gkrs_centralizer_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
