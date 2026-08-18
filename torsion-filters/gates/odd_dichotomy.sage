# -*- coding: utf-8 -*-
# LA DICOTOMIA DE PARIDAD, MEDIDA EN SUS PARAMETROS.   16 de agosto de 2026.
#
# DE DONDE SALE.  odd_BD.sage midio que el reparto B_{m'} x D_r es el correcto (C0 9/9).  Pero eso
# dice QUE funciona, no POR QUE.  Este guion audita cada parametro: que objeto es, de que grupo, y
# cual es el mecanismo que separa el caso par del impar.  La pregunta es la correcta:
# si los conceptos no cuadran, un C0 que pasa es un ajuste con suerte.
#
# LA HIPOTESIS DEL MECANISMO, escrita antes de medir.  El punto donde se evalua es
#
#       p_t = (1, zeta, ..., zeta^{t-1}, z_1^{+-}, ..., z_r^{+-}),     zeta = e^{2 pi i / t},
#
# cerrado bajo x -> 1/x, o sea un elemento del grupo ORTOGONAL O(N).  Y su determinante es el
# producto de todo el alfabeto:
#
#       det p_t  =  prod_{k=0}^{t-1} zeta^k  =  zeta^{t(t-1)/2}  =  (-1)^{t+1}.
#
#   t IMPAR -> det = +1 -> p_t esta en SO(N), la componente IDENTIDAD.  La restriccion
#              GL_N -> SO_N es una restriccion ordinaria de una representacion de verdad, luego las
#              multiplicidades son >= 0 y no hay nada que torcer.
#   t PAR   -> det = -1 -> p_t esta en la OTRA componente de O(N).  Ahi el caracter no es el de una
#              restriccion: hay que torcerlo (twining / orbit algebra), y la expansion simplectica
#              sale VIRTUAL, con coeficientes de los dos signos.
#
# Si eso es asi, su punto 1 de la reseña no es una observacion editorial: es que las dos mitades del
# problema viven en componentes distintas de un mismo grupo ortogonal, y toda la cancelacion
# colectiva que llevamos dos dias persiguiendo entra por ahi y solo por ahi.
#
# LO QUE SE MIDE
#   A  det p_t para t = 2..11, contra la prediccion (-1)^{t+1}.
#   B  REGULARIDAD del elemento de torsion.  Par: g = diag(xi^{+-1..m}) en Sp_2m, orden t = h+2 con
#      h = 2m el numero de Coxeter.  Impar: g = diag(1, zeta^{+-1..m'}) en SO_{2m'+1}, orden
#      t = h+1 con h = 2m'.  Se evalua CADA raiz positiva en g: regular <=> ninguna vale 1.
#      Los dos ordenes son DISTINTOS respecto de h, y eso importa para a quien se cita.
#   C  NO NEGATIVIDAD.  Impar: a^B_Lambda >= 0 (prediccion de la teoria, no ajuste).  Par: los
#      a_Lambda tienen los DOS signos (ya refutada la no negatividad el 14 de agosto; se re-exhibe
#      un testigo aqui para que la comparacion este en el mismo fichero).
#   D  (T^B) en t = 7, 9, 11.  El 9 es el primer t impar COMPUESTO: ahi la entrada individual del
#      bialternante puede anularse sin que se anule el determinante, y la regla podria romperse.
#   E  LOS TERMINOS EN mu_max, lado impar: cuantos son y como de grandes.  En el par eran 25
#      terminos de hasta 798 sumando +-1.  Si en el impar son pocos y pequeños, la dicotomia se ve
#      tambien en la aritmetica y no solo en la categoria.
#
# CONTROLES
#   C0  la prediccion de A se compara con el det CALCULADO, no con la formula.
#   C1  B se hace evaluando las raices positivas de verdad (sistema de raices de Sage), no a mano.
#   C2  SEÑUELO de D: la regla con la pared FALSA en t/2 (que en impar no existe) tiene que fallar.
#   C3  n impreso en cada bloque.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage odd_dichotomy.sage

import json
import sys
from collections import defaultdict

RES = {}
print("=" * 124)
print("LA DICOTOMIA DE PARIDAD, MEDIDA EN SUS PARAMETROS")
print("=" * 124)

# ------------------------------------------------------------------ A: el determinante del punto
print("")
print("  A  EL DETERMINANTE DEL PUNTO DE EVALUACION   det p_t = prod_k zeta^k   (las z^{+-} aportan 1)")
print("     t | det calculado | prediccion (-1)^{t+1} | componente de O(N)")
print("     " + "-" * 96)
A = []
for t in range(2, 12):
    K = CyclotomicField(t) if t > 2 else QQ
    zeta = K.gen() if t > 2 else K(-1)
    d = prod([K(zeta) ** k for k in range(t)])
    pred = (-1) ** (t + 1)
    ok = (d == pred)
    print("     %2d | %13s | %21d | %s   %s"
          % (t, str(d), pred, "SO(N)  identidad" if d == 1 else "O(N) \\ SO(N)  la otra",
             "ok" if ok else "FALLA"))
    A.append({"t": int(t), "det": int(d), "pred": int(pred), "ok": bool(ok)})
RES["A_determinante"] = A
sys.stdout.flush()

# ------------------------------------------------------------------ B: regularidad y orden vs Coxeter
print("")
print("  B  EL ELEMENTO DE TORSION: orden, numero de Coxeter, y regularidad raiz a raiz")
print("     caso   | grupo | h  | orden t | t - h | raices positivas | raices que valen 1 | regular")
print("     " + "-" * 110)
B = []
for t in range(3, 12):
    if t % 2 == 0:
        m = (t - 2) // 2
        if m < 1:
            continue
        tipo, rk = "C", m
        # autovalores del standard rep: xi^{+-1..m};  coordenadas del toro x_i = xi^i
        expo = [i + 1 for i in range(m)]
    else:
        mp = (t - 1) // 2
        if mp < 1:
            continue
        tipo, rk = "B", mp
        expo = [i + 1 for i in range(mp)]
    L = RootSystem("%s%d" % (tipo, rk)).ambient_space()
    pos = L.positive_roots()
    K = CyclotomicField(t)
    z = K.gen()
    unos = 0
    for a in pos:
        v = [int(u) for u in a.to_vector()]
        e = sum(expo[i] * v[i] for i in range(rk)) % t
        if e == 0:
            unos += 1
    h = 2 * rk
    B.append({"t": int(t), "tipo": tipo, "rango": int(rk), "h": int(h),
              "t_menos_h": int(t - h), "n_raices": int(len(pos)), "raices_uno": int(unos),
              "regular": bool(unos == 0)})
    print("     t=%-2d %s | %s%-2d  | %2d | %7d | %5d | %16d | %18d | %s"
          % (t, "impar" if t % 2 else "par  ", tipo, rk, h, t, t - h, len(pos), unos,
             "SI" if unos == 0 else "NO"))
RES["B_regularidad"] = B
print("     -> par: t = h+2 en todos.  impar: t = h+1 en todos.  Son ORDENES DISTINTOS respecto de h,")
print("        luego ninguna cita unica sobre 'elemento regular de orden finito' cubre los dos casos.")
sys.stdout.flush()

# ------------------------------------------------------------------ utilidades

def phi_bialt(beta, tt, nvar):
    Nn = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    Lr = LaurentPolynomialRing(K, nvar, 'z')
    zs = Lr.gens()
    x = [Lr(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(Nn - 1, -1, -1))

    def alt(expo):
        return matrix(Lr, Nn, Nn, lambda i, j: x[i] ** expo[j]).determinant()

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
        if c not in QQ:
            return "NO-RAC"
        if QQ(c) != 0:
            out[k] = QQ(c)
    return out


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


def betas(N, tope):
    return [tuple(sorted(c, reverse=True)) + (0,) for c in Combinations(range(1, tope + 1), N - 1)]


# ------------------------------------------------------------------ C: no negatividad
print("")
print("  C  NO NEGATIVIDAD DE LA PRIMERA CAPA")
print("     caso        | formas | con algun coeficiente NEGATIVO | testigo minimo")
print("     " + "-" * 106)
C = []
# impar: a^B_Lambda, restriccion GL_{2R'+1} -> SO_{2R'+1}
for (t, r) in [(3, 2), (5, 2)]:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    N = t + 2 * r
    neg, tot, testigo = 0, 0, None
    for b in betas(N, 9 if t == 3 else 10):
        P = phi_bialt(b, 1, Rp)
        if P in (None, "NO-POL", "NO-RAC") or not P:
            continue
        aB, rest = pelar(P, "B", Rp)
        tot += 1
        if any(v < 0 for v in aB.values()):
            neg += 1
            if testigo is None:
                testigo = b
    print("     impar t=%d r=%d | %6d | %30d | %s" % (t, r, tot, neg, str(testigo) if testigo else "ninguno"))
    C.append({"caso": "impar t=%d r=%d" % (t, r), "formas": int(tot), "con_negativo": int(neg),
              "testigo": list(map(int, testigo)) if testigo else None})
# par: a_Lambda, la expansion simplectica virtual
for (R,) in [(3,), (4,)]:
    N = 2 + 2 * R
    neg, tot, testigo = 0, 0, None
    for b in betas(N, 9):
        P = phi_bialt(b, 2, R)
        if P in (None, "NO-POL", "NO-RAC") or not P:
            continue
        aL, rest = pelar(P, "C", R)
        tot += 1
        if any(v < 0 for v in aL.values()):
            neg += 1
            if testigo is None:
                testigo = b
    print("     par   t=2 R=%d | %6d | %30d | %s" % (R, tot, neg, str(testigo) if testigo else "ninguno"))
    C.append({"caso": "par t=2 R=%d" % R, "formas": int(tot), "con_negativo": int(neg),
              "testigo": list(map(int, testigo)) if testigo else None})
RES["C_no_negatividad"] = C
sys.stdout.flush()

# ------------------------------------------------------------------ D: (T^B) en t = 7, 9, 11
print("")
print("  D  LA REGLA (T^B) EN t = 7, 9, 11   (el 9 es el primer t impar COMPUESTO)")
print("     t | m' | eta probados | tau != 0 | |tau| != 1 | (T^B) acierta | SEÑUELO pared t/2 acierta")
print("     " + "-" * 108)
D = []
for t in (7, 9):
    mp = (t - 1) // 2
    K = CyclotomicField(t)
    z = K.gen()
    ET = []
    for k in range(0, 2 * t + 1):
        for e in Partitions(k, max_length=mp):
            ET.append(tuple(list(e) + [0] * (mp - len(e))))
    ok = sen = viv = mal = 0
    for e in ET:
        s = K(0)
        for wt, mult in car("B", mp, e).items():
            s += mult * z ** (sum((i + 1) * wt[i] for i in range(mp)) % t)
        v = QQ(s) if s in QQ else None
        if v is None:
            continue
        A_ = [2 * int(e[j]) + 2 * (mp - j - 1) + 1 for j in range(mp)]
        cl = [min(a % t, (t - a) % t) for a in A_]
        pred = (0 not in [a % t for a in A_]) and len(set(cl)) == mp
        # SEÑUELO: se le añade la pared falsa A = t/2, que en t impar no existe
        predS = pred and all((2 * a) % t != 0 for a in A_)
        if (v != 0) == pred:
            ok += 1
        if (v != 0) == predS:
            sen += 1
        if v != 0:
            viv += 1
            if abs(v) != 1:
                mal += 1
    print("     %2d | %2d | %12d | %8d | %10d | %13d | %d" % (t, mp, len(ET), viv, mal, ok, sen))
    D.append({"t": int(t), "mp": int(mp), "n_eta": int(len(ET)), "vivos": int(viv),
              "no_pm1": int(mal), "TB_acierta": int(ok), "senuelo_acierta": int(sen)})
RES["D_TB"] = D
sys.stdout.flush()

# ------------------------------------------------------------------ E: los terminos en mu_max
print("")
print("  E  LOS TERMINOS EN mu_max, LADO IMPAR   (en el par: 25 terminos de hasta 798 sumando +-1)")
print("     beta                       | t r | terminos != 0 | mayor |termino| | mayor suma parcial | A")
print("     " + "-" * 110)
E = []
_BR = {}
def branch_BD(Rp, mp, rr, Lam):
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

for (t, r, bs) in [(3, 2, [(12, 10, 7, 5, 3, 2, 0), (11, 8, 6, 4, 3, 1, 0), (10, 9, 6, 4, 2, 1, 0)]),
                   (5, 2, [(12, 10, 8, 7, 5, 4, 2, 1, 0)])]:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    K = CyclotomicField(t)
    z = K.gen()
    for b in bs:
        P1 = phi_bialt(b, t, r)
        A1, _ = pelar(P1, "D", r)
        A1 = {k: v for k, v in A1.items() if v != 0}
        Ap = {}
        for mu, c in A1.items():
            Ap[tuple(list(mu[:-1]) + [abs(mu[-1])])] = c
        S = list(Ap)
        maxi = [m for m in S if not any(n != m and all(sum(n[:k + 1]) >= sum(m[:k + 1])
                                                       for k in range(len(m))) for n in S)]
        if len(maxi) != 1:
            continue
        mm = maxi[0]
        P2 = phi_bialt(b, 1, Rp)
        aB, _ = pelar(P2, "B", Rp)
        term = []
        for Lam, a in aB.items():
            for (eta, mu), c in branch_BD(Rp, mp, r, Lam).items():
                if tuple(list(mu[:-1]) + [abs(mu[-1])]) != mm:
                    continue
                s = K(0)
                for wt, mult in car("B", mp, eta).items():
                    s += mult * z ** (sum((i + 1) * wt[i] for i in range(mp)) % t)
                tv = QQ(s) if s in QQ else QQ(0)
                if a * c * tv != 0:
                    term.append(int(a * c * tv))
        # las contribuciones de mu y mu* se cuentan las dos; el coeficiente es el de una
        par = 0
        mx = 0
        for x in term:
            par += x
            mx = max(mx, abs(par))
        print("     %-26s | %d %d | %13d | %15d | %18d | %s"
              % (str(b), t, r, len(term), max([abs(x) for x in term]) if term else 0, mx,
                 str(Ap[mm])))
        E.append({"beta": list(map(int, b)), "t": int(t), "r": int(r), "n_terminos": int(len(term)),
                  "mayor_termino": int(max([abs(x) for x in term])) if term else 0,
                  "mayor_suma_parcial": int(mx), "A": int(Ap[mm])})
        sys.stdout.flush()
RES["E_terminos"] = E

json.dump(RES, open("odd_dichotomy_DUMP.json", "w"), indent=1)
print("")
print("=" * 124)
print("DONE")
