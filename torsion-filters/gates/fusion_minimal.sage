# -*- coding: utf-8 -*-
# LOS FILTROS COMO COCIENTES DE FUSION DE NIVEL MINIMO.   16 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 24.  Andersen-Stroppel presentan el anillo de fusion de nivel k como
#
#     Z[chi(omega_1), ..., chi(omega_n)] / < chi(k omega_1 + omega_i) >,
#
# y con nuestros parametros el nivel sale MINIMO en las dos paridades:
#
#     t par  = 2m+2,  tipo C_m:  alcove = { sum m_i omega_i : sum m_i < l/2 - n } con l/2-n = 1,
#              o sea SOLO el peso 0, y nivel k = 0 -> el ideal lo generan los chi(omega_i).
#     t impar = 2m'+1, tipo B_m': alcove  2m_1+...+2m_{n-1}+m_n <= l-2n = 1, o sea 0 y omega_n;
#              pero omega_n es el peso ESPINORIAL y no esta en el retículo TENSORIAL de SO_{2m'+1},
#              luego en nuestro sector queda solo el 0.  Nivel k = 1 -> generadores chi(omega_i),
#              i < n,  y chi(2 omega_n).
#
# Si los generadores del ideal se anulan en nuestro elemento de torsion, la evaluacion FACTORIZA por
# el cociente de fusion minimo, que como anillo es Z.  Y eso explica de golpe el {0,+-1}, el
# "singular muere / regular pliega con signo", y que todos los supervivientes caigan en un punto.
#
# LO QUE SE MIDE
#   F1  los caracteres fundamentales en el elemento:  chi(omega_i) para todo i, las dos paridades.
#       Predicho: CERO todos.
#   F2  el mecanismo:  prod_{alpha en el espectro} (1 + alpha u)  tiene que valer  1 - (-1)^t u^t,
#       o sea todos los e_i intermedios se anulan.  Se calcula el polinomio entero.
#   F3  en el impar, el peso espinorial: se comprueba que omega_n NO es un peso de ninguna
#       representacion TENSORIAL, evaluando la integralidad en el retículo de SO.
#   F4  el generador extra del nivel 1 en tipo B: chi(2 omega_n) = el caracter de Lambda^n(V).
#
# CONTROLES
#   C0  los caracteres se calculan por Freudenthal, sin usar e_i ni el producto.
#   C1  SEÑUELO: los mismos caracteres fundamentales en un elemento de orden t+1.  NO tienen que
#       anularse todos, o el fenomeno no seria del elemento.
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage fusion_minimal.sage

import json

print("=" * 112)
print("LOS GENERADORES DEL IDEAL DE FUSION, EVALUADOS EN NUESTRO ELEMENTO DE TORSION")
print("=" * 112)


def car_eval(typ, rk, mu, t, salto=1):
    """el caracter irreducible de peso mu evaluado en el toro x_i = xi^{salto*i}."""
    K = CyclotomicField(t)
    z = K.gen()
    W = WeylCharacterRing("%s%d" % (typ, rk))
    el = W(W.space().from_vector(vector([QQ(v) for v in mu])))
    s = K(0)
    for wt, mult in el.weight_multiplicities().items():
        v = [QQ(u) for u in wt.to_vector()]
        e = sum(salto * (i + 1) * v[i] for i in range(rk))
        if e not in ZZ:
            return "NO-ENTERO"
        s += int(mult) * z ** (int(e) % t)
    return QQ(s) if s in QQ else s


RES = []
for t in range(3, 12):
    if t % 2 == 0:
        typ, rk = "C", (t - 2) // 2
        nivel = t // 2 - rk - 1
    else:
        typ, rk = "B", (t - 1) // 2
        nivel = t - 2 * rk
    if rk < 1:
        continue
    L = RootSystem("%s%d" % (typ, rk)).ambient_space()

    # los pesos fundamentales, en coordenadas del ambiente
    fund = []
    for i in L.index_set():
        w = L.fundamental_weight(i)
        fund.append((i, [QQ(v) for v in w.to_vector()]))

    # ESPINORIAL o no lo decide el RETICULO, no la evaluacion: omega_n de B_n es (1/2,...,1/2) y
    # simplemente NO esta en el retículo de pesos de SO_{2n+1}.  La primera version de este guion lo
    # decidia mirando si el exponente caia en Z, que es otra cosa y daba 'espinorial' unas veces si y
    # otras no para el mismo peso.
    vals, senu, espin = [], [], []
    for (i, w) in fund:
        es_espin = any(QQ(v).denominator() != 1 for v in w)
        if es_espin:
            espin.append(int(i))
            continue
        v = car_eval(typ, rk, w, t)
        vals.append((i, v))
        senu.append((i, car_eval(typ, rk, w, t + 1)))
    # generador extra del nivel 1 en tipo B: chi(2 omega_n) = Lambda^n(V), que SI es tensorial
    extra = None
    if typ == "B":
        wn = [2 * QQ(v) for v in fund[-1][1]]
        extra = car_eval(typ, rk, wn, t)
    ceros = sum(1 for _, v in vals if v == 0)
    ceros_s = sum(1 for _, s in senu if s == 0)
    print("")
    print("  t=%2d  %s%-2d  nivel k = %d" % (t, typ, rk, nivel))
    print("     chi(omega_i) TENSORIALES en el elemento : %s"
          % str([(int(i), str(v)) for i, v in vals]))
    print("     -> cero en %d de %d   |  pesos espinoriales excluidos del retículo de SO: %s"
          % (ceros, len(vals), str(espin) if espin else "ninguno"))
    if extra is not None:
        print("     generador extra del nivel 1,  chi(2 omega_%d) = Lambda^%d(V) : %s"
              % (rk, rk, str(extra)))
    print("     SEÑUELO, los mismos en un elemento de orden t+1: cero en %d de %d"
          % (ceros_s, len(senu)))
    RES.append({"t": int(t), "tipo": typ, "rango": int(rk), "nivel": int(nivel),
                "n_tensoriales": int(len(vals)), "ceros": int(ceros),
                "indices_espinoriales": espin,
                "chi_2omega_n": str(extra) if extra is not None else None,
                "senuelo_ceros": int(ceros_s)})

print("")
print("=" * 112)
print("F2  EL MECANISMO:  prod_{alpha en el espectro} (1 + alpha u)")
print("=" * 112)
for t in range(3, 12):
    K = CyclotomicField(t)
    z = K.gen()
    P = PolynomialRing(K, 'u')
    u = P.gen()
    if t % 2 == 0:
        m = (t - 2) // 2
        espectro = [z ** (j + 1) for j in range(m)] + [z ** (-(j + 1)) for j in range(m)]
        etiqueta = "C_%d: xi^{+-1..%d}" % (m, m)
    else:
        mp = (t - 1) // 2
        espectro = [K(1)] + [z ** (j + 1) for j in range(mp)] + [z ** (-(j + 1)) for j in range(mp)]
        etiqueta = "B_%d: 1, xi^{+-1..%d}" % (mp, mp)
    pol = prod([1 + a * u for a in espectro])
    print("  t=%2d  %-22s ->  %s" % (t, etiqueta, str(pol)))

print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr: si TODOS los chi(omega_i) tensoriales se anulan y el")
print("  señuelo de orden t+1 no, la evaluacion factoriza por el cociente de fusion minimo y el")
print("  {0,+-1} deja de ser un calculo para ser una consecuencia.")
json.dump(RES, open("fusion_minimal_DUMP.json", "w"), indent=1)
print("=" * 112)
print("DONE")
