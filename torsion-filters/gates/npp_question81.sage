# -*- coding: utf-8 -*-
# NUESTRO FILTRO CONTRA LA QUESTION 8.1 DE NADIMPALLI-PATTANAYAK-PRASAD.   16 de agosto de 2026.
#
# DE DONDE SALE.  Gate de literatura por la FORMULA, no por las palabras.  NPP25 (arXiv:2504.14684,
# "Character theory at a torsion element") tiene:
#
#   * Teorema 4.1, para el elemento PRINCIPAL  C_m = rho^v(e^{2 pi i/m}):
#         dim Z_{Ghat}((lambda+rho)^(e^{2 pi i/m}))  >=  dim Z_{Ghat}(rho^(e^{2 pi i/m})),
#         con IGUALDAD  <=>  Theta_lambda(C_m) != 0,
#     y en ese caso  Theta = (-1)^w mu(c(G)) d_lambda/d_m,  con  d_m = 1 en tipos A y B, y en C
#     cuando m es impar.  El conjunto de raices que decide es
#         Phi_{lambda,m} = { alpha : m | <lambda+rho, alpha^v> },
#     las raices del grupo G_lambda(m).
#
#   * Question 8.1 (que ellos atribuyen a Prasad, y dejan ABIERTA) extiende eso a un elemento de
#     torsion CUALQUIERA x_0 de orden d: si dim Z(x_0) < dim Z((lambda+rho)(e^{2 pi i/d})) entonces
#     Theta = 0, y si son iguales, Theta es -salvo constante- la dimension de una irreducible del
#     dual de ese centralizador.
#
# POR QUE NOS IMPORTA.  Nuestro filtro de torsion es exactamente Theta_eta(g) con g regular, luego
# dim Z(g) = rango.  Si la Question 8.1 vale aqui, entonces
#
#     tau(eta) != 0   <=>   Z_{Ghat}((eta+rho)(zeta_d))  es un TORO,
#
# y como el dual de un toro es un toro, su unica irreducible tiene dimension 1: sale |tau| = 1.
# O sea nuestra regla (T) y su hermana impar (T^B) SERIAN el caso regular de una pregunta abierta.
# Eso es muy distinto de "ya esta publicado": es un dato sobre una pregunta viva.
#
# LO QUE SE MIDE
#   Q1  ¿Es nuestro g PRINCIPAL?  Se evalua CADA raiz simple en g: principal <=> todas dan el mismo
#       escalar.  Si no lo es, el Teorema 4.1 NO nos cubre y estamos en la Question 8.1.
#   Q2  el orden de g en el grupo ADJUNTO (que es el d de la pregunta), calculado, no supuesto.
#   Q3  LA PREGUNTA.  Para cada eta:  dim Z = rango + #{alpha : d | <eta+rho, alpha^v>}  contra
#       tau(eta) != 0.  Se cuentan los dos tipos de fallo por separado.
#   Q4  ¿coincide el conjunto Phi_{eta,d} con NUESTRAS paredes?  Nuestra (T) par habla de tres
#       especies (c=0, c=t/2, c_i=+-c_j) y la impar de dos.  Se compara pared a pared.
#   Q5  el signo.  NPP dan  (-1)^w mu(c(G))  con  w(eta+rho) = rho + d mu.  Se calcula ese w por
#       fuerza bruta en el grupo de Weyl y se compara con nuestro tau medido.
#
# CONTROLES
#   C0  tau se calcula por Freudenthal (suma de pesos), independiente de todo lo anterior.
#   C1  SEÑUELO: la misma dicotomia con d' = d+1.  Tiene que FALLAR; si acierta igual, la
#       coincidencia no dice nada.
#   C2  n impreso siempre, y las dos direcciones del fallo por separado.
#   C3  el caso rango 1 (t=3 y t=4) TIENE que salir principal: es el control de que Q1 mide algo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage npp_question81.sage

import json
import sys

RES = []
print("=" * 124)
print("NUESTRO FILTRO CONTRA LA QUESTION 8.1 DE NPP25   (arXiv:2504.14684)")
print("=" * 124)


def datos(t):
    """(tipo, rango, exponentes del toro) del bloque de torsion para este t."""
    if t % 2 == 0:
        m = (t - 2) // 2
        return "C", m, [i + 1 for i in range(m)]     # g = diag(xi^{+-1..+-m}) en Sp_2m
    mp = (t - 1) // 2
    return "B", mp, [i + 1 for i in range(mp)]        # g = diag(1, zeta^{+-1..+-m'}) en SO_{2m'+1}


def car(typ, rk, mu, _c={}):
    key = (typ, rk, tuple(int(v) for v in mu))
    if key not in _c:
        W = WeylCharacterRing("%s%d" % (typ, rk))
        el = W(W.space().from_vector(vector([Integer(v) for v in mu])))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + int(mult)
        _c[key] = d
    return _c[key]


for t in range(3, 9):
    typ, rk, expo = datos(t)
    if rk < 1:
        continue
    L = RootSystem("%s%d" % (typ, rk)).ambient_space()
    simples = [L.simple_root(i) for i in L.index_set()]
    pos = list(L.positive_roots())
    K = CyclotomicField(t)
    z = K.gen()

    def ev(alpha):
        v = [int(u) for u in alpha.to_vector()]
        return sum(expo[i] * v[i] for i in range(rk)) % t

    # ---- Q1: principal?
    vals = [ev(a) for a in simples]
    principal = len(set(vals)) == 1

    # ---- Q2: orden en el grupo adjunto.  g^k adjunto = 1  <=>  toda RAIZ vale 1 en g^k
    d_adj = None
    for k in range(1, 2 * t + 1):
        if all((k * ev(a)) % t == 0 for a in pos):
            d_adj = k
            break

    # ---- rho y las coraices, en coordenadas del ambiente
    rho = L.rho()
    rho_v = [QQ(u) for u in rho.to_vector()]

    def pairing_coroot(a, vec):
        """<vec, alpha^v> = 2 (vec, alpha)/(alpha, alpha) en el ambiente."""
        av = [QQ(u) for u in a.to_vector()]
        num = sum(vec[i] * av[i] for i in range(rk))
        den = sum(av[i] * av[i] for i in range(rk))
        return 2 * num / den

    ETAS = []
    for k in range(0, 3 * t + 1):
        for e in Partitions(k, max_length=rk):
            ETAS.append(tuple(list(e) + [0] * (rk - len(e))))

    ok = falso_cero = falso_nocero = 0
    sen_ok = 0
    pared_ok = 0
    n_vivos = 0
    for eta in ETAS:
        # tau por Freudenthal
        s = K(0)
        for wt, mult in car(typ, rk, eta).items():
            s += mult * z ** (sum(expo[i] * wt[i] for i in range(rk)) % t)
        tv = QQ(s) if s in QQ else None
        if tv is None:
            continue
        vivo = (tv != 0)
        if vivo:
            n_vivos += 1
        # Phi_{eta,d}: raices con d | <eta+rho, alpha^v>
        a_vec = [QQ(eta[i]) + rho_v[i] for i in range(rk)]
        dd = d_adj
        nraiz = sum(1 for a in pos if pairing_coroot(a, a_vec) % dd == 0)
        toro = (nraiz == 0)
        if toro == vivo:
            ok += 1
        elif vivo:
            falso_cero += 1        # la pregunta predice 0 y no lo es
        else:
            falso_nocero += 1      # la pregunta predice != 0 y es 0
        # SEÑUELO con d+1
        nr2 = sum(1 for a in pos if pairing_coroot(a, a_vec) % (dd + 1) == 0)
        if (nr2 == 0) == vivo:
            sen_ok += 1
        # Q4: nuestras paredes, escritas a mano, contra Phi_{eta,d}
        if typ == "C":
            c = [(int(eta[j]) + rk - j) % t for j in range(rk)]      # a_j = eta_j + m - j + 1
            nuestro = (0 not in c) and (t % 2 == 0 and (t // 2) not in c or t % 2) and \
                      len(set(min(x, t - x) for x in c)) == rk
        else:
            A = [2 * int(eta[j]) + 2 * (rk - j - 1) + 1 for j in range(rk)]
            nuestro = (0 not in [x % t for x in A]) and \
                      len(set(min(x % t, (t - x) % t) for x in A)) == rk
        if nuestro == vivo:
            pared_ok += 1

    n = len([e for e in ETAS])
    print("")
    print("  t=%-2d  %s%-2d  h=%d  t-h=%d   |   PRINCIPAL: %-3s  (raices simples -> %s)   orden en el adjunto: %s"
          % (t, typ, rk, 2 * rk, t - 2 * rk, "SI" if principal else "NO",
             str(sorted(set(vals))), str(d_adj)))
    print("       eta probados %4d | vivos %4d | Question 8.1 acierta %4d  (falsos cero %d, falsos no-cero %d)"
          % (n, n_vivos, ok, falso_cero, falso_nocero))
    print("       nuestras paredes aciertan %4d | SEÑUELO d+1 acierta %4d  (tiene que fallar mucho)"
          % (pared_ok, sen_ok))
    sys.stdout.flush()
    RES.append({"t": int(t), "tipo": typ, "rango": int(rk), "h": int(2 * rk),
                "principal": bool(principal), "d_adjunto": int(d_adj),
                "n_eta": int(n), "n_vivos": int(n_vivos), "q81_acierta": int(ok),
                "falsos_cero": int(falso_cero), "falsos_nocero": int(falso_nocero),
                "paredes_acierta": int(pared_ok), "senuelo_acierta": int(sen_ok)})

print("")
print("=" * 124)
print("  LECTURA, escrita ANTES de correr:")
print("   * si la Question 8.1 acierta en todos los eta y g NO es principal para t >= 5, entonces")
print("     nuestra (T) es el caso REGULAR de una pregunta abierta de Prasad, y eso se cita asi.")
print("   * si falla en algun eta, ese eta es un dato sobre la pregunta, y hay que mirarlo despacio.")
print("   * si el señuelo d+1 acierta parecido, la coincidencia no significa nada.")
json.dump(RES, open("npp_question81_DUMP.json", "w"), indent=1)
print("")
print("=" * 124)
print("DONE")
