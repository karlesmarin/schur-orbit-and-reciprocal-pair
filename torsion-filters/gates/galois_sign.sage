# -*- coding: utf-8 -*-
# EL SIGNO DE GALOIS DEL FILTRO.   16 de agosto de 2026.   (vuelta 28, punto 6)
#
# DE DONDE SALE.  La Observacion 3.9 dice ahora "una evaluacion por orbita y el filtro queda conocido
# en toda la orbita".  Un testigo de dos lineas la mata:  t=4, m=1,  tau(a=1)=+1
# pero tau(a=3)=-1,  y 3 es una unidad.  Lo INVARIANTE en la orbita es el SOPORTE; el valor no.
#
# Y propone el arreglo, que es mejor que retirar la frase.  Multiplicar los residuos por una unidad
# k induce sobre las CLASES PLEGADAS una permutacion con signo  w_k;  poniendo
#
#        gamma_t(k) := det(w_k)  en  {+-1},
#
# el mismo argumento de columnas del bialternante da     tau_t(k.a) = gamma_t(k) . tau_t(a),
# y gamma_t es un caracter  (Z/t)^x -> {+-1}.  Eso recupera "una evaluacion por orbita" CON los signos.
#
# NOTA.  Para nosotros esto no es una conjetura: nuestro Lema 3.1 ya da  tau = sgn(sigma) prod eps.
# Componer con w_k es entonces inmediato.  Lo que aqui se hace es MEDIRLO con un tau independiente
# (evaluacion del caracter, no la formula cerrada), que es la unica forma de que el control valga.
#
# LO QUE SE MIDE
#   G1  gamma_t es homomorfismo:  gamma(k1 k2) = gamma(k1) gamma(k2)  para TODOS los pares.  FATAL.
#   G2  tau_t(k.a) = gamma_t(k) tau_t(a)  sobre una caja de eta y TODAS las unidades, con tau
#       calculado evaluando el caracter (independiente de la formula cerrada).  FATAL.
#   G3  .quien es gamma_t?  se compara con el simbolo de Jacobi (k/t), con el trivial, y con el
#       caracter  k -> (-1)^{(k-1)/2}.  (Zolotarev vive aqui cerca.)
#   G4  los DOS testigos que fijan el enunciado.
#   G5  R^B = {x_i != 0, x_i != +-x_j}  contra  R^C = {2x_i != 0, x_i != +-x_j}:  ambos estables
#       bajo unidades, e iguales si y solo si t es impar.
#
# CONTROLES
#   C0  G1 y G2 son fatales.
#   C1  SENUELO: se prueba tambien con k NO unidad; tiene que fallar.
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage galois_sign.sage

import json
import sys
import itertools
from collections import Counter

IMPARES = [3, 5, 7, 9, 11]
PARES = [4, 6, 8, 10, 12]


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


def tau_eval(typ, rk, eta, t):
    """tau evaluando el CARACTER en el punto de torsion.  Independiente de la formula cerrada."""
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car(typ, rk, eta).items():
        s += mult * z ** (sum((i + 1) * int(wt[i]) for i in range(rk)) % t)
    return int(QQ(s)) if s in QQ else None


def plegar(v, t):
    """clase plegada y signo:  v -> (min(v, t-v), +-1).  Devuelve (clase, eps)."""
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)          # autopareada: la pared de mas del caso par
    return (v, 1) if v < t - v else (t - v, -1)


def w_de_unidad(k, t, n):
    """la permutacion con signo que k induce sobre las clases plegadas {1..n}, y su det.
       n = (t-1)/2 si t impar;  n = t/2 - 1 si t par  (la clase t/2 es fija y se excluye)."""
    cl, eps = [], []
    for j in range(1, n + 1):
        c, e = plegar(k * j, t)
        cl.append(c)
        eps.append(e)
    if sorted(cl) != list(range(1, n + 1)):
        return None
    # det = sgn(perm) . prod eps
    perm = [cl[i] - 1 for i in range(n)]
    s = 1
    visto = [False] * n
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
    for e in eps:
        s *= e
    return int(s)


def delta_de(a, rk, t):
    """det de la permutacion con signo que las clases plegadas de a definen; 0 si no es permutacion."""
    cl, ep = [], []
    for v in a:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, rk + 1)):
        return 0
    perm = [cl[i] - 1 for i in range(rk)]
    sg = 1
    visto = [False] * rk
    for i in range(rk):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            sg = -sg
    for e in ep:
        sg *= e
    return int(sg)


def residuos_a_eta(res, typ, rk, t):
    """dado un vector de residuos mod t, devuelve (eta dominante con esos residuos, signo de la
       permutacion que hubo que aplicar), o None si dos residuos coinciden.  Usa SOLO la
       periodicidad (tau depende del residuo) y que tau es alternante; no la formula cerrada."""
    # tau solo depende del residuo (periodicidad), luego se puede subir todo por un multiplo comun
    # de t.  Se usa el representante MINIMO y se sube solo si hace falta: un eta grande dispara
    # Freudenthal y la corrida se cuelga en B_4 / B_5.
    if typ == "C":
        rep = [int(v) % t for v in res]
        paso = t
    else:
        rep = [(int(v) % t) if (int(v) % t) % 2 == 1 else (int(v) % t) + t for v in res]
        paso = 2 * t
    if len(set(rep)) != rk:
        return None
    for _ in range(4):
        ap0 = sorted(rep, reverse=True)
        if typ == "C":
            e0 = [ap0[j] - (rk - j) for j in range(rk)]
        else:
            e0 = [(ap0[j] - 2 * (rk - j) + 1) // 2 for j in range(rk)]
        if min(e0) >= 0 and all(e0[j] >= e0[j + 1] for j in range(rk - 1)):
            break
        rep = [v + paso for v in rep]
    idx = sorted(range(rk), key=lambda i: -rep[i])
    s = 1
    visto = [False] * rk
    for i in range(rk):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = idx[j]
            L += 1
        if L % 2 == 0:
            s = -s
    ap = [rep[i] for i in idx]
    # rho_{C_n} = (n, ..., 1)  ->  a_j = eta_j + n - j   (j de 0 a n-1)
    # rho_{B_n} = (n-1/2, ..., 1/2)  ->  A_j = 2 eta_j + 2(n-j) - 1
    if typ == "C":
        eta = [ap[j] - (rk - j) for j in range(rk)]
    else:
        eta = [(ap[j] - 2 * (rk - j) + 1) // 2 for j in range(rk)]
    if min(eta) < 0 or any(eta[j] < eta[j + 1] for j in range(rk - 1)):
        return None
    return (tuple(eta), int(s))


def etas(rk, cota):
    def rec(k, tope):
        if k == 0:
            yield ()
            return
        for a in range(tope, -1, -1):
            for resto in rec(k - 1, a):
                yield (a,) + resto
    return list(rec(rk, cota))


print("=" * 118)
print("EL SIGNO DE GALOIS DEL FILTRO:  tau_t(k.a) = gamma_t(k) tau_t(a)")
print("=" * 118)
sys.stdout.flush()

RES = []
for (t, typ) in [(x, "B") for x in IMPARES] + [(x, "C") for x in PARES]:
    if typ == "B":
        rk = (t - 1) // 2                  # B_{m'},  t = 2m'+1
        nfold = (t - 1) // 2
    else:
        rk = (t - 2) // 2                  # C_m,     t = 2m+2
        nfold = t // 2 - 1
    if rk < 1:
        continue
    unidades = [k for k in range(1, t) if gcd(k, t) == 1]
    gamma = {}
    for k in unidades:
        gamma[k] = w_de_unidad(k, t, nfold)

    # G1  homomorfismo
    g1_ok = g1_n = 0
    for k1 in unidades:
        for k2 in unidades:
            if gamma[k1] is None or gamma[k2] is None:
                continue
            k3 = (k1 * k2) % t
            if gamma.get(k3) is None:
                continue
            g1_n += 1
            g1_ok += 1 if gamma[k3] == gamma[k1] * gamma[k2] else 0

    # G2  la covarianza, con tau independiente
    # C3  CONTROL DEL INSTRUMENTO, primero: la formula cerrada  tau = sgn(sigma) prod eps  tiene que
    #     coincidir con la evaluacion del caracter.  Si esto falla, los vectores desplazados estan
    #     mal y todo lo de abajo mide otra cosa.
    #     OJO: en tipo B la formula cerrada solo vale SALVO UNA CONSTANTE  epsilon_t = tau(rho)/delta(rho),
    #     que el paper no enuncia (prop:oddfilter solo da |tau|=1).  Se calibra con eta = 0 y se
    #     IMPRIME, porque esa constante es un dato: vale +1 en t=3,5 y -1 en t=7.
    cota = 3 if rk <= 3 else 2
    c3_ok = c3_n = 0
    c3_fallo = None
    eps_t = None
    # PRIMERA PASADA: calibrar eps_t.  Puntuar en la misma pasada haria fallar los pesos que van
    # ANTES de la calibracion, que es lo que paso en la primera version de este guion.
    for eta in etas(rk, cota):
        tv0 = tau_eval(typ, rk, eta, t)
        if not tv0:
            continue
        if typ == "C":
            aa = [int(eta[j]) + rk - j for j in range(rk)]
        else:
            aa = [2 * int(eta[j]) + 2 * (rk - j) - 1 for j in range(rk)]
        d0 = delta_de(aa, rk, t)
        if d0:
            eps_t = int(tv0) * int(d0)
            break

    for eta in etas(rk, cota):
        tv = tau_eval(typ, rk, eta, t)
        if tv is None:
            continue
        if typ == "C":
            a0 = [int(eta[j]) + rk - j for j in range(rk)]
        else:
            a0 = [2 * int(eta[j]) + 2 * (rk - j) - 1 for j in range(rk)]
        pred = delta_de(a0, rk, t)
        c3_n += 1
        if eps_t is not None and pred * eps_t == tv:
            c3_ok += 1
        elif c3_fallo is None:
            c3_fallo = {"eta": [int(v) for v in eta], "a": [int(v) for v in a0],
                        "delta": int(pred), "eps_t": (int(eps_t) if eps_t is not None else None),
                        "caracter": int(tv)}

    # el caso  tau = 0  cumple la identidad trivialmente, asi que se puntua APARTE el subconjunto
    # con  tau != 0, que es el unico donde la afirmacion tiene contenido.
    g2_ok = g2_n = g2_vivos = g2_ok_vivos = g2_skip = 0
    fallo = None
    for eta in etas(rk, cota):
        tv = tau_eval(typ, rk, eta, t)
        if tv is None:
            continue
        # residuos desplazados:  tipo C  a_j = eta_j + rk - j ;  tipo B  A_j = 2eta_j + 2(rk-j) - 1
        if typ == "C":
            a = [int(eta[j]) + rk - j for j in range(rk)]
        else:
            a = [2 * int(eta[j]) + 2 * (rk - j) - 1 for j in range(rk)]
        for k in unidades:
            gk = gamma[k]
            if gk is None:
                continue
            ka = [(k * v) % t for v in a]
            r2 = residuos_a_eta(ka, typ, rk, t)
            if r2 is None:
                g2_skip += 1
                continue
            eta2, sperm = r2
            tv2 = tau_eval(typ, rk, eta2, t)
            if tv2 is None:
                g2_skip += 1
                continue
            # tau es alternante en los residuos:  T(k.a) = sperm . tau(eta2)
            Tka = sperm * tv2
            g2_n += 1
            if tv != 0:
                g2_vivos += 1
            if Tka == gk * tv:
                g2_ok += 1
                if tv != 0:
                    g2_ok_vivos += 1
            elif fallo is None:
                fallo = {"eta": [int(v) for v in eta], "k": int(k), "gamma": int(gk),
                         "tau": int(tv), "eta_k": [int(v) for v in eta2],
                         "sperm": int(sperm), "T_ka": int(Tka)}

    # G3  .quien es gamma?
    # gamma = sgn(sigma_k) . prod eps_k ;  por el lema de Gauss  prod eps_k = simbolo de Jacobi (k/t).
    # Se separan los dos factores y se identifica cada uno.
    sgn_s, prod_e = {}, {}
    for k in unidades:
        cl, ep = [], []
        for j in range(1, nfold + 1):
            cc_, ee_ = plegar(k * j, t)
            cl.append(cc_)
            ep.append(ee_)
        if sorted(cl) != list(range(1, nfold + 1)):
            sgn_s[k] = prod_e[k] = None
            continue
        perm = [cl[i] - 1 for i in range(nfold)]
        sg = 1
        visto = [False] * nfold
        for i in range(nfold):
            if visto[i]:
                continue
            j, L = i, 0
            while not visto[j]:
                visto[j] = True
                j = perm[j]
                L += 1
            if L % 2 == 0:
                sg = -sg
        sgn_s[k] = int(sg)
        pe = 1
        for e in ep:
            pe *= e
        prod_e[k] = int(pe)
    ident = []
    if all(gamma[k] is not None for k in unidades):
        if all(gamma[k] == 1 for k in unidades):
            ident.append("trivial")
        if t % 2 == 1 and all(gamma[k] == jacobi_symbol(k, t) for k in unidades):
            ident.append("Jacobi (k/t)")
        if t % 2 == 1 and all(gamma[k] == -jacobi_symbol(k, t) for k in unidades):
            ident.append("-Jacobi (k/t)")
        if all(gamma[k] == (-1) ** ((k - 1) // 2) for k in unidades):
            ident.append("(-1)^{(k-1)/2}")
    gauss = None
    if t % 2 == 1 and all(prod_e[k] is not None for k in unidades):
        gauss = all(prod_e[k] == jacobi_symbol(k, t) for k in unidades)
    print("")
    print("  t=%2d  tipo %s  rango %d  (clases plegadas 1..%d)   unidades: %s"
          % (t, typ, rk, nfold, unidades))
    print("     gamma_t : %s" % {int(k): (int(v) if v is not None else None) for k, v in gamma.items()})
    print("     C3  INSTRUMENTO: tau == eps_t . delta        : %d de %d   (eps_t = %s)"
          % (c3_ok, c3_n, eps_t))
    if c3_fallo:
        print("        !! el instrumento esta mal, lo de abajo NO vale: %s" % json.dumps(c3_fallo))
    print("     G1  homomorfismo (Z/t)^x -> {+-1}          : %d de %d" % (g1_ok, g1_n))
    print("     G2  tau(k.a) = gamma(k) tau(a)             : %d de %d   |  con tau != 0: %d de %d"
          " (los saltados por el guardia: %d)" % (g2_ok, g2_n, g2_ok_vivos, g2_vivos, g2_skip))
    print("     G3  gamma_t se identifica con              : %s" % (ident if ident else "ninguno de los probados"))
    print("         gamma = sgn(sigma_k) . prod eps_k      : sgn=%s  prod_eps=%s"
          % ({int(k): sgn_s[k] for k in unidades}, {int(k): prod_e[k] for k in unidades}))
    print("         prod eps_k == Jacobi (k/t)  (lema de Gauss) : %s" % gauss)
    if fallo:
        print("     !! primer fallo de G2: %s" % json.dumps(fallo))
    sys.stdout.flush()
    RES.append({"t": int(t), "tipo": typ, "rango": int(rk), "nfold": int(nfold),
                "gamma": {str(k): (int(v) if v is not None else None) for k, v in gamma.items()},
                "G1": [int(g1_ok), int(g1_n)],
                "G2": [int(g2_ok), int(g2_n), int(g2_ok_vivos), int(g2_vivos), int(g2_skip)],
                "C3": [int(c3_ok), int(c3_n)], "C3_fallo": c3_fallo,
                "eps_t": (int(eps_t) if eps_t is not None else None),
                "sgn_sigma": {str(k): sgn_s[k] for k in unidades},
                "prod_eps": {str(k): prod_e[k] for k in unidades},
                "gauss": (bool(gauss) if gauss is not None else None),
                "identificacion": ident, "fallo": fallo})


# --------------------------------------------------------------- G4  los dos testigos de la vuelta 28
print("")
print("=" * 118)
print("G4  LOS DOS TESTIGOS QUE FIJAN EL ENUNCIADO")
print("=" * 118)
G4 = {}
# (i)  t=4, m=1:  tau(a=1) = +1  pero  tau(a=3) = -1,  y 3 es unidad
#      a = eta_1 + 1  ->  a=1 es eta=(0);  a=3 es eta=(2)
w1 = tau_eval("C", 1, (0,), 4)
w2 = tau_eval("C", 1, (2,), 4)
print("  (i)  t=4, C_1 :  tau(a=1) = %s   tau(a=3) = %s   -> el SOPORTE es invariante, el VALOR no"
      % (w1, w2))
G4["t4_C1"] = {"tau_a1": w1, "tau_a3": w2, "gamma_4_de_3": w_de_unidad(3, 4, 1)}
# (ii) t=5, n=2, eta=(0,0):  tau^B = +1  pero la traslacion (t-1)/2 . 1 = (2,2) da tau~^C = -1
b = tau_eval("B", 2, (0, 0), 5)
cc = tau_eval("C", 2, (2, 2), 5)
print("  (ii) t=5      :  tau^B_5(0,0) = %s   tau~^C_5(2,2) = %s   -> los LOCUS coinciden, los valores no"
      % (b, cc))
G4["t5_traslacion"] = {"tauB_00": b, "tauC_22": cc}

# --------------------------------------------------------------- G5  R^B contra R^C
print("")
print("=" * 118)
print("G5  R^B = {x_i != 0, x_i != +-x_j}   contra   R^C = {2x_i != 0, x_i != +-x_j}")
print("=" * 118)
G5 = []
for t in [3, 4, 5, 6, 7, 8, 9]:
    for n in [1, 2, 3]:
        if n >= t:
            continue
        pts = list(itertools.product(range(t), repeat=n))
        def regB(x):
            return all(v % t != 0 for v in x) and all((x[i] - x[j]) % t != 0 and (x[i] + x[j]) % t != 0
                                                      for i in range(n) for j in range(i + 1, n))
        def regC(x):
            return all((2 * v) % t != 0 for v in x) and all((x[i] - x[j]) % t != 0 and (x[i] + x[j]) % t != 0
                                                            for i in range(n) for j in range(i + 1, n))
        RB = set(x for x in pts if regB(x))
        RC = set(x for x in pts if regC(x))
        uds = [k for k in range(1, t) if gcd(k, t) == 1]
        noud = [k for k in range(1, t) if gcd(k, t) != 1]
        estB = all(set(tuple((k * v) % t for v in x) for x in RB) == RB for k in uds)
        estC = all(set(tuple((k * v) % t for v in x) for x in RC) == RC for k in uds)
        # SENUELO: las no unidades no pueden preservarlo (salvo que RB sea vacio)
        malB = sum(1 for k in noud if RB and set(tuple((k * v) % t for v in x) for x in RB) == RB)
        print("  t=%d n=%d : |R^B|=%3d |R^C|=%3d  iguales=%-5s  unidades estables B/C = %s/%s  "
              "  senuelo no-unidades que preservan R^B: %d de %d"
              % (t, n, len(RB), len(RC), str(RB == RC), estB, estC, malB, len(noud)))
        G5.append({"t": int(t), "n": int(n), "RB": len(RB), "RC": len(RC), "iguales": bool(RB == RC),
                   "estableB": bool(estB), "estableC": bool(estC),
                   "senuelo_no_unidades": [int(malB), int(len(noud))]})
print("")
print("  LECTURA: R^B = R^C  debe darse EXACTAMENTE cuando t es impar (2 invertible mod t).")

json.dump({"por_t": RES, "G4": G4, "G5": G5}, open("galois_sign_DUMP.json", "w"), indent=1)
print("")
print("=" * 118)
print("DONE")
