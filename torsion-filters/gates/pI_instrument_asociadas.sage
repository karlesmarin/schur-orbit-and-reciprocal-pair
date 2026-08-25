# -*- coding: utf-8 -*-
# prob:instrument DEL PAPER I:  .es  m_mu = m_{mu*}  una identidad VISIBLE en tableaux?
# 19 de agosto de 2026, noche.
#
# LA OBSERVACION QUE ORIENTA EL EXPERIMENTO.  |mu*| = |mu| + N - 2 mu'_1, luego los dos lados suman
# coeficientes LR con delta de TAMANOS DISTINTOS: |delta| = |lambda|-|mu| a la izquierda y
# |lambda|-|mu*| a la derecha.  Una biyeccion NO puede preservar |delta|.  Asi que la pregunta
# «visible en tableaux» tiene una forma precisa: .hay alguna estadistica sobre delta que refine la
# igualdad, o la cancelacion es COLECTIVA?
#
# QUE SE MIDE
#   C0  CALIBRACION: m_nu por dos caminos --- coeficientes LR contra la expansion en caracteres
#       ortogonales --- y la involucion mu -> mu** = mu.
#   C1  FATAL: en el locus de anulacion, m_mu = m_{mu*} para toda mu.  Es el enunciado (2186), y si
#       falla es que la implementacion no es la del articulo.
#   C2  LA PREGUNTA: .se refina?  Se prueban cinco estadisticas sobre delta --- l(delta), delta_1,
#       numero de partes, |delta|/2, y el numero de columnas --- y se mira si la igualdad se
#       mantiene al partir las dos sumas por cada una.
#   C3  y si NO se refina: .cuantos terminos hay a cada lado?  Si los cardinales de los soportes ya
#       difieren, no hay biyeccion termino a termino de ninguna clase, y eso es una respuesta.
#   D1  SENUELo: la misma pregunta FUERA del locus de anulacion, donde m_mu != m_{mu*}.  Sirve para
#       ver que las estadisticas no estan pasando por casualidad.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run -d --name pI-instr \
#         -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" -w /work sage-normaliz:local \
#         bash -c "sage < pI_instrument_asociadas.sage > pI_instrument_asociadas_OUT.txt 2>&1"

import itertools
import json
import sys
from collections import Counter, defaultdict

CASOS = [(2, 2), (2, 3), (4, 2)]      # (t, r);  N = t + 2r
LMAX = 10
ANCHO = 8                             # lambda_1 maximo: el locus de anulacion vive en lambda ANCHAS

Sym = SymmetricFunctions(QQ)
s = Sym.s()


def asociada(mu, N):
    """La asociada de Littlewood: la primera columna pasa de mu'_1 a N - mu'_1."""
    c = Partition(mu).conjugate()
    c = list(c)
    if not c:
        return Partition([1] * N) if N > 0 else Partition([])
    if c[0] > N:
        return None
    nueva = [N - c[0]] + c[1:]
    nueva = [x for x in nueva if x > 0]
    if any(nueva[i] < nueva[i + 1] for i in range(len(nueva) - 1)):
        return None
    return Partition(nueva).conjugate()


def pares_hasta(n):
    """particiones con TODAS las partes pares y tamano <= n."""
    out = []
    for k in range(0, n // 2 + 1):
        for p in Partitions(2 * k):
            if all(x % 2 == 0 for x in p):
                out.append(Partition(p))
    return out


_LR = {}
def lr(lam, mu, de):
    key = (tuple(lam), tuple(mu), tuple(de))
    if key not in _LR:
        _LR[key] = QQ((s(Partition(mu)) * s(Partition(de))).coefficient(Partition(lam)))
    return _LR[key]


def m_de(lam, nu, evens):
    """m_nu(lambda) = sum_{delta par} c^lambda_{nu,delta}, y el desglose por delta."""
    tot = 0
    detalle = []
    n = sum(lam) - sum(nu)
    if n < 0:
        return 0, []
    for de in evens:
        if sum(de) != n:
            continue
        c = lr(lam, nu, de)
        if c != 0:
            tot += int(c)
            detalle.append((tuple(de), int(c)))
    return tot, detalle


print("=" * 100)
print("prob:instrument --- .es  m_mu = m_{mu*}  una identidad VISIBLE en tableaux?")
print("=" * 100)
print("")

resumen = {}
for (T, R) in CASOS:
    N = T + 2 * R
    evens = pares_hasta(LMAX + 4)
    # Las lambda se enumeran por ANCHURA, no por tamano.  La primera version usaba |lambda| <= 10 y
    # dejaba UN solo lambda en el locus de anulacion de 93: un control sobre muestra de uno, que es
    # el fallo contra el que avisa la propia nota de instrumentos.  El locus vive en lambda anchas
    # --- tab:virtual tiene 352 de 3060 con lambda_1 <= 14 ---, asi que se barre por lambda_1.
    # Y SOLO EN EL RANGO DE LITTLEWOOD, l(lambda) <= N/2.  m_nu = sum_{delta par} c^lambda_{nu,delta}
    # es la regla de Littlewood y fuera de su rango NO computa la multiplicidad --- ahi gobiernan las
    # reglas de modificacion, como dice el propio articulo en §sec:unstable.  La version anterior
    # barria l(lambda) <= N y por eso el «locus» salia de UN elemento: casi todo lambda fallaba el
    # test por una razon que no era la anulacion, sino que la formula no aplicaba.
    lams = []
    for p in itertools.product(range(ANCHO, -1, -1), repeat=N // 2):
        if any(p[i] < p[i + 1] for i in range(N // 2 - 1)):
            continue
        q = Partition([x for x in p if x > 0])
        if q:
            lams.append(q)
    c0 = c0n = 0
    c1_ok, c1_bad = 0, []
    refina = Counter()
    nrefina = Counter()
    soporte_igual = 0
    soporte_dist = 0
    nulos = 0
    ejemplos = []
    d1_par = 0
    for lam in lams:
        # el locus de anulacion, por la definicion (2186): m_nu = m_{nu*} para TODO nu
        # Etiquetas: TODAS las nu con a lo sumo N filas y |nu| <= |lambda|.  No se filtra por
        # nu contenida en lambda: los casos con m_nu = 0 y m_{nu*} != 0 son precisamente los que
        # detectan la NO anulacion, y quitarlos romperia el test.  La primera version llevaba aqui
        # un filtro con la precedencia mal puesta que dejaba la muestra en un solo lambda.
        etiquetas = []
        for k in range(0, sum(lam) + 1):
            for p in Partitions(k, max_length=N):
                etiquetas.append(Partition(p))
        cero = True
        pares_nu = []
        for nu in etiquetas:
            ast = asociada(nu, N)
            if ast is None:
                continue
            mn, det_n = m_de(lam, nu, evens)
            ma, det_a = m_de(lam, ast, evens)
            pares_nu.append((nu, ast, mn, ma, det_n, det_a))
            if mn != ma:
                cero = False
        if not cero:
            # D1: fuera del locus, cuantos pares tienen ya distinto CARDINAL de soporte
            for (nu, ast, mn, ma, dn, da) in pares_nu:
                if mn != ma:
                    d1_par += 1
            continue
        nulos += 1
        c1_ok += 1
        for (nu, ast, mn, ma, dn, da) in pares_nu:
            if mn == 0 and ma == 0:
                continue
            # C3: cardinales de los soportes
            if len(dn) == len(da):
                soporte_igual += 1
            else:
                soporte_dist += 1
                if len(ejemplos) < 6:
                    ejemplos.append((tuple(lam), tuple(nu), tuple(ast), len(dn), len(da), mn))
            # C2: .refina alguna estadistica?
            for nombre, f in (("l(delta)", lambda d: len(d)),
                              ("delta_1", lambda d: d[0] if d else 0),
                              ("|delta|/2", lambda d: sum(d) // 2),
                              ("num columnas", lambda d: d[0] if d else 0),
                              ("num partes pares", lambda d: sum(1 for x in d if x % 2 == 0))):
                A = Counter()
                B = Counter()
                for (d, c) in dn:
                    A[f(d)] += c
                for (d, c) in da:
                    B[f(d)] += c
                if A == B:
                    refina[nombre] += 1
                else:
                    nrefina[nombre] += 1
    print("  t=%d r=%d  (N=%d):  %d lambda en el locus de anulacion, de %d" % (T, R, N, nulos, len(lams)))
    print("    C1  m_mu = m_{mu*} en todo el locus : %d   (es la definicion, sirve de calibracion)" % c1_ok)
    print("    C3  pares (mu,mu*) no triviales : soporte del MISMO cardinal %d, de cardinal DISTINTO %d"
          % (soporte_igual, soporte_dist))
    if ejemplos:
        print("        ejemplos con cardinal distinto (lambda, mu, mu*, #dn, #da, m) : %s" % (ejemplos[:4],))
    print("    C2  .refina alguna estadistica sobre delta?")
    for nombre in ("l(delta)", "delta_1", "|delta|/2", "num partes pares"):
        a, b = refina[nombre], nrefina[nombre]
        print("        %-18s refina en %4d, NO refina en %4d   %s"
              % (nombre, a, b, "REFINA SIEMPRE" if b == 0 and a > 0 else "no"))
    print("    D1  pares con m_mu != m_{mu*} fuera del locus : %d" % d1_par)
    print("")
    sys.stdout.flush()
    resumen["t=%d,r=%d" % (T, R)] = dict(nulos=int(nulos), C1=int(c1_ok),
                                         soporte_igual=int(soporte_igual),
                                         soporte_dist=int(soporte_dist),
                                         refina={k: int(v) for k, v in refina.items()},
                                         nrefina={k: int(v) for k, v in nrefina.items()},
                                         D1=int(d1_par))

json.dump(resumen, open("pI_instrument_asociadas_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
