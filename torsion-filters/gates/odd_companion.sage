# -*- coding: utf-8 -*-
# M6 -- EL COMPAÑERO IMPAR:  la reduccion a t=1, y un filtro con UNA PARED MENOS.
# 15 de agosto de 2026.
#
# LA IDEA, Y ES NUESTRA.  Todo lo de hoy es t PAR:  mu_t = {+1,-1} u {xi^{+-j}}_{j=1..m},  o sea DOS
# puntos fijos de x -> 1/x mas m pares reciprocos, y de ahi la reduccion a t=2 con m=(t-2)/2.
#
# Para t IMPAR hay UN SOLO punto fijo:
#
#     mu_t = {1} u {xi^{+-j}}_{j=1..m'},        m' = (t-1)/2,
#
# luego la reduccion analoga no va a t=2 sino a  t=1,  con  R' = r + m'.  Y N cuadra:
# 1 + 2R' = 1 + 2r + (t-1) = t + 2r.  Mismo N, mismo beta.
#
# POR QUE ESTO PUEDE SER MEJOR QUE EL CASO PAR.  El filtro de torsion es
#
#     tau_t(eta) = sp_eta(xi, xi^2, ..., xi^{m'}),      a_j = eta_j + m' - j + 1
#
# y con t IMPAR la entrada xi^{i a} - xi^{-i a} se anula solo si 2 i a = 0 (mod t); como t es impar,
# 2 es invertible, luego la UNICA pared es  a = 0 (mod t).  **La pared a = t/2 no existe.**
# Quedan DOS especies de pared (a=0, y a_i = +-a_j) en vez de TRES.  El caso que no hemos mirado es
# estrictamente mas simple que el que llevamos todo el dia atacando.
#
# LO QUE SE MIDE
#   A  (T_impar):  tau != 0  <=>  a_1..a_{m'} ocupan una vez cada clase {+-1},...,{+-m'} mod t,
#      y entonces tau = +-1, con el signo cerrado.  Tres rutas independientes, como en el par.
#   B  la reduccion  Phi_{t,r} = Phi_{1,R'}|_{y=(xi,...,xi^{m'})}, y sobre todo la CUANTIFICACION:
#      que fraccion de los ceros es heredada de t=1 y que fraccion la crea la especializacion.
#      En el par salio 3-6 % heredado / 94-97 % creado.  Si en el impar sale muy distinto, la
#      dificultad no esta donde creemos.
#
# CONTROLES, y todos pueden fallar
#   C0  FATAL.  ruta caracter (Freudenthal) == ruta bialternante, en todos los eta.
#   C1  (T_impar) contra la verdad en LAS DOS direcciones, por separado.  n impreso siempre.
#   C2  SEÑUELO A: "distintos mod t y no nulos", SIN plegado +-.  Tiene que FALLAR.
#   C3  SEÑUELO B: se le añade una pared falsa en (t+1)/2, que en impar NO es pared.  Tiene que
#       FALLAR -- y si no fallara, es que esa pared si estaba actuando y me he equivocado.
#   C4  el signo predicho contra el azar (semilla fija).  El azar ~50 %.
#   C5  FATAL, parte B.  La reduccion, monomio a monomio, contra un SEÑUELO de especializacion
#       (y = (xi^2, ..., ) desplazada) que tiene que discrepar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage odd_companion.sage

import itertools, sys, json
from collections import defaultdict

# ================================================================== PARTE A =====================
def pesos(mm, tope):
    out = []
    def rec(pref, k, cap):
        if k == mm:
            out.append(tuple(pref)); return
        for v in range(cap, -1, -1):
            rec(pref + [v], k + 1, v)
    rec([], 0, tope)
    return out

_C = {}
def tau_caracter(eta, tt, mm):
    K = CyclotomicField(tt)
    xi = K.gen()
    key = (tuple(eta), mm)
    if key not in _C:
        W = WeylCharacterRing("C%d" % mm)
        el = W(W.space().from_vector(vector(list(eta))))
        _C[key] = [(tuple(int(v) for v in wt.to_vector()), mult)
                   for wt, mult in el.weight_multiplicities().items()]
    s = K(0)
    for wt, mult in _C[key]:
        s += mult * xi ** (sum(wt[i] * (i + 1) for i in range(mm)) % tt)
    return s

def tau_bialternante(eta, tt, mm):
    K = CyclotomicField(tt)
    xi = K.gen()
    a = [eta[j] + (mm - (j + 1) + 1) for j in range(mm)]
    d = [mm - (j + 1) + 1 for j in range(mm)]
    def dd(ex):
        return matrix(K, mm, mm,
                      lambda i, j: xi ** ((i + 1) * ex[j]) - xi ** (-(i + 1) * ex[j])).determinant()
    den = dd(d)
    if den == 0:
        return "DENOMINADOR NULO"
    return dd(a) / den

def regla_impar(eta, tt, mm):
    """(vive?, valor).  UNICA pared de tipo a=0; no hay t/2 porque t es impar."""
    a = [eta[j] + (mm - (j + 1) + 1) for j in range(mm)]
    cl, sg = [], 1
    for v in a:
        c = v % tt
        if c == 0:
            return False, 0
        if c <= mm:
            cl.append(c)
        else:
            cl.append(tt - c); sg *= -1
    if len(set(cl)) != mm:
        return False, 0
    perm = [mm - cl[j] for j in range(mm)]
    inv = sum(1 for i in range(mm) for j in range(i + 1, mm) if perm[i] > perm[j])
    return True, sg * (-1) ** inv

def senuelo_A(eta, tt, mm):
    a = [(eta[j] + (mm - (j + 1) + 1)) % tt for j in range(mm)]
    if any(c == 0 for c in a):
        return False
    return len(set(a)) == mm

def senuelo_B(eta, tt, mm):
    """añade una pared FALSA en (t+1)/2, que en t impar no es pared."""
    a = [(eta[j] + (mm - (j + 1) + 1)) % tt for j in range(mm)]
    if any(c == 0 or c == (tt + 1) // 2 for c in a):
        return False
    cl = [min(c, tt - c) for c in a]
    return len(set(cl)) == mm

print("=" * 122)
print("M6-A   EL FILTRO IMPAR  tau_t(eta) = sp_eta(xi,...,xi^{(t-1)/2})   --   UNA PARED MENOS")
print("=" * 122)
print("")
print("  t  m' |    n | C0 caract==bialt | (T) cero-falso / nocero-falso | (T') |tau|=1 | signo | SEÑ.A yerra | SEÑ.B yerra | azar")
print("  " + "-" * 148)
sys.stdout.flush()
set_random_seed(20260815)
TOT = defaultdict(int)
for (tt, tope) in [(3, 14), (5, 10), (7, 7), (9, 5), (11, 4)]:
    mm = (tt - 1) // 2
    n = c0 = fc = fnc = tp = sg = sa = sb = az = nz = 0
    for eta in pesos(mm, tope):
        n += 1
        v1 = tau_caracter(eta, tt, mm)
        v2 = tau_bialternante(eta, tt, mm)
        if v2 == "DENOMINADOR NULO" or v1 != v2:
            c0 += 1
        real = (v1 != 0)
        vive, pred = regla_impar(eta, tt, mm)
        fc += (real and not vive)
        fnc += (vive and not real)
        if real:
            nz += 1
            tp += (v1 not in (1, -1))
            sg += (vive and v1 == pred)
            az += (ZZ.random_element(0, 2) * 2 - 1 == (1 if v1 == 1 else -1))
        sa += (senuelo_A(eta, tt, mm) != real)
        sb += (senuelo_B(eta, tt, mm) != real)
    for k, v in [('n', n), ('c0', c0), ('fc', fc), ('fnc', fnc), ('tp', tp), ('sg', sg),
                 ('sa', sa), ('sb', sb), ('az', az), ('nz', nz)]:
        TOT[k] += v
    print("  %2d  %2d | %4d | %16s | %10d / %-18d | %14s | %2d/%-3d | %11d | %11d | %d/%d"
          % (tt, mm, n, "OK" if c0 == 0 else "*** %d ***" % c0, fc, fnc,
             "OK" if tp == 0 else "*** %d ***" % tp, sg, nz, sa, sb, az, nz))
    sys.stdout.flush()
print("  " + "-" * 148)
print("  TOTAL n=%d | C0 %d | (T) %d + %d | (T') %d | signo %d/%d | SEÑ.A yerra %d | SEÑ.B yerra %d | azar %d/%d"
      % (TOT['n'], TOT['c0'], TOT['fc'], TOT['fnc'], TOT['tp'], TOT['sg'], TOT['nz'],
         TOT['sa'], TOT['sb'], TOT['az'], TOT['nz']))
print("")
print("  SEÑ.B es el control de la idea entera: mete una pared en (t+1)/2, que en t IMPAR no existe.")
print("  Si yerra mucho, confirma que solo hay DOS especies de pared y el impar es mas simple.")
print("  Si NO yerra, esa pared si actuaba y la idea de M6 es falsa.")
print("")

# --------------------------------------------------------- A-bis: la estructura Z/2 del impar ----
# DE DONDE SALE.  Leyendo el §4 de Andersen-Stroppel por las formulas: para ell IMPAR en tipo C el
# alcove es  A = { m_1 + 2m_2 + ... < ell - 2n + 1 }, y con nuestros numeros (n = m', ell = t = 2m'+1)
# eso da  A = {0, omega_1}:  RANGO 2, o sea Z[J]/(J^2-1), no Z como en el caso par.
# ADVERTENCIA, y va antes que el dato: A-S avisan (Warning 2.1.1) que para ell impar su grupo de Weyl
# afin es el del sistema DUAL, y el amigo ya avisó en la vuelta 12 de que las normalizaciones B/C
# cambian.  Asi que esto es CONJETURA pendiente de normalizacion, no hecho.
#
# La evaluacion en g manda J -> -1, luego tau pierde la etiqueta del alcove.  Aqui se mide si algun
# invariante barato la recupera.  El candidato obvio, tau = (-1)^{|eta|}, YA se refuto a mano en
# t=5, eta=(2,1): tau=+1 y |eta|=3.  Se mide su tasa de acierto igualmente, para no fiarnos del
# calculo a mano, y se tabula la distribucion conjunta.
print("=" * 122)
print("M6-A bis   ¿RECUPERA ALGUN INVARIANTE BARATO LA ETIQUETA DEL ALCOVE?  (conjetura de rango 2)")
print("=" * 122)
print("")
print("   t | supervivientes | tau=(-1)^{|eta|} acierta | distribucion conjunta (|eta| mod 2, tau)")
print("   " + "-" * 110)
for (tt, tope) in [(3, 14), (5, 10), (7, 7), (9, 5)]:
    mm = (tt - 1) // 2
    conj = defaultdict(int)
    acierta = nz = 0
    for eta in pesos(mm, tope):
        vive, val = regla_impar(eta, tt, mm)
        if not vive:
            continue
        nz += 1
        par = sum(eta) % 2
        conj[(par, val)] += 1
        acierta += (val == (-1) ** sum(eta))
    print("   %2d | %14d | %10d / %-11d | %s"
          % (tt, nz, acierta, nz, dict(sorted(conj.items()))))
print("")
print("   Si la columna del medio diera nz/nz, la paridad de |eta| SERIA la etiqueta.  Si no, hace")
print("   falta el plegado afin de verdad (en la normalizacion DUAL) y no un atajo.")
print("")
sys.stdout.flush()

# ================================================================== PARTE B =====================
def phi_de(beta, tt, nvar):
    """Phi con bloque congelado = raices tt-esimas de la unidad, y nvar pares libres."""
    N = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(1 if tt == 1 else -1)
    L = LaurentPolynomialRing(K, nvar, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(N - 1, -1, -1))
    def alt(ex):
        return matrix(L, N, N, lambda i, j: x[i] ** ex[j]).determinant()
    den = alt(delta)
    if den == 0:
        return None
    q = alt(list(beta)) / den
    try:
        q = L(q)
    except Exception:
        return "NO-POL"
    return {tuple(e) if hasattr(e, '__iter__') else (e,): c
            for e, c in zip(q.exponents(), q.coefficients()) if c != 0}

def especializar(P, tt, mm, r, desplazar=0):
    """y = (xi^{1+d}, ..., xi^{m'+d}) sobre las PRIMERAS m' variables; deja r libres."""
    K = CyclotomicField(tt)
    xi = K.gen()
    L = LaurentPolynomialRing(K, r, 'w')
    ws = L.gens()
    out = {}
    for e, c in P.items():
        val = K(c) * prod(xi ** ((j + 1 + desplazar) * e[j]) for j in range(mm))
        k = tuple(e[mm:])
        out[k] = out.get(k, 0) + val
    return {k: v for k, v in out.items() if v != 0}

print("=" * 122)
print("M6-B   LA REDUCCION AL CASO t=1,  y  QUE FRACCION DE LOS CEROS LA CREA LA ESPECIALIZACION")
print("=" * 122)
print("")
# LO QUE LA PARTE B MIDE Y LO QUE NO, corregido tras la primera corrida.
#   La identidad  Phi_{t,r} = Phi_{1,R'}|_{y=(xi,...,xi^{m'})}  es un REAGRUPAMIENTO del mismo
#   conjunto de puntos: mu_t = {1} u {xi^{+-j}}.  Luego cualquier implementacion fiel TIENE que dar
#   260/260, y NO es un descubrimiento: es un control de implementacion.  El primer señuelo que puse
#   -- desplazar la especializacion a (xi^2,...) -- empato 0/260 porque para m'=1 eso es xi^{-1}, la
#   MISMA especializacion vista por la simetria reciproca.  Un señuelo que no puede fallar no mide
#   ([[a-decoy-that-ties-means-untested]]).
#   El señuelo bueno es CRUZADO EN beta: comparar esp(beta) contra Phi_{t,r}(beta') de otra forma.
#   Si eso coincidiera, la comparacion entera seria vacua.
#   Y el contenido de la parte B no es la identidad: es la CUANTIFICACION heredados / creados.
print("   t  r  N |  n formas | reduccion exacta | SEÑUELO cruzado en beta | Phi=0 | de esos: heredados de t=1 | CREADOS por la torsion")
print("   " + "-" * 148)
sys.stdout.flush()

for (tt, r, W, TOPE) in [(3, 2, 11, 260), (5, 2, 12, 160), (7, 2, 13, 90)]:
    mm = (tt - 1) // 2
    Rp = r + mm
    N = tt + 2 * r
    n = ok = sen_dif = ceros = hered = 0
    prev_Pt = None
    for comb in itertools.combinations(range(W + 1), N):
        beta = tuple(sorted(comb, reverse=True))
        if n >= TOPE:
            break
        P1 = phi_de(beta, 1, Rp)                 # el objeto de t=1, con R' pares
        if P1 is None or P1 == "NO-POL":
            continue
        Pt = phi_de(beta, tt, r)                 # el objeto real
        if Pt is None or Pt == "NO-POL":
            continue
        n += 1
        esp = especializar(P1, tt, mm, r)
        ok += (esp == Pt)
        # SEÑUELO CRUZADO: esp(beta) contra el Phi de la beta ANTERIOR.  Tiene que discrepar; si
        # coincidiera, la comparacion de la columna anterior no distinguiria formas.
        if prev_Pt is not None:
            sen_dif += (esp != prev_Pt)
        prev_Pt = Pt
        if not Pt:
            ceros += 1
            hered += (not P1)                    # ya era 0 antes de especializar
    if not n:
        print("   %2d %2d %2d | POBLACION VACIA" % (tt, r, N))
        continue
    creados = ceros - hered
    print("   %2d %2d %2d | %9d | %16s | %23s | %5d | %13d (%4.1f %%) | %d (%4.1f %%)"
          % (tt, r, N, n, "%d/%d" % (ok, n), "%d/%d discrepa" % (sen_dif, n - 1), ceros,
             hered, 100.0 * hered / max(ceros, 1), creados, 100.0 * creados / max(ceros, 1)))
    sys.stdout.flush()

print("")
print("   En el caso PAR salio: 3-6 % heredado, 94-97 % creado por la especializacion.")
print("   Si el impar reparte muy distinto, la dificultad no esta donde creemos.")
print("")
print("=" * 122)
print("DONE")
