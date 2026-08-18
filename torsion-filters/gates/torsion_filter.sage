# -*- coding: utf-8 -*-
# EL FILTRO DE TORSION  tau_t(eta) = sp_eta(xi, xi^2, ..., xi^m).   15 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 12 de la consulta externa.  Afirma DOS cosas, sin conjetura segun el:
#
#   (T)   tau_t(eta) != 0  <=>  a_1..a_m ocupan una vez cada clase {+-1},...,{+-m} mod t,
#                              con  a_j = eta_j + m - j + 1,   t = 2m+2,  xi = e^{2 pi i / t}
#   (T')  y cuando ocurre,  tau_t(eta) = +-1.
#
# Su lectura: pared => 0, fuera de paredes => plegado con signo.  Es el mecanismo de Kac-Walton /
# Andersen-Stroppel, verificado en su texto (12_IN, tabla de verificacion).
#
# POR QUE ESTO IMPORTA.  Nuestra reduccion dice  Phi_{t,r} = Phi_{2,R}|_{y=(xi,...,xi^m)}  con
# R = r+m, y el 94-97% de los ceros los crea esa especializacion.  Si (T) vale, el bloque de raices
# de la unidad esta RESUELTO y lo que queda es como interactua con el branching B_{eta,mu}.
#
# TRES RUTAS INDEPENDIENTES, no dos:
#   R1  el caracter de Sage:  sum sobre weight_multiplicities de WeylCharacterRing("C_m"), evaluado
#       en xi^i.  Viene de Freudenthal, NO del bialternante.  Es la definicion de sp_eta.
#   R2  el bialternante de tipo C_m, que es de donde EL deduce la regla.
#   R3  su regla (T)+(T'), y ademas el SIGNO exacto predicho -- no solo |tau|=1.
#
# CONTROLES
#   C0  FATAL.  R1 == R2 en todos los eta.  Si no, una de las dos maquinarias esta mal.
#   C1  la regla contra la verdad, en LAS DOS DIRECCIONES: falsos ceros y falsos no-ceros, por
#       separado.  Se imprime n siempre.
#   C2  SEÑUELO A -- "distintos mod t y no nulos" SIN el plegado +-.  Tiene que FALLAR.
#   C3  SEÑUELO B -- pared solo en 0, sin t/2.  Tiene que FALLAR.
#   C4  SEÑUELO C -- el signo predicho contra un signo al azar (semilla fija).  El azar tiene que
#       acertar ~50%.  Si el azar tambien acierta 100%, la columna no mide nada.
#   C5  el caso t=4 cerrado a mano:  tau_4(k) = (-1)^{k/2} si k par, 0 si impar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage torsion_filter.sage

import itertools

# ---------------------------------------------------------------- enumeracion de pesos ----------
def pesos(m, tope):
    """todos los eta dominantes de C_m con eta_1 <= tope:  eta_1 >= ... >= eta_m >= 0."""
    out = []
    def rec(pref, k, cap):
        if k == m:
            out.append(tuple(pref)); return
        for v in range(cap, -1, -1):
            rec(pref + [v], k + 1, v)
    rec([], 0, tope)
    return out


# ---------------------------------------------------------------- R1: el caracter de Sage -------
_CACHE = {}
def tau_caracter(eta, t, m):
    """sp_eta(xi,...,xi^m) sumando las multiplicidades de peso.  Freudenthal, no bialternante."""
    K = CyclotomicField(t)
    xi = K.gen()
    key = (tuple(eta), m)
    if key not in _CACHE:
        W = WeylCharacterRing("C%d" % m)
        el = W(W.space().from_vector(vector(list(eta))))
        _CACHE[key] = [(tuple(wt.to_vector()), mult)
                       for wt, mult in el.weight_multiplicities().items()]
    s = K(0)
    for wt, mult in _CACHE[key]:
        e = sum(int(wt[i]) * (i + 1) for i in range(m)) % t
        s += mult * xi ** e
    return s


# ---------------------------------------------------------------- R2: el bialternante -----------
def tau_bialternante(eta, t, m):
    """sp_eta = det(x_i^{a_j} - x_i^{-a_j}) / det(x_i^{d_j} - x_i^{-d_j})  en x_i = xi^i."""
    K = CyclotomicField(t)
    xi = K.gen()
    a = [eta[j] + (m - (j + 1) + 1) for j in range(m)]   # a_j = eta_j + m - j + 1, con j desde 1
    d = [m - (j + 1) + 1 for j in range(m)]
    def det_de(expos):
        return matrix(K, m, m,
                      lambda i, j: xi ** ((i + 1) * expos[j]) - xi ** (-(i + 1) * expos[j])
                      ).determinant()
    den = det_de(d)
    if den == 0:
        return "DENOMINADOR NULO"
    return det_de(a) / den


# ---------------------------------------------------------------- R3: su regla ------------------
def clase(a, t, m):
    """(clase en 1..m, signo) o None si cae en pared."""
    c = a % t
    if c == 0 or 2 * c == t:
        return None
    return (c, +1) if c <= m else (t - c, -1)

def regla(eta, t, m):
    """(sobrevive?, valor predicho).  Valor = sgn(permutacion) * producto de signos del plegado."""
    a = [eta[j] + (m - (j + 1) + 1) for j in range(m)]
    cs = [clase(v, t, m) for v in a]
    if any(c is None for c in cs):
        return False, 0
    clases = [c for c, _ in cs]
    if len(set(clases)) != m:
        return False, 0
    signo = 1
    for _, s in cs:
        signo *= s
    # signo de la permutacion  j -> clases[j]  (llevar las columnas al orden del denominador m..1)
    perm = [m - clases[j] for j in range(m)]         # denominador tiene d_j = m-j+1, o sea clases m..1
    inv = sum(1 for i in range(m) for j in range(i + 1, m) if perm[i] > perm[j])
    signo *= (-1) ** inv
    return True, signo


# ---------------------------------------------------------------- SEÑUELOS ----------------------
def senuelo_A(eta, t, m):
    """SIN plegado +-: 'a_j distintos mod t y ninguno 0 ni t/2'."""
    a = [(eta[j] + (m - (j + 1) + 1)) % t for j in range(m)]
    if any(c == 0 or 2 * c == t for c in a):
        return False
    return len(set(a)) == m

def senuelo_B(eta, t, m):
    """pared solo en 0 (sin t/2), con plegado."""
    a = [(eta[j] + (m - (j + 1) + 1)) % t for j in range(m)]
    if any(c == 0 for c in a):
        return False
    cl = [min(c, t - c) for c in a]
    return len(set(cl)) == m


# ================================================================= C5 ===========================
print("=" * 104)
print("C5  EL CASO t=4 A MANO:  su formula  tau_4(k) = (-1)^{k/2} si k par, 0 si impar")
print("=" * 104)
print("")
K4 = CyclotomicField(4); i4 = K4.gen()
malo5 = 0
for k in range(0, 13):
    real = tau_caracter((k,), 4, 1)
    suyo = K4((-1) ** (k // 2)) if k % 2 == 0 else K4(0)
    ok = (real == suyo)
    malo5 += (not ok)
    print("    k=%2d   tau_4(k) = %-6s   el dice %-6s   %s" % (k, real, suyo, "ok" if ok else "*** NO ***"))
print("")
print("    C5 %s   (n = 13)" % ("PASA" if not malo5 else "*** FALLA en %d de 13 ***" % malo5))

# ================================================================= C0, C1, C2, C3, C4 ===========
CONFIG = [(4, 14), (6, 10), (8, 7), (10, 5), (12, 4)]     # (t, tope de eta_1)
set_random_seed(20260815)

print("")
print("=" * 104)
print("(T) y (T') MEDIDAS.   n se imprime SIEMPRE.")
print("=" * 104)
print("")
print("  t   m |    n | C0 R1==R2 | (T) fallos: cero-falso / nocero-falso | (T') |tau|=1 | signo exacto | SEÑUELO A | SEÑUELO B | azar")
print("  " + "-" * 148)

TOT = dict(n=0, c0=0, t_fc=0, t_fnc=0, tp=0, sg=0, sa=0, sb=0, az=0, nz=0)
DETALLE = []
for (t, tope) in CONFIG:
    m = (t - 2) // 2
    P = pesos(m, tope)
    n = c0 = fc = fnc = tp = sg = sa = sb = az = nz = 0
    for eta in P:
        n += 1
        v1 = tau_caracter(eta, t, m)
        v2 = tau_bialternante(eta, t, m)
        if v2 == "DENOMINADOR NULO" or v1 != v2:
            c0 += 1
        vive_real = (v1 != 0)
        vive, pred = regla(eta, t, m)
        if vive_real and not vive:
            fc += 1                       # la regla dice 0 y no lo es
        if vive and not vive_real:
            fnc += 1                      # la regla dice no-cero y es 0
        if vive_real:
            nz += 1
            if v1 not in (1, -1):
                tp += 1                   # (T') falla
            if vive and v1 == pred:
                sg += 1                   # el signo exacto acierta
            if ZZ.random_element(0, 2) * 2 - 1 == (1 if v1 == 1 else -1):
                az += 1                   # el azar, control C4
        if senuelo_A(eta, t, m) != vive_real:
            sa += 1
        if senuelo_B(eta, t, m) != vive_real:
            sb += 1
    for k, v in [('n', n), ('c0', c0), ('t_fc', fc), ('t_fnc', fnc), ('tp', tp), ('sg', sg),
                 ('sa', sa), ('sb', sb), ('az', az), ('nz', nz)]:
        TOT[k] += v
    DETALLE.append((t, m, n, nz))
    print("  %2d  %2d | %4d | %9s | %10d / %-25d | %14s | %5d/%-6d | %9d | %9d | %d/%d"
          % (t, m, n, "OK" if c0 == 0 else "*** %d ***" % c0, fc, fnc,
             "OK" if tp == 0 else "*** %d ***" % tp, sg, nz, sa, sb, az, nz))

print("  " + "-" * 148)
print("  TOTAL  | %4d | C0 fallos %d | (T) fallos %d + %d | (T') fallos %d | signo %d/%d | señuelo A yerra %d | señuelo B yerra %d | azar %d/%d"
      % (TOT['n'], TOT['c0'], TOT['t_fc'], TOT['t_fnc'], TOT['tp'], TOT['sg'], TOT['nz'],
         TOT['sa'], TOT['sb'], TOT['az'], TOT['nz']))
print("")
print("  LECTURA DE LOS SEÑUELOS.  A y B tienen que dar un numero GRANDE de errores.  Si dieran 0,")
print("  la condicion de (T) no estaria diciendo nada que no dijera una version mas floja, y no")
print("  habriamos medido su contenido.  El azar tiene que rondar el 50% del signo.")
print("")
print("  n por configuracion: " + ", ".join("t=%d: %d pesos (%d no nulos)" % (t, m and n or n, nz)
                                            for (t, m, n, nz) in DETALLE))
print("")
print("=" * 104)
print("DONE")
