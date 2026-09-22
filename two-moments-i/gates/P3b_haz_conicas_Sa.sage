# -*- coding: utf-8 -*-
# P3b_haz_conicas_Sa.sage -- el haz de conicas de la superficie cubica S_a, verificado paso a paso.
#
# Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
# What: verifica paso a paso, sobre Q y sobre F_p, la derivacion del recuento
#       #S_a(F_p) = p^2 + [3 + 3*chi_p(D_a) + chi_p(E_a)]*p + 1  para el haz de conicas de
#       S_a = {sum w_i x_i = 0} cap {sum w_i x_i^3 = 0} en P^4 con w = (a,1,1,1,1), o sea
#       a^2(x1^3+x2^3+x3^3+x4^3) - (x1+x2+x3+x4)^3 = 0 en P^3.
# Uso:  sage P3b_haz_conicas_Sa.sage > P3b_haz_conicas_Sa_OUT.txt 2>&1
#       (desde el directorio gates/; pide SageMath)
#
# LOS PASOS QUE SE VERIFICAN, cada uno por separado:
#   1. s=x1+x2, u=x1-x2, t=x3+x4, v=x3-x4  =>  a^2(s^3+t^3+3su^2+3tv^2) - 4(s+t)^3 = 0.
#   2. la recta L : s = t = 0 esta en S_a.
#   3. los planos t = r s dan conicas residuales 3a^2 u^2 + 3a^2 r v^2 + A_a(r) s^2 = 0 con
#      A_a(r) = (r+1)[(a^2-4)r^2 - (a^2+8)r + (a^2-4)].
#   4. las cinco fibras degeneradas son r = 0, infinito, -1, y las dos raices r+- .
#   5. D_a = -3(a^2-4), E_a = (a^2-4)(a^2-16); las raices que hacen falta son
#      ( sqrt(E_a) +- a sqrt(D_a) ) / (2(a^2-4)) y generan exactamente Q(sqrt D_a, sqrt E_a).
#   6. el recuento fibra a fibra da p^2 + [3 + 3 chi(D_a) + chi(E_a)] p + 1 si p no divide
#      6 a (a^2-4)(a^2-16).
#
# El veredicto es la ultima linea (comprobaciones y fallos), no el codigo de salida; y la salida
# no debe contener ningun Traceback.
OK = [0]
FALLOS = [0]

def chk(nombre, cond, extra=""):
    if cond:
        OK[0] += 1
        print("OK    " + nombre + ("  " + extra if extra else ""))
    else:
        FALLOS[0] += 1
        print("FALLO " + nombre + ("  " + extra if extra else ""))

def info(txt):
    print("      " + txt)

print("=== P3b_haz_conicas_Sa : el haz de conicas de S_a ===")
AS = [1, 3, 5, 6, 7, 9]
print("valores de a que se usan en las partes simbolicas y de recuento: " + str(AS))

# ---------------------------------------------------------------- PASO 1, 2, 3 simbolicos
R = PolynomialRing(QQ, ['a', 's', 'u', 't', 'v', 'r'])
a, s, u, t, v, r = R.gens()
x1 = (s + u) / 2
x2 = (s - u) / 2
x3 = (t + v) / 2
x4 = (t - v) / 2
F_orig = a**2 * (x1**3 + x2**3 + x3**3 + x4**3) - (x1 + x2 + x3 + x4)**3
F4 = a**2 * (s**3 + t**3 + 3 * s * u**2 + 3 * t * v**2) - 4 * (s + t)**3
print("--- paso 1: cambio de variables (simbolico en a) ---")
chk("paso1 simbolico: 4*F(x) = a^2(s^3+t^3+3su^2+3tv^2)-4(s+t)^3", R(4 * F_orig) == F4)

def paso1_en_a(av):
    Rx = PolynomialRing(QQ, ['s', 'u', 't', 'v'])
    ss, uu, tt, vv = Rx.gens()
    y1 = (ss + uu) / 2
    y2 = (ss - uu) / 2
    y3 = (tt + vv) / 2
    y4 = (tt - vv) / 2
    G = av**2 * (y1**3 + y2**3 + y3**3 + y4**3) - (y1 + y2 + y3 + y4)**3
    H = av**2 * (ss**3 + tt**3 + 3 * ss * uu**2 + 3 * tt * vv**2) - 4 * (ss + tt)**3
    return Rx(4 * G) == H

for av in AS:
    chk("paso1 a=" + str(av), paso1_en_a(av))

print("--- paso 2: la recta L : s = t = 0 esta contenida en S_a ---")
chk("paso2 simbolico: F4|_{s=t=0} = 0", F4.subs({s: R(0), t: R(0)}) == R(0))
for av in AS:
    Fa = F4.subs({a: R(av)})
    chk("paso2 a=" + str(av), Fa.subs({s: R(0), t: R(0)}) == R(0))

print("--- paso 3: planos t = r s, conica residual y factorizacion de A_a(r) ---")
A_decl = (r + 1) * ((a**2 - 4) * r**2 - (a**2 + 8) * r + (a**2 - 4))
A_nat = a**2 * (1 + r**3) - 4 * (1 + r)**3
chk("paso3 simbolico: A_a(r) = a^2(1+r^3)-4(1+r)^3 = (r+1)[(a^2-4)r^2-(a^2+8)r+(a^2-4)]", A_nat == A_decl)
G_sub = F4.subs({t: r * s})
G_decl = s * (3 * a**2 * u**2 + 3 * a**2 * r * v**2 + A_decl * s**2)
chk("paso3 simbolico: F4(s,u,rs,v) = s*(3a^2u^2+3a^2 r v^2+A_a(r)s^2)", G_sub == G_decl)
for av in AS:
    chk("paso3 a=" + str(av), G_sub.subs({a: R(av)}) == G_decl.subs({a: R(av)}))

print("--- paso 3bis: la fibra r = infinito (plano s = 0) ---")
Ginf = F4.subs({s: R(0)})
chk("paso3bis simbolico: F4(0,u,t,v) = t*((a^2-4)t^2+3a^2 v^2)", Ginf == t * ((a**2 - 4) * t**2 + 3 * a**2 * v**2))

# ---------------------------------------------------------------- PASO 4 simbolico
print("--- paso 4: el determinante del haz se anula exactamente en 0, -1, r+- (y infinito) ---")
det_haz = A_decl * (3 * a**2) * (3 * a**2 * r)
chk("paso4 simbolico: det conica = 9a^4 * r * A_a(r)", det_haz == 9 * a**4 * r * A_decl)
q_sym = (a**2 - 4) * r**2 - (a**2 + 8) * r + (a**2 - 4)
chk("paso4 simbolico: q(0) = a^2-4 y q(-1) = 3a^2, luego r+- distintas de 0 y de -1 si a(a^2-4) != 0",
    q_sym.subs({r: R(0)}) == a**2 - 4 and q_sym.subs({r: R(-1)}) == 3 * a**2)
chk("paso4 simbolico: las CINCO fibras malas son distintas salvo si a(a^2-4)(a^2-16) = 0 (disc = 3a^2(16-a^2))",
    q_sym.discriminant(r) == 3 * a**2 * (16 - a**2))
Rr = PolynomialRing(QQ, 'z')
z = Rr.gen()
for av in AS:
    q = (av**2 - 4) * z**2 - (av**2 + 8) * z + (av**2 - 4)
    A_z = (z + 1) * q
    pol = z * A_z
    natural = all(A_z(k) == av**2 * (1 + k**3) - 4 * (1 + k)**3 for k in range(-5, 6))
    chk("paso4 a=" + str(av) + ": det/(9a^4) = r*A_a(r) de grado 4, raices 0,-1,r+-; disc(q)=3a^2(16-a^2)",
        natural and pol.degree() == 4 and q.degree() == 2 and q.discriminant() == 3 * av**2 * (16 - av**2) and A_z(-1) == 0,
        "disc=" + str(q.discriminant()) + " raices_q_en_Q=" + str([ro for ro, m in q.roots()]))

# ---------------------------------------------------------------- PASO 5
print("--- paso 5: discriminante, campos Q(sqrt D_a), Q(sqrt E_a) y el compuesto ---")
Da_sym = -3 * (a**2 - 4)
Ea_sym = (a**2 - 4) * (a**2 - 16)
disc_sym = (a**2 + 8)**2 - 4 * (a**2 - 4)**2
chk("paso5 simbolico: disc = 3a^2(16-a^2)", disc_sym == 3 * a**2 * (16 - a**2))
chk("paso5 simbolico: disc*(a^2-4)^2 = a^2*D_a*E_a  (luego Q(sqrt disc) = Q(sqrt(D_a E_a)))",
    disc_sym * (a**2 - 4)**2 == a**2 * Da_sym * Ea_sym)

def paso5_en_a(av):
    D = -3 * (av**2 - 4)
    E = (av**2 - 4) * (av**2 - 16)
    sD = QQbar(D).sqrt()
    sE = QQbar(E).sqrt()
    den = 2 * (av**2 - 4)
    Xp = (sE + av * sD) / den
    Xm = (sE - av * sD) / den
    Pz = PolynomialRing(QQbar, 'w')
    w = Pz.gen()
    q = (av**2 - 4) * w**2 - (av**2 + 8) * w + (av**2 - 4)
    raices = [ro for ro, m in q.roots()]
    c1 = any(Xp**2 == -ro for ro in raices)
    c2 = any(Xm**2 == -ro for ro in raices)
    c3 = (Xp**2 != Xm**2) or (len(raices) == 1)
    c4 = (Xp * Xm == QQbar(1))
    c5 = ((Xp + Xm) * (av**2 - 4) == sE)
    c6 = ((Xp - Xm) * (av**2 - 4) == av * sD)
    return (D, E, c1, c2, c3, c4, c5, c6)

for av in AS:
    D, E, c1, c2, c3, c4, c5, c6 = paso5_en_a(av)
    chk("paso5 a=" + str(av) + ": X+^2 y X-^2 son -r-, -r+", c1 and c2 and c3,
        "D=" + str(D) + " E=" + str(E) + " sqfD=" + str(QQ(D).squarefree_part()) + " sqfE=" + str(QQ(E).squarefree_part()))
    chk("paso5 a=" + str(av) + ": X+ X- = 1, (X+ + X-)(a^2-4) = sqrt E, (X+ - X-)(a^2-4) = a sqrt D  [=> Q(X+,X-) = Q(sqrt D, sqrt E)]",
        c4 and c5 and c6)

def grado_compuesto(av):
    D = QQ(-3 * (av**2 - 4)).squarefree_part()
    E = QQ((av**2 - 4) * (av**2 - 16)).squarefree_part()
    gens = set([1])
    for g in [D, E]:
        nuevos = set()
        for h in gens:
            nuevos.add(QQ(h * g).squarefree_part())
        gens = gens | nuevos
    return len(gens)

for av in AS:
    info("a=" + str(av) + ": [Q(sqrt D_a, sqrt E_a):Q] = " + str(grado_compuesto(av)))

# ---------------------------------------------------------------- PASO 6: recuento real
print("--- paso 6: recuento sobre F_p (fuerza bruta en P^3, en P^4, y fibra a fibra) ---")

def pts_P3(p):
    for y in range(p):
        for zz in range(p):
            for w in range(p):
                yield (1, y, zz, w)
    for zz in range(p):
        for w in range(p):
            yield (0, 1, zz, w)
    for w in range(p):
        yield (0, 0, 1, w)
    yield (0, 0, 0, 1)

def pts_P2(p):
    for y in range(p):
        for zz in range(p):
            yield (1, y, zz)
    for zz in range(p):
        yield (0, 1, zz)
    yield (0, 0, 1)

def pts_P4(p):
    for y in range(p):
        for zz in range(p):
            for w in range(p):
                for q in range(p):
                    yield (1, y, zz, w, q)
    for zz in range(p):
        for w in range(p):
            for q in range(p):
                yield (0, 1, zz, w, q)
    for w in range(p):
        for q in range(p):
            yield (0, 0, 1, w, q)
    for q in range(p):
        yield (0, 0, 0, 1, q)
    yield (0, 0, 0, 0, 1)

def cuenta_P3(av, p):
    p = int(p)
    av = int(av)
    cu = [(x * x * x) % p for x in range(p)]
    a2 = (av * av) % p
    tot = 0
    for pt in pts_P3(p):
        y1, y2, y3, y4 = pt
        ssum = (y1 + y2 + y3 + y4) % p
        val = (a2 * (cu[y1] + cu[y2] + cu[y3] + cu[y4]) - cu[ssum]) % p
        if val == 0:
            tot += 1
    return tot

def cuenta_P4(av, p):
    p = int(p)
    av = int(av)
    cu = [(x * x * x) % p for x in range(p)]
    tot = 0
    for pt in pts_P4(p):
        y0, y1, y2, y3, y4 = pt
        if (av * y0 + y1 + y2 + y3 + y4) % p != 0:
            continue
        if (av * cu[y0] + cu[y1] + cu[y2] + cu[y3] + cu[y4]) % p != 0:
            continue
        tot += 1
    return tot

def cuenta_conica(c1, c2, c3, p):
    p = int(p)
    c1 = int(c1) % p
    c2 = int(c2) % p
    c3 = int(c3) % p
    sq = [(x * x) % p for x in range(p)]
    tot = 0
    for pt in pts_P2(p):
        X, Y, Z = pt
        if (c1 * sq[X] + c2 * sq[Y] + c3 * sq[Z]) % p == 0:
            tot += 1
    return tot

def A_de(av, rr):
    return (rr + 1) * ((av**2 - 4) * rr**2 - (av**2 + 8) * rr + (av**2 - 4))

def chi(n, p):
    return kronecker_symbol(Integer(n) % p, p)

def raices_cuad(av, p):
    Pz = PolynomialRing(GF(p), 'w')
    w = Pz.gen()
    q = (av**2 - 4) * w**2 - (av**2 + 8) * w + (av**2 - 4)
    return sorted([Integer(ro) for ro, m in q.roots()])

def analisis_fibras(av, p):
    D = -3 * (av**2 - 4)
    E = (av**2 - 4) * (av**2 - 16)
    rr_malas = raices_cuad(av, p)
    total = 0
    malas_detectadas = []
    fallos_fibra = []
    for rr in range(p):
        c_s = Integer(A_de(av, rr)) % p
        c_u = (3 * av**2) % p
        c_v = (3 * av**2 * rr) % p
        det = (c_s * c_u * c_v) % p
        n = cuenta_conica(c_s, c_u, c_v, p)
        total += n
        if det == 0:
            malas_detectadas.append(rr)
        if det != 0:
            pred = p + 1
        elif rr == 0:
            pred = p + 1 + chi(D, p) * p
        elif rr == p - 1:
            pred = 2 * p + 1
        else:
            pred = p + 1 + chi(-rr, p) * p
        if n != pred:
            fallos_fibra.append(("r=" + str(rr), n, pred))
    n_inf = cuenta_conica(0, (av**2 - 4) % p, (3 * av**2) % p, p)
    total += n_inf
    pred_inf = p + 1 + chi(D, p) * p
    if n_inf != pred_inf:
        fallos_fibra.append(("r=inf", n_inf, pred_inf))
    malas_detectadas.append("inf")
    return (total, malas_detectadas, fallos_fibra, rr_malas)

PRIMOS = [5, 7, 11, 13, 17, 19, 23, 29, 31]

def malos_primos(av):
    N = 6 * av * (av**2 - 4) * (av**2 - 16)
    return [q for q in PRIMOS if N % q == 0]

for av in AS:
    D = -3 * (av**2 - 4)
    E = (av**2 - 4) * (av**2 - 16)
    excl = malos_primos(av)
    info("a=" + str(av) + ": D_a=" + str(D) + " E_a=" + str(E) + " primos excluidos en la lista: " + str(excl))
    for p in PRIMOS:
        if p in excl:
            continue
        nP3 = cuenta_P3(av, p)
        formula = p**2 + (3 + 3 * chi(D, p) + chi(E, p)) * p + 1
        total, malas, fallos_fibra, rr_malas = analisis_fibras(av, p)
        esperadas = sorted(set([0, p - 1] + rr_malas))
        chk("recuento a=" + str(av) + " p=" + str(p) + ": bruto en P^3 = formula",
            nP3 == formula,
            "bruto=" + str(nP3) + " formula=" + str(formula) + " chi(D)=" + str(chi(D, p)) + " chi(E)=" + str(chi(E, p)))
        chk("fibras a=" + str(av) + " p=" + str(p) + ": suma de fibras = bruto en P^3", total == nP3,
            "suma=" + str(total))
        chk("fibras a=" + str(av) + " p=" + str(p) + ": cada fibra mala aporta lo declarado",
            len(fallos_fibra) == 0, str(fallos_fibra))
        chk("fibras a=" + str(av) + " p=" + str(p) + ": las fibras degeneradas son exactamente {0,-1,inf,r+-}",
            sorted([m for m in malas if m != "inf"]) == esperadas,
            "detectadas=" + str(malas) + " esperadas=" + str(esperadas + ["inf"]))
        c_pareja = (chi(D, p) + chi(E, p))
        if len(rr_malas) == 2:
            aporte = chi(-rr_malas[0], p) + chi(-rr_malas[1], p)
        elif len(rr_malas) == 0:
            aporte = 0
        else:
            aporte = None
        chk("paso6 a=" + str(av) + " p=" + str(p) + ": aporte de r+- = chi(D)+chi(E)",
            aporte is not None and aporte == c_pareja,
            "aporte=" + str(aporte) + " chi(D)+chi(E)=" + str(c_pareja) + " r+-=" + str(rr_malas))

print("--- segundo metodo independiente: recuento directo en P^4 con los dos pesos ---")
for av in [1, 3, 5, 7]:
    for p in [5, 7, 11, 13]:
        if p in malos_primos(av):
            continue
        n4 = cuenta_P4(av, p)
        n3 = cuenta_P3(av, p)
        chk("P^4 vs P^3 a=" + str(av) + " p=" + str(p), n4 == n3, "P4=" + str(n4) + " P3=" + str(n3))

print("--- control externo: a = 1 es la cubica de Clebsch, 27 rectas sobre Q(sqrt 5) ---")
for p in [7, 11, 13, 17, 19, 23, 29, 31]:
    n = cuenta_P3(1, p)
    esperado = p**2 + 7 * p + 1 if p % 5 in [1, 4] else p**2 + 5 * p + 1
    chk("Clebsch p=" + str(p), n == esperado, "n=" + str(n) + " esperado=" + str(esperado))

# ---------------------------------------------------------------- CONTROL QUE PUEDE FALLAR
print("--- control (4): a = 2 y a = 4 estan excluidos; se comprueba que la construccion se rompe ---")
Rx4 = PolynomialRing(QQ, ['x1', 'x2', 'x3', 'x4'])
xx = Rx4.gens()

def es_lisa(av):
    F = av**2 * sum(g**3 for g in xx) - (sum(xx))**3
    I = Rx4.ideal([F] + [F.derivative(g) for g in xx])
    return I.dimension() <= 0

for av in [1, 2, 3, 4, 5, 6, 7, 9]:
    lisa = es_lisa(av)
    if av in [2, 4]:
        chk("control a=" + str(av) + ": la superficie ES SINGULAR", not lisa)
    else:
        chk("control a=" + str(av) + ": la superficie es lisa", lisa)

def evalua(av, pt):
    F = av**2 * sum(g**3 for g in xx) - (sum(xx))**3
    d = dict(zip(xx, [Rx4(c) for c in pt]))
    vals = [F.subs(d)] + [F.derivative(g).subs(d) for g in xx]
    return all(w == 0 for w in vals)

chk("control a=2: (1:1:1:-1) es punto singular (nodo)", evalua(2, [1, 1, 1, -1]))
chk("control a=2: tambien (1:1:-1:1), (1:-1:1:1), (-1:1:1:1) -> 4 nodos (cubica de Cayley)",
    evalua(2, [1, 1, -1, 1]) and evalua(2, [1, -1, 1, 1]) and evalua(2, [-1, 1, 1, 1]))
chk("control a=4: (1:1:1:1) es punto singular", evalua(4, [1, 1, 1, 1]))
chk("control a=2: el coeficiente lider (a^2-4) de la cuadratica se anula -> NO hay dos raices r+-",
    (2**2 - 4) == 0)
info("a=2: A_2(r) = " + str(factor(A_decl.subs({a: R(2)}))) + "  -> solo r=0 y r=-1")
info("a=2: la conica en r=0 es (A_2(0), 3a^2, 0) = (0, 12, 0): RANGO 1, recta doble u=0; idem en r=inf")
chk("control a=2: A_2(0) = a^2-4 = 0, fibra de rango 1 (recta doble) en r=0",
    A_decl.subs({a: R(2), r: R(0)}) == R(0))
chk("control a=4: el discriminante 3a^2(16-a^2) se anula -> r+ = r- = 1 (raiz doble)",
    (3 * 4**2 * (16 - 4**2)) == 0 and ((4**2 - 4) * 1 - (4**2 + 8) + (4**2 - 4)) == 0)
print("      a=2 y a=4: la construccion se rompe en el haz, no en la aritmetica. Se mide aparte si")
print("      la FORMULA de recuento sigue cuadrando por accidente:")
for av in [2, 4]:
    D = -3 * (av**2 - 4)
    E = (av**2 - 4) * (av**2 - 16)
    for p in [7, 11, 13, 17, 19]:
        n = cuenta_P3(av, p)
        f = p**2 + (3 + 3 * chi(D, p) + chi(E, p)) * p + 1
        info("a=" + str(av) + " p=" + str(p) + ": bruto=" + str(n) + " formula(con chi(0)=0)=" + str(f) + "  " + ("coincide" if n == f else "NO coincide"))

print("--- control (extra): si p | 6a(a^2-4)(a^2-16) la formula puede fallar ---")
fallos_malos = 0
casos_malos = 0
for av in AS:
    for p in malos_primos(av) + [3]:
        if p == 3 and av % 3 == 0:
            continue
        n = cuenta_P3(av, p)
        D = -3 * (av**2 - 4)
        E = (av**2 - 4) * (av**2 - 16)
        f = p**2 + (3 + 3 * chi(D, p) + chi(E, p)) * p + 1
        casos_malos += 1
        if n != f:
            fallos_malos += 1
        info("a=" + str(av) + " p=" + str(p) + " (EXCLUIDO): bruto=" + str(n) + " formula=" + str(f) + "  " + ("coincide" if n == f else "NO coincide"))
info("primos excluidos en los que la formula NO cuadra: " + str(fallos_malos) + " de " + str(casos_malos) + " (informativo: la derivacion no afirma nada ahi)")

print("")
print("comprobaciones: " + str(OK[0] + FALLOS[0]) + " ; fallos: " + str(FALLOS[0]))
if FALLOS[0] == 0 and OK[0] >= 200:
    print("VEREDICTO: la derivacion (pasos 1-6) SE SOSTIENE. Verificada simbolicamente sobre Q y")
    print("VEREDICTO: por recuento real fibra a fibra sobre F_p, con a = 2 y a = 4 rotos como se dice.")
else:
    print("VEREDICTO: NO se sostiene, o el numero de comprobaciones es insuficiente. Ver los FALLO.")

print("fin")
