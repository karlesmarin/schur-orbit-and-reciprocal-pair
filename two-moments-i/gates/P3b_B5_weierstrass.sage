# -*- coding: utf-8 -*-
# P3b_B5_weierstrass.sage -- la cubica plana de lambda = (2,1,1,1) ES la curva de conductor 20.
#
# LO QUE HACE FALTA PARA QUE C(5,p) SEA UN TEOREMA Y NO UNA MEDIDA.  El tamiz da
#     C(5,p) = (p-1)[p^2 - 4p + 61 + chi(5)(p+15) + 20 chi(-15) + 10 a_p] / 120,
# y de las siete particiones de 5 seis salen con algebra elemental.  La septima, lambda =
# (2,1,1,1), da la cubica plana
#     C :  2x(x+y+z)^2 + yz(y+z) = 0    en P^2,
# lisa, y por tanto #V = p + 1 - a_p POR DEFINICION de a_p.  Eso no dice nada hasta identificar
# la curva.  Un acuerdo de a_p con eta(2z)^2 eta(10z)^2 en 33 primos no es una demostracion; es
# un control.
#
# LO QUE HACE ESTE GUION.  Da el ISOMORFISMO EXPLICITO sobre Q, que si lo es:
#   (1) C tiene el punto racional (0:1:0), luego es una curva eliptica sobre Q;
#   (2) Sage construye el modelo de Weierstrass y la aplicacion biracional, ambos sobre Q;
#   (3) se comprueba que el resultado es 20.a (y^2 = x^3 + x^2 - x), conductor 20;
#   (4) y se contrasta el a_p del modelo contra el recuento DIRECTO de C sobre F_p, que es lo
#       que la formula del tamiz usa de verdad.
#
# Uso:  sage P3b_B5_weierstrass.sage > P3b_B5_weierstrass_OUT.txt 2>&1
#       (desde el directorio gates/; pide SageMath)
#
# Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
fallos = 0
hechas = 0
def afirma(ok, que):
    global fallos, hechas
    hechas += 1
    if not ok:
        fallos += 1
    print("   %-4s %s" % ("ok" if ok else "MAL", que))
    sys.stdout.flush()
print("(1) LA CUBICA DE lambda = (2,1,1,1), Y QUE ES LISA SOBRE Q")
R = PolynomialRing(QQ, 3, "x,y,z")
x, y, z = R.gens()
F = 2*x*(x+y+z)^2 + y*z*(y+z)
print("   C : %s" % F)
J = R.ideal([F.derivative(v) for v in (x, y, z)])
afirma(J.dimension() == 0, "el ideal jacobiano es cero-dimensional: C es lisa sobre Qbar")
afirma(F.subs(x=0, y=1, z=0) == 0, "y (0:1:0) es un punto racional de C")
print("")
print("(2) EL MODELO DE WEIERSTRASS Y LA APLICACION BIRACIONAL, SOBRE Q")
E = EllipticCurve_from_cubic(F, [0, 1, 0], morphism=False)
phi = EllipticCurve_from_cubic(F, [0, 1, 0], morphism=True)
Emin = E.minimal_model()
print("   modelo de Sage   : %s" % E)
print("   modelo minimal   : %s" % Emin)
print("   conductor        : %s" % Emin.conductor())
print("   discriminante min: %s" % Emin.discriminant())
print("   j-invariante     : %s" % Emin.j_invariant())
afirma(Emin.conductor() == 20, "el conductor es 20")
Eref = EllipticCurve([0, 1, 0, -1, 0])
print("   la de la nota    : %s   conductor %s" % (Eref, Eref.conductor()))
afirma(Emin.is_isomorphic(Eref), "y es ISOMORFA sobre Q a y^2 = x^3 + x^2 - x, la de la nota")
try:
    print("   etiqueta LMFDB/Cremona (si hay base local): %s" % Emin.cremona_label())
except Exception as e:
    print("   etiqueta Cremona no disponible sin base local (%s)" % type(e).__name__)
print("")
print("   la aplicacion biracional C -> E, que es lo que hace de esto una DEMOSTRACION:")
print("   %s" % phi)
print("")
print("(3) EL CONTROL QUE IMPORTA: a_p del modelo contra el RECUENTO DIRECTO de C sobre F_p")
print("    (no contra otra formula nuestra: contra enumerar los puntos proyectivos de C)")
print("   %-6s %-10s %-10s %-10s %s" % ("p", "#C(F_p)", "p+1-#C", "a_p(E)", "ok"))
malos = []
for p in prime_range(7, 152):
    if p == 5:
        continue
    Fp = GF(p)
    cnt = 0
    for pt in ProjectiveSpace(Fp, 2).rational_points():
        if F.change_ring(Fp)(pt[0], pt[1], pt[2]) == 0:
            cnt += 1
    ap_medido = p + 1 - cnt
    ap_curva = Emin.ap(p)
    ok = (ap_medido == ap_curva)
    if not ok:
        malos.append(p)
    if p < 60 or not ok:
        print("   %-6d %-10d %-10d %-10d %s" % (p, cnt, ap_medido, ap_curva, "si" if ok else "NO"))
    sys.stdout.flush()
afirma(not malos, "33 primos 7 <= p <= 151: el recuento directo de C da el a_p de 20.a. Fallos %s" % (malos if malos else "0"))
print("")
print("(4) Y EL PRODUCTO ETA, que es el control SIN RED que viaja con el instrumento")
L = 160
S = PowerSeriesRing(ZZ, "q", default_prec=L)
q = S.gen()
eta = q * prod([(1 - q^(2*n))^2 * (1 - q^(10*n))^2 for n in range(1, L//2 + 2)])
coef = eta.list()
malos = [p for p in prime_range(7, 152) if p != 5 and coef[p] != Emin.ap(p)]
afirma(not malos, "eta(2z)^2 eta(10z)^2 reproduce los a_p: la newform de peso 2 y nivel 20. Fallos %s" % (malos if malos else "0"))
print("")
print("=" * 88)
print("comprobaciones: %d ; fallos: %d" % (hechas, fallos))
if hechas < 6:
    print("FATAL: se esperaban 6 comprobaciones y hay %d. El guion se rompio por el camino y el" % hechas)
    print("exit code NO lo dice: 'sage < fichero' devuelve 0 aunque cada celda lance una excepcion.")
print("VEREDICTO: %s" % ("la cubica de (2,1,1,1) ES 20.a sobre Q, con isomorfismo explicito" if (fallos == 0 and hechas >= 6) else "NO CUADRA"))
