# -*- coding: utf-8 -*-
# ¿DIBUJA LA FIGURA 2 LA REGLA EQUIVOCADA EN EL PANEL IMPAR?   16 de agosto de 2026.
#
# Su vuelta 25 dice que el panel t=5 usa la regla SIMPLECTICA vieja (a_j = eta_j + rk - j + 1) en vez
# de la ortogonal correcta (A_j = 2 eta_j + 2(m'-j) + 1), que por casualidad da los MISMOS 49
# supervivientes pero difiere en 66 puntos, 33 falsos en cada direccion.  Se comprueba.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python check_fig2.py
t, rk, TOPE = 5, 2, 15


def vive_C(e1, e2):
    a = [e1 + rk - 0, e2 + rk - 1]            # a_j = eta_j + rk - j + 1  con j = 1,2
    c = [x % t for x in a]
    if 0 in c:
        return False
    cl = [min(x, t - x) for x in c]
    return len(set(cl)) == rk


def vive_B(e1, e2):
    A = [2 * e1 + 2 * (rk - 1) + 1, 2 * e2 + 2 * (rk - 2) + 1]   # 2 eta_j + 2(m'-j) + 1
    c = [x % t for x in A]
    if 0 in c:
        return False
    cl = [min(x, t - x) for x in c]
    return len(set(cl)) == rk


pts = [(a, b) for a in range(TOPE + 1) for b in range(a + 1)]
vc = [p for p in pts if vive_C(*p)]
vb = [p for p in pts if vive_B(*p)]
sc, sb = set(vc), set(vb)
print("parejas dominantes en la caja eta_1 <= %d : %d" % (TOPE, len(pts)))
print("sobreviven con la regla C (la dibujada) : %d" % len(vc))
print("sobreviven con la regla B (la correcta) : %d" % len(vb))
print("puntos donde las dos reglas DIFIEREN     : %d" % len(sc ^ sb))
print("   falsos supervivientes de C (viven en C y no en B) : %d" % len(sc - sb))
print("   supervivientes que C se pierde  (viven en B y no en C) : %d" % len(sb - sc))
print("")
print("los seis primeros de cada lado:")
print("   C y no B :", sorted(sc - sb)[:6])
print("   B y no C :", sorted(sb - sc)[:6])
