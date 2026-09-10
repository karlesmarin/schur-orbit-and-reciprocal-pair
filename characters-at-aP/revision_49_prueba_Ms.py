# -*- coding: utf-8 -*-
# LA DEMOSTRACION DE M_s = M_0, PIEZA A PIEZA.   25 de agosto de 2026.
#
# POR QUE.  La nota asciende 3.2 a teorema pero deja la Prop. 3.3 (M_s = M_0) solo comprobada
# numericamente.  Hay demostracion elemental; se verifica cada paso antes de escribirla.
#
# ⚠ AVISO DE DEFINICION.  M_s es el producto de los MODULOS de los factores que no se anulan.
#    Como numero complejo NO es independiente de s: las fases cambian.  Lo que se necesita, y lo
#    unico que se afirma, es la independencia del modulo.
#
# LOS PASOS
#   S1  prod_{x in U} |a - x| = q,   con U = mu_q \ {a}.   (Es |(X^q-1)'| en X=a.)
#   S2  A_s = prod_{x<y in U} |x - y| = q^{(q-2)/2}.
#       ⚠ CORRECCION: la nota dice "the Vandermonde part alone is q^{q-2}".  Eso es su CUADRADO.
#   S3  C_s = prod ORDENADO sobre (x,y) con xy != 1  =  q^{q-2}|1-a^2| si a^2 != 1,  q^{q-2} si no.
#   S4  D_s = D/|1-a^2| si a^2 != 1,  D si a^2 = 1,  con D = prod_{x in mu_q, x^2 != 1}|1-x^2|.
#   S5  D = q si q impar,  (q/2)^2 si q par.
#   S6  C_s D_s = q^{q-2} D,  y como C_s = B_s^2 D_s se sigue (B_s D_s)^2 = q^{q-2} D.
#   S7  M_s = A_s B_s D_s  =>  M_s = M_0.
#   S8  y el paso que permite mover ell por Weyl:  |1 - z^{-r}| = |1 - z^r| para |z| = 1.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python revision_49_prueba_Ms.py > revision_49_prueba_Ms_OUT.txt 2>&1

import cmath
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []
TOL = 1e-9


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


def mu(q):
    return [cmath.exp(2j * cmath.pi * r / q) for r in range(q)]


def U(q, s):
    return [x for r, x in enumerate(mu(q)) if r != s % q]


print(SEP)
print("S1 -- prod_{x in U} |a - x| = q")
peor = 0.0
for q in range(3, 12):
    for s in range(q):
        a = mu(q)[s % q]
        p = 1.0
        for x in U(q, s):
            p *= abs(a - x)
        peor = max(peor, abs(p - q) / q)
print(f"   peor error relativo sobre q=3..11 y todo s: {peor:.2e}")
ok(peor < 1e-9, "S1: el producto vale q")

print(SEP)
print("S2 -- A_s = prod_{x<y} |x-y| = q^{(q-2)/2}   (la nota decia q^{q-2}: es el CUADRADO)")
print(f"{'q':>3} {'A_0 medido':>16} {'q^{(q-2)/2}':>16} {'q^{q-2}':>16}")
peor = 0.0
for q in range(3, 11):
    Us = U(q, 0)
    A = 1.0
    for i in range(len(Us)):
        for j in range(i + 1, len(Us)):
            A *= abs(Us[i] - Us[j])
    t = float(q) ** ((q - 2) / 2)
    print(f"{q:>3} {A:>16.4f} {t:>16.4f} {float(q)**(q-2):>16.4f}")
    for s in range(q):
        Us = U(q, s)
        As = 1.0
        for i in range(len(Us)):
            for j in range(i + 1, len(Us)):
                As *= abs(Us[i] - Us[j])
        peor = max(peor, abs(As - t) / t)
ok(peor < 1e-8, "S2: A_s = q^{(q-2)/2} para todo s", f"peor {peor:.1e}")

print(SEP)
print("S3/S4/S5/S6 -- C_s (ORDENADO), D_s, y su producto")


def CDB(q, s):
    Us = U(q, s)
    C = 1.0                                  # ordenado: (x,y) con xy != 1
    for x in Us:
        for y in Us:
            if abs(x * y - 1) > TOL:
                C *= abs(1 - x * y)
    B = 1.0                                  # no ordenado, x<y
    for i in range(len(Us)):
        for j in range(i + 1, len(Us)):
            if abs(Us[i] * Us[j] - 1) > TOL:
                B *= abs(1 - Us[i] * Us[j])
    D = 1.0
    for x in Us:
        if abs(x * x - 1) > TOL:
            D *= abs(1 - x * x)
    return C, B, D


def Dfull(q):
    p = 1.0
    for x in mu(q):
        if abs(x * x - 1) > TOL:
            p *= abs(1 - x * x)
    return p


print(f"{'q':>3} {'s':>3} {'C_s medido':>14} {'formula':>14} {'D_s medido':>12} {'formula':>12} "
      f"{'C_s=B^2 D?':>11}")
p3 = p4 = p6 = 0.0
for q in range(3, 10):
    Dq = Dfull(q)
    for s in range(q):
        a = mu(q)[s]
        C, B, D = CDB(q, s)
        a2 = abs(a * a - 1) <= TOL
        Cf = float(q) ** (q - 2) * (1.0 if a2 else abs(1 - a * a))
        Df = Dq if a2 else Dq / abs(1 - a * a)
        p3 = max(p3, abs(C - Cf) / Cf)
        p4 = max(p4, abs(D - Df) / max(Df, 1e-12))
        p6 = max(p6, abs(C - B * B * D) / C)
        if s < 2 and q < 8:
            print(f"{q:>3} {s:>3} {C:>14.4f} {Cf:>14.4f} {D:>12.4f} {Df:>12.4f} "
                  f"{abs(C-B*B*D)/C:>11.1e}")
ok(p3 < 1e-8, "S3: la formula de C_s (ordenado)", f"peor {p3:.1e}")
ok(p4 < 1e-8, "S4: la formula de D_s", f"peor {p4:.1e}")
ok(p6 < 1e-8, "S6a: C_s = B_s^2 D_s", f"peor {p6:.1e}")

print()
print(f"{'q':>3} {'D medido':>12} {'q o (q/2)^2':>14} {'C_s D_s':>16} {'q^{q-2} D':>16}")
p5 = p6b = 0.0
for q in range(3, 10):
    Dq = Dfull(q)
    t = float(q) if q % 2 else (q / 2.0) ** 2
    p5 = max(p5, abs(Dq - t) / t)
    obj = float(q) ** (q - 2) * Dq
    C, B, D = CDB(q, 0)
    for s in range(q):
        C, B, D = CDB(q, s)
        p6b = max(p6b, abs(C * D - obj) / obj)
    print(f"{q:>3} {Dq:>12.4f} {t:>14.4f} {C*D:>16.4f} {obj:>16.4f}")
ok(p5 < 1e-8, "S5: D = q (q impar), (q/2)^2 (q par)", f"peor {p5:.1e}")
ok(p6b < 1e-8, "S6b: C_s D_s = q^{q-2} D, independiente de s", f"peor {p6b:.1e}")

print(SEP)
print("S7 -- M_s = A_s B_s D_s es independiente de s")
peor = 0.0
for q in range(3, 11):
    vals = []
    for s in range(q):
        Us = U(q, s)
        A = 1.0
        for i in range(len(Us)):
            for j in range(i + 1, len(Us)):
                A *= abs(Us[i] - Us[j])
        C, B, D = CDB(q, s)
        vals.append(A * B * D)
    peor = max(peor, max(abs(v / vals[0] - 1) for v in vals))
ok(peor < 1e-8, "S7: M_s = M_0", f"peor {peor:.1e}")

print(SEP)
print("S8 -- |1 - z^{-r}| = |1 - z^r| para |z| = 1   (el paso que permite mover ell por Weyl)")
import random
random.seed(11)
peor = 0.0
for _ in range(200):
    th = random.uniform(0, 2 * cmath.pi)
    z = cmath.exp(1j * th)
    r = random.randint(1, 30)
    peor = max(peor, abs(abs(1 - z ** (-r)) - abs(1 - z ** r)))
ok(peor < 1e-12, "S8: la identidad de modulos", f"peor {peor:.1e}")

print(SEP)
print("S9 -- SENUELO: M_s como numero COMPLEJO (sin modulos) NO es independiente de s")
dep = 0
for q in range(4, 10):
    vals = []
    for s in range(q):
        Us = U(q, s)
        p = 1 + 0j
        for i in range(len(Us)):
            for j in range(i + 1, len(Us)):
                p *= (1 - Us[i] / Us[j])
                if abs(Us[i] * Us[j] - 1) > TOL:
                    p *= (1 - Us[i] * Us[j])
        for x in Us:
            if abs(x * x - 1) > TOL:
                p *= (1 - x * x)
        vals.append(p)
    d = max(abs(v - vals[0]) for v in vals) / max(abs(vals[0]), 1e-12)
    dep += (d > 1e-6)
    print(f"   q={q}: max |M_s - M_0| / |M_0| (complejo) = {d:.2e}")
ok(dep >= 5, "S9: en complejo SI depende de s --- por eso M_s se define con MODULOS",
   f"{dep} de 6")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: los nueve pasos verificados.  La Proposicion 3.3 tiene demostracion.")
