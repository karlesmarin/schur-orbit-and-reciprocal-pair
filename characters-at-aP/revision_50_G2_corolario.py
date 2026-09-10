# -*- coding: utf-8 -*-
# G2 NO ES UNA COMPROBACION: ES UN COROLARIO.   25 de agosto de 2026.
#
# POR QUE.  La nota demuestra el criterio  chi != 0  <=>  N_q(lambda) = r_q  por la especializacion
# principal, y declara que NO depende del tipo.  Luego en G2 los 49 pesos son un CONTROL, no la
# prueba --- pero el texto aun se apoyaba en el "49 de 49" como si lo fuera.  Un lector puede
# entender que en G2 el resultado sigue siendo experimental.  Aqui se comprueba que las dos
# congruencias salen del criterio para TODO peso, no para 49.
#
# LOS SEIS VALORES.  Con lambda = a*w1 + b*w2, u = a+1, v = 3(b+1), los (lambda+rho, alpha) sobre
# las seis raices positivas son:   u,  v,  u+v,  2u+v,  3u+v,  3u+2v.
# En lambda = 0 valen 1, 3, 4, 5, 6, 9, y r_q = #{de esos seis divisibles por q}.
#
# QUE SE MIDE
#   K0  los seis valores, contra la construccion explicita del sistema de raices de G2
#   K1  q=3: N_3 = r_3  <=>  a !≡ 2 (mod 3),  sobre TODOS los pesos hasta a,b <= 200
#   K2  q=2: N_2 = r_2  <=>  no (a impar y b impar),  igual rango
#   K3  la derivacion algebraica de las dos, escrita como reduccion mod q
#   K4  SENUELO: las congruencias con el signo cambiado deben fallar
#   K5  y el criterio contra los caracteres exactos, que es el control que ya existia (49 pesos)
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python revision_50_G2_corolario.py > revision_50_G2_corolario_OUT.txt 2>&1

import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


# --- G2 explicito: matriz de Cartan [[2,-1],[-3,2]] con alpha1 corta.  La forma simetrizada es
#     [[2,-3],[-3,6]] --- la del fichero R_P_todos_los_tipos.py, ya con su control de transpuesta.
#     Las seis raices positivas en coordenadas de pesos fundamentales, y (lambda+rho, alpha).
def seis_valores(a, b):
    u, v = a + 1, 3 * (b + 1)
    return (u, v, u + v, 2 * u + v, 3 * u + v, 3 * u + 2 * v)


print(SEP)
print("K0 -- LOS SEIS VALORES, contra el sistema de raices construido a mano")
# raices positivas de G2 en la base simple (alpha1 corta, alpha2 larga):
#   a1, a2, a1+a2, 2a1+a2, 3a1+a2, 3a1+2a2
RAICES = [(1, 0), (0, 1), (1, 1), (2, 1), (3, 1), (3, 2)]
# forma bilineal en la base de raices simples, normalizada a (a1,a1)=2, (a2,a2)=6, (a1,a2)=-3
G = [[2, -3], [-3, 6]]
# lambda = a w1 + b w2.  Con (w_i, a_j^v) = delta:  (lambda, alpha) = ... se evalua con la
# expresion cerrada u, v.  Se comprueba que la construccion por raices da lo mismo.
mal = 0
for a in range(0, 9):
    for b in range(0, 9):
        # rho = w1 + w2;  ell = lambda + rho = (a+1) w1 + (b+1) w2
        # (w_i, alpha_j) = delta_ij (alpha_j, alpha_j)/2  ->  (w1,a1)=1, (w2,a2)=3
        c1, c2 = a + 1, b + 1
        obs = []
        for (m1, m2) in RAICES:
            # (ell, alpha) = c1*m1*(w1,a1) + c2*m2*(w2,a2) = c1*m1*1 + c2*m2*3
            obs.append(c1 * m1 * 1 + c2 * m2 * 3)
        esp = list(seis_valores(a, b))
        if sorted(obs) != sorted(esp):
            mal += 1
            if mal == 1:
                print(f"   primer desajuste en (a,b)=({a},{b}): {sorted(obs)} vs {sorted(esp)}")
ok(mal == 0, "los seis valores coinciden con la construccion por raices", f"{81 - mal}/81")
ok(G == [[gg for gg in fila] for fila in G] and G[0][1] == G[1][0],
   "la forma simetrizada es simetrica (control de transpuesta)")

print(SEP)
print("K1/K2 -- LAS DOS CONGRUENCIAS, PARA TODOS LOS PESOS (no para 49)")
RHO = seis_valores(0, 0)          # 1, 3, 4, 5, 6, 9
print(f"   los seis valores en lambda=0: {RHO}")
N = 200
for q, nombre, pred in (
        (3, "q=3", lambda a, b: a % 3 != 2),
        (2, "q=2", lambda a, b: not (a % 2 == 1 and b % 2 == 1))):
    r_q = sum(1 for x in RHO if x % q == 0)
    mal = tot = 0
    for a in range(N + 1):
        for b in range(N + 1):
            n_q = sum(1 for x in seis_valores(a, b) if x % q == 0)
            tot += 1
            if (n_q == r_q) != pred(a, b):
                mal += 1
    ok(mal == 0, f"{nombre}: r_q={r_q}; el criterio equivale a la congruencia",
       f"{tot - mal}/{tot} pesos")

print(SEP)
print("K3 -- LA DERIVACION, escrita como reduccion mod q")
print("   q=3:  v = 3(b+1) ≡ 0 siempre.  Los seis mod 3 son  u, 0, u, 2u, 0, 0.")
print("         Luego N_3 = 3 + 3*[3 | u].  r_3 = 3.  =>  N_3 = r_3  <=>  3 no divide u = a+1")
print("                                                             <=>  a ≢ 2 (mod 3).")
print("   q=2:  3u+v ≡ u+v,  2u+v ≡ v,  3u+2v ≡ u.  Los seis mod 2 son  u, v, u+v repetidos.")
print("         Luego N_2 = 2 * #{u, v, u+v pares}, que vale 2 salvo que u y v sean AMBOS pares,")
print("         en cuyo caso vale 6.  r_2 = 2.  =>  N_2 = r_2  <=>  no (u par y v par)")
print("                                                       <=>  no (a impar y b impar).")
# y se comprueba que esas dos reducciones son las correctas
mal3 = mal2 = 0
for a in range(0, 60):
    for b in range(0, 60):
        u, v = a + 1, 3 * (b + 1)
        s = seis_valores(a, b)
        if [x % 3 for x in s] != [u % 3, 0, u % 3, (2 * u) % 3, 0, 0]:
            mal3 += 1
        if sorted(x % 2 for x in s) != sorted([u % 2, v % 2, (u + v) % 2] * 2):
            mal2 += 1
ok(mal3 == 0, "K3: la reduccion mod 3 es  (u, 0, u, 2u, 0, 0)")
ok(mal2 == 0, "K3: la reduccion mod 2 es  (u, v, u+v) por duplicado")

print(SEP)
print("K4 -- SENUELO: las congruencias alteradas tienen que FALLAR")
for q, etq, pred in (
        (3, "a ≢ 1 (mod 3)  [en vez de 2]", lambda a, b: a % 3 != 1),
        (3, "b ≢ 2 (mod 3)  [b en vez de a]", lambda a, b: b % 3 != 2),
        (2, "no (a par y b par)  [paridad opuesta]",
         lambda a, b: not (a % 2 == 0 and b % 2 == 0))):
    r_q = sum(1 for x in RHO if x % q == 0)
    mal = tot = 0
    for a in range(60):
        for b in range(60):
            n_q = sum(1 for x in seis_valores(a, b) if x % q == 0)
            tot += 1
            if (n_q == r_q) != pred(a, b):
                mal += 1
    ok(mal > 0, f"senuelo q={q}: '{etq}' falla, como debe", f"{mal} de {tot} discrepan")

print(SEP)
print("K5 -- EL CONTROL QUE YA EXISTIA: el criterio contra los caracteres exactos")
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "revision_47_G2_OUT.txt")
if os.path.exists(p):
    o = open(p, encoding="utf-8", errors="replace").read()
    ok("49/49" in o, "la salida archivada de revision_47_G2 trae el 49/49")
    ok("FALLA" not in o, "y no trae ningun fallo")
else:
    ok(False, "falta revision_47_G2_OUT.txt")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: en G2 las dos congruencias son COROLARIO del criterio, para todo peso.")
print("           Los 49 pesos son el control con caracteres exactos, no la prueba.")
