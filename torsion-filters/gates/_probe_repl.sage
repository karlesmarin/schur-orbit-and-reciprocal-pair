# -*- coding: utf-8 -*-
# .QUE MATA EXACTAMENTE UN .sage ALIMENTADO AL REPL?   19 de agosto de 2026.
# Tres casos minimos, cada uno aislado.  Sin adivinar: se mira que imprime cada uno.
import sys
print("A_INICIO")
sys.stdout.flush()
x = 0
if x == 1:
    print("A_IF")
else:
    print("A_ELSE_SI_SALE_EL_ELSE_SOBREVIVE")
sys.stdout.flush()
print("B_INICIO")
def f_con_blanco(n):
    y = n + 1

    return y * 10
print("B_f_con_blanco(1) =", f_con_blanco(1), " <-- None significa que el blanco corta el cuerpo")
sys.stdout.flush()
print("C_INICIO")
tot = 0
for k in range(3):
    tot += k

    tot += 100
print("C_tot =", tot, " <-- 3 si el blanco corto el for tras la primera linea; 303 si no")
sys.stdout.flush()
print("D_FIN")
