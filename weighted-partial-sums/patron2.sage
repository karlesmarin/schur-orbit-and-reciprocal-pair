# -*- coding: utf-8 -*-
# EL PATRON, PUESTO A PRUEBA.
#
# Del barrido anterior sale una sospecha: para el caracter cuadratico IMPAR (f = 3 mod 4), la suma
# parcial hasta la mitad del periodo se anula justo cuando f = 7 (mod 8), o sea cuando 2 es resto
# cuadratico:
#
#     T((f-1)/2, f) = sum_{c=1}^{(f-1)/2} c (c|f)  =  0   <=>   (2|f) = +1.
#
# Tiene por que: en la formula clasica del numero de clases para f = 3 (mod 4) aparece el factor
# (2 - (2|f)), que vale 1 cuando 2 es resto y 3 cuando no lo es.  Se comprueba la equivalencia en
# las cuatro celdas, y ademas se mira el cero TRIVIAL del caso par: con f = 1 (mod 4) el caracter es
# par y la suma sobre el periodo entero se anula sola, emparejando c con f-c.
#
# Authors: Carles Marin, Claude (AI assistant).

def T(m, f):
    return sum(c * kronecker_symbol(c, f) for c in range(1, m + 1))


print("EQUIVALENCIA:  f = 3 (mod 4) squarefree.   T((f-1)/2, f) = 0  <=>  f = 7 (mod 8) ?")
tab = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
raros = []
for f in range(3, 400):
    if f % 4 != 3 or not is_squarefree(f):
        continue
    cero = (T((f - 1) // 2, f) == 0)
    siete = (f % 8 == 7)
    tab[(cero, siete)] += 1
    if cero != siete:
        raros.append((f, T((f - 1) // 2, f), kronecker_symbol(2, f)))
print("   cero/7mod8   %3d      cero/3mod8   %3d" % (tab[(True, True)], tab[(True, False)]))
print("   no cero/7    %3d      no cero/3    %3d" % (tab[(False, True)], tab[(False, False)]))
if not raros:
    print("   *** EQUIVALENCIA EXACTA sobre %d conductores ***" % sum(tab.values()))
else:
    print("   %d discrepancias: %s" % (len(raros), raros[:8]))

print()
print("EL CERO TRIVIAL:  f = 1 (mod 4) => caracter PAR => T(f-1, f) = 0 siempre ?")
mal = 0
n = 0
for f in range(5, 300):
    if f % 4 != 1 or not is_squarefree(f):
        continue
    n += 1
    if T(f - 1, f) != 0:
        mal += 1
print("   conductores: %d | fallos: %d  (emparejar c con f-c da f*sum(chi) = 0)" % (n, mal))

print()
print("¿Y CUANTO VALE cuando no se anula?  f = 3 (mod 8), con h(-f) al lado")
print("%-6s %-14s %-10s %-12s %s" % ("f", "T((f-1)/2)", "h(-f)", "T / h(-f)", "(2|f)"))
for f in [3, 11, 19, 43, 59, 67, 83, 107, 131, 139]:
    if f % 8 != 3:
        continue
    t = T((f - 1) // 2, f)
    h = QuadraticField(-f).class_number()
    print("%-6d %-14d %-10d %-12s %d" % (f, t, h, t / h if h else "-", kronecker_symbol(2, f)))
