# -*- coding: utf-8 -*-
# LA SUMA PARCIAL ES EL NUMERO DE CLASES.
#
# Sospecha del barrido anterior:  para f = 3 (mod 4) libre de cuadrados,
#
#     T((f-1)/2, f) = sum_{c=1}^{(f-1)/2} c (c|f)  =  f * h(-f) * 2/w(-f),
#
# con w el numero de unidades del cuerpo cuadratico imaginario (w = 6 para f = 3, w = 2 el resto).
# En particular se anula exactamente cuando... NO: se anula cuando f = 7 (mod 8), y entonces la
# formula de arriba no puede valer porque h nunca es cero.  Asi que hay DOS regimenes y hay que
# separarlos, que es justo lo que distingue (2|f) = +1 de -1.
#
# Se comprueba:
#   A  f = 3 (mod 8):  T = f h(-f) 2/w  exactamente?
#   B  f = 7 (mod 8):  T = 0 siempre?  (ya visto en 42 conductores; se amplia)
#   C  ¿y la suma SIN peso, sum (c|f), que es la de la formula clasica de Dirichlet?
#
# Authors: Carles Marin, Claude (AI assistant).

def T(m, f):
    return sum(c * kronecker_symbol(c, f) for c in range(1, m + 1))


def U(m, f):
    return sum(kronecker_symbol(c, f) for c in range(1, m + 1))


print("A  f = 3 (mod 8):  T((f-1)/2) = f h(-f) 2/w ?")
print("%-6s %-10s %-6s %-4s %-14s %-14s %s" % ("f", "T", "h(-f)", "w", "f*h*2/w", "coincide", "sin peso U"))
okA = malA = 0
for f in range(3, 320):
    if f % 8 != 3 or not is_squarefree(f):
        continue
    K = QuadraticField(-f)
    h = K.class_number()
    w = len(K.unit_group().torsion_generator().multiplicative_order() * [0]) if False else \
        K.number_of_roots_of_unity()
    t = T((f - 1) // 2, f)
    pred = f * h * 2 / w
    ok = (t == pred)
    okA += ok
    malA += (not ok)
    if f < 130 or not ok:
        print("%-6d %-10d %-6d %-4d %-14s %-14s %d" % (f, t, h, w, pred, "si" if ok else "*** NO ***",
                                                       U((f - 1) // 2, f)))
print("   aciertos %d, fallos %d" % (okA, malA))

print()
print("B  f = 7 (mod 8):  T((f-1)/2) = 0 ?")
okB = malB = 0
for f in range(7, 400):
    if f % 8 != 7 or not is_squarefree(f):
        continue
    if T((f - 1) // 2, f) == 0:
        okB += 1
    else:
        malB += 1
        print("   *** f=%d da %d" % (f, T((f - 1) // 2, f)))
print("   se anula en %d conductores, falla en %d" % (okB, malB))

print()
print("C  la formula clasica de Dirichlet, para comparar:  h(-f) = U((f-1)/2) / (2 - (2|f))")
print("%-6s %-8s %-8s %-8s %-10s %s" % ("f", "U", "(2|f)", "2-(2|f)", "U/(2-(2|f))", "h(-f)"))
okC = malC = 0
for f in range(7, 200):
    if f % 4 != 3 or not is_squarefree(f) or f == 3:
        continue
    u = U((f - 1) // 2, f)
    k = 2 - kronecker_symbol(2, f)
    h = QuadraticField(-f).class_number()
    ok = (u == k * h)
    okC += ok
    malC += (not ok)
    if f < 60:
        print("%-6d %-8d %-8d %-8d %-10s %d" % (f, u, kronecker_symbol(2, f), k, u / k, h))
print("   la clasica acierta %d, falla %d  (control de que el montaje es el correcto)" % (okC, malC))
