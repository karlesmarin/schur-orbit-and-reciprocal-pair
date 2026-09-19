# -*- coding: utf-8 -*-
# ¿HAY PATRON EN LAS CASILLAS DONDE MUERE EL CARACTER CUADRATICO?
#
# La singularidad sobre Q es la anulacion de la componente de un caracter real e IMPAR.  Para el
# segmento centrado {-m..m} y chi el caracter cuadratico de conductor f | q', eso es
#
#     T(m, f) = sum_{c=1}^{m} c * chi(c) = 0,
#
# una suma parcial de caracteres.  Aqui se calcula T sobre una rejilla ancha y se buscan
# estructuras: donde se anula, con que frecuencia, si depende de m modulo algo, si el q' esta
# ligado a m, y si el conductor que muere es primo o compuesto.
#
# Se prueban ademas dos lecturas clasicas:
#   L1  chi impar exige f = 3 (mod 4).  ¿Se anula alguna vez con f = 1 (mod 4)?
#   L2  la suma COMPLETA sum_{c<f} c chi(c) es -f h(-f) por la formula de clase; ¿la parcial en
#       m = (f-1)/2 tiene que ver con h(-f)?
#
# Authors: Carles Marin, Claude (AI assistant).

def T(m, f):
    return sum(c * kronecker_symbol(c, f) for c in range(1, m + 1))


print("=" * 90)
print("L1  ¿solo con f = 3 (mod 4)?  (el cuadratico es impar solo entonces)")
print("=" * 90)
cuatro1 = cuatro3 = 0
ceros = []
for f in range(3, 200):
    if not is_squarefree(f) or f % 2 == 0:
        continue
    for m in range(1, min(f, 60)):
        if T(m, f) == 0:
            ceros.append((m, f))
            if f % 4 == 1:
                cuatro1 += 1
            else:
                cuatro3 += 1
print("   ceros hallados: %d | con f = 1 (mod 4): %d | con f = 3 (mod 4): %d" % (
    len(ceros), cuatro1, cuatro3))

print()
print("=" * 90)
print("LOS CEROS, por f")
print("=" * 90)
porf = {}
for (m, f) in ceros:
    porf.setdefault(f, []).append(m)
for f in sorted(porf):
    print("   f=%-4d (%s, %s) : m = %s" % (
        f, "1 mod 4" if f % 4 == 1 else "3 mod 4",
        "primo" if is_prime(f) else "compuesto", porf[f]))

print()
print("=" * 90)
print("¿Y la m mas pequena de cada f?  ¿guarda relacion con f?")
print("=" * 90)
print("%-6s %-9s %-8s %-10s %-10s %s" % ("f", "m minimo", "m/f", "h(-f)", "(2|f)", "(3|f)"))
for f in sorted(porf)[:24]:
    m0 = min(porf[f])
    try:
        h = QuadraticField(-f).class_number() if f % 4 == 3 else QuadraticField(-f).class_number()
    except Exception:
        h = "?"
    print("%-6d %-9d %-8.3f %-10s %-10d %d" % (
        f, m0, float(m0) / f, h, kronecker_symbol(2, f), kronecker_symbol(3, f)))

print()
print("=" * 90)
print("L2  la suma COMPLETA hasta (f-1)/2 y el numero de clases")
print("=" * 90)
print("%-6s %-12s %-12s %-10s %s" % ("f", "T((f-1)/2)", "sum (c|f)", "h(-f)", "se anula T en algun m?"))
for f in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
    med = (f - 1) // 2
    s1 = T(med, f)
    s2 = sum(kronecker_symbol(c, f) for c in range(1, med + 1))
    try:
        h = QuadraticField(-f).class_number()
    except Exception:
        h = "?"
    print("%-6d %-12d %-12d %-10s %s" % (f, s1, s2, h, sorted(porf.get(f, []))[:6]))
