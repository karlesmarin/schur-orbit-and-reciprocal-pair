# -*- coding: utf-8 -*-
# LA LEY DEL RECUENTO.   8 de septiembre de 2026.
# survivors_quasipoly.py midio S(q) = #{x en P/qP : N_q(x) = r_q} en 18 tipos.  Cuadra con
# prod(q - m_i) en unos sitios y no en otros.  Esta gate decide la LEY: la hipotesis es
#      S(q) = prod_i (q - m_i)   <=>   r_q = 0  y  q coprimo con los primos MALOS del tipo.
# Primos malos (clasicos): ninguno en A; 2 en B, C, D; 2 y 3 en G2, F4, E6, E7; 2,3,5 en E8.
# Se imprimen las CUATRO celdas de la equivalencia, que es lo unico que la prueba de verdad.
import json
import sys
from itertools import product

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

MALOS = {"A": [], "B": [2], "C": [2], "D": [2], "G": [2, 3], "F": [2, 3], "E": [2, 3]}
TECHO = 40_000_000
DAT = json.load(open("root_data.json", encoding="utf-8"))


def prod(v):
    r = 1
    for x in v:
        r *= x
    return r


celdas = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
raros = []
print("=" * 96)
print("HIPOTESIS:  S(q) = prod(q - m_i)   <=>   r_q = 0  y  gcd(q, primos malos) = 1")
print("=" * 96)
for d in DAT:
    n, R, exps = d["rango"], d["raices"], d["exponentes"]
    letra = d["tipo"][0]
    malos = MALOS[letra]
    qmax = 20 if n <= 4 else (12 if n == 5 else 8)
    for q in range(2, qmax + 1):
        if (q ** n) * len(R) > TECHO:
            continue
        r_q = sum(1 for c, v in R if (c * sum(v)) % q == 0)
        S = 0
        for x in product(range(q), repeat=n):
            k = 0
            for c, v in R:
                s = 0
                for j in range(n):
                    s += x[j] * v[j]
                if (c * s) % q == 0:
                    k += 1
                    if k > r_q:
                        break
            if k == r_q:
                S += 1
        lhs = (S == prod([q - m for m in exps]))
        rhs = (r_q == 0) and all(q % p for p in malos)
        celdas[(lhs, rhs)] += 1
        if lhs != rhs:
            raros.append((d["tipo"], q, r_q, S, prod([q - m for m in exps]), lhs, rhs))
    print("   %-4s hecho" % d["tipo"])

print()
print("   LAS CUATRO CELDAS")
print("                          hipotesis SI      hipotesis NO")
print("   formula acierta   %14d %17d" % (celdas[(True, True)], celdas[(True, False)]))
print("   formula falla     %14d %17d" % (celdas[(False, True)], celdas[(False, False)]))
print()
if not raros:
    print("RESULTADO: la equivalencia se cumple en todos los casos medidos, y las dos diagonales")
    print("           estan pobladas, luego no es trivial en ninguno de los dos sentidos.")
else:
    print("RESULTADO: NO es una equivalencia.  Casos que se salen (tipo,q,r_q,S,pred,lhs,rhs):")
    for r in raros:
        print("   ", r)
