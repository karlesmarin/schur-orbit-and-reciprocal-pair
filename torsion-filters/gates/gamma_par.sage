# -*- coding: utf-8 -*-
# .SE CIERRA gamma_t TAMBIEN EN LA RAMA PAR?   16 de agosto de 2026.  (vuelta 30, punto 7)
#
# cor:galoisquad cierra gamma_t para t impar con  D_t = det(zeta^{ij} - zeta^{-ij})_{1<=i,j<=n},
# n = (t-1)/2,  via  D_t^2 = (-t)^n.  Se observa que el mismo determinante parece dar
# lo mismo en la rama C con  n = t/2 - 1,  y que entonces el enunciado uniforme seria
#
#      gamma_t = 1  si n es par,      gamma_t = caracter de Q(sqrt(-t))  si n es impar,
#
#      con  n = (t-1)/2  (B, t impar)  y  n = t/2 - 1  (C, t par),
#
# o sea: trivial para  t = 1, 2 (mod 4)  y cuadratico para  t = 3, 0 (mod 4).
#
# LO QUE SE MIDE, en la rama PAR
#   P1  FATAL: D_t^2 == (-t)^n  exacto.
#   P2  sigma_k(D_t) == gamma_t(k) D_t  para toda unidad.
#   P3  gamma_t trivial  <=>  n par.
#   P4  y la clasificacion final por t mod 4, juntando las dos ramas.
#
# CONTROL
#   C0  P1 es fatal.  Si D^2 no es (-t)^n, el argumento del subcuerpo no existe en el par y el
#       enunciado uniforme se cae.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage gamma_par.sage

import json
import sys


def plegar(v, t):
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def sgn_perm(perm):
    n, s, visto = len(perm), 1, [False] * len(perm)
    for i in range(n):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            s = -s
    return int(s)


def gamma_de(k, t, n):
    cl, ep = [], []
    for j in range(1, n + 1):
        c_, e_ = plegar(k * j, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, n + 1)):
        return None
    s = sgn_perm([c - 1 for c in cl])
    for e in ep:
        s *= e
    return int(s)


print("=" * 104)
print("gamma_t EN LAS DOS RAMAS:  n = (t-1)/2 en B,  n = t/2 - 1 en C")
print("=" * 104)
print("")
print("   t  rama   n   n par   D^2 == (-t)^n   sigma_k(D)=gamma D   gamma trivial   t mod 4")
print("   " + "-" * 90)

RES = []
p1 = p1n = p2 = p2n = p3 = p3n = 0
for t in list(range(3, 24)):
    if t < 4 and t % 2 == 0:
        continue
    impar = (t % 2 == 1)
    n = (t - 1) // 2 if impar else t // 2 - 1
    if n < 1:
        continue
    K = CyclotomicField(t)
    z = K.gen()
    M = matrix(K, n, n, lambda i, j: z ** ((i + 1) * (j + 1)) - z ** (-(i + 1) * (j + 1)))
    D = M.determinant()
    ok1 = (D ** 2 == K((-t) ** n))
    unidades = [k for k in range(1, t) if gcd(k, t) == 1]
    ok2 = 0
    triv = True
    for k in unidades:
        g = gamma_de(k, t, n)
        if g is None:
            continue
        sig = K.hom([z ** k])
        ok2 += 1 if sig(D) == g * D else 0
        if g != 1:
            triv = False
    ok3 = (triv == (n % 2 == 0))
    p1 += 1 if ok1 else 0
    p1n += 1
    p2 += ok2
    p2n += len(unidades)
    p3 += 1 if ok3 else 0
    p3n += 1
    print("  %2d   %s   %3d   %-5s   %-5s           %2d/%2d              %-5s        %d"
          % (t, "B" if impar else "C", n, str(n % 2 == 0), str(ok1), ok2, len(unidades),
             str(triv), t % 4))
    sys.stdout.flush()
    RES.append({"t": int(t), "rama": "B" if impar else "C", "n": int(n),
                "n_par": bool(n % 2 == 0), "D2": bool(ok1),
                "sigma": [int(ok2), int(len(unidades))], "trivial": bool(triv),
                "t_mod_4": int(t % 4), "P3": bool(ok3)})

print("")
print("  P1  D^2 == (-t)^n            : %d de %d   <== FATAL" % (p1, p1n))
print("  P2  sigma_k(D) = gamma_t(k) D: %d de %d" % (p2, p2n))
print("  P3  gamma trivial <=> n par  : %d de %d" % (p3, p3n))
print("")
print("  clasificacion final por t mod 4:")
for m in (0, 1, 2, 3):
    sub = [x for x in RES if x["t_mod_4"] == m]
    if not sub:
        continue
    tri = set(x["trivial"] for x in sub)
    print("     t = %d (mod 4) : t = %-26s  gamma trivial = %s"
          % (m, str([x["t"] for x in sub])[:26], sorted(tri)))

json.dump(RES, open("gamma_par_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
