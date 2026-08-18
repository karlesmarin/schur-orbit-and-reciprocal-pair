# comprobacion puntual del testigo: Sp_4, t=6, lambda=(1,0).  Tres rutas.
K = CyclotomicField(6)
xi = K.gen()
m, t = 2, 6
rho = [2, 1]
eta = (1, 0)
a = [eta[j] + rho[j] for j in range(m)]

# R1 Freudenthal
W = WeylCharacterRing("C2")
el = W(W.space().from_vector(vector([Integer(1), Integer(0)])))
s = K(0)
for wt, mult in el.weight_multiplicities().items():
    v = [int(u) for u in wt.to_vector()]
    s += int(mult) * xi ** (sum((i + 1) * v[i] for i in range(m)) % t)
print("R1 Freudenthal      :", s, "   (dim de la irreducible:", el.degree(), ")")

# R2 bialternante
num = matrix(K, m, m, lambda i, j: xi ** ((i + 1) * a[j] % t) - xi ** ((-(i + 1) * a[j]) % t))
den = matrix(K, m, m, lambda i, j: xi ** ((i + 1) * rho[j] % t) - xi ** ((-(i + 1) * rho[j]) % t))
print("R2 bialternante     :", num.determinant() / den.determinant(), "  (denominador:", den.determinant(), ")")

# R3 a mano: la estandar de Sp_4 tiene autovalores xi^{+-1}, xi^{+-2}
print("R3 a mano           :", xi + xi ** 5 + xi ** 2 + xi ** 4)

# el centralizador: ninguna coraiz de C_2 cumple 6 | <a, alpha^v>
L = RootSystem("C2").ambient_space()
print("a = lambda + rho    :", a, "   a mod 6:", [x % t for x in a])
for al in L.positive_roots():
    av = [QQ(u) for u in al.to_vector()]
    p = 2 * sum(QQ(a[i]) * av[i] for i in range(m)) / sum(av[i] * av[i] for i in range(m))
    print("   raiz %-12s coraiz-pairing <a,alpha^v> = %-4s   divisible por 6: %s"
          % (str(av), str(p), p % t == 0))
print("")
print("y la pared que SI lo mata:  2*a_1 =", 2 * a[0], "= 0 mod 6 ->", (2 * a[0]) % t == 0)
print("DONE")
