# -*- coding: utf-8 -*-
# EL TESTIGO MINIMO, POR TRES RUTAS INDEPENDIENTES.   16 de agosto de 2026.
#
# QUE SE AFIRMA.  Para  G = Sp_4  (tipo C_2, numero de Coxeter h = 4) y el elemento
#
#     g = diag(xi, xi^2, xi^{-1}, xi^{-2}),   xi = e^{2 pi i / 6},
#
# que es REGULAR y de orden 6 en el grupo adjunto, la representacion irreducible de peso maximo
# lambda = (1,0) -- la estandar, de dimension 4 -- tiene caracter CERO en g, y sin embargo el
# centralizador de (lambda+rho)(zeta_6) en el dual es un TORO, igual que el de g.
#
# O sea: la dicotomia de centralizadores NO decide la anulacion fuera del elemento principal.  Y g
# no es principal: npp_principal_check.sage lo mide (NINGUN C_m tiene su multiconjunto de valores).
#
# TRES RUTAS, y ninguna comparte codigo con las otras:
#   R1  Freudenthal: suma de los pesos de la irreducible, con multiplicidad, evaluados en g.
#   R2  bialternante: det(x_i^{a_j} - x_i^{-a_j}) / det(x_i^{rho_j} - x_i^{-rho_j}).
#   R3  a mano: la representacion estandar de Sp_4 tiene autovalores xi^{+-1}, xi^{+-2}; su traza es
#       (xi + xi^5) + (xi^2 + xi^4) = 2cos(60) + 2cos(120) = 1 - 1 = 0.
#
# Y EL DIAGNOSTICO, que es lo que hay que escribir en el paper: la pared que se pierde es
#     2 a_i = 0  (mod t),  o sea  a_i = t/2,
# que es la pared AFIN del nivel, no una raiz del sistema finito.  Existe solo si t es par -- y por
# eso el caso impar, donde el elemento SI es principal, no tiene el problema.
#
# CONTROL
#   C0  el mismo circuito sobre un eta que SI sobrevive tiene que dar |tau| = 1 por las tres rutas.
#   C1  se listan TODOS los testigos de la caja, no solo el minimo, y se cuenta cuantos hay.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage npp_witness.sage

import json

print("=" * 118)
print("EL TESTIGO MINIMO DE QUE LA DICOTOMIA DE CENTRALIZADORES NO DECIDE FUERA DEL PRINCIPAL")
print("=" * 118)

RES = []
for (t, m) in [(6, 2), (8, 3)]:
    K = CyclotomicField(t)
    xi = K.gen()
    L = RootSystem("C%d" % m).ambient_space()
    pos = list(L.positive_roots())
    rho = [m - j for j in range(m)]          # rho_{C_m} = (m, m-1, ..., 1)
    W = WeylCharacterRing("C%d" % m)

    def freudenthal(eta):
        el = W(W.space().from_vector(vector([Integer(v) for v in eta])))
        s = K(0)
        for wt, mult in el.weight_multiplicities().items():
            v = [int(u) for u in wt.to_vector()]
            s += int(mult) * xi ** (sum((i + 1) * v[i] for i in range(m)) % t)
        return s

    def bialternante(eta):
        a = [int(eta[j]) + rho[j] for j in range(m)]
        num = matrix(K, m, m, lambda i, j: xi ** ((i + 1) * a[j] % t) - xi ** ((-(i + 1) * a[j]) % t))
        den = matrix(K, m, m, lambda i, j: xi ** ((i + 1) * rho[j] % t) - xi ** ((-(i + 1) * rho[j]) % t))
        dd = den.determinant()
        if dd == 0:
            return None
        return num.determinant() / dd

    def pairing_coroot(alpha, vec):
        av = [QQ(u) for u in alpha.to_vector()]
        return 2 * sum(vec[i] * av[i] for i in range(m)) / sum(av[i] * av[i] for i in range(m))

    testigos = []
    ok_tres = 0
    n = 0
    for k in range(0, 2 * t + 1):
        for e in Partitions(k, max_length=m):
            eta = tuple(list(e) + [0] * (m - len(e)))
            n += 1
            f = freudenthal(eta)
            b = bialternante(eta)
            if b is None:
                continue
            if f != b:
                print("   !! R1 y R2 DISCREPAN en eta=%s : %s vs %s" % (str(eta), str(f), str(b)))
                continue
            ok_tres += 1
            a = [QQ(eta[j] + rho[j]) for j in range(m)]
            nraiz = sum(1 for al in pos if pairing_coroot(al, a) % t == 0)
            if nraiz == 0 and f == 0:
                testigos.append((tuple(int(v) for v in eta), [int(x) for x in a], [int(x) % int(t) for x in a]))

    print("")
    print("  t=%d  C_%d  h=%d :  eta probados %d | R1 == R2 en %d | TESTIGOS (centralizador toro y tau=0): %d"
          % (t, m, 2 * m, n, ok_tres, len(testigos)))
    if testigos:
        eta, a, c = testigos[0]
        print("     el minimo:  eta=%s   a = eta+rho = %s   a mod %d = %s" % (str(eta), str(a), t, str(c)))
        print("     tau por Freudenthal = %s | por bialternante = %s"
              % (str(freudenthal(eta)), str(bialternante(eta))))
        print("     ninguna coraiz cumple %d | <a, alpha^v>, y sin embargo el caracter se anula." % t)
        print("     la pared que lo mata es  2*a_i = 0 (mod %d), o sea a_i = %d : %s"
              % (t, t // 2, str([x for x in c if (2 * x) % t == 0 and x != 0])))
        print("     los %d primeros testigos: %s" % (min(6, len(testigos)),
                                                     str([x[0] for x in testigos[:6]])))
    RES.append({"t": int(t), "m": int(m), "n_eta": int(n), "R1_eq_R2": int(ok_tres),
                "n_testigos": int(len(testigos)),
                "testigo_minimo": {"eta": [int(x) for x in testigos[0][0]],
                                   "a": testigos[0][1], "a_mod_t": testigos[0][2]} if testigos else None,
                "testigos": [[int(x) for x in w[0]] for w in testigos[:20]]})

print("")
print("=" * 118)
print("  CONTROL C0: un eta que SI sobrevive, por las tres rutas, en t=6")
K = CyclotomicField(6); xi = K.gen()
W = WeylCharacterRing("C2")
el = W(W.space().from_vector(vector([Integer(2), Integer(0)])))
s = K(0)
for wt, mult in el.weight_multiplicities().items():
    v = [int(u) for u in wt.to_vector()]
    s += int(mult) * xi ** (sum((i + 1) * v[i] for i in range(2)) % 6)
print("    eta=(2,0):  a = (4,1),  a mod 6 = (4,1)  -> ninguna pared  ->  tau = %s" % str(s))
json.dump(RES, open("npp_witness_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
