# -*- coding: utf-8 -*-
# LA DIAGONAL, PARA TODO CARACTER IMPAR:  ¿T = 0  <=>  chi(2) = 1 ?
#
# El borrador cierra la diagonal m = (f-1)/2 solo para el caracter CUADRATICO, via la formula de
# numero de clases de Dirichlet, y deja "falta un criterio uniforme".  El barrido de
# caracteres_orden_mayor.py encuentra seis casillas con un caracter impar de orden > 2 de
# componente nula, y las SEIS estan en esa misma diagonal.  En las seis, chi(2) = 1.
#
# Eso sugiere que la diagonal no es un fenomeno cuadratico sino uno de <2>:
#
#     T(chi) = sum_{c=1}^{(f-1)/2} c chi(c) = (f/2) B_{1,chi} (chi(2)^{-1} - 1)
#
# y como B_{1,chi} != 0 para chi impar, T = 0 <=> chi(2) = 1.  Para el cuadratico eso es
# (2|f) = 1, o sea f = 7 (mod 8): exactamente la dicotomia del borrador.
#
# Se mide lo UNO y lo OTRO: el criterio y la forma cerrada, en aritmetica ciclotomica exacta, para
# f impar compuesto tambien, y para caracteres NO primitivos.  Con dos senuelos que deben fallar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  sh sg.sh diagonal_chi2.sage

print("=" * 100)
print("f impar, TODO caracter impar mod f:  T = sum_{c<=(f-1)/2} c chi(c)")
print("=" * 100)
crit_ok = crit_mal = forma_ok = forma_mal = 0
prim_ok = prim_mal = 0
nulos = 0
ordenes_nulos = {}
prim_no = 0
for f in range(3, 122, 2):
    G = DirichletGroup(f)
    for chi in G:
        if chi(-1) != -1:
            continue
        T = sum(c * chi(c) for c in range(1, (f - 1) // 2 + 1))
        B1 = sum(a * chi(a) for a in range(1, f)) / f
        c2 = chi(2)
        # el criterio.  Se separa PRIMITIVO de imprimitivo: si chi viene inducido de conductor f*
        # menor, B_{1,chi} = B_{1,chi*} prod_{p | f} (1 - chi*(p)) se ANULA en cuanto algun p | f
        # tiene chi*(p) = 1, y entonces T = 0 sin que chi(2) = 1.  Ahi el criterio no puede valer.
        prim = (chi.conductor() == f)
        if (T == 0) == (c2 == 1):
            crit_ok += 1
            if prim:
                prim_ok += 1
        else:
            crit_mal += 1
            if prim:
                prim_mal += 1
            print("   *** CRITERIO FALLA: f=%d orden=%d chi(2)=%s T=%s B_1=%s primitivo=%s ***" % (
                f, chi.order(), c2, T, B1, prim))
        # la forma cerrada
        pred = (f / 2) * B1 * (c2 ** (-1) - 1)
        if T == pred:
            forma_ok += 1
        else:
            forma_mal += 1
            print("   *** FORMA FALLA: f=%d orden=%d  T=%s  pred=%s ***" % (
                f, chi.order(), T, pred))
        if T == 0:
            nulos += 1
            o = chi.order()
            ordenes_nulos[o] = ordenes_nulos.get(o, 0) + 1
            if chi.conductor() != f:
                prim_no += 1
print()
print("caracteres impares examinados: %d" % (crit_ok + crit_mal))
print("criterio  T = 0 <=> chi(2) = 1 :  %d si, %d no" % (crit_ok, crit_mal))
print("   restringido a caracteres PRIMITIVOS mod f: %d si, %d no" % (prim_ok, prim_mal))
print("forma cerrada T = (f/2) B_1 (chi(2)^-1 - 1) :  %d si, %d no" % (forma_ok, forma_mal))
print("con T = 0: %d  | por orden del caracter: %s" % (nulos, sorted(ordenes_nulos.items())))
print("de ellos, NO primitivos mod f: %d" % prim_no)

print()
print("=" * 100)
print("SENUELOS (tienen que fallar; si pasan, el criterio no dice nada)")
print("=" * 100)
# S1: fuera de la diagonal, m = (f-1)/2 - 1, el criterio no debe valer
s1_ok = s1_mal = 0
for f in range(5, 60, 2):
    G = DirichletGroup(f)
    for chi in G:
        if chi(-1) != -1:
            continue
        m = (f - 1) // 2 - 1
        if m < 1:
            continue
        T = sum(c * chi(c) for c in range(1, m + 1))
        if (T == 0) == (chi(2) == 1):
            s1_ok += 1
        else:
            s1_mal += 1
print("S1  el mismo criterio en m = (f-1)/2 - 1: %d coinciden, %d NO  (debe haber muchos NO)" % (
    s1_ok, s1_mal))
# S2: caracteres PARES en la diagonal
s2_ok = s2_mal = 0
for f in range(5, 60, 2):
    G = DirichletGroup(f)
    for chi in G:
        if chi(-1) != 1 or chi.is_trivial():
            continue
        T = sum(c * chi(c) for c in range(1, (f - 1) // 2 + 1))
        if (T == 0) == (chi(2) == 1):
            s2_ok += 1
        else:
            s2_mal += 1
print("S2  el criterio para caracteres PARES: %d coinciden, %d NO  (debe haber NO)" % (
    s2_ok, s2_mal))

print()
print("=" * 100)
print("EL CUADRATICO COMO COROLARIO:  chi(2) = 1  <=>  f = +-1 (mod 8)")
print("=" * 100)
print("%-6s %-10s %-10s %-10s %s" % ("f", "f mod 8", "chi(2)", "T", "prop:dirichlet"))
for f in [7, 11, 15, 19, 23, 31, 35, 39, 43, 47, 51, 55, 59]:
    if f % 4 != 3:
        continue
    G = DirichletGroup(f)
    chis = [chi for chi in G if chi.order() == 2 and chi(-1) == -1 and chi.conductor() == f]
    if not chis:
        continue
    chi = chis[0]
    T = sum(c * chi(c) for c in range(1, (f - 1) // 2 + 1))
    print("%-6d %-10d %-10s %-10s %s" % (
        f, f % 8, chi(2), T, "0 si f=7(8)" if f % 8 == 7 else "!=0 si f=3(8)"))
