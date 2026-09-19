# -*- coding: utf-8 -*-
# ¿SOLO MUERE EL CARACTER CUADRATICO?   El borrador decia que si, en su rango.  Aqui se barre el
# MISMO rango que dibuja fig_ceros -- conductores primos f < 160, todo m < f -- y se pregunta por
# TODOS los caracteres impares, no solo los de orden 2.
#
# ARITMETICA EXACTA, sin flotantes.  Para f primo el grupo es ciclico: fijado un generador g, el
# caracter de orden EXACTAMENTE d es chi(c) = zeta_d^{ind(c)}, con ind el logaritmo discreto en
# base g.  La suma  T = sum_{c<=m} c chi(c)  se acumula como vector entero en Z[x]/(x^d - 1)
# plegando los indices modulo d, y se anula si y solo si Phi_d divide a ese polinomio.  Phi_d se
# construye por divisiones enteras exactas desde x^d - 1.
#
# Los caracteres de orden exactamente d son chi^j con gcd(j,d) = 1, que son los CONJUGADOS de
# Galois de chi: si uno tiene componente nula, todos.  Basta un representante por d.
#
# El caracter es impar si y solo si chi(-1) = zeta_d^{(f-1)/2} = -1, o sea d par y
# (f-1)/2 = d/2 (mod d).
#
# CONTROL C1.  Para d = 2 la suma es sum_{c<=m} c (c|f) con el simbolo de Kronecker, que es otro
# camino de codigo: se calculan los dos y se comparan casilla a casilla.  Un PASA sobre 0 casillas
# no es un PASA, asi que se imprime el numero de casillas comparadas.
#
# PRIMERA VERSION EQUIVOCADA, registrada para no repetirla: se indexaba chi(c) = zeta_d^{e ind(c)}
# con e = (f-1)/d, que es el caracter de orden d/gcd(d,e), NO de orden d.  Daba 10 casillas de
# orden > 2, y ninguna lo era.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python caracteres_orden_mayor.py > caracteres_orden_mayor_OUT.txt


def primos(n):
    crib = [True] * (n + 1)
    crib[0] = crib[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if crib[i]:
            for j in range(i * i, n + 1, i):
                crib[j] = False
    return [i for i in range(n + 1) if crib[i]]


def divide(a, b):
    """division exacta de polinomios enteros a/b, listas de coeficientes de menor a mayor."""
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] // b[-1]
        q[i] = c
        for j, bb in enumerate(b):
            a[i + j] -= c * bb
    assert all(x == 0 for x in a), "division no exacta"
    return q


CICL = {}


def ciclotomico(n):
    """Phi_n como lista de coeficientes, por division entera desde x^n - 1."""
    if n in CICL:
        return CICL[n]
    p = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            p = divide(p, ciclotomico(d))
    CICL[n] = p
    return p


def divisible(a, b):
    """¿b divide a a, sobre Z?  a y b con b monico."""
    a = a[:]
    while len(a) >= len(b):
        if a[-1] == 0:
            a.pop()
            continue
        c = a[-1]
        if c % b[-1]:
            return False
        k = len(a) - len(b)
        for j, bb in enumerate(b):
            a[k + j] -= (c // b[-1]) * bb
        while a and a[-1] == 0:
            a.pop()
    return all(x == 0 for x in a)


def raiz_primitiva(f):
    for g in range(2, f):
        vis, x = set(), 1
        for _ in range(f - 1):
            x = x * g % f
            vis.add(x)
        if len(vis) == f - 1:
            return g
    raise RuntimeError("sin raiz primitiva")


def kron(a, n):
    """simbolo de Kronecker (a|n) para n impar, el otro camino de codigo del control C1."""
    if n == 1:
        return 1
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


print("CONDUCTORES PRIMOS f < 160, TODO m < f: ¿que caracteres impares tienen componente nula?")
print("%-5s %-5s %-6s %s" % ("f", "m", "orden", "T = sum_{c<=m} c chi(c)"))
total = orden2 = mayor = 0
casos = []
ctrl_ok = ctrl_mal = 0
for f in primos(160):
    if f < 5:
        continue
    g = raiz_primitiva(f)
    ind = [0] * f
    x = 1
    for k in range(f - 1):
        ind[x] = k
        x = x * g % f
    ordenes = [d for d in range(2, f) if (f - 1) % d == 0 and d % 2 == 0]
    for m in range(1, f):
        for d in ordenes:
            # chi de orden EXACTAMENTE d:  chi(c) = zeta_d^{ind(c)};  impar <=> chi(-1) = -1
            if ((f - 1) // 2) % d != d // 2:
                continue
            vec = [0] * d
            for c in range(1, m + 1):
                vec[ind[c] % d] += c
            nulo = divisible(vec[:], ciclotomico(d))
            if d == 2:                       # C1: el mismo cero, por Kronecker
                nulo2 = (sum(c * kron(c, f) for c in range(1, m + 1)) == 0)
                if nulo == nulo2:
                    ctrl_ok += 1
                else:
                    ctrl_mal += 1
                    print("   *** C1 DISCREPA en f=%d m=%d: %s vs %s ***" % (f, m, nulo, nulo2))
            if not nulo:
                continue
            total += 1
            if d == 2:
                orden2 += 1
            else:
                mayor += 1
                casos.append((f, m, d))
                print("%-5d %-5d %-6d %s" % (f, m, d, "0  *** ORDEN > 2 ***"))
print()
print("C1 (orden 2 por dos caminos): %d casillas comparadas, %d coinciden, %d discrepan" % (
    ctrl_ok + ctrl_mal, ctrl_ok, ctrl_mal))
print("componentes impares nulas: %d | de orden 2: %d | de orden > 2: %d" % (total, orden2, mayor))
print("casillas (f, m, orden) con orden > 2:", casos)
