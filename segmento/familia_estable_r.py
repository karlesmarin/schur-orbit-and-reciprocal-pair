# -*- coding: utf-8 -*-
# prop:fam CON EL MODULO r:  ¿ r | n-1, n, n+1  =>  el segmento es estable en q' ?
#
# El enunciado del borrador usaba q' = M/p^v.  Al corregirlo a r = q/p^{v_p(q)} hay que rehacer la
# comprobacion en el caso PAR, donde q' = 2r y el segmento son los IMPARES modulo 2r: alli la clase
# singleton no es el centro 0 sino c = r, que es la unica impar con 2c = 0 (mod 2r).
#
# Se mide la implicacion en los dos sentidos y se separa por paridad de n y de r, porque el
# argumento es distinto en cada rama.  Tambien se cuenta el reciproco, que es FALSO y debe verse.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python familia_estable_r.py > familia_estable_r_OUT.txt


def factoriza(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def mcd(a, b):
    while b:
        a, b = b, a % b
    return a


def exps_su(n):
    if n % 2 == 1:
        return [(n + 1 - 2 * j) // 2 for j in range(1, n + 1)]
    return [n + 1 - 2 * j for j in range(1, n + 1)]


def estable(exps, qp):
    base = sorted(e % qp for e in exps)
    for u in range(1, qp):
        if mcd(u, qp) == 1 and sorted((u * e) % qp for e in exps) != base:
            return False
    return True


print("prop:fam con r = q/p^{v_p(q)}:  r | n-1,n,n+1  =>  estable en q' = M/p^{v_p(M)}")
print("%-4s %-4s %-4s %-5s %-5s %-4s %-8s %-8s %s" % (
    "n", "q", "p", "M", "q'", "r", "familia", "estable", "veredicto"))
ok = mal = 0
recip_si = recip_no = 0
fuera_estables = []
ramas = {}
for n in range(3, 26):
    for q in range(3, 61):
        M = q if n % 2 == 1 else 2 * q
        ex = exps_su(n)
        for p in factoriza(M):
            v = 0
            t = M
            while t % p == 0:
                t //= p
                v += 1
            qp = M // p ** v
            if qp < 3:
                continue
            r = q
            while r % p == 0:
                r //= p
            fam = any(x % r == 0 for x in (n - 1, n, n + 1)) if r > 1 else True
            est = estable(ex, qp)
            if fam:
                if est:
                    ok += 1
                else:
                    mal += 1
                    print("%-4d %-4d %-4d %-5d %-5d %-4d %-8s %-8s *** FALLA ***" % (
                        n, q, p, M, qp, r, "si", "NO"))
                rama = ("n impar" if n % 2 else "n par") + ", " + ("r impar" if r % 2 else "r par")
                ramas[rama] = ramas.get(rama, 0) + 1
            else:
                if est:
                    recip_si += 1      # estable FUERA de familia: el reciproco es falso
                    # y ESTO es lo que importa para prop:sop -- si alguna de estas tuviera
                    # d >= 2 y E >= 2, el Teorema 3(d) forzaria defecto con un primo FUERA de
                    # las familias, y prop:sop seria falsa.  Se imprimen todas con su d y su E.
                    E = 1
                    if v > 0:
                        E = p ** (v - 1) * (p - 1)
                    d = 1
                    cnt = sum(1 for c in range(1, qp) if mcd(c, qp) == 1)
                    d = cnt // 2
                    fuera_estables.append((n, q, p, qp, r, d, E))
                else:
                    recip_no += 1
print()
print("en familia: %d casillas | estables: %d | NO estables: %d" % (ok + mal, ok, mal))
print("por rama:", sorted(ramas.items()))
print("fuera de familia: %d estables (el reciproco es FALSO) y %d inestables" % (
    recip_si, recip_no))
print()
print("LAS ESTABLES FUERA DE FAMILIA, UNA A UNA:  si alguna tuviera d >= 2 y E >= 2, el Teorema")
print("3(d) forzaria defecto y prop:sop seria FALSA.  (n, q, p, q', r, d, E):")
for c in fuera_estables:
    print("   n=%-3d q=%-3d p=%-3d q'=%-3d r=%-3d d=%-3d E=%-3d %s" % (
        c + ("*** d>=2 y E>=2: MIRAR ***" if c[5] >= 2 and c[6] >= 2 else "",)))
print()
print("Y LA CLASE SINGLETON: las c con 2c = 0 (mod q') son c = 0 y, si q' es par, c = q'/2.")
print("En el caso par el segmento es impar modulo q' = 2r, asi que la unica disponible es c = r.")
for (n, q) in [(4, 15), (6, 35), (8, 21), (10, 33)]:
    M = 2 * q
    for p in factoriza(M):
        if p == 2:
            continue
        v = 0
        t = M
        while t % p == 0:
            t //= p
            v += 1
        qp = M // p ** v
        r = q
        while r % p == 0:
            r //= p
        if qp != 2 * r or r < 3:
            continue
        ex = exps_su(n)
        N = [0] * qp
        for e in ex:
            N[e % qp] += 1
        raras = [c for c in range(qp) if N[c] not in (min(x for x in N if x) if any(N) else 0,
                                                     max(N))]
        print("   n=%d q=%d p=%d q'=%d r=%d  N = %s  clase c=r vale %d" % (
            n, q, p, qp, r, N, N[r % qp]))
