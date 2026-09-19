# -*- coding: utf-8 -*-
# (1) EL CASO n PAR, LOS DOS REGIMENES.   (2) EL HUECO DEL CRITERIO DE SOPORTE.
#
# (1)  Para n = 2m el segmento son los impares {+-1, +-3, ..., +-(2m-1)} sobre M = 2q.  Localmente:
#        p = 2     => v_2(M) = v_2(q)+1 y q' es IMPAR.  2 es unidad: multiplicar por 2^{-1} lleva
#                     los impares a un segmento CONSECUTIVO centrado en q'/2 -- medio periodo --
#                     que es la segunda familia del paper.  Como multiplicar por una unidad es
#                     Galois, el rango de la matriz de momentos no cambia.
#        p impar   => q' = 2r con r la parte prima-a-p de q, o sea q' PAR.  Ahi 2 no es unidad y
#                     entra el Lema de medio periodo y duplicacion: L_{2r} = 2 L_r, y para p impar
#                     M_{2r} y M_r tienen el MISMO rango modulo p.
#      Se comprueban las dos reducciones contra W_1 medido con el instrumento de capas.
#
# (2)  El criterio de soporte tiene dos huecos declarados; aqui se mide el segundo:
#        "estable pero FUERA de familia  =>  W_1 = d".
#      Si sale entero, el criterio de familias queda medido de punta a punta y solo falta probarlo.
#
# Authors: Carles Marin, Claude (AI assistant).
load("inst_tipoA.sage")


def S_de(exps, qp):
    S = [0] * qp
    for e in exps:
        S[e % qp] += e
    return S


def rango_S(exps, qp, p):
    S = S_de(exps, qp)
    units = [u for u in range(1, qp) if gcd(u, qp) == 1 and u < qp - u]
    if not units:
        return 0
    filas = [[S[(inverse_mod(u, qp) * c) % qp] for c in range(qp)] for u in units]
    return matrix(GF(p), filas).rank()


def estable(exps, qp):
    base = sorted([e % qp for e in exps])
    for u in range(1, qp):
        if gcd(u, qp) == 1 and sorted([(u * e) % qp for e in exps]) != base:
            return False
    return True


print("=" * 98)
print("(1a) p = 2, q' IMPAR: ¿el segmento x 2^{-1} es un consecutivo centrado en q'/2?")
print("=" * 98)
print("%-3s %-5s %-4s %-30s %-24s %s" % ("n", "q", "q'", "impares mod q'", "x 2^{-1} ordenado", "consecutivo?"))
bien = mal = 0
for n in range(4, 13, 2):
    for q in range(3, 40, 2):        # q impar => q' = q, p = 2 con v = 1
        qp = q
        if qp < 3:
            continue
        ex = exps_su(n)
        img = sorted([(inverse_mod(2, qp) * e) % qp for e in ex])
        # consecutivo (ciclicamente) si los valores distintos forman un bloque
        dis = sorted(set(img))
        cons = all((dis[i + 1] - dis[i]) == 1 for i in range(len(dis) - 1))
        centro = (sum(dis) / len(dis)) if dis else 0
        if cons:
            bien += 1
        else:
            mal += 1
        if n <= 8 and q <= 21:
            print("%-3d %-5d %-4d %-30s %-24s %s  centro=%s" % (
                n, q, qp, str(sorted([e % qp for e in ex]))[:30], str(img)[:24],
                "si" if cons else "NO", centro))
print("   consecutivos: %d | no: %d" % (bien, mal))

print()
print("=" * 98)
print("(1b) p IMPAR, q' = 2r PAR: ¿rango modulo p del segmento impar mod 2r == rango mod r?")
print("=" * 98)
print("%-3s %-5s %-3s %-5s %-4s %-7s %-9s %-9s %s" % (
    "n", "q", "p", "q'=2r", "r", "W_1", "rango 2r", "rango r", "veredicto"))
ok = malo = 0
for n in range(4, 11, 2):
    for q in range(3, 46):
        M = 2 * q
        if euler_phi(M) > 32:
            continue
        for p, _ in factor(M):
            if p == 2:
                continue
            v = valuation(M, p)
            qp = M // p ** v
            if qp < 3 or euler_phi(p ** v) < 2 or qp % 2 != 0:
                continue
            r = qp // 2
            if r < 3:
                continue
            ex = exps_su(n)
            try:
                I = ImgSeg(M, p, ex, n - 1, 2)
            except Exception:
                continue
            W = I.W()
            if len(W) < 2:
                continue
            r2r = rango_S(ex, qp, p)
            rr = rango_S(ex, r, p)
            coin = (r2r == rr)
            if coin:
                ok += 1
            else:
                malo += 1
            print("%-3d %-5d %-3d %-5d %-4d %-7d %-9d %-9d %s" % (
                n, q, p, qp, r, W[1], r2r, rr, "ok" if coin else "*** DISCREPA ***"))
print("   coinciden: %d | discrepan: %d" % (ok, malo))
if malo == 0 and ok:
    print("   *** la duplicacion vale: el segmento impar mod 2r mide lo mismo que mod r ***")

print()
print("=" * 98)
print("(2) HUECO DEL SOPORTE: estable pero FUERA de familia => W_1 = d ?")
print("=" * 98)
print("%-3s %-5s %-3s %-5s %-4s %-6s %-9s %s" % ("n", "q", "p", "q'", "d", "W_1", "W_1 = d?", "familia"))
llena = corta = 0
for n in range(3, 11):
    for q in range(3, 61):
        M = q if n % 2 == 1 else 2 * q
        if euler_phi(M) > 32:
            continue
        for p, _ in factor(M):
            v = valuation(M, p)
            qp = M // p ** v
            if qp < 3 or euler_phi(p ** v) < 2:
                continue
            ex = exps_su(n)
            if not estable(ex, qp):
                continue
            qq = q // p ** valuation(q, p)
            enfam = any(x % qq == 0 for x in (n - 1, n, n + 1)) if qq >= 3 else True
            if enfam:
                continue
            try:
                I = ImgSeg(M, p, ex, n - 1, 2)
            except Exception:
                continue
            W = I.W()
            if len(W) < 2:
                continue
            d = euler_phi(qp) / 2
            if W[1] == d:
                llena += 1
            else:
                corta += 1
            print("%-3d %-5d %-3d %-5d %-4s %-6d %-9s %s" % (
                n, q, p, qp, d, W[1], "si" if W[1] == d else "*** NO ***", "fuera"))
print("   casillas estables fuera de familia: %d | con W_1 = d: %d | cortas: %d" % (
    llena + corta, llena, corta))
if corta == 0 and llena:
    print("   *** el criterio de familias queda medido de punta a punta ***")
elif corta:
    print("   hay %d casillas estables fuera de familia con capa CORTA: ahi podria haber indice" % corta)
