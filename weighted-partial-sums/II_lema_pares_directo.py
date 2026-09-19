# -*- coding: utf-8 -*-
r"""II_lema_pares_directo.py -- el lema de pares, contra las sumas de subconjuntos CALCULADAS.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

El lema de pares (Lemma intervalos (b) de la nota II) promete: con los pares de II_certificados_A.py, las sumas de subconjuntos de Q
cubren todo [L, P_Q - L].  Aqui se calculan las sumas de subconjuntos de verdad (bitset en un entero
de Python) en los primeros bloques y se comprueba la inclusion.  Se informa tambien el hueco REAL
(el menor L' con [L', P_Q - L'] cubierto) para ver cuanto regala el lema.
Control negativo: con Q sin su menor primo q0 y sin los pares, el mismo L deberia fallar en algun
bloque (si nunca falla, el test no discrimina).  MEDIDO: no falla -- el hueco real (~L/3.5) es muy
inferior al del lema, asi que este control no discrimina; lo que cuenta es la columna L(real) < L(lema).
Uso: python II_lema_pares_directo.py [NBLOQUES]
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 12
N = 200000
criba = bytearray([1]) * (N + 1)
criba[0:2] = b"\x00\x00"
for i in range(2, int(N ** .5) + 1):
    if criba[i]:
        criba[i * i::i] = bytearray(len(range(i * i, N + 1, i)))


def pares_de(Q, a, b):
    Qs = set(Q); usados = {Q[0]}; q0 = Q[0]; H = U = 0; pares = []
    while 2 * H - q0 + 2 < b:
        hecho = False
        for d in range(min(H + 1, (a - Q[0]) // 2), 0, -1):
            for u in Q:
                v = u + 2 * d
                if v > a:
                    break
                if u in usados or v in usados or v not in Qs:
                    continue
                pares.append((u, v)); usados.update((u, v)); H += d; U += u; hecho = True
                break
            if hecho:
                break
        assert hecho
    return pares, U + q0 - 1


def sumas(Q):
    bs = 1
    for q in Q:
        bs |= bs << q
    return bs


def hueco_real(bs, P):
    Lr = 0
    while Lr <= P // 2:
        ok = all((bs >> x) & 1 for x in (Lr,))  # candidato
        # menor L' con [L', P-L'] cubierto: buscar el mayor x <= P/2 no cubierto
        break
    fuera = [x for x in range(0, P // 2 + 1) if not (bs >> x) & 1]
    return (max(fuera) + 1) if fuera else 0


a = 5000; nb = 0; fallos = 0; ctrl_falla = 0
while nb < NB:
    b = 11 * a // 10
    Q = [p for p in range(b // 2 + 1, a + 1) if criba[p]]
    P = sum(Q)
    pares, L = pares_de(Q, a, b)
    bs = sumas(Q)
    cubre = all((bs >> x) & 1 for x in range(L, P - L + 1))
    Lr = hueco_real(bs, P)
    fallos += not cubre
    # control: quitar q0 y los primos de los pares; mismo L
    quitar = {Q[0]} | {x for pr in pares for x in pr}
    Qc = [q for q in Q if q not in quitar]
    bc = sumas(Qc); Pc = sum(Qc)
    cc = all((bc >> x) & 1 for x in range(L, Pc - L + 1))
    ctrl_falla += not cc
    print("bloque %2d [%6d,%6d] |Q|=%4d  L(lema)=%7d  L(real)=%5d  cubre=%s  control(sin pares) cubre=%s"
          % (nb + 1, a, b, len(Q), L, Lr, cubre, cc)); sys.stdout.flush()
    nb += 1; a = b + 1
print("bloques: %d, fallos del lema: %d, controles que fallan: %d" % (nb, fallos, ctrl_falla))
