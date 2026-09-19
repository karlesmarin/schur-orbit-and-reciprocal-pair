# -*- coding: utf-8 -*-
r"""II_certificados_A.py -- certificados PROPIOS de la linea A para 5000 <= m < 10^8.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Prueba (construccion propia de los pares y comprobacion en enteros):
 * Segundo momento: para Q subconjunto de primos de (m/2, m], con signos aleatorios en los demas
   primos, E R_Q^2 = sum_{d libre de cuadrados} b_m(d)^2 - sum_{q in Q} q^2 <= m^3 (log m + 13)/9,
   b_m(d) = d sum_{j <= sqrt(m/d)} j^2.  Luego existe una eleccion con |R_Q| <= B*(m),
   B*(m)^2 = m^3 (log m + 13)/9.
 * Lema de pares: primos distintos q0, (u_j, v_j) de Q con d_j = (v_j - u_j)/2, 1 <= d_j <= H_{j-1}+1,
   H = sum d_j, U = sum u_j, L = U + q0 - 1; si 2H - q0 + 2 >= b, las sumas con signo de Q cubren
   todos los enteros de la paridad de P_Q en [-D, D], D = P_Q - 2L.
 * Un bloque [a, b] con Q = primos de (b/2, a] sirve para todo m en [a, b]; basta 90 D^2 > b^3 (130 + 7e),
   e = ceil(log2 b)  (log b < 7e/10), en enteros.
Bloques a_1 = 5000, b_j = min(floor(11 a_j / 10), 10^8 - 1), a_{j+1} = b_j + 1.
Construccion de los pares: voraz, d_j lo mayor posible <= H_{j-1} + 1, hasta 2H - q0 + 2 >= b.
Control: se comprueba ademas, con la desigualdad, que el lema FALLA si se quita el ultimo par (la
condicion (5) deja de cumplirse) -- el certificado no es vacuo.
Uso: python II_certificados_A.py
"""
import sys, math, time
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
N = 10 ** 8
t0 = time.time()
criba = bytearray([1]) * (N + 1)
criba[0:2] = b"\x00\x00"
for i in range(2, int(N ** .5) + 1):
    if criba[i]:
        criba[i * i::i] = bytearray(len(range(i * i, N + 1, i)))
print("criba hasta 10^8: %.0f s" % (time.time() - t0)); sys.stdout.flush()


def primos_en(lo, hi):  # primos en (lo, hi]
    return [p for p in range(lo + 1, hi + 1) if criba[p]]


a = 5000
bloques = pares_tot = 0
peor = 0.0
pref16 = True
filas_csv = ["a,b,N,pares,L,D,cota_sobre_D2"]
while a <= N - 1:
    b = min(11 * a // 10, N - 1)
    Q = primos_en(b // 2, a)
    Qs = set(Q)
    PQ = sum(Q)
    usados = set()
    q0 = Q[0]
    usados.add(q0)
    H = U = 0
    pares = []
    while 2 * H - q0 + 2 < b:
        hecho = False
        for d in range(min(H + 1, (a - Q[0]) // 2), 0, -1):
            for u in Q:
                v = u + 2 * d
                if v > a:
                    break
                if u in usados or v in usados or v not in Qs:
                    continue
                pares.append((u, v)); usados.update((u, v))
                H += d; U += u
                hecho = True
                break
            if hecho:
                break
        if not hecho:
            sys.exit("SIN PAR en el bloque [%d, %d] con H=%d" % (a, b, H))
    L = U + q0 - 1
    D = PQ - 2 * L
    e = (b - 1).bit_length()  # ceil(log2 b) para b no potencia de 2; cota superior en general
    if (1 << e) < b:
        e += 1
    ok = 90 * D * D > b ** 3 * (130 + 7 * e)
    # comprobaciones del lema, en enteros
    Hc = 0
    for (u, v) in pares:
        dj = (v - u) // 2
        assert (v - u) % 2 == 0 and 1 <= dj <= Hc + 1 and b // 2 < u < v <= a
        Hc += dj
    assert Hc == H and 2 * H - q0 + 2 >= b
    # control: sin el ultimo par la condicion (5) no se cumple
    ult = pares[-1]
    assert 2 * (H - (ult[1] - ult[0]) // 2) - q0 + 2 < b
    if not ok:
        sys.exit("FALLA la desigualdad en [%d, %d]: D=%d" % (a, b, D))
    # prefijos en {2,3,5,7}: factor 2^4 = 16 en el segundo momento, para b >= 60000
    if b >= 60000:
        pref16 &= 90 * D * D > 16 * b ** 3 * (130 + 7 * e)
    bloques += 1; pares_tot += len(pares)
    cociente = (b ** 3 * (130 + 7 * e)) / (90 * D * D)
    filas_csv.append("%d,%d,%d,%d,%d,%d,%.6e" % (a, b, len(Q), len(pares), L, D, cociente))
    peor = max(peor, cociente)
    if bloques <= 3 or bloques % 20 == 0 or b == N - 1:
        print("bloque %3d [%9d, %9d] |Q|=%7d pares=%2d H=%9d q0=%9d L=%12d D=%15d  B*^2/D^2 <= %.4f"
              % (bloques, a, b, len(Q), len(pares), H, q0, L, D, cociente))
        sys.stdout.flush()
    a = b + 1
print("bloques: %d (cubren 5000..%d), pares: %d, max (cota)/(ventana)^2 = %.4f  (debe ser < 1)"
      % (bloques, N - 1, pares_tot, peor))
print("prefijos en {2,3,5,7} (factor 16) en todos los bloques con b >= 60000: %s" % pref16)
import pathlib
(pathlib.Path(__file__).resolve().parent / "II_certificados_A_bloques.csv").write_text(
    "\n".join(filas_csv) + "\n", encoding="utf-8")
print("tiempo: %.0f s" % (time.time() - t0))
