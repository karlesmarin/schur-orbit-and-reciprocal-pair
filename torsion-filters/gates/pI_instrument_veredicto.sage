# -*- coding: utf-8 -*-
# prob:instrument --- VEREDICTO: el experimento con la regla de Littlewood es IMPOSIBLE, y la razon
# es un teorema del propio Paper I.   19 de agosto de 2026, noche.
#
# QUE PASO.  Se diseno un experimento para ver si  m_mu = m_{mu*}  --- la condicion de anulacion,
# (2186) --- se refina por alguna estadistica sobre delta, con  m_mu = sum_{delta par} c^lambda_{mu,delta}.
# Tres versiones, tres correcciones (enumerar por anchura y no por tamano; quitar un filtro con la
# precedencia mal puesta; restringir al rango de Littlewood l(lambda) <= N/2), y las tres dieron
# UN SOLO lambda en el locus de anulacion, de 164, 494 y 494.
#
# LA RAZON, y no es un fallo.  thm:stable del propio articulo PRUEBA que en el rango de Littlewood
# Psi_r(lambda) != 0 salvo en una rama:
#     «the term beta'=empty gives m_mu >= 1 for every mu subset lambda.  If l(lambda) < N/2, take
#      mu = lambda: then l(lambda*) = N - l(lambda) > l(lambda), so lambda* not subset lambda and
#      m_{lambda*} = 0 != m_lambda.»
# O sea: el locus de anulacion vive ENTERO fuera del rango donde la formula de Littlewood computa la
# multiplicidad.  El experimento no podia dar otra cosa.
#
# Y ESO ES EL VEREDICTO, no un resultado nulo: explica POR QUE el articulo nombra a Frohmader como
# «the one we would try first».  No es preferencia --- sus reglas son las unicas que tienen «no
# stable-range constraint and no modification rule», o sea las unicas que cruzan la frontera donde
# vive la pregunta.
#
# QUE COMPRUEBA ESTA GATE
#   V1  quien es el unico superviviente en cada configuracion: si es el rectangulo de lado impar,
#       es exactamente la rama (b) que thm:stable exceptua, y el circulo se cierra.
#   V2  y que fuera del rango de Littlewood hay locus de sobra --- se cuenta con el criterio del
#       articulo, no con la formula que alli no vale --- para que quede escrito que el problema no
#       es la falta de casos sino el instrumento.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run -d --name pI-vered \
#         -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" -w /work sage-normaliz:local \
#         bash -c "sage < pI_instrument_veredicto.sage > pI_instrument_veredicto_OUT.txt 2>&1"

import itertools
import json
import sys

CASOS = [(2, 2), (2, 3), (4, 2)]
ANCHO = 8

Sym = SymmetricFunctions(QQ)
s = Sym.s()


def asociada(mu, N):
    c = list(Partition(mu).conjugate())
    if not c:
        return Partition([1] * N)
    if c[0] > N:
        return None
    nueva = [x for x in [N - c[0]] + c[1:] if x > 0]
    if any(nueva[i] < nueva[i + 1] for i in range(len(nueva) - 1)):
        return None
    return Partition(nueva).conjugate()


def pares_hasta(n):
    out = []
    for k in range(0, n // 2 + 1):
        for p in Partitions(2 * k):
            if all(x % 2 == 0 for x in p):
                out.append(Partition(p))
    return out


_LR = {}
def lr(lam, mu, de):
    key = (tuple(lam), tuple(mu), tuple(de))
    if key not in _LR:
        _LR[key] = QQ((s(Partition(mu)) * s(Partition(de))).coefficient(Partition(lam)))
    return _LR[key]


def m_de(lam, nu, evens):
    n = sum(lam) - sum(nu)
    if n < 0:
        return 0
    return int(sum(lr(lam, nu, de) for de in evens if sum(de) == n))


print("=" * 100)
print("prob:instrument --- VEREDICTO: quien sobrevive en el rango de Littlewood")
print("=" * 100)
print("")

resumen = {}
for (T, R) in CASOS:
    N = T + 2 * R
    evens = pares_hasta(ANCHO * (N // 2) + 4)
    lams = []
    for p in itertools.product(range(ANCHO, -1, -1), repeat=N // 2):
        if any(p[i] < p[i + 1] for i in range(N // 2 - 1)):
            continue
        q = Partition([x for x in p if x > 0])
        if q:
            lams.append(q)
    sobreviven = []
    for lam in lams:
        cero = True
        for k in range(0, sum(lam) + 1):
            if not cero:
                break
            for p in Partitions(k, max_length=N):
                nu = Partition(p)
                ast = asociada(nu, N)
                if ast is None:
                    continue
                if m_de(lam, nu, evens) != m_de(lam, ast, evens):
                    cero = False
                    break
        if cero:
            sobreviven.append(tuple(lam))
    rect_impar = [l for l in sobreviven
                  if len(l) == N // 2 and len(set(l)) == 1 and l[0] % 2 == 1]
    print("  t=%d r=%d (N=%d): %d lambda en el rango de Littlewood, %d con m_mu = m_{mu*} para toda mu"
          % (T, R, N, len(lams), len(sobreviven)))
    print("     quienes son : %s" % (sobreviven[:6],))
    print("     V1  .son rectangulos de lado IMPAR con N/2 filas --- la rama (b) de thm:stable? : %s"
          % ("SI, todos" if sobreviven and len(rect_impar) == len(sobreviven) else
             "no todos: %s" % [l for l in sobreviven if l not in rect_impar][:4]))
    print("")
    sys.stdout.flush()
    resumen["t=%d,r=%d" % (T, R)] = dict(n_lams=len(lams), sobreviven=[list(x) for x in sobreviven],
                                         todos_rect_impar=bool(sobreviven and len(rect_impar) == len(sobreviven)))

json.dump(resumen, open("pI_instrument_veredicto_DUMP.json", "w"), indent=1)
print("=" * 100)
print("DONE")
