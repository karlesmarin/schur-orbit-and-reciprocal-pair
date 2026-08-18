# -*- coding: utf-8 -*-
# EL FILTRO ES UNA FUNCION SOBRE (Z/t)^rango, NO SOBRE PARTICIONES.   16 de agosto de 2026.
#
# DE DONDE SALE.  Cruzar (R) con la formula del signo del Lema 3.1.
#
#   (R) dice:  tau_t(eta) != 0  <=>  (eta+rho)(xi) regular,  y esa condicion se escribe entera con
#       los residuos  c_j = a_j mod t,  a = eta + rho.
#   El Lema 3.1 da el SIGNO tambien en terminos de los c_j: los epsilon_j = +-1 segun c_j <= m, y la
#       permutacion sigma que ordena min(c_j, t-c_j).
#
# Componiendo:  tau_t(eta)  deberia depender SOLO de  a mod t.  O sea el filtro no es una funcion
# sobre particiones sino sobre  (Z/t)^rango,  un conjunto FINITO.  Consecuencias, si vale:
#
#   * tau es invariante por  eta_j -> eta_j + t,  luego el paso natural en el lado del BRANCHING es
#     t -- el mismo tamaño que el t-ribbon del lado de la orbita.  Eso convierte la observacion
#     "el tamaño del paso no esta ajustado" en una consecuencia, no en una coincidencia.
#   * E^{(4)} se construye a partir de una tabla finita, y por eso sale matriz de transferencia.
#
# CUIDADO CON UN DETALLE: sumar t a eta_j puede romper que eta siga siendo una PARTICION (hay que
# mantener eta_1 >= eta_2 >= ...).  Se testa solo sobre los desplazamientos que dejan particion, y
# se dice cuantos son.
#
# LO QUE SE MIDE
#   P1  tau(eta) contra tau(eta') siempre que a = a' (mod t), con eta, eta' particiones distintas.
#       Se comparan VALOR y no solo anulacion.
#   P2  el numero de clases de residuo distintas que aparecen, y cuantos eta caen en cada una.
#   P3  el desplazamiento eta_j -> eta_j + t explicitamente: tau tiene que ser el mismo.
#
# CONTROLES
#   C0  tau por Freudenthal, que no sabe nada de residuos.
#   C1  SEÑUELO: comparar eta con a = a' (mod t-1).  Tiene que FALLAR, o el modulo no es t.
#   C2  n impreso siempre; y si alguna clase tiene un solo representante no prueba nada, asi que se
#       cuenta cuantas clases tienen 2 o mas.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage filter_is_periodic.sage

import json
import sys
from collections import defaultdict

_CH = {}
def car(typ, rk, mu):
    key = (typ, rk, tuple(int(v) for v in mu))
    if key not in _CH:
        W = WeylCharacterRing("%s%d" % (typ, rk))
        el = W(W.space().from_vector(vector([Integer(v) for v in mu])))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + int(mult)
        _CH[key] = d
    return _CH[key]


def tau(typ, rk, eta, t):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car(typ, rk, eta).items():
        s += mult * z ** (sum((i + 1) * wt[i] for i in range(rk)) % t)
    return QQ(s) if s in QQ else None


print("=" * 116)
print("¿ES EL FILTRO UNA FUNCION DE  a = eta + rho  MODULO t?")
print("=" * 116)
print("   t | G    | eta probados | clases mod t | clases con >=2 | tau constante en la clase |"
      " SEÑUELO mod t-1")
print("   " + "-" * 110)

RES = []
for t in range(3, 10):
    if t % 2 == 0:
        typ, rk = "C", (t - 2) // 2
        rho = [rk - j for j in range(rk)]
    else:
        typ, rk = "B", (t - 1) // 2
        rho = None                       # semientero: se trabaja con A = 2a
    if rk < 1:
        continue
    ETAS = []
    for k in range(0, 4 * t + 1):
        for e in Partitions(k, max_length=rk):
            ETAS.append(tuple(list(e) + [0] * (rk - len(e))))
    clases = defaultdict(list)
    clases_sen = defaultdict(list)
    val = {}
    for eta in ETAS:
        v = tau(typ, rk, eta, t)
        if v is None:
            continue
        val[eta] = v
        if typ == "C":
            a = tuple((int(eta[j]) + rho[j]) % t for j in range(rk))
            asen = tuple((int(eta[j]) + rk - j) % (t - 1) for j in range(rk))
        else:
            a = tuple((2 * int(eta[j]) + 2 * (rk - j - 1) + 1) % t for j in range(rk))
            asen = tuple((2 * int(eta[j]) + 2 * (rk - j - 1) + 1) % (t - 1) for j in range(rk))
        clases[a].append(eta)
        clases_sen[asen].append(eta)

    con2 = [c for c, L in clases.items() if len(L) >= 2]
    malas = [c for c in con2 if len(set(val[e] for e in clases[c])) > 1]
    con2s = [c for c, L in clases_sen.items() if len(L) >= 2]
    malas_s = [c for c in con2s if len(set(val[e] for e in clases_sen[c])) > 1]
    print("   %2d | %s%-2d  | %12d | %12d | %14d | %-25s | %s"
          % (t, typ, rk, len(val), len(clases), len(con2),
             "%d de %d" % (len(con2) - len(malas), len(con2)),
             "falla en %d de %d" % (len(malas_s), len(con2s)) if con2s else "sin clases"))
    if malas:
        c = malas[0]
        print("        CONTRAEJEMPLO en la clase %s: %s con tau = %s"
              % (str(c), str(clases[c][:3]), str([val[e] for e in clases[c][:3]])))
    sys.stdout.flush()
    RES.append({"t": int(t), "tipo": typ, "rango": int(rk), "n_eta": int(len(val)),
                "n_clases": int(len(clases)), "clases_con_2": int(len(con2)),
                "clases_constantes": int(len(con2) - len(malas)),
                "senuelo_clases": int(len(con2s)), "senuelo_malas": int(len(malas_s))})

# P3: el desplazamiento explicito eta_j -> eta_j + t
print("")
print("   P3  el desplazamiento  eta_j -> eta_j + t,  sobre los que siguen siendo particion:")
for t in (4, 6, 5, 7):
    if t % 2 == 0:
        typ, rk = "C", (t - 2) // 2
    else:
        typ, rk = "B", (t - 1) // 2
    if rk < 1:
        continue
    ok = tot = 0
    for k in range(0, 3 * t + 1):
        for e in Partitions(k, max_length=rk):
            eta = tuple(list(e) + [0] * (rk - len(e)))
            for j in range(rk):
                nu = list(eta)
                nu[j] += t
                if nu != sorted(nu, reverse=True):
                    continue
                tot += 1
                if tau(typ, rk, eta, t) == tau(typ, rk, tuple(nu), t):
                    ok += 1
    print("       t=%2d  %s%d :  %d de %d desplazamientos conservan tau" % (t, typ, rk, ok, tot))
    sys.stdout.flush()

json.dump(RES, open("filter_is_periodic_DUMP.json", "w"), indent=1)
print("")
print("=" * 116)
print("  LECTURA, escrita ANTES de correr: si tau es constante en cada clase mod t y el señuelo")
print("  mod t-1 falla, el filtro es una funcion sobre un conjunto FINITO, y el paso t del lado del")
print("  branching deja de ser una coincidencia con el t-ribbon: es la misma periodicidad.")
print("=" * 116)
print("DONE")
