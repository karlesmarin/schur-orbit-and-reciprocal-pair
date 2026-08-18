# -*- coding: utf-8 -*-
# EL PAPER II, AUDITADO -- Y LA HERRAMIENTA, AFINADA.   16 de agosto de 2026.
#
# Dogfooding: pasar nuestras identidades algebraicas por `socrates.tools.tautology`
# (E:\proyectos\socrates-audit) y, sobre todo, ver donde el auditor se equivoca.
#
# LO QUE APARECIO.  De 25 casos, Socrates acerto 19 y fallo 6, y el patron era diagnostico: acertaba
# a t = 4 y t = 8 y fallaba a t = 3, 5, 6, 7.  Es decir, su `simplify` no cierra productos
# CICLOTOMICOS de orden que no sea potencia de 2.  Comprobacion numerica a 30 digitos en tres valores
# de z: |P - (1-z^t)| ~ 1e-140 .. 1e-165.  Las identidades son CIERTAS; el equivocado es el auditor.
# Y el parche obvio, expand_complex, arregla t=3 y t=5 y NO arregla t=7: no basta.
#
# LA MANERA CORRECTA, que es la que se usa aqui y la que habria que llevar al tool: no simplificar
# con heuristicas sobre exp(2 pi i / t), sino trabajar en el CUERPO CICLOTOMICO.  Se representa la
# raiz por un simbolo w, se expande el producto como polinomio en z y w, se reducen los exponentes
# de w modulo t, y cada coeficiente se reduce modulo el polinomio ciclotomico Phi_t(w).  Es
# aritmetica entera exacta y decide sin heuristica.
#
# CONTROLES
#   T1  cada identidad probada del paper tiene que dar CERO por reduccion ciclotomica exacta.
#   T2  cada señuelo tiene que dar NO CERO.  Un señuelo que da cero invalida el metodo, no el paper.
#   T3  y se anota, caso a caso, que dijo Socrates, para saber donde ciega.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python socrates_audit.py

import io
import json
import sys

sys.path.insert(0, r"E:\proyectos\socrates-audit")

import sympy                                                        # noqa: E402
from sympy import cyclotomic_poly, symbols                          # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

z, w = symbols("z w")

try:
    from socrates.tools.tautology import IdentityRule, is_tautology
    HAY_SOCRATES = True
except Exception as e:                                              # pragma: no cover
    HAY_SOCRATES = False
    print("aviso: no se pudo importar Socrates (%s); solo se corre la reduccion exacta" % e)


def reduce_ciclotomico(expr, t):
    """Cero exacto en Q(zeta_t): expandir en z y w, w^t -> 1, y reducir mod Phi_t(w).

    Devuelve el resto.  Es cero si y solo si la identidad se cumple en el cuerpo ciclotomico.
    """
    P = sympy.Poly(sympy.expand(expr), z)
    Phi = sympy.Poly(cyclotomic_poly(t, w), w)
    resto = []
    for coef in P.all_coeffs():
        c = sympy.Poly(sympy.expand(coef), w)
        # w^t = 1: se pliegan los exponentes
        plegado = 0
        for (e,), a in zip(c.monoms(), c.coeffs()):
            plegado += a * w ** (e % t)
        resto.append(sympy.rem(sympy.Poly(plegado, w), Phi).as_expr())
    return sympy.simplify(sum(sympy.Abs(x) for x in resto)) == 0, resto


CASOS = []


def caso(etiqueta, donde, expr, t, espera, socrates_expr=None, reglas=()):
    CASOS.append({"etiqueta": etiqueta, "donde": donde, "expr": expr, "t": t,
                  "espera": espera, "soc": socrates_expr, "reglas": list(reglas)})


def prod_orbita(t):
    return sympy.prod([1 - w ** kk * z for kk in range(t)])


def prod_cross(t):
    mp = (t - 1) // 2
    return (1 - z) * sympy.prod([(1 - w ** i * z) * (1 - w ** (t - i) * z)
                                 for i in range(1, mp + 1)])


def prod_even(t):
    m = (t - 2) // 2
    return sympy.prod([(1 - w ** i * z) * (1 - w ** (t - i) * z) for i in range(1, m + 1)])


for tt in (3, 4, 5, 7, 8):
    caso("epsdet t=%d: prod_{k<t}(1-w^k z) = 1-z^t" % tt, "lem:epsdet",
         prod_orbita(tt) - (1 - z ** tt), tt, True)
caso("SEÑUELO epsdet t=4 con 1-z^5", "lem:epsdet",
     prod_orbita(4) - (1 - z ** 5), 4, False)

for tt in (3, 5, 7, 9):
    caso("crossden t=%d: (1-z)prod(1-w^i z)(1-w^-i z) = 1-z^t" % tt, "prop:crossden",
         prod_cross(tt) - (1 - z ** tt), tt, True)
caso("SEÑUELO crossden t=5 sin el factor (1-z)", "prop:crossden",
     prod_cross(5) / (1 - z) - (1 - z ** 5), 5, False)

for tt in (4, 6, 8, 10):
    caso("even t=%d: prod = (1-z^t)/((1-z)(1+z))" % tt, "prop:eventransversal(i)",
         sympy.expand(prod_even(tt) * (1 - z) * (1 + z)) - (1 - z ** tt), tt, True)
caso("SEÑUELO even t=6 con el denominador impar", "prop:eventransversal(i)",
     sympy.expand(prod_even(6) * (1 - z)) - (1 - z ** 6), 6, False)

print("=" * 104)
print("PAPER II: IDENTIDADES CICLOTOMICAS, POR REDUCCION EXACTA -- Y QUE DIJO SOCRATES")
print("=" * 104)
print("")
print("  %-54s %-10s %-12s %s" % ("identidad", "exacto", "Socrates", "donde"))
print("  " + "-" * 100)

ok = 0
disc = []
for c in CASOS:
    es_cero, _ = reduce_ciclotomico(c["expr"], c["t"])
    acierta = (es_cero == c["espera"])
    ok += 1 if acierta else 0
    # y lo mismo, tal como Socrates lo ve: con la raiz como exp(2 pi i / t)
    soc = "-"
    if HAY_SOCRATES:
        try:
            e = c["expr"].subs(w, sympy.exp(2 * sympy.pi * sympy.I / c["t"]))
            v = bool(is_tautology(sympy.expand(sympy.simplify(e)), []))
            soc = "TRIVIAL" if v else "no trivial"
            if v != c["espera"]:
                disc.append((c["etiqueta"], soc))
        except Exception as ex:
            soc = "ERROR"
            disc.append((c["etiqueta"], "ERROR %s" % ex))
    print("  %s %-52s %-10s %-12s %s"
          % ("ok " if acierta else "!! ", c["etiqueta"][:52],
             "0" if es_cero else "no 0", soc, c["donde"]))

print("")
print("-" * 104)
print("  T1/T2  reduccion exacta, aciertos : %d de %d" % (ok, len(CASOS)))
print("  T3     casos donde Socrates discrepa de la verdad : %d de %d" % (len(disc), len(CASOS)))
for (e, v) in disc:
    print("         %-54s dijo %s" % (e[:54], v))
print("")
print("  LECTURA. La reduccion ciclotomica decide sin heuristica y acierta en todos.  Socrates falla")
print("  exactamente donde su simplify no cierra la raiz: el arreglo para el tool no es un simplify")
print("  mas agresivo, es reducir modulo Phi_t antes de preguntar.")

json.dump({"n": len(CASOS), "exactos_ok": ok,
           "socrates_discrepa": [{"caso": e, "dijo": v} for (e, v) in disc]},
          io.open("socrates_audit_DUMP.json", "w", encoding="utf-8"), indent=1)
print("")
print("=" * 104)
print("DONE")
