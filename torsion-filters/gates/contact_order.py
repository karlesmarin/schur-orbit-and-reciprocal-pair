# -*- coding: utf-8 -*-
# LA PROFUNDIDAD ES ORDEN DE CONTACTO CON EL LUGAR CONCENTRICO.  13 de agosto de 2026.
#
# DE DONDE SALE.  witness_family.py saco los testigos de profundidad maxima de t=4 r=2 en las anchuras
# 18..38, y todos tienen la misma anatomia: de los N/2 = 4 pares  beta_i + beta_{N-1-i},  TRES valen
# exactamente C = W y UNO falla -- y falla SIEMPRE POR 4, en las cinco anchuras.  O sea el defecto de
# concentricidad se queda FIJO mientras la escala W crece.
#
#     (38, 37, 23, 18, 16, 15,  1, 0)    37+1=38  23+15=38  18+16=34   defecto 4
#     (34, 33, 23, 16, 14, 11,  1, 0)    33+1=34  23+11=34  16+14=30   defecto 4
#     (18, 17, 11,  8,  7,  6,  1, 0)    17+1=18  11+ 6=17   8+ 7=15   defecto 1+3 = 4
#
# LA LECTURA.  Phi_t se anula EXACTAMENTE sobre el lugar concentrico ((i) y (ii), teorema del 12).  Un
# superviviente profundo es un punto que se ACERCA a ese lugar sin estar en el, y la profundidad es el
# ORDEN DE CONTACTO.  Con defecto fijo y escala creciente el punto es relativamente cada vez mas
# concentrico, luego el contacto sube y NO PUEDE haber cota.  Eso no mata el enunciado: lo reubica.
# La pregunta no era "que profundidad" sino "a que distancia del lugar concentrico".
#
# LO QUE SE MIDE AQUI, y las tres cosas pueden salir que no:
#   N1  LA FAMILIA.  Los testigos de anchuras consecutivas se generan sumando UN VECTOR FIJO
#       v = (4,2,2,2,2,2,2,0):  beta_{j+1} = beta_j + v.  Verificado a ojo en 18 -> 22 -> 26 -> 30.
#       Aqui se itera hasta j = 24 (anchura 114), MUY por encima de cualquier barrido exhaustivo
#       posible, y se comprueba que cada beta_j SIGUE siendo superviviente y que prof = 6 + 2j.
#       Si aguanta, "para todo K hay un testigo de profundidad > K" deja de ser una medida y pasa a
#       ser una FAMILIA EXHIBIDA -- que es lo que se puede intentar demostrar en cerrado.
#   N2  EL DEFECTO.  Para cada miembro: el multiconjunto de defectos de par, su suma, y prof.  La
#       hipotesis es defecto CONSTANTE con prof creciente.  Si el defecto creciera con j, la lectura
#       "orden de contacto" se cae y hay que decirlo.
#   N3  LA OTRA FAMILIA, la que cancela mas: (26,17,15,12,11,10,9,0) -> (30,19,17,14,13,12,11,0) ->
#       (38,23,21,18,17,16,15,0) tenia 3, 3, 4 peldanos que CANCELAN.  Si tambien se extiende, la cota
#       en la escalera CON SOPORTE muere tambien -- mas despacio, pero muere -- y se dice.
#
# EL CRITERIO NO SE REESCRIBE: probe() sale de ejecutar el preambulo de k_vs_m.py (mismos bytes).
# C0 lo re-firma contra scan() de survivors_wide.py, y ademas re-mide los testigos publicados por
# witness_family.py: si probe() aqui no reprodujera SUS numeros, nada de esto valdria.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python contact_order.py

import itertools
import json
import os
import sys
from collections import Counter

from survivors_wide import scan

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
assert "def probe(" in _head, "k_vs_m.py cambio de forma"
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]

T_, R_ = 4, 2
N_ = T_ + 2 * R_
V = (4, 2, 2, 2, 2, 2, 2, 0)          # el generador, leido de los testigos
JMAX = 24

# los testigos publicados en witness_family_OUT.txt, copiados a mano para re-medirlos
TESTIGOS_PUB = [
    ((18, 17, 11, 8, 7, 6, 1, 0), 6, 2),
    ((22, 19, 13, 10, 9, 8, 3, 0), 8, 2),
    ((26, 25, 15, 12, 11, 10, 1, 0), 10, 2),
    ((26, 17, 15, 12, 11, 10, 9, 0), 10, 3),
    ((30, 29, 19, 14, 12, 11, 1, 0), 12, 2),
    ((30, 19, 17, 14, 13, 12, 11, 0), 12, 3),
    ((34, 33, 23, 16, 14, 11, 1, 0), 14, 2),
    ((38, 23, 21, 18, 17, 16, 15, 0), 16, 4),
]


def defecto(beta):
    """(lista de defectos beta_i + beta_{N-1-i} - C, con C = max+min ; su suma en valor absoluto)."""
    C = beta[0] + beta[-1]
    d = [beta[i] + beta[len(beta) - 1 - i] - C for i in range(len(beta) // 2)]
    return d, sum(abs(x) for x in d)


def add(beta, j):
    return tuple(b + j * v for b, v in zip(beta, V))


# ===================================================================== C0 ========================
print("=" * 108)
print("C0  ACEPTACION -- probe() es el de k_vs_m.py byte a byte, y aqui se le exige DOS cosas")
print("=" * 108)
print("")
bad = 0
print("  C0a  reproduce scan() de survivors_wide.py")
for (t, r, M) in [(4, 2, 15), (6, 3, 18)]:
    n_ref, cont_ref, sv_ref = scan(t, r, M)
    mine, mb = Counter(), []
    for comb in itertools.combinations(range(M + 1), t + 2 * r):
        beta = tuple(sorted(comb, reverse=True))
        rec = probe(beta, t, r, deep=False)
        if rec is None:
            continue
        mine[(rec['e'] == t, rec['surv'])] += 1
        if rec['surv']:
            mb.append(beta)
    ok = (sum(mine.values()) == n_ref and mine == cont_ref
          and sorted(mb) == sorted(x['beta'] for x in sv_ref))
    bad += not ok
    print("       t=%d r=%d M=%d : objetivo %d/%d, betas %d/%d   %s"
          % (t, r, M, sum(mine.values()), n_ref, len(mb), len(sv_ref), "ok" if ok else "*** FALLA ***"))

print("")
print("  C0b  re-mide los TESTIGOS que publico witness_family.py: prof y peldanos que cancelan")
for (beta, p_pub, c_pub) in TESTIGOS_PUB:
    rec = probe(beta, T_, R_)
    ok = rec is not None and rec['surv'] and rec['prof'] == p_pub and rec['vac_cancelan'] == c_pub
    bad += not ok
    print("       %-38s prof %s/%d, cancelan %s/%d   %s"
          % (str(beta), rec['prof'] if rec else "-", p_pub,
             rec['vac_cancelan'] if rec else "-", c_pub, "ok" if ok else "*** FALLA ***"))
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")

# ===================================================================== N1, N2 ====================
print("")
print("=" * 108)
print("N1  LA FAMILIA  beta_j = beta_0 + j*(4,2,2,2,2,2,2,0)   |   N2  su defecto de concentricidad")
print("=" * 108)
SEMILLAS = [
    # A y B son los DOS testigos de W=18 con beta_7 = 0; B es la rama que mas cancela, y se obtiene
    # restando 2v al testigo (26,17,15,12,11,10,9,0), o sea la familia ya estaba en la tabla.
    ("A  testigo de W=18, rama alta", (18, 17, 11, 8, 7, 6, 1, 0)),
    ("B  testigo de W=18, rama que mas cancela", (18, 13, 11, 8, 7, 6, 5, 0)),
    # C es el reflejo de A bajo x -> C - x: los testigos salen en parejas reflejadas, y hay que
    # comprobar que la familia no privilegia una de las dos mitades.
    ("C  el REFLEJO de A bajo x -> W - x", (18, 17, 12, 11, 10, 7, 1, 0)),
]
RES = {}
for (nombre, beta0) in SEMILLAS:
    print("")
    print("  %s   semilla %s" % (nombre, str(beta0)))
    print("     j    W   beta                                          objetivo surv  prof  cancelan  sinsop  defecto")
    print("  " + "-" * 118)
    filas = []
    roto = None
    for j in range(JMAX + 1):
        beta = add(beta0, j)
        W = beta[0] - beta[-1]
        rec = probe(beta, T_, R_)
        d, dt = defecto(beta)
        if rec is None:
            estado = "NO objetivo"
            print("  %4d %4d   %-44s %s" % (j, W, str(beta), estado))
            roto = roto or (j, "deja de ser poblacion objetivo")
            break
        if not rec['surv'] or rec['prof'] is None:
            print("  %4d %4d   %-44s si       NO   (Delta != 0 o Phi_t == 0)" % (j, W, str(beta)))
            roto = roto or (j, "deja de sobrevivir")
            break
        print("  %4d %4d   %-44s si       si  %5d %9d %7d   %s  (suma %d)"
              % (j, W, str(beta), rec['prof'], rec['vac_cancelan'], rec['vac_sin_soporte'],
                 str(d), dt))
        filas.append(dict(j=j, W=W, beta=list(beta), prof=rec['prof'],
                          cancelan=rec['vac_cancelan'], sinsop=rec['vac_sin_soporte'],
                          defecto=d, defecto_total=dt))
        sys.stdout.flush()
    RES[nombre] = dict(semilla=list(beta0), filas=filas, roto=roto)
    if filas:
        profs = [f['prof'] for f in filas]
        cans = [f['cancelan'] for f in filas]
        defs = {f['defecto_total'] for f in filas}
        print("     -> %d miembros, prof de %d a %d, cancelan de %d a %d, defecto total %s"
              % (len(filas), profs[0], profs[-1], min(cans), max(cans),
                 "CONSTANTE = %d" % defs.pop() if len(defs) == 1 else "VARIA: %s" % sorted(defs)))
        pred = [profs[0] + 2 * i for i in range(len(profs))]
        print("     -> prof = %d + 2j : %s" % (profs[0], "ACIERTA los %d" % len(profs)
                                               if profs == pred else "FALLA -> %s" % profs))
    if roto:
        print("     -> la familia SE ROMPE en j = %d (%s), y se dice" % roto)

json.dump(RES, open("contact_order_FAMILY.json", "w"), indent=1)

# ===================================================================== N3 ========================
print("")
print("=" * 108)
print("N3  VEREDICTO")
print("=" * 108)
print("")
for nombre in RES:
    f = RES[nombre]['filas']
    if not f:
        print("  %s : sin miembros" % nombre)
        continue
    print("  %s" % nombre)
    print("     anchura %d -> %d, profundidad %d -> %d, peldanos que CANCELAN %d -> %d"
          % (f[0]['W'], f[-1]['W'], f[0]['prof'], f[-1]['prof'], f[0]['cancelan'], f[-1]['cancelan']))
    print("     defecto de concentricidad: %s"
          % ("CONSTANTE = %d en los %d miembros" % (f[0]['defecto_total'], len(f))
             if len({x['defecto_total'] for x in f}) == 1 else "VARIA"))
    print("     %s" % (RES[nombre]['roto'] and "SE ROMPE en j=%d (%s)" % RES[nombre]['roto']
                       or "aguanta los %d miembros SIN ROMPERSE" % len(f)))
print("")
print("  LO QUE ESTO ES, y lo que NO es:")
print("     Si una familia aguanta con prof = prof_0 + 2j, entonces 'para todo K existe un beta con")
print("     profundidad > K' esta EXHIBIDO, no medido: no depende de ningun tope de barrido.  Y si el")
print("     defecto de concentricidad se queda fijo mientras W crece, la profundidad es el ORDEN DE")
print("     CONTACTO con el lugar donde Phi_t se anula, no una cantidad combinatoria suelta.")
print("     Lo que NO es: una prueba.  Es un testigo explicito, que es la materia prima de una.")
print("")
print("DONE")
