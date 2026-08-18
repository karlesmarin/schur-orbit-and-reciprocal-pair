# -*- coding: utf-8 -*-
# ¿ES LA SATURACION PARTE DE LA HISTORIA DE conj:gcom, O ES OTRO PAPER? -- 15 de agosto de 2026.
#
# LA SOSPECHA.  La saturacion se mide sobre formas con prof FINITO, o sea Phi != 0: el estrato de
# arriba cancela pero el objeto NO se anula.  La Conjetura 8.44 vive en el otro lado, las formas que
# SI se anulan.  Son poblaciones DISJUNTAS, y los dos puentes propuestos hoy para cruzarlas estan
# muertos: el pelado con centro conservado (nuestro) y CANCEL gobernado por g_com (de la consulta).
#
# Asi que antes de invertir una semana en "que invariante es 4/24/18/42" hay que contestar:
#
#       ¿se mueve algo de lo que le importa a conj:gcom cuando CANCEL sube?
#
# A conj:gcom le importan exactamente tres cosas: C = min S + max S, el valor de empate tau, y cuanto
# de S es simetrico respecto del centro.  Si CANCEL sube de 2 a 4 mientras esas tres se quedan
# quietas, la saturacion es INDEPENDIENTE de la anulacion y es otro trabajo.
#
# COLUMNAS
#   CANCEL   niveles soportados que cancelan por delante
#   C, tau   las dos constantes; C==tau es la conclusion del Corolario 8.31 bajo anulacion
#   simC     fraccion de S con su reflejo C-v tambien en S    (1.0 = S entero simetrico)
#   simTau   lo mismo respecto de tau
#   cero     Phi == 0
#
# CONTROL.  Se corre tambien una trayectoria que ARRANQUE en una forma que se anula, para ver si al
# separar los extremos sigue anulandose o sale del locus.  Si la anulacion se pierde de inmediato,
# las dos poblaciones no solo son disjuntas: no se tocan ni por continuidad.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python cancel_vs_vanishing.py

import os

from second_stratum import setup
from peel_zero import phi_zero, C_and_tau

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]


def S_of(beta, t):
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    return sorted({v for k in E for v in Cd[k]})


def sim(S, c):
    if c is None or not S:
        return 0.0
    s = set(S)
    return sum(1 for v in S if (c - v) in s) / float(len(S))


CASOS = [(4, 2, 26, (18, 17, 11, 8, 7, 6, 1, 0), 10),
         (4, 2, 28, (18, 17, 11, 8, 7, 6, 1, 0), 10)]

print("=" * 104)
print("¿SE MUEVE ALGO DE conj:gcom CUANDO CANCEL SUBE?")
print("=" * 104)
for (t, r, s, seed, J) in CASOS:
    S0 = S_of(seed, t)
    hi, lo = S0[-1], S0[0]
    print("\n--- t=%d r=%d  paso s=%d  semilla %s" % (t, r, s, seed))
    print("      j |    W | CANCEL |    C |  tau | C==tau | simC | simTau | cero")
    filas = []
    for j in range(J + 1):
        b = tuple(sorted([(x + s * j if x == hi else (x - s * j if x == lo else x)) for x in seed],
                         reverse=True))
        rec = probe(b, t, r)
        if rec is None or not rec['surv'] or rec['prof'] is None:
            break
        S = S_of(b, t)
        ct = C_and_tau(b, t, r)
        C, tau = (ct if ct else (None, None))
        z = phi_zero(b, t, r)
        filas.append((rec['vac_cancelan'], C, tau, C == tau, round(sim(S, C), 3),
                      round(sim(S, tau), 3), z))
        print("    %3d | %4d | %6s | %4s | %4s | %6s | %4.2f | %6.2f | %s"
              % (j, b[0] - b[-1], rec['vac_cancelan'], C, tau, C == tau,
                 sim(S, C), sim(S, tau), z))
    if filas:
        cans = [f[0] for f in filas]
        movio = len(set(cans)) > 1
        quietos = []
        for k, nom in ((3, 'C==tau'), (4, 'simC'), (5, 'simTau'), (6, 'cero')):
            quietos.append((nom, len(set(f[k] for f in filas)) == 1))
        print("    -> CANCEL %d..%d (%s).  Quietos: %s"
              % (min(cans), max(cans), 'SE MUEVE' if movio else 'quieto',
                 ', '.join(n for n, q in quietos if q) or 'ninguno'))
        if movio and all(q for _, q in quietos):
            print("       VEREDICTO: CANCEL se mueve y NADA de conj:gcom se mueve"
                  "  ->  poblaciones INDEPENDIENTES")

print()
print("=" * 104)
print("  Si CANCEL sube mientras C, tau, la simetria y la anulacion no se mueven, la saturacion")
print("  no es parte de la historia de conj:gcom: es otro trabajo.")
print("=" * 104)
