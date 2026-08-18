# -*- coding: utf-8 -*-
# FALSEAR LA LEY DEL EXTREMO FUGITIVO.  13 de agosto de 2026.
#
# LA LEY, tal como queda escrita despues de defect_cone.py (donde mi prediccion del CONO fallo: de 118
# vectores concentricos solo 3 generan familia, y los 3 son el mismo mas una constante, o sea UNO mod
# traslacion).  Restandole la traslacion, el generador se desnuda:
#
#     LEY DEL EXTREMO FUGITIVO.  Si beta es superviviente, tambien lo es  beta + j*(2,0,...,0,-2)
#     -- los N-2 puntos interiores CONGELADOS, los dos extremos separandose 2 por lado y por paso --
#     y  prof(beta_j) = prof(beta) + 2j,  mientras el numero de estratos que CANCELAN se satura.
#
# Medido hasta ahora en DOS semillas de t=4 r=2.  Dos semillas no es una ley: es una anecdota con dos
# puntos.  Este guion existe para romperla, y tiene cuatro maneras de hacerlo:
#
#   N1  UNIVERSALIDAD.  Se aplica a TODOS los supervivientes que salen del barrido exhaustivo de
#       t=4 r=2 hasta W=30 -- no a los dos elegidos.  Si la ley es de los supervivientes, ninguno
#       falla; si falla alguno, la condicion que lo separa es el objeto y hay que enseñarla.
#   N2  EL PASO.  Por que 2 y no 1, 3, 4?  Y por que simetrico?  Se prueban los pasos s=1..4 y las
#       variantes asimetricas (solo arriba, solo abajo).  Si s=1 tambien valiera, la anchura crecería
#       de 2 en 2 y la pendiente de la recta NO seria 1/2: la pendiente saldria del paso, no del objeto.
#   N3  SEÑUELO.  Mover un punto INTERIOR en vez de los extremos.  Si tambien generase familias, no
#       seria "el extremo fugitivo" sino cualquier deformacion, y la ley perderia su contenido.
#   N4  OTRAS (t, r) -- la pregunta original.  Los supervivientes de t=6 r=2, t=6 r=3, t=8 r=3, t=6 r=4
#       se someten al mismo generador.  La PENDIENTE dK/dW de cada configuracion sale de aqui sin
#       barrido exhaustivo: es 2/(2s) = 1/s con s el paso minimo que funciona.
#
# UNA PREDICCION, escrita antes de correr: los saltos de K medidos por k_vs_m.py estan espaciados 4 en
# t=4 r=2 y t=6 r=3, y espaciados 2 en t=8 r=3 y t=6 r=4.  Si la pendiente es 1/s, eso dice s=2 en las
# dos primeras y s=1 en las dos ultimas -- y t=2r se cumple exactamente en las dos primeras.
# PREDIGO: s = 2 si t = 2r, s = 1 si no.  N4 la confirma o la mata.
#
# EL CRITERIO NO SE REESCRIBE: probe() sale de ejecutar el preambulo de k_vs_m.py, mismos bytes.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python falsify_law.py

import itertools
import json
import os
import sys
import time
from collections import Counter, defaultdict

from survivors_wide import scan

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
assert "def probe(" in _head and "def shapes_of_width(" in _head, "k_vs_m.py cambio de forma"
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]
shapes_of_width = _ns["shapes_of_width"]

JMAX = 10
OUT_JSON = "falsify_law_RESULT.json"


def gen(N, s, arriba=True, abajo=True, pos=None):
    """el vector generador.  pos != None mueve ESE indice interior en vez de los extremos (senuelo)."""
    v = [0] * N
    if pos is not None:
        v[pos] = s
        return tuple(v)
    if arriba:
        v[0] = s
    if abajo:
        v[N - 1] = -s
    return tuple(v)


def itera(beta0, v, t, r, jmax=JMAX):
    filas = []
    for j in range(jmax + 1):
        beta = tuple(b + j * x for b, x in zip(beta0, v))
        if any(beta[i] <= beta[i + 1] for i in range(len(beta) - 1)):
            break
        rec = probe(beta, t, r)
        if rec is None or not rec['surv'] or rec['prof'] is None:
            break
        filas.append((j, beta[0] - beta[-1], rec['prof'], rec['vac_cancelan']))
    return filas


def ley_ok(filas, jmax=JMAX):
    """la ley exige: llega a jmax Y prof crece exactamente +2 por paso."""
    if len(filas) < jmax + 1:
        return False, "se rompe en j=%d" % len(filas)
    p = [f[2] for f in filas]
    d = [p[i] - p[i - 1] for i in range(1, len(p))]
    if set(d) != {2}:
        return False, "prof no crece +2: %s" % p
    return True, "ok"


def supervivientes(t, r, Wmax, Wmin=None):
    """todos los betas supervivientes con anchura en [Wmin, Wmax], barrido exhaustivo."""
    N = t + 2 * r
    out = []
    for W in range(Wmin or (N - 1), Wmax + 1):
        for beta in shapes_of_width(W, N):
            rec = probe(beta, t, r)
            if rec is None or not rec['surv'] or rec['prof'] is None:
                continue
            out.append((beta, rec['prof'], rec['vac_cancelan']))
    return out


# ===================================================================== C0 ========================
print("=" * 112)
print("C0  ACEPTACION -- fatal")
print("=" * 112)
print("")
bad = 0
for (t, r, M) in [(4, 2, 15), (6, 2, 17), (6, 3, 18)]:
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
    print("  C0a  probe() == scan()  t=%d r=%d M=%d : %s" % (t, r, M, "ok" if ok else "*** FALLA ***"))
f = itera((18, 17, 11, 8, 7, 6, 1, 0), gen(8, 2), 4, 2)
o, m = ley_ok(f)
bad += not o
print("  C0b  la familia minima ya publicada se re-mide: %d miembros, %s" % (len(f), m))
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")
RES = {}

# ===================================================================== N1 ========================
print("")
print("=" * 112)
print("N1  UNIVERSALIDAD -- la ley contra TODOS los supervivientes de t=4 r=2 hasta W=30")
print("=" * 112)
print("")
t0 = time.time()
SV = supervivientes(4, 2, 30)
print("  %d supervivientes en W = 7..30  (%.0f s)" % (len(SV), time.time() - t0))
buenos, malos = [], []
for (beta, p, c) in SV:
    filas = itera(beta, gen(8, 2), 4, 2)
    o, m = ley_ok(filas)
    (buenos if o else malos).append((beta, p, c, m, len(filas)))
print("  la ley (prof +2 por paso, %d pasos) se cumple en %d de %d  (%.1f%%)"
      % (JMAX, len(buenos), len(SV), 100.0 * len(buenos) / len(SV)))
if malos:
    print("  LOS QUE FALLAN, uno a uno (hasta 20):")
    for (beta, p, c, m, n) in malos[:20]:
        print("     %-42s prof=%-3d canc=%-2d  %s" % (str(beta), p, c, m))
    if len(malos) > 20:
        print("     ... y %d mas" % (len(malos) - 20))
    dW = Counter(b[0][0] - b[0][-1] for b in malos)
    print("  los fallos por anchura: %s" % dict(sorted(dW.items())))
    dP = Counter(b[1] for b in malos)
    print("  los fallos por profundidad de partida: %s" % dict(sorted(dP.items())))
else:
    print("  NINGUN superviviente falla.  La ley no es de dos semillas.")
RES['N1'] = dict(total=len(SV), ok=len(buenos), fallan=len(malos),
                 ejemplos_fallo=[[list(b[0]), b[1], b[2], b[3]] for b in malos[:40]])
sys.stdout.flush()

# ===================================================================== N2 ========================
print("")
print("=" * 112)
print("N2  EL PASO -- por que 2, y por que simetrico")
print("=" * 112)
print("")
print("     generador                        de %d semillas, cuantas dan familia de %d pasos" % (min(len(SV), 60), JMAX))
print("  " + "-" * 92)
muestra = SV[:60]
N2 = {}
variantes = ([("simetrico  s=%d  (%d,0..0,-%d)" % (s, s, s), gen(8, s)) for s in (1, 2, 3, 4)]
             + [("solo ARRIBA s=%d" % s, gen(8, s, abajo=False)) for s in (1, 2)]
             + [("solo ABAJO  s=%d" % s, gen(8, s, arriba=False)) for s in (1, 2)])
for (nombre, v) in variantes:
    n = 0
    for (beta, p, c) in muestra:
        o, _ = ley_ok(itera(beta, v, 4, 2))
        n += o
    N2[nombre] = n
    print("     %-34s %d de %d   %s" % (nombre, n, len(muestra),
                                        "<-- funciona" if n else ""))
    sys.stdout.flush()
RES['N2'] = N2

# ===================================================================== N3 ========================
print("")
print("=" * 112)
print("N3  SEÑUELO -- mover un punto INTERIOR en vez de los extremos")
print("=" * 112)
print("")
N3 = {}
for pos in range(1, 7):
    for s in (2, -2):
        v = gen(8, s, pos=pos)
        n = sum(ley_ok(itera(beta, v, 4, 2))[0] for (beta, p, c) in muestra)
        N3["pos%d_%+d" % (pos, s)] = n
        print("     mover beta[%d] en %+d : %d de %d %s"
              % (pos, s, n, len(muestra), "*** TAMBIEN FUNCIONA -- la ley pierde contenido ***" if n else ""))
    sys.stdout.flush()
RES['N3'] = N3
if max(N3.values()) == 0:
    print("     ningun senuelo interior genera familia.  Son LOS EXTREMOS, no una deformacion cualquiera.")

# ===================================================================== N4 ========================
print("")
print("=" * 112)
print("N4  LAS OTRAS (t, r) -- la pregunta original, sin barrido exhaustivo largo")
print("=" * 112)
print("")
print("     PREDICHO en la cabecera, antes de correr:  s = 2 si t = 2r,  s = 1 si no.")
print("")
CFG = [(4, 2, 26), (6, 2, 26), (6, 3, 21), (8, 3, 22), (6, 4, 21), (8, 2, 24), (10, 2, 23)]
print("     t   r   N | supervivientes | s=1  s=2  s=3  s=4 | s minimo | pendiente 1/s | t=2r?")
print("  " + "-" * 104)
N4 = {}
for (t, r, Wmax) in CFG:
    N = t + 2 * r
    t0 = time.time()
    sv = supervivientes(t, r, Wmax)
    if not sv:
        print("    %2d  %2d  %2d |  %13s | -    -    -    -   | -        | -             | %s"
              % (t, r, N, "0 (hasta W=%d)" % Wmax, "si" if t == 2 * r else "no"))
        N4["%d_%d" % (t, r)] = dict(n=0)
        continue
    m = sv[:60]
    cuenta = {}
    for s in (1, 2, 3, 4):
        cuenta[s] = sum(ley_ok(itera(beta, gen(N, s), t, r))[0] for (beta, p, c) in m)
    smin = min([s for s in cuenta if cuenta[s] > 0], default=None)
    print("    %2d  %2d  %2d |  %5d (%3d us) | %-4d %-4d %-4d %-4d | %-8s | %-13s | %s"
          % (t, r, N, len(sv), len(m), cuenta[1], cuenta[2], cuenta[3], cuenta[4],
             smin if smin else "ninguno", ("1/%d" % smin) if smin else "-",
             "SI" if t == 2 * r else "no"))
    N4["%d_%d" % (t, r)] = dict(n=len(sv), muestra=len(m), cuenta=cuenta, smin=smin,
                                t_eq_2r=(t == 2 * r))
    sys.stdout.flush()
RES['N4'] = N4

# ===================================================================== N5 ========================
print("")
print("=" * 112)
print("N5  VEREDICTO")
print("=" * 112)
print("")
print("  N1  universalidad : %d de %d supervivientes de t=4 r=2 cumplen la ley"
      % (RES['N1']['ok'], RES['N1']['total']))
print("  N2  el paso       : %s" % ", ".join("%s=%d" % (k.split()[0] + k.split()[1], v)
                                             for k, v in RES['N2'].items() if v))
print("  N3  senuelo       : %s" % ("0 de todos -- son los extremos"
                                    if max(RES['N3'].values()) == 0 else "ALGUNO FUNCIONA, mirar N3"))
vivos = {k: v for k, v in N4.items() if v.get('smin')}
print("  N4  otras (t,r)   : familia encontrada en %d de %d configuraciones" % (len(vivos), len(CFG)))
pred_ok = all((v['smin'] == 2) == v['t_eq_2r'] for v in vivos.values())
print("      la prediccion 's=2 si y solo si t=2r' : %s"
      % ("ACIERTA en las %d" % len(vivos) if pred_ok else "FALLA -- y se dice"))
for k, v in sorted(vivos.items()):
    print("        t=%s r=%s : s = %d, pendiente 1/%d %s"
          % (tuple(k.split("_")) + (v['smin'], v['smin'],
                                    "(t=2r)" if v['t_eq_2r'] else "")))
print("")
print("  ALCANCE: los supervivientes salen de barridos exhaustivos con tope de RELOJ, asi que")
print("  'ninguno falla' significa 'ninguno de los que hay hasta ese tope'.  Lo que NO depende del")
print("  tope es la iteracion: cada familia llega a anchuras que ningun barrido alcanza.")
print("")
json.dump(RES, open(OUT_JSON, "w"), indent=1)
print("DONE")
