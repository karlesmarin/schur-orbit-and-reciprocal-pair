# -*- coding: utf-8 -*-
# ============================================================================================
#  C = tau, REDUCIDO A UN ENUNCIADO SOBRE DOS ELEMENTOS.  14 de agosto de 2026.
#
#  DE DONDE VIENE.  c_eq_tau.py dejo esto (su _OUT.txt): sigma_V(K) = K esta PROBADO (P1), luego
#  min K + max K = V; y C == tau se cumple 2522 de 2522, pero H* ("min S y max S estan en K") FALLA
#  en 919.  El N5 de aquel guion reformulo la brecha como "los dos recortes se compensan", medido
#  0 pares con recorte arriba != recorte abajo.  Eso es CIERTO pero es una REESCRITURA de la tesis,
#  no una reduccion: por eso su Contingencia salio "SIN VALOR, la hipotesis se cumple SIEMPRE".
#
#  LA REDUCCION DE VERDAD, y es de dos lineas.  sigma_V(v) = V - v invierte el orden.  Entonces:
#
#      sigma_V(max S) en S   =>   sigma_V(max S) >= min S   =>   max S + min S <= V
#      sigma_V(min S) en S   =>   sigma_V(min S) <= max S   =>   max S + min S >= V
#      luego  C = max S + min S = V = tau.
#
#  Y las dos hipotesis son UNA SOLA, porque  S = K  u  g_com  u  {x1,x2,y1,y2}  con
#  sigma_V(K) = K (PROBADO) y sigma_V{x1,x2} = {x1,x2}, sigma_V{y1,y2} = {y1,y2} (el empate,
#  PROBADO: x1+x2 = y1+y2 = V).  El unico trozo de S que sigma_V puede sacar fuera es g_com.
#  Por tanto:
#
#      (A1)   max S no esta en g_com   Y   min S no esta en g_com        =>   C = tau
#
#  o sea: NINGUN EXTREMO DE S LO DESCARTAN LOS DOS MAXIMIZADORES A LA VEZ.  Eso SI es un enunciado
#  sobre dos elementos, y es sobre la OPTIMALIDAD de la transversal, no sobre la simetria.
#
#  LO QUE ESTE GUION MIDE, y en particular donde se rompe el intento de prueba.
#  Intento (caso facil, y esta cerrado): si m = max S esta en g_com, su clase k tiene >= 2
#  elementos, todos los demas estan en T.  Si ALGUNO de ellos, m', esta en H, re-elegir el pick de
#  la clase k como m' da  deg' - deg = m - m' > 0: contradiccion.  Luego, si A1 falla, TODA la
#  clase de max S menos max S vive en L.  El caso que queda abierto da
#      deg' - deg = m + m' - 2*h_r,
#  cuyo signo NO es obvio.  N3 mide si ese caso tiene habitantes: si no los tiene, A1 esta PROBADA.
#
#  COLUMNAS
#    C0  fatal: re-monta la anatomia de c_eq_tau.py y comprueba sus dos numeros archivados
#        (2522 formas con [Phi]_top = 0, y H* en 1603 de ellas).  Si no cuadran, el resto no vale.
#    N1  DONDE VIVEN max S y min S: en K, en el empate {x1,x2,y1,y2}, o en g_com.  Tabla cruzada.
#    N2  A1 medida, y la cadena A1 => C = tau, con Implicacion (que cuenta el antecedente).
#    N3  EL HUECO DE LA PRUEBA: de las formas con max S en g_com (si las hay), en cuantas hay algun
#        elemento de su clase en H.  Esas estan probadas.  El resto es el hueco exacto.
#    N4  SEÑUELO QUE PUEDE FALLAR: para un elemento NO extremo de S, sigma_V(x) en S.  Si eso
#        tambien se cumpliera siempre, la propiedad no seria de los extremos y N2 no diria nada.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python extremes_of_S.py   (desde gates/)
# ============================================================================================

import json
import os
import sys
from collections import Counter

from _control import Contingencia, Implicacion

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "c_eq_tau.py")
_src = open(SRC, encoding="utf-8").read()
# OJO: c_eq_tau.py contiene ESA MISMA cadena como literal (hace el mismo truco con k_vs_m.py), y
# el literal aparece ANTES de anatomia().  Hay que cortar por el marcador CON su relleno de '='.
_head = _src.split("# ===================================================================== C0 =")[0]
assert "def anatomia(" in _head, "c_eq_tau.py cambio de forma"
_ns = {"__name__": "c_eq_tau_preamble", "__file__": SRC}
exec(compile(_head, SRC, "exec"), _ns)
anatomia = _ns["anatomia"]
shapes_of_width = _ns["shapes_of_width"]
CFG = _ns["CFG"]

OUT_JSON = "extremes_of_S_RESULT.json"


def poblacion():
    """todas las formas con [Phi]_top == 0, con su anatomia, en las 7 configuraciones de c_eq_tau."""
    for (t, r, Wmax) in CFG:
        N = t + 2 * r
        for W in range(N - 1, Wmax + 1):
            for b in shapes_of_width(W, N):
                a = anatomia(b, t, r)
                if a is None or not a.get('top0') or a.get('tau') is None:
                    continue
                yield (t, r, b, a)


print("=" * 116)
print("C0  ACEPTACION -- fatal.  Se re-montan los dos numeros ARCHIVADOS de c_eq_tau_OUT.txt.")
print("=" * 116)
print("")
POB = list(poblacion())
n_top0 = len(POB)
n_hstar = sum(1 for (t, r, b, a) in POB if a['S'][0] in a['K'] and a['S'][-1] in a['K'])
n_ctau = sum(1 for (t, r, b, a) in POB if a['C'] == a['tau'])
print("     formas con [Phi]_top == 0 : %d   (archivado: 2522)  %s"
      % (n_top0, "ok" if n_top0 == 2522 else "*** NO CUADRA ***"))
print("     H* (extremos en K)        : %d   (archivado: 1603)  %s"
      % (n_hstar, "ok" if n_hstar == 1603 else "*** NO CUADRA ***"))
print("     C == tau                  : %d   (archivado: 2522)  %s"
      % (n_ctau, "ok" if n_ctau == 2522 else "*** NO CUADRA ***"))
print("")
if n_top0 != 2522 or n_hstar != 1603 or n_ctau != 2522:
    print("  C0 FALLA -- la anatomia no reproduce lo archivado.  Veredicto SUSPENDIDO.")
    print("DONE (suspendido)")
    raise SystemExit(1)
print("  C0 PASA")
sys.stdout.flush()


def donde(v, a):
    if v in a['K']:
        return "K"
    if v in a['dif']:
        return "empate"
    if v in a['gcom']:
        return "g_com"
    return "???"


# ===================================================================== N1 ========================
print("")
print("=" * 116)
print("N1  DONDE VIVEN LOS DOS EXTREMOS DE S")
print("=" * 116)
print("")
tab = Counter()
for (t, r, b, a) in POB:
    tab[(donde(a['S'][-1], a), donde(a['S'][0], a))] += 1
print("     (donde max S, donde min S) : n")
for k in sorted(tab, key=lambda x: -tab[x]):
    print("        %-22s : %5d" % (str(k), tab[k]))
print("")
en_gcom = sum(v for k, v in tab.items() if "g_com" in k)
print("     con ALGUN extremo en g_com : %d" % en_gcom)
print("")

# ===================================================================== N2 ========================
print("=" * 116)
print("N2  A1  =>  C = tau,  y A1 misma")
print("=" * 116)
print("")
im = Implicacion("A1 (ningun extremo de S en g_com)", "C == tau")
cont = Contingencia("A1", "sigma_V(max S) en S  Y  sigma_V(min S) en S")
for (t, r, b, a) in POB:
    S = set(a['S'])
    V = a['tau']
    A1 = (a['S'][0] not in a['gcom']) and (a['S'][-1] not in a['gcom'])
    im.add(A1, a['C'] == a['tau'], b)
    cont.add(A1, (V - a['S'][-1]) in S and (V - a['S'][0]) in S, b)
im.informe(indent="     ")
print("")
print("     y la equivalencia algebraica A1 <=> sigma_V manda los dos extremos dentro de S")
print("     (deberia ser EXACTA: si no lo es, hay un error de codigo en alguna de las dos)")
cont.informe(indent="     ")
print("")

# ===================================================================== N3 ========================
print("=" * 116)
print("N3  EL HUECO DE LA PRUEBA: si max S estuviera en g_com, donde vive su clase")
print("=" * 116)
print("")
casos = [(t, r, b, a) for (t, r, b, a) in POB
         if a['S'][-1] in a['gcom'] or a['S'][0] in a['gcom']]
print("     formas con algun extremo en g_com : %d" % len(casos))
if not casos:
    print("")
    print("     -> A1 NO TIENE CONTRAEJEMPLO en esta poblacion.  El intento de prueba se cierra en")
    print("        su caso facil siempre que la clase del extremo tenga algun elemento en H; N3b")
    print("        mide cuantas veces hace falta el caso duro.")
else:
    ct = Counter()
    for (t, r, b, a) in casos[:20]:
        print("        t=%d r=%d beta=%-34s max S=%d min S=%d" % (t, r, str(b), a['S'][-1], a['S'][0]))
print("")
print("  N3b  el caso FACIL de la prueba: la clase de max S tiene algun elemento en H(T_A).")
print("       Si sale que SIEMPRE lo tiene, el argumento de re-eleccion cierra A1 sin el caso duro.")
facil = duro = 0
ej_duro = []
for (t, r, b, a) in POB:
    m = a['S'][-1]
    k = m % t
    clase = [v for v in a['S'] if v % t == k]
    H = a['TA'][-r:]                      # TA viene ORDENADO creciente de anatomia()
    otros = [v for v in clase if v != m]
    if any(v in H for v in otros):
        facil += 1
    else:
        duro += 1
        if len(ej_duro) < 6:
            ej_duro.append((t, r, b, m, clase, H))
print("       clase de max S con algun otro elemento en H : %d" % facil)
print("       sin ninguno (el caso DURO)                  : %d" % duro)
for e in ej_duro:
    print("          t=%d r=%d beta=%-32s max S=%d clase=%s H=%s" % e)
print("")

# ===================================================================== N4 ========================
print("=" * 116)
print("N4  SEÑUELO -- la propiedad tiene que ser de los EXTREMOS, no de todo S")
print("=" * 116)
print("")
print("     Para cada forma se toma un elemento NO extremo de S (el de en medio) y se pregunta lo")
print("     mismo: sigma_V(x) en S.  Si esto tambien saliera siempre, N2 no diria nada de nada.")
se = Contingencia("x = un elemento NO extremo de S", "sigma_V(x) en S")
n_med = n_ok = 0
for (t, r, b, a) in POB:
    S = a['S']
    if len(S) < 3:
        continue
    x = S[len(S) // 2]
    if x == S[0] or x == S[-1]:
        continue
    n_med += 1
    n_ok += ((a['tau'] - x) in set(S))
print("     elementos de en medio probados : %d ;  con sigma_V(x) en S : %d  (%.1f%%)"
      % (n_med, n_ok, 100.0 * n_ok / n_med if n_med else 0.0))
print("     -> %s" % ("el señuelo FALLA a menudo: la propiedad SI es de los extremos"
                      if n_ok < n_med else
                      "*** el señuelo tambien se cumple SIEMPRE: N2 no separa nada ***"))
print("")

# ===================================================================== N5 ========================
print("=" * 116)
print("N5  EL INTERCAMBIO -- la medida que A1 necesita para dejar de ser medida")
print("=" * 116)
print("")
print("     Sea m = max S y supongamos que una transversal T lo DESCARTA.  Toda su clase menos m")
print("     esta en T; sea m' el mayor de ellos.  Re-elegir el pick de esa clase como m' da")
print("         T' = T u {m} \\ {m'},   y como m > todo T, m entra en H.")
print("       - si m' esta en H :  Delta = m - m' > 0.  CERRADO.")
print("       - si m' esta en L :  Delta = m + m' - 2*h_r,  de signo NO obvio.  Esto es el hueco.")
print("     Y basta Delta >= 0, NO Delta > 0: con |G| <= 2 PROBADO, un TERCER maximizador ya es")
print("     contradiccion.  Aqui se mide Delta sobre TODAS las transversales que descartan max S.")
print("")
all_transversals = _ns["all_transversals"]
setup = _ns["setup"]
deg_of = _ns["deg_of"] if "deg_of" in _ns else (lambda T, r: sum(T[:r]) - sum(T[r:]))
cnt5 = Counter()
peor = []
form_mala = set()
disc_formula = 0
for (t, r, b, a) in POB:
    st = setup(b, t)
    if st is None:
        continue
    cl, E, Cd = st
    m = a['S'][-1]
    k = m % t
    clase = sorted((v for v in b if v % t == k), reverse=True)
    if len(clase) < 2 or clase[0] != m:
        continue
    mp = clase[1]
    for (sel, T, w, d) in all_transversals(b, cl, r, t):
        if m in T:
            continue
        Tp = tuple(sorted(list(T) + [m], reverse=True))
        Tp = tuple(x for x in Tp if x != mp) if mp in Tp else None
        if Tp is None or len(Tp) != 2 * r:
            continue
        delta = deg_of(Tp, r) - deg_of(T, r)
        H = T[:r]
        cnt5["m' en H" if mp in H else "m' en L"] += 1
        pred = (m - mp) if mp in H else (m + mp - 2 * T[r - 1])
        disc_formula += (pred != delta)
        cnt5["Delta > 0"] += (delta > 0)
        cnt5["Delta == 0"] += (delta == 0)
        cnt5["Delta < 0"] += (delta < 0)
        if delta <= 0:
            form_mala.add((t, r, b))
            if len(peor) < 6:
                peor.append((t, r, b, m, mp, delta, "H" if mp in H else "L"))
n5 = cnt5["Delta > 0"] + cnt5["Delta == 0"] + cnt5["Delta < 0"]
print("     transversales que descartan max S : %d" % n5)
print("        con m' en H : %d      con m' en L : %d" % (cnt5["m' en H"], cnt5["m' en L"]))
print("        Delta > 0   : %d" % cnt5["Delta > 0"])
print("        Delta == 0  : %d" % cnt5["Delta == 0"])
print("        Delta < 0   : %d   <-- si esto es > 0, el intercambio NO cierra A1" % cnt5["Delta < 0"])
print("        discrepancias formula-vs-calculo directo : %d %s"
      % (disc_formula, "" if not disc_formula else "*** LA FORMULA DE LA CABECERA ESTA MAL ***"))
for e in peor:
    print("        Delta<=0: t=%d r=%d beta=%-30s m=%d m'=%d Delta=%d (m' en %s)" % e)
print("        formas distintas con algun Delta <= 0 : %d de %d" % (len(form_mala), n_top0))
print("")

# ===================================================================== VEREDICTO =================
print("=" * 116)
print("VEREDICTO")
print("=" * 116)
print("")
print("     A1 (ningun extremo de S en g_com) : %d de %d" % (n_top0 - en_gcom, n_top0))
print("     A1 => C = tau : implicacion de DOS LINEAS, no medida -- el algebra esta en la cabecera.")
print("     Lo que queda por PROBAR es solo A1, y su caso facil cubre %d de %d." % (facil, n_top0))
print("")
print("     ALCANCE: la misma poblacion de c_eq_tau.py, o sea barrido con tope de reloj.")
json.dump({"n_top0": n_top0, "n_hstar": n_hstar, "A1": n_top0 - en_gcom,
           "facil": facil, "duro": duro, "senuelo_ok": n_ok, "senuelo_n": n_med,
           "tabla": {str(k): v for k, v in tab.items()}}, open(OUT_JSON, "w"), indent=1)
print("DONE")
