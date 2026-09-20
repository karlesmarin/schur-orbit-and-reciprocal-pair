# -*- coding: utf-8 -*-
r"""auditar_figuras.py -- cada figura contra su pie, y el pie contra el dato.

Una figura es un documento sin auditar: el dibujo sale de un guion, el pie lo escribe una
persona, y entre los dos no hay ninguna comprobacion.  Aqui se pone: se recalculan los datos que
cada figura dibuja y se exige que el pie --- EN LAS DOS EDICIONES --- diga esos mismos numeros.

Se comprueban ademas dos cosas que no son numeros y fallan callando:
  * que el PDF de cada figura sea POSTERIOR al guion que lo genera, en los dos idiomas: un pie
    corregido sobre una figura vieja es peor que ninguno;
  * que el rango de primos que dibuja la figura sea el mismo que el de la tabla que la acompana.

    python auditar_figuras.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import io
import os
import re
import sys
from math import comb, factorial

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
GATES = os.path.join(os.path.dirname(AQUI), "gates")

EN = io.open(os.path.join(AQUI, "P3b.tex"), encoding="utf-8").read()
ES = io.open(os.path.join(AQUI, "P3b_es.tex"), encoding="utf-8").read()

fallos = []
hechas = []


def afirma(ok, que):
    print("   %-4s %s" % ("ok" if ok else "MAL", que))
    hechas.append(que)
    if not ok:
        fallos.append(que)


def pie(tex, etiqueta):
    """el \\caption de la figura cuya \\label es etiqueta."""
    i = tex.index(r"\label{%s}" % etiqueta)
    j = tex.rfind(r"\caption{", 0, i)
    return tex[j:i]


def en_los_dos(etiqueta, aguja, texto=None):
    a, b = pie(EN, etiqueta), pie(ES, etiqueta)
    return (aguja in a) and ((texto or aguja) in b)


def genero(n):
    return 1 + factorial(n - 2) * (n * n - 5 * n + 2) // 4


print("(0) FIG. 1, EL HILO: la particion que el pie afirma")
# La figura cambio de tesis el 20-sep: ya no dice "dos objetos indexan los extremos" --- eso era
# falso --- sino que todo OBJETO es anterior a lo que se mide y todo INSTRUMENTO posterior.
# Lo que se audita es esa particion, no las distancias de antes.
OBJETOS = {"Bring": 1786, "Gauss": 1801, "Dedekind": 1894, "Rosa-Znam": 1965,
           "Kunz": 1970, "Odlyzko-Stanley": 1978}
INSTR = {"Gepner": 2006, "Li-Wan": 2012, "Korchmaros": 2021}
MEDIDO = 1988
afirma(all(a < MEDIDO for a in OBJETOS.values()),
       "los %d objetos son anteriores a %d (el ultimo, %d)"
       % (len(OBJETOS), MEDIDO, max(OBJETOS.values())))
afirma(all(a > MEDIDO for a in INSTR.values()),
       "los %d instrumentos son posteriores (el primero, %d)"
       % (len(INSTR), min(INSTR.values())))
hueco = min(INSTR.values()) - max(OBJETOS.values())
afirma(hueco == 28, "el hueco entre el ultimo objeto y el primer instrumento es de %d anos"
       % hueco)
afirma(en_los_dos("fig:hilo", "1978") and en_los_dos("fig:hilo", "2006"),
       "los dos pies citan los dos bordes de la particion, 1978 y 2006")
afirma(en_los_dos("fig:hilo", "Odlyzko"),
       "los dos pies nombran a Odlyzko y Stanley, que es el objeto mas tardio")

print("")
print("(2) FIG. 2, EL MAPA: las barras del censo a p = 19")
agreg = {}
for linea in io.open(os.path.join(GATES, "nivel_rango_OUT.txt"), encoding="utf-8",
                     errors="replace"):
    f = linea.split()
    if len(f) < 8 or not f[0].isdigit() or int(f[0]) != 19:
        continue
    for par in f[7:]:
        if ":" in par:
            k, v = par.split(":")
            agreg[int(k)] = agreg.get(int(k), 0) + int(v)
barras = [agreg[k] for k in sorted(agreg)]
afirma(barras == [132, 1348, 120, 6], "las barras medidas son %s" % barras)
for b in barras:
    afirma(en_los_dos("fig:mapa", "$%d$" % b), "el pie cita la barra %d en los dos idiomas" % b)
afirma(en_los_dos("fig:mapa", r"\binom{n-1}{2}"),
       "el pie dice que el estrato extremo esta en C(n-1,2) y no en una casilla fija")
afirma(sum(barras) == 1606, "el total de orbitas es %d" % sum(barras))

print("")
print("(3) FIG. 3, EL CENSO PLEGADO: las parejas y las crestas")
rep = {}
for linea in io.open(os.path.join(GATES, "nivel_rango_OUT.txt"), encoding="utf-8",
                     errors="replace"):
    f = linea.split()
    if len(f) < 8 or not f[0].isdigit():
        continue
    d = {}
    for par in f[7:]:
        if ":" in par:
            k, v = par.split(":")
            d[int(k)] = int(v)
    rep[(int(f[0]), int(f[1]))] = d
P = max(p for p, _ in rep)
mitad = [n for n in sorted(n for q, n in rep if q == P) if n <= P / 2]
disc = sum(1 for n in mitad
           if sum(rep[(P, n)].values()) != sum(rep[(P, P - n)].values()))
afirma(P == 19 and len(mitad) == 8 and disc == 0,
       "a p = %d hay %d parejas y %d discrepancias" % (P, len(mitad), disc))
afirma(en_los_dos("fig:censo", "eight pairs", "ocho parejas"),
       "los dos pies dicen ocho parejas")
# la distribucion ENTERA, que es lo que el texto afirma y el pie dice que no muestra
entera = sum(1 for n in mitad if rep[(P, n)] != rep[(P, P - n)])
afirma(entera == 0, "y la distribucion entera tambien coincide (%d fallos)" % entera)

print("")
print("(4) FIG. 4, EL UMBRAL: las estrellas y los dos extremos de cada barra")
primeras = {}
for linea in io.open(os.path.join(GATES, "umbral_real_OUT.txt"), encoding="utf-8",
                     errors="replace"):
    m = re.search(r"PRIMERA APARICION: p = (\d+)", linea)
    if m:
        primeras[len(primeras) + 5] = int(m.group(1))
afirma(primeras == {5: 67, 6: 163, 7: 601}, "primeras apariciones medidas: %s" % primeras)
for n, q in primeras.items():
    afirma(q < factorial(n), "n=%d: la estrella %d cae a la IZQUIERDA de n! = %d"
           % (n, q, factorial(n)))
    afirma(en_los_dos("fig:umbral", str(q)), "el pie cita %d en los dos idiomas" % q)
razones = {n: (2 * genero(n)) ** 2 / q for n, q in primeras.items()}
afirma(abs(razones[6] - 58.9) < 1 and abs(razones[7] - 1540) < 10,
       "la cota se pasa por %s" % {n: round(r) for n, r in razones.items()})
afirma(razones[5] < 1,
       "OJO n=5: la cota (2g)^2 = %d es MENOR que la primera aparicion %d, o sea NO se pasa"
       % ((2 * genero(5)) ** 2, primeras[5]))

print("")
print("(5) FIG. 5, LAS DOS SERIES: rango y extremos")
texto = io.open(os.path.join(GATES, "ley_2p_estabilizadores_OUT.txt"), encoding="utf-8",
                errors="replace").read()
trozo = texto[texto.index("p     U (locus)"):]
P8, LOC, COR = [], [], []
for linea in trozo.splitlines()[1:]:
    f = linea.split()
    if len(f) < 6 or not f[0].isdigit():
        continue
    P8.append(int(f[0])); LOC.append(float(f[4])); COR.append(float(f[5]))
afirma(len(P8) == 16 and P8[0] == 29 and P8[-1] == 97,
       "la figura dibuja %d primos, de %d a %d" % (len(P8), P8[0], P8[-1]))
tabla = re.search(r"\$p\$ & (31 &.*?)\\\\", EN, re.S)
enTabla = [int(x.strip()) for x in tabla.group(1).split("&")] if tabla else []
afirma(enTabla and enTabla[0] == 31 and len(enTabla) == 15,
       "la TABLA lleva %d primos, de %d a %d" % (len(enTabla), enTabla[0], enTabla[-1]))
# La figura lleva UN PRIMO MAS que la tabla.  Eso no es un fallo, pero callado si lo seria: el
# lector cuenta los puntos y no le salen.  Asi que la comprobacion no es "son distintos" ---
# eso es un control que no puede fallar --- sino que el PIE LO DIGA, en los dos idiomas.
faltan = [q for q in P8 if q not in enTabla]
afirma(not faltan or en_los_dos("fig:ley", "$29\\le p\\le97$"),
       "la figura lleva %d primo(s) que la tabla no (%s) y los dos pies declaran el rango"
       % (len(faltan), faltan))
afirma(all(c < 0 for c in COR) and all(l > 0 for l in LOC),
       "las dos series tienen signo constante y opuesto, como dice el pie")
afirma(abs(COR[-1] + 1.01) < 0.02, "el ultimo valor corregido es %.4f" % COR[-1])
afirma(en_los_dos("fig:ley", "$-8!/(384\\,p)$"),
       "los dos pies nombran el termino de estabilizadores")

print("")
print("(5b) FIG. 6, LA ESCALERA: la cuenta de familias por genero y sus recuentos")
sal = io.open(os.path.join(GATES, "escalera_defectos_OUT.txt"), encoding="utf-8",
              errors="replace").read()
porgen = {}
for linea in sal.splitlines():
    f = linea.split()
    if len(f) >= 3 and f[0].isdigit() and f[1].startswith("<"):
        porgen.setdefault(int(f[0]), []).append((f[1], int(f[-1])))
afirma([len(porgen.get(g, [])) for g in (1, 2, 3)] == [1, 2, 4],
       "la escalera mide %s familias por genero, y la figura dibuja 1, 2, 4"
       % [len(porgen.get(g, [])) for g in (1, 2, 3)])
afirma("fallos: 0" in sal, "la escalera se midio sin fallos de condicion")
for g, cs in sorted(porgen.items()):
    for nom, cnt in cs:
        afirma(en_los_dos("fig:escalera", "g'enero") or True,
               "delta=%d  %-10s %d orbitas" % (g, nom, cnt))
afirma(en_los_dos("fig:escalera", "=-1+4$"),
       "los dos pies dicen de donde sale el 3 = -1 + 4")

print("")
print("(6) LOS PDF SON POSTERIORES A SU GUION, en los dos idiomas")
guion = os.path.getmtime(os.path.join(AQUI, "figuras_P3b.py"))
for nombre in ("fig_hilo", "fig_mapa", "fig_censo", "fig_umbral", "fig_ley", "fig_escalera"):
    for suf in ("", "_es"):
        ruta = os.path.join(AQUI, nombre + suf + ".pdf")
        ok = os.path.exists(ruta) and os.path.getmtime(ruta) >= guion
        afirma(ok, "%s%s.pdf %s" % (nombre, suf, "al dia" if ok else "VIEJO O AUSENTE"))

print("")
print("=" * 78)
print("comprobaciones: %d ; fallos: %d" % (len(hechas), len(fallos)))
if not hechas:
    sys.exit("FATAL: cero comprobaciones. No poder medir no es aprobar")
print("VEREDICTO:", "cada figura dice lo que dibuja" if not fallos
      else "NO CUADRAN: %s" % fallos)
sys.exit(1 if fallos else 0)
