# -*- coding: utf-8 -*-
# LA CASTELLANA CONTRA LA INGLESA, EN SUSTANCIA Y NO EN CIFRAS.   25 de agosto de 2026.
#
# POR QUE.  El bloque N9 de check_numeros compara CIFRAS, y por eso no vio que la castellana
# estuviera tres revisiones por detras.  Se le anadio la comparacion de etiquetas y de entornos.
# Sigue faltando lo demas: ecuaciones en display, marcadores \ver, figuras, secciones.  Si la
# castellana pierde una ecuacion o una cita de guion, Carles revisa un documento distinto del que
# viaja --- y una copia de revision atrasada es peor que ninguna.
# [[a-numeric-guard-cannot-see-a-stale-translation]]
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_paridad_es.py > check_paridad_es_OUT.txt 2>&1

import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


def cuerpo(p):
    s = io.open(p, encoding="utf-8").read()
    s = s.split("\\begin{document}")[-1]
    return s.split("\\begin{thebibliography}")[0], s


EN, EN_full = cuerpo("characters_at_aP.tex")
ES, ES_full = cuerpo("characters_at_aP_es.tex")

print(SEP)
print("P1 -- INVENTARIO ESTRUCTURAL")
MEDIDAS = [
    ("ecuaciones en display", r"\\\["),
    ("ecuaciones numeradas", r"\\begin\{equation\}"),
    ("teoremas", r"\\begin\{theorem\}"),
    ("proposiciones", r"\\begin\{proposition\}"),
    ("corolarios", r"\\begin\{corollary\}"),
    ("observaciones", r"\\begin\{observation\}"),
    ("demostraciones", r"\\begin\{proof\}"),
    ("fichas (keybox)", r"\\begin\{keybox\}"),
    ("figuras", r"\\includegraphics"),
    ("secciones", r"\\section\{"),
    ("subsecciones", r"\\subsection\*\{"),
    ("etiquetas", r"\\label\{"),
]
print(f"{'':>26} {'EN':>5} {'ES':>5}")
for etq, pat in MEDIDAS:
    a, b = len(re.findall(pat, EN)), len(re.findall(pat, ES))
    print(f"{etq:>26} {a:>5} {b:>5}" + ("" if a == b else "   <<< DIFIERE"))
    ok(a == b, f"P1: mismo numero de {etq}", f"EN {a}, ES {b}")

print(SEP)
print("P2 -- LOS MARCADORES DE GUION \\ver{}: las dos versiones tienen que citar los MISMOS")
va = sorted(re.findall(r"\\ver\{([^}]*)\}", EN))
vb = sorted(re.findall(r"\\ver\{([^}]*)\}", ES))
sa, sb = set(va), set(vb)
ok(sa == sb, "P2: el mismo conjunto de guiones citados",
   f"solo EN: {sorted(sa - sb)}   solo ES: {sorted(sb - sa)}" if sa != sb else f"{len(sa)} guiones")
ok(len(va) == len(vb), "P2: y el mismo numero de citas a guion", f"EN {len(va)}, ES {len(vb)}")

print(SEP)
print("P2b -- TODO \\ver{} TIENE QUE NOMBRAR UN GUION DE VERDAD")
# ⚠ La nota llevaba un \ver{\,} en la pagina de atribucion --- "Every number marked [] comes from
#   a script" --- que renderizaba como un CORCHETE VACIO: la frase que explica los marcadores no
#   ensenaba ninguno.  Solo se vio MIRANDO LA PAGINA.  Ahora se exige que el argumento sea un .py.
import os
GATES = os.path.join("..", "gates")


def guion(marcador):
    """El nombre del .py que hay dentro de un \\ver{}.  Algunos llevan sufijo legitimo, como
    'prasad2016_vs_nuestro.py, \\S W', que nombra el guion Y el bloque: se admite el sufijo,
    pero tiene que haber un .py delante y el fichero tiene que existir."""
    t = marcador.replace("\\_", "_").strip()
    m = re.match(r"([A-Za-z0-9_\-.]+\.py)", t)
    return m.group(1) if m else None


for etq, lst in (("EN", va), ("ES", vb)):
    malos = [x for x in lst if guion(x) is None]
    ok(not malos, f"P2b: en {etq} todo \\ver{{}} nombra un guion .py", str(malos))
falta = sorted({guion(x) for x in va if guion(x)}
               - {f for f in os.listdir(GATES) if f.endswith(".py")})
ok(not falta, "P2b: y todos existen en gates/", str(falta))
ok(guion("\\,") is None and guion("aP\\_potencias.py") == "aP_potencias.py",
   "P2b: SENUELO, un marcador vacio SI se detecta y uno bueno no")

print(SEP)
print("P3 -- LAS ETIQUETAS Y LAS CITAS")
la, lb = set(re.findall(r"\\label\{([^}]*)\}", EN_full)), set(re.findall(r"\\label\{([^}]*)\}", ES_full))
ok(la == lb, "P3: las mismas etiquetas",
   f"solo EN: {sorted(la - lb)}   solo ES: {sorted(lb - la)}" if la != lb else f"{len(la)}")


def citas(t):
    return {k.strip() for g in re.findall(r"\\cite\{([^}]+)\}", t) for k in g.split(",")}


ca, cb = citas(EN_full), citas(ES_full)
ok(ca == cb, "P3: las mismas fuentes citadas",
   f"solo EN: {sorted(ca - cb)}   solo ES: {sorted(cb - ca)}" if ca != cb else f"{len(ca)}")

print(SEP)
print("P3b -- REFERENCIAS DEICTICAS: las que LaTeX no puede comprobar")
print("      ⚠ Al mover la seccion de G2 delante de los avisos, se rompieron TRES referencias que")
print("        no son \\ref y por tanto no dan 'undefined': el bloque de cierre se quedo DENTRO de")
print("        la seccion de G2, 'the next section' paso a apuntar a los avisos, y 'two paragraphs")
print("        below' dejo de contar bien.  Un documento con secciones moviles no puede apuntar")
print("        por posicion.")
DEICTICAS = ["the next section", "the last section", "two paragraphs below",
             "the previous section", "the section above", "la sección siguiente",
             "la última sección", "dos párrafos más abajo", "la sección anterior"]
for etq, t in (("EN", EN_full), ("ES", ES_full)):
    malas = [d for d in DEICTICAS if d in t]
    ok(not malas, f"P3b: en {etq} no hay referencias por posicion de seccion", str(malas))

# y el orden: la frase de cierre tiene que ser lo ULTIMO del cuerpo.
# 8-sep-2026: la frase que anclaba esto ("measure whatever you would like measured") era una oferta
# dirigida a un lector concreto y se cayo con la despersonalizacion.  El ancla pasa a ser la frase
# que de verdad cierra el cuerpo, y que no va dirigida a nadie.
for etq, t, frase, sec in (("EN", EN_full, "What I would most like to know, in order", r"\section{What this leaves"),
                           ("ES", ES_full, "Lo que más me gustaría saber, por orden", r"\section{Lo que esto deja")):
    if frase in t and sec in t:
        ok(t.index(frase) > t.index(sec),
           f"P3b: en {etq} la frase de cierre esta DENTRO de la seccion de cierre")
    else:
        ok(False, f"P3b: en {etq} no encuentro la frase de cierre o su seccion")

print(SEP)
print("P3c -- EL RECUENTO DE PARRAFOS DE PROSA")
print("      ⚠ La castellana llevaba un parrafo entero DESPUES de su frase de cierre, sin")
print("        equivalente en ingles, que resucitaba LA OBSTRUCCION FALSA --- 'la diferencia")
print("        segunda no sale del determinante como un producto' --- desmentida por el lema de")
print("        filas tres parrafos antes.  Comparar entornos y ecuaciones no lo veia: era prosa.")


def parrafos(t):
    cuerpo = t.split("\\begin{document}")[-1].split("\\begin{thebibliography}")[0]
    # fuera entornos, displays y tablas: quedan los parrafos de texto
    cuerpo = re.sub(r"\\begin\{(figure|tabular|center|itemize|equation)\}.*?"
                    r"\\end\{\1\}", " ", cuerpo, flags=re.S)
    cuerpo = re.sub(r"\\\[.*?\\\]", " ", cuerpo, flags=re.S)
    trozos = [p.strip() for p in re.split(r"\n\s*\n", cuerpo)]
    return [p for p in trozos if len(re.sub(r"\\[a-zA-Z]+|[{}$\\]", "", p).strip()) > 90]


pa, pb = parrafos(EN_full), parrafos(ES_full)
ok(abs(len(pa) - len(pb)) <= 1,
   "P3c: las dos versiones tienen (casi) los mismos parrafos de prosa",
   f"EN {len(pa)}, ES {len(pb)}")
# Senuelo: si se le anade un parrafo a una version, se tiene que notar.
# ⚠ hay que inyectarlo ANTES de \begin{thebibliography}, que es donde parrafos() corta: pegado al
#   final del fichero no lo veia nadie, y el control no podia fallar.
_inj = ES_full.replace("\\begin{thebibliography}",
                       "\n\nEste es un parrafo de mas, suficientemente largo como para que el "
                       "filtro de longitud lo cuente sin dudarlo, y asi el senuelo pueda "
                       "fallar de verdad.\n\n\\begin{thebibliography}", 1)
ok(len(parrafos(_inj)) == len(pb) + 1,
   "P3c: SENUELO, un parrafo de mas SI se detecta",
   f"{len(pb)} -> {len(parrafos(_inj))}")

print(SEP)
print("P4 -- SENUELO: el control tiene que ver una ecuacion que falte")
_falso = ES.replace("\\[", "", 1)
ok(len(re.findall(r"\\\[", _falso)) != len(re.findall(r"\\\[", EN)),
   "P4: quitar una ecuacion a la castellana SI se detecta")

print(SEP)
print("P5 -- ORTOGRAFIA: palabras que SIEMPRE llevan tilde, escritas sin ella")
# ⚠ 9-sep-2026.  Esta compuerta compara ESTRUCTURA, y por eso no vio lo siguiente: la
#   demostracion castellana del Teorema 5.1, escrita ese dia, salio ENTERA SIN ACENTOS.  Compilo
#   limpia, la paridad cuadro, las siete compuertas dieron verde, y el PDF imprimia "Ahi ponemos
#   ... un nucleo de Dirichlet".  Lo vio Carles leyendo.  Diecisiete palabras.
#   La lista lleva SOLO palabras que llevan tilde siempre: un primer intento incluyo "como",
#   "aun", "cardinal" y "criterio", que no la llevan, y dio 90 falsos positivos.  Una lista mal
#   hecha no mide nada.

SIEMPRE_CON_TILDE = """ahi alli aqui armonico asi caracter cientifico codigo composicion
computacion condicion conjugacion construccion contraccion cuestion demostracion descomposicion
dimension direccion division ecuacion eleccion especifico estadistico evaluacion expresion
factorizacion formula funcion generacion geometrico grafico identificacion informacion
interpretacion iteracion logica matematico maximo metodo minimo modulo multiplicacion notacion
numerico nucleo numero operacion parametro particion proposicion proyeccion publicacion rapido
razon reduccion region relacion representacion restriccion seccion segun simetrico simplectico
sintesis situacion solucion tambien teorico topologico traduccion transformacion ultimo union
unico unicos version""".split()


def prosa(s):
    """El texto que se IMPRIME: sin bibliografia, sin matematicas, sin comandos ni etiquetas.
    Las etiquetas importan: \\ref{sec:modulo} y \\ref{fig:union3d} traen 'modulo' y 'union'
    dentro, y son los diez falsos positivos que salieron al probar esto."""
    # ES ya viene de cuerpo(): sin preambulo y sin bibliografia.
    s = re.sub(r"(?s)\\\[.*?\\\]", " ", s)
    s = re.sub(r"\$[^$]*\$", " ", s)
    s = re.sub(r"\\[A-Za-z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?", " ", s)   # comando + su argumento
    return s


def sin_tilde(texto):
    fuera = []
    for p in sorted(set(SIEMPRE_CON_TILDE)):
        # \u26a0 SIN re.I esto no cazaba "Ahi" ni "Asi" a principio de frase, que son DOS DE LAS
        #   PALABRAS QUE FALLARON el 9-sep.  Un control que no ve el caso que lo motivo no sirve.
        for m in re.finditer(r"(?<![A-Za-z\u00c0-\u017f])" + p + r"(?![A-Za-z\u00c0-\u017f])",
                             texto, re.I):
            fuera.append((p, re.sub(r"\s+", " ", texto[max(0, m.start() - 45):m.end() + 20])))
    return fuera


_pr = prosa(ES)
_malas = sin_tilde(_pr)
ok(not _malas, f"P5: ninguna de las {len(set(SIEMPRE_CON_TILDE))} palabras vigiladas va sin tilde",
   "" if not _malas else f"{len(_malas)}: {[m[0] for m in _malas[:6]]}")
for _p, _c in _malas[:6]:
    print(f"           {_p}: ...{_c}...")

print(SEP)
print("P5b -- SENUELOS: los dos casos reales del 9-sep tienen que saltar")
_inj2 = _pr.replace("funci\u00f3n", "funcion", 1)
ok(len(sin_tilde(_inj2)) > len(_malas),
   "P5b: una tilde quitada en minusculas SI se detecta",
   f"{len(_malas)} -> {len(sin_tilde(_inj2))}")
# el caso REAL: "Ah\u00ed" al principio de frase, con mayuscula.  Sin re.I esto pasaba de largo.
_inj3 = _pr.replace("Ah\u00ed ", "Ahi ", 1)
ok(len(sin_tilde(_inj3)) > len(_malas),
   "P5b: y con MAYUSCULA inicial tambien --- 'Ahi', el caso que ocurrio",
   f"{len(_malas)} -> {len(sin_tilde(_inj3))}")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: la copia castellana tiene la misma estructura que la que viaja,")
print("           y ninguna palabra de la lista vigilada le falta la tilde.")
