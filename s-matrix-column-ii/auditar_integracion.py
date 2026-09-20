# -*- coding: utf-8 -*-
r"""auditar_integracion.py -- la nota como UN texto, no como seis sueltos.

Las otras tres auditorias miran piezas: los numeros (auditar_datos), las figuras contra su pie
(auditar_figuras), las citas y las etiquetas (auditar_refs).  Ninguna mira si el conjunto ENCAJA.
Eso es lo que falla cuando se mueve una seccion o se quita un apendice: el mapa de la seccion 1
promete algo que ya no esta donde decia, un teorema se queda sin nadie que lo invoque, un simbolo
se usa antes de definirse, o el resumen describe una version anterior.

Seis cosas, y cada una mirada en las DOS ediciones:

  (1) EL MAPA.  La seccion 1 anuncia que hace cada seccion.  Se exige que cada promesa nombre una
      seccion que existe y que esa seccion contenga lo prometido.
  (2) LAS COSTURAS.  Cada seccion cierra apuntando a la siguiente, que es la forma de la casa.
      Una seccion que acaba sin costura es una que se quedo suelta al reordenar.
  (3) ENUNCIADOS HUERFANOS.  Un teorema que nadie cita no esta integrado: o se usa, o sobra.
  (4) DIRECCION DE LAS REFERENCIAS.  Una referencia hacia delante DENTRO de una demostracion es
      una dependencia circular en el tiempo del lector.  Fuera de una demostracion es un anuncio,
      y eso si vale.
  (5) EL RESUMEN CONTRA EL CUERPO.  Cada cifra del resumen tiene que aparecer en el cuerpo.
  (6) NUMEROS QUE VIVEN EN VARIOS SITIOS.  El mismo dato dicho dos veces son dos claims.

    python auditar_integracion.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))

fallos = []
hechas = []


def afirma(ok, que):
    print("   %-4s %s" % ("ok" if ok else "MAL", que))
    hechas.append(que)
    if not ok:
        fallos.append(que)


def cuerpo_de(ruta):
    """el texto sin los comentarios de LaTeX: lo que el lector ve, no nuestras notas."""
    s = io.open(ruta, encoding="utf-8").read()
    return "\n".join(l for l in s.splitlines() if not l.lstrip().startswith("%"))


def secciones(s):
    """[(posicion, etiqueta, titulo, fin)] en orden de aparicion."""
    ms = list(re.finditer(r"\\section\{([^}]*)\}\\label\{(sec:\w+)\}", s))
    out = []
    for i, m in enumerate(ms):
        fin = ms[i + 1].start() if i + 1 < len(ms) else len(s)
        out.append((m.start(), m.group(2), m.group(1), fin))
    return out


def donde(secs, pos):
    n = "(cabecera)"
    for p, e, _, _ in secs:
        if p <= pos:
            n = e
    return n


# lo que el mapa de la seccion 1 promete de cada seccion, y como se comprueba que esta
PROMESAS = {
    "sec:defecto": ("lem:delta3", "la clasificacion de los defectos pequenos"),
    "sec:ley": ("thm:ley", "la ley 2/p - 1/p^2"),
    "sec:bring": ("eq:umbral", "el estrato maximal y el umbral"),
    "sec:centro": ("prop:semisistema", "los semisistemas"),
}

for fich, idioma in (("P3b.tex", "EN"), ("P3b_es.tex", "ES")):
    s = cuerpo_de(os.path.join(AQUI, fich))
    secs = secciones(s)
    orden = {e: i for i, (_, e, _, _) in enumerate(secs)}
    intro = s[secs[0][0]:secs[0][3]]

    print("")
    print("=" * 94)
    print("%s  --  %d secciones: %s" % (fich, len(secs), " -> ".join(e[4:] for _, e, _, _ in secs)))

    print("")
    print("(1) EL MAPA DE LA SECCION 1 PROMETE LO QUE CADA SECCION ENTREGA")
    for et, (etiqueta, que) in PROMESAS.items():
        prometida = ("\\ref{%s}" % et) in intro
        cuerpo_sec = next((s[p:f] for p, e, _, f in secs if e == et), "")
        entregada = ("\\label{%s}" % etiqueta) in cuerpo_sec
        afirma(prometida and entregada,
               "%s: el mapa la anuncia (%s) y ella contiene %s (%s)"
               % (et[4:], "si" if prometida else "NO", que, "si" if entregada else "NO"))

    print("")
    print("(2) NINGUNA SECCION ACABA EN UN DATO")
    # El criterio NO es "lleva una \ref hacia delante": la seccion de la ley cierra nombrando a
    # Bring y a Gauss en prosa, sin referencia, y eso es mejor prosa y no un defecto.  Lo que si
    # es binario, y lo que el reordenado rompio dos veces, es acabar dentro de un bloque \medido:
    # una seccion que termina en un numero se quedo sin entregar el testigo a la siguiente.
    for i, (p, et, _, f) in enumerate(secs[:-1]):
        cola = s[f - 400:f].rstrip()
        acaba_en_dato = cola.endswith("}") and cola.rfind("\\medido{") > cola.rfind("\n\n")
        afirma(not acaba_en_dato, "%s no termina dentro de un bloque medido" % et[4:])

    print("")
    print("(3) NINGUN ENUNCIADO SE QUEDA SIN QUE NADIE LO USE")
    for m in re.finditer(r"\\begin\{(?:theorem|proposition|lemma|corollary|observation)\}"
                         r"\\label\{([\w:]+)\}", s):
        et = m.group(1)
        citas = len(re.findall(r"\\ref\{%s\}" % re.escape(et), s))
        afirma(citas >= 1, "%-18s se invoca %d vez/veces" % (et, citas))

    print("")
    print("(4) NINGUNA DEMOSTRACION SE APOYA EN ALGO POSTERIOR")
    malas = 0
    for m in re.finditer(r"\\begin\{proof\}(?:\[[^\]]*\])?(.*?)\\end\{proof\}", s, re.S):
        desde = donde(secs, m.start())
        for r in re.finditer(r"\\ref\{([\w:]+)\}", m.group(1)):
            d = re.search(r"\\label\{%s\}" % re.escape(r.group(1)), s)
            if d and orden.get(donde(secs, d.start()), -1) > orden.get(desde, -1):
                print("        %s citado desde una demostracion de %s" % (r.group(1), desde))
                malas += 1
    afirma(malas == 0, "demostraciones que dependen de una seccion posterior: %d" % malas)

    print("")
    print("(5) CADA CIFRA DEL RESUMEN APARECE TAMBIEN EN EL CUERPO")
    res = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", s, re.S).group(1)
    resto = s[s.index("\\end{abstract}"):]
    faltan = []
    for n in sorted(set(re.findall(r"(?<![\w\\{])(\d{2,4})(?![\w}])", res))):
        if n not in resto:
            faltan.append(n)
    afirma(not faltan, "cifras del resumen ausentes del cuerpo: %s" % (faltan or "ninguna"))

    print("")
    print("(6) LOS DATOS QUE VIVEN EN VARIOS SITIOS DICEN LO MISMO")
    # Las tres primeras apariciones de defecto maximo: el dato que ya nos mordio una vez, y que
    # vive en tres sitios --- resumen, pie de figura y cuerpo.  Se cuentan los sitios donde los
    # tres numeros aparecen juntos, sin exigir que esten pegados: el display los separa.
    trios = len(re.findall(r"\b67\b(?:[^.]|\.\d){0,120}?\b163\b(?:[^.]|\.\d){0,120}?\b601\b", s))
    afirma(trios >= 3, "el trio 67/163/601 sale entero en %d sitios y coincide en todos" % trios)
    # los valores falsos, como NUMERO y no como subcadena: 101607 es un numero de articulo
    falsos = [v for v in ("233", "607")
              if re.search(r"(?<![\d.])%s(?![\d])" % v, re.sub(r"\\bibitem.*", "", s))]
    afirma(not falsos, "no queda ningun %s suelto, que eran los valores falsos"
           % (falsos or "233/607"))
    # las barras del censo viven SOLO en el pie de la figura 2, y es ahi donde auditar_figuras
    # las contrasta con la compuerta; aqui basta exigir que no se hayan duplicado a mano
    barras = {b: s.count(b) for b in ("132", "1348")}
    afirma(all(v == 1 for v in barras.values()),
           "las barras del censo viven en un solo sitio, el pie: %s" % barras)

print("")
print("=" * 94)
print("comprobaciones: %d ; fallos: %d" % (len(hechas), len(fallos)))
if not hechas:
    sys.exit("FATAL: cero comprobaciones. No poder medir no es aprobar")
print("VEREDICTO:", "la nota encaja consigo misma" if not fallos
      else "NO ENCAJA:\n   - " + "\n   - ".join(fallos))
sys.exit(1 if fallos else 0)
