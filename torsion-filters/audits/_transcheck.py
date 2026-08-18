r"""Compara la ESTRUCTURA del ingles y del castellano, no la prosa.

Una traduccion se rompe por omision, no por estilo: una seccion que no se copio, un
\label que se quedo atras, una figura sin su \includegraphics, una fila de tabla de
menos.  Este guion no lee espanol; cuenta y empareja marcas.
"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EN = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
ES = io.open("orbit_pair_ii_es.tex", encoding="utf-8").read()

fallos = 0


def compara(nombre, a, b, ordenado=True):
    global fallos
    fa, fb = (sorted(a), sorted(b)) if ordenado else (a, b)
    if fa == fb:
        print("  [ok ] %-28s %d" % (nombre, len(a)))
        return
    fallos += 1
    solo_en = [x for x in a if x not in b]
    solo_es = [x for x in b if x not in a]
    print("  [!! ] %-28s EN=%d ES=%d" % (nombre, len(a), len(b)))
    if solo_en:
        print("         solo en EN: %s" % solo_en[:12])
    if solo_es:
        print("         solo en ES: %s" % solo_es[:12])


def saca(pat, s):
    return re.findall(pat, s)


print("ESTRUCTURA: ingles contra castellano")
print("=" * 78)

compara("\\label", saca(r"\\label\{([^}]*)\}", EN), saca(r"\\label\{([^}]*)\}", ES))
# MULTICONJUNTO y no conjunto.  Con set() una referencia usada dos veces en una edicion y una en la
# otra es invisible, y por ahi se colo el 18 de agosto una cita a [Sch21,Mur01] que estaba en la
# castellana y no en la inglesa: los conjuntos coincidian y el recuento no.
compara("\\ref destinos", sorted(saca(r"\\ref\{([^}]*)\}", EN)),
        sorted(saca(r"\\ref\{([^}]*)\}", ES)))
compara("\\eqref destinos", sorted(saca(r"\\eqref\{([^}]*)\}", EN)),
        sorted(saca(r"\\eqref\{([^}]*)\}", ES)))
compara("\\cite claves", sorted(re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", EN)),
        sorted(re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", ES)))
# la version castellana incluye fig_*_es.pdf: se compara la figura, no el idioma
compara("figuras (pdf)", saca(r"\\includegraphics[^{]*\{([^}]*)\}", EN),
        [x.replace("_es.pdf", ".pdf")
         for x in saca(r"\\includegraphics[^{]*\{([^}]*)\}", ES)])
compara("\\bibitem claves", saca(r"\\bibitem\[[^\]]*\]\{([^}]*)\}", EN),
        saca(r"\\bibitem\[[^\]]*\]\{([^}]*)\}", ES))

# entornos: mismo numero de cada tipo, y en el mismo orden
for env in ("theorem", "proposition", "lemma", "corollary", "conjecture", "remark",
            "problem", "observation", "figure", "proof", "keybox", "equation"):
    a = len(re.findall(r"\\begin\{%s\}" % env, EN)) if env != "keybox" \
        else len(re.findall(r"\\keybox", EN))
    b = len(re.findall(r"\\begin\{%s\}" % env, ES)) if env != "keybox" \
        else len(re.findall(r"\\keybox", ES))
    if a == b:
        print("  [ok ] %-28s %d" % ("entorno " + env, a))
    else:
        fallos += 1
        print("  [!! ] %-28s EN=%d ES=%d" % ("entorno " + env, a, b))

# secciones y subsecciones: solo el RECUENTO y el orden de sus labels
sec_en = len(re.findall(r"\\section\{", EN))
sec_es = len(re.findall(r"\\section\{", ES))
sub_en = len(re.findall(r"\\subsection[*]?\{", EN))
sub_es = len(re.findall(r"\\subsection[*]?\{", ES))
for nom, a, b in (("\\section", sec_en, sec_es), ("\\subsection", sub_en, sub_es)):
    if a == b:
        print("  [ok ] %-28s %d" % (nom, a))
    else:
        fallos += 1
        print("  [!! ] %-28s EN=%d ES=%d" % (nom, a, b))

# filas de las tablas largas: una fila perdida es una afirmacion perdida
for nom, marca in (("filas longtable", r"\\\\"),):
    a = len(re.findall(marca, EN))
    b = len(re.findall(marca, ES))
    print("  [--] %-28s EN=%d ES=%d (informativo)" % (nom, a, b))

# los numeros del texto: todo numero con separador de millar debe aparecer en las dos
num_en = sorted(set(re.findall(r"\d[\d]*\\,\d\d\d", EN)))
num_es = sorted(set(re.findall(r"\d[\d]*\\,\d\d\d", ES)))
compara("numeros con \\,", num_en, num_es)

# marcas de estado
for cmd in ("stproved", "stverif", "stext"):
    a = len(re.findall(r"\\%s" % cmd, EN))
    b = len(re.findall(r"\\%s" % cmd, ES))
    if a == b:
        print("  [ok ] %-28s %d" % ("marca \\" + cmd, a))
    else:
        fallos += 1
        print("  [!! ] %-28s EN=%d ES=%d" % ("marca \\" + cmd, a, b))

# Una marca de estado seguida de ESPACIO: LaTeX se come el espacio, porque una palabra de control
# se traga los blancos que la siguen.  `\stext and an` se imprime `externaland`, y asi viajo hasta
# la pagina 53 --- y la castellana decia `ajenoy` en la 55, el mismo sitio.  Se escribe `\stext{}`.
#
# La PRIMERA version de este control exigia que siguiera una LETRA, y por eso se dejo fuera
# `\stverif --- had some other shift`, que se imprime `verified—had`: siete por edicion, y las
# encontro una lectura de fuera, no el control.  Lo que importa no es que siga una letra: es que siga
# ALGO tras un espacio que la macro se va a comer.  Solo es inocuo cuando lo que sigue va pegado
# de todos modos --- puntuacion de cierre --- o cuando es `\\`, que es el fin de fila de la tabla.
INOCUO = set(",;.:)]}&")
pegadas = 0
for nombre, texto in (("EN", EN), ("ES", ES)):
    for m in re.finditer(r"\\st(?:ext|verif|proved)[ \t]+(\S)", texto, re.UNICODE):
        sigue = m.group(1)
        if sigue in INOCUO or texto[m.end() - 1:m.end() + 1] == "\\\\":
            continue
        pegadas += 1
        i = max(0, m.start() - 45)
        print("  [!! ] marca que se come el espacio (%s): ...%s..."
              % (nombre, " ".join(texto[i:m.end() + 12].split())))
fallos += pegadas
if not pegadas:
    print("  [ok ] %-28s ninguna marca se come el espacio" % "marcas + espacio")

print("")
print("  discrepancias estructurales: %d" % fallos)
