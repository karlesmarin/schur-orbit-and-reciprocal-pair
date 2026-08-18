r"""Audita las figuras castellanas: ingles residual, y texto fuera del lienzo.

Dos cosas que el guion que las dibuja no puede ver por si solo:

  (1) una cadena que se colo en ingles porque nadie la mando a la capa de texto
      que se intercepto --- por ejemplo un rotulo puesto por una leyenda o por un
      objeto que se construyo antes del parcheo;
  (2) un rotulo que, al alargarse en castellano, se sale del lienzo de la figura.
      Eso no rompe nada ni avisa en el log: sale cortado en el PDF y ya.

Se mide sobre los PDF, no sobre el codigo.
"""
import io
import os
import sys

import fitz

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# palabras que en un rotulo solo pueden ser ingles.  Se dejan fuera las que se
# escriben IGUAL en castellano --- regular, ideal, diagonal, transversal, minimal,
# mod, log --- porque marcarlas fabrica un fallo donde no lo hay: la primera
# version de esta lista denuncio cinco figuras correctas por eso.
INGLES = {
    "the", "and", "with", "that", "this", "from", "when", "does", "what",
    "where", "which", "point", "points", "forms", "form", "weight", "weights",
    "survive", "survives", "surviving", "vanish", "vanishing", "odd",
    "cores", "empty", "frozen", "count", "signed", "proved",
    "cited", "single", "fusion", "group", "condition",
    "step", "matrix", "explicit", "left", "still", "lattice",
    "unhit", "quotient", "numerator", "landing",
    "inset", "magnified", "magnitude", "value",
    "absolute", "partial", "terms", "added", "decreasing", "dominant",
    "vertex", "exponents", "vertices", "subtracted", "denominator", "read",
    "coordinates", "locus", "equal", "unequal", "tested",
    "cases", "orbits", "saving", "colour", "evaluation", "walls",
    "families", "middle", "pair", "absent", "collapse", "available", "dies",
    "generator", "torsion", "interlacing", "parity",
    "killed", "cone", "sublattice", "edge", "planes", "alternate", "survivors",
    "adding", "always", "leaves", "meets", "support", "once",
    "twice", "cancel", "finite", "difference", "filled", "ringed", "blue",
    "arrows", "shifts", "larger", "than", "coefficient",
    "axes", "shapes", "wanders", "lands", "every", "path",
}

def desborde(nombre):
    """Cuanto se sale el texto del lienzo, en puntos, y cuantas palabras."""
    doc = fitz.open(nombre)
    pg = doc[0]
    caja = pg.rect
    pal = pg.get_text("words")
    peor, n = 0.0, 0
    for w in pal:
        d = max(caja.x0 - w[0], w[2] - caja.x1, caja.y0 - w[1], w[3] - caja.y1)
        if d > 0.5:
            n += 1
            peor = max(peor, d)
    doc.close()
    return len(pal), n, peor


fallos = 0
print("Figuras castellanas: ingles residual, y desborde CONTRA la inglesa")
print("=" * 78)
print("  (un desborde que ya existe en la inglesa es diseno de la figura, no")
print("   dano de la traduccion; lo que importa es que el castellano no lo empeore)")
print("")

figs = sorted(f for f in os.listdir(".") if f.startswith("fig_") and f.endswith("_es.pdf"))
for f in figs:
    doc = fitz.open(f)
    pg = doc[0]
    palabras = pg.get_text("words")
    sospechosas = sorted({w[4] for w in palabras if w[4].lower().strip(".,;:()") in INGLES})
    doc.close()

    _, n_es, peor_es = desborde(f)
    en = f.replace("_es.pdf", ".pdf")
    if os.path.exists(en):
        _, n_en, peor_en = desborde(en)
    else:
        n_en, peor_en = -1, -1.0

    empeora = peor_es > peor_en + 1.0
    mal = bool(sospechosas) or empeora
    fallos += 1 if mal else 0
    print("  [%s] %-26s palabras=%3d  ingles=%d  desborde EN=%.1fpt (%d) ES=%.1fpt (%d)"
          % ("!! " if mal else "ok ", f, len(palabras), len(sospechosas),
             peor_en, n_en, peor_es, n_es))
    if sospechosas:
        print("         %s" % sospechosas[:14])

# y el articulo compuesto: ninguna pagina anterior a la bibliografia debe llevar
# esas palabras.  La bibliografia SI las lleva, y con razon: son titulos citados.
doc = fitz.open("orbit_pair_ii_es.pdf")
bib = None
for i in range(doc.page_count):
    if "Referencias" in doc[i].get_text() or "References" in doc[i].get_text():
        bib = i
        break
if bib is None:
    bib = doc.page_count
malas = {}
for i in range(bib):
    for w in doc[i].get_text().split():
        t = w.lower().strip(".,;:()$")
        # ojo: `t not in ("value")` seria pertenencia a la CADENA, no al conjunto
        if t in INGLES and t not in {"value", "element"}:
            malas.setdefault(t, []).append(i + 1)
print("")
print("  articulo compuesto: %d paginas antes de la bibliografia (pag. %s)"
      % (bib, bib + 1))
if malas:
    print("  [!! ] palabras inglesas en el cuerpo:")
    for t, ps in sorted(malas.items())[:20]:
        print("         %-14s paginas %s" % (t, sorted(set(ps))[:8]))
    fallos += 1
else:
    print("  [ok ] ninguna palabra inglesa de la lista en el cuerpo")
doc.close()

print("")
print("  figuras auditadas: %d   fallos: %d" % (len(figs), fallos))
