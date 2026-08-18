r"""Dibuja las doce figuras en castellano SIN duplicar los doce guiones.

El problema: los rotulos internos de las figuras estan en ingles, y la version
castellana del articulo los necesita traducidos.  Copiar cada guion seria crear
doce ficheros que se desincronizan del original en la primera correccion.

En su lugar se intercepta la capa de texto de matplotlib --- `text`, `annotate`,
`set_xlabel`, `set_ylabel`, `set_title`, `legend`, `set_xticklabels`, ... --- y se
pasa cada cadena por un diccionario antes de dibujarla; y se reescribe el nombre
del fichero de salida a `*_es.pdf`.  Los guiones se ejecutan tal cual, de modo
que TODOS SUS CONTROLES SIGUEN CORRIENDO: si un guion se niega a dibujar porque
sus datos no cuadran, tambien se niega aqui.

Dos modos:

    python figs_es.py --recoge     escribe figs_es_CADENAS.txt con lo que se dibuja
    python figs_es.py              dibuja, traduciendo por FIGLANG

Lo que no se traduce, a proposito: la matematica (`$...$` sin palabras), los
nombres propios, y las cadenas que no aparecen en el diccionario --- esas se
listan al final de la corrida para que no pasen desapercibidas.
"""
import io
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RECOGE = "--recoge" in sys.argv
VISTAS = []          # cadenas que llegaron a la capa de texto, en orden
SIN_TRADUCIR = set()  # las que no estaban en el diccionario

from figlang import FIGLANG, NOTACION  # noqa: E402


# Texto puesto DENTRO del modo matematico.  `\mathrm{alphabet}` es una palabra inglesa aunque
# viva entre $...$, y los comandos de LaTeX --- \prod, \lambda, \sum --- no lo son.
_TEXTO_EN_MATE = re.compile(
    r"\\(?:math(?:rm|it|sf|tt|bf|cal)|text(?:rm|it|bf|sf|tt|normal)?|operatorname\*?)"
    r"\s*\{([^{}]*)\}"
)


def _solo_matematica(s):
    """True si la cadena es matematica pura: nada que traducir dentro.

    Lo de entre `$...$` es matematica, PERO no todo lo que hay ahi lo es: `\\mathrm{...}`,
    `\\text{...}` y `\\operatorname{...}` son TEXTO colocado en modo matematico.  La primera
    version de esta regla borraba el modo matematico entero y miraba solo lo que quedaba fuera,
    de modo que

        $\\prod(\\mathrm{alphabet})=(-1)^{t-1}=-1$

    no dejaba nada fuera, se declaraba matematica pura, y se dibujo en INGLES en los tres paneles
    de `fig_alphabet_es.pdf` del Paper I --- sin traducir y **sin aparecer en la lista de sin
    traducir**, que es lo que la hizo invisible durante semanas.  Ahora el contenido de esos
    comandos se extrae y se juzga como texto.  Portado aqui el 18 de agosto de 2026.
    """
    fuera = re.sub(r"\$[^$]*\$", "", s)
    dentro = " ".join(w for w in _TEXTO_EN_MATE.findall(s) if w.strip() not in NOTACION)
    return not re.search(r"[A-Za-z]{3,}", fuera + " " + dentro)


YA_TRADUCIDO = set(FIGLANG.values())


def tr(s):
    if not isinstance(s, str) or not s.strip():
        return s
    VISTAS.append(s)
    if s in FIGLANG:
        return FIGLANG[s]
    if s in YA_TRADUCIDO:
        # algunos guiones releen un rotulo y lo vuelven a poner: ya esta en castellano
        return s
    if not _solo_matematica(s):
        SIN_TRADUCIR.add(s)
    return s


# --- parcheo de la capa de texto -------------------------------------------
_orig = {}


def envuelve(cls, nombre, pos):
    """Envuelve cls.nombre traduciendo el argumento posicional `pos`."""
    f = getattr(cls, nombre, None)
    if f is None:
        return
    _orig[(cls, nombre)] = f

    def envuelto(self, *a, **k):
        a = list(a)
        if len(a) > pos:
            if isinstance(a[pos], str):
                a[pos] = tr(a[pos])
            elif isinstance(a[pos], (list, tuple)):
                a[pos] = [tr(x) for x in a[pos]]
        for clave in ("label", "title", "s", "text", "xlabel", "ylabel"):
            if isinstance(k.get(clave), str):
                k[clave] = tr(k[clave])
        return f(self, *a, **k)

    setattr(cls, nombre, envuelto)


for nombre, pos in (("text", 2), ("set_xlabel", 0), ("set_ylabel", 0),
                    ("set_zlabel", 0), ("set_title", 0), ("annotate", 0),
                    ("set_xticklabels", 0), ("set_yticklabels", 0),
                    ("set_zticklabels", 0), ("bar", 99), ("plot", 99),
                    ("scatter", 99), ("axhline", 99), ("axvline", 99),
                    ("fill_between", 99), ("stairs", 99), ("step", 99)):
    envuelve(Axes, nombre, pos)
try:
    from mpl_toolkits.mplot3d import Axes3D

    for nombre, pos in (("text", 3), ("set_xlabel", 0), ("set_ylabel", 0),
                        ("set_zlabel", 0), ("set_title", 0), ("text2D", 2)):
        envuelve(Axes3D, nombre, pos)
except Exception:
    pass

envuelve(Figure, "text", 2)
envuelve(Figure, "suptitle", 0)

# El gancho que de verdad cierra la puerta.  Envolver las funciones de alto nivel
# deja fuera todo rotulo que llegue por otra via --- `set_xticks(pos, labels)`, una
# leyenda construida a mano, un `Text` creado directamente --- y la primera pasada
# dejo asi seis figuras con ingles dentro.  TODO texto de matplotlib termina en
# `Text.set_text`, asi que se traduce ahi y no hay puerta de atras.
from matplotlib.text import Text  # noqa: E402

_set_text = Text.set_text


def _set_text_env(self, s):
    return _set_text(self, tr(s) if isinstance(s, str) else s)


Text.set_text = _set_text_env

_text_init = Text.__init__


def _text_init_env(self, x=0, y=0, text="", *a, **k):
    return _text_init(self, x, y, tr(text) if isinstance(text, str) else text, *a, **k)


Text.__init__ = _text_init_env

_legend_axes = Axes.legend


def _legend_env(self, *a, **k):
    if a and isinstance(a[0], (list, tuple)) and a and all(isinstance(x, str) for x in a[0]):
        a = ([tr(x) for x in a[0]],) + tuple(a[1:])
    if len(a) > 1 and isinstance(a[1], (list, tuple)):
        a = (a[0], [tr(x) for x in a[1]]) + tuple(a[2:])
    if isinstance(k.get("title"), str):
        k["title"] = tr(k["title"])
    return _legend_axes(self, *a, **k)


Axes.legend = _legend_env

# --- redireccion de la salida ----------------------------------------------
SALIDAS = []
_savefig_fig = Figure.savefig


def _savefig_env(self, fname, *a, **k):
    if isinstance(fname, str) and fname.endswith(".pdf") and not fname.endswith("_es.pdf"):
        fname = fname[:-4] + "_es.pdf"
        SALIDAS.append(fname)
    return _savefig_fig(self, fname, *a, **k)


Figure.savefig = _savefig_env
plt.savefig = lambda *a, **k: _savefig_env(plt.gcf(), *a, **k)

# --- ejecucion de los doce guiones -----------------------------------------
GUIONES = sorted(f for f in os.listdir(".") if f.startswith("fig_") and f.endswith(".py"))

print("Figuras en castellano: %d guiones" % len(GUIONES))
print("=" * 78)
fallos = 0
for g in GUIONES:
    src = io.open(g, encoding="utf-8").read()
    entorno = {"__name__": "__main__", "__file__": os.path.abspath(g)}
    try:
        exec(compile(src, g, "exec"), entorno)
        print("  [ok ] %s" % g)
    except SystemExit as ex:
        # Un SystemExit es como los guiones de figura AVISAN de que un control ha disparado --- el
        # de pixeles, por ejemplo, que en castellano salta donde en ingles no, porque el castellano
        # es mas largo.  Contarlo como "salida limpia" dejaba la figura sin generar y sin ruido.
        fallos += 1
        print("  [!! ] %s -> control disparado: %s" % (g, str(ex).strip().splitlines()[0]
                                                       if str(ex).strip() else "SystemExit"))
    except Exception as ex:
        fallos += 1
        print("  [!! ] %s -> %s: %s" % (g, type(ex).__name__, ex))
    plt.close("all")

if RECOGE:
    # una cadena por linea, en repr: varios rotulos llevan \n DENTRO, y escribirlos
    # en crudo los parte en dos y hace creer que son dos claves distintas.
    io.open("figs_es_CADENAS.txt", "w", encoding="utf-8").write(
        "\n".join(repr(s) for s in sorted(set(VISTAS))))
    print("")
    print("  cadenas dibujadas, sin repetir: %d -> figs_es_CADENAS.txt"
          % len(set(VISTAS)))
else:
    print("")
    print("  figuras escritas: %d" % len(SALIDAS))
    # CONTROL: cada guion tiene que haber dejado su PDF castellano.  Comparar el numero de guiones
    # con el de salidas no basta --- un guion puede escribir dos figuras y otro ninguna ---, asi que
    # se comprueba fichero a fichero y se nombran los que faltan.
    faltan = [g for g in GUIONES
              if g not in ("figs_es.py", "figlang.py")
              and not os.path.exists(g[:-3] + "_es.pdf")]
    if faltan:
        fallos += len(faltan)
        print("  SIN PDF CASTELLANO (%d): %s" % (len(faltan), ", ".join(faltan)))
    else:
        print("  todos los guiones dejaron su *_es.pdf")
    if SIN_TRADUCIR:
        print("  SIN TRADUCIR (%d):" % len(SIN_TRADUCIR))
        for s in sorted(SIN_TRADUCIR):
            print("     %r" % s)
    else:
        print("  ninguna cadena con palabras quedo sin traducir.")
print("  guiones que fallaron: %d" % fallos)
