# -*- coding: utf-8 -*-
r"""Prueba de la regla `_solo_matematica` de `figs_es.py`, con el caso que la rompio.

El 18 de agosto de 2026 la formula

    $\prod(\mathrm{alphabet})=(-1)^{t-1}=-1$

se dibujo en INGLES en los tres paneles de `fig_alphabet_es.pdf` y NO aparecio en la lista de
"sin traducir": la regla borraba el modo matematico entero y miraba solo lo de fuera, asi que una
cadena enteramente matematica se declaraba intraducible por construccion.  Estaba asi en la copia
publicada, arXiv:2608.09619v2 pagina 3.

Esta prueba fija el comportamiento en las dos direcciones: que el texto en modo matematico se
detecte, y que la matematica de verdad se siga dejando pasar --- que es lo que hace util la regla.

    python test_figs_es_heuristic.py

Autores: Carles Marin + Claude (AI assistant).
"""
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

_TEXTO_EN_MATE = re.compile(
    r"\\(?:math(?:rm|it|sf|tt|bf|cal)|text(?:rm|it|bf|sf|tt|normal)?|operatorname\*?)"
    r"\s*\{([^{}]*)\}"
)


def _solo_matematica(s):
    fuera = re.sub(r"\$[^$]*\$", "", s)
    dentro = " ".join(_TEXTO_EN_MATE.findall(s))
    return not re.search(r"[A-Za-z]{3,}", fuera + " " + dentro)


# (cadena, esperado, por que)
CASOS = [
    # --- lo que hay que DETECTAR (esperado False = hay algo que traducir) ---------------------
    (r"$\prod(\mathrm{alphabet})=(-1)^{t-1}=-1$", False, "EL CASO: la que se colo tres veces"),
    (r"$\operatorname{rank}(M)=r$", False, "operatorname con palabra"),
    (r"$\text{odd width}$", False, "text con palabra"),
    (r"$\mathrm{alfabeto}$", False, "tambien el castellano: la regla mira palabras, no idioma"),
    (r"the frozen orbit $\mu_t$", False, "palabra fuera del modo matematico"),
    (r"$C-\mathcal{S}=\mathcal{S}$ and some $\Delta_i(k)=C$", False, "palabra entre dos formulas"),

    # --- lo que hay que DEJAR PASAR (esperado True = matematica pura) -------------------------
    (r"$\lambda+\rho$", True, "comandos de LaTeX no son palabras"),
    (r"$x_1$", True, "trivial"),
    (r"$\mathcal{S}$", True, "mathcal de una letra sigue siendo notacion"),
    (r"$\sum_{i=1}^{n}\binom{n}{i}$", True, "sum y binom son comandos, no texto"),
    (r"$(d_1,d_2,d_3)$", True, "una terna"),
    (r"$O(6)^{-}$", True, "un grupo"),
]

fallos = 0
for s, esperado, nota in CASOS:
    got = _solo_matematica(s)
    ok = got == esperado
    fallos += 0 if ok else 1
    print("  %-5s solo_matematica=%-5s esperado=%-5s  %s" % ("ok" if ok else "FALLA", got, esperado, nota))

print()
if fallos:
    print("FALLAN %d de %d." % (fallos, len(CASOS)))
    sys.exit(1)
print("PASA: %d de %d." % (len(CASOS), len(CASOS)))
