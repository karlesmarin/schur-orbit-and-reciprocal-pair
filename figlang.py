# -*- coding: utf-8 -*-
r"""Diccionario de los rotulos dibujados dentro de las figuras, ingles -> castellano.

Lo usa `figs_es.py`, que intercepta la capa de texto de matplotlib y pasa por aqui cada cadena
antes de dibujarla.  Los guiones no se tocan ni se duplican, de modo que todos sus controles
siguen corriendo.

Reglas seguidas:

  * La matematica no se traduce.  Una cadena que entre `$...$` no tiene tres letras seguidas
    fuera pasa tal cual, y `figs_es.py` la deja pasar sin avisar.
  * El vocabulario es el de `orbit_pair_es.tex`, no el que a uno le salga: *lugar de anulacion*,
    *maximizador*, *transversal*, *autocomplementaria*, *etiqueta*, *rachas*, *input externo*.
  * Lo que queda fuera del diccionario y lleva palabras se imprime al final de la corrida, para
    que un rotulo nuevo no se cuele en ingles sin que nadie lo vea.

Autores: Carles Marin + Claude (AI assistant).
"""

FIGLANG = {
    # --- fig_thread: el hilo de siete preguntas ---------------------------------------------
    "What does a Schur polynomial\ndo on a full orbit $\\mu_t$?":
        "¿Qué hace un polinomio de Schur\nen una órbita completa $\\mu_t$?",
    "Littlewood: it vanishes exactly\nwhen the $t$-core is nonempty":
        "Littlewood: se anula exactamente\ncuando el $t$-núcleo no es vacío",
    "And if one free reciprocal\npair joins that orbit?":
        "¿Y si un par recíproco libre\nse une a esa órbita?",
    "then it factorises: the value is\na triple $(d_1,d_2,d_3)$ and a sign":
        "entonces factoriza: el valor es\nuna terna $(d_1,d_2,d_3)$ y un signo",
    "Is that triple all the\nshape leaves behind?":
        "¿Es esa terna todo lo que\nla forma deja atrás?",
    # El hilo es la figura mas apretada del paper: en ingles la respuesta de la perla 3 ya
    # termina justo donde empieza la pregunta de la perla 1.  "de cualquier tamano" son nueve
    # caracteres mas y la cruzaba.  "cualquier particion" dice lo mismo y cabe.
    "yes — a partition of any size\nis compressed to three integers":
        "sí — cualquier partición\nse comprime a tres enteros",
    "Does the factorising\nsurvive a second pair?":
        "¿Sobrevive la factorización\na un segundo par?",
    "no. But at $t=2$ the set where\nit vanishes does, for every $r$":
        "no. Pero en $t=2$ sí el conjunto\ndonde se anula, para todo $r$",
    "At $t=2$, which $\\lambda$\nvanish for every $r$?":
        "En $t=2$, ¿qué $\\lambda$\nse anulan para todo $r$?",
    "two families, and nothing else":
        "dos familias, y nada más",
    "And for every $t$\nat once?":
        "¿Y para todo $t$\na la vez?",
    "one reflection of the excess part\nforces it — for all $t$ and $r$":
        "una reflexión de la parte de exceso\nlo fuerza — para todo $t$ y $r$",
    "So what is\nstill missing?":
        "¿Y qué\nfalta todavía?",
    "one set: is the common omitted\npart $g_{\\mathrm{com}}$ symmetric too?":
        "un conjunto: ¿es simétrica también\nla parte omitida común $g_{\\mathrm{com}}$?",
    "Why did one pair\nbehave so differently?":
        "¿Por qué un solo par\nse portaba tan distinto?",
    "there the folded expansion is\n$\\pm$ genuine; from two on it need not be":
        "allí la expansión plegada es\n$\\pm$ genuina; de dos en adelante no tiene por qué",
    "classical": "clásico",
    "proved here": "demostrado aquí",
    "open": "abierto",

    # --- fig_alphabet -----------------------------------------------------------------------
    "(a)  $t=4$, one pair, $|z|=1$": "(a)  $t=4$, un par, $|z|=1$",
    "on the circle, and free of the orbit": "sobre la circunferencia, y libre de la órbita",
    "(b)  $t=4$, one pair, $z$ real": "(b)  $t=4$, un par, $z$ real",
    "$z$ outside, $z^{-1}$ inside: a free direction":
        "$z$ fuera, $z^{-1}$ dentro: una dirección libre",
    "(c)  $t=2$, two pairs": "(c)  $t=2$, dos pares",
    "an $r$-dimensional free torus, inside $O(6)^{-}$":
        "un toro libre de dimensión $r$, dentro de $O(6)^{-}$",
    "the frozen orbit $\\mu_t$": "la órbita congelada $\\mu_t$",
    "a free reciprocal pair $z,\\,z^{-1}$": "un par recíproco libre $z,\\,z^{-1}$",
    # Esta se colaba en ingles en los TRES paneles de la edicion castellana, y sin aviso: la
    # cadena entera va entre $...$, y la regla de `_solo_matematica` miraba solo FUERA del modo
    # matematico, asi que `\mathrm{alphabet}` era invisible.  Corregida la regla el 18 de agosto,
    # ahora tambien mira dentro de \mathrm{} y \text{}.  Se listan las dos ramas del signo aunque
    # esta figura solo produzca la de -1, para que un cambio de parametros no la vuelva a colar.
    "$\\prod(\\mathrm{alphabet})=(-1)^{t-1}=-1$":
        "$\\prod(\\mathrm{alfabeto})=(-1)^{t-1}=-1$",
    "$\\prod(\\mathrm{alphabet})=(-1)^{t-1}=+1$":
        "$\\prod(\\mathrm{alfabeto})=(-1)^{t-1}=+1$",
    # NOTACION, no idioma: `inv` es el numero de inversiones y se escribe igual en castellano.
    # Va aqui, con su cadena identica, porque la regla nueva la detecta --- correctamente, es una
    # palabra en modo matematico --- y sin esta entrada apareceria en cada corrida en la lista de
    # "sin traducir", que es la lista que hay que mantener vacia para que sirva de algo.
    "$M_S=(-1)^{\\mathrm{inv}(b_S)}V$": "$M_S=(-1)^{\\mathrm{inv}(b_S)}V$",

    # --- fig_beta ---------------------------------------------------------------------------
    "(b) $\\beta_j=\\lambda_j+N-j$, strictly decreasing":
        "(b) $\\beta_j=\\lambda_j+N-j$, estrictamente decreciente",
    "(c) the classes mod $t=4$, and the excess":
        "(c) las clases módulo $t=4$, y el exceso",
    "the dashed classes are the excess $\\mathcal{E}$; $|\\mathcal{E}|=e=2$":
        "las clases a trazos son el exceso $\\mathcal{E}$; $|\\mathcal{E}|=e=2$",
    "residue $0$": "residuo $0$",
    "residue $1$": "residuo $1$",
    "residue $2$": "residuo $2$",
    "residue $3$": "residuo $3$",

    # --- fig_plane --------------------------------------------------------------------------
    "$t$  (size of the frozen orbit $\\mu_t$)":
        "$t$  (tamaño de la órbita congelada $\\mu_t$)",
    "$r$  (free pairs)": "$r$  (pares libres)",
    "the value itself (Thm 3.1)": "el valor mismo (Teo. 3.1)",
    "criterion, no external input (Cor 8.22)": "criterio, sin input externo (Cor. 8.22)",
    "criterion, modulo one input (Thm 8.6)": "criterio, módulo un input (Teo. 8.6)",
    "one implication proved (Thm 8.35);\nconverse conjectural (Conj 8.43)":
        "una implicación demostrada (Teo. 8.35);\nrecíproco conjetural (Conj. 8.43)",
    "where the zero locus is known, and how":
        "dónde se conoce el lugar de anulación, y cómo",
    "columns of odd $t$ are settled outright;\nthe row $r=1$ is the only place the value,\n"
    "and not merely its vanishing, is known":
        "las columnas de $t$ impar quedan resueltas;\nla fila $r=1$ es el único sitio donde se\n"
        "conoce el valor y no sólo su anulación",

    # --- fig_laplace ------------------------------------------------------------------------
    "Laplace along the $t$\nfrozen rows": "Laplace por las $t$\nfilas congeladas",
    "Lemma 4.1: zero unless the $t$": "Lema 4.1: cero salvo que las $t$",
    "residues on $S$ are distinct": "residuos sobre $S$ sean distintos",
    "the free $2\\times2$ block": "el bloque libre $2\\times2$",
    # Este va MEDIDO, no traducido a ojo: es el rotulo de la izquierda de la fila superior de
    # fig_laplace, y el de la derecha empieza donde el ingles termina.  La traduccion literal
    # ---"por cada eleccion de las dos columnas libres"--- es 13 caracteres mas larga y se comia
    # al vecino.  "por cada par de columnas libres" dice lo mismo en la anchura del ingles.
    "one term per choice of the two free columns $U=\\{j_1<j_2\\}$":
        "un término por cada par de columnas libres $U=\\{j_1<j_2\\}$",
    "and every choice falls into one of three types":
        "y toda elección cae en uno de tres tipos",
    "profile $(2,2,1,1)$": "perfil $(2,2,1,1)$",
    "profile $(3,1,1,1)$": "perfil $(3,1,1,1)$",
    "profile $(3,2,1,0)$": "perfil $(3,2,1,0)$",
    "two classes of size two": "dos clases de tamaño dos",
    "one class of size three": "una clase de tamaño tres",
    "a residue class is empty": "una clase de residuos es vacía",
    "terms": "términos",
    "value is zero": "el valor es cero",

    # --- fig_image / fig_fibres -------------------------------------------------------------
    "$t=3$:   860 partitions $\\longrightarrow$ 96 points":
        "$t=3$:   860 particiones $\\longrightarrow$ 96 puntos",
    "$t=4$:   883 partitions $\\longrightarrow$ 73 points":
        "$t=4$:   883 particiones $\\longrightarrow$ 73 puntos",
    "$t=3$:   160 fibres over $|\\lambda|\\leq 20$":
        "$t=3$:   160 fibras sobre $|\\lambda|\\leq 20$",
    "$t=4$:   121 fibres over $|\\lambda|\\leq 20$":
        "$t=4$:   121 fibras sobre $|\\lambda|\\leq 20$",
    "floor grid spacing $= t = 3$": "paso de la rejilla del suelo $= t = 3$",
    "floor grid spacing $= t = 4$": "paso de la rejilla del suelo $= t = 4$",
    "partitions collapsing onto one point": "particiones que colapsan en un punto",
    "number of partitions": "número de particiones",
    "the coupling term $d_3$": "el término de acoplamiento $d_3$",
    "accent: $I_{3}=(3,3,2,+1)$, 26 partitions of sizes 0–20":
        "acento: $I_{3}=(3,3,2,+1)$, 26 particiones de tamaños 0–20",
    "accent: $I_{4}=(4,4,2,+1)$, 36 partitions of sizes 0–20":
        "acento: $I_{4}=(4,4,2,+1)$, 36 particiones de tamaños 0–20",

    # --- fig_locus --------------------------------------------------------------------------
    "$t$-cores  ($6$)": "$t$-núcleos  ($6$)",
    "$t$-cores  ($10$)": "$t$-núcleos  ($10$)",
    "$t$-cores  ($21$)": "$t$-núcleos  ($21$)",
    "extra family  ($0$)": "familia adicional  ($0$)",
    "extra family  ($2$)": "familia adicional  ($2$)",
    "extra family  ($3$)": "familia adicional  ($3$)",
    "neither": "ninguna",
    "$t=3$   (odd: no extra family)": "$t=3$   (impar: sin familia adicional)",

    # --- fig_signed -------------------------------------------------------------------------
    "columns $c$ of the shape $\\lambda=(c^{\\,a})$":
        "columnas $c$ de la forma $\\lambda=(c^{\\,a})$",
    "rows $a$": "filas $a$",
    "signed count": "recuento con signo",

    # --- fig_runs ---------------------------------------------------------------------------
    "every $k$ even: the runs alternate": "todo $k$ par: las rachas alternan",
    "every $k$ odd: no two are adjacent": "todo $k$ impar: no hay dos adyacentes",

    # --- fig_zeros --------------------------------------------------------------------------
    "branch (a)": "rama (a)",
    "branch (b)": "rama (b)",
    "branch (a): $\\beta$ of constant parity": "rama (a): $\\beta$ de paridad constante",
    "branch (b): self-complementary, odd width":
        "rama (b): autocomplementaria, de anchura impar",
    "branch (b) only at $|\\lambda|=w(r{+}1)$, $w$ odd":
        "rama (b) sólo en $|\\lambda|=w(r{+}1)$, $w$ impar",
    "all shapes with at most $N$ parts": "todas las formas con a lo sumo $N$ partes",
    "$57$ zeros of $15226$ shapes": "$57$ ceros de $15226$ formas",
    "$118$ zeros of $8547$ shapes": "$118$ ceros de $8547$ formas",
    "$319$ zeros of $2724$ shapes": "$319$ ceros de $2724$ formas",

    # --- fig_involution ---------------------------------------------------------------------
    "$\\lambda=(3,3,2,1)$,  $\\beta=(8,7,5,3,1,0)$,  $c=8$ --- self-complementary, "
    "$w=3$ odd: the character vanishes":
        "$\\lambda=(3,3,2,1)$,  $\\beta=(8,7,5,3,1,0)$,  $c=8$ --- autocomplementaria, "
        "$w=3$ impar: el carácter se anula",
    "$\\lambda=(3,3,2,2)$,  $\\beta=(8,7,5,4,1,0)$,  $c=8$ --- one box away: no pairing, "
    "and it does not":
        "$\\lambda=(3,3,2,2)$,  $\\beta=(8,7,5,4,1,0)$,  $c=8$ --- a una casilla: sin "
        "emparejamiento, y no se anula",
    "$x\\mapsto c-x$  closes on $\\beta$": "$x\\mapsto c-x$  cierra sobre $\\beta$",
    "$x\\mapsto c-x$  misses $\\beta$": "$x\\mapsto c-x$  se sale de $\\beta$",
    "$x\\mapsto c-x$, landing in $\\beta$": "$x\\mapsto c-x$, cayendo en $\\beta$",
    "landing outside $\\beta$": "cayendo fuera de $\\beta$",
    "$\\beta_j$ even": "$\\beta_j$ par",
    "$\\beta_j$ odd": "$\\beta_j$ impar",

    # --- fig_reduction ----------------------------------------------------------------------
    "number of labels $\\nu$": "número de etiquetas $\\nu$",
    "standard,\n$\\nu'_1<N/2$": "estándar,\n$\\nu'_1<N/2$",
    "standard,\n$\\nu'_1=N/2$": "estándar,\n$\\nu'_1=N/2$",
    "standard,\n$\\nu'_1>N/2$": "estándar,\n$\\nu'_1>N/2$",
    "non-standard,\n$\\nu'_1+\\nu'_2>N$": "no estándar,\n$\\nu'_1+\\nu'_2>N$",
    "survives: an element of the basis": "sobrevive: un elemento de la base",
    "folds onto the basis, as $\\pm\\,o_\\mu(A)$": "se pliega sobre la base, como $\\pm\\,o_\\mu(A)$",
    "vanishes identically": "se anula idénticamente",
    "$r=1$,  $N=4$:  97 labels, basis of 10": "$r=1$,  $N=4$:  97 etiquetas, base de 10",
    "$r=2$,  $N=6$:  67 labels, basis of 25": "$r=2$,  $N=6$:  67 etiquetas, base de 25",

    # --- fig_phase --------------------------------------------------------------------------
    "vanishes at the endpoint only": "se anula sólo en el extremo",
    "$r=1$,  $N=4$:  no such shape in 1292": "$r=1$,  $N=4$:  ninguna forma así entre 1292",
    "$r=2$,  $N=6$:  4 of 1513 shapes": "$r=2$,  $N=6$:  4 de 1513 formas",
    "$r=3$,  $N=8$:  3 of 795 shapes": "$r=3$,  $N=8$:  3 de 795 formas",

    # --- fig_increments ---------------------------------------------------------------------
    "increment $\\Delta_i(k)=c_{i,k}+c_{i,k+1}$": "incremento $\\Delta_i(k)=c_{i,k}+c_{i,k+1}$",
    "tie: $\\tau=8$": "empate: $\\tau=8$",
    "$|G|\\leq2$, tie in\nclasses $i,\\ i{+}t/2$": "$|G|\\leq2$, empate en las\nclases $i,\\ i{+}t/2$",
    "increments and the\n2-torsion of $\\mathbb{Z}/t$":
        "los incrementos y la\n2-torsión de $\\mathbb{Z}/t$",
    "the two maximisers, class by class:\nthey differ in the tied classes only":
        "los dos maximizadores, clase a clase:\nsólo difieren en las clases empatadas",
    "top half $H$": "mitad superior $H$",
    "bottom half $L$": "mitad inferior $L$",
    "omitted $g_i$": "$g_i$ omitido",

    # --- fig_reflection ---------------------------------------------------------------------
    "(a)  $[\\det]_{D_1}=0$:  $T_{\\mathrm{b}}=\\tau-T_{\\mathrm{a}}$, and both extremes of "
    "$\\mathcal{S}$ avoid $g_{\\mathrm{com}}$":
        "(a)  $[\\det]_{D_1}=0$:  $T_{\\mathrm{b}}=\\tau-T_{\\mathrm{a}}$, y los dos extremos de "
        "$\\mathcal{S}$ evitan $g_{\\mathrm{com}}$",
    "(b)  $|G|=2$ but $[\\det]_{D_1}\\neq0$:  no reflection, and $\\min\\mathcal{S}$ falls in "
    "$g_{\\mathrm{com}}$":
        "(b)  $|G|=2$ pero $[\\det]_{D_1}\\neq0$:  sin reflexión, y $\\min\\mathcal{S}$ cae en "
        "$g_{\\mathrm{com}}$",
    "$T_{\\mathrm{b}}=\\tau-T_{\\mathrm{a}}$:\nthe two reflect":
        "$T_{\\mathrm{b}}=\\tau-T_{\\mathrm{a}}$:\nlos dos se reflejan",
    "$v\\leftrightarrow\\tau-v$ inside $\\mathcal{S}$":
        "$v\\leftrightarrow\\tau-v$ dentro de $\\mathcal{S}$",
    "$g_{\\mathrm{com}}$ (omitted by both)": "$g_{\\mathrm{com}}$ (omitido por los dos)",
    "the tie $\\{p_1,p_2,q_1,q_2\\}$": "el empate $\\{p_1,p_2,q_1,q_2\\}$",
    "the extremes of $\\mathcal{S}$\navoid $g_{\\mathrm{com}}$":
        "los extremos de $\\mathcal{S}$\nevitan $g_{\\mathrm{com}}$",

    # --- fig_pairing ------------------------------------------------------------------------
    "(a) both clauses hold": "(a) se cumplen las dos cláusulas",
    "(b) symmetric, but no increment $=C$": "(b) simétrico, pero ningún incremento $=C$",
    # Los dos bloques de rotulos de fig_pairing van uno a cada lado del panel, y en ingles se
    # rozan sin tocarse.  La segunda linea es la que manda la anchura: en castellano hay que
    # dejarla en la del ingles o invade el bloque de enfrente.
    "16 transversals, 0 fixed\nall pair up with opposite sign":
        "16 transversales, 0 fijas\nse emparejan con signo opuesto",
    "12 transversals, 2 fixed\na fixed one cannot cancel":
        "12 transversales, 2 fijas\nuna fija no se cancela",
    "class $0$: 2 elements (even)\nclass $2$: 2 elements (even)":
        "clase $0$: 2 elementos (par)\nclase $2$: 2 elementos (par)",
    "class $0$: 3 elements (odd)": "clase $0$: 3 elementos (impar)",
    "the reflection acts on\nevery transversal": "la reflexión actúa sobre\ntoda transversal",
    "$w(g^\\dagger)=-w(g)$\nand $A(C-T)=A(T)$": "$w(g^\\dagger)=-w(g)$\ny $A(C-T)=A(T)$",
    "$C-\\mathcal{S}=\\mathcal{S}$ and\nsome $\\Delta_i(k)=C$":
        "$C-\\mathcal{S}=\\mathcal{S}$ y\nalgún $\\Delta_i(k)=C$",

    # --- fig_virtual ------------------------------------------------------------------------
    "(a)  $\\lambda=(9,8,8,6,0,0)$: a genuine character, all $a_{\\lambda\\mu}>0$":
        "(a)  $\\lambda=(9,8,8,6,0,0)$: un carácter genuino, todos los $a_{\\lambda\\mu}>0$",
    "(b)  $\\lambda=(9,9,8,6,0,0)$: properly virtual":
        "(b)  $\\lambda=(9,9,8,6,0,0)$: propiamente virtual",
    "sign footprint, on the floor": "huella del signo, sobre el suelo",
    "the cut at rank $r=2$": "el corte en el rango $r=2$",

    # --- fig_map ----------------------------------------------------------------------------
    "three routes to the zero locus, and what each one reaches":
        "tres rutas al lugar de anulación, y hasta dónde llega cada una",
    "the criterion, every $r$\n(Theorem 8.6)": "el criterio, todo $r$\n(Teorema 8.6)",
    "vanishing, every $t$\nand $r$ (Thm 8.35)": "anulación, todo $t$\ny $r$ (Teo. 8.35)",
    "Littlewood's reduction\nto the $C_\\mu$": "la reducción de Littlewood\na los $C_\\mu$",
    "the dictionary:\ntwo Schur factors": "el diccionario:\ndos factores de Schur",
    "external:\nrigidity of $s_\\lambda s_\\mu$": "externo:\nrigidez de $s_\\lambda s_\\mu$",
    "external input": "input externo",
    "a certificate for it\n(refuted)": "un certificado para ello\n(refutado)",
    "a counterexample, so proved": "un contraejemplo, luego demostrado",
    "OPEN: the converse\n(Conjecture 8.43)": "ABIERTO: el recíproco\n(Conjetura 8.43)",
    "OPEN: explain the\n$C_\\mu$ (Problem 10.5)": "ABIERTO: explicar los\n$C_\\mu$ (Problema 10.5)",
}
