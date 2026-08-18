r"""Diccionario ingles -> castellano de los rotulos de las figuras.

Las claves son EXACTAMENTE las cadenas que los guiones `fig_*.py` mandan a la capa
de texto de matplotlib; se recogieron ejecutandolos, no leyendo el codigo
(`python figs_es.py --recoge`).  Lo que no esta aqui no se traduce: la matematica
pura y los numeros pasan tal cual, y cualquier cadena con palabras que falte se
imprime al final de la corrida en lugar de colarse en ingles.

Varios rotulos son UNA cadena con un salto de linea dentro, no dos cadenas.  El
primer intento los partio por ese salto y quedaron sin traducir a pesar de estar
en el diccionario; van aqui con su salto conservado, porque la figura coloca cada
mitad en su propio renglon.

Las traducciones del hilo se mantienen cortas a proposito: `fig_thread` mide el
solape de las etiquetas en pixeles y se niega a dibujar si una se sale de su
perla, y el castellano es mas largo que el ingles casi siempre.
"""

# NOTACION, no idioma.  Son las palabras que aparecen dentro de `\mathrm{...}` en modo
# matematico y que el CUERPO del articulo castellano escribe igual --- `\mathrm{top}\Newt(...)`
# esta asi en la version espanola, y `mod`, `Newt` y `regular` tambien ---.  Se declaran aqui
# para que la lista de "sin traducir" de `figs_es.py` signifique algo: si queda vacia, no hay
# nada pendiente; una palabra nueva dentro de modo matematico aparecera ahi y habra que decidir.
NOTACION = {"regular", "Newt", "top", "mod"}

FIGLANG = {
    # ---- fig_galois -------------------------------------------------------
    r"$24$ points, $4$ orbits: one evaluation per colour, a saving of $\varphi(t)=6$":
        r"$24$ puntos, $4$ órbitas: una evaluación por color, un ahorro de $\varphi(t)=6$",
    r"$t=7$: the locus is $4$ Galois orbits":
        r"$t=7$: el lugar son $4$ órbitas de Galois",
    r"$t=5$: type $B$ moved by $\frac{t-1}{2}$ is type $C$":
        r"$t=5$: el tipo $B$ movido $\frac{t-1}{2}$ es el tipo $C$",
    r"every arrow lands on a circle; $\times$ marks a jump across the torus":
        r"toda flecha cae en un círculo; $\times$ marca un salto por el toro",
    r"size of the locus, rank $2$": r"tamaño del lugar, rango $2$",
    r"equal for odd $t$; unequal in all tested even cases":
        r"iguales para $t$ impar; desiguales en todos los casos pares probados",
    r"weight coordinates: $t=4$ is the one even case with nothing to compare, both sets empty":
        r"coordenadas de pesos: $t=4$ es el único caso par sin nada que comparar, los dos conjuntos vacíos",

    # ---- fig_cone3d -------------------------------------------------------
    r"$819$ lattice points satisfy the interlacing; $252$ survive the parity ($128$ with $+1$, $124$ with $-1$), and $567$ are killed":
        r"$819$ puntos del retículo cumplen el entrelazado; $252$ sobreviven a la paridad ($128$ con $+1$, $124$ con $-1$) y $567$ mueren",
    r"the cone, and the surviving sublattice": r"el cono, y el subretículo superviviente",
    r"edge on: the surviving planes alternate": r"de canto: los planos supervivientes alternan",
    r"the survivors lie on alternate planes $|\Lambda|+|\mu|$ even: adding a single box always leaves the set":
        r"los supervivientes están en planos alternos de $|\Lambda|+|\mu|$ par: añadir una casilla siempre sale del conjunto",

    # ---- fig_walls3d ------------------------------------------------------
    r"$C_3$ at $t=8$: the wall $2a_i\equiv0$ is invisible in the dual":
        r"$C_3$ en $t=8$: la pared $2a_i\equiv0$ es invisible en el dual",
    r"dies on $2a_i\equiv0$; the dual reads it as regular":
        r"muere en $2a_i\equiv0$; el dual lo lee como regular",
    r"on a wall both readings see": r"en una pared que ven las dos lecturas",

    # ---- fig_alcove -------------------------------------------------------
    r"even, $t=6$: $C_2$ at level $0$": r"par, $t=6$: $C_2$ en nivel $0$",
    r"odd, $t=5$: $B_2$ at level $1$": r"impar, $t=5$: $B_2$ en nivel $1$",
    r"1 point available": r"1 punto disponible",
    r"every generator of the fusion ideal dies": r"todo generador del ideal de fusión muere",
    r"value at the torsion element": r"valor en el elemento de torsión",
    "spin $\\omega_2$: not in the\nweight lattice of $SO_5$":
        "espín $\\omega_2$: fuera del\nretículo de pesos de $SO_5$",

    # ---- fig_filter -------------------------------------------------------
    r"even, $t=6$   (three wall families)": r"par, $t=6$   (tres familias de paredes)",
    r"odd, $t=5$   (two: the middle pair is absent)": r"impar, $t=5$   (dos: falta el par central)",
    r"49 of 136 survive": r"49 de 136 sobreviven",
    r"36 of 136 survive": r"36 de 136 sobreviven",
    r"surviving weights (%)": r"pesos supervivientes (%)",
    r"rank $m$ of the frozen factor": r"rango $m$ del factor congelado",
    r"the collapse with the rank": r"el colapso con el rango",
    r"box $\eta_1\leq 9$": r"caja $\eta_1\leq 9$",

    # ---- fig_transversal --------------------------------------------------
    r"a transversal: one index per nonzero class, none from class $0$":
        r"una transversal: un índice por clase no nula, ninguno de la clase $0$",
    r"one class unhit, and the whole numerator is zero":
        r"una clase sin tocar, y el numerador entero es cero",
    r"$n_j=1,2,1$, so $|\mathrm{supp}\,\nu|=2\prod_j n_j=4$.  Filled: the transversal $S_{\min}$ of lowest $V$, which carries the top weight.":
        r"$n_j=1,2,1$, luego $|\mathrm{supp}\,\nu|=2\prod_j n_j=4$.  En relleno: la transversal $S_{\min}$ de $V$ más bajo, que lleva el peso superior.",
    "$\\nu\\equiv0$: no transversal exists, so nothing is left to divide.\nRead off the residues of $V=2(\\Lambda+\\rho)$; no character is evaluated.":
        "$\\nu\\equiv0$: no existe transversal, así que no queda nada que dividir.\nSe lee de los residuos de $V=2(\\Lambda+\\rho)$; no se evalúa ningún carácter.",
    r"class $0$": r"clase $0$",

    # ---- fig_determinant --------------------------------------------------
    r"the fibre meets the support once": r"la fibra corta el soporte una vez",
    r"it meets it twice, and they cancel": r"lo corta dos veces, y se cancelan",

    # ---- fig_division -----------------------------------------------------
    r"$t=3$, $r=2$, $\Lambda=(6,6,6)$: $\nu$ is the finite difference of $c$":
        r"$t=3$, $r=2$, $\Lambda=(6,6,6)$: $\nu$ es la diferencia finita de $c$",
    r"filled: $\mathrm{supp}\,c$ ($25$ points).   ringed: $\mathrm{supp}\,\nu$ ($4$).   blue $+1$, red $-1$.   arrows: the $2^r$ shifts of one point":
        r"relleno: $\mathrm{supp}\,c$ ($25$ puntos).   anillado: $\mathrm{supp}\,\nu$ ($4$).   azul $+1$, rojo $-1$.   flechas: los $2^r$ desplazamientos de un punto",
    r"the quotient is far larger than the numerator, and still a unit":
        r"el cociente es mucho mayor que el numerador, y sigue siendo una unidad",
    r"every coefficient on both axes is $\pm1$, on all $186$ shapes; the diagonal is $y=x$":
        r"todo coeficiente en los dos ejes es $\pm1$, en las $186$ formas; la diagonal es $y=x$",
    r"$|\mathrm{supp}\,\nu|$   (numerator: proved $\{0,\pm1\}$)":
        r"$|\mathrm{supp}\,\nu|$   (numerador: probado $\{0,\pm1\}$)",
    r"$|\mathrm{supp}\,c|$   (quotient)": r"$|\mathrm{supp}\,c|$   (cociente)",

    # ---- fig_collapse -----------------------------------------------------
    r"$n=16$ forms;  largest single term $|B\tau|=798$":
        r"$n=16$ formas;  mayor término $|B\tau|=798$",
    r"every path wanders, every path lands on $\pm1$":
        r"todo recorrido divaga, todo recorrido aterriza en $\pm1$",
    "the band $|A_{\\mu_{\\max}}|\\leq1$ is the thin strip at zero;\nthe inset is the same strip, magnified":
        "la banda $|A_{\\mu_{\\max}}|\\leq1$ es la franja fina en cero;\nel recuadro es la misma franja, ampliada",
    r"partial sum of $B_{\eta,\mu_{\max}}\,\tau_t(\eta)$":
        r"suma parcial de $B_{\eta,\mu_{\max}}\,\tau_t(\eta)$",
    r"terms added, in decreasing order of $\eta$":
        r"términos añadidos, en orden decreciente de $\eta$",
    r"an order of magnitude, then one value": r"un orden de magnitud, y luego un valor",
    r"$\max$ partial sum": r"suma parcial $\max$",
    r"absolute value (log)": r"valor absoluto (log)",
    r"the landing": r"el aterrizaje",

    # ---- fig_law3d --------------------------------------------------------
    r"dominant vertex $(14,10,6)$": r"vértice dominante $(14,10,6)$",
    r"dominant vertex $(5,3,1)$": r"vértice dominante $(5,3,1)$",
    r"dominant vertex $(9,7,5)$": r"vértice dominante $(9,7,5)$",
    r"1080 exponents, 48 vertices": r"1080 exponentes, 48 vértices",
    r"2168 exponents, 48 vertices": r"2168 exponentes, 48 vértices",
    r"466 exponents, 48 vertices": r"466 exponentes, 48 vértices",
    r"and the subtracted vertex is $\sigma_r=(9,7,5)=2\rho_{C_r}+(t-1)$, the dominant vertex of the denominator, read in exponent coordinates":
        r"y el vértice restado es $\sigma_r=(9,7,5)=2\rho_{C_r}+(t-1)$, el vértice dominante del denominador, leído en coordenadas de exponentes",

    # ---- fig_residue ------------------------------------------------------
    r"all forms": r"todas las formas",
    r"vanishing": r"se anulan",
    r"not seen by the $t$-core": r"no lo ve el $t$-núcleo",
    r"how much is left": r"cuánto queda",
    r"forms ($t=6$, $r=2$, $\beta_i\leq13$)": r"formas ($t=6$, $r=2$, $\beta_i\leq13$)",
    r"forms": r"formas",
    r"the 36 $t$-cores of the occupied population":
        r"los 36 $t$-núcleos de la población ocupada",
    r"the vanishing concentrates on three cores":
        r"la anulación se concentra en tres núcleos",
    r"33 cores with no vanishing form at all":
        r"33 núcleos sin ninguna forma que se anule",
    r"core $(1,1,1)$": r"núcleo $(1,1,1)$",
    r"core $(2,1)$": r"núcleo $(2,1)$",
    r"$\lambda^{(0)}\neq\varnothing$: perfect, then not":
        r"$\lambda^{(0)}\neq\varnothing$: perfecto, y luego no",
    r"empty": r"vacío",

    # ---- fig_thread: cada rotulo es UNA cadena con su salto de linea -------
    "when does\nit vanish?": "¿cuándo se\nanula?",
    "reduce to $t\\in\\{1,2\\}$\nplus a specialisation":
        "reducir a $t\\in\\{1,2\\}$\ny una especialización",
    "what does the\nspecialisation do?": "¿qué hace la\nespecialización?",
    "it is a regularity\ncondition": "es una condición\nde regularidad",
    "regular in which\ngroup?  ($t$ even)": "¿regular en qué\ngrupo?  ($t$ par)",
    "regular in which\ngroup?  ($t$ odd)": "¿regular en qué\ngrupo?  ($t$ impar)",
    "twining: the input\ngoes virtual": "entrelazado: entra\nun virtual",
    "$\\tau^C_t$ with its sign,\nproved here": "$\\tau^C_t$ con su signo,\nprobado aquí",
    "$B_{R'}\\!\\downarrow B_{m'}\\!\\times\\! D_r$:\nordinary, genuine":
        "$B_{R'}\\!\\downarrow B_{m'}\\!\\times\\! D_r$:\nordinaria, genuino",
    "principal element:\ncited": "elemento principal:\ncitado",
    "where do the two\nfilters live?": "¿dónde viven los\ndos filtros?",
    # La clave de la izquierda cambio en fig_thread.py y esta entrada se quedo apuntando a la
    # cadena vieja: el nodo viajaba en INGLES dentro de la figura castellana, en la pagina 2.  Lo
    # vio _figs_es_audit.py y no figs_es.py, cuya heuristica de "matematica pura" es la que se
    # arreglo aparte.  Se deja la vieja por si algun guion la recupera.
    "minimal fusion:\na single point": "fusión mínima:\nun solo punto",
    "minimal fusion: one\ntensor-sector point": "fusión mínima: un punto\ndel sector tensorial",
    "one point: which\nweight is on top?": "un punto: ¿qué\npeso va arriba?",
    "how big is that\ncoefficient?": "¿cuánto vale ese\ncoeficiente?",
    "primitive:\n$q_t[M]=\\pm[0]$": "primitiva:\n$q_t[M]=\\pm[0]$",
    "and for $t$ odd?": "¿y para $t$ impar?",
    "and its filter?": "¿y su filtro?",
    "and the division\nby $\\Delta_t$?": "¿y la división\npor $\\Delta_t$?",
    "a sum along a\nprogression of step $2t$": "suma sobre una\nprogresión de paso $2t$",
    "and that sum?": "¿y esa suma?",
    "$c=\\pm\\epsilon_t\\det M$,\nan explicit $0/{\\pm}1$ matrix":
        "$c=\\pm\\epsilon_t\\det M$,\nmatriz $0/{\\pm}1$ explícita",
    "so what is left?": "¿qué queda?",
    "$M$ unimodular;\n(L2),(L3) still open": "$M$ unimodular;\n(L2),(L3) abiertas",
    "$\\nu$: a signed\ntransversal count": "$\\nu$: transversales\ncon signo",
    r"frozen": r"congelado",
    r"not dominant": r"no dominante",
    r"end": r"fin",
    r"even $t=2m+2$": r"$t$ par $=2m+2$",
    r"odd $t=2m+1$": r"$t$ impar $=2m+1$",
    r"$t$ even": r"$t$ par",
    r"$t$ odd": r"$t$ impar",

    # ---- etiquetas de leyenda ---------------------------------------------
    # No llegaban por ninguna de las funciones de alto nivel: matplotlib guarda
    # el `label=` y lo dibuja despues.  Aparecieron al enganchar `Text.set_text`.
    r"edges of the interlacing cone": r"aristas del cono del entrelazado",
    r"in the cone, killed by the parity": r"en el cono, muertos por la paridad",
    r"false positives of the discriminant": r"falsos positivos del discriminante",
    r"vanishing forms": r"formas que se anulan",
    r"the rest of the population": r"el resto de la población",
    r"lands on $A_{\mu_{\max}}=+1$": r"aterriza en $A_{\mu_{\max}}=+1$",
    r"lands on $A_{\mu_{\max}}=-1$": r"aterriza en $A_{\mu_{\max}}=-1$",
    r"survives, $\tau_t=+1$": r"sobrevive, $\tau_t=+1$",
    r"survives, $\tau_t=-1$": r"sobrevive, $\tau_t=-1$",
    r"on a wall, $\tau_t=0$": r"en una pared, $\tau_t=0$",
    r"wall": r"pared",
    r"wall absent for odd $t$": r"pared ausente para $t$ impar",
    r"type $B$ keeps": r"el tipo $B$ conserva",
    r"type $C$ keeps": r"el tipo $C$ conserva",
    r"proved here": r"probado aquí",
    r"not ours: cited": r"no es nuestro: citado",
    r"verified by computation": r"verificado por cálculo",
    r"open": r"abierto",

    # ---- fig_problems -----------------------------------------------------
    # Los titulos de los seis problemas los lee el guion del .aux INGLES, asi que aqui van sus
    # equivalentes castellanos, que son los que el lector de la edicion en castellano vera en el
    # cuerpo.  Si se renombra un problema, la cadena deja de casar y el informe lo lista.
    # Cadenas NO crudas: el salto tiene que ser un salto de verdad, no una barra y una ene.
    "what is left open, and what each open problem sits on":
        "qué queda abierto, y sobre qué se apoya cada problema",
    "even $t$": "$t$ par",
    "odd $t$": "$t$ impar",
    "external input": "aporte ajeno",
    "answered by measurement": "respondido midiendo",
    "the reduction to\n$t\\in\\{1,2\\}$": "la reducción a\n$t\\in\\{1,2\\}$",
    "both filters are minimal\nfusion projections":
        "proyecciones de fusión\nde nivel mínimo",
    "$\\mathrm{Newt}(N_\\delta)$ is the\n$C_r$ zonotope":
        "$\\mathrm{Newt}(N_\\delta)$ es el\nzonotopo de $C_r$",
    "OPEN: extremal\nprimitivity (Conj. 8.2)": "ABIERTO: primitividad\nextremal (Conj. 8.2)",
    "$\\det p_t=-1$: a twining\nto $C_m\\times C_r$":
        "$\\det p_t=-1$: entrelazado\na $C_m\\times C_r$",
    "$\\tau^C$ with its sign\n(the bialternant)": "$\\tau^C$ con su signo\n(el bialternante)",
    "a transversal count,\ntwo forbidden classes":
        "recuento transversal,\ndos clases prohibidas",
    "no per-$\\Lambda$ reduction:\n$a_\\Lambda$ has both signs":
        "sin reducción por $\\Lambda$:\n$a_\\Lambda$ toma los dos signos",
    "$\\det p_t=+1$: ordinary\n$B_{R'}\\to B_{m'}\\times D_r$":
        "$\\det p_t=+1$: ordinaria\n$B_{R'}\\to B_{m'}\\times D_r$",
    "$\\tau^B$ at a principal\nelement": "$\\tau^B$ en un elemento\nprincipal",
    "$\\nu$ is a signed\ntransversal count": "$\\nu$ es un recuento con\nsigno de transversales",
    "the division inverts:\n$\\Delta_t=\\psi^t(\\Delta_1)$":
        "la división se invierte:\n$\\Delta_t=\\psi^t(\\Delta_1)$",
    "$c=\\pm\\epsilon_t\\det M$,\nthen OPEN (L1)": "$c=\\pm\\epsilon_t\\det M$,\ny ABIERTO (L1)",
    "why primitive": "por qué primitiva",
    "where two filters agree": "dónde coinciden\ndos filtros",
    "a discriminant for the\nresidue": "un discriminante para\nel residuo",
    "is $G_t$ their induction\noperator?": "¿es $G_t$ su operador\nde inducción?",
    "the other classical types": "los otros tipos clásicos",
    "one diagram for both\nparities": "un diagrama para las\ndos paridades",
    "Problem 17.1": "Problema 17.1",
    "Problem 17.2": "Problema 17.2",
    "Problem 17.3": "Problema 17.3",
    "Problem 17.4": "Problema 17.4",
    "Problem 17.5": "Problema 17.5",
    "Problem 17.6": "Problema 17.6",

    # ---- identicas en los dos idiomas, declaradas para que no las liste ----
    # el informe de "sin traducir" como si se hubieran escapado.
    r"$\eta_1$ mod $t$": r"$\eta_1$ mod $t$",
    r"$\eta_2$ mod $t$": r"$\eta_2$ mod $t$",
    r"$a_1$ mod $t$": r"$a_1$ mod $t$",
    r"$a_2$ mod $t$": r"$a_2$ mod $t$",
}
