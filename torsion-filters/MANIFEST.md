# Ancillary files for *Torsion filters at a root-of-unity orbit*

Everything the paper appeals to, in three directories.

* `gates/` --- the computations behind the claims of the verification section. Each script
  is self-contained and prints its own controls; where it produced a machine-readable
  result the archived run is beside it as `*_OUT.txt` and `*_DUMP.json`.
* `figures/` --- the scripts that draw the figures. Each recomputes its data and refuses to
  draw if its controls fail, so a figure in this paper cannot disagree with its own data.
* `audits/` --- the checks run on the manuscript itself: citations, cross-table
  consistency, formula typography, ambiguous terms, and the numbers in the text against
  the dumps.

Scripts ending in `.sage` need SageMath; the rest are plain Python 3 and need only the
standard library, except the figure scripts, which need matplotlib.

The header of every script states what it tests and, where it applies, the *decoy* that was
run against it --- a deliberately wrong variant which had to fail. That discipline is the
reason the verification tables can distinguish a result from an agreement.

## Gates (241)

| script | what it does | archived output |
|---|---|---|
| `A1_is_purbhoo_vw.sage` | A1 NO ES NUESTRA: es el Teorema 2.5 de Purbhoo-van Willigenburg (2008).  Aqui esta el diccionario, | yes |
| `E4_operator.sage` | E^{(4)} COMO OPERADOR, independiente de beta.   15 de agosto de 2026. | yes |
| `G2_size.py` | EL TAMANO DE G2, que es la PRECONDICION de transplantar la seccion 5 un piso mas abajo. | yes |
| `G3_size.py` | EL TAMANO DEL TERCER ESTRATO.  Se puede controlar como se controlo /G/ <= 2? | yes |
| `Gbot_anatomy.sage` | LA ANATOMIA DE Gbot.  Arriba el maximizador era casi unico; abajo hay hasta 16.  Que los ordena? | yes |
| `L1_subset_formula.sage` | (L1) POR SUBCONJUNTOS: la formula que la prueba, puesta a prueba.   16 de agosto de 2026. | yes |
| `L3_caja_grande.sage` | (L3) DONDE PUEDE FALLAR.   16 de agosto de 2026. | yes |
| `L3_libre_de_multiplicidad.sage` | .ES LIBRE DE MULTIPLICIDAD LA RESTRICCION IMPAR?   16 de agosto de 2026. | yes |
| `L3_multiplicidad.sage` | (L3): EL Lambda QUE APORTA TIENE MULTIPLICIDAD UNO.   16 de agosto de 2026. | yes |
| `_control.py` | _control.py -- LIBRERIA, no gate.  Controles que se niegan a certificar cuando no pueden fallar. | -- |
| `_pf_gkrs.sage` | GKRS PARA (L1).   16 de agosto de 2026. | yes |
| `_preflight.sage` | CONTROL EXTERNO del criterio del 14 de agosto, por VIA INDEPENDIENTE, con cascada y checkpoint. | -- |
| `_probe2.sage` | (sin titulo) | -- |
| `_probe_313.py` | .DE DONDE SALEN LOS 313 PARES?   17 de agosto de 2026. | yes |
| `_probe_afin.py` | .ES EL PLEGADO AFIN?  DESCOMPOSICION DEL SIGNO EN LOS PARES QUE SE CANCELAN. | yes |
| `_probe_bd.sage` | sonda: ¿existe la regla de ramificacion B_R' -> B_m' x D_r en esta imagen de Sage? | -- |
| `_probe_branch.sage` | sonda: ¿existe la regla Sp_8 -> Sp_4 x Sp_4 en esta imagen? | -- |
| `_probe_cancelacion.py` | .POR QUE SE CANCELAN?  LA ESTRUCTURA DE LAS PROGRESIONES CON VARIOS TERMINOS. | yes |
| `_probe_cero_uno.py` | .ES M EQUIVALENTE A UNA MATRIZ 0/1? | yes |
| `_probe_conceptos.sage` | DOS CONCEPTOS QUE EL PAPER USA SIN COMPROBAR.   16 de agosto de 2026. | -- |
| `_probe_determinante.py` | LA DESCOMPOSICION QUE FALTABA: c(X) COMO UN DETERMINANTE. | yes |
| `_probe_dossignos.py` | .PUEDEN LOS DOS SIGNOS REALIZAR UNA MISMA RANURA LIBRE? | yes |
| `_probe_even_top.py` | .Por que falla E7 en el par?  El argumento del impar no usa la paridad: v es estrictamente | yes |
| `_probe_fix_tautology.py` | .QUE ARREGLO SIRVE PARA SOCRATES?  UNA SOLA PRUEBA, TODOS LOS CANDIDATOS. | yes |
| `_probe_gammaC.sage` | LAS DOS PROPUESTAS DE LA VUELTA 29, verificadas antes de escribirlas.   16 de agosto de 2026. | -- |
| `_probe_grupo.py` | OTRA DESCOMPOSICION DE W^1: LEER EL ELEMENTO DE GRUPO, NO ADIVINAR EL FACTOR. | yes |
| `_probe_involucion.py` | LA INVOLUCION QUE CANCELA.  .Que mueve un termino de la progresion al siguiente? | yes |
| `_probe_multiplete.py` | .ES NUESTRO nu UN MULTIPLETE GKRS EN EL SENTIDO DE LANDWEBER-SJAMAAR? | yes |
| `_probe_quiralidad.py` | INTERROGATORIO DE LA QUIRALIDAD.  .Que es simetrico y que antisimetrico? | yes |
| `_probe_tauC.sage` | .Tiene el filtro de tipo C una constante de normalizacion, como el de tipo B? | -- |
| `_probe_timing.sage` | (sin titulo) | -- |
| `_probe_weyl_doble.py` | .CUANTA CANCELACION ES SOLO ANTISIMETRIA DE WEYL? | yes |
| `_reprint_square.py` | Reimprime la tabla del cuadrado conmutativo DESDE commuting_square_DUMP.json. | -- |
| `_witness_quick.sage` | comprobacion puntual del testigo: Sp_4, t=6, lambda=(1,0).  Tres rutas. | -- |
| `a1_proof.py` | ============================================================================================ | yes |
| `affine_A_check.py` | ¿SE CUMPLE LA HIPOTESIS (A) DE LA CONSULTA EXTERNA?  -- 15 de agosto de 2026. | yes |
| `all_A_mu.sage` | ¿VALE /A_mu/ <= 1 PARA TODO mu, O SOLO PARA EL MAXIMO?   15 de agosto de 2026. | yes |
| `atoms_distance.sage` | LOS ATOMOS (nu,Q) Y SUS DISTANCIAS.   15 de agosto de 2026. | yes |
| `ayyer_component.sage` | ayyer_component.sage | yes |
| `band_law.py` | LA BANDA +2s: de donde sale el valor de saturacion.  13 de agosto de 2026. | yes |
| `bottom_anatomy.sage` | EL ESTRATO DE ABAJO ES MAS SIMPLE QUE EL DE ARRIBA: no hace falta rigidez de Schur, solo signos. | yes |
| `bottom_sees_gcom.py` | EL ESTRATO DE ABAJO ES JUSTO LO QUE VE g_com.  Hipotesis del 12 de agosto de 2026. | yes |
| `bottom_stratum.sage` | LA FILTRACION POR GRADOS TIENE DOS EXTREMOS.  Solo habiamos mirado uno. | yes |
| `branch_filter.sage` | LA RUTA DEL BRANCHING FILTRADO.   Phi_{t,r} = sum_mu ( sum_eta B_{eta,mu} tau_t(eta) ) sp_mu. | yes |
| `branching_B.sage` | The general-t (eq:Cmu): the freezing is a branching matrix onto the Sp(2r) basis. | yes |
| `c_eq_tau.py` | C = tau: EL HUECO DEL 12, y su reduccion a H*.  13 de agosto de 2026. | yes |
| `cancel_traj_gcom.py` | ¿ESTA LA CANCELACION GENUINA CONTROLADA POR g_com?  El test decisivo -- 15 de agosto de 2026. | yes |
| `cancel_vs_gcom.py` | LA SATURACION CONTRA g_com -- 15 de agosto de 2026. | yes |
| `cancel_vs_vanishing.py` | ¿ES LA SATURACION PARTE DE LA HISTORIA DE conj:gcom, O ES OTRO PAPER? -- 15 de agosto de 2026. | yes |
| `check_24in.sage` | LAS DOS AFIRMACIONES DE SU VUELTA 24, VERIFICADAS ANTES DE TOCAR NADA.   16 de agosto de 2026. | yes |
| `check_e6.py` | SU CONTRAEJEMPLO A "e <= 4", CON NUESTRA PROPIA ANATOMIA.   16 de agosto de 2026. | yes |
| `check_fig2.py` | ¿DIBUJA LA FIGURA 2 LA REGLA EQUIVOCADA EN EL PANEL IMPAR?   16 de agosto de 2026. | yes |
| `clases_t3r2.py` | LAS CLASES DE MATRICES A t=3, r=2, IMPRESAS. | yes |
| `closed_form_prof.py` | LA PROFUNDIDAD EN FORMA CERRADA, y la hipotesis exacta.  13 de agosto de 2026. | yes |
| `collision_graph.py` | EL GRAFO DE COLISIONES: ¿la cancelacion de Phi tiene PAREJA UNICA?  15 de agosto de 2026. | yes |
| `collision_graph2.py` | EL MONOMIO TESTIGO: anatomia de las fibras AISLADAS.  15 de agosto de 2026. | yes |
| `commuting_square.sage` | EL CUADRADO CONMUTATIVO:  las dos factorizaciones del mismo objeto.   15 de agosto de 2026. | yes |
| `conj_crit_t2.sage` | t = 2: EL ESTRATO DE ARRIBA BASTA, Y ESO CIERRA conj:crit.  Con la prueba, y con el control que | yes |
| `construct.py` | EL CONSTRUCTOR: supervivientes de profundidad A LA CARTA, sin barrer.  13 de agosto de 2026. | yes |
| `contact_order.py` | LA PROFUNDIDAD ES ORDEN DE CONTACTO CON EL LUGAR CONCENTRICO.  13 de agosto de 2026. | yes |
| `contradiction_hunt.py` | Authors: Carles Marin, Claude (AI assistant). | -- |
| `core_conditional.sage` | ¿ES EL t-CORE EL DISCRIMINANTE?   Condicionar, que es lo que convierte una señal en criterio. | yes |
| `criterion_S.py` | EL CRITERIO SOBRE S.  14 de agosto de 2026. | yes |
| `criterion_control.py` | CONTROL INDEPENDIENTE del criterio del 14 de agosto, en Python exacto y sin Sage. | yes |
| `criterion_sage_check.sage` | CONTROL EXTERNO del criterio del 14 de agosto, por VIA INDEPENDIENTE, con cascada y checkpoint. | yes |
| `cross_formulas.sage` | Simplify the day's formulas to canonical form and cross them. | yes |
| `defect_cone.py` | LA LEY DEL DEFECTO INVARIANTE, y el cono de generadores.  13 de agosto de 2026. | yes |
| `depth.py` | EL ENUNCIADO DE PROFUNDIDAD.  Baja algun estrato de D1 - 4? | yes |
| `depth_histogram.py` | EL HISTOGRAMA DE PROFUNDIDADES SOBRE LOS 124 SUPERVIVIENTES.  13 de agosto de 2026. | yes |
| `descent.py` | EL DESCENSO Y LA FRONTERA: la ruta VERTICAL para t >= 4.  13 de agosto de 2026. | yes |
| `dim_certificate.py` | EL CERTIFICADO DE DIMENSION, y con el un barrido que llega donde el polinomio no llegaba. | yes |
| `divided_differences.py` | LA DIVISION POR Delta_t COMO OPERADOR DE DIFERENCIAS DIVIDIDAS.   16 de agosto de 2026. | yes |
| `dominant_vector.py` | EL VECTOR DOMINANTE v(T), Y SI SU MAXIMO DECIDE.  15 de agosto de 2026. | yes |
| `dump_polytopes.sage` | LOS TRES POLITOPOS DE LA LEY, calculados para dibujarlos.   15 de agosto de 2026. | -- |
| `e_equals_t.py` | LA CONTINGENCIA e = t, ARCHIVADA.  Es el recibo que faltaba. | yes |
| `even_gkrs_identity.sage` | E6: LA IDENTIDAD GKRS PAR, COMPLETA.   16 de agosto de 2026. | yes |
| `even_transversal.sage` | EL ANALOGO PAR DEL TRANSVERSAL.   16 de agosto de 2026. | yes |
| `excess4_sp_support.sage` | Excess 4: expand Psi_2 in symplectic characters and look at the SUPPORT. | yes |
| `excess_collapse.sage` | Does the value collapse, and does the collapse survive above excess 2? | yes |
| `excess_fibre_anatomy.sage` | What do the members of a large fibre at excess 4 share? | yes |
| `excess_fibre_growth.sage` | Does the fibre GROW?  The range-free form of the question. | yes |
| `excess_invariant.sage` | What IS the invariant above excess 2?  Measured against the true fibres, not guessed. | yes |
| `excess_s3_falsify.sage` | Trying to break the previous conclusion, before believing it. | yes |
| `excess_s3_symmetry.sage` | The formula is symmetric in d1, d2, d3.  Is anything else? | yes |
| `excess_t2_block.sage` | t = 2: what blocks a permutation, when the arithmetic mod t allows it? | yes |
| `excess_t2_minsize.sage` | Is "missing" at t = 2 a block, or just a bigger partition than I swept? | yes |
| `extra_by_kernel.sage` | extra_by_kernel.sage | yes |
| `extremes_of_S.py` | ============================================================================================ | yes |
| `falsify_law.py` | FALSEAR LA LEY DEL EXTREMO FUGITIVO.  13 de agosto de 2026. | yes |
| `fibres_1007.py` | Problem 10.7 -- the fibres of the evaluation invariant. | yes |
| `fibres_1007b.py` | Problem 10.7, second pass: the fibre as a LATTICE, with the paper's own (A,B) convention. | yes |
| `filter_is_periodic.sage` | EL FILTRO ES UNA FUNCION SOBRE (Z/t)^rango, NO SOBRE PARTICIONES.   16 de agosto de 2026. | yes |
| `filter_translate.py` | ¿SON LOS DOS FILTROS LA MISMA FUNCION, TRASLADADA?   16 de agosto de 2026. | yes |
| `flag.py` | LA BANDERA: beta no es una configuracion con un numero, es una CADENA.  13 de agosto de 2026. | yes |
| `folding_sage_check.sage` | ============================================================================================ | yes |
| `folding_t2.py` | ============================================================================================ | yes |
| `folding_t2_dump.py` | ============================================================================================ | yes |
| `folding_t2_refino.py` | LA PREDICCION CORREGIDA: no UN caracter simplectico, sino una COMBINACION ENTERA NO NEGATIVA | -- |
| `formulas_from_scratch.sage` | Independent check of the paper's central displayed formulas, written from the PRINTED statements | yes |
| `fs_indicator.py` | .ES EL FILTRO UN INDICADOR DE FROBENIUS-SCHUR DISFRAZADO?   17 de agosto de 2026. | yes |
| `fusion_minimal.sage` | LOS FILTROS COMO COCIENTES DE FUSION DE NIVEL MINIMO.   16 de agosto de 2026. | yes |
| `galois_sign.sage` | EL SIGNO DE GALOIS DEL FILTRO.   16 de agosto de 2026.   (vuelta 28, punto 6) | yes |
| `gamma_par.sage` | .SE CIERRA gamma_t TAMBIEN EN LA RAMA PAR?   16 de agosto de 2026.  (vuelta 30, punto 7) | yes |
| `gap_inequality.py` | LA DESIGUALDAD  g_ext < g_int,  ahora con prueba y no solo medida. | yes |
| `gate_folded_lemma.sage` | LITERATURE GATE on the folded lemma.  Is it Nadimpalli-Pattanayak-Prasad in type C? | yes |
| `gate_hall.sage` | GATE on the Hall / matching framing.  Is it a matching condition, or just a rank drop in a | yes |
| `gate_selfcomp.sage` | THE GATE THAT CAN KILL IT: is the criterion just branch (b) in disguise? | yes |
| `gate_specialization.sage` | The interior is the t = 2 problem with pairs specialized to roots of unity. | yes |
| `gate_specialization2.sage` | The specialization reading, with the two defects of gate_specialization.sage fixed. | yes |
| `gcom_branching.py` | g_com EN EL LENGUAJE DEL BRANCHING.   16 de agosto de 2026. | yes |
| `gkrs_L1.sage` | GKRS PARA (L1).   16 de agosto de 2026. | yes |
| `gkrs_centralizer.sage` | .EN QUE BASE VIVE (L1)?  EL CENTRALIZADOR DEL ELEMENTO DE TORSION.   16 de agosto de 2026. | yes |
| `global_check.sage` | DOS COMPROBACIONES GLOBALES que la mirada de conjunto exige antes de opinar de nada. | yes |
| `k_vs_m.py` | SE MUEVE K CON M?  El barrido que hasta ahora no se podia hacer.  13 de agosto de 2026. | yes |
| `kernel_residual.sage` | The second brick: the residual kernel of the specialization. | yes |
| `kernel_residual2.sage` | The residual kernel, second attempt: measured against ratios, with no denominator to get wrong. | yes |
| `kt_calibrate.sage` | kt_calibrate.sage -- DERIVE the Koike-Terada determinant instead of recalling it. | yes |
| `kumari_618.sage` | kumari_618.sage -- does the SUPER version, Theorem 6.18 of the thesis (= Theorem 4.5 of the | yes |
| `kumari_618_t15_t21.sage` | kumari_618_t15_t21.sage -- the supertableaux theorem at the next two odd composites. | yes |
| `kumari_gate.sage` | kumari_gate.sage -- two independent gates on Kumari's thesis, both run from the definitions. | yes |
| `kumari_repair.sage` | kumari_repair.sage -- does the REPAIRED theorem still do the job the original was written for? | yes |
| `kumari_rescue.sage` | kumari_rescue.sage | yes |
| `kumari_root.sage` | kumari_root.sage -- WHY does Theorem 4.4 of the withdrawn arXiv:2211.14093 fail exactly at the | yes |
| `kumari_where.sage` | kumari_where.sage -- WHERE does the negative multiplicity sit? | -- |
| `layer_condition.py` | LA CAPA: que le falta a la condicion necesaria para ser criterio.  14 de agosto de 2026. | -- |
| `lemma_V_eq_C.sage` | LEMA: si la reflexion intercambia los dos maximizadores, el valor de empate es EXACTAMENTE C. | yes |
| `lemstep_audit.py` | Authors: Carles Marin, Claude (AI assistant). | yes |
| `levi_frame.py` | EL MARCO DE LEVI, verificado a maquina. | yes |
| `link_P_separates.sage` | EL ESLABON QUE QUEDA: [Phi]_top = 0 => T_B = C - T_A.  Se reduce a "P separa", y eso se testea. | yes |
| `link_closed.sage` | LA CADENA CERRADA: Phi_t = 0 => (ii), salvo UN solo enunciado, y ya no es nuestro. | yes |
| `local_parity.sage` | SU LEMA (10):  los r_i SIN reordenar, y la anatomia de la cancelacion.  15 de agosto de 2026. | yes |
| `master_family.py` | LA FAMILIA MAESTRA: ¿es CANCEL eventualmente CONSTANTE o solo eventualmente PERIODICA? | yes |
| `mechanism_rn.sage` | The mechanism at general r, checked before it is believed. | yes |
| `middle_block.sage` | THE MIDDLE BLOCK: a PROVED sufficient condition for non-vanishing, and where it stops. | yes |
| `minimal_audit.py` | Authors: Carles Marin, Claude (AI assistant). | yes |
| `moves_table.py` | LA TABLA DE CUATRO CASOS DEL LEMA DE LOS DOS MOVIMIENTOS, testada fila por fila. | yes |
| `mu_max_law.py` | LA LEY DEL PESO SIMPLECTICO MAXIMO.  15 de agosto de 2026. | yes |
| `mumax_odd.sage` | LA LEY DEL PESO SUPERIOR EN EL IMPAR: ¿que rho lleva el desplazamiento?   16 de agosto de 2026. | yes |
| `mumax_transversal.py` | .ES EL TRANSVERSAL DE mu_max EL DE COORDENADAS MAS PEQUENAS?   16 de agosto de 2026. | yes |
| `mumax_vs_transversal.sage` | EL TOP TRANSVERSAL CONTRA LA LEY DE mu_max.   16 de agosto de 2026. | yes |
| `newt_denominator.sage` | EL VERTICE DOMINANTE DEL DENOMINADOR, POR SUMAS DE MINKOWSKI.   16 de agosto de 2026. | yes |
| `newt_zonotope.sage` | Newt(N_delta) POR SUMAS DE MINKOWSKI, ENTERO.   16 de agosto de 2026. | yes |
| `newton_vertices.py` | LOS VERTICES DEL POLIGONO DE NEWTON SON ATOMOS SUELTOS.  15 de agosto de 2026. | yes |
| `npp_minus.sage` | NPP on the det = -1 component: the count runs over the LONG roots. | yes |
| `npp_minus2.sage` | NPP on the det = -1 component, second attempt: the count runs over the FOLDED residues. | yes |
| `npp_principal_check.sage` | ¿ES NUESTRO ELEMENTO DE TORSION CONJUGADO A UN ELEMENTO PRINCIPAL?   16 de agosto de 2026. | yes |
| `npp_question81.sage` | NUESTRO FILTRO CONTRA LA QUESTION 8.1 DE NADIMPALLI-PATTANAYAK-PRASAD.   16 de agosto de 2026. | yes |
| `npp_witness.sage` | EL TESTIGO MINIMO, POR TRES RUTAS INDEPENDIENTES.   16 de agosto de 2026. | yes |
| `nu_transversal.py` | EL NUMERADOR GKRS ES UN RECUENTO DE TRANSVERSALES.   16 de agosto de 2026. | yes |
| `odd_BD.sage` | EL IMPAR EN SU GRUPO:  B_{R'} -> B_{m'} x D_r,  la propuesta 3 de su reseña.   16 de agosto de 2026. | yes |
| `odd_atoms.sage` | LOS ETA QUE SOBREVIVEN EN EL PESO SUPERIOR, LADO IMPAR.   16 de agosto de 2026. | yes |
| `odd_basis.sage` | ¿EN QUE BASE VIVE EL OBJETO IMPAR?   15 de agosto de 2026. | yes |
| `odd_companion.sage` | M6 -- EL COMPAÑERO IMPAR:  la reduccion a t=1, y un filtro con UNA PARED MENOS. | yes |
| `odd_dichotomy.sage` | LA DICOTOMIA DE PARIDAD, MEDIDA EN SUS PARAMETROS.   16 de agosto de 2026. | yes |
| `odd_extremal.sage` | LA ESTRUCTURA EXTREMAL DEL LADO IMPAR: ¿por que sale +-1?   16 de agosto de 2026. | yes |
| `odd_per_lambda.sage` | LA PRIMITIVIDAD IMPAR, LOCALIZADA POR Lambda.   16 de agosto de 2026. | yes |
| `odd_sign_formula.sage` | EL SIGNO DEL FILTRO IMPAR, EN FORMA CERRADA.   16 de agosto de 2026. | yes |
| `oddsign_convention.py` | LA CONSTANTE DEL FILTRO IMPAR, EN EL CONVENIO DEL PAPER.   16 de agosto de 2026. | yes |
| `ordertype_split.sage` | ¿LA CANCELACION OCURRE DENTRO DE CADA ORDER TYPE O ENTRE ELLOS?   15 de agosto de 2026. | yes |
| `parity2.sage` | The Laplace parity, rewritten from scratch, with the printed case as an acceptance test. | yes |
| `parity_alternates.sage` | Does the pairwise cancellation alternate with the parity of r? | yes |
| `parity_anchor_r1.sage` | *** THE CLOSING PROSE OF L3 IS WRONG AND IS LEFT STANDING.  It predicts that at r = 1 the two | yes |
| `parity_closed.sage` | The reflection sign q, in closed form, and the criterion on a big sample. | yes |
| `parity_r3_wide.sage` | Does the alternation survive thirty zeros? | yes |
| `parity_selfpaired.sage` | A closed form for the reflection sign q, and the criterion at a second t. | yes |
| `parity_terms.sage` | The term signs of the Laplace expansion, computed instead of derived. | yes |
| `peel_gcom.py` | EL PELADO CONTRA g_com -- 15 de agosto de 2026. | yes |
| `peel_identity.py` | LA IDENTIDAD DE PELADO -- y con ella el pelado deja de ser medida y pasa a ser lema. | -- |
| `peel_zero.py` | EL PELADO DE UN CERO: 14 de agosto de 2026. | yes |
| `pieri_population.sage` | LA RECURSION DE PIERI SOBRE POBLACION, Y DONDE SE ROMPE.   16 de agosto de 2026. | yes |
| `pieri_recursion.sage` | UNA RECURSION DE PIERI PARA LOS A_mu, Y QUE LE HACE AL PESO SUPERIOR.   16 de agosto de 2026. | yes |
| `pob_helper.py` | Envoltorio minimo para que un guion de Sage use las poblaciones de peel_zero sin que el | -- |
| `proof_E_even.sage` | PROOF that /E/ is even, with every step machine-checked, and the control that would make it vacuous. | yes |
| `proof_G_le_2.sage` | /G/ <= 2, DEMOSTRADO -- y con ello t impar no tiene segunda rama, por una ruta nueva. | yes |
| `quadrant_big.sage` | The same question with 30-40 zeros instead of 8.  Needs a different engine. | yes |
| `quadrant_criterion.sage` | The criterion, stated on beta alone. | yes |
| `quadrant_debug.sage` | Two of my scripts disagree.  One case, printed term by term, decides which. | yes |
| `quadrant_mechanism.sage` | *** SU CONTEO ES FALSO. NO CITAR EL 30/30 NI EL 0/635. *** | yes |
| `quadrant_necessity.sage` | The attack on my own sweep: I filtered by concentric BEFORE computing. | yes |
| `quadrant_parity.sage` | *** AVISO RETIRADO: ESTE GUION ESTABA BIEN. *** | -- |
| `quadrant_separator.sage` | t = 4, r = 2: what separates the vanishing concentric shapes from the rest? | yes |
| `quadrant_which5.sage` | In the interior, concentric is necessary and no longer sufficient.  Which concentric shapes vanish? | yes |
| `quadrant_zeros.sage` | The interior of the (t, r) map: where does s_lambda(mu_t, z_1^{+-1}, ..., z_r^{+-1}) vanish? | yes |
| `quadrant_zeros2.sage` | The same quadrant, with the control the first run was missing, and a longer range. | yes |
| `quotient_determines.py` | Authors: Carles Marin, Claude (AI assistant). | yes |
| `quotient_split.sage` | DENTRO DE UN SOLO CORE:  6 se anulan y 6 no.  ¿QUE LAS SEPARA?   15 de agosto de 2026. | yes |
| `reduction_general_t.sage` | Towards the converse: the general-t analogue of (eq:Cmu), and a defect of mine corrected. | yes |
| `reduction_to_t2.py` | LA REDUCCION A t=2, Y DOS CONTROLES QUE MATAN DOS IDEAS.  15 de agosto de 2026. | yes |
| `reflection_pairs.py` | LA REFLEXION EMPAREJA LAS CLASES DE EXCESO -- verificacion de la Proposicion probada el 14. | yes |
| `regular_in_G.py` | (R_C): EL FILTRO ES REGULARIDAD EN EL GRUPO ORIGINAL, NO EN EL DUAL.   16 de agosto de 2026. | yes |
| `rename_collisions.py` | Renombra las colisiones de notacion de note_t2 antes de integrarla en el paper. | -- |
| `root.py` | LA RAIZ: por que cancela la banda +2s -- causa medida, y mecanismo REFUTADO.  13 de agosto de 2026. | yes |
| `saturation.py` | LA SATURACION -- y la muerte de mi propio "K_sop <= 4".  13 de agosto de 2026. | yes |
| `second_stratum.py` | EL SEGUNDO ESTRATO.  Arriba esta probado; el de abajo NO es separable (H7 murio en Gbot_anatomy). | yes |
| `second_stratum_refine.py` | CIERRE DE LOS DOS CABOS QUE ABRIO second_stratum.py. | yes |
| `second_vanishes.py` | SE ANULA DE VERDAD EL SEGUNDO ESTRATO?  La pregunta que decide el alcance de la seccion 7. | yes |
| `selection_rule.sage` | UNA REGLA DE SELECCION QUE SALE DE CRUZAR DOS FORMULAS.   16 de agosto de 2026. | yes |
| `sieving_1008.sage` | sieving_1008.sage | yes |
| `sigma_involution.py` | LA INVOLUCION sigma -- el mecanismo de la direccion SUFICIENTE del criterio, para todo t y r. | yes |
| `sign_layer.py` | LA CAPA DE SIGNOS DEL FILTRO: .que modulo ve?   16 de agosto de 2026. | yes |
| `sign_lemma.py` | EL LEMA DE SIGNO w(sigma P) = -w(P), y el ataque que lo pudo tumbar. | yes |
| `sign_ratio.sage` | EL ESLABON, REDUCIDO A UN SIGNO: w_B = -w_A  <=>  K = V - K.  Y aqui esta su forma cerrada. | yes |
| `sixteen_betas.py` | LAS 16 FORMAS DE  t=6, r=2  CON  Delta != 0.   15 de agosto de 2026. | yes |
| `socrates_audit.py` | EL PAPER II, AUDITADO -- Y LA HERRAMIENTA, AFINADA.   16 de agosto de 2026. | yes |
| `sp_expansion.sage` | LA EXPANSION SIMPLECTICA  Phi_{t,r} = sum_mu A_mu · sp_mu.  15 de agosto de 2026. | -- |
| `sp_law.sage` | LA LEY DEL PESO MAXIMO: ¿se lee mu_max del defecto de simetria de g_com?  15 de agosto de 2026. | yes |
| `sp_law_fit.py` | ¿DE QUE DEPENDE mu_max?  El ajuste, con su control de vacuidad.  15 de agosto de 2026. | -- |
| `step_law.py` | LA LEY DE PASO -- el Lema 4.4 del paper, para clases de tamano cualquiera. | yes |
| `stratum_sign.sage` | THE TOP STRATUM HAS ITS OWN (i)+(ii): closure is the concentricity, and the missing half is a SIGN. | yes |
| `support_ladder.py` | LA ESCALERA CON SOPORTE, EN TODAS LAS (t, r).  13 de agosto de 2026. | yes |
| `survivors_are_weyl.sage` | LOS SUPERVIVIENTES SON EL GRUPO DE WEYL, Y tau ES SU CARACTER SIGNO.   16 de agosto de 2026. | yes |
| `survivors_wide.py` | RE-MEDIDA DEL "e = t" DE LOS 12 SUPERVIVIENTES, CON M GRANDE.  13 de agosto de 2026. | yes |
| `thm31_audit.py` | Authors: Carles Marin, Claude (AI assistant). | yes |
| `threshold_alcove.sage` | .ES EL EMPATE DEL SENUELO UN UMBRAL DE NIVEL?   prob:threshold, 17 de agosto de 2026. | yes |
| `top_stratum_blind.py` | EL ESTRATO DE ARRIBA NO PUEDE DECIDIR g_com -- y aqui esta el testigo que lo prueba. | yes |
| `topdeg_gate.sage` | THE TOP-DEGREE COMPONENT, which factorises into two GL(r) alternants -- and the anatomy of the | yes |
| `topstratum.sage` | WHY THE TOP-DEGREE PART MISSES EXACTLY 20 SHAPES: the stratum is a smaller copy of the problem. | yes |
| `topterm_gate.sage` | THE CONVERSE, attacked by its extremal term.  Does the leading monomial survive? | yes |
| `torsion_filter.sage` | EL FILTRO DE TORSION  tau_t(eta) = sp_eta(xi, xi^2, ..., xi^m).   15 de agosto de 2026. | yes |
| `twelve_forms.sage` | LAS DOCE.  El residuo del problema abierto, y es pequeño.   15 de agosto de 2026. | yes |
| `two_strata_audit.py` | AUDITORIA DE MI PROPIA REFUTACION.  12 de agosto de 2026 (noche). | yes |
| `two_strata_depth.py` | LA PREGUNTA QUE ABRE S4.  12 de agosto de 2026 (noche). | yes |
| `two_strata_fail.py` | LAS OCHO QUE ROMPEN EL ENUNCIADO DE DOS ESTRATOS.  Verificacion dura, 12 de agosto de 2026 (noche, tras el com | yes |
| `two_strata_wide.py` | EL BARRIDO ANCHO, y las dos lecturas mias que mata.  12 de agosto de 2026 (noche). | yes |
| `unimodularidad_barrido.py` | BARRIDO ANCHO DE LA UNIMODULARIDAD.   16 de agosto de 2026. | yes |
| `unit_invariance.py` | EL LUGAR DEL FILTRO ES GALOIS-INVARIANTE.   16 de agosto de 2026. | yes |
| `v28_witnesses.sage` | LOS TESTIGOS DE LA VUELTA 28, verbatim.   16 de agosto de 2026. | yes |
| `v29_two_lemmas.sage` | LAS DOS PROPUESTAS DE LA VUELTA 29, verificadas antes de escribirlas.   16 de agosto de 2026. | yes |
| `v2_formulas.py` | ============================================================================================ | yes |
| `vanishing_criterion.sage` | M1 -- POR QUE SE ANULA:  ¿por SOPORTE o por CANCELACION?   15 de agosto de 2026. | yes |
| `verify_rename.py` | Verifica el renombrado en el FUENTE: aplica el mapeo inverso al .tex actual y lo compara con la | -- |
| `wall_table.sage` | LA TABLA DE LOS 16:  que pared mata al candidato, y quien sobrevive.   15 de agosto de 2026. | yes |
| `weyl_population.sage` | ¿ES EL CONJUNTO DE SUPERVIVIENTES UN TORSOR DEL GRUPO DE WEYL?   16 de agosto de 2026. | yes |
| `weyl_structure.sage` | ¿ES EL CONJUNTO DE SUPERVIVIENTES UN TORSOR DEL GRUPO DE WEYL?   16 de agosto de 2026. | yes |
| `why_count_dies.sage` | Why the root count dies at r = 2 and the covering condition does not. | yes |
| `why_delta.py` | POR QUE Delta != 0.  De medir el fenomeno a buscarle el mecanismo. | yes |
| `widen_G.sage` | WIDENING THE ONE SIZE THAT MATTERS: is /G/ <= 2 real, or an artefact of the range? | yes |
| `witness_family.py` | LA RECTA K = W/2 - 3: ¿en que escalera?  Y QUIENES son los testigos.  13 de agosto de 2026. | yes |
| `yacobi_factor.sage` | ¿ES NUESTRO Sp_2 UNO DE LOS FACTORES SL_2 DE YACOBI, O NO?   15 de agosto de 2026. | yes |
| `yacobi_parity.sage` | LA RUTA DE YACOBI, EJECUTADA:  A_mu por paridad de los r_i, sin pasar por los eta. | yes |
| `zeros_propagate.sage` | EL LOCUS DE ANULACION ES UN IDEAL, Y ESO LE PROHIBE COSAS.   16 de agosto de 2026. | yes |

## Figures (15)

| script | what it does | archived output |
|---|---|---|
| `fig_alcove.py` | (sin titulo) | yes |
| `fig_collapse.py` | (sin titulo) | yes |
| `fig_cone3d.py` | (sin titulo) | yes |
| `fig_determinant.py` | LA FIBRA Y EL DETERMINANTE SON LA MISMA SUMA.   16 de agosto de 2026. | yes |
| `fig_division.py` | LA UNICA DIFICULTAD QUE QUEDA EN (L1): LA DIVISION POR Delta_t.   16 de agosto de 2026. | yes |
| `fig_filter.py` | (sin titulo) | yes |
| `fig_galois.py` | (sin titulo) | yes |
| `fig_law3d.py` | (sin titulo) | yes |
| `fig_problems.py` | ============================================================================================ | yes |
| `fig_residue.py` | (sin titulo) | yes |
| `fig_thread.py` | (sin titulo) | yes |
| `fig_transversal.py` | EL NUMERADOR GKRS, DIBUJADO: un transversal de las clases plegadas.   16 de agosto de 2026. | yes |
| `fig_walls3d.py` | (sin titulo) | yes |
| `figlang.py` | (sin titulo) | yes |
| `figs_es.py` | (sin titulo) | yes |

## Audits (27)

| script | what it does | archived output |
|---|---|---|
| `_abslen.py` | Longitud del abstract, contada COMO LA CUENTA arXiv: caracteres LITERALES del texto que se pega | yes |
| `_ambiguity.py` | AUDITORIA DE AMBIGUEDAD.  .Usa el paper la misma palabra para dos cosas distintas? | yes |
| `_attraudit.py` | AUDITORIA DE ATRIBUCION.  .Nos estamos colgando medallas de otros? | yes |
| `_bibaudit.py` | Audita la bibliografia del paper II: toda entrada citada, toda cita con entrada, y ningun nombre | yes |
| `_bibparity.py` | .CITAN LAS DOS EDICIONES LO MISMO, EN LOS MISMOS SITIOS?   18 de agosto de 2026. | yes |
| `_citeaudit.py` | Auditoria de citas: entradas de la bibliografia que nunca se citan, y citas sin entrada. | yes |
| `_claimaudit.py` | Inventario de afirmaciones con marcador de estado, para el repaso de principio a fin. | yes |
| `_ctrlchars.py` | NI UN CARACTER DE CONTROL EN EL FUENTE, Y NADA QUE LA EXTRACCION SE COMA EN UN ENUNCIADO. | yes |
| `_dumpcheck.py` | CADA NUMERO DEL PAPER, CONTRA SU VOLCADO. | yes |
| `_figaudit.py` | Auditoria de figuras: cuales tienen \label y cuales se CITAN en el texto, y el hueco de pagina. | yes |
| `_figs_es_audit.py` | (sin titulo) | yes |
| `_fixbib.py` | Rehace la bibliografia del paper II y engancha los \cite que faltaban. | yes |
| `_fixclaims.py` | Corrige tres sobre-enunciados del manuscrito y actualiza dos cifras.  Cada uno es un error mio, | yes |
| `_formaudit.py` | AUDITORIA TIPOGRAFICA DE LAS FORMULAS, ANTES DE PUBLICAR. | yes |
| `_icollide.py` | LA COLISION DE LA UNIDAD IMAGINARIA. | yes |
| `_insert.py` | Inserta las dos secciones que faltaban para que el paper sea hermano del companion. | yes |
| `_insert2.py` | Inserta las tres secciones nuevas, las filas de verificacion que faltan y las referencias. | yes |
| `_introaudit.py` | Auditoria de la INTRODUCCION contra el cuerpo del paper. | yes |
| `_moveregular.py` | Mueve la subseccion del Lema de regularidad de la seccion del caso impar a la del filtro, | yes |
| `_numaudit.py` | AUDITORIA DE NUMEROS ANTES DE PUBLICAR. | yes |
| `_tablecross.py` | CRUCE DE LAS DOS TABLAS.  La tabla de verificacion dice el estado con un marcador de color; la de | yes |
| `_taurename.py` | Renombra tau_t -> tau^C_t en el TEXTO PAR, dejando el impar con su tau^B_t. | yes |
| `_transcheck.py` | (sin titulo) | yes |
| `_verifaudit.py` | Barre la tabla de verificacion fila a fila: busca filas duplicadas, filas sin estado, y filas cuyo | yes |
| `_vuelta26.py` | Auditoria de los puntos de su vuelta 26 contra el .tex, uno por uno. | yes |
| `_vuelta27.py` | Auditoria de los cinco puntos de su vuelta 27, contra el .tex. | yes |
| `_whatsleft.py` | Inventario de lo que queda: toda fila de Attribution que NO diga "proved", y los problemas abiertos. | yes |
