Ancillary material for

    "Factorization and vanishing of Schur polynomials twisted by roots of unity and reciprocal pairs"
    Carles Marin

Every numerical claim in the paper regenerates from these scripts, and the saved output of each run
is included in outputs/ so that any count in the verification table (Section 9) can be located
without running anything.

The material comes in two groups.

    GROUP 1  (Sections 3-7)   plain Python, needs only `mpmath`
                              (`numpy` and `matplotlib` for the figure scripts).
    GROUP 2  (Section 8)      Sage. These compute in exact Laurent rings over several variables and
                              solve linear systems over Q, which is what Sage is for.

    pip install mpmath numpy matplotlib          # group 1
    rm -f script.sage.py && sage script.sage     # group 2  -- see the warning below

Run them from this directory.

WARNING, and it cost us a wrong answer once. `sage script.sage` preparses the source into
`script.sage.py` and will happily re-run that cached copy: an edited .sage can produce, silently,
the output of the code it replaced -- identical byte for byte, so nothing looks wrong. Delete the
.sage.py first, every time. None of them ship with this paper for the same reason.

Section and statement numbers below refer to the paper these files accompany. What is collected here
is the scripts and their saved output, not the paper itself. The published rank-one case at t=2 is
Part IV, doi:10.5281/zenodo.21463000.

--------------------------------------------------------------------------------------------------
GROUP 1 -- THE EVALUATION AND ITS SHARPNESS  (plain Python)
--------------------------------------------------------------------------------------------------

theorem_full.py         Theorem 3.1 with its exact sign, and Lemmas 4.1-4.5 separately.
                        The bialternant is evaluated from scratch and compared against the closed
                        form; the sign is compared too, not only the magnitude.
                        -> 959 exact matches, 476 zeros, 0 failures; L4 724/724.   (Section 9)

law_control.py          An independent second implementation of the same check, written from the
                        statement of the theorem alone rather than from theorem_full.py. This is
                        the "two implementations along different code paths" of Section 9. It also
                        supplies the shared determinant and partition routines.
                        -> 959 / 476 / 0.

thm_main_independent.py A THIRD implementation of Theorem 3.1, written from the printed statement
                        alone: it reads the beta set, builds the bialternant at 60 digits and
                        compares against the closed form, sign included. It shares no code with
                        theorem_full.py or law_control.py. Its decoy forces eps=+1 and must fail.
                        -> 749 shapes over t=2..6, 0 failures; decoy fails 35 of 133. (Section 9)

invariant_separates.py  Proposition 3.10 the same way, in both directions: same datum => same
                        value, and same value => same datum. Its decoy drops the sign from the
                        datum and must then collide.
                        -> 676 non-vanishing shapes, 0 collisions; decoy 104.      (Section 9)

falsify.py              The controls (D3) and (D4) of Section 7: the same test applied to the coset
                        alphabet and to a free (non-reciprocal) pair.
                        -> orbit 600/0, coset 383/217, free pair 200/400.          (Sections 7, 9)

d_from_quotient.py      Proposition 3.4: the three arguments read off the t-quotient.
                        -> 2970/2970 in the two-class profile.                     (Section 3.1)

sign_ayyer_idiom.py     Proposition 3.11: the short form of the sign in the notation of [AK25].
                        Prints a CONTROL first -- the same comparison with the two blocks in
                        arbitrary order, which fails 592/904 and is meant to -- then the ordered
                        comparison, then the cell test. The control is what makes the ordering
                        hypothesis load-bearing.
                        -> control 592/904; ordered 1496/1496; 112 cells, 0 mixed.  (Section 3)

AUDIT_FORMULAS.py       The displayed FORMULAS rather than the counts. Every other script here checks
                        a claim about numbers; this one evaluates both sides of each displayed
                        identity at numeric points, straight from the definitions -- s_lambda is
                        always the bialternant and h_j always comes off the generating function --
                        so a wrong sign or a stray constant cannot survive by being shared with the
                        machinery that produced the tables. Covers Theorem 2.1 with its sign, the
                        interval reading, Lemmas 4.2 and 4.5, the splitting of the alphabet, the
                        complementation identity, d'Ocagne, and the scalar caveat of Section 10.
                        -> 3297 evaluations over 9 formulas, 0 failures.  (all sections)

alternation_proof.py    Proposition 6.4, one line per step of the PROOF rather than its conclusion:
                        the beta-form M_ij = [beta_i >= c_j][beta_i = c_j mod 2] of Jacobi-Trudi;
                        that sigma != 0 exactly when each parity class carries as many columns as
                        rows with segment lengths 1..size, which is why |sigma| <= 1; that k odd
                        admits no two consecutive nonzeros; and that consecutive nonzeros exchange
                        the two segment lengths, which is what leaves both block determinants alone
                        so that only the shuffle sign flips.
                        -> 4823 / 4823 / 1341 / 217, 0 failures.  (Section 6)

involution_runs.py      Proposition 6.4, second half: the cancellation across nu, at t = 2. For fixed
                        k = nu1 - nu2 the admissible nu are the single family (m+k, m), so the whole
                        question is the word sigma_m(k) = s_{lambda/(m+k,m)}(1,-1), computed by
                        Jacobi-Trudi with h_j(1,-1) = [j even]. Reports that |sigma| <= 1 (a set,
                        not a multiset), that the sign alternates inside every maximal run of
                        consecutive nonzeros, and -- the control -- that the survivors, one per
                        odd-length run, rebuild the closed form coefficient by coefficient.
                        -> max 1; 0 of 1406 runs break it; 333/333.  (Section 6)

concentric_locus.py     Proposition 3.5: the concentric branch d3 = 0 in quotient coordinates. Checks
                        that it is empty for t ODD -- with the even-t count printed beside it, 43,
                        because a script that finds nothing anywhere would report the same zero --
                        that for t EVEN it is exactly rB - rA = t/2 with |lam^(rA)| = |lam^(rB)| + 1,
                        both directions, and that d3 is NOT a function of the residue profile, so
                        the locus is an arrangement in the quotient sizes and not in the core
                        lattice.
                        -> odd 0, even 43; 1331/1331; 116 of 156 profiles split.  (Section 3)

extra_sign.py           The sign of Theorem 5.2, on both of its families, so that the theorem is an
                        equality and not an equality up to sign.  On the t-cores that is (4) itself.
                        On the extra family it is the inversion count (11): with t = 2m and
                        lambda_2 = m + j the residue word of beta is (j, m+j, 2m-1, ..., 0) and
                        inv(w) = C(2m,2) + j + (m+j) = 2m^2 + 2j, which is even.  Recomputed here
                        from beta, with eps by the LONG formula (7) and, as a route check, by the
                        SHORT one (8) that the proof uses.  Three controls: eps must not be
                        identically +1 (it is -1 on 80 of 3038 two-row shapes in range), the
                        near-miss 2m^2 + j must be refuted (21 of 21 with j > 0), and no shape
                        outside the two families may satisfy the equality (0 of 2835).
                        -> 28 extras + 219 cores, all +1, 0 failures.  (Section 7)

invariant_minimal.py    Proposition 3.10 and identity (10).  First, that the evaluation invariant is
                        MINIMAL and not merely complete: distinct (multiset, sign) give distinct
                        values.  The map tested for injectivity is (multiset, eps) -> eps * prod_i
                        (u^d_i - u^-d_i) with u = z^(1/2), exact integer Laurent arithmetic, over
                        every triple with entries up to a bound rather than only the realised ones.
                        The control grades the same population by a lossy invariant, the SUM instead
                        of the multiset, which collides 75284 times of 75640.  Second, the shift law:
                        beta -> beta + 1 is lambda -> lambda + (1^N), the triple does not move and
                        the value picks up det(A) = (-1)^(t+1), so eps carries the whole reversal.
                        Its control drops the exponent, and fails on every shape at even t.
                        -> 75640 invariants, 0 collisions; 826/826 on the shift law.  (Section 3)

fig_plane.py            Draws fig_plane, the two-parameter family (t, r) that the introduction
                        promises to map and previously described only in prose.  Nothing is placed
                        by hand: each cell's status is computed from the statements themselves, one
                        rule per line of the paper -- r = 1 gives the value (Theorem 3.1), odd t the
                        criterion with no external input (Corollary 8.22), t = 2 the criterion with
                        one (Theorem 8.6), and the rest one implication (Corollary 8.34) with the
                        converse conjectural.  Its own control: the cell counts are printed and must
                        sum to the grid, so a miscoloured or duplicated cell cannot pass silently.
                        (Section 1)

fig_involution2.py      Draws fig_pairing, the mechanism of Corollary 8.34: S on a line with the arcs
                        v <-> C-v and the axis at C/2, a transversal and its reflection with both
                        signs COMPUTED, and the pairing counted.  The second panel is the control and
                        is FOUND BY SEARCH, not chosen: the smallest beta with C-S = S but no
                        increment equal to C.  There the fixed class has odd size, two transversals
                        are their own reflection, and the sum cannot cancel -- which is exactly what
                        the second clause of the hypothesis buys.  Without that panel the figure
                        would teach that symmetry alone suffices, and it does not.  (Section 8)

fig_intro.py            Draws fig_thread (the paper as a chain of seven questions), fig_alphabet (the
                        alphabet in the complex plane) and fig_beta (the beta/residue/profile
                        dictionary carried out on one shape).  Every coordinate is computed from the
                        definitions rather than placed by hand.  (Sections 1-2)

fig_proof.py            Draws fig_laplace (the architecture of the proof of Theorem 3.1, with the
                        three residue profiles and the surviving column-pair counts 0, 4, 3 counted
                        rather than asserted) and fig_map (what the criterion of Section 8 rests on,
                        with statements named rather than numbered so the figure cannot drift when
                        the numbering does).  (Sections 4 and 10)

fig_v2.py               Draws fig_increments (Lemma 8.15 carried out on one shape, tie and prefix
                        form visible), fig_reflection (Proposition 8.28 and Corollary 8.30 on two
                        shapes, one with |G| = 2 and the reflection and one with |G| = 2 and
                        without) and fig_virtual (the sp_mu expansion on two shapes differing in a
                        single box, genuine against properly virtual).  Also computes Table 3.
                        (Section 8)

extra_locus_kernel.py   Problem 10.6.  The extra locus of the independence criterion against the
                        KERNEL of the specialization, over three free parts of size two: the
                        reciprocal pair, the zeta_2-orbit (z,-z), and a genuinely free pair.  The
                        third is Remark 5.3's control and must return the t-cores and nothing else.
                        The reciprocal pair turns out to be the collapsing specialization of
                        SMALLEST kernel, which is why its extra locus is one family and (z,-z)'s is
                        28 solutions at t=2 and 186 at t=8.

fibre_lattice.py        Proposition 3.6, Corollary 3.7 and Remark 3.8.  The two-class stratum as a
                        lattice, the generating function of a fibre, and the fact that stopped the
                        remark from being wrong: the sign sees the t-2 invisible parts from t = 3 on.
                        Controls K1 (wrong denominator) and K2 (setup()'s tie-break instead of the
                        theorem's) must fail, and do.  Output in outputs/fibre_lattice.txt.

sign_one_object.py      Remark 6.5.  The sorting sign of Section 2 -- the one Proposition 3.11 puts
                        inside eps_lambda -- is the sign that makes the runs of Proposition 6.4
                        alternate.  Checks s_{lam/mu}(1,-1) = sgn(sigma_lam) sgn(sigma_mu) over
                        18037 skew pairs, then that parts (ii) and (iii) are one statement about
                        sigma alone.  Four controls, each of which must be refuted.

sign_proof_check.py     Proposition 3.11 again, but its PROOF rather than its conclusion, one step per
                        line so that a disagreement can be localised: sgn(sigma) = (-1)^inv(w); the
                        parity count for the inversions at a letter, tested on random words where it
                        can fail; the resulting formula for inv(w) - inv(b_S) with the blocks
                        ordered; and floor(t/2) + t + C(t+3,2) = 1 mod 2, which closes it. Its
                        control repeats the third step with the ordering dropped and fails 904 of
                        1496 -- that is the single point at which the hypothesis enters.
                        -> 1529 / 4000 / 1496 / 198, 0 failures; control 592/904.  (Section 3)

single_char.py          Lemma 5.1: Phi is a single sl_2 character iff the interval triple contains
                        t twice.
                        -> 2113/2113.                                              (Section 5)

ak53_consistency.py     The comparison with [AK25, Theorem 5.3] on the reciprocal locus, in both
                        directions, and the t-core routine used throughout. The core routine is
                        validated against the example in Remark 5.4 of that paper before use.
                        -> 119 t-cores + 21 extra, no others.                      (Section 5)

extra_locus.py          Theorem 5.2: the extra family, enumerated and matched against the closed
                        description.
                        -> 21 found, 21 predicted, sets identical.                 (Section 5)

extra_structure.py      The core and quotient of the extra family, and the t-quotient routine.
                        -> core (t/2-1+j, j), quotient a single (2) in slot j.     (Section 5)

enumeration.py          Corollary 6.1: the weighted tableau sum, built from Gelfand-Tsetlin chains,
                        against the closed form.
                        -> 14/14.                                                  (Section 6)

rect_degeneracy.py      The degeneracy behind Proposition 6.3: for a rectangle two of the three
                        arguments coincide. This checks the d-triple only.
                        -> c even 18/18, c odd 2/18.                               (Section 6)

rect_boxes.py           Proposition 6.3 read literally over t=2, r<=4, c<=60 (240 cases), and the
                        four box families of the Example (16 cases). The signed count is taken as
                        the z->1 limit of the closed form. Note that case (iv) of the proposition
                        must exclude r=4: at full height the d-triple is (2,2,2) for every c and the
                        value is (-1)^c, which is why r=4 is stated separately.
                        -> 240/240 and boxes 16/16.                                (Section 6)

core_vs_criterion.py    Remark 8.8, the two halves of it. First that the residue profile determines
                        the t-core at fixed N, so that branch (a) IS a core condition -- and which
                        core it is. Then the half that matters: that the core does NOT decide the
                        vanishing, since branch (b) leaves the profile untouched. Exact rational
                        arithmetic; a shape counts as vanishing only if it vanishes at three generic
                        points. Controls in the same run: branch (a) <=> some n_i = 0, and
                        (a) or (b) <=> Psi_r = 0.
                        -> 260 profiles at t=2..5, 0 profile->core clashes; 1338 shapes at r=1,2,3,
                           0 violations of either control, and 5 core classes carrying both
                           behaviours -- the empty core among them.                  (Section 8)

check_nonstd_bound.py   Lemma 8.10 against the archived outputs: every non-standard label printed by
                        the group-2 scripts satisfies |nu| >= 2r+3, and every residue shape satisfies
                        |lambda| >= 2r+3.
                        -> 0 violations.                                           (Section 8.4)

sieve_counts.py         The endpoint sieve quantified, from the same exact data the Section-8 figures
                        use. Runs the nesting control first (a shape vanishing identically must also
                        vanish at the endpoint).
                        -> control 0; flags 143 / 22 / 9 shapes at r=1,2,3, of which 0 / 4 / 3 are
                           spurious, i.e. 0% / 18% / 33%.                          (Section 8)

check_layout.py         Not a mathematical check. Reports pages carrying a large blank band inside
                        the text area, or ending early with a blank tail -- neither of which a
                        clean LaTeX run detects. It takes several files at once, because the one
                        defect it caught was in an edition nobody had pointed it at.
                        Usage: python check_layout.py orbit_pair.pdf [more.pdf ...]
                               python check_layout.py --threshold 0.2 orbit_pair.pdf

check_xref.py           Not a mathematical check either. Rebuilds from the source which environment
                        every \label sits in, then checks each reference against the word that
                        introduces it -- "Theorem \ref{lem:...}" compiles and is wrong -- and lists
                        the statements nothing ever refers to.
                        Usage: python check_xref.py orbit_pair.tex orbit_pair_es.tex
                        -> 0 wrong-kind or dangling references in both editions.

check_refs.py           Not a mathematical check either. Two things LaTeX cannot tell you: whether a
                        cross-reference points at the right KIND of object ("Theorem \ref{lem:...}"
                        compiles and is wrong), and whether any bibliography entry is never cited.
                        Usage: python check_refs.py orbit_pair.tex orbit_pair_es.tex

minimality.py           Why the alphabet of this paper is the SMALLEST one with anything to say. The
                        only extension smaller than a reciprocal pair is mu_t u {-1} with t odd, and
                        it is frozen: over t = 3,5,7 and |lambda| <= 14 the value takes only 0 and
                        +-1, with nothing left to depend on.
                        -> 1092 shapes, values outside {0,+1,-1}: 0.        (Introduction)

fixed_set.py            The fixed set of the involution of [KumT], checked against the thesis rather
                        than assumed: over skew shapes with at most 4 rows and |lambda| <= 10, the
                        surviving fillings carry one sign and the fixed set is a single filling.
                        -> 1134 skew shapes; |sigma| != #coverable: 0; #coverable != 1: 0.
                                                                            (Section 6)

check_parity.py         Not a mathematical check. Compares the two editions statement by statement so
                        that a change made in one and forgotten in the other cannot survive a build.
                        Companions: check_math.py (the displayed formulas), check_prose.py (the
                        paragraph counts, which catches a dropped paragraph that check_parity does
                        not), check_numbers.py (every number against an archived stdout),
                        check_abstract.py (the short abstract against the arXiv character limit),
                        check_floats.py (every figure against the page that refers to it, and the
                        figures nothing refers to at all).

fig_runs.py             Proposition 6.4, drawn: sigma_m(k) over m and k for two shapes, one with
                        every k even and one with every k odd, so that the three parts of the
                        proposition are visible as bar height, gaps and colour alternation.

fig_fibres.py           The fibres of the evaluation invariant -- see the entry above under the
                        figures.

figs_es.py              Draws every figure again with its internal labels in Spanish, for the
figlang.py              Spanish edition of this paper, WITHOUT duplicating the drawing scripts.
                        It intercepts matplotlib's text layer -- and in particular
                        `Text.set_text`, through which every label passes, so that a legend or a
                        hand-built label cannot slip past -- looks each string up in the
                        dictionary `figlang.py`, and redirects the output to `*_es.pdf`. The
                        scripts above run unchanged, so all of their own controls still run.
                        A string with words that is not in the dictionary is listed at the end of
                        the run rather than drawn silently in English.
                        Usage: python figs_es.py --recoge   (collect the strings)
                               python figs_es.py            (draw them)
                        The Spanish edition and its figures are not part of the arXiv submission;
                        the two files are here because they belong beside the scripts they wrap.

--------------------------------------------------------------------------------------------------
GROUP 1b -- THE OPEN PROBLEMS, ATTACKED  (Sage)
--------------------------------------------------------------------------------------------------

These are not part of the proof of any theorem. Each one attacks a statement the paper leaves open,
and two of them killed a conjecture we had believed.

coupling_rank.sage      Proposition 10.3 as first stated, tested where it had never been run: rank one
                        with B LARGER than a single reciprocal pair. A first attempt counted
                        irreducible factors, which is wrong -- chi_k is itself reducible.

coupling_rank_tight.sage  The same, done properly. "A product of characters" is tested by the roots:
                        a product of sl2 characters has every root on the unit circle, which is
                        decidable numerically and cannot be an artefact of factoring over the wrong
                        field. The paper's own theorem is the control.
                        -> control 0 off-circle; z,z^2: 17; z,z^3: 18. Rank one is
                           therefore NOT the condition.                     (Proposition 10.3)

iso_first_row.sage      The isolating-witness mechanism at r = 2, on the shapes where the certificate
                        was expected to be easiest.

iso_counterexample.sage The shape that kills it: lambda = (5,5,2,1) has no isolating mu, because at
                        ell(lambda) = 4 and N = 6 the isolation lemma leaves only one-row mu, and a
                        shape that wide contains every associate.           (Proposition 8.11)

_prove_W_ext.sage       The isolating witness at r = 2, extended to |lambda| <= 14.
                        -> 227 proved witnesses, 11 residue; the true minimum of ell(nu) over
                           non-standard nu with o_nu(A) != 0 is 5, against the 4 the bound gives.

_prove_W_r3.sage        The same at r = 3, which is the control on the explanation: isolation needs
                        ell(mu) <= 7 - ell(lambda) there, so the residue should need a wider lambda
                        than the range provides, and it does.
                        -> |lambda| <= 13: 149 proved witnesses, 0 residue.  (Section 8)

sp_at_A.sage            Problem 10.1, the symplectic column: what is sp_nu(W,1,-1)? Adjoining a
                        letter is the ring map p_k -> p_k + c^k, so the question is an identity in
                        Lambda and needs no length hypothesis.

sp_branching.sage       The answer, via branching: 247 coefficients over 74 distinct skew shapes with
                        0 clashes, values 0 and +-1, and the orthogonal identity of Section 8 still
                        holding as the control (45 hold, 0 fail).           (Problem 10.1)

so_at_A.sage            The odd-orthogonal column of the same problem. Recorded as it stands: Sage
                        carries no so basis, and doubling is not adjunction, which is why that
                        column resists where the symplectic one gave way.

extra_locus_types.sage  Problem 10.6: does Theorem 5.2 have an analogue for sp and o? Both sides of
                        the criterion are specializations of power sums -- the orbit sends p_k to
                        p_k + t when t | k, the reciprocal pair to p_k + z^k + z^-k -- so the locus
                        { lambda : LHS = +- RHS } is computable in Lambda. Type s is the control and
                        must return the t-cores plus exactly the family of Theorem 5.2. Two unrelated
                        evaluation points, since one alone would admit a coincidence.

sec8_derivation.sage    See Group 2.

--------------------------------------------------------------------------------------------------
GROUP 1c -- THE CRITERION OF SECTION 8, AND ITS TWO INDEPENDENT CONTROLS
--------------------------------------------------------------------------------------------------

criterion_base.py       The routines the five files below share: the Laplace expansion over the
                        transversals (perm_sign, setup, deg_of, all_transversals) and the degree
                        filtration (split_sign, alt, stratify, stratum, measure). They were
                        extracted from two larger scripts belonging to a separate line of work not
                        yet published; the bodies are the originals, cut out and checked character
                        by character, not rewritten.

step_law.py             Lemma 8.14, the column move at any class size: w(g')/w(g) = (-1)^{1+B+M}.
                        Its decoy drops the count M and must fail.
                        -> 285600/285600; the one-count decoy fails 149652 of 285600.

sign_lemma.py           The closed form of the sign behind Theorem 8.35, w(sigma P) = -w(P), with
                        the attack that could have brought it down: the ratio is -1 exactly when
                        eta = 2, and the criterion has to forbid eta = 1.
                        -> 438/438 over thirteen configurations; eta = 1 occurs 0 times.

criterion_S.py          The criterion itself, swept over the WHOLE population rather than the
                        restricted one: Phi = 0 iff a class is missing, or C - S = S. Its decoy
                        asks for the symmetry of the whole beta set instead of S.
                        -> 0 failures; 0 disagreements against Theorem 3.1 at r = 1; the decoy
                           fails in 24 of 24 configurations.

criterion_control.py    The same question by a route sharing no line with the above: Jacobi-Trudi
                        in exact rational arithmetic, no Laplace, no strata, no greedy. It exists
                        because a fault in the one instrument would turn 24 clean configurations
                        into 24 agreements with itself.
                        -> 43010 shapes; VEREDICTO: CONFIRMA.

criterion_sage_check.sage   And a third time, in Sage, through the generating function of the h's
                        over the whole alphabet. Caution: the header of an earlier version said
                        "power sums"; what runs is Jacobi-Trudi, and the label had been left
                        behind -- see the docstring of value().
                        -> 12937 shapes, 14 configurations; VEREDICTO: CONFIRMA; the decoy fails
                           in 14 of 14; the sieve makes 0 false discards.

--------------------------------------------------------------------------------------------------
GROUP 2 -- THE ZERO LOCUS FOR EVERY r  (Sage)
--------------------------------------------------------------------------------------------------

AUDIT_ALL.sage          Rebuilds every load-bearing claim of Section 8 from the definitions in a
                        single pass, independently of the scripts that first produced them.
                        -> 3414 checks, 0 failures. Includes branch (a) 16, branch (b) 24, the
                           even-width control 56, complementation 474, d'Ocagne 496 and 39.

selfcomp_law.sage       Theorem 8.1 in both directions: for r=1,2,3 every shape in range is tested
                        for exact vanishing and against the criterion, with false negatives and
                        false positives counted separately, plus the even-width control.
                        -> scanned 2157 + 4542 + 3262 = 9961 shapes, 242 + 55 + 21 = 318 vanishers,
                           0 false negatives and 0 false positives at every r.

close_X_r1.sage         Theorem 8.3, the converse at one pair, as a complete finite case analysis
                        over the parity patterns of the beta set.
                        -> the ten buckets sum to 14950 beta sets; 0 violations.
                        CAUTION: this script also runs an exploratory sub-check ("the 4x4 determinant
                        IS the signed sum of single S's"), which reports about 4180 mismatches out of
                        6655. That sub-claim is NOT used in the paper and its failure is not a
                        failure of anything stated there. Read only the case-analysis block.

associates_witness.sage Theorem 8.4, the converse inside Littlewood's range: an explicit associate
                        witness is constructed for every shape that is not an odd rectangle of full
                        height, and the exceptions are checked to be exactly those.
                        -> 76 + 144 + 152 = 372 witnesses, 0 failures.

prove_W.sage            Proposition 8.11, the isolating-mu statement, at r=2 and r=3, and the count of
                        non-standard labels that survive on this alphabet.
                        -> 55 + 25 = 80 proved isolating witnesses, 0 residue at both r.
                        CAUTION: the three runs use different ranges -- run(1,12), run(2,10),
                        run(3,9) -- so the counts 33 / 14 / 0 of non-standard labels are NOT
                        comparable across r. Non-standard labels need |nu| >= 2r+3 (Lemma 8.10), so
                        the usable slack above that threshold is 7 / 3 / 0: the zero at r=3 is the
                        range stopping at the threshold, not the obstruction disappearing.

typeD_rule.sage         The type-D table on this alphabet, obtained rather than recalled: every label
                        is placed in one of four classes and its value computed.
                        -> 139 labels at r=1 and 67 at r=2; 18 and 3 of them have no associate in
                           range.

typeD_residue.sage      Those 21 residual labels, reduced exactly.
                        -> 18/18 at r=1 and 3/3 at r=2, with integer coefficients.

unstable_closed.sage    The finite criterion of Section 8 -- the vanishing of every C_mu -- against
                        the exact object, in both the stable and the unstable range.
                        -> 94 + 90 = 184 shapes, 0 disagreements.

sec8_formulas.sage      Every displayed formula of Section 8 re-derived from its definitions rather
                        than quoted: the recurrence, and the sign of the column operation that turns
                        one character type into the next.
                        -> 0 failures on both.

rank_one_law.sage       The rank-one law: at r = 1 the object is +- a genuine character, so a single
                        order-two specialization decides whether it vanishes. Reports the shapes the
                        law covers and the identity holding on them.
                        -> 36/36 shapes vanish as predicted; 272 shapes, 0 mismatches.

nonstandard_survive.sage  The control that keeps Lemma 8.10 honest: non-standard labels do NOT all
                        die, so the reduction cannot be waved through. Prints the smallest survivor.

karmakar_caseB_min.sage The minimum of case (B), obtained rather than asserted: the staircase
                        (2r-1, ..., 1), of sizes 1, 6, 15, 28, 45, so r(2r-1) = 15 at r = 3.

graph_datum.sage        Whether a graph datum with an orientation separates the shapes that the sign
                        alone does not: sgn(a1-b1) by itself leaves 3 collisions at t = 2, and the
                        oriented datum leaves 0.

--------------------------------------------------------------------------------------------------
FIGURES
--------------------------------------------------------------------------------------------------

fig_data.py             Shared data generator for the image of the compression map.
fig_image.py            The image of the compression map.
fig_locus.py            The independence locus. Recomputes it rather than reading it off a table:
                        both sides are evaluated for every two-row partition in range.
fig_signed.py           The signed counts of the Example. Computes eps * d1 d2 d3 from scratch, so
                        the picture is data rather than assertion; prints a cross-check (16/16).
fig_data_new.sage       Exact data for the four Section-8 figures -> fig_data_new.json.
fig_zeros.py            The zero locus counted by |lambda|.
fig_phase.py            Where one reciprocal pair stops being typical. Each panel names its smallest
                        endpoint-only shape; the rest are listed in the accompanying remark.
fig_involution.py       The reflection involution, on a self-complementary shape of odd width and on
                        the same shape with one box moved.
fig_reduction.py        The type-D reduction as stacked bars. Bar heights are label counts over the
                        range each panel was computed on, which is smaller than the range of the
                        corresponding row of Section 9; the classification is what the picture is
                        about.
fig_fibres.py           The fibres of the evaluation invariant, the complement of fig_image: the
                        same range, with |lambda| on the vertical axis, so a column is a set of
                        partitions of different sizes carrying one value. Carries its own controls:
                        the bialternant is evaluated on all 860 (t=3) and 883 (t=4) partitions and
                        must be constant on every fibre, and the count of distinct VALUES is
                        compared against the count of invariants -- at t=4, 121 invariants take
                        only 110 values, because the closed form is symmetric in d1, d2, d3.

--------------------------------------------------------------------------------------------------
outputs/
--------------------------------------------------------------------------------------------------

The full stdout of every script above, as run on 2026-07-30. Every count in the verification table
appears in one of these files; a few are sums of printed lines rather than a single printed total,
and where that is so the summands are the per-r lines of the same file (9961 = 2157 + 4542 + 3262,
318 = 242 + 55 + 21, 372 = 76 + 144 + 152, 14950 = the ten buckets, 184 = 94 + 90, 236 = 139 + 97,
80 = 55 + 25, 21 = 18 + 3).

Regenerate with:

    python _audit_group1.py                  # group 1 -> _out/, then outputs/
    sage _audit_table.sh                     # group 2, inside the Sage container
    bash _save_outputs.sh                    # consolidate into outputs/

--------------------------------------------------------------------------------------------------
NOTES
--------------------------------------------------------------------------------------------------

Conventions follow the paper: t is the order of the root of unity and z is the free reciprocal
variable. (Much of the surrounding literature uses t for the free variable; the scripts do not.)

Two routines carry a caution worth repeating here. In `theorem_full.setup` the two distinguished
classes are ordered by column, which is what the proof of Theorem 3.1 needs; the short sign of
Proposition 3.11 instead requires them ordered by residue, and `sign_ayyer_idiom.py` reorders them for
that reason. Using the wrong order makes the short form fail -- which is the control that script
prints first. Note that `lambda11` is antisymmetric under exchanging the two blocks, so it must be
recomputed after any reordering; reusing it is a silent error.

Determinants are computed by Gaussian elimination with partial pivoting rather than by mpmath's
`det`, because the singular case is precisely the one that has to be scored rather than raised.

Author: Carles Marin, with Claude (Anthropic) as an AI research assistant.
