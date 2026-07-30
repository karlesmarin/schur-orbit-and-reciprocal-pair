Ancillary material for

    "Factorization of Schur polynomials twisted by roots of unity and a reciprocal pair"
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
    sage script.sage                             # group 2

Run them from this directory.

Section and statement numbers below refer to the preprint, which is forthcoming; this repository
carries the scripts and their saved output, not the paper. The published rank-one case at t=2 is
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

falsify.py              The controls (D3) and (D4) of Section 7: the same test applied to the coset
                        alphabet and to a free (non-reciprocal) pair.
                        -> orbit 600/0, coset 383/217, free pair 200/400.          (Sections 7, 9)

d_from_quotient.py      Proposition 3.4: the three arguments read off the t-quotient.
                        -> 2970/2970 in the two-class profile.                     (Section 3.1)

sign_ayyer_idiom.py     Question 3.7: the short form of the sign in the notation of [AK25]. Prints a
                        CONTROL first -- the same comparison with the two blocks in arbitrary order,
                        which fails 592/904 and is meant to -- then the ordered comparison, then the
                        cell test. The control is what makes the ordering hypothesis load-bearing.
                        -> control 592/904; ordered 1496/1496; 112 cells, 0 mixed.  (Section 3)

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

check_nonstd_bound.py   Lemma 8.6 against the archived outputs: every non-standard label printed by
                        the group-2 scripts satisfies |nu| >= 2r+3, and every residue shape satisfies
                        |lambda| >= 2r+3.
                        -> 0 violations.                                           (Section 8.4)

sieve_counts.py         The endpoint sieve quantified, from the same exact data the Section-8 figures
                        use. Runs the nesting control first (a shape vanishing identically must also
                        vanish at the endpoint).
                        -> control 0; flags 143 / 22 / 9 shapes at r=1,2,3, of which 0 / 4 / 3 are
                           spurious, i.e. 0% / 18% / 33%.                          (Section 8)

check_layout.py         Not a mathematical check. Reports pages carrying a large blank band inside
                        the text area, which a clean LaTeX run does not detect.
                        Usage: python check_layout.py orbit_pair.pdf

check_refs.py           Not a mathematical check either. Two things LaTeX cannot tell you: whether a
                        cross-reference points at the right KIND of object ("Theorem \ref{lem:...}"
                        compiles and is wrong), and whether any bibliography entry is never cited.
                        Usage: python check_refs.py orbit_pair.tex orbit_pair_es.tex

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

prove_W.sage            Conjecture 8.7, the isolating-mu statement, at r=2 and r=3, and the count of
                        non-standard labels that survive on this alphabet.
                        -> 55 + 25 = 80 proved isolating witnesses, 0 residue at both r.
                        CAUTION: the three runs use different ranges -- run(1,12), run(2,10),
                        run(3,9) -- so the counts 33 / 14 / 0 of non-standard labels are NOT
                        comparable across r. Non-standard labels need |nu| >= 2r+3 (Lemma 8.6), so
                        the usable slack above that threshold is 7 / 3 / 0: the zero at r=3 is the
                        range stopping at the threshold, not the obstruction disappearing.

typeD_rule.sage         The type-D table on this alphabet, obtained rather than recalled: every label
                        is placed in one of four classes and its value computed.
                        -> 139 labels at r=1 and 67 at r=2; 18 and 3 of them have no associate in
                           range.

typeD_residue.sage      Those 21 residual labels, reduced exactly.
                        -> 18/18 at r=1 and 3/3 at r=2, with integer coefficients.

unstable_closed.sage    The criterion (13) against the exact object, in both the stable and the
                        unstable range.
                        -> 94 + 90 = 184 shapes, 0 disagreements.

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

--------------------------------------------------------------------------------------------------
outputs/
--------------------------------------------------------------------------------------------------

The full stdout of every script above, as run on 2026-07-30. Every count in the verification table
appears in one of these files; a few are sums of printed lines rather than a single printed total,
and where that is so the summands are the per-r lines of the same file (9961 = 2157 + 4542 + 3262,
318 = 242 + 55 + 21, 372 = 76 + 144 + 152, 14950 = the ten buckets, 184 = 94 + 90, 206 = 139 + 67,
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
Question 3.7 instead requires them ordered by residue, and `sign_ayyer_idiom.py` reorders them for
that reason. Using the wrong order makes the short form fail -- which is the control that script
prints first. Note that `lambda11` is antisymmetric under exchanging the two blocks, so it must be
recomputed after any reordering; reusing it is a silent error.

Determinants are computed by Gaussian elimination with partial pivoting rather than by mpmath's
`det`, because the singular case is precisely the one that has to be scored rather than raised.

Author: Carles Marin, with Claude (Anthropic) as an AI research assistant.
