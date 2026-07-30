# -*- coding: utf-8 -*-
"""Do the two language editions still say the same thing?

Corrections get applied to one edition and not propagated to the other, and nothing in a LaTeX build
notices. This compares the two sources on things that must match regardless of language:

  * every number of two or more digits, with multiplicity -- a count that appears in one edition and
    not the other means a sentence was added, dropped or left stale;
  * the number of numbered environments of each kind, and of labels -- a missing theorem or a missing
    paragraph carrying one shows up here;
  * the set of \\label keys, which must be identical for the bilingual build to work.

It found two real defects on 2026-07-30: the Spanish edition was missing the paragraph explaining the
verification table's shaded block together with the "3414 checks, 0 failures" audit, and its open
problem on the other classical types still said "we have not attempted this" after the English had
been corrected to report the measurement that made that sentence false.

    python check_parity.py orbit_pair.tex orbit_pair_es.tex

Authors: Carles Marin, Claude (AI assistant)."""
import collections, io, re, sys

ENVS = ["theorem", "lemma", "corollary", "proposition", "conjecture", "remark", "question",
        "problem", "example", "figure", "equation"]


def load(path):
    src = io.open(path, encoding="utf-8").read().split("\\begin{thebibliography}")[0]
    src = re.sub(r"(?<!\\)%.*", "", src)          # drop comments, keep \%
    return src


def main(a_path, b_path):
    A, B = load(a_path), load(b_path)
    bad = 0

    na = collections.Counter(re.findall(r"(?<![\d.])(\d{2,})(?![\d.])", A))
    nb = collections.Counter(re.findall(r"(?<![\d.])(\d{2,})(?![\d.])", B))
    diff = sorted((k for k in set(na) | set(nb) if na[k] != nb[k]),
                  key=lambda k: -abs(na[k] - nb[k]))
    print("numbers of two or more digits, by multiplicity")
    if not diff:
        print("   identical")
    for k in diff:
        bad += 1
        print("   %-10s %s=%d   %s=%d" % (k, a_path, na[k], b_path, nb[k]))

    print("\nnumbered environments")
    for e in ENVS:
        ca = len(re.findall(r"\\begin\{%s\}" % e, A))
        cb = len(re.findall(r"\\begin\{%s\}" % e, B))
        if ca != cb:
            bad += 1
            print("   %-12s %s=%d   %s=%d   MISMATCH" % (e, a_path, ca, b_path, cb))
    else:
        pass
    if all(len(re.findall(r"\\begin\{%s\}" % e, A)) == len(re.findall(r"\\begin\{%s\}" % e, B))
           for e in ENVS):
        print("   all counts equal")

    la = set(re.findall(r"\\label\{([^}]+)\}", A))
    lb = set(re.findall(r"\\label\{([^}]+)\}", B))
    print("\nlabels")
    if la == lb:
        print("   identical (%d)" % len(la))
    else:
        for k in sorted(la - lb):
            bad += 1
            print("   only in %s: %s" % (a_path, k))
        for k in sorted(lb - la):
            bad += 1
            print("   only in %s: %s" % (b_path, k))

    print("\nTOTAL divergences: %d" % bad)
    return bad


if __name__ == "__main__":
    args = sys.argv[1:] or ["orbit_pair.tex", "orbit_pair_es.tex"]
    sys.exit(0 if main(args[0], args[1]) == 0 else 1)
