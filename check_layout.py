# -*- coding: utf-8 -*-
"""Layout audit: finds pages with a large blank band inside the text area.

A clean LaTeX compile reports overfull boxes but says nothing about half-empty pages caused by float
congestion, which is what a reader notices first. For each page this measures the biggest vertical
gap between consecutive content blocks, ignoring the running head, and reports any gap above a
threshold as a fraction of the text height.

TWO shapes of blank, because for a while it only saw one. A hole BETWEEN two blocks is the float
congestion above. A blank TAIL -- content that stops early with nothing after it -- has no hole at
all, so the between-blocks measure returns zero and the page passes. That was page 54 of the arXiv
edition sitting at 18 % of its height with nine guards green: a forty-row tabular is an indivisible
box, so LaTeX moved it whole to the next page rather than splitting it, and left the section heading
alone on its own. The fix there was longtable; the fix here is to measure the tail as well. The
control lives in the scratchpad as ctl.tex -- the same source with that one table reverted -- and it
reports page 54 at 90 %, which is how we know this check can fail.

    python check_layout.py                    # every edition built by build_all.sh
    python check_layout.py orbit_pair.pdf orbit_pair_Z.pdf ...
    python check_layout.py --threshold 0.2 orbit_pair.pdf

It takes SEVERAL files, and with no argument it audits them all, because the one defect it has
caught was in an edition nobody ran it on. The long-abstract editions for Zenodo carried a first
page with the title and nothing else -- 73 % of the column blank -- from the day they were first
built: amsart sets the abstract in a single unbreakable box, so an abstract that does not fit under
the title moves whole to page 2. This script would have said so in one line at any point. It was
called with its default argument, which is the arXiv edition, and the other four were never asked.
A checker that is only ever pointed at one of five artifacts is a checker for one artifact.

Exits non-zero if anything is reported, so a build script can gate on it.

Authors: Carles Marin, Claude (AI assistant)."""
import glob
import sys

import fitz

DEFAULT = ["orbit_pair.pdf", "orbit_pair_es.pdf", "orbit_pair_Z.pdf", "orbit_pair_Z_es.pdf",
           "orbit_pair_bilingual.pdf"]

args = sys.argv[1:]
THR = 0.14
if "--threshold" in args:
    k = args.index("--threshold")
    THR = float(args[k + 1])
    del args[k:k + 2]
paths = args or [p for p in DEFAULT if glob.glob(p)]

bad = 0
for path in paths:
    doc = fitz.open(path)
    worst = []
    for i, page in enumerate(doc):
        H = page.rect.height
        head, foot = 0.10 * H, 0.93 * H          # drop the running head; keep the folio
        blocks = []
        for b in page.get_text("blocks"):
            x0, y0, x1, y1 = b[:4]
            text = b[4].strip() if len(b) > 4 else ""
            if not text:
                continue
            if y1 < head or y0 > foot:
                continue
            blocks.append((y0, y1))
        # images and drawings count as content too
        for d in page.get_image_info():
            bb = d.get("bbox")
            if bb and bb[3] > head and bb[1] < foot:
                blocks.append((bb[1], bb[3]))
        if len(blocks) < 2:
            continue
        blocks.sort()
        # merge, then find the largest hole
        merged = [list(blocks[0])]
        for y0, y1 in blocks[1:]:
            if y0 <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], y1)
            else:
                merged.append([y0, y1])
        gap, at = 0.0, None
        for a, b in zip(merged, merged[1:]):
            g = b[0] - a[1]
            if g > gap:
                gap, at = g, a[1]
        frac = gap / (foot - head)
        if frac >= THR:
            worst.append((frac, i + 1, gap, at, "band"))

        # And the blank that is NOT between two blocks: the one after the last one.  A page whose
        # content stops early has no internal hole at all, so the loop above sees nothing, and that
        # is precisely how p.54 sat at 18 % of its height through nine green guards until a reader
        # noticed.  The cause there was a forty-row tabular, an indivisible box that LaTeX moves
        # whole rather than splitting.  The last page of the document is exempt: ending short is
        # what a last page does.
        if i + 1 < doc.page_count:
            tail = foot - merged[-1][1]
            tfrac = tail / (foot - head)
            if tfrac >= THR:
                worst.append((tfrac, i + 1, tail, merged[-1][1], "tail"))

    worst.sort(reverse=True)
    print("%s : %d pages, threshold %.0f%% of the text height"
          % (path, doc.page_count, THR * 100))
    if not worst:
        print("  no page has a blank band or a short tail above the threshold.")
    for frac, pg, gap, at, kind in worst:
        what = ("blank band" if kind == "band" else
                "content stops early, blank tail")
        print("  page %3d : %s %.0f%% of the text height (%.0fpt, starting at y=%.0f)"
              % (pg, what, frac * 100, gap, at))
    bad += len(worst)

print("TOTAL pages with a blank band or a short tail: %d  (over %d file%s)"
      % (bad, len(paths), "" if len(paths) == 1 else "s"))
sys.exit(1 if bad else 0)
