#!/bin/bash
# Consolidates every backing script's saved stdout into anc/outputs/, which is the durable
# evidence layer for the paper's verification table: every number in that table must be
# greppable here.  Regenerate with _audit_group1.py (pure Python) and, in the Sage container,
# ayyer_bridge/_audit_table.sh + _run_provew.sh.
# Carles Marin + Claude (AI assistant).
set -e
ANC="E:/proyectos/Curiosity/research/orbit-pair/paper/anc"
BR="E:/proyectos/Curiosity/research/smeft_formalization/ayyer_bridge"
OUT="$ANC/outputs"
mkdir -p "$OUT"

# group 1 -- pure Python, ancillary to the paper
for f in "$ANC"/_out/*.txt; do [ -f "$f" ] && cp "$f" "$OUT/$(basename "$f")"; done
# group 2 -- Sage, the zero-locus section
for f in "$BR"/_out/*.txt; do [ -f "$f" ] && cp "$f" "$OUT/$(basename "$f")"; done

echo "saved:"; ls -1 "$OUT" | sed 's/^/   /'
echo "total files: $(ls -1 "$OUT" | wc -l)"
