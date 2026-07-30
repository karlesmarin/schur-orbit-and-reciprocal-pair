#!/bin/bash
# Re-runs every script that backs the second group of the paper's verification table, saving the
# FULL output of each so the table can be checked line by line against what the scripts print.
# Carles Marin + Claude (AI assistant).
mkdir -p _out
for f in selfcomp_law close_X_r1 associates_witness typeD_rule typeD_residue unstable_closed AUDIT_ALL; do
  if [ -f "$f.sage" ]; then
    echo "running $f"
    timeout 900 sage "$f.sage" > "_out/$f.txt" 2>&1
    echo "  exit=$? lines=$(wc -l < _out/$f.txt)"
  else
    echo "MISSING: $f.sage"
  fi
done
