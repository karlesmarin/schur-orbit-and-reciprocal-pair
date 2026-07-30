#!/bin/bash
# Runs prove_W.sage (the isolating-witness / unstable-residue count) and saves its full output.
# Carles Marin + Claude (AI assistant).
mkdir -p _out
timeout 900 sage prove_W.sage > _out/prove_W.txt 2>&1
echo "exit=$?"
tail -20 _out/prove_W.txt
