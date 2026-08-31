#!/bin/bash
# check-omark-db.sh
#
# Standalone integrity check for an OMArk .h5 database (e.g. LUCA.h5).
# NOT part of the per-taxon 03-BUSCO-and-OMArk.sh pipeline — walking every
# leaf table and sampling a real data read is too slow to repeat on every
# taxon run. Run this by hand instead, whenever there's reason to doubt the
# file (e.g. right after a fresh/resumed download).
#
# A bare "can it be opened" check isn't enough here: a truncated download
# can still open fine, since PyTables' table-of-contents lives near the
# file's header/pre-allocated space, not necessarily at the very end. This
# forces an actual read from every leaf table, which is what catches a
# missing tail.
#
# Usage: bash check-omark-db.sh [path-to-h5]
#   Defaults to /bank/ncbi/gnathostome/OMArk/LUCA.h5

set -euo pipefail

DB_PATH="${1:-/bank/ncbi/gnathostome/OMArk/LUCA.h5}"

if [[ ! -s "$DB_PATH" ]]; then
    echo "[FAIL] ${DB_PATH} does not exist or is empty" >&2
    exit 1
fi

echo "[INFO] Checking ${DB_PATH} ($(du -h "$DB_PATH" | cut -f1))..." >&2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate omark

set +e
python3 - "$DB_PATH" << 'PYEOF'
import sys
import tables

path = sys.argv[1]

try:
    f = tables.open_file(path, mode="r")
except Exception as e:
    print(f"[FAIL] could not open {path}: {e}")
    sys.exit(1)

bad = []
n_leaves = 0
for node in f.walk_nodes("/", classname="Leaf"):
    n_leaves += 1
    try:
        if hasattr(node, "shape") and node.shape != ():
            _ = node[tuple(0 for _ in node.shape)]
        else:
            _ = node.read()
    except Exception as e:
        bad.append((node._v_pathname, str(e)))

f.close()

print(f"Walked {n_leaves} leaf node(s), sampled a data read from each.")
if bad:
    print(f"[FAIL] {len(bad)} broken node(s):")
    for name, err in bad[:10]:
        print(f"  {name}: {err}")
    sys.exit(1)

print(f"[PASS] {path} looks structurally intact ({n_leaves} nodes checked)")
PYEOF
status=$?
set -e

conda deactivate
exit "$status"
