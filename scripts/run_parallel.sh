#!/bin/sh
# Shard the family across N workers (Xeon offload box) and merge results.
# Usage: sh scripts/run_parallel.sh [NWORKERS]   (default 16)
# Writes tmp/results_full.json (merged, families.json order) + per-shard logs.
set -e
cd "$(dirname "$0")/.."
NW=${1:-16}
mkdir -p tmp/chunks

python3 - "$NW" <<'EOF'
import json, sys
nw = int(sys.argv[1])
fam = json.load(open("data/families.json", encoding="utf-8"))
labels = list(fam["curves"])
# round-robin balances heavy conductors across shards
for i in range(nw):
    open(f"tmp/chunks/chunk_{i}.txt", "w").write(" ".join(labels[i::nw]))
print(f"{len(labels)} labels -> {nw} chunks")
EOF

for i in $(seq 0 $((NW - 1))); do
    python3 -u src/bsd_family.py --label $(cat tmp/chunks/chunk_$i.txt) \
        --out tmp/res_$i.json > tmp/log_$i.txt 2>&1 &
done
wait

python3 - <<'EOF'
import json, glob
fam = json.load(open("data/families.json", encoding="utf-8"))
rows = {}
fails = 0
for fn in glob.glob("tmp/res_*.json"):
    for r in json.load(open(fn, encoding="utf-8")):
        rows[r["label"]] = r
        ok = (r["conductor_rederived"] and r["sha_is_square"]
              and r["fe_selfcheck"] and r["torsion_agrees"]
              and r["mismatches"] == "none")
        fails += 0 if ok else 1
ordered = [rows[l] for l in fam["curves"] if l in rows]
json.dump(ordered, open("tmp/results_full.json", "w"), indent=1)
print(f"merged {len(ordered)}/{len(fam['curves'])} rows; failures: {fails}")
missing = [l for l in fam["curves"] if l not in rows]
if missing:
    print("MISSING:", missing)
EOF

grep -h "FAIL" tmp/log_*.txt || echo "no FAIL lines"
