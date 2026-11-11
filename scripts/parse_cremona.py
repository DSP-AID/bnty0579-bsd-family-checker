"""Build data/families.json from live-fetched Cremona database heads (ecdata repo).

Input: raw eclib table fragments downloaded by scripts/fetch_cremona.sh
  allbsd.00000-09999 head : "N class num [a1,a2,a3,a4,a6] R T CP OM L1 REG SHA"
  aplist.00000-09999 head : "N class AP25(+-)"   (+/- = W_p eigenvalue at bad p)
  count.00000-09999       : "N NNF"  (isogeny classes per conductor)
Output: families.json with every rank-0 curve of conductor <= MAXN, plus the
per-conductor class counts and a_p tables used as independent cross-checks.
"""
from __future__ import annotations

import json
import os
import re
import sys

MAXN = 100
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def lines(name):
    with open(os.path.join(DATA, name), encoding="utf-8", errors="replace") as f:
        for ln in f:
            ln = ln.strip()
            if ln:
                yield ln


def main():
    curves, classes = {}, {}
    for ln in lines("allbsd_head.txt"):
        parts = ln.split()
        n = int(parts[0])
        if n > MAXN:
            break
        ncl, num = parts[1], int(parts[2])
        ai = [int(x) for x in re.findall(r"-?\d+", parts[3])]
        assert len(ai) == 5, ln
        r, t, cp = int(parts[4]), int(parts[5]), int(parts[6])
        om, l1, _reg, sha = parts[7], parts[8], parts[9], parts[10]
        lab = f"{n}{ncl}{num}"
        classes.setdefault(n, set()).add(ncl)
        curves[lab] = {
            "label": lab,
            "lmfdb_label": f"{n}.{ncl}{num}",
            "conductor": n,
            "ainvs": ai,
            "rank": r,
            "torsion": t,
            "tamagawa_product": cp,
            "real_period": om,
            "l1_over": l1,  # L^{(r)}(E,1)/r! (Cremona L1 column)
            "sha": sha,
        }
    rank0 = {k: v for k, v in curves.items() if v["rank"] == 0}

    ap = {}
    for ln in lines("aplist_head.txt"):
        parts = ln.split()
        n = int(parts[0])
        if n > MAXN:
            break
        ap[f"{n}{parts[1]}"] = parts[2:]

    counts = {}
    for ln in lines("count.txt"):
        n, nnf = ln.split()
        if int(n) <= MAXN:
            counts[int(n)] = int(nnf)

    fam = {
        "source": "Cremona ecdata (JohnCremona/ecdata@master, eclib tables), "
        "fetched live from raw.githubusercontent.com; LMFDB is downstream of this DB",
        "max_conductor": MAXN,
        "n_curves_all_ranks": len(curves),
        "n_rank0_curves": len(rank0),
        "classes_per_conductor_lmfdb_count": counts,
        "classes_seen_in_allbsd": {str(k): sorted(v) for k, v in sorted(classes.items())},
        "ap25": ap,
        "curves": rank0,
    }
    with open(os.path.join(DATA, "families.json"), "w", encoding="utf-8") as f:
        json.dump(fam, f, indent=1)
    print(
        f"rank-0 curves N<=100: {len(rank0)} in "
        f"{sum(1 for _ in classes)} conductors listed; total curves all ranks {len(curves)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
