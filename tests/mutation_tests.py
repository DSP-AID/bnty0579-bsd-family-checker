"""Mutation suite: each mutation of a stored cross-check value must make the
checker FAIL (the checks bite), and the unmutated curve must PASS.

Run: python3 tests/mutation_tests.py
"""
from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAM = os.path.join(ROOT, "data", "families.json")
CHECKER = os.path.join(ROOT, "src", "bsd_family.py")


def mutations():
    """(name, label, apply-fn) - each breaks one cross-checked quantity."""

    def setf(field, val):
        def f(d, lab):
            d["curves"][lab][field] = val

        return f

    def mul(field, factor):
        def f(d, lab):
            v = d["curves"][lab][field]
            d["curves"][lab][field] = str(float(v) * factor)

        return f

    def perturb(field, rel):
        def f(d, lab):
            v = float(d["curves"][lab][field])
            d["curves"][lab][field] = repr(v * (1 + rel))

        return f

    def a4_shift(k):
        def f(d, lab):
            d["curves"][lab]["ainvs"][3] += k

        return f

    def flip_ap():
        def f(d, lab=None):
            tab = d["ap25"]["11a"]
            for i, tok in enumerate(tab):
                if tok.lstrip("-").isdigit():
                    tab[i] = str(-int(tok))
                    break

        return f

    return [
        ("tamagawa+2", "11a1", setf("tamagawa_product", 7)),
        ("torsion=10", "11a1", setf("torsion", 10)),
        ("sha=4", "11a1", setf("sha", 4)),
        ("omega x2", "11a1", mul("real_period", 2.0)),
        ("L1 +1e-6", "11a1", perturb("l1_over", 1e-6)),
        ("a4 -9", "11a1", a4_shift(1)),
        ("N=13", "11a1", setf("conductor", 13)),
        ("ap25 flip", "11a1", flip_ap()),
        ("cprod 6->5", "14a2", setf("tamagawa_product", 5)),
        ("omega/2 (two-comp)", "14a2", mul("real_period", 0.5)),
        # regression coverage for the two bugs found in the family run:
        ("torsion 8->4 (Z/2xZ/4)", "15a1", setf("torsion", 4)),
        ("cprod 2->1 (III at p=2)", "32a2", setf("tamagawa_product", 1)),
    ]


def run(label, data_dir=None):
    cmd = [sys.executable, CHECKER, "--label", label]
    if data_dir:
        cmd += ["--data", data_dir]
    return subprocess.run(cmd, capture_output=True, text=True, check=False).returncode


def main():
    with open(FAM, encoding="utf-8") as f:
        fam = json.load(f)
    bad = 0
    for name, label, apply_fn in mutations():
        d = copy.deepcopy(fam)
        d["curves"] = {label: d["curves"][label]}
        apply_fn(d, label)
        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, "families.json"), "w", encoding="utf-8") as f:
                json.dump(d, f)
            rc = run(label, tmp)
        ok = rc != 0
        bad += 0 if ok else 1
        print(f"{'bites' if ok else 'SILENT'}   {name:<20} ({label}) rc={rc}")
    # positive controls: the unmutated curves must pass (incl. additive-type and
    # non-cyclic-torsion regressions)
    for label in ("11a1", "14a2", "32a2", "15a1"):
        rc = run(label)
        ok = rc == 0
        bad += 0 if ok else 1
        print(f"{'clean' if ok else 'BROKEN'} {label:<20} control rc={rc}")
    print(f"\n{len(mutations()) + 4} tests, silent/broken: {bad}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
