#!/bin/sh
# Generalized BSD leading-term verifier for every rank-0 elliptic curve over Q
# of conductor <= 100 (generalization of the 32a2-only checker from bnty_0573).
#
#   ./verify.sh              check the whole family (offline Cremona snapshot)
#   ./verify.sh 32a2         one curve by Cremona label, with a live LMFDB call
#   ./verify.sh --refresh    re-download the family data, then check everything
#
# Exit code 0 = every computed invariant agrees with the cross-check sources
# and the BSD prediction is an integer square matching the stored sha.
set -e
cd "$(dirname "$0")"

if [ "$1" = "--refresh" ]; then
    sh scripts/fetch_cremona.sh
    python3 scripts/parse_cremona.py
    shift
fi

if [ -n "$1" ]; then
    exec python3 src/bsd_family.py --label "$1" --live
fi
exec python3 src/bsd_family.py
