#!/bin/sh
# Fetch the family data live from the Cremona eclib database (the source DB of
# LMFDB's elliptic-curve pages; CC-BY-SA / Artistic-2.0 as marked upstream).
# Only the head of each conductor-range file is needed (conductor <= 100).
set -e
cd "$(dirname "$0")/../data"
BASE=https://raw.githubusercontent.com/JohnCremona/ecdata/master
curl -s -r 0-499999 -o allbsd_head.txt $BASE/allbsd/allbsd.00000-09999
curl -s -r 0-499999 -o aplist_head.txt $BASE/aplist/aplist.00000-09999
curl -s -r 0-199999 -o count.txt $BASE/count/count.00000-09999
curl -s -r 0-499999 -o alllabels_head.txt $BASE/alllabels/alllabels.00000-09999
{ date -u +"fetched %Y-%m-%dT%H:%M:%SZ from $BASE"; sha256sum allbsd_head.txt aplist_head.txt count.txt; } > PROVENANCE.txt
wc -l allbsd_head.txt
