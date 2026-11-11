#!/bin/sh
# Mutation suite entry point: proves the cross-checks bite.
set -e
cd "$(dirname "$0")/.."
python3 tests/mutation_tests.py
