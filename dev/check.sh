#!/usr/bin/env bash
# Full check: compile the Java modules, run the Clojure suite against the
# nibabel corpus, then have nibabel re-read everything nicloj wrote.
#
# Needs `clojure` and `javac` on PATH, plus a Python env with nibabel and
# nilearn for the corpus steps.
set -euo pipefail
cd "$(dirname "$0")/.."

CLOJURE=${CLOJURE:-clojure}
# Prefer the project venv, so this works the same on the host and inside a
# container where `python` may not have nibabel.
if [ -z "${PYTHON:-}" ] && [ -x .venv/bin/python ]; then
  PYTHON=$PWD/.venv/bin/python
fi
PYTHON=${PYTHON:-python}

echo "== javac =="
$CLOJURE -T:build javac

if [ ! -f test-data/manifest.edn ]; then
  echo "== generating test corpus =="
  $PYTHON scripts/gen_testdata.py
fi

echo "== clojure tests =="
$CLOJURE -M:test

echo "== nibabel verification =="
$PYTHON scripts/verify_roundtrip.py
