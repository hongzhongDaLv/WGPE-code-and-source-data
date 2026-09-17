#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python "$ROOT/tests/smoke_test.py"
python "$ROOT/tests/test_september_update.py"
