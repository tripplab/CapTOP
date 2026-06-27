#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${CAPTOP_BIN:-$ROOT/build/captop}"
WORK="$ROOT/tests/stage6_persistence"
mkdir -p "$WORK"
if ! "$BIN" --version | grep -q 'GUDHI support: enabled'; then
  echo "Stage 6 persistence tests require a CAPTOP binary built with GUDHI support." >&2
  exit 0
fi

echo "Stage 6 test harness is ready. Add/refresh synthetic meshes in $WORK and run captop persist cases."
