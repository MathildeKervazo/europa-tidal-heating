#!/bin/sh
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=${1:-$ROOT/work}
JOBS=${JOBS:-1}
export KEEP_ALL=${KEEP_ALL:-0}
make -C "$ROOT/src" FC="${FC:-gfortran}" > /dev/null
BIN="$ROOT/src/Forced_tide"
export BIN
mkdir -p "$WORK"
cp -R "$ROOT/input/." "$WORK/"
find "$WORK" -name input.nml -print0 | xargs -0 -n 1 -P "$JOBS" sh -c '
  d=$(dirname "$1")
  cd "$d" && mkdir -p Melt Htide_Output/Tot && "$BIN" > Forced_tide.log || { echo "FAILED: $d" >&2; exit 1; }
  if [ "$KEEP_ALL" != 1 ]; then rm -f yi-real.dat yi-imag.dat Hmu_htide.txt k2_Q_Ptide.txt Forced_tide.log; fi
  echo "done: $d"
' sh
echo "all calculations finished in $WORK"
