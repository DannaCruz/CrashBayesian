#!/usr/bin/env bash
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
out="$repo/model-output"
mkdir -p "$out"

if ! command -v make >/dev/null; then
  echo "make is required" >&2
  exit 1
fi

: "${CMDSTAN:?Set CMDSTAN to the CmdStan installation directory}"

make -C "$CMDSTAN" "$repo/models/nb_nonspatial_article" "$repo/models/nb_icar_disconnected_article"

run_four_chains() {
  local executable="$1" data="$2" name="$3" seed="$4" delta="$5"
  mkdir -p "$out/$name"
  local pids=()
  for chain in 1 2 3 4; do
    "$executable" sample num_warmup=500 num_samples=500 thin=1 \
      adapt delta="$delta" algorithm=hmc engine=nuts max_depth=10 stepsize=0.01 \
      data file="$data" init=0 random seed=$((seed + chain)) id="$chain" \
      output file="$out/$name/chain_${chain}.csv" refresh=200 \
      >"$out/$name/chain_${chain}.log" 2>&1 &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do wait "$pid"; done
}

run_four_chains "$repo/models/nb_nonspatial_article" "$out/data/intersection_nonspatial.json" intersection_nonspatial 202608240 0.95
run_four_chains "$repo/models/nb_nonspatial_article" "$out/data/segment_nonspatial.json" segment_nonspatial 202608250 0.95
run_four_chains "$repo/models/nb_icar_disconnected_article" "$out/data/intersection_car.json" intersection_car 202608260 0.97
run_four_chains "$repo/models/nb_icar_disconnected_article" "$out/data/segment_car.json" segment_car 202608270 0.97

echo "All models completed."

