#!/usr/bin/env bash
set -euo pipefail

PDOS_METHOD="eq"
OUTPUT="R2SCAN_QHA_phonons_exp.png"
labels=(eq exp_1 exp_2 exp_3 exp_5 exp_6 exp_7)
args=()
shopt -s nullglob

for label in "${labels[@]}"; do
  if [[ $label == eq ]]; then
    base_dir=comp_0; tag=eq
  else
    base_dir=$label; tag=$label
  fi

  ph_dirs=( "$base_dir"/ph_* )
  if (( ${#ph_dirs[@]} != 1 )); then
    echo "Warning: skipping '$base_dir' (found ${#ph_dirs[@]} ph_* dirs)" >&2
    continue
  fi

  args+=( "${tag}:${ph_dirs[0]}/" )
done

shopt -u nullglob

if (( ${#args[@]} == 0 )); then
  echo "Error: no phonon dirs found in any of the specified labels" >&2
  exit 1
fi

python ~/scripts/plot_phonons_multiple.py \
  --base-paths "${args[@]}" \
  --pdos-method "$PDOS_METHOD" \
  --output "$OUTPUT"
