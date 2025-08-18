#!/usr/bin/env bash
set -euo pipefail

PDOS_METHOD=comp_6
OUTPUT="QHA_phonons_pdos_comp.png"
labels=(eq comp_3 comp_5 comp_6)
args=()
shopt -s nullglob

for label in "${labels[@]}"; do
    comp_dir=$([[ $label == eq ]] && echo "comp_0" || echo "$label")
    ph_dirs=( "$comp_dir"/ph_* )
    if (( ${#ph_dirs[@]} != 1 )); then
        echo "Warning: skipping '$comp_dir' (found ${#ph_dirs[@]} dirs)" >&2
        continue
    fi
    args+=( "${label}:${ph_dirs[0]}/" )
done

shopt -u nullglob

if (( ${#args[@]} == 0 )); then
    echo "Error: no dirs found" >&2
    exit 1
fi

python ~/scripts/plot_phonons_multiple_dos.py \
    --base-paths "${args[@]}" \
    --pdos-method "$PDOS_METHOD" \
    --output "$OUTPUT"
