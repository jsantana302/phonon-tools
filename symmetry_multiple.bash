#!/bin/bash

for dir in comp_* exp_*; do
    if [ -d "$dir" ]; then
        (cd "$dir" && ~/scripts/symmetry.bash)
    fi
done

