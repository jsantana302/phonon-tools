#!/bin/bash
# Create a temporary file
tmpfile=$(mktemp)

# Extract first 2000 lines
head -n 2000 band.yaml > "$tmpfile"

# Extract and print frequency values using grep
grep -oP 'frequency:\s*\K[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?' "$tmpfile"

# Optionally remove the temporary file
rm "$tmpfile"

