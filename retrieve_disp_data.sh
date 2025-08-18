#!/bin/bash

# Create a directory called plot_data if it doesn't exist
mkdir -p plot_data

# Loop through directories named dispersion_XXX
for dir in dispersion_*; do
  if [ -d "$dir" ]; then
    # Extract the temperature value from the directory name
    temp="${dir#dispersion_}"
    # Create a directory for the temperature inside plot_data
    mkdir -p "plot_data/dispersion_${temp}"
    # Copy the phonopy.yaml and FORCE_CONSTANTS files to the new directory
    cp "${dir}/phonopy.yaml" "plot_data/dispersion_${temp}/"
    cp "${dir}/FORCE_CONSTANTS" "plot_data/dispersion_${temp}/"
  fi
done

echo "Files copied successfully."

