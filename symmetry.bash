#!/bin/bash
mkdir crys 2>/dev/null
cd crys 2>/dev/null

# Copy and prepare the input file (suppress errors if files are missing)
cp ../PBE-1.restart ./PBE.inp 2>/dev/null
cp ../PBE0-1.restart ./PBE.inp 2>/dev/null
cp ../R2SCAN-1.restart ./PBE.inp 2>/dev/null

sed -i '/IGNORE_CONVERGENCE_FAILURE/d' PBE.inp
sed -i '/&AUXILIARY_DENSITY_MATRIX_METHOD/,/&END AUXILIARY_DENSITY_MATRIX_METHOD/d' PBE.inp
sed -i '/MIN_PAIR_LIST_RADIUS -1/d' PBE.inp

# Define an array of tolerances in reverse order
tolerances=(
  0.1
  0.05
  0.01
  0.005
  0.001
  0.0001
  0.0005
  0.00001
  0.00005
)

# Loop through the tolerances and execute phonopy
for tolerance in "${tolerances[@]}"; do
  output_file="sym_${tolerance}.out"
  phonopy --cp2k -c PBE.inp --symmetry --tolerance ${tolerance} > "${output_file}"
  mv "Punitcell.inp" "Punit_${tolerance}.inp" 2>/dev/null
  mv "Bunitcell.inp" "Bunit_${tolerance}.inp" 2>/dev/null
done

~/scripts/file_conversion/cp2k2xyz.py > conversion.out 
mkdir -p Bunit 2>/dev/null

# Move files, suppressing errors if source/destination issue occurs
mv Bunit* Bunit 2>/dev/null
