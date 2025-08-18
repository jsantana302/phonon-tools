#!/bin/tcsh

num=$(echo $1 + $2 | bc)

echo "Cs3Sb" > POSCAR-01
echo "1.00000000000000" >> POSCAR-01
echo -e $num >> POSCAR-01

grep -v $1 POSCAR-00 | grep -v Cs3Sb | grep -v 1.00000000000000 >> POSCAR-01 #busca las líneas que no contienen NSW o EDIFF o LMAXMIX en el INCAR de OPT_PRIM y las copia en un nuevo INCAR

