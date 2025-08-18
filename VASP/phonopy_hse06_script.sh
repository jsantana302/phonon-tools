#!/bin/tcsh
set echo

# Activate the conda environment
conda activate ~/miniconda3/envs/phonopyenv

# Ask for user input to determine if it's pre or post processing
echo "pre or post?"
set phonopy = $<

# Perform operations for pre processing
if ($phonopy == "pre") then

   # Copy CONTCAR to POSCAR
   cp ../CONTCAR ./POSCAR

   # Perform phonopy symmetry operations
   phonopy --symmetry --tolerance 0.001
   sleep 4 # Wait for phonopy to complete

   # Ask user for cell type
   echo "prim or conv?"
   set cell = $<

   # Choose POSCAR file based on cell type
   if ($cell == "prim") then
      cp PPOSCAR POSCAR
      rm BPOSCAR
   else if ($cell == "conv") then
      cp BPOSCAR POSCAR
      rm PPOSCAR
   else
      echo "I don't understand"
      exit
   endif
   
   # Ask for supercell size
   echo "supercell size"
   set s1 = $<:q
   set s2 = $<:q
   set s3 = $<:q

   # Run phonopy with provided supercell size
   phonopy -d --dim="$s1 $s2 $s3"
   sleep 4
   
   # Count POSCAR files and create a directory for each one
   set nPoscar = `ls POSCAR-* | wc -l`
   mkdir 0
   echo 0
   mv SPOSCAR 0/POSCAR

   # Move POSCAR files to respective directories
   foreach i ( `seq 1 $nPoscar` )
      mkdir $i
      echo $i
      if ($i < 10) then
         set prefix="00"
      else if ($i < 100) then
         set prefix="0"
      else
         set prefix=""
      endif
      mv "POSCAR-$prefix$i" "$i/POSCAR"
   end  

   # Copy necessary files for VASP calculations
   cp ../../files_VASP/KPOINTS . 
   cp ../../files_VASP/INCAR.* .
   cp ../../files_VASP/vasp.job .
   cp ../../files_VASP/POTCAR .

   # Create symbolic links in each directory and submit jobs
   foreach dir ( `seq 0 $nPoscar` )
      cd $dir
      ln -sf ../KPOINTS . 
      ln -sf ../INCAR.1 ./INCAR
      ln -sf ../POTCAR .
      ln -sf ../vasp.job .   
      sbatch -J $dir vasp.job
      cd ..
   end

# Perform operations for post processing
else if ($phonopy == "post") then
   
   echo "How many directories?"
   set num = $<
   
   ## Set the variable for the command
   set cmd = "phonopy -f"
   ## Append the folder names to the end of the command
   foreach i ( `seq 1 $num` )
      set cmd = "$cmd $i/vasprun.xml"
      end
   echo "$cmd"
   $cmd
   sleep 5
   
   cp home/nipj5103/scripts/band-pdos-hp.conf .

   phonopy -p -s band-pdos-hp.conf
else
   echo "I don't understand"
endif
