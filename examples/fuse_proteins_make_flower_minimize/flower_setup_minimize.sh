#!/bin/bash

# Example run:
#bash setup_flowers_3.sh DDDDK_DDDDK 2hiu_a_ddddk_b_ddddk_1rr.pdb 1brs_comp_1rr_1rr.pdb "A:F,B:C" 45 5 45 5

# This script creates a flower with a center (stigma) structure and a rotated (petal) structure.
# This scripts works for linking two chains between two protein structure. A more general script for linking between trimers or
# structures with higher number of monomers is not programmed yet.
if [ -z "$10" ];
then
   echo "ERROR: This script takes ten parameters. Follow usage below."
   echo "usage: $0 Folder insulin_pdb_path lue_zip_pdb_path A:A,B:B 45 30 45 5 terminals regions"
   exit
else
   script_cmd cd flower_pot
   dir=$1
   if [ ! -d $dir ];
   then
      #rm -f -r $dir It is too risky to remove. It might delete a lot of information. Remove manully if needed. 
      echo "Warning: Folder for linker must exist. Script will exit without linking the two PDBs and creating a flower."
      exit
   fi
   script_cmd cd $dir
   # Relative from where the script is run.
   relative_path="../../"
   centerPDB=$2
   rotatePDB=$3
   chn_lnk1=$4
   angle=$5
   initdist=$6
   distance=$7
   interval=$8
   
   flwr1=${chn_lnk1//[:]/_}
   flwr1=${flwr1//[,]/}
   rnm_chn1=${flwr1: 1:1}
   rnm_chn2=${flwr1: 4:1}
   centerHeader=`basename ${centerPDB%.*b}`
   rotateHeader=`basename ${rotatePDB%.*b}`
   for (( i=initdist; i <= distance; i+=interval ))
   do
      current_flwr=$centerHeader"_"$rotateHeader"_"$flwr1"_"$i
      current_flwr=${current_flwr,,}
      if [ -d $current_flwr ];
      then
         echo "Warning: Folder for center structure already exist. Script will exit without any changes to the folder to avoid losing data. Delete Manually if needed."
         exit
      fi
      script_cmd mkdir $current_flwr
      script_cmd cd $current_flwr
      script_cmd flower.py --center $relative_path$centerPDB --rotate $relative_path$rotatePDB --angle $angle --distance $i --map yes --link $chn_lnk1
      for file in $(ls s_*);
      do  
        folder=${file%.*b}
        folder=${folder,,}
        script_cmd mkdir $folder
        script_cmd grep ATOM $file > $folder"/"${file,,}
        script_cmd rm $file
        script_cmd cd $folder
        
        script_cmd pdb_cif.py prepare --input $file --terminals $9
        script_cmd minimization.py --input $folder"r.crd" --info ${10}

        script_cmd cd ..
        script_cmd mv $folder "d_"$folder
        exit
      cd ..
   done
   cd ..
fi
