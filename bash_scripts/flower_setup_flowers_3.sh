#!/bin/bash

# Example run:
#bash setup_flowers_3.sh DDDDK_DDDDK 2hiu_a_ddddk_b_ddddk_1rr.pdb 1brs_comp_1rr_1rr.pdb "A:F,B:C" 45 5 45 5

# This script creates a flower with a center (stigma) structure and a rotated (petal) structure.
# This scripts works for linking two chains between two protein structure. A more general script for linking between trimers or
# structures with higher number of monomers is not programmed yet.
if [ -z "$1" ];
then
   echo usage: $0 Linker_Directory PDB_center_structure PDB_rotate_structure first_link_chains distance angle
   exit
fi
if [ -z "$2" ];
then
   echo usage: $0 Linker_Directory PDB_center_structure PDB_rotate_structure first_link_chains distance angle
   exit
fi
if [ -z "$3" ];
then
   echo usage: $0 Linker_Directory PDB_center_structure PDB_rotate_structure first_link_chains distance angle
   exit
fi
if [ -z "$4" ];
then
   echo usage: $0 Linker_Directory PDB_center_structure PDB_rotate_structure first_link_chains distance angle
   exit
fi
if [ -z "$5" ];
then
   echo usage: $0 Linker_Directory PDB_center_structure PDB_rotate_structure first_link_chains distance angle
   exit
fi
if [ -z "$6" ];
then
   echo usage: $0 Linker_Directory PDB_center_structure PDB_rotate_structure first_link_chains distance angle
   exit
else
   dir=$1
   if [ ! -d $dir ];
   then
      #rm -f -r $dir It is too risky to remove. It might delete a lot of information. Remove manully if needed. 
      echo "Warning: Folder for linker must exist. Script will exit without linking the two PDBs and creating a flower."
      exit
   fi
   cd $dir
   centerPDB=$2
   rotatePDB=$3
   chn_lnk1=$4
   angle=$5
   initdist=$6
   distance=$7
   interval=$8
   
   flwr1=${chn_lnk1//[,]/_}
   flwr1=${flwr1//[:]/}
   rnm_chn1=${flwr1: 1:1}
   rnm_chn2=${flwr1: 4:1}
   for (( i=initdist; i <= distance; i+=interval ))
   do
      if [ -d $flwr1"_"$i ];
      then
         echo "Warning: Folder for center structure already exist. Script will exit without any changes to the folder to avoid losing data. Delete Manually if needed."
         exit
      fi
      flwr1=${flwr1,,}
      mkdir $flwr1"_"$i
      cd $flwr1"_"$i
      flower.py" --center "../"$centerPDB --rotate "../"$rotatePDB --angle $angle --distance $i --map yes --link $chn_lnk1
    
      for file in $(ls s_*);
      do  
        folder=${file%.*b}
        folder=${folder,,}
        mkdir $folder
        #converting file to lower case in new directory
        grep ATOM $file > $folder"/"${file,,}
        #cat $file > $folder"/"$file
        rm $file
        cd $folder
        cp $charmm"setup_one.inp" .
        cp $charmm"min_enr_mmgbsa_rmsd_flower.inp" .
        reduce -HIS -FLIP -OH -ROTEXOH -BUILD -OCC0.0 -H2OOCC0.0 -H2OB1000 $folder".pdb" > $folder"r.pdb"
        pdb.py --prepare --pdbin1 $folder"r.pdb" --crdout $folder"r.crd" --seqfix yes
        perl -pi -e 's/generate A first none last none setup/generate '$rnm_chn1' first ACE last CTER setup/g' setup_one.inp
        perl -pi -e 's/generate B first none last none setup/generate '$rnm_chn2' first ACE last CTER setup/g' setup_one.inp
        perl -pi -e 's/A.SEQ/'$rnm_chn1'.SEQ/g' setup_one.inp
        perl -pi -e 's/B.SEQ/'$rnm_chn2'.SEQ/g' setup_one.inp
        perl -pi -e 's/A_FIXRES.INP/'$rnm_chn1'_FIXRES.INP/g' setup_one.inp
        perl -pi -e 's/B_FIXRES.INP/'$rnm_chn2'_FIXRES.INP/g' setup_one.inp
        perl -pi -e 's/INFILE/'$folder'r/g' setup_one.inp
        perl -pi -e 's/OUTFILE/'$folder'rr/g' setup_one.inp
        charmm_xxlarge < setup_one.inp > setup_one.out
        perl -pi -e 's/INPUT/'$folder'rr/g' min_enr_mmgbsa_rmsd_flower.inp
        perl -pi -e 's/OUTPUT/'$folder'rr/g' min_enr_mmgbsa_rmsd_flower.inp
        perl -pi -e 's/FILEID/'$folder'/g' min_enr_mmgbsa_rmsd_flower.inp
        # The remaining substitutions in this loop should be read from the command prompt because it could be different for other systems.
        perl -pi -e 's/segid A/segid C/g' min_enr_mmgbsa_rmsd_flower.inp
        perl -pi -e 's/segid B/segid F/g' min_enr_mmgbsa_rmsd_flower.inp
        rotate_a="1:110"
        rotate_b="1:89"
        linker_a="111:114"
        linker_b="90:93"
        center_a="115:135"
        center_b="94:123"
        perl -pi -e 's/ROTA/'$rotate_a'/g' min_enr_mmgbsa_rmsd_flower.inp
        perl -pi -e 's/ROTB/'$rotate_b'/g' min_enr_mmgbsa_rmsd_flower.inp
        perl -pi -e 's/LIKA/'$linker_a'/g' min_enr_mmgbsa_rmsd_flower.inp
        perl -pi -e 's/LIKB/'$linker_b'/g' min_enr_mmgbsa_rmsd_flower.inp
        perl -pi -e 's/CENA/'$center_a'/g' min_enr_mmgbsa_rmsd_flower.inp
        perl -pi -e 's/CENB/'$center_b'/g' min_enr_mmgbsa_rmsd_flower.inp
        charmm_xxlarge < min_enr_mmgbsa_rmsd_flower.inp > min_enr_mmgbsa_rmsd_flower.out
        rm fort*
        rm *.inp
        rm *.out
        rm *.SEQ
        rm *.INP
        cd ..
        mv $folder "d_"$folder
      done
      cd ..
   done
fi
