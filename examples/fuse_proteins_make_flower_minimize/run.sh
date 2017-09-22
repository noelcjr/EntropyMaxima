#!/bin/bash

# WARNING: To run this program, you must have run the fuse protein exmples as it needs
#          to files inside the flower_pot/ directory for input here

if [ -d `DDDDK_DDDDK`];
then
   #rm -f -r $dir It is too risky to remove. It might delete a lot of information. Remove manully if needed. 
   echo "Warning: Folder for flowers already exist becasue it was probably ran already."
   echo "         This is to avoid running the script on top of already ran data which"
   echo "         will result in unknow behavior/result, or at least waste resources."
   echo "         If you still want to run the script again, delete DDDDK_DDDDK at your"
   echo "         own risk and run it."
   exit
fi

mkdir DDDDK_DDDDK
insulin_path="../add_linker_prepare_proteins_for_fusion/setup/2hiu_1r.pdb"
lue_zip_path="../add_linker_prepare_proteins_for_fusion/setup/2zta_1r.pdb" 
terminals="A,ACE,CTER:B,ACE,CTER"
regions="A,1-34,35-39,40-60:B,1-34,35-39,40-69"

time ./flower_setup_minimize.sh DDDDK_DDDDK $insulin_path $lue_zip_path "A:A,B:B" 45 30 45 5 $terminals $regions &> DDDDK_DDDDK/2hiu_1r_2zta_1r_AA_BB_45_30_45_5.out &
time ./flower_setup_minimize.sh DDDDK_DDDDK $insulin_path $lue_zip_path "A:B,B:A" 45 30 45 5 $terminals $regions &> DDDDK_DDDDK/2hiu_1r_2zta_1r_AB_BA_45_30_45_5.out &
time ./flower_setup_minimize.sh DDDDK_DDDDK $lue_zip_path $insulin_path "A:A,B:B" 45 30 45 5 $terminals $regions &> DDDDK_DDDDK/2zta_1r_2hiu_1r_AA_BB_45_30_45_5.out &
time ./flower_setup_minimize.sh DDDDK_DDDDK $lue_zip_path $insulin_path "A:B,B:A" 45 30 45 5 $terminals $regions &> DDDDK_DDDDK/2zta_1r_2hiu_1r_AB_BA_45_30_45_5.out &
