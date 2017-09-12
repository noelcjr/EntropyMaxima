#!/bin/bash

yourfile="s_90_315_18_aa_bb"

cp "../Lecture_4/"$yourfile".pdb" .

mkdir $yourfile
cp $yourfile".pdb" $yourfile
cp /home/Programs/EntropyMaxima/charmm_templates/setup_one.inp $yourfile"/setup_one.inp"
cp /home/Programs/EntropyMaxima/charmm_templates/minimize.inp $yourfile"/minimize.inp"
cd $yourfile

perl -pi -e "s/code/home\/noel\/Projects\/Protein_design\/EntropyMaxima/g" setup_one.inp
perl -pi -e "s/first none last none/first ACE last CTER/g" setup_one.inp
perl -pi -e "s/code/home\/noel\/Projects\/Protein_design\/EntropyMaxima/g" minimize.inp



# Warning setup_one.inp has a path that will not be found for par and top
perl -pi -e 's/INFILE/'$yourfile'/g' setup_one.inp
perl -pi -e 's/OUTFILE/'$yourfile'r/g' setup_one.inp

pdb_cif.py prepare --input $yourfile".pdb" --crdout $yourfile".crd" --seqfix yes

charmm < setup_one.inp > setup_one.out

perl -pi -e 's/INPUT/'$yourfile'r/g' minimize.inp
perl -pi -e 's/OUTPUT/'$yourfile'r/g' minimize.inp

charmm < minimize.inp > minimize.out

# Now, copy and paste the pdb files to your computer and load them in VMD

mkdir box_setup
cd box_setup
cp /home/Programs/EntropyMaxima/charmm_templates/waterbox.inp .

perl -pi -e "s/PATH/\/home\/CCL_Lectures\/Lecture_5\/"$yourfile"/g" waterbox.inp
perl -pi -e "s/INFILE/..\/"$yourfile"r/g" waterbox.inp
perl -pi -e "s/OUTFILE/"$yourfile"r_min3_box/g" waterbox.inp

charmm < waterbox.inp > waterbox.out

cd ..
mkdir add_ions
cd add_ions
cp /home/Programs/EntropyMaxima/charmm_templates/add_NaCl.inp .

perl -pi -e "s/PATH/\/home\/CCL_Lectures\/Lecture_5\/"$yourfile"box_setup/g" add_NaCl.inp
perl -pi -e "s/INPUT/"$yourfile"r_min3_box/g" add_NaCl.inp
perl -pi -e "s/OUTPUT/"$yourfile"r_min3_box_ions/g" add_NaCl.inp

#resid 23294 .or. resid 7732 .or. resid 23188 .or. resid 18521 .or. resid 13532 .or. -
#resid 26196 .or. resid 11922 .or. resid 2797 .or. resid 31311 .or. resid 22681 .or. -
#resid 30947 .or. resid 27474 .or. resid 28245 .or. resid 2507 .or. resid 1439 .or. -
#resid 3133 .or. resid 27125 .or. resid 6735 .or. resid 12568 .or. resid 16326 .or. -
#resid 4105 .or. resid 17784 .or. resid 20966 .or. resid 28651 .or. resid 22830 .or. -
#resid 5108 .or. resid 5796 .or. resid 21058 .or. resid 22218 .or. resid 24112 .or. -
#resid 11472 

#resid 9420 .or. resid 27150 .or. resid 16208 .or. resid 29689 .or. resid 11628 .or. -
#resid 16053 .or. resid 16674 .or. resid 29483 .or. resid 15944 .or. resid 9225 .or. -
#resid 18921 .or. resid 21926 .or. resid 4376 .or. resid 12951 .or. resid 11388 .or. -
#resid 389 .or. resid 24358 .or. resid 4103 .or. resid 16723 .or. resid 32512 .or. -
#resid 8475 .or. resid 33687 .or. resid 19381 .or. resid 31323 .or. resid 18149 .or.

charmm < add_NaCl.inp > add_NaCl.out

cd ..
mkdir NAMDsim
cd NAMDsim
cp /home/Programs/EntropyMaxima/charmm_templates/NAMD.conf $yourfile".conf"

perl -pi -e "s/PATH/\/home\/CCL_Lectures\/Lecture_5\/"$yourfile"\/add_ions/g" $yourfile".conf"
perl -pi -e "s/INPUT/"$yourfile"r_min3_box_ions/g" $yourfile".conf"
perl -pi -e "s/OUTPUT/"$yourfile"r_min3_box_ions_namd/g" $yourfile".conf"

pdb_cif.py maxmin --input ../add_ions/s_0_0_0_aa_bbr_min3_box_ions.pdb 
#('Number of Atoms ', 33357)
#('Xmin - Xmax', -54.646999, ' - ', 54.738998)
#('Ymin - Ymax', -30.139999, ' - ', 30.367001)
#('Zmin - Zmax', -25.535, ' - ', 25.527)
#--------------------------------
#('Celbasis X, Y, Z:', 109.386, 60.507, 51.062)

perl -pi -e "s/CELL_X/109.386/g" $yourfile".conf"
perl -pi -e "s/CELL_Y/60.507/g" $yourfile".conf"
perl -pi -e "s/CELL_Z/51.062/g" $yourfile".conf"

cp /home/Programs/NAMD/NAMD_2.11_Linux-x86_64-multicore/namd2 /usr/local/bin/namd2_multicore

namd2_multicore +p8 s_0_0_0_aa_bb.conf &> s_0_0_0_aa_bb_1.out &
