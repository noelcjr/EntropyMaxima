#!/bin/bash

yourfile="s_90_315_18_aa_bb"
mypath=$(echo "`pwd`" | sed 's/\//\\\//g')

echo "WARNING: This script will run a molecular dynamics simulations. It will take"
echo "         about 26 minutes using only one core of your computer. If you have more"
echo "         than one core, and you want to run this script faster,"
echo "         enter a number for the number of cores you want to use:"

read cores

echo The NAMD simulatiosn will be run using $cores cores.

# ...
# call it
pause 'Press [Enter] key to continue...'

mkdir $yourfile
cp $yourfile".pdb" $yourfile
cp /home/Programs/EntropyMaxima/charmm_templates/setup_one.inp $yourfile"/setup_one.inp"
cp /home/Programs/EntropyMaxima/charmm_templates/minimize.inp $yourfile"/minimize.inp"
cd $yourfile

perl -pi -e "s/code/home\/Programs\/EntropyMaxima/g" setup_one.inp
perl -pi -e "s/first none last none/first ACE last CTER/g" setup_one.inp
perl -pi -e "s/code/home\/Programs\/EntropyMaxima/g" minimize.inp

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

perl -pi -e "s/code/home\/Programs\/EntropyMaxima/g" waterbox.inp
perl -pi -e "s/PATH/"$mypath"\/"$yourfile"/g" waterbox.inp
perl -pi -e "s/INFILE/"$yourfile"r/g" waterbox.inp
perl -pi -e "s/OUTFILE/"$yourfile"r_min3_box/g" waterbox.inp

charmm < waterbox.inp > waterbox.out

cd ..
mkdir add_ions
cd add_ions
cp /home/Programs/EntropyMaxima/charmm_templates/add_NaCl.inp .

perl -pi -e "s/code/home\/Programs\/EntropyMaxima/g" add_NaCl.inp

perl -pi -e "s/PATH/"$mypath"/g" add_NaCl.inp
perl -pi -e "s/INPUT/"$yourfile"\/box_setup\/"$yourfile"r_min3_box/g" add_NaCl.inp
perl -pi -e "s/OUTPUT/"$yourfile"r_min3_box_ions/g" add_NaCl.inp

perl -pi -e "s/SOD/resid 7396 .or. resid 3578 .or. resid 9121 .or. resid 4861 .or. resid 8441 .or. -\nresid 7003 .or. resid 9138 .or. resid 3376 .or. resid 7666 .or. resid 56 .or. -\nresid 1875 .or. resid 2241 .or. resid 1947 .or. resid 2124 .or. resid 4390 .or. -\nresid 6670 .or. resid 3750 .or. resid 3084 .or. resid 3987 .or. resid 8163 .or. -\nresid 2099 .or. resid 6717 .or. resid 7759 .or. resid 7148 .or. resid 3498 .or. -\nresid 1334 .or. resid 967/g" add_NaCl.inp

perl -pi -e "s/CLA/resid 8466 .or. resid 4115 .or. resid 8082 .or. resid 8530 .or. resid 2587 .or. -\nresid 7260 .or. resid 9543 .or. resid 8080 .or. resid 7897 .or. resid 5145 .or. -\nresid 5966 .or. resid 6796 .or. resid 5363 .or. resid 7311 .or. resid 8720 .or. -\nresid 2715 .or. resid 1375 .or. resid 7256 .or. resid 4421 .or. resid 3185 .or. -\nresid 1943/g" add_NaCl.inp


charmm < add_NaCl.inp > add_NaCl.out

cd ..
mkdir NAMDsim
cd NAMDsim
cp /home/Programs/EntropyMaxima/charmm_templates/NAMD.conf $yourfile".conf"

perl -pi -e "s/code/home\/Programs\/EntropyMaxima/g" $yourfile".conf"
perl -pi -e "s/PATH/"$mypath"/g" $yourfile".conf"
perl -pi -e "s/INPUT/"$yourfile"\/add_ions\/"$yourfile"r_min3_box_ions/g" $yourfile".conf"
perl -pi -e "s/OUTPUT/"$yourfile"r_min3_box_ions_namd/g" $yourfile".conf"

pdb_cif.py maxmin --input ../add_ions/"$yourfile"r_min3_box_ions.pdb 
#('Results from ', '../add_ions/s_90_315_18_aa_bbr_min3_box_ions.pdb')
#('Number of Atoms ', 29567)
#('Xmin - Xmax', -53.646999, ' - ', 53.597)
#('Ymin - Ymax', -30.639999, ' - ', 30.701)
#('Zmin - Zmax', -22.627001, ' - ', 22.811001)
#--------------------------------
#('Celbasis X, Y, Z:', 107.244, 61.341, 45.438004)

perl -pi -e "s/CELL_X/107.244/g" $yourfile".conf"
perl -pi -e "s/CELL_Y/61.341/g" $yourfile".conf"
perl -pi -e "s/CELL_Z/45.438/g" $yourfile".conf"

cp /home/Programs/NAMD/NAMD_2.11_Linux-x86_64-multicore/namd2 /usr/local/bin/namd2_multicore
cp ../add_ions/s_90_315_18_aa_bbr_min3_box_ions.psf .
cp ../../NAMD_simulation.vmd .

echo "After running the simulation, copy the NAMDsim/ folder to your computer and open the "
echo "NAMD_simulation.vmd file using vmd by opening vmd and go to main -> Load Visualisation state, and click on NAMD_simulation.vmd"

# I do not comment out commands without good reason.
namd2_multicore +p$cores s_90_315_18_aa_bb.conf &> s_90_315_18_aa_bb.out &
