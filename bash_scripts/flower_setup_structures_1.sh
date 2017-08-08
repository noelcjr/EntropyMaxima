#!/bin/bash

ex="/home/noel/Projects/Protein_design/EntropyMaxima/src/"
py="/usr/bin/python"
charmm="/home/noel/Projects/Protein_design/EntropyMaxima/charmm_templates/"

# We will get structures that have already been setup and prepared.
# Insulin 2hiu has 10 models from NMR. We will pick the first one only.
# in /home/noel/Projects/Protein_design/EntropyMaxima_Tests_Examples/examples/Linker_minimization
# there is an example of inuline where all 10 models of 2hiu are considered.
# I found this not necessary because insulin is going to be modified to link to Barnase-barstar
# and when linkers are added and minimized, the original structure is so changed that is
# not helpful to consider the other 9 different conformations. Differences in the conformations
# are mostly from fluctuations in the N anc C termials of the A and B chains.

#mkdir 1brs_barnase_barstar
#cp /home/noel/Projects/Protein_design/Barnase_barstar/Struct_Prep/reduce_prep_charmm/1brs_comp_1rr* 1brs_barnase_barstar/
cd 1brs_barnase_barstar 
$py $ex"gen_csv.py" --frompsfcrd --crd 1brs_comp_1rr.crd --psf 1brs_comp_1rr_xplo.psf
$py $ex"del_residue.py" --rem "1,1,C,ACE" --inp 1brs_comp_1rr.csv --out 1brs_comp_1rr.csv
$py $ex"del_residue.py" --rem "1,2,F,ACE" --inp 1brs_comp_1rr.csv --out 1brs_comp_1rr.csv
$py $ex"del_residue.py" --rem "110,1,C,CTER" --inp 1brs_comp_1rr.csv --out 1brs_comp_1rr.csv
$py $ex"del_residue.py" --rem "89,2,F,CTER" --inp 1brs_comp_1rr.csv --out 1brs_comp_1rr.csv
cd ..

#mkdir 2hiu_insulin_no_disulfides
#cp ../../Insulin/struct_prep/2hiu/init_setup/2hiu_1rr_no_disulf* 2hiu_insulin_no_disulfides
cd 2hiu_insulin_no_disulfides
$py $ex"gen_csv.py" --frompsfcrd --crd 2hiu_1rr_no_disulf.crd --psf 2hiu_1rr_no_disulf_xplo.psf 
$py $ex"del_residue.py" --rem "1,1,A,ACE" --inp 2hiu_1rr_no_disulf.csv --out 2hiu_1rr_no_disulf.csv
$py $ex"del_residue.py" --rem "1,2,B,ACE" --inp 2hiu_1rr_no_disulf.csv --out 2hiu_1rr_no_disulf.csv
$py $ex"del_residue.py" --rem "21,1,A,CTER" --inp 2hiu_1rr_no_disulf.csv --out 2hiu_1rr_no_disulf.csv
$py $ex"del_residue.py" --rem "30,2,B,CTER" --inp 2hiu_1rr_no_disulf.csv --out 2hiu_1rr_no_disulf.csv
cd ..

#mkdir 2hiu_insulin_disulfides
#cp ../../Insulin/struct_prep/2hiu/init_setup/2hiu_1rr.* 2hiu_insulin_disulfides
#cp ../../Insulin/struct_prep/2hiu/init_setup/2hiu_1rr_xplo.psf 2hiu_insulin_disulfides
cd 2hiu_insulin_disulfides
$py $ex"gen_csv.py" --frompsfcrd --crd 2hiu_1rr.crd --psf 2hiu_1rr_xplo.psf
$py $ex"del_residue.py" --rem "1,1,A,ACE" --inp 2hiu_1rr.csv --out 2hiu_1rr.csv
$py $ex"del_residue.py" --rem "1,2,B,ACE" --inp 2hiu_1rr.csv --out 2hiu_1rr.csv
$py $ex"del_residue.py" --rem "21,1,A,CTER" --inp 2hiu_1rr.csv --out 2hiu_1rr.csv
$py $ex"del_residue.py" --rem "30,2,B,CTER" --inp 2hiu_1rr.csv --out 2hiu_1rr.csv
cd ..
