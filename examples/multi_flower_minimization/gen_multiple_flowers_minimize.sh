#!/usr/bin/env bash

#TITLE: Making a fusion peptide by joining Insulin and Leucine Zipper with a linker.

# REQUIREMENTS: Entropy and VMD. A VMD state for visualization is created in the last
#               two lines of this script. This VMD file needs to be copied outside the
#               conteiner to be opened in VMD.
               
# Extending a protein structure by adding a linker to the N-terminal and fusing this
# to another protein.  

linker1="DDDDK"
linker2="DDDDK"

mkdir $linker1"_"$linker2
cd $linker1"_"$linker2

mkdir "2hiu_2zta_AA_"$linker1s"_BB_"$linker2s
cd "2hiu_2zta_AA_"$linker1s"_BB_"$linker2s
script_cmd flower.py --center ../2hiu_1r.pdb --rotate ../2zta_1r.pdb --angle 45 --distance 45 --map yes --link "A:A,B:B"


cd ../..
