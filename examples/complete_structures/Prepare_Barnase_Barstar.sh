#!/bin/bash

# Title: Extracting structure's chains from within a CIF or PDB file

# Requirements: 1. Install Entropy Maxima.

#  Introduction:
#  This script downloads 1BRS.cif from the protein data-bank.
#  1BRS has 6 chains, A,B,C,D,E and F. AD, BE and CF form dimers.
#  It is a crystal structure with three dimers, and we only need one.
#  How do I know that the crystal structure has three separate dimers 
#  of the same protein instead of one oligomer formed by three pairs of 
#  proteins? The simplest way to find out is by reading the crystalographic 
#  paper that is in the CIF headers with meta data information.
   
#  Chains belonging to the same dimers are not contiguous to each other
#  in alphatical order by chain identifiers in this crystal structure, but it is
#  a common practice in other crystal structures. You would expect to have A and B
#  chains forming a dimer, as well as CD and EF, but for unknown reasons the 
#  crystalographers chose to label chains differently. This could only be observed 
#  after opening the structure on VMD outside the container.
#  From this observation, we write the first command. 

#1. Get barnase-barstar crystal structure from rcsb.org
mkdir Barnase_Barstar
cd Barnase_Barstar 
wget https://files.rcsb.org/view/1BRS.cif

#1. Extract chains corresponding to dimers in the PDB file. For most protein
#   modeling exercises, we only need one dimer that we need to extract, or separate,
#   from the other two dimers. The following command separates chains from structures
#   by AD,BE,CF and generates outputs in PDB format.

pdb_cif.py extract --input 1BRS.cif --chains --groups AD,BE,CF

#   As always, confirm that what you did this right by copying the generated PDB 
#   files outside the container for visualization with VMD. We will pick one of
#   the three dimers for more analysis.

#2. It is difficult to detect all the gaps in structure by visual inspection.
#   Gaps in the middle of the sequence are easy to miss, and gaps at the
#   begining or end of each chain of the proteins are even harder to detect visually.
#   The following command does a gap analysis that gives the amino acids of the protein
#   that are missing coordinates. It detecks gaps at N-terminal if structure begins
#   at a number higher than 1, in the middle if residue numbers are not present, but
#   it will not detect gaps at the C-terminal. That information could only be 
#   detected from sequence information that the gap command ignores.

pdb_cif.py gaps --input 1BRS_0_AD.pdb
pdb_cif.py gaps --input 1BRS_0_BE.pdb
pdb_cif.py gaps --input 1BRS_0_CF.pdb

#   The output gave a Model number and a Chain identifier. SOme PDBs have more than
#   one structure of the same protein, and they are usally number starting at 0.
 
#   The three numbers inside the parenthesis correpond to the gap number and the
#   amino acid number ranges. For example: (1,2,6) is read as gap 1 that corresponds
#   to missing amino acids 2,3,4,5,6. When there are no gaps, the output is (0,0,0).

#3. It is now possible to complete missing regions in some dimers, from regions
#   present in other dimers. Completing resideus 1 and 2 missing in chain C in one 
#   dimer from chain B of another dimer using the following command.

pdb_cif.py align --input 1BRS_0_BE.pdb --refatoms CA,B,3,6 --fit 1BRS_0_CF.pdb --fitatoms CA,C,3,6 --out 1BRS_complete_CD_1.pdb --addatoms B,1,2:C,1,2

#   The --addatoms section is optional. Excluding this option allows alignments 
#   that are required for comparing structural features without changing them. 
#   See the Preparing_Heterotrimeric_G_Protein.sh for an example of this.

#4. Copy 1BRS_0_BE.pdb, 1BRS_0_CF.pdb and 1BRS_complete_CD_1.pdb to your operating system 
#   and open them with VMD to check that atoms were added from one structure to
#   the other after aligning the structures.

cd ..
#CONCLUSION: It is important to know what is inside a CIF structure file to find out
#            if there are mutiple repeated structures of the same protein. We then can
#extract only one dimer from the original crystal structure and use the other dimers
#`to fill gaps in the selected dimer. The selected dimer could be prepared for simulation.
