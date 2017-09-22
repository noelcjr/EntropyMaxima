#!/usr/bin/env bash

#TITLE: Making a fusion peptide by joining Insulin and Leucine Zipper with a linker.

# Extending a protein structure by adding a linker to the N-terminal and fusing this
# to another protein's C-terminal.  

mkdir flower_pot
cd flower_pot

c_echo GREEN "1. Download insulin (2HIU) and Leucine zipper (2ZTA)."

script_cmd wget https://files.rcsb.org/view/2HIU.cif
script_cmd wget https://files.rcsb.org/view/2ZTA.cif

c_echo GREEN "2. Check for gaps in the structures."

script_cmd pdb_cif.py gaps --input 2HIU.cif
script_cmd pdb_cif.py gaps --input 2ZTA.cif

#   The insulin structure (2HIU.cif) has 10 models. Models are repeated versions of the
#   same protein inside a structure. If you open a PDB structure in VMD, you will see only one
#   structure at the time. If you click the arrows at the bottom of the Main window in vmd, 
#   the models will be shown one after the other, and it appears as if it is moving. This looks
#   like a simulation, but it is not. 

#   The insuline structure was obtained with Nuclear Magnetic Resonance (NMR),
#   and not by x-ray diffractions of crystalized structures. NMR structures can give multiple
#   snapshots, or models, of the structure to identify regions that move more. NMR techniques
#   are restricted to small peptides. For more info:
#   https://en.wikipedia.org/wiki/Nuclear_magnetic_resonance_spectroscopy_of_proteins

#   In one of the example where structures are completed, we saw that 1BRS.pdb has one structure 
#   with the protein of interest repeated three times in the same model. For 2HIU has multiple 
#   models instead of one. This is not a formating inconsistency because 1BRS formed a crystal
#   latice by the oligomerization of three identical proteins. The 2HIU.pdb structure was obtain
#   by NMR in solution, so the repeated structures are just representative of the range of motions
#   of a single protein, and they are superimposed and not forming a lattice. This is even more
#   evident by the number of chaiins in bothe structures. There are 6 chains in 1BRS, A,B,C,D,E and
#   F and only one model stricture. Insulin (2HIU) has only two chains repeated across 10 models.    

#3. The last two steps found no gaps in the two structures, the next steps considers sequence 
#   information for detecting gaps, and it found two missing residues in Leucine zipper's (2ZTA)
#   C-term.

script_cmd gen_csv.py --fromcif --cif 2HIU.cif --out1 2HIU.csv
script_cmd gen_csv.py --fromcif --cif 2ZTA.cif --out1 2ZTA.csv

#4. The previous steps added N and C terminals to cap the ends of the proteins.

script_cmd del_residue.py --rem "1,1,A,ACE" --inp 2HIU.csv --out 2HIU.csv
script_cmd del_residue.py --rem "1,2,B,ACE" --inp 2HIU.csv --out 2HIU.csv

script_cmd del_residue.py --rem "34,1,A,CTER" --inp 2ZTA.csv --out 2ZTA.csv
script_cmd del_residue.py --rem "34,1,B,CTER" --inp 2ZTA.csv --out 2ZTA.csv

#  Terminals are protein modifications that signal the begining and end of a chain. They have bioological and
#  experimental relevance.

#5. Now let's add a string of amino acids to insulin. These amino acids will work as a linker to Leucine zipper.
#   They are added at the end of insuline chains where the N-terminal had been removed in the previous step, and it
#   also specifies the direction the addition is made. Ndir adds the amino acids in the direction of where the
#   N-terminal would be, and Cdir in the opositite directions. The two possible directions would add amino acids
#   in oposites directions using atoms on the amino acid that is attached to as reference point.

linker1="ASP,ASP,ASP,ASP,LYS"
linker1s="DDDDK"
linker2="ASP,ASP,ASP,ASP,LYS"
linker2s="DDDDK"
script_cmd add_residues.py --apn "1,1,A,Ndir" --res $linker1 --inp 2HIU.csv --out 2HIU.csv
script_cmd add_residues.py --apn "1,2,B,Ndir" --res $linker2 --inp 2HIU.csv --out 2HIU.csv

#6. Clean up. We will work on only one insulin structure, so we will delete some other files to
#   make things neat.

script_cmd rm 2HIU_2.pdb 2HIU_3.pdb 2HIU_4.pdb 2HIU_5.pdb 2HIU_6.pdb 2HIU_7.pdb 2HIU_8.pdb 2HIU_9.pdb 2HIU_10.pdb

#7. For insulin and Leucine Zipper, reduce adds hydrogens to histidines, pdb_cif.py --prepare prepares files for
#   charmm, and then charmm is ran.
script_cmd pdb_cif.py prepare --input 2HIU_1.pdb --terminals "A,none,CTER:B,none,CTER"
script_cmd pdb_cif.py prepare --input 2ZTA_1.pdb --terminals "A,ACE,none:B,ACE,none"

#9. Fix PDBs from Charmm output. For some reason charmm places the chain identifier in a different column that
#   Biopyhton's parser can't detect.

script_cmd pdb_cif.py fixpdb --input 2hiu_1r.pdb
script_cmd pdb_cif.py fixpdb --input 2zta_1r.pdb

# CONCLUSION: In this example we fused 2hiu (insulin) with a linker. The structures were properly prepared
#             by adding missing residues and atoms, and adding terminals. So far this script has not fused the 
#             extended 2hiu (insulin) to 2zta (leucine zipper). The following example, make_flower_and_minimize,
#             the structures prepared in this example (2hiu_1r.pdb and 2zta_1r.pdb) will be joined to be one
#             protein. The joining of these two structures will be done from every possible angle of approach
#             and at different distances. This will create 26 structures of the same joined proteins in differen
#             orientational approaches that when open all at once for vizualization, it resembles a dandalion, or
#             a flower with insulin at the middle and leucine zipper around resembling sepals.  

cd ../..
