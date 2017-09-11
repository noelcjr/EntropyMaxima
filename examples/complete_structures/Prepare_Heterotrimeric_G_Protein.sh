#!/bin/bash

# TITLE: Completing the heterotrimeric G-protein Alpha subunit, an insight into why some structural
#       information is missing, and how some biologically relevant features must be kept.

# In the other example in this folder, we completed a Barnase-Barstar structure missing structure in
# one dimer from another sequence-identical dimer structure from within the same crystal structure. 

# In this tutorial, you will complete a region of one structure (1GDD.cif) from another structure
# (1GIA.pdb). The two structures are identical in amino acid sequence, but 1GDD.cif has a hydrolyzed
# Guanosine diphosphate (GDP) molecule (i.e. only two phosphate groups), and 1GIA.pdb has a 
# un-hydrolyzed Guanosine triphosphate (GTP) (i.e. three phosphate groups). These two molecules bind 
# to the same protein in an area near a section that that is missing in 1GDD.cif. 

# Why is there missing regions in some parts of a protein and not others. The reason is that even a 
# the low temperatures that proteins are crystallize there is thermal motion that prevents some part
# of proteins to be photographed with x-rays. It is something like taking a picture of a car moving
# fast, and just getting a blurred image with low resolution. Likewise, some parts of the protein 
# move faster than others.

# Completing a second structure gives us a glimpse of the lack of consistency in naming protein
# chains and molecules bound to the protein. This is important to realize because you have 
# to be familiar with the organization of chains in this molecules and their biological relevance
# to model the protein the proper way. You need to read about the function of the protein
# to know what to complete, keep or rename without losing relevant structural information.

# We will also learn a different way to complete missing amino acids when we do not have
# another protein structure to copy and paste the missing structural information.

#1. if you are working on the same container as in Lecture 1, you should still
#   have the wget program and the following command should run:

wget https://files.rcsb.org/view/1GDD.cif
wget https://files.rcsb.org/view/1GIA.cif

#2. Check for gaps test pased for both 1GDD and 1GIA.
   
pdb_cif.py gaps --input 1GDD.cif
pdb_cif.py gaps --input 1GIA.cif

#   This shows that amino acid numbers 202 to 217 are missing in 1GDD and present in 1GIA.

#3. First, we will just align 1GDD.pdb to 1GIA.pdb using atoms present and near the missing
#   region, but we will not cut and paste the missing atoms from one structure to another.
#   Alignments are almost always done on alpha carbons (CA) in the back bone of the crystal
#   structure. Side chains have more mobility than back bones in proteins, and using
#   other atoms can produce bad alignments. 

pdb_cif.py align --input 1GIA.cif --refatoms CA,A,198,201:CA,A,218,221 --fit 1GDD.cif --fitatoms CA,A,198,201:CA,A,218,221 --out 1GDD_to_1GIA_just_fit.pdb

#   This is to show you an added feature of the align command in pdb_cif.py. The alignment
#   is done, and a Root mean square deviation (RMSD) of the atoms used for alignment is given.
#   This measure of closeness in alignment is important to analyse simulations later.

#4. We now do the same alignment, but this time we add the missing region. When adding regions,
#   every atom from the amino acid is added to the structure with the gap. This could potentially
#   cause unwanted clashes between amino acids. We will deal with those later, for now it is just
#   important to make you aware that added regions in one protein might no be able to fit well in
#   the pocket of the structure with missing amino acids.  

pdb_cif.py align --input 1GIA.cif --refatoms CA,A,198,201:CA,A,218,221 --fit 1GDD.cif --fitatoms CA,A,198,201:CA,A,218,221 --out 1GDD_to_1GIA_fit_and_add.pdb --addatoms A,202,217:A,202,217

#   Note: If you open the structures aligned by the last command, there is a slight misalignment
#   in placing one region of one protein on another because is not a perfect fit. The added region 
#   will not show up if you use ribbons to display the backbone in VMD, but the atoms are there. 
#   This is to remind you that the alignments are rigid and do not always fit perfectly, and VMD 
#   picks this slight misalignment in atoms structures, and chooses not to draw ribbons. But the
#   alignments worked, and the added atoms will show up if rendered as atoms, lines or bonds
#   rendering modes.

#5. Something to think about. If you open 1GIA.pdb and 1GDD.pdb with your favorite text editor 
#   and go to where the GDP and GTP are in the file, you will notice that while these two molecules 
#   are surrounded by the protein, they are a separate molecule, yet they have the same chain 
#   identifier. 1GIA even has a additional Magnesium ion with the same chain identifier. This 
#   could lead to problems down the road. For now, just keep in mind that while the formatting 
#   of PDBs where improved in the CIF files, certain inconsistencies still are present because of 
#   decisions made by crystalographers. This only highlights the important of fully understanding
#   the biological relevance of the protein we are modeling and what we want to keep, separate, align
#   or complete. In this case, we would want to make GDP, GTP and Mg ion to a different chain
#   identifier. This can be done in a text editor.  

#6. Now, There is another way to complete a structure when there is no 'twin' structure
#   from which we can get missing amino acids parts. The following scripts will detect
#   gaps, including those in the C-terminal, and will add amino acids in a straight line.
#   In this case the missing parts extend out into the solvent, and not towards the
#   interior of the protein. This could be a problem for modeling other proteins with this
#   method. This method takes a CIF file, and it outputs the structure in a comma separated 
#   file as well as a PDB file for each model in the structure.

gen_csv.py --fromcif --cif 1GDD.cif --out1 1GDD_0.csv

#   CONLCUSION: In this lecture we saw that two identical proteins have a region that is
#                missing when the protein is bound to GDP, and present when bound to GTP.
#   This is crystallographic evidence of the hydrolization process of GTPase proteins. The
#   missing region suggests that there is a role by the dynamics of the protein in biological
#   function. This is at the heart of why we want to model proteins. In future lectures
#   we will learn to do molecular dynamics simulations and understand energies. Simulations
#   are models that enhance what is observed in crystal structures, and they could be used,
#   for example, in understanding the roles of the amino acids involved in the hydrolization
#   of GTP so that we can do computational alanine scanning mutagenesis, and measure what
#   amino acids are more important to this function. This would reduce the amount of work
#   in the laboratory when doing alanine scanning mutagenesis by guiding the process in a
#   more informed way. 

#   We also learn that we cannot trust crystal structures in their naming conventions of 
#   chains, or even their numbering. This should easily be solved by a careful modeller of
#   proteins, but it requires a lot of reading to be sure that what we model is sound.
