#!/usr/bin/python

"""
Created on Sat Apr  9 12:52:52 2016

@author: noel
"""
import os
import sys
# For debugging in spyder
sys_path = '/home/noel/Projects/Protein_design/EntropyMaxima/src'
charmmdir = '/home/noel/Projects/Protein_design/EntropyMaxima/params/charmm27.ff/'
is_in_path = False
for i in sys.path:
    if i == sys_path:
        is_in_path = True
if not is_in_path:
    sys.path.append(sys_path)
import numpy as np
import CHARMM_Parser as CP
import Molecular_Rigid_Manipulations as MRM
import Molecular_Descriptors as md
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import *
import Bio.PDB.Superimposer as sip
#####        INPUT PARAMETERS    ####################################
#outer = "/home/noel/Projects/Protein_design/Insulin/Struct_Prep/pdbs/init_setup_bridge/out/2zta_3.pdb"
center = "/home/noel/Projects/Protein_design/EntropyMaxima/examples/Super_Structure/2hiu_2.pdb"
linker = "/home/noel/Projects/Protein_design/EntropyMaxima/charmm_templates/peptides.pdb"
####################################################################
'''ALL LINKING is done from NTERM of centered protein to CTERM of outer-rotated protein only.
   Linker is generated with the CTERM as the highest residue number.
   Linker must have GLYs at each end that are not part of the linker, are used for
   alignments and need to be deleted after connecting the structure.
   center and outer structures have capping atoms that need to be deleted.
   Linker CTERM residue is aligned to NTERM of linkers using GLs.'''
######  INPUTS REQUIRED #####################################
# List description. 
# [[Center info],[Outer info],[Linker info]]
# Inner list description
# Description [[PATH,chain,residue number],......]
bridge1 = [[center,'A',1],[outer,'A',31],[linker,'C',24]]
bridge2 = [[center,'B',1],[outer,'B',31],[linker,'C',24]]
outer_displacement = [45,0,0]
## NEEDS TO PROVIDE A WAY TO GIVE VECTOR ORIENTATIONS OF center AND outer
## AND DISPLAY IT FOR GUIDANCE IN CHOOSING THIS VALUES
## This option should give, for example:
##    1. orientations for outer and center closest by distance between two termini. 
##    2. orientations for outer and ceter closest by average distance between all termini.
##    3. orientation for outer only by one termini
## FOR NOW THIS ORIENTATIONS ARE SET TO DEFAULTS:
outer_orientation = [-1,0,0]
center_orientation = [0,0,0]
#############################################################
#rig = MRM.Molecular_Rigid_Manipulation(charmmdir)
pdb_parser = PDBParser()
params = CP.read_charmm_FF(charmmdir)
mnp = MRM.Molecular_Rigid_Manipulation(params)
cmc = md.CenterOfMassCalculator()
#############################################################
cen = pdb_parser.get_structure('Centered', center)
mnp.fix_pdb_from_CHARMM(center,17,72)
'''Puts centered molecule's center of masss at (0,0,0)'''
mnp.translate_molecule(cen,mnp.center_molecule(cmc.calculate_center_of_mass(cen)))
out = pdb_parser.get_structure('Rotated', outer)
mnp.fix_pdb_from_CHARMM(outer,17,72)
'''Puts outer molecule's center of masss at (0,0,0)'''     
mnp.translate_molecule(out,mnp.center_molecule(cmc.calculate_center_of_mass(out)))
ccc = md.ChargeCalculator(charmmdir)
'''The center of charge is not a point, but a vector. It indicates the direction
   or dipole moment of the molecule'''
cm = ccc.calculate_center_of_charge(out)
cm = cm/np.linalg.norm(cm)
RM = mnp.alignVectors([1.0,0.0,0.0],cm)
for j in out.get_atoms():
    v2 = [j.get_coord()[0],j.get_coord()[1],j.get_coord()[2]]
    jj = np.dot(v2,RM)
    j.set_coord(jj) 
# FLIP Molecule along y-axis while molecule at centered by center of mass
m = rotaxis(np.pi, Vector(0, 1, 0))
for j2 in out.get_atoms():
    v2 = Vector([j2.get_coord()[0],j2.get_coord()[1],j2.get_coord()[2]]) 
    v3 = v2.left_multiply(m)
    j2.set_coord(v3.get_array())
mnp.translate_molecule(out,outer_displacement)

lin = pdb_parser.get_structure('Linker', linker)
mnp.fix_pdb_from_CHARMM(linker,17,72) 
# Aligns linker and centreal protein by ref(cen) and mob(lin) atoms
ref_atom = [cen[0]['A'][1]['N'],cen[0]['A'][1]['CA'],cen[0]['A'][1]['C']]
mob_atom = [lin[0]['C'][24]['N'],lin[0]['C'][24]['CA'],lin[0]['C'][24]['C']]
super_imposer = Superimposer()
super_imposer.set_atoms(ref_atom, mob_atom)
super_imposer.apply(lin.get_atoms())
''' Fix this to merge all three structures onto one.
for i in lin.get_residues():
    cen.add(i)
    
for i in out.get_residues():
    cen.add(i)
'''
io = PDBIO()
io.set_structure(cen)
io.save('/home/noel/Projects/Protein_design/Insulin/Struct_Prep/pdbs/cen_out.pdb')

io = PDBIO()
io.set_structure(lin)
io.save('/home/noel/Projects/Protein_design/Insulin/Struct_Prep/pdbs/lin_out.pdb')
        
io = PDBIO()
io.set_structure(out)
io.save('/home/noel/Projects/Protein_design/Insulin/Struct_Prep/pdbs/out_out.pdb')
