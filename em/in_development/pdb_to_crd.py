#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 10:47:42 2018

@author: noel
"""
import em.code.input_output as IO
import Bio.PDB as struct

dir_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1B27/brs_A41C_A83C/structures/' 
#'/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1B3S/brn_H102A_brs_Y29F_C41A_C83A/structure/'
for i in range(0,990):
    structure = struct.PDBParser(QUIET=True).get_structure('pdb', dir_path+str(i)+'.pdb')
    line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
    
    count_atoms = 1
    count_resid = 1
    lines = []
    for j in structure.get_atoms():
        if count_atoms == 1:
            current_resid = j.get_parent().id[1]
        if current_resid != j.get_parent().id[1]:
            current_resid = j.get_parent().id[1]
            count_resid += 1
        lines.append(line.format(count_atoms, count_resid, j.get_parent().resname, j.get_full_id()[4][0], j.get_coord()[0], 
              j.get_coord()[1], j.get_coord()[2], j.get_full_id()[2], j.get_parent().id[1], 0.00))
        count_atoms += 1
    count_atoms -= 1    
    len(str(count_atoms))
    lines.insert(0,(5-len(str(count_atoms)))*' '+str(count_atoms))
    lines.insert(0,'*')
    lines.insert(0,'* Date: XX/XX/XXXX by Noel Carrascal')
    lines.insert(0,'* Converting PDB to CRD')
    
    outFile = open(dir_path+str(i)+'.crd', 'w')
    for j in lines:
        outFile.write(j+'\n')
    outFile.close()


