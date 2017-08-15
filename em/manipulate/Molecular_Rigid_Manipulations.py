# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 19:27:04 2016
@author: Noel Carrascal
@email noelcarrascal@gmail.com
"""
import numpy as np
import Bio.PDB.Vector as vec
import CHARMM_Parser as CP
from Bio.PDB import *
# FIXIT: this function is out Molecular_Rigid_manipulation
# because it does not need charmdir to run and Molecular_rigid_manipulation does.
# Do we need to load charmmdir for Molecular_rigid_manipulation in the constructor?
class Histidin_Correction(object):
    # HSE: HE1, HE2,    , HD2
    # HSD: HE1,    , HD1, HD2
    # HSP: HE1, HE2, HD1, HD2
    def rename_HIS(self, resid):
        HE1 = HE2 = HD1 = HD2 = False
        for i in resid.get_atom():
            if i.id == 'HE1':
                HE1 = True
            elif i.id == 'HE2':
                HE2 = True
            elif i.id == 'HD1':
                HD1 = True
            elif i.id == 'HD2':
                HD2 = True
        
        if HE1 and HD2:
            if HE2 and HD1:
                return 'HSP'
            else:
                if HE2:
                    return 'HSE'
                elif HD1:
                    return 'HSD'
                else:
                    print('WARNNING: Unusual number of hydrogen in Histidine. HSD is default.')
                    return 'HSD'
        else:
            print('WARNNING: Unusual number of hydrogen in Histidine. HSD is default.')
            return 'HSD'

class dihedral(object):
    def __init__(self, chain, res, resnum, strng, val):
        self.chain = chain
        self.res = res
        self.resnum = resnum
        self.strng = strng
        self.val = val
        self.atoms = []
        
    def add_atoms(self, atm):
        self.atoms.append(atm)
        
class dihedral_manipulator(object):
    def __init__(self):
        self.resi = []
        
    def add_residue(self, res):
        self.resi.append(res)

class Molecular_Rigid_Manipulation(object):
    # http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
    SC_ROT = {'ch1':{'ARG':{'N CA CB CG':['HB1','HB2','CG','HG1','HG2','CD','HD1','HD2','NE','HE','CZ','NH1','HH11','HH12','NH2','HH22','HH21']},\
                     'ASN':{'N CA CB CG':['HB1','HB2','CG','OD1','ND2','HD21','HD22']},\
                     'ASP':{'N CA CB CG':['HB1','HB2','CG','OD1','OD2']},\
                     'CYS':{'N CA CB SG':['HB1','HB2','SG','HG1']},\
                     'GLN':{'N CA CB CG':['HB1','HB2','CG','HG1','HG2','CD','OE1','NE2','HE21','HE22']},\
                     'GLU':{'N CA CB CG':['HB1','HB2','CG','HG1','HG2','CD','OE1','OE2']},\
                     'HSD':{'N CA CB CG':['HB1','HB2','CG','ND1','HD1','CD2','HD2','CE1','HE1','NE2']},\
                     'HSE':{'N CA CB CG':['HB1','HB2','CG','ND1','CD2','HD2','CE1','HE1','NE2','HE2']},\
                     'HSP':{'N CA CB CG':['HB1','HB2','CG','ND1','HD1','CD2','HD2','CE1','HE1','NE2','HE2']},\
                     'ILE':{'N CA CB CG1':['HB','CG1','HG11','HG12','CG2','HG21','HG22','HG23','CD','HD1','HD2','HD3']},\
                     'LEU':{'N CA CB CG':['HB1','HB2','CG','HG','CD1','HD11','HD12','HD13','CD2','HD21','HD22','HD23']},\
                     'LYS':{'N CA CB CG':['HB1','HB2','CG','HG1','HG2','CD','HD1','HD2','CE','HE1','HE2','NZ','HZ1','HZ2','HZ3']},\
                     'MET':{'N CA CB CG':['HB1','HB2','CG','HG1','HG2','SD','CE','HE1','HE2','HE3']},\
                     'PHE':{'N CA CB CG':['HB1','HB2' 'CG','CD1','HD1','CD2','HD2','CE1','HE1','CE2','HE2','CZ','HZ']},\
                     'PRO':{'N CA CB CG':['HB1','HB2','CG','HG1','HG2','CD','HD1','HD2']},\
                     'SER':{'N CA CB OG':['HB1','HB2','OG','HG1']},\
                     'THR':{'N CA CB OG1':['HB','OG1','HG1','CG2','HG21','HG22','HG23']},\
                     'TRP':{'N CA CB CG':['HB1','HB2','CG','CD1','HD1','CD2','CE2','NE1','HE1','CZ2','HZ2','CE3','HE3','CH2','HH2','CZ3','HZ3']},\
                     'TYR':{'N CA CB CG':['HB1','HB2','CG','CD1','HD1','CD2','HD2','CE1','HE1','CE2','HE2','CZ','OH','HH']},\
                     'VAL':{'N CA CB CG1':['HB','CG1','HG11','HG12','HG13','CG2','HG21','HG22','HG23']}},
              'ch2':{'ARG':{'CA CB CG CD':['HG1','HG2','CD','HD1','HD2','NE','HE','CZ','NH1','HH11','HH12','NH2','HH22','HH21']},\
                     'ASN':{'CA CB CG OD1':['OD1','ND2','HD21','HD22']},\
                     'ASP':{'CA CB CG OD1':['OD1','OD2']},\
                     'GLN':{'CA CB CG CD':['HG1','HG2','CD','OE1','NE2','HE21','HE22']},\
                     'GLU':{'CA CB CG CD':['HG1','HG2','CD','OE1','OE2']},\
                     'HSD':{'CA CB CG ND1':['ND1','HD1','CD2','HD2','CE1','HE1','NE2']},\
                     'HSE':{'CA CB CG ND1':['ND1','CD2','HD2','CE1','HE1','NE2','HE2']},\
                     'HSP':{'CA CB CG ND1':['ND1','HD1','CD2','HD2','CE1','HE1','NE2','HE2']},\
                     'ILE':{'CA CB CG1 CD':['HG11','HG12','CD','HD1','HD2','HD3']},\
                     'LEU':{'CA CB CG CD1':['HG','CD1','HD11','HD12','HD13','CD2','HD21','HD22','HD23']},\
                     'LYS':{'CA CB CG CD':['HG1','HG2','CD','HD1','HD2','CE','HE1','HE2','NZ','HZ1','HZ2','HZ3']},\
                     'MET':{'CA CB CG SD':['HG1','HG2','SD','CE','HE1','HE2','HE3']},\
                     'PHE':{'CA CB CG CD1':['CD1','HD1','CD2','HD2','CE1','HE1','CE2','HE2','CZ','HZ']},\
                     'PRO':{'CA CB CG CD':['HG1','HG2','CD','HD1','HD2']},\
                     'TRP':{'CA CB CG CD1':['CG','CD1','HD1','CD2','CE2','NE1','HE1','CZ2','HZ2','CE3','HE3','CH2','HH2','CZ3','HZ3']},\
                     'TYR':{'CA CB CG CD1':['CG','CD1','HD1','CD2','HD2','CE1','HE1','CE2','HE2','CZ','OH','H']}},              
              'ch3':{'ARG':{'CB CG CD NE':['HD1','HD2','NE','HE','CZ','NH1','HH11','HH12','NH2','HH22','HH21']},\
                     'GLN':{'CB CG CD OE1':['OE1','NE2','HE21','HE22']},\
                     'GLU':{'CB CG CD OE1':['OE1','OE2']},\
                     'LYS':{'CB CG CD CE':['HD1','HD2','CE','HE1','HE2','NZ','HZ1','HZ2','HZ3']},\
                     'MET':{'CB CG SD CE':['CE','HE1','HE2','HE3']}},
              'ch4':{'ARG':{'CG CD NE CZ':['HE','CZ','NH1','HH11','HH12','NH2','HH22','HH21']},\
                     'LYS':{'CG CD CE NZ':['HE1','HE2','NZ','HZ1','HZ2','HZ3']}},  
              'ch5':{'ARG':{'CD NE CZ NH1':['NH1','HH11','HH12','NH2','HH22','HH21']}}}
    def __init__(self,parameter_object):
        self.parms = parameter_object
        self.dm = dihedral_manipulator()

    def prep_dihedrals(self, strctr, model): 
        """ Add dihedral atoms from topology information excluding impropers from IC entries """
        for h in strctr.get_models():
            if h.get_id() == model:
                for i in h.get_chains():
                    for j in i.get_residues():
                        res = []
                        for k in self.parms.AA[j.get_resname()].ic:
                            if ' '.join(k[:4]).find('*') == -1:
                                res.append(dihedral(i.get_id(),j.get_resname(),j.get_id()[1],k[:4],k[6]))
                        self.dm.add_residue(res)
        # Considers residues only between 2 and n-1. Eliminares - + and assigns 
        # proper resnum. This is becase residue 1 and n are used for alignment.
        last = {}
        for i in self.dm.resi:
            for j in i:
                last[j.chain] = j.resnum
        for i in self.dm.resi:
            for j in i:
                if j.resnum != 1 and j.resnum != last[j.chain]:
                    j.atoms = []
                    for k in j.strng:
                        if k[0] == '-':
                            j.atoms.append(strctr[model][j.chain][j.resnum-1][k[1:]])
                        elif k[0] == '+':
                            j.atoms.append(strctr[model][j.chain][j.resnum+1][k[1:]])
                        else:
                            j.atoms.append(strctr[model][j.chain][j.resnum][k[:]])
        '''OLD: for i in self.dm.resi:
            for j in i:
                if j.resnum != 1 and j.resnum != last[j.chain]:
                    for k in j.strng:
                        if k[0] == '-':
                            j.atoms.append(strctr[model][j.chain][j.resnum-1][k[1:]])
                        elif k[0] == '+':
                            j.atoms.append(strctr[model][j.chain][j.resnum+1][k[1:]])
                        else:
                            j.atoms.append(strctr[model][j.chain][j.resnum][k[:]])'''

    def translate_molecule(self, strctr, model, direction):
        for i in strctr.get_models():
            if i.get_id() == model:
                for j in i.get_atoms():
                    j.set_coord(j.get_coord() + direction)

    def translate_molecule_model(self, model, direction):
            for j in model.get_atoms():                
                j.set_coord(j.get_coord() + direction)

    def center_molecule(self, center):
        return [x * -1 for x in center]
    '''Deprecated and Moved input_output.py
    def fix_pdb_from_CHARMM(self,path,a,b):
        """ This fixes the problem of CHARMM not generating a chain id.
        But needs to be tested in more systems """
        file_lines = []
        count = 0
        infile = open(path, 'r')
        off1 = a + 5
        off2 = a + 3
        with infile as f:
            for line in f:
                if line[0:4] == "ATOM":
                    #get chain from column count location
                    if count == 0:
                        print('Check the following lines to ensure right changes:')
                        print(line)
                    line = line.replace(line[a:off1],line[a:off2]+' '+line[b],1)
                    file_lines.append(line)
                    if count == 0:
                        print(line)
                    count += 1
                else:
                    file_lines.append(line)
        infile.close()
        outfile = open(path,'w')
        for i in file_lines:
            outfile.write(i)
        outfile.close()'''
    
    def print_residue_atoms(self, struct, s_id, chain, res_num):
        for i in struct[s_id]:
            if i.get_id() == chain:
                for j in i:
                    if j.get_id()[1] == res_num:
                        for k in j:
                            print k

    def delete_atom(self, struct, s_id, chain, res_num, atom_id):
        for i in struct[s_id]:
            if i.get_id() == chain:
                for j in i:
                    if j.get_id()[1] == res_num:
                        j.detach_child(atom_id)
    
    def get_vector_list(self, strctr, mdl, chain, resnum, atom_list):
        vect_list = []
        for i in strctr.get_models():
            if i.get_id() == mdl:
                for j in i.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_residues():
                            if k.get_id()[1] == resnum:
                                for l in k:
                                    if l.id in atom_list:
                                        vect_list.append(l)
        return vect_list

    def vector_lists_rmsd(self,a,b):
        if len(a) == len(b):
            temp = 0
            for i in range(0,len(a)):
                temp = temp + np.sqrt(sum(np.power(a[i].get_coord()-b[i].get_coord(),2)))
            return temp/len(a)
        else:
            print("ERROR: the two input lists must be equal in length. Returning -1")
            return -1

    # resi Amino Acid three leter id. e.g LYS, ALA etc. that are found in top file.
    # group_type == amino, carbonyl, side, all.
    # value 0 = name1, 1 = name2, 2 = partial charge
    def get_residue_component(self, resi, group_type):
        resi =  resi.upper()
        group_type = group_type.upper()
        count_group = 0
        atms = []
        if group_type == 'AMINO':
            for i in self.parms.AA[resi].atoms:
                if count_group == 0:
                    for j in i:
                        atms.append(j)
                count_group += 1
        elif group_type == 'CARBONYL':
            for i in self.parms.AA[resi].atoms:
                if count_group == self.parms.AA[resi].current_groups:
                    for j in i:
                        atms.append(j)
                count_group += 1
        elif group_type == 'SIDE':
            for i in self.parms.AA[resi].atoms:
                if count_group > 0 and count_group < self.parms.AA[resi].current_groups:
                    for j in i:
                        atms.append(j)
                count_group += 1
        elif group_type == 'ALL':
            for i in self.parms.AA[resi].atoms:
                for j in i:
                    atms.append(j)
                count_group += 1
        return atms
    
    # The Following function works for baackbone dihedraal angles only
    def select_atoms_rotate(self, strctr, model, chain, res, dihed, terminus):
        # TODO I changed the algorith so that it works
        num_atoms = 0
        for i in strctr.get_models():
            if i.get_id() == model:
                for j in i.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_atoms():
                            num_atoms += 1
        rot = np.zeros((num_atoms,1), dtype=bool)
        # PSI angel
        if dihed == 'N CA C +N':
            if terminus == 'NTERM':
                current_atom = -1
                for i in strctr.get_models():
                    if i.get_id() == model:
                        for j in i.get_chains():
                            if j.get_id() == chain:
                                for k in j.get_residues():
                                    for l in k.get_atom():
                                        current_atom += 1
                                        if k.get_id()[1] == res:
                                            if l.id in self.get_residue_component(k.resname,'amino')+\
                                                       self.get_residue_component(k.resname,'side'):
                                                if k.resname != 'PRO':
                                                    if l.id not in ['CA']:
                                                        rot[current_atom] = True
                                        elif k.get_id()[1] > res:
                                            pass
                                        elif k.get_id()[1] < res:
                                            rot[current_atom] = True
            elif terminus == 'CTERM':
                # TODO: CTERM needs testing
                print('WARNING: CTERM is not tested sufficiently.')
                sys.exit(1)
                current_atom = -1
                for i in strctr.get_models():
                    if i.get_id() == model:
                        for j in i.get_chains():
                            if j.get_id() == chain:
                                for k in j.get_residues():
                                    for l in k.get_atom():
                                        current_atom += 1
                                        if k.get_id()[1] == res:
                                            if l.id in self.get_residue_component(k.resname,'carbonyl'):
                                                rot[current_atom] = True
                                        elif k.get_id()[1] > res:
                                            rot[current_atom] = True
            else:
                print("ERROR: Only CTERM or NTERM are valid options for terminus.")
            return rot
        # PHI angle
        elif dihed == '-C N CA C':
            if terminus == 'NTERM':
                current_atom = -1
                for i in strctr.get_models():
                    if i.get_id() == model:
                        for j in i.get_chains():
                            if j.get_id() == chain:
                                for k in j.get_residues():
                                    if k.get_resname() == 'PRO':
                                        pass
                                    else:
                                        for l in k.get_atom():
                                            current_atom += 1
                                            if k.get_id()[1] == res:
                                                # N should also be true but it would be an unnecesary rotation.
                                                if l.id == 'HN':
                                                    rot[current_atom] = True
                                            elif k.get_id()[1] > res:
                                                pass
                                            elif k.get_id()[1] < res:
                                                rot[current_atom] = True
            elif terminus == 'CTERM':
                # TODO: CTERM needs testing
                print('WARNING: CTERM is not tested sufficiently.')
                sys.exit(1)
                current_atom = -1
                for i in strctr.get_models():
                    if i.get_id() == model:
                        for j in i.get_chains():
                            if j.get_id() == chain:
                                for k in j.get_residues():
                                    if k.get_resname() == 'PRO':
                                        pass
                                    else:
                                        for l in k.get_atom():
                                            current_atom += 1
                                            if k.get_id()[1] == res:
                                                if l.id in ['HA']+self.get_residue_component(k.resname,'side')+self.get_residue_component(k.resname,'carbonyl'):
                                                    rot[current_atom] = True
                                            elif k.get_id()[1] > res:
                                                rot[current_atom] = True
                                            elif k.get_id()[1] < res:
                                                pass
            else:
                print("ERROR: Only CTERM or NTERM are valid options for terminus.")
            return rot
        else:
            print('ERROR: ONLY phy psy angles can be adjusted for now.')
            return rot

    def select_atoms_rotate_model(self, model, chain, res, aa, dname, dihed, terminus):
        # TODO I changed the algorith so that it works
        num_atoms = 0
        for j in model.get_chains():
            if j.get_id() == chain:
                for k in j.get_atoms():
                    num_atoms += 1

        rot = np.zeros((num_atoms,1), dtype=bool)
        # PSI angel
        if dihed == 'N CA C +N':
            if terminus == 'NTERM':
                current_atom = -1
                for j in model.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_residues():
                            for l in k.get_atom():
                                current_atom += 1
                                if k.get_id()[1] == res:
                                    if l.id in self.get_residue_component(k.resname,'amino')+\
                                               self.get_residue_component(k.resname,'side'):
                                        if k.resname != 'PRO':
                                            if l.id not in ['CA']:
                                                rot[current_atom] = True
                                elif k.get_id()[1] > res:
                                    pass
                                elif k.get_id()[1] < res:
                                    rot[current_atom] = True
            elif terminus == 'CTERM':
                # TODO: CTERM needs testing
                print('WARNING: CTERM is not tested sufficiently.')
                sys.exit(1)
                current_atom = -1
                for j in model.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_residues():
                            for l in k.get_atom():
                                current_atom += 1
                                if k.get_id()[1] == res:
                                    if l.id in self.get_residue_component(k.resname,'carbonyl'):
                                        rot[current_atom] = True
                                elif k.get_id()[1] > res:
                                    rot[current_atom] = True
            else:
                print("ERROR: Only CTERM or NTERM are valid options for terminus.")
            return rot
        # PHI angle
        elif dihed == '-C N CA C':
            if terminus == 'NTERM':
                current_atom = -1
                for j in model.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_residues():
                            if k.get_resname() == 'PRO':
                                pass
                            else:
                                for l in k.get_atom():
                                    current_atom += 1
                                    if k.get_id()[1] == res:
                                        # N should also be true but it would be an unnecesary rotation.
                                        if l.id == 'HN':
                                            rot[current_atom] = True
                                    elif k.get_id()[1] > res:
                                        pass
                                    elif k.get_id()[1] < res:
                                        rot[current_atom] = True
            elif terminus == 'CTERM':
                # TODO: CTERM needs testing
                print('WARNING: CTERM is not tested sufficiently.')
                sys.exit(1)
                current_atom = -1
                for j in model.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_residues():
                            if k.get_resname() == 'PRO':
                                pass
                            else:
                                for l in k.get_atom():
                                    current_atom += 1
                                    if k.get_id()[1] == res:
                                        if l.id in ['HA']+self.get_residue_component(k.resname,'side')+self.get_residue_component(k.resname,'carbonyl'):
                                            rot[current_atom] = True
                                    elif k.get_id()[1] > res:
                                        rot[current_atom] = True
                                    elif k.get_id()[1] < res:
                                        pass
            else:
                print("ERROR: Only CTERM or NTERM are valid options for terminus.")
            return rot
        # CHI 1
        else:
            current_atom = -1
            for j in model.get_chains():
                if j.get_id() == chain:
                    for k in j.get_residues():
                        for l in k.get_atom():
                            current_atom += 1
                            if k.get_id()[1] == res:
                                if k.resname != 'PRO':
                                    if l.id in self.SC_ROT[dname][aa][dihed]:
                                        rot[current_atom] = True
                                else:
                                    # TODO could probably figure out a way move PRO a little bit
                                    pass
            return rot
    def get_init_dih_ang(self,chn,resnm,dhda):
        last = 0
        for i in self.dm.resi:
            for j in i:
                if j.chain == chn:
                    last = j.resnum
        for i in self.dm.resi:
            for j in i:
                if j.chain == chn:
                    if j.resnum != 1 and j.resnum != last:
                        if j.resnum == resnm:
                            if ' '.join(j.strng) == dhda:
                                a1 = j.atoms[0].get_vector()
                                a2 = j.atoms[1].get_vector()
                                a3 = j.atoms[2].get_vector()
                                a4 = j.atoms[3].get_vector()
                                return calc_dihedral(a1, a2, a3, a4)
        return np.nan

    def pivot_and_orient(self, chn, resnm, dname, dhda):
        last = 0
        for i in self.dm.resi:
            for j in i:
                if j.chain == chn:
                    last = j.resnum
        for i in self.dm.resi:
            for j in i:
                if j.chain == chn:
                    if j.resnum != 1 and j.resnum != last:
                        if j.resnum == resnm:
                            if ' '.join(j.strng) == dhda:
                                # TODO for NTER only...
                                if (dname == 'psi') or (dname == 'phi'):
                                    pivot = j.atoms[2].get_vector()
                                    pivot_name = j.atoms[2].get_name()
                                    pivot_segid = j.atoms[2].get_full_id()[3][1]
                                    n1 = j.atoms[1].get_vector() - j.atoms[2].get_vector()
                                    rot_name = j.atoms[1].get_name()
                                    rot_segid = j.atoms[1].get_full_id()[3][1]
                                    return [pivot, n1, pivot_name,pivot_segid,rot_name,rot_segid]
                                elif dname[0:2] == 'ch':
                                    pivot = j.atoms[1].get_vector()
                                    pivot_name = j.atoms[1].get_name()
                                    pivot_segid = j.atoms[1].get_full_id()[3][1]
                                    n1 = j.atoms[2].get_vector() - j.atoms[1].get_vector()
                                    rot_name = j.atoms[2].get_name()
                                    rot_segid = j.atoms[2].get_full_id()[3][1]
                                    return [pivot, n1, pivot_name,pivot_segid,rot_name,rot_segid]
        sys.exit(1)
        return [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
    ''' Deprecated it was imporved by pivot_and_orient
    def pivot_and_orient2(self, chn, resnm, dhda):
        last = 0
        for i in self.dm.resi:
            for j in i:
                if j.chain == chn:
                    last = j.resnum
        for i in self.dm.resi:
            count_dih = 0
            for j in i:
                if j.chain == chn:
                    if j.resnum != 1 and j.resnum != last:
                        if count_dih < 2:
                            if j.resnum == resnm:
                                if ' '.join(j.strng) == dhda:
                                    pivot = j.atoms[2].get_vector()
                                    pivot_name = j.atoms[2].get_name()
                                    pivot_segid = j.atoms[2].get_full_id()[3][1]
                                    n1 = j.atoms[1].get_vector() - j.atoms[2].get_vector()
                                    rot_name = j.atoms[1].get_name()
                                    rot_segid = j.atoms[1].get_full_id()[3][1]
                        count_dih += 1
        return [pivot, n1, pivot_name,pivot_segid,rot_name,rot_segid]
    '''
    def dihedral_rotation(self, strctr, model, chain, rsnum, aa, dname, dih_ang, term, angl):
        m1 = rotaxis(angl, Vector(1, 0, 0))
        pivot_r1 = self.pivot_and_orient(chain,rsnum,dih_ang)
        # Centers the molecule to the 2nd atom in the dihedral
        self.translate_molecule(strctr,model,self.center_molecule(list(pivot_r1[0])))
        r1 = rotmat(pivot_r1[1], Vector([1,0,0]))
        for i in strctr.get_models():
            if i.get_id() == model:
                for j in i.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_atoms():
                            vv = Vector(k.get_coord()).left_multiply(r1)
                            k.set_coord(vv.get_array())
        rot = self.select_atoms_rotate(strctr,model,chain,rsnum,aa,dname,dih_ang,term)
        current_atom = 0
        for i in strctr.get_models():
            if i.get_id() == model:
                for j in i.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_atoms():
                            if rot[current_atom]:
                                vv = Vector(k.get_coord()).left_multiply(m1)
                                k.set_coord(vv.get_array())
                            current_atom += 1
        r1 = rotmat(Vector([1,0,0]),pivot_r1[1])
        for i in strctr.get_models():
            if i.get_id() == model:
                for j in i.get_chains():
                    if j.get_id() == chain:
                        for k in j.get_atoms():
                            vv = Vector(k.get_coord()).left_multiply(r1)
                            k.set_coord(vv.get_array())
        self.translate_molecule(strctr,model,list(pivot_r1[0]))

    def dihedral_rotation_model(self, model, chain, rsnum, aa, dname, dih_ang, term, angl):
        m1 = rotaxis(angl, Vector(1, 0, 0))
        pivot_r1 = self.pivot_and_orient(chain,rsnum,dname,dih_ang)
        # Centers the molecule to the 2nd atom in the dihedral
        self.translate_molecule_model(model,self.center_molecule(list(pivot_r1[0])))
        r1 = rotmat(pivot_r1[1], Vector([1,0,0]))
        for j in model.get_chains():
            if j.get_id() == chain:
                for k in j.get_atoms():
                    vv = Vector(k.get_coord()).left_multiply(r1)
                    k.set_coord(vv.get_array())
        rot = self.select_atoms_rotate_model(model,chain,rsnum,aa,dname,dih_ang,term)
        current_atom = 0
        for j in model.get_chains():
            if j.get_id() == chain:
                for k in j.get_atoms():
                    if rot[current_atom]:
                        vv = Vector(k.get_coord()).left_multiply(m1)
                        k.set_coord(vv.get_array())
                    current_atom += 1

        r1 = rotmat(Vector([1,0,0]),pivot_r1[1])
        for j in model.get_chains():
            if j.get_id() == chain:
                for k in j.get_atoms():
                    vv = Vector(k.get_coord()).left_multiply(r1)
                    k.set_coord(vv.get_array())
        self.translate_molecule_model(model,list(pivot_r1[0]))

    def model_energy(self, model,chnA,aaA,resA):
        # energy is assumed to be -1 or favorable, then try to falsify, if so, energy = +1
        for i in model.get_chains():
            if i.get_id() == chnA:
                for j in i.get_residues():
                    if j.get_id()[1] == resA:
                        for k in j.get_atom():
                            atm1 = self.parms.NONBONDED[self.parms.AA[aaA].atom_type[k.get_id()]][2]
                            for l in model.get_residues():
                                if ((l.get_id()[1]+1) == resA) and (l.get_parent().get_id() == chnA):
                                    for m in l.get_atom():
                                        if m.get_id() != 'C':
                                            if m.get_id() == 'OT1' or m.get_id() == 'OT2':
                                                atm2 = self.parms.NONBONDED[self.parms.AA['CTER'].atom_type[m.get_id()]][2]
                                                sum_radii = float(atm1) + float(atm2)
                                                dist = np.linalg.norm(k.get_coord()-m.get_coord())
                                                if dist < sum_radii:
                                                    #print(aaA,k.get_id(),l.get_resname(),m.get_id(),dist,sum_radii)
                                                    return 1
                                            else:
                                                atm2 = self.parms.NONBONDED[self.parms.AA[l.get_resname()].atom_type[m.get_id()]][2]
                                                sum_radii = float(atm1) + float(atm2)
                                                dist = np.linalg.norm(k.get_coord()-m.get_coord())
                                                if dist < sum_radii:
                                                    #print(aaA,k.get_id(),l.get_resname(),m.get_id(),dist,sum_radii)
                                                    return 1
                                elif ((l.get_id()[1]-1) == resA)  and (l.get_parent().get_id() == chnA):
                                    for m in l.get_atom():
                                        if m.get_id() != 'N':
                                            # TODO do the same for NTER
                                            if m.get_id() == 'OT1' or m.get_id() == 'OT2':
                                                atm2 = self.parms.NONBONDED[self.parms.AA['CTER'].atom_type[m.get_id()]][2]
                                                sum_radii = float(atm1) + float(atm2)
                                                dist = np.linalg.norm(k.get_coord()-m.get_coord())
                                                if dist < sum_radii:
                                                    #print(aaA,k.get_id(),l.get_resname(),m.get_id(),dist,sum_radii)
                                                    return 1
                                            else:
                                                atm2 = self.parms.NONBONDED[self.parms.AA[l.get_resname()].atom_type[m.get_id()]][2]
                                                sum_radii = float(atm1) + float(atm2)
                                                dist = np.linalg.norm(k.get_coord()-m.get_coord())
                                                if dist < sum_radii:
                                                    #print(aaA,k.get_id(),l.get_resname(),m.get_id(),dist,sum_radii)
                                                    return 1
                                elif (l.get_id()[1] == resA) and (l.get_parent().get_id() == chnA):
                                    pass
                                else:
                                    for m in l.get_atom():
                                        if m.get_id() == 'OT1' or m.get_id() == 'OT2':
                                            atm2 = self.parms.NONBONDED[self.parms.AA['CTER'].atom_type[m.get_id()]][2]
                                            sum_radii = float(atm1) + float(atm2)
                                            dist = np.linalg.norm(k.get_coord()-m.get_coord())
                                            if dist < sum_radii:
                                                #print(aaA,k.get_id(),l.get_resname(),m.get_id(),dist,sum_radii)
                                                return 1
                                        else:
                                            atm2 = self.parms.NONBONDED[self.parms.AA[l.get_resname()].atom_type[m.get_id()]][2]
                                            sum_radii = float(atm1) + float(atm2)
                                            dist = np.linalg.norm(k.get_coord()-m.get_coord())
                                            if dist < sum_radii:
                                                #print(aaA,k.get_id(),l.get_resname(),m.get_id(),dist,sum_radii)
                                                return 1
        return -1
        
    def genMatrix2(self, angle, axis):
        rcos = np.cos(angle);
        rsin = np.sin(angle)
        matrix = np.zeros((3,3))
        u = axis[0]
        v = axis[1]
        w = axis[2]
        matrix[0][0] =      rcos + u*u*(1-rcos)
        matrix[1][0] =  w * rsin + v*u*(1-rcos)
        matrix[2][0] = -v * rsin + w*u*(1-rcos)
        matrix[0][1] = -w * rsin + u*v*(1-rcos)
        matrix[1][1] =      rcos + v*v*(1-rcos)
        matrix[2][1] =  u * rsin + w*v*(1-rcos)
        matrix[0][2] =  v * rsin + u*w*(1-rcos)
        matrix[1][2] = -u * rsin + v*w*(1-rcos)
        matrix[2][2] =      rcos + w*w*(1-rcos)
        return matrix

    # REMEBER WEBSITE
    #http://www.j3d.org/matrix_faq/matrfaq_latest.html
    # This method works but it deprecated in favor of Bio.Vector class
    def alignVectors(self, v1,v2):
        v1 = v1/np.linalg.norm(v1)
        v2 = v2/np.linalg.norm(v2)
        axis = np.cross(v1,v2)
        # This IF checks if vextors are colinear.
        if axis[0] == 0 and axis[1] == 0 and axis[2] == 0:
            #If identical they are colinear and the same, return [0,0,0] and
            #be sure to catch this return when generating the rotation matrix.
            if v1[0] == v2[0] and v1[1] == v2[1] and v1[2] == v2[2]:
                return np.identity((3,3))
            else:
                # if colinear and oposite, do it in two steps 
                # go from v2 to v3 = [1,0,0] then from v3 to v1
                v3 = np.array([1,0,0])
                # v3 could be colinear with v1 by chance, then
                # go from v2 tp v3 = [0,1,0] and then go from v3 to v1
                if v3[0] == v1[0] and v3[1] == v1[1] and v3[2] == v1[2]:
                    v3 = np.array([0,1,0])
        
                # Now do transform in two steps V2 to v3 followed by v3 to v1
                v3 = v3/np.linalg.norm(v3)
                axis = np.cross(v3,v2)
                # It is possible that when v3 is colinear or identical with v2
                # so we overide v3 and axis to a new intermediary vector.
                if axis[0] == 0 and axis[1] == 0 and axis[2] == 0:
                    v3 = np.array([0,1,0])
                    v3 = v3/np.linalg.norm(v3)
                    axis = np.cross(v3,v2)
        
                # Then we rotate v2 to v3
                angle = np.arccos(np.dot(v3,v2))
                axis = axis/np.linalg.norm(axis)
                v2 = np.dot(v3,self.genMatrix(angle,axis))
                # With v2 aligned to v3, 
                v2 = v2/np.linalg.norm(v2)
                axis = np.cross(v1,v2)
                angle = np.arccos(np.dot(v1,v2))
                axis = axis/np.linalg.norm(axis)
                return self.genMatrix2(angle,axis)
        else:
            angle = np.arccos(np.dot(v1,v2))
            axis = axis/np.linalg.norm(axis)
            return self.genMatrix2(angle,axis)
