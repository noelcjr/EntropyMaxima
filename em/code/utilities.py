# -*- coding: utf-8 -*-
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
aa = 'ARNDCEQGHILKMFPSTWYV'
polar_charged_aa = 'RHKDESTNQ'
non_polar_hydrophobic = 'CGPAVILMFYW'
residueDict1_1 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', \
                  'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', \
                  'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', \
                  'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
residueDict1_2 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', \
                  'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', \
                  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', \
                  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
# TODO This dihedral atom manipulation could be done from parameters alone and IC tables, but I am not sure if or how.
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
