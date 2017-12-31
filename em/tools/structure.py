#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 10:15:22 2017

@author: noel
"""
import os
import sys
import Bio.PDB as struct
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import em.describe.utilities as utilities
import em.tools.sequence as sequence
import em.tools.CHARMM_Parser as CP
import pandas as pd
import numpy as np

class structure(object):
    def __init__(self,path, pdb_id='pdb', cif_id='cif'):
        # copied from input_output.py _read_structure
        self.file_name = os.path.basename(path).split('.')[0]
        self.file_sufix = os.path.basename(path).split('.')[1]
        self.dir_path = os.path.dirname(path)
        self.params = CP.read_charmm_FF()
        self.chains = []
        self.models = {}
        if self.file_sufix == 'pdb':
            self.header = struct.parse_pdb_header(path)
            self.structure = struct.PDBParser(QUIET=True).get_structure(pdb_id, path)
            self.has_sequence = False
        elif self.file_sufix == 'cif':
            self.header = struct.MMCIF2Dict()
            self.structure = struct.MMCIFParser().get_structure(cif_id, path)
            self.has_sequence = True
        else:
            print("ERROR: Unreognized file type " + self.file_sufix + " in " + self.file_name)
            sys.exit(1)
    
    def add_sequence(self, sequence_string_or_file, seq_id='from_string', seq_format='fasta'):
        if seq_format == 'str':
            self.seq = sequence.sequence(sequence_string_or_file,seq_id,'str')
        elif seq_format == 'fasta':
            self.seq = sequence.sequence(sequence_string_or_file,'','fasta')
        else:
            print("ERROR: "+seq_id+" not suported.")
            sys.exit(1)
    
    def get_sequence_from_XYZ(self, model, chain_id):
        xyz_seq = ''
        #for i in self.structure:
        #if model.get_id() == model_id:
        for j in model.get_chains():
            if j.get_id() == chain_id:
                for k in j.get_residues():
                    if k.get_resname().upper() in utilities.residueDict1_2:
                        xyz_seq += utilities.residueDict1_2[k.get_resname().upper()]
        return xyz_seq
        
    def align_structure_and_sequenceinfo(self, model):
        for j in model.get_chains():
            key_id = prot_str.file_name.upper()+":"+j.get_id()
            top_seq = str(prot_str.seq.chain_sequences[key_id])
            bottom_seq = prot_str.get_sequence_from_XYZ(model,j.get_id())
            alignments = pairwise2.align.globalxx(top_seq, bottom_seq)
            # Use format_alignment method to format the alignments in the list
            print('  XYZ-fasta Alignment for Chain', j.get_id())
            for a in alignments:
                print(format_alignment(*a))
                count = 1
                gaps = []
                for i in alignments[0][1]:
                    if i == '-':
                        gaps.append(count)
                    count += 1
            missing_residues_struct = []
            for k in j.get_residues():
                if k.get_id()[1] != 'W' and k.resname != 'HOH':
                    missing_residues_struct.append(k.get_full_id()[3][1])
            gap_report = []
            for i in self.check_gaps(missing_residues_struct):
                gap_report.append(i[1])
            print('Gap Locations for chain '+j.get_id()+' from structural residue numbers:',gap_report)
            print('**********************************************')
        
    def check_gaps(self,data):
        ''' The following code identify missing indexs in an monotonically
        increasing array. Full sequence is obtained from begining of sequence
        to the last amino acid in the sequence. Then we look for missing amino
        acid indexes in the available structural information.
        '''
        #data=[2,3,4,5,8,9]
        count = 0
        inserts = []
        for i in range(1, data[-1] + 1):
            if i not in data:
                inserts.append((count, i))
            count += 1
        data2 = [i for i in data]
        for i in inserts:
            data2.insert(i[0], i[1])
        return inserts
    
    def min_max(self, model):
        #atom_coordinates = self.get_atom_coordinates()
        #for i in self.structure.get_models():
        for j in model.get_chains():
            atoms = []
            for k in j.get_atoms():
                atoms.append(k.get_coord())
            x,y,z = zip(*atoms)
        
            Xmin,Xmax = min(x), max(x)
            Ymin,Ymax = min(y), max(y)
            Zmin,Zmax = min(z), max(z)
        
            print('  Number of Atoms ', len(atoms), 'Chain', j.get_id())
            print('  Xmin - Xmax', Xmin, ' - ', Xmax)
            print('  Ymin - Ymax', Ymin, ' - ', Ymax)
            print('  Zmin - Zmax', Zmin, ' - ', Zmax)
            print('  Celbasis(Side length) X, Y, Z:', (Xmax - Xmin), (Ymax - Ymin), (Zmax - Zmin))
    
    def check_for_aa_atomtypes_in_parameters(self, parameters):
        chain_check_count = {}
        residue_check_count = {}
        atom_check_count = {}
        models_summary = {}
        for i in self.structure.get_models():
            for j in i.get_chains():
                chain_check_count[j.get_id()] = {}
                for k in j.get_residues():
                    if k not in residue_check_count:
                        residue_check_count[k] = 1
                    else:
                        residue_check_count[k] += 1
                for k in j.get_atoms():
                    if k not in atom_check_count[k]:
                        atom_check_count[k] = 1
                    else:
                        atom_check_count[k] += 1
                chain_check_count[j.get_id()]['res'] = residue_check_count
                chain_check_count[j.get_id()]['atm'] = atom_check_count
            models_summary[i.get_id()] = chain_check_count
    
    def load_struct_into_dictionary(self):
        for i in self.structure.get_models():
            self.models[i.get_id()] = {}
            for j in i.get_chains():
                self.chains.append(j.get_id())
                self.models[i.get_id()][j.get_id()] = {}
                residues = []
                atoms = []
                for k in j.get_residues():
                    residues.append(k.get_resname())
                for k in j.get_atoms():
                    atoms.append(k.get_id())
                self.models[i.get_id()][j.get_id()]['res'] = residues
                self.models[i.get_id()][j.get_id()]['atm'] = atoms
    
    def get_count_report(self):
        print("Amino acid and atom type counts per chain.")
        for selection in ['res','atm']:
            unique = []
            for i in self.models:
                for j in sorted(set(self.models[i])):
                    for k in sorted(set(self.models[i][j][selection])):
                        unique.append(k)
            unique = sorted(set(unique))
            df = pd.DataFrame(index=unique)
            for i in self.models:
                aa_model_sum = {}
                for k in unique:
                    aa_model_sum[k] = 0
                for j in sorted(set(self.models[i])):
                    aa_count = {}
                    for k in unique:
                        aa_count[k] = 0
                    for k in self.models[i][j][selection]:
                        aa_count[k] += 1
                    df[(i,j)] = pd.Series(aa_count)
                    for k in unique:
                        aa_model_sum[k] += aa_count[k]
                df[(i,'Tot')] = pd.Series(aa_model_sum)
            print(df)

    def check_for_missing_paramters(self):
        for i in self.structure.get_models():
            print("Model:"+str(i.get_id()))
            print("  Missing atom parameters:")
            for j in i.get_chains():
                missing_aa = {}
                missing_atm = []
                print("  Chain:"+j.get_id())
                for k in j.get_residues():
                    if k.get_resname() in self.params.AA:
                        for l in k.get_atoms():
                            if l.get_id() not in self.params.AA[k.get_resname()].atom_type:
                                missing_atm.append(str((k.get_id()[1],k.get_resname(),l.get_id())))
                    else:
                        if k.get_resname() not in missing_aa:
                            missing_aa[k.get_resname()] = []
                            missing_aa[k.get_resname()].append(k.get_id()[1])
                        else:
                            missing_aa[k.get_resname()].append(k.get_id()[1])
                print("    Missing AA paramaters:\n"+str(missing_aa))           
                print("    Missing Atom paramaters:\n"+str(missing_atm))
                print('--------------------------------')
    
    def rename(self, feature, value_from, value_to, where_ids = {'chain_id':'\0','res_id':'\0'}):
        """ For ease of coding, it will do it in all models. 
            I see no reason for doing it in only one of the models, and not the others.
            If it is needed to treat a model different, separate the structure by models and
            work on renaming components in each model separately.
        """
        possible_features = ['chain','residue','atom']
        if feature not in possible_features:
            print("Error: feature "+feature+" non existent or not allowed to be modified. Exiting now.")
            sys.exit(1)
        if where_ids['res_id'] != '\0':
            resid_range = [int(i) for i in where_ids['res_id'].split(':')]
            if len(resid_range) == 1:
                resid_range.append(resid_range[0])
            elif len(resid_range) != 2:
                print("Error: only one where_ids range for res_id. Exiting now.")
                sys.exit(1)
        if feature == possible_features[0]: # chain
            if (where_ids['chain_id'] != '\0') and (where_ids['res_id'] != '\0'):
                for i in self.structure.get_chains():
                    if i.get_id() == where_ids['chain_id']:
                        for j in i.get_residues():
                            if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                                temp = [k for k in j.get_full_id()]
                                temp[2] = where_ids['chain_id']
                                j.full_id = tuple(temp)
            elif (where_ids['chain_id'] == '\0') and (where_ids['res_id'] != '\0'):
                for i in self.structure.get_chains():
                    if i.get_id() == value_from:
                        for j in i.get_residues():
                            if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                                temp = [k for k in j.get_full_id()]
                                temp[2] = value_to
                                j.full_id = tuple(temp)
            elif (where_ids['chain_id'] != '\0') and (where_ids['res_id'] == '\0'):
                for i in self.structure.get_chains():
                    if i.get_id() == where_ids['chain_id']:
                        for j in i.get_residues():
                            #if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                            temp = [k for k in j.get_full_id()]
                            temp[2] = where_ids['chain_id']
                            j.full_id = tuple(temp)
            elif (where_ids['chain_id'] == '\0') and (where_ids['res_id'] == '\0'):
                for i in self.structure.get_chains():
                    if i.get_id() == value_from:
                        for j in i.get_residues():
                            #if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                            temp = [k for k in j.get_full_id()]
                            temp[2] = value_to
                            j.full_id = tuple(temp)
        elif feature == possible_features[1]: # residue
            print('1.')
            if (where_ids['chain_id'] != '\0') and (where_ids['res_id'] != '\0'):
                for i in self.structure.get_chains():
                    if i.get_id() == where_ids['chain_id']:
                        for j in i.get_residues():
                            if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                                j.resname = value_to
            elif (where_ids['chain_id'] == '\0') and (where_ids['res_id'] != '\0'):
                for i in self.structure.get_chains():
                    for j in i.get_residues():
                        if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                            j.resname = value_to
            elif (where_ids['chain_id'] != '\0') and (where_ids['res_id'] == '\0'):
                for i in self.structure.get_chains():
                    if i.get_id() == where_ids['chain_id']:
                        for j in i.get_residues():
                            if j.resname == value_from:
                                j.resname = value_to
            elif (where_ids['chain_id'] == '\0') and (where_ids['res_id'] == '\0'):
                print('2.')
                for i in self.structure.get_chains():
                    print('3.', i.get_id())
                    for j in i.get_residues():
                        if j.resname == value_from:
                            j.resname = value_to
        elif feature == possible_features[2]: # atoms
            if (where_ids['chain_id'] != '\0') and (where_ids['res_id'] != '\0'):
                for i in self.structure.get_chains():
                    if i.get_id() == where_ids['chain_id']:
                        for j in i.get_residues():
                            if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                                for k in j.get_atoms():
                                    if k.id == value_from:
                                        k.id = value_to
            elif (where_ids['chain_id'] == '\0') and (where_ids['res_id'] != '\0'):
                for i in self.structure.get_chains():
                    #if i.get_id() == where_ids['chain_id']:
                    for j in i.get_residues():
                        if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                            for k in j.get_atoms():
                                if k.id == value_from:
                                    k.id = value_to
            elif (where_ids['chain_id'] != '\0') and (where_ids['res_id'] == '\0'):
                for i in self.structure.get_chains():
                    if i.get_id() == where_ids['chain_id']:
                        for j in i.get_residues():
                            #if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                            for k in j.get_atoms():
                                if k.id == value_from:
                                    k.id = value_to
            elif (where_ids['chain_id'] == '\0') and (where_ids['res_id'] == '\0'):
                for i in self.structure.get_chains():
                    #if i.get_id() == where_ids['chain_id']:
                    for j in i.get_residues():
                        #if j.get_id()[1] >= resid_range[0] and j.get_id()[1] <= resid_range[1]:
                        for k in j.get_atoms():
                            if k.id == value_from:
                                k.id = value_to

    def summary(self, model_id=0):
        print('Results from ', self.dir_path+'/'+self.file_name+'.'+self.file_sufix)
        list_model_ids = []
        for i in self.structure.get_models():
            list_model_ids.append(i)
        print("This structure has "+str(len(list_model_ids))+" models. Only model 0 will be sumaryzed. For summaries")
        print("of other models. enter the model number in the command line.")
        print("Models:"+str(list_model_ids))
        if len(list_model_ids) != 1:
            print("This structure has multiple models. Only one model summary will be output at the time.")
            print("Model 0 will be the default output. For different models, pass a different id after --model.")
            print("This is the list of model ids ="+str(list_model_ids)+". Enter each one separately for a summary.")
        else:
            pass
        print("Model="+str(model_id))
        print('--------------------------------')
        for i in self.structure.get_models():
            if i.get_id() == model_id:
                self.min_max(i)
                print('--------------------------------')
                self.align_structure_and_sequenceinfo(i)
                print('--------------------------------')
                self.load_struct_into_dictionary()
                self.get_count_report()
                self.check_for_missing_paramters()
                print('--------------------------------')
# Test
#strt_path = '/home/noel/Projects/Protein_design/Insulin/struct_prep/2hiu/2hiu.pdb'
#sqnc_path = '/home/noel/Projects/Protein_design/Insulin/rcsb_seq_struct/2HIU/2hiu.fasta.txt'
#prot_str = structure(strt_path)
#prot_str.add_sequence(sqnc_path)
strt_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1BRS/fit_complete/'
sqnc_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Alignments/'
prot_str = structure(strt_path+'1brs.pdb')
prot_str.add_sequence(sqnc_path+'1brs.fasta.txt')
prot_str.summary()
prot_str.rename('residue','HOH','TIP3')
prot_str.summary()
prot_str.rename('atom','O','OH2',{'chain_id':'A','res_id':'111:255'})
prot_str.rename('atom','O','OH2',{'chain_id':'B','res_id':'111:247'})
prot_str.rename('atom','O','OH2',{'chain_id':'C','res_id':'111:160'})
prot_str.rename('atom','O','OH2',{'chain_id':'D','res_id':'90:166'})
prot_str.rename('atom','O','OH2',{'chain_id':'E','res_id':'90:153'})
prot_str.rename('atom','O','OH2',{'chain_id':'F','res_id':'90:129'})
prot_str.summary()
for h in ['A','B','C']:
    for i in [4,25,51,55,76,88,96,109]:
        prot_str.rename('atom','CD1','CD',{'chain_id':h,'res_id':str(i)})
    prot_str.rename('atom','OXT','O',{'chain_id':h,'res_id':str(110)})
for h in ['D','E','F']:
    for i in [5,10,13,84,86,87]:
        prot_str.rename('atom','CD1','CD',{'chain_id':h,'res_id':str(i)})
    if h in ['E','F']:
        prot_str.rename('atom','OXT','O',{'chain_id':h,'res_id':str(89)})
    
 