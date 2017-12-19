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

class structure(object):
    def __init__(self,path, pdb_id='pdb', cif_id='cif'):
        # copied from input_output.py _read_structure
        self.file_name = os.path.basename(path).split('.')[0]
        self.file_sufix = os.path.basename(path).split('.')[1]
        self.dir_path = os.path.dirname(path)
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
            self.seq = sequence(sequence_string_or_file,seq_id,'str')
        elif seq_format == 'fasta':
            self.seq = sequence(sequence_string_or_file,'','fasta')
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
                print(a)
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
    
    def summary(self, model_id=0):
        print('Results from ', self.dir_path+'/'+self.file_name+'.'+self.file_sufix)
        print("Model="+str(model_id)+". For a different model, pass model number as a parameter(default=0) to summary.")
        print('--------------------------------')
        for i in self.structure.get_models():
            if i.get_id() == model_id:
                self.min_max(i)
                print('--------------------------------')
                self.align_structure_and_sequenceinfo(i)
                print('--------------------------------')

#Test
strt_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1BRS/fit_complete/'
sqnc_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Alignments/'
#1brs_D_fasta.txt'
prot_str = structure(strt_path+'1brs.pdb')
prot_str.add_sequence(sqnc_path+'1brs.fasta.txt')
prot_str.summary()

    