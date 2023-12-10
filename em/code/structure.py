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
import em.code.utilities
import em.code.sequence as sqnce
import em.code.structural_expression as SE
import pandas as pd
import numpy as np
import math

class structure(object):
    def __init__(self,str_path):
        self.sequence_path = '\0'
        self.chains = []
        self.models = {}
        self.file_name = os.path.basename(str_path).split('.')[0]
        self.file_sufix = os.path.basename(str_path).split('.')[1]
        self.dir_path = os.path.dirname(str_path)
        if self.file_sufix == 'pdb':
            self.header = struct.parse_pdb_header(str_path)
            self.structure = struct.PDBParser(QUIET=True).get_structure(
                    self.file_name, str_path)
            self.seq = sqnce.sequence()
            self.has_sequence = False
        elif self.file_sufix == 'cif':
            self.structure = struct.MMCIFParser(QUIET=True).get_structure(
                    self.file_sufix, str_path)
            self.seq = sqnce.sequence()
            self.has_sequence = False
        else:
            print("ERROR: Unreognized file type "+self.file_sufix[-1]+" in "+
                  self.file_name[-1])
            sys.exit(1)
               
    def add_sequence(self, seq_format='fasta'):
        """
        Associates a sequence in Fasta format to a structure. This helps find
        missing amino acids at the begining and end of the structure.
        Called by:
            subscripts/gaps.py
            subscripts/summary.py
        """
        if seq_format == 'fasta':
            if self.sequence_path == '\0':
                print("ERROR: "+seq_format+" not found.")
                sys.exit(1)
            self.seq.clear_sequences()
            self.seq.add_Fasta_sequence(self.sequence_path)
            self.has_sequence = True
        else:
            print("ERROR: "+seq_format+" not suported.")
            sys.exit(1)
    
    def gap_report(self, alignment_format=60):
        """
        Checks for gaps in the crystal structure. For greater accuracy, include
        a fasta sequence. Without a fasta sequence, amino acids at the end of
        the structure will be ignored.
        Called by:
            subscripts/gaps.py
        """
        if self.has_sequence:
            self.align_seq_from_structure_and_sequence(alignment_format)
        else:
            print("Gap report for "+self.dir_path+"/"+self.file_name+"."+
                  self.file_sufix+" in tuple format (gap_number, first aa, last aa):")
            first_model = True
            for i in self.structure.get_models():
                if first_model:
                    for j in i.get_chains():
                        aa_present = []
                        for k in j.get_residues():
                            if k.get_id()[1] != 'W' and k.resname != 'HOH':
                                aa_present.append(k.get_full_id()[3][1])
                        print("Gaps in model "+str(i.get_id())+" chain "+j.get_id()+':')
                        self.report_gaps(aa_present)
                    first_model = False

    def report_gaps(self,aa_present):
        """
           Finds gaps. Unsued.
        """
        # TODO This report_gaps can be turned into a more general function that
        # detects breaks in a sequence of integers. found in flower.py?
        inserts = self.check_gaps(aa_present)
        count = 0
        init = True
        first = 0
        last = 0
        report = []
        for l in inserts:
            if init:
                first = l[1]
                last = l[1]
                init = False
                count += 1
            if l[1] == last+1:
                last = l[1]
            else:
                if first != last:
                    report.append((count,first,last))
                    count += 1
                first = l[1]
                last = l[1]
        report.append((count,first,last))
        for c in report:
            print(c)
                      
    def get_sequence_from_XYZ(self, model, chain_id):
        """
           Self explanatory. Unused.
        """
        xyz_seq = ''
        for j in model.get_chains():
            if j.get_id() == chain_id:
                for k in j.get_residues():
                    if k.get_resname().upper() in utilities.residueDict1_2:
                        xyz_seq += utilities.residueDict1_2[k.get_resname().
                                                              upper()]
        return xyz_seq
    
    def align_seq_from_structure_and_sequence(self, alignment_width):
        """
           Self explanatory. Called internally by gap_report. Consider making 
           it private or merge with gap_report function.
        """
        first_model = True
        for i in self.structure.get_models():
            if first_model:
                for j in i.get_chains():
                    key_id = self.file_name.upper()+":"+j.get_id()
                    top_seq = str(self.seq.chain_sequences[key_id])
                    bottom_seq = self.get_sequence_from_XYZ(i,j.get_id())
                    # for alignment parameters:
                    # http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
                    # TODO the points for matching, penalize, gaps and extended 
                    # gaps can be optional input parameters
                    alignments = pairwise2.align.globalms(top_seq, bottom_seq, 2, 
                                                          -1, -.5, -.1)
                    # Use format_alignment method to format alignments in the list
                    print('X-structure-fasta Alignment for Model '+
                          str(i.get_id())+' Chain '+j.get_id())
                    #alignment_count = 0
                    for k in alignments:
                        #print(format_alignment(*k))
                        #print(type(k),type(format_alignment(*k)))
                        #print(k)
                        seq1 = k[0]
                        seq2 = k[1]
                        score = k[2]
                        seq_count = 0
                        num_columns = math.ceil(float(len(seq1)) / 
                                                float(alignment_width))
                        #last_column_width = len(seq1) % alignment_width
                        for kk in range(int(num_columns)):
                            if kk == (num_columns-1):
                                # do the last one
                                print('     '+seq1[seq_count:])
                                print('     '+'|'*(len(seq1)-seq_count)+
                                      ' '*(alignment_width-(len(seq1)-seq_count))+
                                      '   '+str(seq_count+(len(seq1)-seq_count)))
                                print('     '+seq2[seq_count:])
                                print('Alignment Score: '+str(score))
                            else:
                                # print alignment_width at the time
                                print('     '+
                                      seq1[seq_count:(seq_count+alignment_width)])
                                print('     '+'|'*(alignment_width)+'   '+
                                      str(seq_count+alignment_width))
                                print('     '+
                                      seq2[seq_count:(seq_count+alignment_width)]+
                                      '\n')
                                seq_count += alignment_width
                        first_seq_residue = True
                        first_seq_residue_type = ''
                        first_seq_residue_number = 0
                        for l in alignments[0][1]:
                            if first_seq_residue:
                                if l == '-':
                                    first_seq_residue_number += 1
                                else:
                                    first_seq_residue = False
                                    first_seq_residue_type = l
                                    first_seq_residue_number += 1
                        aa_present = []
                        aa_not_present = []
                        count = 0
                        for l in alignments[0][1]:
                            if l == '-':
                                count += 1
                                aa_not_present.append(count)
                            else:
                                count += 1
                                aa_present.append(count)
                        first_X_residue = True
                        first_X_residue_type = ''
                        first_X_residue_number = 0
                        for m in j.get_residues():
                            if first_X_residue:
                               if m.get_resname().upper() in utilities.residueDict1_2:
                                   first_X_residue_type = utilities.residueDict1_2[m.get_resname().upper()]
                               else:
                                   print("    WARNING:First residue in X-structure\
                                         is not a standard amino acid: "+
                                         m.get_resname().upper())
                                   first_X_residue_type = m.get_resname().upper()
                               first_X_residue_number = m.get_full_id()[3][1]
                               first_X_residue = False
                            #if m.get_id()[1] != 'W' and m.resname != 'HOH':
                            # missing_residues_struct.append(m.get_full_id()[3][1])
    
                        if (first_X_residue_number-first_seq_residue_number) == 1:
                            if first_seq_residue_type == first_X_residue_type:
                                print("\nWARNING: The first residue index for "+
                                      "the fasta sequence and the X-structure are"+
                                      " off by one.")
                                print("\nfirst_seq_residue_number and "+
                                      "first_seq_residue_type:  "+
                                      str(first_seq_residue_number)+
                                      " and "+first_seq_residue_type)
                                print("\nfirst_X_residue_number and "+
                                      "first_X_residue_type  :  "+
                                      str(first_X_residue_number)+" and "+
                                      first_X_residue_type)
                        read_sequence = True
                        first_sequence = 0
                        end_sequence = 0
                        count = 1
                        # TODO like report gaps the folowing loop detects breaks 
                        # in integer sequences, and could be turned into a general
                        # function. See report_gaps function.
                        for iii in aa_not_present:
                            if read_sequence:
                                print("    Gaps found based on structure indexing:")
                                first_sequence = iii
                                end_sequence = iii
                                read_sequence = False
                            else:
                                if (end_sequence+1) == iii:
                                    end_sequence = iii
                                else:
                                    print("    ("+str(count)+", "+
                                          str(first_sequence)+
                                          ", "+str(end_sequence)+")")
                                    first_sequence = iii
                                    end_sequence = iii
                                    count += 1
                        if (first_sequence == 0) and (end_sequence == 0):
                            print("     No gaps found")
                        else:
                            print("    ("+str(count)+", "+str(first_sequence)+
                                  ", "+str(end_sequence)+")")
                            print()
                    print("<---------------------------------------------->")
                first_model = False    
        
    def check_gaps(self,data):
        ''' The following code identify missing indexs in an monotonically 
        increasing array. Full sequence is obtained from begining of sequence 
        to the last amino acid in the sequence. Then we look for missing amino
        acid indexes in the available structural information. 
        Ex: data=[2,3,4,5,8,9] would return inserts = [(0, 1), (5, 6), (6, 7)] 
        which provide inserts to complete the array output in 
        data2 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        Unused. 
        '''
        # TODO Move to utilities if this function has more genral applications
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

    def min_max_se(self, optionsse):
        """ This function give the minimum and maximum X, Y, Z coordinates
        of a structure file PDB/CIF. The defailt SE is m[*]c[*]r[*]a[*].
        The coordinates can be used to calculate a box's size around the
        molecule that contains all atoms.
        Called by:
            subscripts/minmax.py
        """
        SE1 = SE.structural_expression(optionsse)
        selected_atoms = SE1.atom_list(self.structure)
        atoms = []
        for i in selected_atoms:
            for j in i:
                atoms.append(j[3].get_coord())
        x,y,z = zip(*atoms)
    
        Xmin,Xmax = min(x), max(x)
        Ymin,Ymax = min(y), max(y)
        Zmin,Zmax = min(z), max(z)
    
        print(self.structure.get_id(),'SE',optionsse,
              'Number of Atoms ', len(atoms))
        print('  Xmin - Xmax', Xmin, ' - ', Xmax)
        print('  Ymin - Ymax', Ymin, ' - ', Ymax)
        print('  Zmin - Zmax', Zmin, ' - ', Zmax)
        print('  Celbasis(Side length) X, Y, Z:', (Xmax - Xmin), 
              (Ymax - Ymin), (Zmax - Zmin))
        
    """
    TODO this function is not used. Consider movig to energy
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
    """
    def load_struct_into_dictionary(self):
        """
        Self Explanatory. Called internally by self.summary(self)
        """
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
        """
        Self Explanatory. Called internally by self.summary(self)
        """
        print("Amino acid and atom type counts per chain.")
        for selection in ['res','atm']:
            unique = []
            for i in self.models:
                for j in sorted(set(self.models[i])):
                    for k in sorted(set(self.models[i][j][selection])):
                        unique.append(k)
            unique = sorted(set(unique))
            df = pd.DataFrame(index=unique)
            first_model = True
            for i in self.models:
                if first_model:
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
                    first_model = False
            with pd.option_context('display.max_rows', None, 
                                   'display.max_columns', None):
                print(df)

    """
    TODO Move this to energy. Parameters are not needed for just structure
         manipulation.
    def check_for_missing_paramters(self):
        self.params = CP.read_charmm_FF()
        first_model = True
        for i in self.structure.get_models():
            if first_model:
                print("Model:"+str(i.get_id()))
                print("  Missing atom parameters:")
                for j in i.get_chains():
                    missing_aa = {}
                    missing_atm = []
                    print("  Chain:"+j.get_id())
                    for k in j.get_residues():
                        if k.get_resname() in self.params.AA:
                            for l in k.get_atoms():
                                pars = self.params.AA[k.get_resname()].atom_type
                                if l.get_id() not in pars:
                                    missing_atm.append(str((k.get_id()[1], 
                                                            k.get_resname(),
                                                            l.get_id())))
                        else:
                            if k.get_resname() not in missing_aa:
                                missing_aa[k.get_resname()] = []
                                missing_aa[k.get_resname()].append(k.get_id()[1])
                            else:
                                missing_aa[k.get_resname()].append(k.get_id()[1])
                    print("    Missing AA paramaters:\n"+str(missing_aa))           
                    print("    Missing Atom paramaters:\n"+str(missing_atm))
                    print('--------------------------------')
                first_model = False
    """
    # TODO rename does not uses structural expressions because it was written 
    # before they were created. A rename is almost finished.
    def rename(self, feature, value_from, value_to, 
               where_ids = {'chain_id':'\0','res_id':'\0'}):
        """ For ease of coding, it will do it in all models. I see no reason 
            for doing it in only one of the models, and not the others.
            If it is needed to treat a model different, separate the structure 
            by models and work on renaming components in each model separately.
        """
        possible_features = ['chain','residue','atom']
        if feature not in possible_features:
            print("Error: feature "+feature+" non existent or not allowed to \
                  be modified. Exiting now.")
            sys.exit(1)
        if where_ids['res_id'] != '\0':
            resid_range = [int(i) for i in where_ids['res_id'].split(':')]
            if len(resid_range) == 1:
                resid_range.append(resid_range[0])
            elif len(resid_range) != 2:
                print("Error: only one where_ids range for res_id. Exiting \
                      now.")
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
        print('Results from ', self.dir_path+'/'+self.file_name+'.'+
              self.file_sufix)
        print("If a structure has multiple models, only the first model will \
be sumaryzed since they are all identical. Only Cell basis will be shown for \
all models because it depends on coordinates.")
        for i in self.structure.get_models():
            if i.get_id() == model_id:
                # TODO mim_max_se takes one optionsse parameter check and find it.
                # self.min_max_se()
                print('--------------------------------')
                self.gap_report()
                print('--------------------------------')
                self.load_struct_into_dictionary()
                self.get_count_report()
                # TODO missing parameter excluded from stctr for Local Maxima.
                #self.check_for_missing_paramters()
                #print('--------------------------------')

    def save_substructure(self, model, path):
        io = struct.PDBIO()
        io.set_structure(model)
        io.save(path)
        
    def save_structure(self, path='\0'):
        io = struct.PDBIO()
        io.set_structure(self.structure)
        if path == '\0':
            path = self.dir_path+'/'+self.file_name+'.pdb'
        io.save(path)
    
    def get_se_output(self, optionsse, optionsout):
        """print pdb directly to avoid structural deepcopies of biopython output
           deepcopy takes a long time. 
        """
        SE1 = SE.structural_expression(optionsse)
        selected_atoms = SE1.atom_list(self.structure)
        group_all = {}
        for i in selected_atoms:
            for j in i:
                if j[0] not in group_all:
                    group_all[j[0]] = {}
                if j[1] not in group_all[j[0]]:
                    group_all[j[0]][j[1]] = {}
                if j[2] not in group_all[j[0]][j[1]]:
                    group_all[j[0]][j[1]][j[2]] = []
                group_all[j[0]][j[1]][j[2]].append(j[3])
        # NOTE This is the official PDB formating line. it requires a special
        #      treatment for the third, atom type, column.
        line = '{:6}{:>5} {:<5}{:<4}{:1}{:>4d}{:>12.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>12}'
        mdl = '{:5}{:9}'
        lines = []
        for i in group_all:
            atm = 1
            lines.append(mdl.format('MODEL',i))
            for j in group_all[i]:
                rsd = 1
                for k in group_all[i][j]:
                    for l in group_all[i][j][k]:
                        atm_typ = l.get_full_id()[4][0]
                        if len(atm_typ) != 4:
                            atm_typ = ' '+atm_typ
                        lines.append(line.format('ATOM',str(atm),\
                                                 atm_typ,\
                                                 l.get_parent().get_resname(),j,\
                                                 int(k),\
                                                 float(l.get_coord()[0]),\
                                                 float(l.get_coord()[1]),\
                                                 float(l.get_coord()[2]),float('0'),\
                                                 float('0'),\
                                                 l.get_full_id()[4][0][0]))
                        atm += 1
                    rsd += 1
            lines.append('ENDMDL')
        if optionsout == "OUT_2_TERMINAL":
            for i in lines:
                print(i)
        else:
            outFile = open(optionsout, 'w')
            for i in lines:
                outFile.write(i + '\n')
            outFile.close()

class crd(object):
    '''
     This fixes when there are multiple chains and the residue number gets
     reset at the beginning of each chain. With this fix, residue numbers
     will be renumbered
      1 - 5 Integer Atom no. sequential
     6 - 10 Integer ires Residue position from file beginning
    11 - 11 Space
    12 - 15 Achar resname Residue name
    16 - 16 Space
    17 - 20 Achar type Atom type, IUPAC name of atom left justified
    21 - 30 Real(10.5) x Orthogonal coordinates for X, Angstroms
    31 - 40 Real(10.5) y Orthogonal coordinates for Y, Angstroms
    41 - 50 Real(10.5) z Orthogonal coordinates for Z, Angstroms
    51 - 51 Space
    52 - 55 Achar segid Segment identifier
    56 - 56 Space
    57 - 60 Achar resid Residue number within the segment
    61 - 70 Real(10.5) Weighting array value
    line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
    '''
    # https://www.charmm.org/charmm/documentation/by-version/c40b1/params/doc/io/#IOFORM
    def __init__(self, crd_file):
        print("A-",crd_file)
        self.CRD_file = crd_file
        print("B-")
        self.crd_lines = self.read_crd_file()
        print("C-")
        self.form_line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
        print("D-")
        self.xyz = np.array([])
        # The following lines gave errors when initialized in __init__, so they were moved to subroutines.
        # self.col_length = [(0,5),(5,11),(11,15),(15,16),(16,20),(20,30),(30,40),(40,50),(50,51),(51,55),(55,56),(56,60),(60,72)]
        # self.col_types = [int,int,str,str,str,float,float,float,str,str,str,int,float]
        # self.col_names = ['atm_num','res_num','res_name','sp1','atm_nam','x','y','z','sp2','chn_id','sp3','res_num2','holder']
        # self.comments = []
        # self.crd_contents = []

    def read_crd_file(self):
        inFile = open(self.CRD_file, 'r')
        contents = inFile.read()
        self.comments = []
        self.crd_contents = []
        self.col_types = [int,int,str,str,str,float,float,float,str,str,str,int,float]
        self.col_length = [(0,5),(5,11),(11,15),(15,16),(16,20),(20,30),(30,40),(40,50),(50,51),(51,55),(55,56),(56,60),(60,72)]
        self.col_names = ['atm_num','res_num','res_name','sp1','atm_nam','x','y','z','sp2','chn_id','sp3','res_num2','holder']
        #  the next loop could probably be simpler, but I made sure it can handle atoms with 
        atom_number = 0
        number_atoms = 0
        for i in contents.split('\n')[:-1]:
            if i[0] == '*':
                self.comments.append(i)
            elif len(i) == 5:
                number_atoms = int(i.strip())
                self.xyz = np.zeros((number_atoms,3))
            else:
                line_contents = []
                field_count = 0
                for j in self.col_length:
                    val = self.col_types[field_count]((i[j[0]:j[1]]).strip())
                    line_contents.append(val)
                    if self.col_names[field_count] == 'x':
                        self.xyz[atom_number,0] = val
                    elif self.col_names[field_count] == 'y':
                        self.xyz[atom_number,1] = val
                    elif self.col_names[field_count] == 'z':
                        self.xyz[atom_number,2] = val
                    field_count += 1
                self.crd_contents.append(line_contents)
                atom_number += 1
        return self.crd_contents
    # TODO check /home/noel/Datasets_nobackup/Constructs/Insulin_n_Barnase_barstar_2/proximity_optimization.py
    # for a corrected line to out put PDB directly from coordinates.
    # also check that lines do not come pre-formatted and that is why they are
    # different from proximity_optimization
    def write_pdb_from_crd(self,basedir,filename,lines):
        ''' TODO: write pdbdirectly with all chains the right way.
        http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        COLUMNS        DATA  TYPE    FIELD        DEFINITION
        -------------------------------------------------------------------------------------
         1 -  6        Record name   "ATOM  "
         7 - 11        Integer       serial       Atom  serial number.
        13 - 16        Atom          name         Atom name.
        17             Character     altLoc       Alternate location indicator.
        18 - 20        Residue name  resName      Residue name.
        22             Character     chainID      Chain identifier.
        23 - 26        Integer       resSeq       Residue sequence number.
        27             AChar         iCode        Code for insertion of residues.
        31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        55 - 60        Real(6.2)     occupancy    Occupancy.
        61 - 66        Real(6.2)     tempFactor   Temperature  factor.
        77 - 78        LString(2)    element      Element symbol, right-justified.
        79 - 80        LString(2)    charge       Charge  on the atom.
        '''
        line = '{:6}{:>4} {:5}{:1}{:3} {:1}{:>4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:2}'
        pdb_converted_lines = []
        for i in lines:
            pdb_converted_lines.append(line.format('ATOM',i[0],i[4],'',i[2],i[9],int(i[1]),'',float(i[5]),float(i[6]),\
                                                   float(i[7]),float('0'),float('0'),i[9],''))
        if basedir == '':
            basedir = '.'
        if basedir[-1] != "/":
            basedir += "/"
        outFile = open(basedir+filename, 'w')
        for i in pdb_converted_lines:
            outFile.write(i + '\n')
        outFile.close()
        
    def write_crd_from_pdb(self,basedir,first,last=''):
        '''
         This fixes when there are multiple chains and the residue number gets
         reset at the beginning of each chain. With this fix, residue numbers
         will be renumbered
          1 - 5 Integer Atom no. sequential
         6 - 10 Integer ires Residue position from file beginning
        11 - 11 Spacevi 
        12 - 15 Achar resname Residue name
        16 - 16 Space
        17 - 20 Achar type Atom type, IUPAC name of atom left justified
        21 - 30 Real(10.5) x Orthogonal coordinates for X, Angstroms
        31 - 40 Real(10.5) y Orthogonal coordinates for Y, Angstroms
        41 - 50 Real(10.5) z Orthogonal coordinates for Z, Angstroms
        51 - 51 Space
        52 - 55 Achar segid Segment identifier
        56 - 56 Space
        57 - 60 Achar resid Residue number within the segment
        61 - 70 Real(10.5) Weighting array value
        '''
        line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
        if last == '':
            file_name = first.split('.')
            if len(file_name) != 2:
                print("Error: A full PDB file name and suffix is required. Example: 0.pdb, protein.pdb")
                print("       Please enter a correct file name. Script will exit now.")
                sys.exit(1)
            elif len(file_name) == 2:
                if file_name[1] != 'pdb':
                    print("Error: A correct PDB file name and suffix is required. Example: 0.pdb, protein.pdb")
                    print("       Please enter a correct file name with lower case suffix. Script will exit now.")
                    sys.exit(1)
                else:
                    structure = struct.PDBParser(QUIET=True).get_structure('pdb', basedir+str(first))
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
                    
                    outFile = open(basedir+str(file_name[0])+'.crd', 'w')
                    for j in lines:
                        outFile.write(j+'\n')
                    outFile.close()
            else:
                print("Error: A correct PDB file name and suffix is required. Example: 0.pdb, protein.pdb")
                print("       Please enter a correct file name. Script will exit now.")
                sys.exit(1)               
        else:
            if (type(int(first)) is int) and (type(int(last)) is int):
                for i in range(0,990):
                    structure = struct.PDBParser(QUIET=True).get_structure('pdb', basedir+str(i)+'.pdb')
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
                    
                    outFile = open(basedir+str(i)+'.crd', 'w')
                    for j in lines:
                        outFile.write(j+'\n')
                    outFile.close()

    # TODO this method coulfd be giving a wrong output.
    def write_crd(self,structure, path_n_filename, model=0):
        line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
        count_atoms = 1
        count_resid = 1
        lines = []
        for a in structure.get_models():
            if a.id == model:
                for j in a.get_atoms():
                    # The following uncomment lines seem to be needed only
                    # if the first residue is not 1. do we need them?
                    #if count_atoms == 1:
                    #    current_resid = j.get_parent().id[1]
                    # TODO could this cause a bug?
                    #if current_resid != j.get_parent().id[1]:
                    #    current_resid = j.get_parent().id[1]
                    #    count_resid += 1
                    lines.append(line.format(count_atoms, count_resid, \
                                             j.get_parent().resname, \
                                             j.get_full_id()[4][0], \
                                             j.get_coord()[0], \
                                             j.get_coord()[1], \
                                             j.get_coord()[2], \
                                             j.get_full_id()[2], \
                                             j.get_parent().id[1], 0.00))
            count_atoms += 1
        count_atoms -= 1    
        len(str(count_atoms))
        lines.insert(0,(5-len(str(count_atoms)))*' '+str(count_atoms))
        lines.insert(0,'*')
        lines.insert(0,'* Date: XX/XX/XXXX by Noel Carrascal')
        lines.insert(0,'* Converting PDB to CRD')
        
        outFile = open(path_n_filename, 'w')
        for j in lines:
            outFile.write(j+'\n')
        outFile.close()
class psf(object):
    # http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win-html/node24.html
    def __init__(self, psf_file):
        self.PSF_file = psf_file
        self.psf_lines = self.read_psf_file()

    def read_psf_file(self):
        inFile = open(self.PSF_file, 'r')
        contents = inFile.read()
        psf_contents = []
        atom_number = 0
        number_atoms = 0
        is_natom = False
        for i in contents.split('\n'):
            ii = i.strip().split()
            if '!NATOM' in ii:
                number_atoms = int(ii[0])
                is_natom = True
                continue
            if is_natom:
                psf_contents.append(ii)
                atom_number += 1
                if atom_number == number_atoms:
                    is_natom = False
        return psf_contents
