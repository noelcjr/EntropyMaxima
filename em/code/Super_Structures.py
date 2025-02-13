# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:18:10 2016
http://mmcif.wwpdb.org/pdbx-mmcif-home-page.html
@author: noel
"""
#from Bio.PDB import *
import sys
import Bio.PDB.MMCIF2Dict as mmcifD
from Bio.PDB.MMCIFParser import *
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.Atom import Atom
from em.code.structure import crd
import pandas as pd
import numpy as np
#from CHARMM_Parser import *
# TODO: Make this class an extended class of Bio.PDB that solves all those
# problems modifying values.

class index_super_structure(object):
    """ Deprecrated. Overkill. Too many indexes. Some of them can be handled
        with structural expressions. 
    """
    def __init__(self, ss):
        self.Sup_Struc = ss
        self.sorted_by_aaid = False
        self.Fixed_Prolines = False
    def generate_indexes_from_Full_Structure(self):
        self.aa_indx = []
        self.aaid_indx = []
        self.aaid_indx2 = []
        self.ent_id_indx = []
        self.ent_id_indx2 = []
        self.chn_indx = []
        self.chn_indx2 = []
        self.comp_indx = []
        self.comp = []
        current_aaid = 0
        current_ent_id = 0
        current_chn = ''
        current_comp = ''
        self.aaid_end = []
        self.ent_id_end = []
        self.chn_end = []
        self.comp_end = []
        aaid_cnt = 0
        ent_id_cnt = 0
        chn_cnt = 0
        comp_cnt = 0
        for i in range(self.Sup_Struc.shape[0]):
            if i == 0:
                self.aa_indx.append(self.Sup_Struc.loc[i,'aa']) 
                self.aaid_indx.append(i)
                self.ent_id_indx.append(i)
                self.ent_id_indx2.append(self.Sup_Struc.loc[i,'ent_id'])
                self.chn_indx.append(i)
                self.chn_indx2.append(self.Sup_Struc.loc[i,'chain'])
                self.comp_indx.append(i)
                self.comp.append(self.Sup_Struc.loc[i,'component'][0:4])
                self.aaid_indx2.append(self.Sup_Struc.loc[i,'aaid'])
                current_aaid = self.Sup_Struc.loc[i,'aaid']
                current_ent_id = self.Sup_Struc.loc[i,'ent_id']
                current_chn = self.Sup_Struc.loc[i,'chain']
                current_comp = self.Sup_Struc.loc[i,'component'][0:4]
                aaid_cnt = 0
                ent_id_cnt = 0
                chn_cnt = 0
                comp_cnt = 0
            else:
                if self.Sup_Struc.loc[i,'aaid'] != current_aaid:
                    self.aaid_indx.append(i)
                    current_aaid = self.Sup_Struc.loc[i,'aaid']
                    self.aaid_end.append(aaid_cnt)
                    aaid_cnt += 1
                else:
                    aaid_cnt += 1
                if self.Sup_Struc.loc[i,'ent_id'] != current_ent_id:
                    self.ent_id_indx.append(i)
                    current_ent_id = self.Sup_Struc.loc[i,'ent_id']
                    self.ent_id_end.append(ent_id_cnt)
                    ent_id_cnt += 1
                else:
                    ent_id_cnt += 1
                if self.Sup_Struc.loc[i,'chain'] != current_chn:
                    self.chn_indx.append(i)
                    current_chn = self.Sup_Struc.loc[i,'chain']
                    self.chn_end.append(chn_cnt)
                    chn_cnt += 1
                else:
                    chn_cnt += 1
                if self.Sup_Struc.loc[i,'component'][0:4] != current_comp:
                    self.comp_indx.append(i)
                    self.comp.append(self.Sup_Struc.loc[i,'component'][0:4])
                    self.aaid_indx2.append(self.Sup_Struc.loc[i,'aaid'])
                    current_comp = self.Sup_Struc.loc[i,'component'][0:4]
                    self.aa_indx.append(self.Sup_Struc.loc[i,'aa'])
                    self.chn_indx2.append(self.Sup_Struc.loc[i,'chain'])
                    self.ent_id_indx2.append(self.Sup_Struc.loc[i,'ent_id'])
                    self.comp_end.append(comp_cnt)
                    comp_cnt += 1
                else:
                    comp_cnt += 1
        self.aaid_end.append(aaid_cnt)
        self.ent_id_end.append(ent_id_cnt)
        self.chn_end.append(chn_cnt)
        self.comp_end.append(comp_cnt)
        for i in range(len(self.aaid_indx)): 
            self.aaid_indx[i] = (self.aaid_indx[i],self.aaid_end[i])
        for i in range(len(self.ent_id_indx)):
            self.ent_id_indx[i] = (self.ent_id_indx[i],self.ent_id_end[i])
        for i in range(len(self.chn_indx)):
            self.chn_indx[i] = (self.chn_indx[i],self.chn_end[i])
        for i in range(len(self.comp_indx)):
            #self.comp_indx[i] = (self.comp_indx[i],self.comp_end[i])
            #self.comp_indx[i] = (self.aa_indx[i],self.aaid_indx[i],self.comp[i],self.comp_indx[i],self.comp_end[i])
            self.comp_indx[i] = (self.chn_indx2[i],self.ent_id_indx2[i],self.aaid_indx2[i],self.aa_indx[i],self.comp[i],self.comp_indx[i],self.comp_end[i])
    # TODO Re arranging atom positions in the Super Structure Full_Structure is something that is done for component
    # analysis only becuase in CHARMM the order of atoms might be an issue, but I need to check.
    def fix_proline_components(self):
        if not self.Fixed_Prolines:
            idx_ss_length = len(self.comp_indx)
            count_idx_ss = 0
            while count_idx_ss < idx_ss_length:
                if self.comp_indx[count_idx_ss][3] == 'PRO':
                    pass
    def sort_comp_index_by_aaid_within_chains(self):
        ''' This to make the output match the output by the Graduate School Perl program
        to verify that the energies calculated are equal
        '''
        if not self.sorted_by_aaid:
            idx_ss_length = len(self.comp_indx)
            for i in range(0,idx_ss_length-1):
                for j in range(i+1,idx_ss_length):
                    if (self.comp_indx[i][2] != self.comp_indx[j][2]) and \
                       (self.comp_indx[i][0] == self.comp_indx[j][0]):
                        if (str(self.comp_indx[i][2]) > str(self.comp_indx[j][2])):
                            temp = self.comp_indx[i]
                            self.comp_indx[i] = self.comp_indx[j]
                            self.comp_indx[j] = temp
            # Bubble sort above just order the tuples by aaid, but it also alters the AMIN-SIDE-CARB order of components
            # which is important for consistency. The following while loop adjusts the order of components.
            count_idx_ss = 0
            while count_idx_ss < idx_ss_length:
                if self.comp_indx[count_idx_ss][3] == 'GLY':
                    if self.comp_indx[count_idx_ss][4] == 'AMIN':
                        pass
                    else:
                        temp = self.comp_indx[count_idx_ss]
                        self.comp_indx[count_idx_ss] = self.comp_indx[count_idx_ss+1]
                        self.comp_indx[count_idx_ss+1] = temp
                    count_idx_ss += 2
                else:
                    if self.comp_indx[count_idx_ss][4] == 'AMIN':
                        if self.comp_indx[count_idx_ss+1][4] == 'SIDE':
                            pass
                        elif self.comp_indx[count_idx_ss+1][4] == 'CARB':
                            temp = self.comp_indx[count_idx_ss+1]
                            self.comp_indx[count_idx_ss+1] = self.comp_indx[count_idx_ss+2]
                            self.comp_indx[count_idx_ss+2] = temp
                    elif self.comp_indx[count_idx_ss][4] == 'CARB':
                        if self.comp_indx[count_idx_ss+1][4] == 'SIDE':
                            temp = self.comp_indx[count_idx_ss]
                            self.comp_indx[count_idx_ss] = self.comp_indx[count_idx_ss+2]
                            self.comp_indx[count_idx_ss+2] = temp
                        elif self.comp_indx[count_idx_ss+1][4] == 'AMIN':
                            temp = self.comp_indx[count_idx_ss]
                            self.comp_indx[count_idx_ss] = self.comp_indx[count_idx_ss+1]
                            self.comp_indx[count_idx_ss+1] = temp
                            temp = self.comp_indx[count_idx_ss+1]
                            self.comp_indx[count_idx_ss+1] = self.comp_indx[count_idx_ss+2]
                            self.comp_indx[count_idx_ss+2] = temp
                    elif self.comp_indx[count_idx_ss][4] == 'SIDE':
                        if self.comp_indx[count_idx_ss+1][4] == 'AMIN':
                            temp = self.comp_indx[count_idx_ss]
                            self.comp_indx[count_idx_ss] = self.comp_indx[count_idx_ss+1]
                            self.comp_indx[count_idx_ss+1] = temp
                        elif self.comp_indx[count_idx_ss+1][4] == 'CARB':
                            temp = self.comp_indx[count_idx_ss]
                            self.comp_indx[count_idx_ss] = self.comp_indx[count_idx_ss+1]
                            self.comp_indx[count_idx_ss+1] = temp
                            temp = self.comp_indx[count_idx_ss]
                            self.comp_indx[count_idx_ss] = self.comp_indx[count_idx_ss+2]
                            self.comp_indx[count_idx_ss+2] = temp
                    count_idx_ss += 3
            self.sorted_by_aaid = True

class Super_Structure(object):
    correct = {'ALA':{'H':'HN'},\
               'ARG':{'H':'HN','HB3':'HB1','HG3':'HG1','HD3':'HD1'},\
               'ASN':{'H':'HN','HB3':'HB1'},\
               'ASP':{'H':'HN'},\
               'CYS':{'H':'HN','HB3':'HB1'},\
               'GLU':{'H':'HN','HB3':'HB1','HG3':'HG1'},\
               'GLN':{'H':'HN','HB3':'HB1','HG3':'HG1'},\
               'GLY':{'H':'HN','HA3':'HA1'},\
               'HSD':{'H':'HN','HB3':'HB1'},\
               'HSE':{'H':'HN'},\
               'HSP':{'H':'HN'},\
               'ILE':{'CD1':'CD','H':'HN','HG13':'HG11','HD11':'HD1','HD12':'HD2',\
                      'HD13':'HD3','HB3':'HB1'},\
               'LEU':{'H':'HN','HB3':'HB1'},\
               'LYS':{'H':'HN','HB3':'HB1','HG3':'HG1','HD3':'HD1','HE3':'HE1'},\
               'MET':{'H':'HN'},\
               'PHE':{'H':'HN','HB3':'HB1'},\
               'PRO':{'H':'HN','HB3':'HB1','HG3':'HG1','HD3':'HD1'},\
               'SER':{'H':'HN','HB3':'HB1','HG':'HG1'},\
               'THR':{'H':'HN'},\
               'TRP':{'H':'HN'},\
               'TYR':{'H':'HN','HB3':'HB1'},\
               'VAL':{'H':'HN'}}
    inv_correct = {'ALA':{'HN':'H'},\
                   'ARG':{'HN':'H','HB1':'HB3','HG1':'HG3','HD1':'HD3','OT1':'O'},\
                   'ASN':{'HN':'H','HB1':'HB3','OT1':'O'},\
                   'ASP':{'HN':'H','OT1':'O'},\
                   'CYS':{'HN':'H','HB1':'HB3','OT1':'O'},\
                   'GLU':{'HN':'H','HB1':'HB3','HG1':'HG3','OT1':'O'},\
                   'GLN':{'HN':'H','HB1':'HB3','HG1':'HG3','OT1':'O'},\
                   'GLY':{'HN':'H','HA1':'HA3','OT1':'O'},\
                   'HSD':{'HN':'H','HB1':'HB3','OT1':'O'},\
                   'HSE':{'HN':'H','OT1':'O'},\
                   'HSP':{'HN':'H','OT1':'O'},\
                   'ILE':{'CD':'CD1','HN':'H','HG11':'HG13','HD1':'HD11','HD2':'HD12',\
                          'HD3':'HD13','HB1':'HB3','OT1':'O'},\
                   'LEU':{'HN':'H','HB1':'HB3','OT1':'O'},\
                   'LYS':{'HN':'H','HB1':'HB3','HG1':'HG3','HD1':'HD3','HE1':'HE3','OT1':'O'},\
                   'MET':{'HN':'H','OT1':'O'},\
                   'PHE':{'HN':'H','HB1':'HB3','OT1':'O'},\
                   'PRO':{'HN':'H','HB1':'HB3','HG1':'HG3','HD1':'HD3','OT1':'O'},\
                   'SER':{'HN':'H','HB1':'HB3','HG1':'HG','OT1':'O'},\
                   'THR':{'HN':'H','OT1':'O'},\
                   'TRP':{'HN':'H','OT1':'O'},\
                   'TYR':{'HN':'H','HB1':'HB3','OT1':'O'},\
                   'VAL':{'HN':'H','OT1':'O'}}
    his = {0:{0:'HSD',1:'HSE',2:'HSE'},1:{0:'HSD',1:'HSD',2:'HSE'},2:{0:'HSD',1:'HSD',2:'HSP'}}
    def __init__(self,parameter_object,structure_path,mode):
        if mode == 'setup':
            # TODO: for insulin 2hiu, the Histadines are not corrected to CHARMM's 
            #       three options.
            self.params = parameter_object
            self.parser = MMCIFParser()
            self.strctr = self.parser.get_structure('Center',structure_path)
            self.header = mmcifD.MMCIF2Dict(structure_path)
            self.Full_Structure = pd.DataFrame()
            # For adding peptides and their anchers
            self.pep = {}
            # NTERM is when building missing residues in the direction of the N-Term
            self.NTERM_anchor = {}
            # CTERM is when building missing residues in the direction of the C-Term
            self.CTERM_anchor = {}
            self.contiguous_LL = []
            self.hist = {}
        elif mode == 'add_linker':
            self.params = parameter_object
            self.Full_Structure = pd.read_csv(structure_path)
            # When csv comes from del_residue use first if.
            # when it comes from another add_residue use second if. better way?
            if "Unnamed: 0" in self.Full_Structure.columns:
                self.Full_Structure.set_index(['Unnamed: 0'],inplace=True)
                self.Full_Structure.index.names = ['indx']
            elif "indx" in self.Full_Structure.columns:
                self.Full_Structure.set_index(['indx'],inplace=True)
            # For adding peptides and their anchers
            self.pep = {}
            # NTERM is when building missing residues in the direction of the N-Term
            self.NTERM_anchor = {}
            # CTERM is when building missing residues in the direction of the C-Term
            self.CTERM_anchor = {}
            self.contiguous_LL = []
            self.hist = {}
        elif mode == 'charmm_input':
            self.params = parameter_object
            # directory is just for holding pth for output of CSV Super Structure.
            self.directory = structure_path
            self.Full_Structure = pd.DataFrame()
        else:
            print('ERROR: Super Structure can only be ran witht the following options:')
            print('       setup       to convert a cif to a csv or pdb file.')
            print('       add_linker  to add a linker and search its conformational space.')
            print('Try again. Exiting Now!')
            sys.exit(1)

    def check_models(self):
        '''CIF files have multiple models sometimes. It is better to check that they have the same number of atoms than
        assume that they do.'''
        for k in range(len(self.models)):
            for l in range(1,len(self.models)-1):
                a = list(self.atom_site_df.label_atom_id[self.atom_site_df.pdbx_PDB_model_num == k])
                b = list(self.atom_site_df.label_atom_id[self.atom_site_df.pdbx_PDB_model_num == l])
                if len(a) != len(b):
                    print('ERROR: Models '+k+' and '+l+' have different number of atoms.')
                    print('       Modify the code to account for this. This structure')
                    print('       is not processed. Exit Now!')
                    sys.exit(1)
                else:
                    for m in range(len(a)):
                        if a[m] != b[m]:
                            print('ERROR: Models '+k+' and '+l+' have differen atoms that are different.')
                            print('              '+a[m]+' '+b[m]+' do not match.')
                            print('       Modify the code to account for this. This structure')
                            print('       is not processed. Exit Now!')
                            sys.exit(1)
        print('All Models have identical atom types, and in the same order.')

    def corrections(self,res,atm):
        if res in self.correct:
            if atm in self.correct[res]: 
                return self.correct[res][atm]
            else:
                return atm
        else:
            return atm
    
    def inv_corrections(self,res,atm):
        if res in self.inv_correct:
            if atm in self.inv_correct[res]: 
                return self.inv_correct[res][atm]
            else:
                return atm
        else:
            return atm

    def create_super_structure_df_from_CRD_PSF(self, crd_file, psf_file):
        # CRD and PSF are list of strings corresponding to every line in the crd and psf files
        aa = []
        aaid = []
        entity_id = []
        chain_id = []
        atmtyp1 = []
        atmtyp2 = []
        charg = []
        component = []
        x = []
        y = []
        z = []
        occupancy = []
        B_iso = []
        ent = 1
        nuc = ['GUA','ADE','CYT','THY','URA']
        pro = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HSE','HSD','HSP','ILE','LEU','LYS','MET','PHE','PRO','SER',\
               'THR','TRP','TYR','VAL','TIP3','SOD','CLA']
         
        if len(crd_file.crd_lines) == len(psf_file.psf_lines):
            for i in range(len(crd_file.crd_lines)):
                crd_lne = crd_file.crd_lines[i]
                psf_lne = psf_file.psf_lines[i]
                if i == 0:
                    ent = 1
                    chn = psf_lne[1]
                else:
                    if chn != psf_lne[1]:
                        ent += 1
                        chn = psf_lne[1]
                aa.append(crd_lne[2])
                #aaid.append(crd_lne[1])
                # TODO check that this change won't affect MMGBSA calculations, or of it fixes the bug.
                # print("vvv",crd_lne[0]," v ",crd_lne[1]," w ",crd_lne[2]," x ",crd_lne[3]," y ",crd_lne[4]," z ",crd_lne[5],crd_lne[6],crd_lne[7])
                aaid.append(crd_lne[8])
                entity_id.append(ent)
                chain_id.append(psf_lne[1])
                atmtyp1.append(psf_lne[4])
                atmtyp2.append(psf_lne[5])
                charg.append(float(psf_lne[6]))
                #component.append('NA')
                x.append(crd_lne[5])
                y.append(crd_lne[6])
                z.append(crd_lne[7])
                occupancy.append(psf_lne[8])
                B_iso.append(psf_lne[9])
            print(crd_lne)
        self.Full_Structure['aa'] = pd.Series(aa)
        self.Full_Structure['aaid'] = pd.Series(aaid)
        self.Full_Structure['ent_id'] = pd.Series(entity_id)
        self.Full_Structure['chain'] = pd.Series(chain_id)
        self.Full_Structure['atmtyp1'] = pd.Series(atmtyp1)
        self.Full_Structure['atmtyp2'] = pd.Series(atmtyp2)
        c_cter = {}
        for j in range(len(self.Full_Structure['atmtyp1'])):
            resi = self.Full_Structure.loc[j,'aa']
            atmo = self.Full_Structure.loc[j,'atmtyp1']
            chan = self.Full_Structure.loc[j,'chain']
            if resi in nuc:
                count = 1
                found_atm = False
                for k in self.params.AA[resi].atoms:
                    if atmo in k:
                        found_atm = True
                        if count == 1:
                            component.append('NUC1')
                        elif count == len(self.params.AA[resi].atoms):
                            component.append('NUC5')
                        else:
                            component.append(('NUC'+str(count)))
                    count += 1
                if not found_atm:
                    if atmo == 'H5T':
                        component.append('NUC1')
                    elif atmo == 'H3T':
                        component.append('NUC5')
            elif resi in pro:
                count = 1
                found_atm = False
                for k in self.params.AA[resi].atoms:
                    if atmo in k:
                        found_atm = True
                        if count == 1:
                            component.append('AMINO')
                        elif count == len(self.params.AA[resi].atoms):
                            component.append('CARBO')
                        else:
                            component.append(('SIDE'+str(count)))
                    count += 1  
                if not found_atm:
                    count = 1
                    for k in self.params.AA['ACE'].atoms:
                        if atmo in k:
                            found_atm = True
                            component.append(('ACE'+str(count)))
                        count += 1
                    if not found_atm:
                        count = 1
                        for k in self.params.AA['NTER'].atoms:
                            if atmo in k:
                                found_atm = True
                                component.append(('NTER'+str(count)))
                            count += 1
                    if not found_atm:
                        count = 1
                        for k in self.params.AA['CTER'].atoms:
                            if atmo in k:
                                found_atm = True
                                component.append(('CTER'+str(count)))
                                if atmo == 'OT1':
                                    c_cter[chan] = j
                            count += 1
            else:
                print('ERROR: Amino Acid '+resi+' not found in parameters. Exit now.')
                sys.exit(1)
        self.Full_Structure['component'] = pd.Series(component)
        # Make C part of CTER1
        for j in c_cter:
            self.Full_Structure.loc[c_cter[j]-1,'component'] = 'CTER1'
        self.Full_Structure['charg'] = pd.Series(charg)
        mass = []
        atmNum = []
        atmtyp3 = []
        epsilon = []
        rmin_half = []
        atminfo = []
        for i in self.Full_Structure['atmtyp2']:
            atmNum.append(self.params.am.MASS[i][0])
            mass.append(self.params.am.MASS[i][1])
            atmtyp3.append(self.params.am.MASS[i][2])
            epsilon.append(self.params.NONBONDED[i][1])
            rmin_half.append(self.params.NONBONDED[i][2])
            atminfo.append(False)
        #print("AAA",atmNum[-1])
        self.Full_Structure['epsilon'] = pd.Series(epsilon)
        self.Full_Structure['rmin_half'] = pd.Series(rmin_half)
        self.Full_Structure['atmtyp3'] = pd.Series(atmtyp3)
        self.Full_Structure['mass'] = pd.Series(mass)
        self.Full_Structure['atmNum'] = pd.Series(atmNum)
        #print("BBB",self.Full_Structure['aaid'])
        #print("BBBB",self.Full_Structure.dtypes)
        self.Full_Structure['aaid'] = pd.to_numeric(self.Full_Structure['aaid'])
        #print("CCC",self.Full_Structure.dtypes)
        self.Full_Structure['ent_id'] = self.Full_Structure['ent_id'].apply(int)
        self.Full_Structure['mass'] = self.Full_Structure['mass'].apply(float)
        #atom_site_df['id'] = atom_site_df['id'].apply(int)
        self.Full_Structure['epsilon'] = self.Full_Structure['epsilon'].apply(float)
        self.Full_Structure['rmin_half'] = self.Full_Structure['rmin_half'].apply(float)
        self.Full_Structure['atmNum'] = self.Full_Structure['atmNum'].apply(int)
        self.Full_Structure['aainfo1'] = pd.Series(atminfo)
        self.Full_Structure['atminfo1'] = pd.Series(atminfo)
        self.Full_Structure['x1'] = pd.Series(x)
        self.Full_Structure['y1'] = pd.Series(y)
        self.Full_Structure['z1'] = pd.Series(z)   
        self.Full_Structure['occupancy'] = pd.Series(occupancy)
        self.Full_Structure['B_factor'] = pd.Series(B_iso)
        # Makes sure that some columns are the right type
        #self.Full_Structure['x1'] = self.Full_Structure['x1'].apply(float)
        #self.Full_Structure['y1'] = self.Full_Structure['y1'].apply(float)
        #self.Full_Structure['z1'] = self.Full_Structure['z1'].apply(float)
        #self.Full_Structure['charg'] = self.Full_Structure['charg'].apply(float)
        #print("BBBB",self.Full_Structure.dtypes)
    
    def create_column_with_ztranslated(self, new_column_name, from_column, chains_translated, distance_translated):
        z_d = []
        displaced_chains = chains_translated.split(',')
        for i in self.Full_Structure.index:
            if self.Full_Structure.loc[i,'chain'] in displaced_chains:
                z_d.append(self.Full_Structure.loc[i,from_column]+distance_translated)
            else:
                z_d.append(self.Full_Structure.loc[i,from_column])
        self.Full_Structure[new_column_name] = pd.Series(z_d)
        self.Full_Structure[new_column_name] = self.Full_Structure[new_column_name].apply(float)
    
    def add_column_to_super_structure_df_from_CRD_PSF(self, column_name, from_file):
        if column_name in ['sa_z0','sa_z500']:
            temp = np.zeros(self.Full_Structure.shape[0])
            for i in open(from_file, 'r').read().split('\n'):
                ii = i.strip().split()
                if len(ii) > 0:
                    temp[int(ii[0])-1] = float(ii[1])
            self.Full_Structure[column_name] = pd.Series(temp)
        else:
            temp = []
            for i in open(from_file, 'r').read().split('\n'):
                ii = i.strip().split()
                if len(ii) > 0:
                    temp.append(ii[1])
            self.Full_Structure[column_name] = pd.Series(temp)
        self.Full_Structure[column_name] = self.Full_Structure[column_name].apply(float)

    def get_histidine_info(self,ent):
        #for i in entity_poly_df.entity_id:
        self.hist = {}
        hist_entity = {}
        for i in ent.id:
            hist_chain = {}
            for chain in str(self.entity_poly_df.pdbx_strand_id[self.entity_poly_df.entity_id == i][int(i)-1]).split(","):
                resid = 1 
                for j in self.entity_poly_seq_df.index[self.entity_poly_seq_df.entity_id == i]:
                    if self.entity_poly_seq_df.loc[j,'mon_id']  == 'HIS':
                        temp = list(self.atom_site_df.label_atom_id[(self.atom_site_df.pdbx_PDB_model_num == '1') & \
                                                                    (self.atom_site_df.label_entity_id == i) & \
                                                                    (self.atom_site_df.label_seq_id == self.entity_poly_seq_df.loc[j,'num'])])
                        HDcount = 0
                        HEcount = 0
                        for k in temp:
                            if k[0] == 'H':
                                if len(k) > 1:
                                    if k[1] == 'E':
                                        HEcount += 1
                                    elif k[1] == 'D':
                                        HDcount += 1
                        hist_chain[resid] = (resid,HDcount,HEcount)
                    resid += 1
                hist_entity[chain] = hist_chain
            self.hist[i] = hist_entity
    def create_super_structure_df(self):
        # Find HIS and analyze for the next step to have the HIS ready to identify
        ent = self.entity_df[['id','pdbx_number_of_molecules']][self.entity_df.type == 'polymer']
        self.get_histidine_info(ent)
        ###########################################################################
        # create a Data Frame from sequence information. Atom types and
        # charges are added to the sequence information.
        aa = []
        aaid = []
        entity_id = []
        chain_id = []
        atmtyp1 = []
        atmtyp2 = []
        charg = []
        component = []
        add_ACE = {}  # TODO this is boolean not Dictionary
        add_CTR = {}  # TODO this is boolean not Dictionary
        beg_seq = {}
        end_seq = {}  
        # It is possible to have more than one chain per entity.
        for entid in ent.id:
            for chain in str(self.entity_poly_df.pdbx_strand_id[self.entity_poly_df.entity_id == entid][int(entid)-1]).split(","):
                beg_seq[chain] = int(self.struct_ref_seq_df.seq_align_beg[self.struct_ref_seq_df.pdbx_strand_id == chain])
                end_seq[chain] = int(self.struct_ref_seq_df.seq_align_end[self.struct_ref_seq_df.pdbx_strand_id == chain])
                add_ACE = True
                add_CTR = False
                # ASSUME ACE AND CTER are always added, unless ACE or CTER are explicitly present.
                # Assign parameters to all atoms with coordinates
                # include atoms without coordinate that could be missing in cystal structure
                # According to residue type and CHARMM top27 force filed.                  
                for j in self.entity_poly_seq_df.index[self.entity_poly_seq_df.entity_id == entid]:
                    # HIS is a special case that have three variants in CHARMM
                    if self.entity_poly_seq_df.loc[j,'mon_id'] == 'HIS':
                        for k in self.hist[entid]:
                            for l in self.hist[entid][k]:
                                HD = self.hist[entid][k][l][1]
                                HE = self.hist[entid][k][l][2]
                                if (HD <= 2) and (HE <= 2):
                                    # the dictionary simplifies the algorithm
                                    chrm = self.his[HD][HE]
                                else:
                                    print('Error: Problem with histadine assignment')
                                    print('#HD =',HD,' #HE=',HE)
                                    print('       Program exit.')
                                    sys.exit(1)
                        pdbx = 'HIS'
                    else:
                        chrm = self.entity_poly_seq_df.ix[j]['mon_id']
                        pdbx = self.entity_poly_seq_df.ix[j]['mon_id']
                    # Set a flag that tells if a CYS is forming a disulfide bond
                    is_disulfide = False
                    if self.entity_poly_seq_df.loc[j,'mon_id'] == 'CYS':
                        for k in self.struct_conn_df.index:
                            if self.struct_conn_df.loc[k,'conn_type_id'] == 'disulf':
                                if self.struct_conn_df.loc[k,'ptnr1_label_asym_id'] == chain:
                                    if self.struct_conn_df.loc[k,'ptnr1_label_comp_id'] == 'CYS':
                                        if self.struct_conn_df.loc[k,'ptnr1_label_seq_id'] == self.entity_poly_seq_df.loc[j,'num']:
                                            is_disulfide = True
                                if self.struct_conn_df.loc[k,'ptnr2_label_asym_id'] == chain:
                                    if self.struct_conn_df.loc[k,'ptnr2_label_comp_id'] == 'CYS':
                                        if self.struct_conn_df.loc[k,'ptnr2_label_seq_id'] == self.entity_poly_seq_df.loc[j,'num']:
                                            is_disulfide = True
                    # TODO, the program is designed to detect disulfide bonds from information in the CIF file.
                    # The line following this comment will ignore the presence of disulfide bonds because they
                    # are dependent of reducing or oxidicing environment. Whether a HG1 atom is present
                    # or not should be based on design desitions for protein. Insulin forms in the cytoplasm
                    # and it should be reduce for accurate calculations.
                    is_disulfide = False
                    snum = self.entity_poly_seq_df.loc[j,'num']
                    ent_id = self.entity_poly_seq_df.loc[j,'entity_id']
                    # print(entid,chain,j,chrm,pdbx,snum,ent_id)
                    # A peptide is a zwitterion at near neutral pH, protonated at low pH 
                    # (NH3+, CO2H), or deprotonated at high pH (NH2, CO2-). For a peptide 
                    # with uncharged ends, common patches to apply (via generate) are ACE 
                    # (NYacetyl) and CT3 (C-terminal methylamine);
                    # We will modify all crytal structure initial and last resude to have
                    # N-terminal ACE and C-terminal CTER.
                    if int(snum) == end_seq[chain]:
                        add_CTR = True
                    if chrm in self.params.AA:
                        if chrm == 'ACE':
                            # THis assuemes that ACE is always first in the sequence
                            for k in self.params.AA['ACE'].atoms:
                                for l in k:
                                    aa.append(list(self.entity_poly_seq_df.mon_id[(self.entity_poly_seq_df.entity_id == entid) & \
                                                                                  (self.entity_poly_seq_df.num == str(beg_seq[chain]))])[0])
                                    aaid.append(snum)
                                    entity_id.append(ent_id)
                                    chain_id.append(chain)
                                    atmtyp1.append(l)
                                    atmtyp2.append(self.params.AA['ACE'].atom_type[l])
                                    charg.append(self.params.AA['ACE'].atom_chrg[l])
                                    component.append('ACETY')
                            add_ACE = False
                        else:
                            if add_ACE:
                                for k in self.params.AA['ACE'].atoms:
                                    for l in k:
                                        aa.append(pdbx)
                                        aaid.append(snum)
                                        entity_id.append(ent_id)
                                        chain_id.append(chain)
                                        atmtyp1.append(l)
                                        atmtyp2.append(self.params.AA['ACE'].atom_type[l])
                                        charg.append(self.params.AA['ACE'].atom_chrg[l])
                                        component.append('ACETY')
                                add_ACE = False
                                comp = 1
                                for k in self.params.AA[chrm].atoms:
                                    for l in k:
                                        if l == 'HG1' and chrm == 'CYS' and is_disulfide:
                                            pass
                                        else:
                                            aa.append(pdbx)
                                            aaid.append(snum)
                                            entity_id.append(ent_id)
                                            chain_id.append(chain)
                                            atmtyp1.append(self.corrections(chrm,l))
                                            atmtyp2.append(self.params.AA[chrm].atom_type[self.corrections(chrm,l)])
                                            charg.append(self.params.AA[chrm].atom_chrg[self.corrections(chrm,l)])
                                            if comp == 1:
                                                component.append('AMINO')
                                            else:
                                                if l in ['C','O']:
                                                    component.append('CARBO')
                                                else:
                                                    component.append(('SIDE'+str(comp)))
                                    comp += 1
                            else:
                                if add_CTR:
                                    comp = 1
                                    for k in self.params.AA[chrm].atoms:
                                        for l in k:
                                            if l == 'HG1' and chrm == 'CYS' and is_disulfide:
                                                pass
                                            else:
                                                aa.append(pdbx)
                                                aaid.append(snum)
                                                entity_id.append(ent_id)
                                                chain_id.append(chain)
                                                if l == 'O':
                                                    atmtyp1.append('OT1')
                                                    atmtyp2.append(self.params.AA['CTER'].atom_type['OT1'])
                                                    charg.append(self.params.AA['CTER'].atom_chrg['OT1'])
                                                    component.append('CTERM')
                                                elif l == 'C':
                                                    atmtyp1.append('C')
                                                    atmtyp2.append(self.params.AA['CTER'].atom_type['C'])
                                                    charg.append(self.params.AA['CTER'].atom_chrg['C'])
                                                    component.append('CTERM')
                                                else:
                                                    atmtyp1.append(self.corrections(chrm,l))
                                                    atmtyp2.append(self.params.AA[chrm].atom_type[self.corrections(chrm,l)])
                                                    charg.append(self.params.AA[chrm].atom_chrg[self.corrections(chrm,l)])
                                                    if comp == 1:
                                                        component.append('AMINO')
                                                    else:
                                                        if l in ['C','O']:
                                                            component.append('CARBO')
                                                        else:
                                                            component.append(('SIDE'+str(comp)))
                                        comp += 1
                                    aa.append(pdbx)
                                    aaid.append(snum)
                                    entity_id.append(ent_id)
                                    chain_id.append(chain)
                                    atmtyp1.append('OT2')
                                    atmtyp2.append(self.params.AA['CTER'].atom_type['OT2'])
                                    charg.append(self.params.AA['CTER'].atom_chrg['OT2'])
                                    component.append('CTERM')  
                                else:
                                    comp = 1
                                    for k in self.params.AA[chrm].atoms:
                                        for l in k:
                                            if l == 'HG1' and chrm == 'CYS' and is_disulfide:
                                                pass
                                            else:
                                                aa.append(pdbx)
                                                aaid.append(snum)
                                                entity_id.append(ent_id)
                                                chain_id.append(chain)
                                                atmtyp1.append(self.corrections(chrm,l))
                                                atmtyp2.append(self.params.AA[chrm].atom_type[self.corrections(chrm,l)])
                                                charg.append(self.params.AA[chrm].atom_chrg[self.corrections(chrm,l)])
                                                if comp == 1:
                                                    component.append('AMINO')
                                                else:
                                                    if l in ['C','O']:
                                                        component.append('CARBO')
                                                    else:
                                                        component.append(('SIDE'+str(comp)))
                                        comp += 1
        self.Full_Structure['aa'] = pd.Series(aa)
        self.Full_Structure['aaid'] = pd.Series(aaid)
        self.Full_Structure['ent_id'] = pd.Series(entity_id)
        self.Full_Structure['chain'] = pd.Series(chain_id)
        self.Full_Structure['atmtyp1'] = pd.Series(atmtyp1)
        self.Full_Structure['atmtyp2'] = pd.Series(atmtyp2)
        self.Full_Structure['component'] = pd.Series(component)
        self.Full_Structure['charg'] = pd.Series(charg)
        ###########################################################################
        # Add atomtyp, masses and atmNumber to each atom type
        mass = []
        atmNum = []
        atmtyp3 = []
        epsilon = []
        rmin_half = []
        atminfo = []
        for i in self.Full_Structure['atmtyp2']:
            atmNum.append(self.params.am.MASS[i][0])
            mass.append(self.params.am.MASS[i][1])
            atmtyp3.append(self.params.am.MASS[i][2])
            epsilon.append(self.params.NONBONDED[i][1])
            rmin_half.append(self.params.NONBONDED[i][2])
            atminfo.append(False)
        
        self.Full_Structure['epsilon'] = pd.Series(epsilon)
        self.Full_Structure['rmin_half'] = pd.Series(rmin_half)
        self.Full_Structure['atmtyp3'] = pd.Series(atmtyp3)
        self.Full_Structure['mass'] = pd.Series(mass)
        self.Full_Structure['atmNum'] = pd.Series(atmNum)
        ###########################################################################
        # DF Type correction.
        self.Full_Structure['aaid'] = self.Full_Structure['aaid'].apply(int)
        self.Full_Structure['ent_id'] = self.Full_Structure['ent_id'].apply(int)
        self.Full_Structure['mass'] = self.Full_Structure['mass'].apply(float)
        self.atom_site_df['id'] = self.atom_site_df['id'].apply(int)
        self.Full_Structure['epsilon'] = self.Full_Structure['epsilon'].apply(float)
        self.Full_Structure['rmin_half'] = self.Full_Structure['rmin_half'].apply(float)
        self.Full_Structure['atmNum'] = self.Full_Structure['atmNum'].apply(int)
        ###########################################################################
        # Distinguish between imputed and available structural information T or F 
        entities = list(set(self.Full_Structure.ent_id))
        for g in self.models:
            x = []
            y = []
            z = []
            occupancy = []
            B_iso = []
            for h in entities:
                for i in self.Full_Structure.index[self.Full_Structure.ent_id == h]:
                    aa1 = self.Full_Structure.ix[i]['aa']
                    type1 = self.inv_corrections(aa1,self.Full_Structure.ix[i]['atmtyp1'])
                    aaid1 = str(self.Full_Structure.ix[i]['aaid'])
                    chn = self.Full_Structure.ix[i]['chain']
                    cx = self.atom_site_df.Cartn_x[(self.atom_site_df.label_atom_id == type1) & \
                                                   (self.atom_site_df.label_comp_id == aa1) & \
                                                   (self.atom_site_df.label_seq_id == aaid1) & \
                                                   (self.atom_site_df.pdbx_PDB_model_num == g) & \
                                                   (self.atom_site_df.label_entity_id == str(h)) & \
                                                   (self.atom_site_df.label_asym_id == chn)]
                    if len(cx) == 1:
                        x.append(float(cx))
                    else:
                        x.append(float('nan'))
                    cy = self.atom_site_df.Cartn_y[(self.atom_site_df.label_atom_id == type1) & \
                                                   (self.atom_site_df.label_comp_id == aa1) & \
                                                   (self.atom_site_df.label_seq_id == aaid1) & \
                                                   (self.atom_site_df.pdbx_PDB_model_num == g) & \
                                                   (self.atom_site_df.label_entity_id == str(h)) & \
                                                   (self.atom_site_df.label_asym_id == chn)]
                    if len(cy) == 1:
                        y.append(float(cy))
                    else:
                        y.append(float('nan'))
                    cz = self.atom_site_df.Cartn_z[(self.atom_site_df.label_atom_id == type1) & \
                                                   (self.atom_site_df.label_comp_id == aa1) & \
                                                   (self.atom_site_df.label_seq_id == aaid1) & \
                                                   (self.atom_site_df.pdbx_PDB_model_num == g) & \
                                                   (self.atom_site_df.label_entity_id == str(h)) & \
                                                   (self.atom_site_df.label_asym_id == chn)]
                    if len(cz) == 1:
                        z.append(float(cz))
                    else:
                        z.append(float('nan'))
                    ocu = self.atom_site_df.occupancy[(self.atom_site_df.label_atom_id == type1) & \
                                                      (self.atom_site_df.label_comp_id == aa1) & \
                                                      (self.atom_site_df.label_seq_id == aaid1) & \
                                                      (self.atom_site_df.pdbx_PDB_model_num == g) & \
                                                      (self.atom_site_df.label_entity_id == str(h)) & \
                                                      (self.atom_site_df.label_asym_id == chn)]
                    if len(ocu) == 1:
                        occupancy.append(float(ocu))
                    else:
                        occupancy.append('nan')
                    bfc = self.atom_site_df.B_iso_or_equiv[(self.atom_site_df.label_atom_id == type1) & \
                                                           (self.atom_site_df.label_comp_id == aa1) & \
                                                           (self.atom_site_df.label_seq_id == aaid1) & \
                                                           (self.atom_site_df.pdbx_PDB_model_num == g) & \
                                                           (self.atom_site_df.label_entity_id == str(h)) & \
                                                           (self.atom_site_df.label_asym_id == chn)]
                    if len(bfc) == 1:
                        B_iso.append(float(bfc))
                    else:
                        B_iso.append(float('nan'))
            self.Full_Structure['aainfo'+g] = pd.Series(atminfo)
            self.Full_Structure['atminfo'+g] = pd.Series(atminfo)
            self.Full_Structure['x'+g] = pd.Series(x)
            self.Full_Structure['y'+g] = pd.Series(y)
            self.Full_Structure['z'+g] = pd.Series(z)
            self.Full_Structure['occupancy'] = pd.Series(occupancy)
            self.Full_Structure['B_factor'] = pd.Series(B_iso)
        ###########################################################################
        #self.Full_Structure['B_factor'] = self.Full_Structure['B_factor'].apply(float)
        #self.Full_Structure['occupancy'] = self.Full_Structure['occupancy'].apply(float)
        ###########################################################################
        # For each models, set to True for atom coordinates with structural info
        # from the CIF file. pd.notnull gives True is not null, false if null
        for h in self.models:
            x = pd.notnull(self.Full_Structure['x'+h])
            y = pd.notnull(self.Full_Structure['y'+h])
            z = pd.notnull(self.Full_Structure['z'+h])
            # This loop checks that an atom location is obtained from x-struct
            # or added from rebuilding missing peptides
            for i in range(self.Full_Structure.shape[0]):
                if x[i] & y[i] & z[i]:
                    self.Full_Structure.loc[i,'atminfo'+h] = True
            # This loop, aa with all atms added from rebuilding missing peptides
            # are considered Reabuilt (false). From x-structure (True)
            for i in self.entity_poly_df.entity_id:
                for j in self.entity_poly_seq_df.num[self.entity_poly_seq_df.entity_id == i]:
                    set_true = False
                    for k in self.Full_Structure['atminfo'+h][(self.Full_Structure.ent_id == int(i)) & \
                                                              (self.Full_Structure.aaid == int(j))]:
                        if k:
                            set_true = True
                    if set_true:
                        for l in self.Full_Structure.index[(self.Full_Structure.ent_id == int(i)) & \
                                                            (self.Full_Structure.aaid == int(j))]:
                            self.Full_Structure.loc[l,'aainfo'+h] = True

    def get_missing_aa_listlist(self,mdl,ent,chn):
        """ Maps amino acids that are present and missing and puts
            contiguous missing aminoacid in a list of lists"""
        self.missing = []
        self.present = []
        for j in self.entity_poly_seq_df.num[(self.entity_poly_seq_df.entity_id == ent)]:
            set_true = False
            for k in self.Full_Structure['atminfo'+mdl][(self.Full_Structure.ent_id == int(ent)) & \
                                                        (self.Full_Structure.aaid == int(j)) & \
                                                        (self.Full_Structure.chain == chn)]:
                if k:
                    set_true = True
            if not set_true:
                self.missing.append(list(self.Full_Structure.aaid[(self.Full_Structure.ent_id == int(ent)) & \
                                                                  (self.Full_Structure.aaid == int(j)) &\
                                                                  (self.Full_Structure.chain == chn)])[0])
            else:
                self.present.append(list(self.Full_Structure.aaid[(self.Full_Structure.ent_id == int(ent)) & \
                                                                  (self.Full_Structure.aaid == int(j)) &\
                                                                  (self.Full_Structure.chain == chn)])[0])
        temp = 0
        contiguous_temp = []
        self.contiguous_LL = []
        for i in range(len(self.missing)):
            if i == 0:
                temp = self.missing[i]
                contiguous_temp.append(temp)
            else:
                if self.missing[i] == (temp+1):
                    contiguous_temp.append(self.missing[i])
                    temp = self.missing[i]
                else:
                    self.contiguous_LL.append(contiguous_temp)
                    temp = self.missing[i]
                    contiguous_temp = []
                    contiguous_temp.append(self.missing[i])
        self.contiguous_LL.append(contiguous_temp)
        
    def build_missing_aa(self):
        for h in self.models:
            model = str(h)
            NTER_ent = {}
            CTER_ent = {}
            for hh in self.entity_poly_df.entity_id:
                entity_id = int(hh)
                NTER_ch = {}
                CTER_ch = {}
                NTER_ent[hh] = NTER_ch
                CTER_ent[hh] = CTER_ch
                # TODO would entity_id-1 be right all the time?
                for hhh in str(self.entity_poly_df.pdbx_strand_id[self.entity_poly_df.entity_id == hh][entity_id-1]).split(","):
                    self.get_missing_aa_listlist(model,hh,hhh)
                    CTER_ent[hh][hhh] = False
                    NTER_ent[hh][hhh] = False
                    # This two ifs check that the list of lists is not empty
                    if((len(self.contiguous_LL) != 1) and (len(self.contiguous_LL[0]) != 0)):
                        for i in range(len(self.contiguous_LL)):
                            if i == 0 and self.contiguous_LL[i][0] == 1:
                                for j in self.contiguous_LL[i][::-1]:
                                    place = self.entity_poly_seq_df['mon_id'][j-1]
                                    if place == 'HIS':
                                        place = self.his[self.hist[hh]['A'][j][1]][self.hist[hh]['A'][j][2]]
                                    if place == 'ACE':
                                        NTER_ent[hh][hhh] = True
                                    if place == 'CTER':
                                        print('ERROR: Structure attempts to place a CTER with a CTER_ancher.')
                                        print('       CTER has no CTERM_ancher by definition because CTER caps')
                                        print('       ,or ends, the chain. With nothing to ancher beyond that.')
                                        sys.exit(1)
                                    self.fit_coordinates('Ndir',j,entity_id,hhh,model,place)
                            else:
                                for j in self.contiguous_LL[i][:]:
                                    place = self.entity_poly_seq_df['mon_id'][j-1]
                                    if place == 'HIS':
                                        place = self.his[self.hist[hh]['A'][j][1]][self.hist[hh]['A'][j][2]]
                                    if place == 'CTER':
                                        CTER_ent[hh][hhh] = True
                                    if place == 'ACE':
                                        print('ERROR: Structure attempts to place a ACE with a NTER_ancher.')
                                        print('       ACE has no NTERM_ancher by definition because ACE caps')
                                        print('       ,or begins, the chain. With nothing to ancher before that.')
                                        sys.exit(1)
                                    self.fit_coordinates('Cdir',j,entity_id,hhh,model,place)
        # TODO This could be added to the above function somehow. I did it separately
        # because the mess of arrays made it hard to make together.
        for h in self.models:
            model = str(h)
            min_aaid = {}
            max_aaid = {}
            for hh in self.entity_poly_df.entity_id:
                entity_id = int(hh)
                min_ch = {}
                max_ch = {}
                min_aaid[hh] = min_ch
                max_aaid[hh] = max_ch
                for hhh in str(self.entity_poly_df.pdbx_strand_id[self.entity_poly_df.entity_id == hh][entity_id-1]).split(","):
                    for j in self.entity_poly_seq_df.num[(self.entity_poly_seq_df.entity_id == hh)]:
                        if int(j) == 1:
                            min_aaid[hh][hhh] = int(j)
                            max_aaid[hh][hhh] = int(j)
                        if j > max_aaid:
                            max_aaid[hh][hhh] = int(j)
            for hh in self.entity_poly_df.entity_id:
                entity_id = int(hh)
                for hhh in str(self.entity_poly_df.pdbx_strand_id[self.entity_poly_df.entity_id == hh][entity_id-1]).split(","):
                    if not NTER_ent[hh][hhh]:
                        self.fit_coordinates('NTER',min_aaid[hh][hhh],entity_id,hhh,model,'ACE')
                    if not CTER_ent[hh][hhh]:
                        self.fit_coordinates('CTER',max_aaid[hh][hhh],entity_id,hhh,model,'CTER') 
                                                
    def fit_coordinates(self,mod,aaid,entity_id,chain,model,place):
        '''mod tells this method how to deal with the residue numbers
        to which added aa are fixed. This residue number are modified by + or - 1'''
        if mod == 'Ndir':
            aaid2 = aaid + 1
        elif mod == 'Cdir':
            aaid2 = aaid - 1
        elif mod == 'NTER':
            aaid2 = aaid
        elif mod == 'CTER':
            aaid2 = aaid
        else:
            print('ERROR: mod variable for fit_coordinates does not exists.')
            print('       Only Ndir, Cdir, NTER and CTER are valid.')
            print('       Exiting program before finished.')
            sys.exit(1)
        fixed1 = []
        moving1 = []
        moving2 = []
        for k in self.Full_Structure.index[(self.Full_Structure.aaid == aaid2) &\
                                           (self.Full_Structure.ent_id == entity_id) &\
                                           (self.Full_Structure.chain == chain)]:
            temp = self.Full_Structure.loc[k,'atmtyp1']
            if temp == 'N' or temp == 'CA' or temp =='C':
                cord = np.array([self.Full_Structure.loc[k,'x'+model],\
                                 self.Full_Structure.loc[k,'y'+model],\
                                 self.Full_Structure.loc[k,'z'+model]],\
                                 dtype=float)
                fixed1.append(Atom(temp,cord,0.0,0.0,1,temp+'k',k))
        if (mod == 'Ndir') | (mod == 'NTER'):
            for l in self.CTERM_anchor[place]:
                if l == 'N':
                    tempN = Atom(self.CTERM_anchor[place][l].get_id(),\
                                 self.CTERM_anchor[place][l].get_coord(),\
                                 0.0,0.0,1,\
                                 self.CTERM_anchor[place][l].get_id(),\
                                 self.CTERM_anchor[place][l].get_id())
                elif l == 'CA':
                    tempCA = Atom(self.CTERM_anchor[place][l].get_id(),\
                                  self.CTERM_anchor[place][l].get_coord(),\
                                  0.0,0.0,1,\
                                  self.CTERM_anchor[place][l].get_id(),\
                                  self.CTERM_anchor[place][l].get_id())
                elif l =='C':
                    tempC = Atom(self.CTERM_anchor[place][l].get_id(),\
                                 self.CTERM_anchor[place][l].get_coord(),\
                                 0.0,0.0,1,\
                                 self.CTERM_anchor[place][l].get_id(),\
                                 self.CTERM_anchor[place][l].get_id())
        elif (mod == 'Cdir') | (mod == 'CTER'):
            for l in self.NTERM_anchor[place]:
                if l == 'N':
                    tempN = Atom(self.NTERM_anchor[place][l].get_id(),\
                                 self.NTERM_anchor[place][l].get_coord(),\
                                 0.0,0.0,1,\
                                 self.NTERM_anchor[place][l].get_id(),\
                                 self.NTERM_anchor[place][l].get_id())
                elif l == 'CA':
                    tempCA = Atom(self.NTERM_anchor[place][l].get_id(),\
                                  self.NTERM_anchor[place][l].get_coord(),\
                                  0.0,0.0,1,\
                                  self.NTERM_anchor[place][l].get_id(),\
                                  self.NTERM_anchor[place][l].get_id())
                elif l =='C':
                    tempC = Atom(self.NTERM_anchor[place][l].get_id(),\
                                 self.NTERM_anchor[place][l].get_coord(),\
                                 0.0,0.0,1,\
                                 self.NTERM_anchor[place][l].get_id(),\
                                 self.NTERM_anchor[place][l].get_id())
        else:
            print('Error: mod not found. Program will exit.')
            sys.exit(1)
        moving1.append(tempN)
        moving1.append(tempCA)
        moving1.append(tempC)
        for m in self.pep[place]:
            moving2.append(Atom(self.pep[place][m].get_id(),\
                                self.pep[place][m].get_coord(),\
                                0.0,0.0,1,\
                                self.pep[place][m].get_id(),\
                                self.pep[place][m].get_id()))
        sup = Superimposer()
        sup.set_atoms(fixed1,moving1)
        sup.apply(moving1)
        sup.apply(moving2)
        # TODO: It is here where we need to get some angles and place the added amino acid in the right conformation
        #       relative to the planar shape of the amino and carbonyl gropus. 
        moving3 = {}
        for k in self.Full_Structure.index[(self.Full_Structure.aaid == aaid) &\
                                           (self.Full_Structure.ent_id == entity_id) &\
                                           (self.Full_Structure.chain == chain)]:
            moving3[self.inv_corrections(place,self.Full_Structure.loc[k,'atmtyp1'])] = k
        for i in moving2:
            self.Full_Structure.loc[moving3[self.inv_corrections(place,i.get_id())],'x'+model] = i.get_coord()[0]
            self.Full_Structure.loc[moving3[self.inv_corrections(place,i.get_id())],'y'+model] = i.get_coord()[1]
            self.Full_Structure.loc[moving3[self.inv_corrections(place,i.get_id())],'z'+model] = i.get_coord()[2]
            self.Full_Structure.loc[moving3[self.inv_corrections(place,i.get_id())],'occupancy'] = 0.0
            self.Full_Structure.loc[moving3[self.inv_corrections(place,i.get_id())],'B_factor'] = 0.0

    def delete_aa(self,aaid,ent_id,chain,aa):
        del_atoms = []
        if aa in self.params.AA:
            for i in range(len(self.params.AA[aa].atoms)):
                if aa == 'CTER':
                    if 'OT2' in self.params.AA[aa].atoms[i]:
                        del_atoms = del_atoms + ['OT2']
                else:
                    del_atoms = del_atoms + self.params.AA[aa].atoms[i]
        else:            
            print('ERROR: Amino Acid or Terminal '+aa+' is not found in parameters for deletion.')
            print('       '+aa+' must be in the parameter filed read.')
            print('       Exiting Now!')
            sys.exit(1)
        for i in del_atoms:
            self.Full_Structure = self.Full_Structure[(self.Full_Structure.aaid != aaid) |\
                                                      (self.Full_Structure.chain != chain) |\
                                                      (self.Full_Structure.atmtyp1 != i)]
        # When the CTER is deleted, the group parmeters have to be change to CARBONYL from CTERM
        if aa == 'CTER':
            indx = self.Full_Structure.index[(self.Full_Structure.aaid == aaid) &\
                                             (self.Full_Structure.atmtyp1 == 'C') &\
                                             (self.Full_Structure.chain == chain)][0]
            self.Full_Structure.loc[indx,['atmtyp2']] = 'C'
            self.Full_Structure.loc[indx,['component']] = 'CARBO'
            self.Full_Structure.loc[indx,['charg']] = 0.51
            self.Full_Structure.loc[indx,['epsilon']] = -0.11
            self.Full_Structure.loc[indx,['rmin_half']] = 2.0
            indx = self.Full_Structure.index[(self.Full_Structure.aaid == aaid) &\
                                             (self.Full_Structure.atmtyp1 == 'OT1') &\
                                             (self.Full_Structure.chain == chain)][0]
            self.Full_Structure.loc[indx,['atmtyp1']] = 'O'
            self.Full_Structure.loc[indx,['atmtyp2']] = 'O'
            self.Full_Structure.loc[indx,['component']] = 'CARBO'
            self.Full_Structure.loc[indx,['charg']] = -0.51
            # NOTE: I removed (self.Full_Structure.ent_id != ent_id) |\ from the aboive resignment. 
            #       deletions must be for all models or the models will be incosistent.
        self.Full_Structure = self.Full_Structure.reset_index(drop=True)

    def build_pep_and_anchers(self,p1):
        """Generate ancher for attaching residues to a protein.
        An ancher is a set of atoms used for alignment. These atoms are
        not the ones being attached, and they are usually backbone atoms
        from the previous or next amino acid in the peptide PDB file.
        """
        pep_order = []
        first = {}
        last = {}
        ACE = {}
        CTER = {}
        h = [g.get_id()[1] for g in p1.get_residues()]
        pep_naa = len(h)
        for j in p1.get_residues():
            if j.get_id()[1] > 1 and j.get_id()[1] < pep_naa:
                atms = {}
                for k in j.get_atom():
                    atms[k.get_name()] = k
                self.pep[j.get_resname()] = atms
                pep_order.append(j.get_resname())
            if j.get_id()[1] == 1:
                for k in j.get_atom():
                    if k.get_name() in ['N','CA','C']:
                        first[k.get_name()] = k
                    if k.get_name() in ['CAY','HY1','HY2','HY3','CY','OY']:
                        ACE[k.get_name()] = k
            if j.get_id()[1] == pep_naa:
                for k in j.get_atom():
                    if k.get_name() in ['N','CA','C']:
                        last[k.get_name()] = k
                    if k.get_name() in ['C','OT1','OT2']:
                        CTER[k.get_name()] = k
        self.pep['ACE'] = ACE
        self.pep['CTER'] = CTER
        for i in range(len(pep_order)):
            if i == 0:
                self.NTERM_anchor[pep_order[i]] = first
                atms = {}
                for k in self.pep[pep_order[i+1]]:
                    if k in ['N','CA','C']:
                        atms[k] = self.pep[pep_order[i+1]][k]
                self.CTERM_anchor[pep_order[i]] = atms
                self.CTERM_anchor['ACE'] = first
            elif i == (len(pep_order)-1):
                self.CTERM_anchor[pep_order[i]] = last
                atms = {}
                for k in self.pep[pep_order[i-1]]:
                    if k in ['N','CA','C']:
                        atms[k] = self.pep[pep_order[i-1]][k]
                self.NTERM_anchor[pep_order[i]] = atms
                self.NTERM_anchor['CTER'] = last 
            else:
                atms = {}
                for k in self.pep[pep_order[i+1]]:
                    if k in ['N','CA','C']:
                        atms[k] = self.pep[pep_order[i+1]][k]
                self.CTERM_anchor[pep_order[i]] = atms
                atms = {}
                for k in self.pep[pep_order[i-1]]:
                    if k in ['N','CA','C']:
                        atms[k] = self.pep[pep_order[i-1]][k]
                self.NTERM_anchor[pep_order[i]] = atms
    
    def check_residue_parameters(self):
        """Checks that aminoacids from structures in the rcsb.org database are
        present int CHARMM parameter file. Use as a filter."""
        parameters_complete = False
        for i in self.strctr.get_models():
            for j in i.get_chains():
                for k in j.get_residues():
                    if k.get_resname() not in self.params.AA:
                        # HIS are excluded to be exchanged by HSD, HSP or HSE
                        if k.get_resname() != 'HIS':
                            print('Missing parameters for residue '+k.get_resname()+\
                            ' in model '+str(i.get_id())+' chain '+str(j.get_id()))
                        else:
                            parameters_complete = False
                    else:
                        parameters_complete = True
        return parameters_complete
                        
    def check_atom_parameters(self):
        """Checks that atoms from structures in the rcsb.org database are
        present int CHARMM parameter file. Use as a filter."""
        for i in self.strctr.get_models():
            for j in i.get_chains():
                for k in j.get_residues():
                    for l in k.get_atom():
                        if l.get_id() not in params.AA[k.get_resname()].atom_type:
                            print('    Missing parameters for atom '+l.get_id()+\
                            ' in model '+str(i.get_id())+' chain '+str(j.get_id())+\
                            ' residue '+k.get_resname()+' '+\
                            str(l.get_parent().get_full_id()[3][1]))

    def write_csv(self,basedir,filename):
        if basedir == '':
            basedir = '.'
        if basedir[-1] != "/":
            basedir += "/"
        self.Full_Structure.to_csv(basedir+filename+'.csv')
        
    def write_crd(self,file_path, model=0):
        print("Printing CRD")
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
        '''
        line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4}{:3}{:>12.5f}'
        atom_count = 0
        resi_count = 0
        lines = []
        for a in sp.get_models():
            if a.id == model:
                for b in a.get_chains():
                    for c in b.get_residues():
                        resi_count += 1
                        for d in c.get_atom():
                            atom_count += 1
                            lines.append(line.format(atom_count,resi_count, \
                            c.get_resname(),d.get_name(), d.get_coord()[0], \
                            d.get_coord()[1],d.get_coord()[2], b.id, \
                            c.get_full_id()[3][1],0))
        
        outFile = open(file_path, 'w')
        outFile.write('* CRD generated with fixed HS \n')
        outFile.write('* Name: NN\n')
        outFile.write('* Title: New Molecule Generated by pdb_preparation.py\n')
        outFile.write('*\n')
        outFile.write('{:>5}'.format(len(lines))+'\n')
        for i in lines:
            outFile.write(i+'\n')
        outFile.close()
        
    def report_sequence_structure_comparison(self):
        pass
    
    def read_dict_into_dataframes(self):
        """Puts structural information inside a DataFrame and checks for missing
        structural information by checking with sequence information available.
        The data below is stracted for programming clarity, but they could be
        obtain directly from the dictionaries."""
        # TODO: Some CIF file might have more that one chain, can be of different 
        # types, so this code may need to accomodate for that. When the coding id finished
        # I do not need to assign all of the variables below, I could get that directly
        # from the CIF parser header. I do this to remember what I got from the CIF
        # file during coding.
        '''Data items in the CELL category record details about the
           crystallographic cell parameters.'''
        col1 = '_cell.'
        col2 = ['entry_id','length_a','length_b','length_c','angle_alpha',\
                'angle_beta','angle_gamma','Z_PDB','pdbx_unique_axis']
        self.cell_df = pd.DataFrame(columns=col2)
        for i in col2:
            self.cell_df[i] = pd.Series(self.header[col1+i])
        '''Data items in the SYMMETRY category record details about the
           space-group symmetry.'''
        col1 = '_symmetry.'
        col2 = ['entry_id','space_group_name_H-M','pdbx_full_space_group_name_H-M',\
                'cell_setting','Int_Tables_number']
        self.symmetry_df = pd.DataFrame(columns=col2)
        for i in col2:
            self.symmetry_df[i] = pd.Series(self.header[col1+i])
        '''Data items in the ENTITY category record details (such as
           chemical composition, name and source) about the molecular
           entities that are present in the crystallographic structure.'''
        col1 = '_entity.'
        col2 = ['id','type','src_method','pdbx_description','formula_weight',\
                'pdbx_number_of_molecules','details']
        self.entity_df = pd.DataFrame(columns=col2)
        for i in col2:
            self.entity_df[i] = pd.Series(self.header[col1+i])
        '''Data items in the ENTITY_POLY category record details about the
           polymer, such as the type of the polymer, the number of
           monomers and whether it has nonstandard features.'''
        col1 = '_entity_poly.'
        col2 = ['entity_id','type','nstd_linkage','nstd_monomer',\
                'pdbx_seq_one_letter_code','pdbx_seq_one_letter_code_can',\
                'pdbx_strand_id']
        self.entity_poly_df = pd.DataFrame(columns=col2)
        for i in col2:
            self.entity_poly_df[i] = pd.Series(self.header[col1+i])
        '''Data items in the ENTITY_POLY_SEQ category specify the sequence
           of monomers in a polymer. Allowance is made for the possibility
           of microheterogeneity in a sample by allowing a given sequence
           number to be correlated with more than one monomer ID. The
           corresponding ATOM_SITE entries should reflect this
           heterogeneity.''' 
        col1 = '_entity_poly_seq.'
        col2 = ['entity_id','num','mon_id','hetero']
        self.entity_poly_seq_df = pd.DataFrame(columns=col2)
        for i in col2:
            self.entity_poly_seq_df[i] = pd.Series(self.header[col1+i])
        '''Data items in the STRUCT_REF category allow the author of a
           data block to relate the entities or biological units
           described in the data block to information archived in external
           databases.'''
        col1 = '_struct_ref.'
        # col2 = ['id','db_name','db_code','entity_id','pdbx_db_accession',\
        #         'pdbx_align_begin','pdbx_seq_one_letter_code','biol_id']
        # col2 = ['id','db_name','db_code','entity_id','pdbx_db_accession',\
        #         'pdbx_align_begin','biol_id']
        col2 = ['id','db_name','db_code','entity_id','pdbx_db_accession',\
                'pdbx_align_begin']
        self.struct_ref_df = pd.DataFrame(columns=col2)
        for i in col2:
            self.struct_ref_df[i] = pd.Series(self.header[col1+i])
        '''Data items in the STRUCT_REF_SEQ category provide a mechanism
           for indicating and annotating a region (or regions) of alignment
           between the sequence of an entity or biological unit described
           in the data block and the sequence in the referenced database
           entry.'''
        col1 = '_struct_ref_seq.'
        #    col2 = ['align_id','ref_id','pdbx_PDB_id_code','pdbx_strand_id',\
        #            'seq_align_beg','pdbx_seq_align_beg_ins_code','seq_align_end',\
        #            'pdbx_seq_align_end_ins_code','pdbx_db_accession','db_align_beg',\
        #            'pdbx_db_align_beg_ins_code','db_align_end',\
        #            'pdbx_db_align_end_ins_code','pdbx_auth_seq_align_beg',\
        #            'pdbx_auth_seq_align_end']
        col2 = ['align_id','ref_id','pdbx_PDB_id_code','pdbx_strand_id',\
                'seq_align_beg','pdbx_seq_align_beg_ins_code','seq_align_end',\
                'pdbx_seq_align_end_ins_code','pdbx_db_accession','db_align_beg',\
                'db_align_end','pdbx_auth_seq_align_beg','pdbx_auth_seq_align_end']
        self.struct_ref_seq_df = pd.DataFrame(columns=col2)
        for i in col2:
            self.struct_ref_seq_df[i] = pd.Series(self.header[col1+i])
        '''Data items in the STRUCT_ASYM category record details about the
           structural elements in the asymmetric unit.'''
        col1 = '_struct_asym.'
        col2 = ['id','pdbx_blank_PDB_chainid_flag','pdbx_modified','entity_id',\
                'details']
        self.struct_asym_df = pd.DataFrame(columns=col2)
        for i in col2:
            self.struct_asym_df[i] = pd.Series(self.header[col1+i])
        '''Data items in the ATOM_SITE category record details about
           the atom sites in a macromolecular crystal structure, such as
           the positional coordinates, atomic displacement parameters,
           magnetic moments and directions.
           The data items for describing anisotropic atomic
           displacement factors are only used if the corresponding items
           are not given in the ATOM_SITE_ANISOTROP category.'''
        col1 = "_atom_site."
        #col2 = ['group_PDB','id','type_symbol','label_atom_id','label_alt_id', \
        #        'label_comp_id','label_asym_id','label_entity_id','label_seq_id', \
        #        'pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z','occupancy', \
        #        'B_iso_or_equiv','Cartn_x_esd','Cartn_y_esd','Cartn_z_esd', \
        #        'occupancy_esd','B_iso_or_equiv_esd','pdbx_formal_charge', \
        #        'auth_seq_id','auth_comp_id','auth_asym_id','auth_atom_id', \
        #        'pdbx_PDB_model_num']
        # Changed du to inconsitencies in formatting
        col2 = ['group_PDB','id','type_symbol','label_atom_id','label_alt_id', \
                'label_comp_id','label_asym_id','label_entity_id','label_seq_id', \
                'pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z', 'occupancy',\
                'B_iso_or_equiv', \
                'pdbx_formal_charge', \
                'auth_seq_id','auth_comp_id','auth_asym_id','auth_atom_id', \
                'pdbx_PDB_model_num']
        self.atom_site_df = pd.DataFrame()
        for i in col2:
            self.atom_site_df[i] = self.header[col1+i]
        ##### So far the aboove fields will throw an exception if not found.
        ##### The following header fields will be checked for their presence
        ##### and if not found program is still allowed to run
        col1 = "_struct_conn."
        col2 = ['id','conn_type_id','ptnr1_label_asym_id','ptnr1_label_comp_id',\
                'ptnr1_label_seq_id','ptnr1_label_atom_id','ptnr1_symmetry',\
                'ptnr2_label_asym_id','ptnr2_label_comp_id','ptnr2_label_seq_id',\
                'ptnr2_label_atom_id','pdbx_dist_value']
        self.struct_conn_df = pd.DataFrame()     
        for i in col2:
            if (col1+i) in self.header:
                self.struct_conn_df[i] = self.header[col1+i]
            else:
                print('WARNING: Field '+(col1+i)+' is not found in CIF header')          
        # These model are obtained because it is a function that is used a lot
        # so it is obtained once to save repetitive iterations. It is also more
        # reliable to obtain this information from structural information directly.
        # Because the CIF format for other fields can vary 
        # TODO: Find a better place for this little one line function call.
        self.models = self.get_models()
    def get_models(self):
        mdls = list(set(self.atom_site_df.pdbx_PDB_model_num))
        mdls = [int(i) for i in mdls]
        mdls = sorted(mdls)
        mdls = [str(i) for i in mdls]
        return mdls
