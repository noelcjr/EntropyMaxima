# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 23:01:32 2017

@author: noel
"""
import os
import sys
folder = 'EntropyMaxima/'
path = '/home/noel/Projects/Protein_design/'+folder
sys_path = path+'/src'
if not os.path.isdir(sys_path):
    print("Error: The parameter directory "+sys_path+" is not found.")
    sys.exit(1)
param_path = path+'params/charmm27.ff/'
if not os.path.isdir(param_path):
    print("Error: The parameter directory "+param_path+" is not found.")
    sys.exit(1)
import CHARMM_Parser as CP
import Energy_Functions as EF
import Super_Structures as SS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from input_output import crd,psf

# This class manages the search for a linker or loop in search of some outcome. It could me moved later if it is too
# large to be here. It uses mmgbsa_ca_analysis.pyc and bash scripts to automate the search for a structure.
# Start by placinf in this class Component_Analysis_linkers.py
class mmgbsa_ca_searcher(object):
    def __init__(self):
        pass

class mmgbsa_ca_analysis(object):
    def __init__(self,oscwd):
        self.dirpath = oscwd
        self.components = {}
        self.df_pair = ''
        self.df_self = ''
    def rank_mmgbsa(self):
        folder = '/home/noel/Projects/Protein_design/EntropyMaxima/examples/CA_analysis/analyze_CA/'
        file_name = '2hiu_1rr_min_mmgbsa.txt'
        mmgbsa = pd.read_csv(folder+file_name)
        mmgbsa_v =mmgbsa.values
        mmgbsa_dict = {}
        for i in range(mmgbsa_v.shape[0]):
            if mmgbsa_v[i][0] in mmgbsa_dict:
                mmgbsa_dict[mmgbsa_v[i][0]] = mmgbsa_dict[mmgbsa_v[i][0]] + mmgbsa_v[i][2:]
            else:
                mmgbsa_dict[mmgbsa_v[i][0]] = mmgbsa_v[i][2:]
            if mmgbsa_v[i][1] in mmgbsa_dict:
                mmgbsa_dict[mmgbsa_v[i][1]] = mmgbsa_dict[mmgbsa_v[i][1]] + mmgbsa_v[i][2:]
            else:
                mmgbsa_dict[mmgbsa_v[i][1]] = mmgbsa_v[i][2:]
    # TODO https://en.wikipedia.org/wiki/Connected-component_labeling
    def connected_component_labeling(self):
        pass
    def mmgbsa_raw_summary(self,options):
        # TODO: NOT USED. For now only working with delta free energies. In the future I plan to use direct free
        #    energies of components in search of paterns or new insights. Direct free energies need 1-4 interactions to
        #    be more meaningful, and the 1-4 interactions code has a bug!
        if not options.range:
            options.range = '(-10,10,200):(-5,5,300)'
        filename_prefix_for_output = options.pair2.split('.')[0][0:-5]
        fullpath = self.dirpath+'/'+filename_prefix_for_output
        self.df_pair = pd.read_csv(self.dirpath+'/'+options.pair2)
        self.df_self = pd.read_csv(self.dirpath+'/'+options.self2)
        self.fix_cols_components()
        r1 = 'pos_ene(>+'+str(options.cut)+'kcal/mol)'
        r2 = 'neg_ene(<-'+str(options.cut)+'kcal/mol)'
        rows = ['min_val','max_val','mean','std','total_components','nan',r1,'pos_eng(%)',r2,'neg_eng(%)','included',\
                'included(%)','excluded','excluded(%)']
        self_energy_terms = 's'+self.df_self.columns[2:]
        pair_energy_terms = 'p'+self.df_pair.columns[2:]
        summary = pd.DataFrame(index=rows,columns=pair_energy_terms.union(self_energy_terms))
        for i in pair_energy_terms:
            summary[i] = summary[i].astype(float)
        for i in self_energy_terms:
            summary[i] = summary[i].astype(float)
        leng_self = len(self.df_self)
        for i in self_energy_terms:
            d = self.df_self.describe()[i[1:]]
            temp = []
            temp.append(d['min'])
            temp.append(d['max'])
            temp.append(d['mean'])
            temp.append(d['std'])
            temp.append(leng_self)
            np_self = self.df_self.loc[:,i[1:]].values
            temp.append(np.isnan(np_self).sum())
            neg_eng = np.where(np_self < -options.cut)
            pos_eng = np.where(np_self > options.cut)
            temp.append(len(pos_eng[0]))
            temp.append(len(pos_eng[0])/float(leng_self))
            temp.append(len(neg_eng[0]))
            temp.append(len(neg_eng[0])/float(leng_self))
            included = 100*(len(pos_eng[0])+len(neg_eng[0]))/float(leng_self)
            excluded = 100-included
            pos_neg_tot = len(pos_eng[0])+len(neg_eng[0])
            temp.append(pos_neg_tot)
            temp.append(included)
            temp.append(leng_self-pos_neg_tot)
            temp.append(excluded)
            summary[i] = np.array(temp)
        leng_pair = len(self.df_pair)
        for i in pair_energy_terms:
            d = self.df_pair.describe()[i[1:]]
            temp = []
            temp.append(d['min'])
            temp.append(d['max'])
            temp.append(d['mean'])
            temp.append(d['std'])
            temp.append(leng_pair)
            np_pairs = self.df_pair.loc[:,i[1:]].values
            temp.append(np.isnan(np_pairs).sum())
            neg_eng = np.where(np_pairs < -options.cut)
            pos_eng = np.where(np_pairs > options.cut)
            temp.append(len(pos_eng[0]))
            temp.append(len(pos_eng[0])/float(leng_pair))
            temp.append(len(neg_eng[0]))
            temp.append(len(neg_eng[0])/float(leng_pair))
            included = 100*(len(pos_eng[0])+len(neg_eng[0]))/float(leng_pair)
            excluded = 100-included
            pos_neg_tot = len(pos_eng[0])+len(neg_eng[0])
            temp.append(pos_neg_tot)
            temp.append(included)
            temp.append(leng_pair-pos_neg_tot)
            temp.append(excluded)
            summary[i] = np.array(temp)
        print(summary)
        summary.to_csv(fullpath+'_rawsum.csv')
        fig, axes = plt.subplots(nrows=2, ncols=2)
        ax0, ax1, ax2, ax3 = axes.flatten()
        pair_labels = ['pGBE', 'pGBE_Z', 'pEE','pEE_Z','pVWE','pVWE_Z']
        self_labels = ['sGBE', 'sGBE_Z', 'sSAE','sSAE_Z']
        x = self.df_pair.loc[:,self.df_pair.columns[2:]]
        y = self.df_self.loc[:,self.df_self.columns[2:]]
        tuples = []
        for i in options.range.split(':'):
            tuples.append(tuple(i[1:-1].split(',')))
        if len(tuples) != 2:
            print("ERROR: The string parameters for histograms must be two tuples separated by columns. For example:\
            \"(-10,10,200):(-5,5,300)\". First two values of the tuples can be real numbers corresponding to a minimum \
            and maimum range for a histogram. The las value is an integer corresponding to the number of bins.")
            sys.exit(1)
        # First histogram's ranges and numbers of bins are selcted automatically to contain the minimum and maximum
        # values of the distribution and 100 bins.
        n_bins = 100
        ax0.hist(x.values, n_bins, histtype='step', stacked=True, fill=False, label=pair_labels)
        ax0.hist(y.values, n_bins, histtype='step', stacked=True, fill=False, label=self_labels)
        ax0.legend(prop={'size': 10})
        ax0.set_title('Full Energy Range. Bins=100')
        # Second histogram uses the options.cut as limits on the x-axis and 10000 bins
        n_bins = 10000
        ax1.hist(x.values, n_bins, histtype='step', stacked=True, fill=False, label=pair_labels)
        ax1.hist(y.values, n_bins, histtype='step', stacked=True, fill=False, label=self_labels)
        ax1.legend(prop={'size': 10})
        ax1.set_title("Energy Range within options.cut. Bins=10000' ")
        ax1.set_xlim((-options.cut, options.cut))
        # Third and Fourth Histograms are defulted two (-10,10,200):(-5,5,300), but can be changed by user.
        n_bins = int(tuples[0][2])
        ax2.hist(x.values, n_bins, histtype='step', stacked=True, fill=False, label=pair_labels)
        ax2.hist(y.values, n_bins, histtype='step', stacked=True, fill=False, label=self_labels)
        ax2.legend(prop={'size': 10})
        ax2.set_title('Energy Range '+str(tuples[0][0:2])+". Bins="+str(tuples[0][2]))
        ax2.set_xlim((float(tuples[0][0]), float(tuples[0][1])))
        #axes.set_ylim((-10, 10))
        n_bins = int(tuples[1][2])
        ax3.hist(x.values, n_bins, histtype='step', stacked=True, fill=False, label=pair_labels)
        ax3.hist(y.values, n_bins, histtype='step', stacked=True, fill=False, label=self_labels)
        ax3.legend(prop={'size': 10})
        ax3.set_title('Energy Range '+str(tuples[1][0:2])+". Bins="+str(tuples[1][2]))
        ax3.set_xlim((float(tuples[1][0]), float(tuples[1][1])))
        #axes.set_ylim((-10, 10))
        fig.set_size_inches(12.0, 12.0)
        fig.savefig(self.dirpath+'/'+filename_prefix_for_output+'_rawhis.png')
    def mmgbsa_binding(self, options):
        filename_prefix_for_output = options.pair1.split('.')[0][0:-5]
        fullpath = self.dirpath+'/'+filename_prefix_for_output
        self.df_pair = pd.read_csv(self.dirpath+'/'+options.pair1)
        self.df_self = pd.read_csv(self.dirpath+'/'+options.self1)
        self.fix_cols_components()

        np_self = self.df_self.loc[:,self.df_self.columns].values
        np_pair = self.df_pair.loc[:,self.df_pair.columns].values
        # Binding between the protein and DNA and RNA
        MMGBSA = np.zeros((np_pair.shape[0],12))
        for i in range(np_pair.shape[0]):
            MMGBSA[i][0] = int(np_pair[i][0])                                           # Component 1
            MMGBSA[i][1] = int(np_pair[i][1])                                           # Component 2
            self1 = np_self[int(MMGBSA[i][0])][2] - np_self[int(MMGBSA[i][0])][3]       # Self 1 GB
            self2 = np_self[int(MMGBSA[i][1])][2] - np_self[int(MMGBSA[i][1])][3]       # Self 2 GB
            MMGBSA[i][2] = self1 + self2 + (np_pair[i][2] - np_pair[i][3])              # pair GB
            MMGBSA[i][3] = np_pair[i][4] - np_pair[i][5]                                # pair EE
            MMGBSA[i][4] = np_pair[i][6] - np_pair[i][7]                                # pair VDW folled by SA
            MMGBSA[i][5] = (np_self[int(MMGBSA[i][0])][4]-np_self[int(MMGBSA[i][0])][5])+(np_self[int(MMGBSA[i][1])][4]-np_self[int(MMGBSA[i][1])][5])
            MMGBSA[i][6] = MMGBSA[i][2] + MMGBSA[i][3]                                  # Polar
            MMGBSA[i][7] = MMGBSA[i][4] + MMGBSA[i][5]                                  # Non-Polar
            MMGBSA[i][8] = MMGBSA[i][6] + MMGBSA[i][7]                                  # MMGBSA
        # First rank unfiltered. Output results and decide later/
        mmgbsa = []
        for i in range(MMGBSA.shape[0]):
            mmgbsa.append([self.components[int(MMGBSA[i][0])],self.components[int(MMGBSA[i][1])],MMGBSA[i][2],\
                           MMGBSA[i][3],MMGBSA[i][4],MMGBSA[i][5],MMGBSA[i][6],MMGBSA[i][7],MMGBSA[i][8]])
        outFile = open(fullpath+'_mmgbsa.txt', 'w')
        outFile.write('comp1,comp2,dGB,dEE,dVDW,dSA,Polar,Nonpolar,MMGBSA\n')
        for i in mmgbsa:
            temp = ''
            temp += i[0]+','+i[1]+','
            for j in i[2:]:
                temp += "{:.3e}".format(j)+','
            outFile.write(temp[0:-1]+'\n')
        outFile.close()
    def mmgbsa_CA_bindingMatrix(self,options):
        ####################################      READ FILES     #######################################################
        #directory, filename = os.path.split()
        crd_file = crd(self.dirpath+'/'+options.crd)
        psf_file = psf(self.dirpath+'/'+options.psf)
        out_file = os.path.basename(options.crd).split('.')[0]
        ###################### After reading files, Generate and Index and Super Structure  ##############################
        params = CP.read_charmm_FF(param_path)
        insu = SS.Super_Structure(params, self.dirpath,'charmm_input')
        # At this point, a XPLOR psf could only have been creted from a complete structure, so no worries of gaps.
        insu.create_super_structure_df_from_CRD_PSF(crd_file,psf_file)
        ################################################################################################################
        insu.add_column_to_super_structure_df_from_CRD_PSF('gb_z0', self.dirpath+'/'+options.gb)
        insu.add_column_to_super_structure_df_from_CRD_PSF('gb_z500', self.dirpath+'/'+options.gbz)
        insu.add_column_to_super_structure_df_from_CRD_PSF('sa_z0', self.dirpath+'/'+options.sa)
        insu.add_column_to_super_structure_df_from_CRD_PSF('sa_z500', self.dirpath+'/'+options.saz)
        insu.create_column_with_ztranslated('z1_d', 'z1', "A",float(500))
        ###################### After reading files, Generate and Index a Super Structure  ##############################
        idx_ss = SS.index_super_structure(insu.Full_Structure)
        idx_ss.generate_indexes_from_Full_Structure()
        idx_ss.sort_comp_index_by_aaid_within_chains()
        new_index = []
        for i in range(len(idx_ss.comp_indx)):
            if idx_ss.comp_indx[i][3] == 'PRO':
                if idx_ss.comp_indx[i][4] == 'AMIN':
                    temp_list = range(idx_ss.comp_indx[i][5],idx_ss.comp_indx[i][6]+1)
                    temp_list = [temp_list[j] for j in [0,4,5,1,2,3]]
                    for j in temp_list[3:]:
                        insu.Full_Structure.loc[j,'component'] = 'SIDE2'
                    new_index = new_index + temp_list
                elif idx_ss.comp_indx[i][4] == 'SIDE':
                    new_index = new_index + range(idx_ss.comp_indx[i][5],idx_ss.comp_indx[i][6]+1)
                elif idx_ss.comp_indx[i][4] == 'CARB':
                    new_index = new_index + range(idx_ss.comp_indx[i][5],idx_ss.comp_indx[i][6]+1)
            else:
                new_index = new_index + range(idx_ss.comp_indx[i][5],idx_ss.comp_indx[i][6]+1)
        insu.Full_Structure = insu.Full_Structure.loc[new_index]
        insu.Full_Structure.component[(insu.Full_Structure.aa == 'PRO')&(insu.Full_Structure.atmtyp1 == 'CB')]= 'SIDE3'
        insu.Full_Structure.component[(insu.Full_Structure.aa == 'PRO')&(insu.Full_Structure.atmtyp1 == 'HB1')]= 'SIDE3'
        insu.Full_Structure.component[(insu.Full_Structure.aa == 'PRO')&(insu.Full_Structure.atmtyp1 == 'HB2')]= 'SIDE3'
        insu.Full_Structure.component[(insu.Full_Structure.aa == 'PRO')&(insu.Full_Structure.atmtyp1 == 'CG')]= 'SIDE4'
        insu.Full_Structure.component[(insu.Full_Structure.aa == 'PRO')&(insu.Full_Structure.atmtyp1 == 'HG1')]= 'SIDE4'
        insu.Full_Structure.component[(insu.Full_Structure.aa == 'PRO')&(insu.Full_Structure.atmtyp1 == 'HG2')]= 'SIDE4'
        insu.Full_Structure = insu.Full_Structure.reset_index(drop=True)
        idx_ss = SS.index_super_structure(insu.Full_Structure)
        idx_ss.generate_indexes_from_Full_Structure()
        idx_ss.sort_comp_index_by_aaid_within_chains()
        # CHARMM_Test/gbsw_ab_nbxmod5.dat 
        # ELECB,ELECU, GBENB, GBENU, VDWAB, VDWAU, ASPB, ASPU
        # -379.654,-332.641,-495.252, -595.916, 5672.62, 2338.65, 37.9186, 48.284
        ch_AB = insu.Full_Structure.loc[idx_ss.chn_indx[0][0]:idx_ss.chn_indx[-1][-1],['charg','epsilon',\
                      'rmin_half','x1','y1','z1','gb_z0','gb_z500','z1_d','sa_z0','sa_z500','atmtyp1','chain',\
                      'component','aa','aaid']].values
        EE_K = -327.90
        EE_K2 = -163.95
        #Surface tension coefficient        (sgamma) =   0.010 [kcal/mol/Angs**2]
        SA_K = 0.010
        l = len(idx_ss.comp_indx)
        idxFile = open(self.dirpath+'/'+out_file+'_MMGBSA.idx', 'w')
        if options.ver:
            outFile = open(self.dirpath+'/'+out_file+'_MMGBSA_verbose.out', 'w')
        else:
            outFile = open(self.dirpath+'/'+out_file+'_MMGBSA.out', 'w')
        tempGB = 0
        tempGB_Z = 0
        tempSA = 0
        tempSA_Z = 0
        tempEE = 0
        tempEE_Z = 0
        tempVDW = 0
        tempVDW_Z = 0
        pairGB = 0
        pairGB_Z = 0
        DGB = 0.0
        DGB_Z = 0.0
        DEE = 0.0
        DEE_Z = 0.0
        DVW = 0.0
        DVW_Z = 0.0
        DSA = 0.0
        DSA_Z = 0.0
        for i in range(0,l):
            idxFile.write(str(i)+" "+str(idx_ss.comp_indx[i][2])+" "+idx_ss.comp_indx[i][0]+" "+idx_ss.comp_indx[i][3]+\
                          " "+idx_ss.comp_indx[i][4]+"\n")
            for k in range(idx_ss.comp_indx[i][5],idx_ss.comp_indx[i][6]):
                temp5 = EE_K2*ch_AB[k][0]*ch_AB[k][0]/ch_AB[k][6]
                tempGB += temp5
                temp5_Z = EE_K2*ch_AB[k][0]*ch_AB[k][0]/ch_AB[k][7]
                tempGB_Z += temp5_Z
                tempSA += SA_K*ch_AB[k][9]
                tempSA_Z += SA_K*ch_AB[k][10]
                for kk in range(k+1,idx_ss.comp_indx[i][6]+1):
                    last = kk
                    d1 = np.linalg.norm(ch_AB[k][3:6]-ch_AB[kk][3:6])
                    r2 = d1*d1
                    d1_Z = np.linalg.norm(ch_AB[k][3:6]-ch_AB[kk][3:6])
                    r2_Z = d1_Z*d1_Z
                    temp3 = ch_AB[k][6]*ch_AB[kk][6]
                    temp4 = (-1*r2)/(4*temp3)
                    temp3_Z = ch_AB[k][7]*ch_AB[kk][7]
                    temp4_Z = (-1*r2_Z)/(4*temp3_Z)
                    denm = np.sqrt(r2 + (temp3*np.exp(temp4)))
                    denm_Z = np.sqrt(r2_Z + (temp3_Z*np.exp(temp4_Z)))
                    Kqq = EE_K*ch_AB[k][0]*ch_AB[kk][0]
                    temp5 = Kqq/denm
                    temp5_Z = Kqq/denm_Z
                    tempGB += temp5
                    tempGB_Z += temp5_Z
            tempSA += SA_K*ch_AB[last][9]
            tempSA_Z += SA_K*ch_AB[last][10]
            temp5 = EE_K2*ch_AB[last][0]*ch_AB[last][0]/ch_AB[last][6]
            tempGB += temp5
            temp5_Z = EE_K2*ch_AB[last][0]*ch_AB[last][0]/ch_AB[last][7]
            tempGB_Z += temp5_Z
            DGB += tempGB
            DGB_Z += tempGB_Z
            DSA += tempSA
            DSA_Z += tempSA_Z
            if options.ver:
                # Output for comparison with GB_Comp_Analy.prl
                outFile.write("SELF "+str(i+1)+" "+str(i+1)+" "+idx_ss.comp_indx[i][0]+" "+str(idx_ss.comp_indx[i][2])\
                              +" "+str(idx_ss.comp_indx[i][3])+" "+str(idx_ss.comp_indx[i][4])+"\n      GB "+\
                              str(tempGB)+"\n      GBZ "+str(tempGB_Z)+"\n      SA "+str(tempSA)+"\n      SAZ "+\
                              str(tempSA_Z)+"\n")
            else:
                outFile.write(str(i)+" "+str(i)+" "+str(tempGB-tempGB_Z)+" "+str(tempSA-tempSA_Z)+"\n")
            tempSA = 0.0
            tempSA_Z = 0.0
            for j in range(i+1,l):
                for m in range(idx_ss.comp_indx[i][5],idx_ss.comp_indx[i][6]+1):
                    for n in range(idx_ss.comp_indx[j][5],idx_ss.comp_indx[j][6]+1):
                        r = np.linalg.norm(ch_AB[m][3:6]-ch_AB[n][3:6])
                        r2 = r*r
                        r_Z = np.linalg.norm(np.array((ch_AB[m][3],ch_AB[m][4],ch_AB[m][8]))-\
                              np.array((ch_AB[n][3],ch_AB[n][4],ch_AB[n][8])))
                        r2_Z = r_Z*r_Z
                        temp = (332.06*ch_AB[m][0]*ch_AB[n][0])/r
                        temp_Z = (332.06*ch_AB[m][0]*ch_AB[n][0])/r_Z
                        tempEE += temp
                        tempEE_Z += temp_Z
                        Eps = np.sqrt(ch_AB[m][1]*ch_AB[n][1])
                        Rmin = ch_AB[m][2] + ch_AB[n][2]
                        A = Rmin/r
                        A_Z = Rmin/r_Z
                        A2 = A*A
                        A2_Z = A_Z*A_Z
                        A6 = A2*A2*A2
                        A6_Z = A2_Z*A2_Z*A2_Z
                        A12 = A6*A6
                        A12_Z = A6_Z*A6_Z
                        tempVDW += Eps*(A12-(2*A6))
                        tempVDW_Z += Eps*(A12_Z-(2*A6_Z))
                        temp3 = ch_AB[m][6]*ch_AB[n][6]
                        temp4 = (-1*r2)/(4*temp3)
                        temp3_Z = ch_AB[m][7]*ch_AB[n][7]
                        temp4_Z = (-1*r2_Z)/(4*temp3_Z)
                        denm = np.sqrt(r2 + (temp3*np.exp(temp4)))
                        denm_Z = np.sqrt(r2_Z + (temp3_Z*np.exp(temp4_Z)))
                        temp5 = (EE_K*ch_AB[m][0]*ch_AB[n][0])/denm
                        temp5_Z = (EE_K*ch_AB[m][0]*ch_AB[n][0])/denm_Z
                        pairGB += temp5
                        pairGB_Z += temp5_Z
                DGB += pairGB
                DGB_Z += pairGB_Z
                DEE += tempEE
                DEE_Z += tempEE_Z
                DVW += tempVDW
                DVW_Z += tempVDW_Z
                if options.ver:
                    # Output for comparison with GB_Comp_Analy.prl
                    outFile.write("PAIR "+str(i+1)+" "+str(j+1)+" ")
                    outFile.write(idx_ss.comp_indx[i][0]+" "+str(idx_ss.comp_indx[i][2])+" "+str(idx_ss.comp_indx[i][3])\
                                  +" "+str(idx_ss.comp_indx[i][4])+" - "+idx_ss.comp_indx[j][0]+" "+\
                                  str(idx_ss.comp_indx[j][2])+" "+str(idx_ss.comp_indx[j][3])+" "+\
                                  str(idx_ss.comp_indx[j][4])+"\n      EE "+str(tempEE)+"\n      EEZ "+str(tempEE_Z)+\
                                  "\n      VDW "+str(tempVDW)+"\n      VDWZ "+str(tempVDW_Z)+"\n      pGB "+str(pairGB)\
                                  +"\n      pGBZ "+str(pairGB_Z)+"\n")
                else:
                    outFile.write(str(i)+" "+str(j)+" "+str(tempEE-tempEE_Z)+" "+str(tempVDW-tempVDW_Z)+" "+\
                                  str(pairGB-pairGB_Z)+"\n")
                pairGB = 0
                pairGB_Z = 0
                tempEE = 0
                tempEE_Z = 0
                tempVDW = 0
                tempVDW_Z = 0
            tempGB = 0
            tempGB_Z = 0
            tempSA = 0
            tempSA_Z = 0

    def mmgbsa_CA_bindingMatrix_14exclusions(self,options):
        # TODO: Not working. Some problems with 1-4exclusions. Some missing exclusion? not enough testing.
        ef = EF.Energy_Functions(insu, idx_ss)
        ef.get_exclusions_1_4()
        count = 0
        exclu = []
        for i in ef.exclusions_1_4:
            for j in ef.exclusions_1_4[i]:
                if not len(ef.exclusions_1_4[i][j]) == 0:
                    for k in ef.exclusions_1_4[i][j]:
                        exclu.append(k)
                        count += 1
        ELECB = 0.0
        ELECU = 0.0
        count = 0
        for k in range(ch_AB.shape[0]-1):
            for kk in range(k+1,ch_AB.shape[0]):
                if not (k,kk) in exclu:
                    Eqq = 332.06*ch_AB[k][0]*ch_AB[kk][0]
                    d1 = np.linalg.norm(ch_AB[k][3:6]-ch_AB[kk][3:6])
                    ELECB += Eqq/d1
                    d1_z = np.linalg.norm(np.array((ch_AB[k][3],ch_AB[k][4],ch_AB[k][8]))-\
                                          np.array((ch_AB[kk][3],ch_AB[kk][4],ch_AB[kk][8])))
                    ELECU += Eqq/d1_z
                    count += 1
        print("Pair calculations=",count)
        print("ELECB=",ELECB)
        print("ELECU=",ELECU)

    def fix_cols_components(self):
        for i in self.df_self.columns:
            self.df_self.rename(columns={i:i.replace(' ','')}, inplace=True)
        for i in self.df_pair.columns:
            self.df_pair.rename(columns={i:i.replace(' ','')}, inplace=True)
        for i in self.df_self.index:
            self.components[self.df_self.loc[i,'i']] = self.df_self.loc[i,'Chain_entity_aa_aaid_comp']
