#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 22:35:08 2017

@author: noel
"""
import os
import sys
#import Super_Structures as SS
import pandas as pd
import numpy as np
#from CHARMM_Parser import *
# TODO: Make this class an extended class of Bio.PDB that solves all those
# problems modifying values.
class Energy_Functions(object):
    def __init__(self, superstructure, idx_superstructure):
        self.ss = superstructure
        self.idx_ss = idx_superstructure
        self.exclusions_1_4 = {}
        self.component_self_enrg = []
        self.component_pair_enrg = []
    def MMGBSA_CA2(self):
        ch_AB = self.ss.Full_Structure.loc[self.idx_ss.chn_indx[0][0]:self.idx_ss.chn_indx[-1][-1],['charg','epsilon',\
                      'rmin_half','x1','y1','z1','gb_z0','gb_z500','z1_d','sa_z0','sa_z500','atmtyp1','charg']].values
        l = len(self.idx_ss.comp_indx)
        self.component_self_enrg.append('i, Chain_entity_aa_aaid_comp, GBE, GBE_Z, SAE, SAE_Z')
        EE_K = -327.90
        EE_K2 = -163.95
        #Surface tension coefficient        (sgamma) =   0.010 [kcal/mol/Angs**2]
        #SA_K = 0.010
        SA_K = 1.0
        path='/home/noel/Projects/Protein_design/EntropyMaxima_Tests_Examples/examples/CA_analysis/SB_MMGBSA_CA/'
        outFile = open(path+'GBCA_python.out', 'w')
        for i in range(0,l):
            GBE = 0
            SAE = 0
            EEE = 0
            VWE = 0
            GBE_Z = 0
            SAE_Z = 0
            EEE_Z = 0
            VWE_Z = 0
            for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                GBE += EE_K2*ch_AB[k][0]*ch_AB[k][0]/ch_AB[k][6]
                GBE_Z += EE_K2*ch_AB[k][0]*ch_AB[k][0]/ch_AB[k][7]
                SAE += SA_K*ch_AB[k][9]
                SAE_Z += SA_K*ch_AB[k][10]
                #for kk in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                for kk in range(k,self.idx_ss.comp_indx[i][6]+1):
                    #if not (k,kk) in self.exclusions_1_4[i][i]:
                    d1 = np.linalg.norm(ch_AB[k][3:6]-ch_AB[kk][3:6])
                    Kqq = EE_K*ch_AB[k][0]*ch_AB[kk][0]
                    rr = np.sqrt(ch_AB[k][6]*ch_AB[kk][6])
                    dd1 = d1*d1
                    GBE += Kqq/np.sqrt((dd1)+(rr*np.exp(-1*dd1/4*rr)))
                    rr_z = np.sqrt(ch_AB[k][7]*ch_AB[kk][7])
                    GBE_Z += Kqq/np.sqrt((dd1)+(rr_z*np.exp(-1*dd1/4*rr_z)))
            temp = "SELF "+str(i+1)+' '+self.idx_ss.comp_indx[i][0]+' '+str(self.idx_ss.comp_indx[i][2])+' '
            temp += self.idx_ss.comp_indx[i][3]+' '+self.idx_ss.comp_indx[i][4]+' '+str(GBE)+' '+str(GBE_Z)+' '
            temp += str(SAE)+' '+str(SAE_Z)
            outFile.write(temp+"\n")
            # FORMAT "{:.3e}".format(GBE)
            #self.component_self_enrg.append(temp)
            #self.component_pair_enrg.append('i, j, GBE, GBE_Z, EEE, EEE_Z, VWE, VWE_Z')
            #for i in range(0,l-1):
            for j in range(i+1,l):
                GBE = 0.0
                SAE = 0.0
                EEE = 0.0
                VWE = 0.0
                GBE_Z = 0.0
                SAE_Z = 0.0
                EEE_Z = 0.0
                VWE_Z = 0.0
                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                        #if not (k,kk) in self.exclusions_1_4[i][j]:
                        d1 = np.linalg.norm(ch_AB[k][3:6]-ch_AB[kk][3:6])
                        qq = ch_AB[k][0]*ch_AB[kk][0]
                        Eqq = 332.06*qq
                        EEE += Eqq/d1
                        Rmin = ch_AB[k][2] + ch_AB[kk][2]
                        EPS = np.sqrt(ch_AB[k][1]*ch_AB[kk][1])
                        A = Rmin/d1
                        A2 = A*A
                        A6 = A2*A2*A2
                        A12 = A6*A6
                        VWE += EPS*(A12-(2*A6))
                        Gqq = EE_K*qq
                        rr = np.sqrt(ch_AB[k][6]*ch_AB[kk][6])
                        dd1 = d1*d1
                        GBE += Gqq/np.sqrt((dd1)+(rr*np.exp(-1*dd1/4*rr)))
                        d1_z = np.linalg.norm(np.array((ch_AB[k][3],ch_AB[k][4],ch_AB[k][8]))-\
                                              np.array((ch_AB[kk][3],ch_AB[kk][4],ch_AB[kk][8])))
                        EEE_Z += Eqq/d1_z
                        A = Rmin/d1_z
                        A2 = A*A
                        A6 = A2*A2*A2
                        A12 = A6*A6
                        VWE_Z += EPS*(A12-(2*A6))
                        rr_z = np.sqrt(ch_AB[k][7]*ch_AB[kk][7])
                        dd1_z = d1_z*d1_z
                        GBE_Z += Gqq/np.sqrt((dd1_z)+(rr_z*np.exp(-1*dd1_z/4*rr_z)))
                temp = str(i)+','+str(j)
                temp += ','+"{:.3e}".format(GBE)+','+"{:.3e}".format(GBE_Z)
                temp += ','+"{:.3e}".format(EEE)+','+"{:.3e}".format(EEE_Z)
                temp += ','+"{:.3e}".format(VWE)+','+"{:.3e}".format(VWE_Z)
                self.component_pair_enrg.append(temp)
    def MMGBSA_CA(self):
        ch_AB = self.ss.Full_Structure.loc[self.idx_ss.chn_indx[0][0]:self.idx_ss.chn_indx[-1][-1],['charg','epsilon',\
                      'rmin_half','x1','y1','z1','gb_z0','gb_z500','z1_d','sa_z0','sa_z500','atmtyp1','charg']].values
        l = len(self.idx_ss.comp_indx)
        self.component_self_enrg.append('i, Chain_entity_aa_aaid_comp, GBE, GBE_Z, SAE, SAE_Z')
        EE_K = -327.90
        EE_K2 = -163.95
        #Surface tension coefficient        (sgamma) =   0.010 [kcal/mol/Angs**2]
        SA_K = 0.010
        for i in range(0,l):
            GBE = 0
            SAE = 0
            EEE = 0
            VWE = 0
            GBE_Z = 0
            SAE_Z = 0
            EEE_Z = 0
            VWE_Z = 0
            for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                GBE += EE_K2*ch_AB[k][0]*ch_AB[k][0]/ch_AB[k][6]
                GBE_Z += EE_K2*ch_AB[k][0]*ch_AB[k][0]/ch_AB[k][7]
                SAE += SA_K*ch_AB[k][9]
                SAE_Z += SA_K*ch_AB[k][10]
                #for kk in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                for kk in range(k,self.idx_ss.comp_indx[i][6]+1):
                    if not (k,kk) in self.exclusions_1_4[i][i]:
                        d1 = np.linalg.norm(ch_AB[k][3:6]-ch_AB[kk][3:6])
                        Kqq = EE_K*ch_AB[k][0]*ch_AB[kk][0]
                        rr = np.sqrt(ch_AB[k][6]*ch_AB[kk][6])
                        dd1 = d1*d1
                        GBE += Kqq/np.sqrt((dd1)+(rr*np.exp(-1*dd1/4*rr)))
                        rr_z = np.sqrt(ch_AB[k][7]*ch_AB[kk][7])
                        GBE_Z += Kqq/np.sqrt((dd1)+(rr_z*np.exp(-1*dd1/4*rr_z)))
            temp = str(i)+','+self.idx_ss.comp_indx[i][0]+'_'+str(self.idx_ss.comp_indx[i][1])+'_'+str(self.idx_ss.comp_indx[i][2])
            temp += '_'+self.idx_ss.comp_indx[i][3]+'_'+self.idx_ss.comp_indx[i][4]+','+"{:.3e}".format(GBE)+','+"{:.3e}".format(GBE_Z)
            temp += ','+"{:.3e}".format(SAE)+','+"{:.3e}".format(SAE_Z)
            self.component_self_enrg.append(temp)
        self.component_pair_enrg.append('i, j, GBE, GBE_Z, EEE, EEE_Z, VWE, VWE_Z')
        for i in range(0,l-1):
            for j in range(i+1,l):
                GBE = 0.0
                SAE = 0.0
                EEE = 0.0
                VWE = 0.0
                GBE_Z = 0.0
                SAE_Z = 0.0
                EEE_Z = 0.0
                VWE_Z = 0.0
                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                        if not (k,kk) in self.exclusions_1_4[i][j]:
                            d1 = np.linalg.norm(ch_AB[k][3:6]-ch_AB[kk][3:6])
                            qq = ch_AB[k][0]*ch_AB[kk][0]
                            Eqq = 332.06*qq
                            EEE += Eqq/d1
                            Rmin = ch_AB[k][2] + ch_AB[kk][2]
                            EPS = np.sqrt(ch_AB[k][1]*ch_AB[kk][1])
                            A = Rmin/d1
                            A2 = A*A
                            A6 = A2*A2*A2
                            A12 = A6*A6
                            VWE += EPS*(A12-(2*A6))
                            Gqq = EE_K*qq
                            rr = np.sqrt(ch_AB[k][6]*ch_AB[kk][6])
                            dd1 = d1*d1
                            GBE += Gqq/np.sqrt((dd1)+(rr*np.exp(-1*dd1/4*rr)))
                            d1_z = np.linalg.norm(np.array((ch_AB[k][3],ch_AB[k][4],ch_AB[k][8]))-\
                                                  np.array((ch_AB[kk][3],ch_AB[kk][4],ch_AB[kk][8])))
                            EEE_Z += Eqq/d1_z
                            A = Rmin/d1_z
                            A2 = A*A
                            A6 = A2*A2*A2
                            A12 = A6*A6
                            VWE_Z += EPS*(A12-(2*A6))
                            rr_z = np.sqrt(ch_AB[k][7]*ch_AB[kk][7])
                            dd1_z = d1_z*d1_z
                            GBE_Z += Gqq/np.sqrt((dd1_z)+(rr_z*np.exp(-1*dd1_z/4*rr_z)))
                temp = str(i)+','+str(j)
                temp += ','+"{:.3e}".format(GBE)+','+"{:.3e}".format(GBE_Z)
                temp += ','+"{:.3e}".format(EEE)+','+"{:.3e}".format(EEE_Z)
                temp += ','+"{:.3e}".format(VWE)+','+"{:.3e}".format(VWE_Z)
                self.component_pair_enrg.append(temp)
    def get_exclusions_1_4(self):
        for i in range(0,len(self.idx_ss.aaid_indx2)):
            self.exclusions_1_4[i] = {}
            for j in range(i,len(self.idx_ss.aaid_indx2)):
                self.exclusions_1_4[i][j] = []
                if self.idx_ss.comp_indx[i][0] == self.idx_ss.comp_indx[j][0]:
                    if self.idx_ss.comp_indx[i][2] == self.idx_ss.comp_indx[j][2]:
                        if (self.idx_ss.comp_indx[i][4] == 'ACE1'):
                            if (self.idx_ss.comp_indx[j][4] == 'ACE1'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        self.exclusions_1_4[i][j].append((k,kk))
                            elif (self.idx_ss.comp_indx[j][4] == 'ACE2'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if (self.idx_ss.comp_indx[i][5] == k) or (self.idx_ss.comp_indx[j][5] == kk):
                                            self.exclusions_1_4[i][j].append((k,kk))
                            elif (self.idx_ss.comp_indx[j][4] == 'AMIN'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if (self.idx_ss.comp_indx[i][5] == k) and (self.idx_ss.comp_indx[j][5] == kk):
                                            self.exclusions_1_4[i][j].append((k,kk))
                            elif self.idx_ss.comp_indx[j][4] in ['SIDE','CARB','CTER']:
                                pass
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'ACE2'):
                            if (self.idx_ss.comp_indx[j][4] == 'ACE2'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        self.exclusions_1_4[i][j].append((k,kk))
                            elif (self.idx_ss.comp_indx[j][4] == 'AMIN'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    if (self.idx_ss.comp_indx[i][5] == k):
                                        for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                            if kk > (self.idx_ss.comp_indx[j][5]+1):
                                                if self.ss.Full_Structure.loc[kk,['atmtyp1']][0] != 'CA':
                                                    pass
                                                else:
                                                    self.exclusions_1_4[i][j].append((k,kk))
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                    else:
                                        for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                            if kk != self.idx_ss.comp_indx[j][5]:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                            elif self.idx_ss.comp_indx[j][4] in ['SIDE','CARB','CTER']:
                                pass
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'AMIN'):
                            if (self.idx_ss.comp_indx[j][4] == 'AMIN'):
                                base_indx = self.idx_ss.comp_indx[i][5]
                                if self.idx_ss.comp_indx[j][3] == 'PRO':
                                    l = {0:'N', 1:'CD', 2:'HD1', 3:'HD2', 4:'CA', 5:'HA'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(3,3),\
                                         (4,4),(4,5),(5,5)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'GLY':
                                    l = {0:'N', 1:'HN', 2:'CA', 3:'HA1', 4:'HA2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(1,1),(1,2),(2,2),(2,3),(3,3),(2,4),(3,3),(3,4),\
                                         (4,4)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                else:
                                    l = {0:'N', 1:'HN', 2:'CA', 3:'HA'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(1,1),(1,2),(2,2),(2,3),(3,3)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                            elif (self.idx_ss.comp_indx[j][4] == 'SIDE'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'N':
                                            if self.idx_ss.comp_indx[j][3] == 'PRO':
                                                if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['CB','CG']:
                                                    pass
                                                else:
                                                    self.exclusions_1_4[i][j].append((k,kk))
                                            else:
                                                if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['CB']:
                                                    pass
                                                else:
                                                    self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'HN':
                                            pass
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'HA':
                                            if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['CB']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'CA':
                                            if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['CG','CB','HB1','HB2','HB3','HB','SG','CG1','CG2','OG','OG1']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'CD':
                                            if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['CB','CG','HG1','HG2']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'HD1':
                                            if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['CD','HD2','CG']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'HD2':
                                            if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['CD','HD1','CG']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        else:
                                            self.exclusions_1_4[i][j].append((k,kk))
                            elif (self.idx_ss.comp_indx[j][4] == 'CARB'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'N':
                                            if self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['O']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] in ['HA','HA1','HA2']:
                                            if self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['O']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'CA':
                                            self.exclusions_1_4[i][j].append((k,kk))
                                        else:
                                            pass
                            elif (self.idx_ss.comp_indx[j][4] == 'CTER'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'N':
                                            if self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['OT1','OT2']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] in ['HA','HA1','HA2']:
                                            if self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['OT1','OT2']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'CA':
                                            self.exclusions_1_4[i][j].append((k,kk))
                                        else:
                                            pass
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'SIDE'):
                            if (self.idx_ss.comp_indx[j][4] == 'SIDE'):
                                base_indx = self.idx_ss.comp_indx[i][5]
                                if self.idx_ss.comp_indx[j][3] == 'ALA':
                                    l = {0:'CB',1:'HB1',2:'HB2',3:'HB3'}
                                    e =[(0,0),(0,1),(0,2),(0,3),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'ARG':
                                    l = {0:'CB',1:'HB1',2:'HB2',3:'CG',4:'HG1',5:'HG2',6:'CD',7:'HD1',8:'HD2',9:'NE',\
                                         10:'HE',11:'CZ',12:'NH1',13:'NH11',14:'NH12',15:'NH2',16:'HH21',17:'HH22'}
                                    e =[(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),\
                                        (3,4),(3,5),(3,6),(3,7),(3,8),(3,9),(4,4),(4,5),(4,6),(5,5),(5,6),(6,6),(6,7),\
                                        (6,8),(6,9),(6,10),(6,11),(7,7),(7,8),(7,9),(8,8),(8,9),(9,9),(9,10),(9,11),\
                                        (9,12),(9,15),(10,10),(10,11),(11,11),(11,12),(11,13),(11,14),(11,15),(11,16),\
                                        (11,17),(12,12),(12,13),(12,14),(12,15),(13,13),(13,14),(14,14),(15,15),(15,16),\
                                        (15,17),(16,16),(16,17),(17,17)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'ASN':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'OD1', 5:'ND2', 6:'HD21', 7:'HD22'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),(3,4),\
                                         (3,5),(3,6),(3,7),(4,4),(4,5),(5,5),(5,6),(5,7),(6,6),(6,7),(7,7)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'ASP':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'OD1', 5:'OD2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),(3,4),\
                                         (3,5),(4,4),(4,5),(5,5)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'CYS':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'SG', 4:'HG1'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),(3,4)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                    #exclusions CB to itself, CB to HB1, CB to HB2, CB to SG, CB to HG1, HB1 to itself
                                    # HB1 to HB2, HB1 to SG, HB2 to itself, HB2 to SG, SG to itself, SG to HG1
                                elif self.idx_ss.comp_indx[j][3] == 'GLN':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'HG1', 5:'HG2', 6:'CD', 7:'OE1', 8:'NE2',\
                                         9:'HE21', 10:'HE22'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),\
                                         (3,4),(3,5),(3,6),(3,7),(3,8),(4,4),(4,5),(4,6),(5,5),(5,6),(6,6),(6,7),(6,8),\
                                         (6,9),(6,10),(7,7),(7,8),(7,9),(7,10),(8,8),(8,9),(8,10),(9,9),(9,10),(10,10)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'GLU':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'HG1', 5:'HG2', 6:'CD', 7:'OE1', 8:'OE2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),\
                                         (3,4),(3,5),(3,6),(3,7),(3,8),(4,4),(4,5),(4,6),(5,5),(5,6),(6,6),(6,7),(6,8),\
                                         (7,7),(7,8),(8,8)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'HSD':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'ND1', 4:'HD1', 5:'CG', 6:'CE1', 7:'HE1', 8:'NE2',\
                                         9:'CD2',10:'HD2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,5),(0,9),(1,1),(1,2),(1,5),(2,2),(2,5),(3,3),\
                                         (3,4),(3,5),(3,6),(3,7),(3,8),(3,9),(4,4),(4,5),(4,6),(5,5),(5,6),(5,8),(5,9),
                                         (5,10),(6,6),(6,7),(6,8),(6,9),(7,7),(7,8),(8,8),(8,9),(8,10),(9,9),(9,10),(10,10)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'HSE':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'ND1', 4:'CG', 5:'CE1', 6:'HE1', 7:'NE2', 8:'HE2',\
                                         9:'CD2',10:'HD2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,9),(1,1),(1,2),(1,4),(2,2),(2,4),(3,3),(3,4),\
                                         (3,5),(3,6),(3,7),(3,9),(4,4),(4,5),(4,7),(4,9),(4,10),(5,5),(5,6),(5,7),(5,8),
                                         (5,9),(6,6),(6,7),(7,7),(7,8),(7,9),(7,10),(8,8),(8,9),(9,9),(9,10),(10,10)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'HSP':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CD2', 4:'HD2', 5:'CG', 6:'NE2', 7:'HE2', 8:'ND1',\
                                         9:'HD1',10:'CE1',11:'HE1'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,5),(0,8),(1,1),(1,2),(1,5),(2,2),(2,5),(3,3),(3,4),\
                                         (3,5),(3,6),(3,7),(3,8),(3,10),(4,4),(4,5),(4,6),(5,5),(5,6),(5,8),(5,9),(5,10),\
                                         (6,6),(6,7),(6,8),(6,10),(6,11),(7,7),(7,10),(8,8),(8,9),(8,10),(8,11),(9,9),\
                                         (9,10),(10,10),(10,11),(11,11)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'ILE':
                                    l = {0:'CB', 1:'HB', 2:'CG2', 3:'HG21', 4:'HG22', 5:'HG23', 6:'CG1', 7:'HG11',\
                                         8:'HG12',9:'CD',11:'HD2',12:'HD3'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(0,8),(0,9),(1,1),(1,2),(1,6),\
                                         (2,2),(2,3),(3,3),(2,4),(2,5),(2,6),(3,3),(3,4),(3,5),(4,4),(4,5),(5,5),(6,6),\
                                         (6,7),(6,8),(6,9),(6,10),(6,11),(6,12),(7,7),(7,8),(7,9),(8,8),(8,9),(9,9),\
                                         (9,10),(9,11),(9,12),(10,10),(10,11),(10,12),(11,11),(11,12),(12,12)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'LEU':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'HG', 5:'CD1', 6:'HD11', 7:'HD12',\
                                         8:'HD13',9:'CD2',11:'HD21',12:'HD22',13:'HD23'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,9),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),\
                                         (3,4),(3,5),(3,6),(3,7),(3,8),(3,9),(3,10),(3,11),(3,12),(4,4),(4,5),(4,9),\
                                         (5,5),(5,6),(5,7),(5,8),(5,9),(6,6),(6,7),(6,8),(7,7),(7,8),(8,8),(9,9),(9,10),\
                                         (9,11),(9,12),(10,10),(10,11),(10,12),(11,11),(11,12),(12,12)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'LYS':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'HG1', 5:'HG2', 6:'CD', 7:'HD1',\
                                         8:'HD2',9:'CE',10:'HE1',11:'HE2',12:'HZ',13:'HZ1',14:'HZ2',15:'HZ3'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),\
                                         (3,4),(3,5),(3,6),(3,7),(3,8),(3,9),(4,4),(4,5),(4,6),(5,5),(5,6),(6,6),(6,7),\
                                         (6,8),(6,9),(6,10),(6,11),(6,12),(7,7),(7,8),(7,9),(8,8),(8,9),(9,9),(9,10),\
                                         (9,11),(9,12),(9,13),(9,14),(9,15),(10,10),(10,11),(10,12),(11,11),(11,12),\
                                         (12,13),(12,14),(12,15),(12,12),(13,13),(13,14),(13,15),(14,14),(14,15),(15,15)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'MET':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'HG1', 5:'HG2', 6:'SD', 7:'CE', 8:'HE1',\
                                         9:'HE2',10:'HE3'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),\
                                         (3,4),(3,5),(3,6),(3,7),(4,4),(4,5),(4,6),(5,5),(5,6),(6,6),(6,7),(6,8),(6,9),\
                                         (6,10),(7,7),(7,8),(7,9),(7,10),(8,8),(8,9),(8,10),(9,9),(9,10),(10,10)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'PHE':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'CD1', 5:'HD1', 6:'CE1', 7:'HE1',\
                                         8:'CZ',9:'HZ',10:'CD2',11:'HD2',12:'CE2',13:'HE2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,10),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),\
                                         (3,4),(3,5),(3,6),(3,10),(3,11),(3,12),(4,4),(4,5),(4,6),(4,7),(4,8),(4,10),\
                                         (5,5),(5,6),(6,6),(6,7),(6,8),(6,9),(6,12),(7,7),(7,8),(8,8),(8,9),(8,10),\
                                         (8,12),(8,13),(9,9),(9,12),(10,10),(10,11),(10,12),(10,13),(11,11),(11,12),\
                                         (12,12),(12,13),(13,13)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'PRO':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'HG1', 5:'HG2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),(3,4),\
                                         (3,5),(4,4),(4,5),(5,5)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'SER':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'OG', 4:'HG1'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),(3,4),(4,4)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'THR':
                                    l = {0:'CB', 1:'HB', 2:'OG1', 3:'HG1', 4:'CG2', 5:'HG21', 6:'HG22', 7:'HG23'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(1,1),(1,2),(1,4),(2,2),(2,3),\
                                         (2,4),(3,3),(4,4),(4,5),(4,6),(4,7),(5,5),(5,6),(5,7),(6,6),(6,7),(7,7)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'TRP':
                                    l = {0:'CB',1:'HB1',2:'HB2',3:'CG',4:'CD1',5:'HD1',6:'NE1',7:'HE1',8:'CE2',9:'CD2',\
                                         10:'CE3',11:'HE3',12:'CZ3',13:'HZ3',14:'CZ2',15:'HZ2',16:'CH2',17:'HH2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,9),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),(3,4),\
                                         (3,5),(3,6),(3,8),(3,9),(3,10),(4,4),(4,5),(4,6),(4,7),(4,8),(4,9),(5,5),(5,6),\
                                         (6,6),(6,7),(6,8),(6,9),(6,14),(7,7),(7,8),(8,8),(8,9),(8,10),(8,14),(8,15),\
                                         (8,16),(9,9),(9,10),(9,11),(9,12),(9,14),(10,10),(10,11),(10,12),(10,13),(10,16),\
                                         (11,11),(11,12),(12,12),(12,13),(12,14),(12,16),(12,17),(13,13),(13,16),(14,14),
                                         (14,15),(14,16),(14,17),(15,15),(15,16),(16,16),(16,17),(17,17)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'TYR':
                                    l = {0:'CB', 1:'HB1', 2:'HB2', 3:'CG', 4:'CD1', 5:'HD1', 6:'CE1', 7:'HE1',\
                                         8:'CZ',9:'OH',10:'HH',12:'HD2',13:'CE2',14:'HE2'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,11),(1,1),(1,2),(1,3),(2,2),(2,3),(3,3),\
                                         (3,4),(3,5),(3,6),(3,11),(3,12),(3,13),(4,4),(4,5),(4,6),(4,7),(4,8),(4,11),\
                                         (5,5),(5,6),(6,6),(6,7),(6,8),(6,9),(6,13),(7,7),(7,8),(8,8),(8,9),(8,10),\
                                         (8,11),(8,13),(8,14),(9,9),(9,10),(9,13),(10,10),(11,11),(11,12),(11,13),\
                                         (11,14),(12,12),(12,13),(13,13),(13,14),(14,14)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'VAL':
                                    l = {0:'CB', 1:'HB', 2:'CG1', 3:'HG11', 4:'HG12', 5:'HG13', 6:'CG2', 7:'HG21',\
                                         8:'HG22', 9:'HG23'}
                                    e = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(0,8),(0,9),(1,1),(1,2),(1,6),\
                                         (2,2),(2,3),(2,4),(2,5),(2,6),(3,3),(3,4),(3,5),(4,4),(4,5),(5,5),(6,6),(6,7),\
                                         (6,8),(6,9),(7,7),(7,8),(7,9),(8,8),(8,9),(9,9)]
                                    for nn in e:
                                        self.exclusions_1_4[i][j].append((base_indx+nn[0],base_indx+nn[1]))
                                elif self.idx_ss.comp_indx[j][3] == 'GUA':
                                    pass
                                elif self.idx_ss.comp_indx[j][3] == 'ADE':
                                    pass
                                elif self.idx_ss.comp_indx[j][3] == 'CYT':
                                    pass
                                elif self.idx_ss.comp_indx[j][3] == 'THY':
                                    pass
                                elif self.idx_ss.comp_indx[j][3] == 'URA':
                                    pass
                                else:
                                    print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][3]+". Exit now.")
                                    sys.exit(1)
                            elif (self.idx_ss.comp_indx[j][4] == 'CARB'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'CB':
                                            if self.ss.Full_Structure.loc[kk,['atmtyp1']][0] != 'C':
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        else:
                                            pass
                            elif (self.idx_ss.comp_indx[j][4] == 'CTER'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'CB':
                                            if self.ss.Full_Structure.loc[kk,['atmtyp1']][0] != 'C':
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        else:
                                            pass
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'CARB'):
                            if (self.idx_ss.comp_indx[j][4] == 'CARB'):
                                base_indx = self.idx_ss.comp_indx[i][5]
                                self.exclusions_1_4[i][j].append((base_indx,base_indx))      # C to itself
                                self.exclusions_1_4[i][j].append((base_indx,base_indx+1))    # C to O
                                self.exclusions_1_4[i][j].append((base_indx+1,base_indx+1))  # O to itself
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'CTER'):
                            if (self.idx_ss.comp_indx[j][4] == 'CTER'):
                                base_indx = self.idx_ss.comp_indx[i][5]
                                self.exclusions_1_4[i][j].append((base_indx,base_indx))      # C to itself
                                self.exclusions_1_4[i][j].append((base_indx,base_indx+1))    # C to OT1
                                self.exclusions_1_4[i][j].append((base_indx,base_indx+2))    # C to OT2
                                self.exclusions_1_4[i][j].append((base_indx+1,base_indx+1))  # OT1 to itself
                                self.exclusions_1_4[i][j].append((base_indx+1,base_indx+2))  # OT1 to OT2
                                self.exclusions_1_4[i][j].append((base_indx+2,base_indx+2))  # OT2 to itself
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                    elif self.idx_ss.comp_indx[i][2] == (self.idx_ss.comp_indx[j][2]-1):
                        if (self.idx_ss.comp_indx[i][4] == 'ACE1'):
                            if self.idx_ss.comp_indx[j][4] in ['ACE1','ACE2','AMIN','SIDE','CARB','CTER']:
                                pass
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'ACE2'):
                            if self.idx_ss.comp_indx[j][4] in ['ACE1','ACE2','AMIN','SIDE','CARB','CTER']:
                                pass
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'AMIN'):
                            if self.idx_ss.comp_indx[j][4] in ['ACE1','ACE2','SIDE','CARB','CTER']:
                                pass
                            elif (self.idx_ss.comp_indx[j][4] == 'AMIN'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'CA':
                                            if self.ss.Full_Structure.loc[kk,['atmtyp1']][0] == 'N':
                                                self.exclusions_1_4[i][j].append((k,kk))
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'SIDE'):
                            if self.idx_ss.comp_indx[j][4] in ['ACE1','ACE2','AMIN','SIDE','CARB','CTER']:
                                pass
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'CARB'):
                            if self.idx_ss.comp_indx[j][4] in ['ACE1','ACE2','SIDE','CARB','CTER']:
                                pass
                            elif (self.idx_ss.comp_indx[j][4] == 'AMIN'):
                                for k in range(self.idx_ss.comp_indx[i][5],self.idx_ss.comp_indx[i][6]+1):
                                    for kk in range(self.idx_ss.comp_indx[j][5],self.idx_ss.comp_indx[j][6]+1):
                                        if self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'C':
                                            if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['N','HN','CA']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        elif self.ss.Full_Structure.loc[k,['atmtyp1']][0] == 'O':
                                            if not self.ss.Full_Structure.loc[kk,['atmtyp1']][0] in ['N']:
                                                pass
                                            else:
                                                self.exclusions_1_4[i][j].append((k,kk))
                                        else:
                                            self.exclusions_1_4[i][j].append((k,kk))
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                        elif (self.idx_ss.comp_indx[i][4] == 'CTER'):
                            if self.idx_ss.comp_indx[j][4] in ['ACE1','ACE2','AMIN','SIDE','CARB','CTER']:
                                pass
                            else:
                                print("Error: Unknown or Impossible component type "+self.idx_ss.comp_indx[j][4]+". Exit now.")
                                sys.exit(1)
                    else:
                        pass
                else:
                    pass
