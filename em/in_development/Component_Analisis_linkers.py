#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 13:09:31 2017

@author: noel
"""
import numpy as np
import pandas as pd
from os import walk
import matplotlib.pyplot as plt

linker_path = "/home/noel/Projects/Protein_design/EntropyMaxima/examples/Linker_minimization/"
linkers = ["DDDDK","SGDDDDK","SGSGDDDDK"]
gbsa_AB = "gbsa"

MMGBSA_df = pd.DataFrame(columns=["dir","linker","c","r","ELECB","ELECU","GBENB","GBENU","VDWAB","VDWAU","ASPB","ASPU"])
MMGBSA_df.columns = ["dir","linker","c","r","ELECB","ELECU","GBENB","GBENU","VDWAB","VDWAU","ASPB","ASPU"]
MMGBSA_df_index = 0
for h in linkers:
    for i in range(1,(10+1)):
        for j in range(1,(1+1)):
            f = []
            for (dirpath, dirnames, filenames) in walk(linker_path+h+"/c_"+str(i)+"/r_"+str(j)+"/"):
                f.extend(filenames)
                break
            for k in dirnames:
                temp_df = pd.read_csv(dirpath+k+"/gbsw_ab.dat")
                temp_df.columns = [l.strip() for l in temp_df.columns]
                temp_df["dir"] = k
                temp_df["linker"] = h
                temp_df["c"] = i
                temp_df["r"] = j
                temp_df = temp_df[MMGBSA_df.columns]
                MMGBSA_df.loc[MMGBSA_df_index] = [m for m in temp_df.loc[0]]
                MMGBSA_df_index += 1

MMGBSA_df['dEE'] =  MMGBSA_df['ELECB'] - MMGBSA_df['ELECU']
MMGBSA_df['dGB'] =  MMGBSA_df['GBENB'] - MMGBSA_df['GBENU']
MMGBSA_df['dVW'] =  MMGBSA_df['VDWAB'] - MMGBSA_df['VDWAU']
MMGBSA_df['dSA'] =  MMGBSA_df['ASPB'] - MMGBSA_df['ASPU']
MMGBSA_df['dPolar'] =  MMGBSA_df['dEE'] + MMGBSA_df['dGB']
MMGBSA_df['dnonPolar'] =  MMGBSA_df['dVW'] + MMGBSA_df['dSA']
MMGBSA_df['MMGBSA'] = MMGBSA_df['dPolar'] + MMGBSA_df['dnonPolar']
############################################################################################
ddddk = MMGBSA_df.MMGBSA[MMGBSA_df.linker == "DDDDK"]
sgddddk = MMGBSA_df.MMGBSA[MMGBSA_df.linker == "SGDDDDK"]
sgsgddddk = MMGBSA_df.MMGBSA[MMGBSA_df.linker == "SGSGDDDDK"]
# from http://matplotlib.org/examples/statistics/boxplot_color_demo.html
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# rectangular box plot
bplot1 = axes.boxplot((ddddk,sgddddk,sgsgddddk),vert=True,patch_artist=True)   # vertical box aligmnent fill with color
plt.title('MMGBSA Binding Energy boxplots for three linkers')
plt.ylim([-200,200])
plt.setp(axes, xticks=[1,2,3],xticklabels=['DDDDK', 'SGDDDDK', 'SGSGDDDDK'])
plt.show()
##########################################################################################
Pddddk = MMGBSA_df.dPolar[MMGBSA_df.linker == "DDDDK"]
Psgddddk = MMGBSA_df.dPolar[MMGBSA_df.linker == "SGDDDDK"]
Psgsgddddk = MMGBSA_df.dPolar[MMGBSA_df.linker == "SGSGDDDDK"]
# from http://matplotlib.org/examples/statistics/boxplot_color_demo.html
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# rectangular box plot
bplot1 = axes.boxplot((Pddddk,Psgddddk,Psgsgddddk),vert=True, patch_artist=True)   # vertical box aligmnent fill with color
plt.title('Polar Binding Energy boxplots for three linkers')
#plt.ylim([-200,200])
plt.setp(axes, xticks=[1,2,3],xticklabels=['DDDDK', 'SGDDDDK', 'SGSGDDDDK'])
plt.show()
############################################################################################
NPddddk = MMGBSA_df.dnonPolar[MMGBSA_df.linker == "DDDDK"]
NPsgddddk = MMGBSA_df.dnonPolar[MMGBSA_df.linker == "SGDDDDK"]
NPsgsgddddk = MMGBSA_df.dnonPolar[MMGBSA_df.linker == "SGSGDDDDK"]
# from http://matplotlib.org/examples/statistics/boxplot_color_demo.html
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# rectangular box plot
bplot1 = axes.boxplot((NPddddk,NPsgddddk,NPsgsgddddk),vert=True,patch_artist=True)   # vertical box aligmnent fill with color
plt.title('NonPolar Binding Energy boxplots for three linkers')
#plt.ylim([-200,200])
plt.setp(axes, xticks=[1,2,3],
         xticklabels=['DDDDK', 'SGDDDDK', 'SGSGDDDDK'])
plt.show()
######################################################################################
### Find most favorable energy for each linker
min_linker_dict = {}
for i in linkers:
    min_linker_dict[i] = [-1,1000000,-1,-1,""]
for i in linkers:
    for j in MMGBSA_df.index[MMGBSA_df.linker == i]:
        if MMGBSA_df.loc[j,'MMGBSA'] < min_linker_dict[i][1]:
            min_linker_dict[i][0] = j
            min_linker_dict[i][1] = MMGBSA_df.loc[j,'MMGBSA']
            min_linker_dict[i][2] = MMGBSA_df.loc[j,'c']
            min_linker_dict[i][3] = MMGBSA_df.loc[j,'r']
            min_linker_dict[i][4] = MMGBSA_df.loc[j,'dir']
    print("Minimum Energy for Linker "+i+" = "+str(min_linker_dict[i][1])+" index = "+str(min_linker_dict[i][0]))
    print("        center structure = "+str(min_linker_dict[i][2])+" rotated structure = "+str(min_linker_dict[i][3]))
    print("               directory = "+min_linker_dict[i][4])
