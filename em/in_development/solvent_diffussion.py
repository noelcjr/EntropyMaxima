#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 16:18:10 2018

@author: noel
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import Bio.PDB as struct
import em.describe.Molecular_Descriptors as md
import em.tools.CHARMM_Parser as CP
import em.manipulate.Molecular_Rigid_Manipulations as MRM

base_dir = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/'
input_files = ['1B3S/brn_H102A_brs_Y29F/stripped_structures/input_gbsw_ab.dat',
               '1B3S/brn_H102A_brs_Y29F_C41A_C83A/stripped_structures/input_gbsw_ab.dat',
               '1B27/as_is/stripped_structures/input_gbsw_ab.dat',
               '1B27/brs_A41C_A83C/stripped_structures/input_gbsw_ab.dat']
               
# Simple data to display in various forms
cols = ['indx','Eng','ee','gb','vdw','sa']
MMGBSA_enrg = pd.DataFrame(columns=['1B3S_2','1B3S_4','1B27_0','1B27_2'])
df1B3S_brn_H102A_brs_Y29F = pd.read_csv(base_dir+input_files[0], names = cols)
df1B3S_brn_H102A_brs_Y29F = df1B3S_brn_H102A_brs_Y29F.set_index('indx')
df1B3S_brn_H102A_brs_Y29F_C41A_C83A = pd.read_csv(base_dir+input_files[1], names = cols)
df1B3S_brn_H102A_brs_Y29F_C41A_C83A = df1B3S_brn_H102A_brs_Y29F_C41A_C83A.set_index('indx')
df1B27_as_is = pd.read_csv(base_dir+input_files[2], names = cols)
df1B27_as_is = df1B27_as_is.set_index('indx')
df1B27_brs_A41C_A83C = pd.read_csv(base_dir+input_files[3], names = cols)
df1B27_brs_A41C_A83C = df1B27_brs_A41C_A83C.set_index('indx')

MMGBSA_enrg['1B3S_2'] = df1B3S_brn_H102A_brs_Y29F['Eng']
MMGBSA_enrg['1B3S_4'] = df1B3S_brn_H102A_brs_Y29F_C41A_C83A['Eng']
MMGBSA_enrg['1B27_0'] = df1B27_as_is['Eng']
MMGBSA_enrg['1B27_2'] = df1B27_brs_A41C_A83C['Eng']

MMGBSA_enrg.describe()
MMGBSA_enrg.plot()
MMGBSA_enrg.boxplot()
MMGBSA_enrg.loc[0:500].boxplot()
MMGBSA_enrg.loc[500:].boxplot()

#dir_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1B3S/brn_H102A_brs_Y29F/structures/'
dir_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1B27/brs_A41C_A83C/structures/'
residue_list = []
count_for_diffusion = 0
i = 0
structure = struct.PDBParser(QUIET=True).get_structure('pdb', dir_path+str(i)+'.pdb') 
for j in structure.get_residues():
    if j.get_resname() in ['TIP','SOD','CLA']:
        if j.get_resname() == 'CLA':
            print(j.get_full_id())
        count_for_diffusion += 1
        residue_list.append(j.get_full_id()[2]+str(j.get_full_id()[3][1]))
# TODO Need to rename all CLA chains to something different from I(SOD)
diffusion = pd.DataFrame(index=residue_list)
count_frames = 0
for i in range(0,990):
    vector_list = []
    structure1 = struct.PDBParser(QUIET=True).get_structure('pdb', dir_path+str(i)+'.pdb')
    for j in structure1.get_residues():
        if j.get_resname() in ['TIP','SOD','CLA']:
            #indx = j.get_full_id()[2]+str(j.get_full_id()[3][1])
            for k in j.get_atoms():
                if k.get_id() in ['OH2','SOD','CLA']:
                    vector_list.append(k.get_coord())
    diffusion[str(count_frames)] = list(np.array(vector_list))
    print(count_frames)
    count_frames += 1
   
length = len(diffusion.ix['B111'])
# Delta of translation without correction for image transform
mean_diff = []
std_diff = []
for i in diffusion.index:
    from_ix = 0
    to_ix = length
    diffusion_per_10ps = []
    for j in range(len(diffusion.columns)-1):
        diffusion_per_10ps.append(np.linalg.norm(diffusion.ix[i][str(j+1)]-diffusion.ix[i][str(j)]))
    print(i)
    mean_diff.append(np.mean(diffusion_per_10ps))
    std_diff.append(np.std(diffusion_per_10ps))
count = 0
for i in diffusion.index:
    diffusion.ix[i,'mean'] = mean_diff[count]
    diffusion.ix[i,'std'] = std_diff[count]
    count += 1
diffusion.to_csv(dir_path+'diffusion.csv')
# get waters with mean diffusions of less than
diffusion = pd.read_csv(dir_path+'diffusion.csv')
###############################################################
# Delta of translation with correction for image transform
protein_XYZ_length = []
for i in range(0,990):
    atoms = []
    structure = struct.PDBParser(QUIET=True).get_structure('pdb', dir_path+str(i)+'.pdb') 
    for k in structure.get_atoms():
        atoms.append(k.get_coord())
    x,y,z = zip(*atoms)
    X = max(x) - min(x)
    Y = max(y) - min(y)
    Z = max(z) - min(z)
    protein_XYZ_length.append([(max(x),min(x)),(max(y),min(y)),(max(z),min(z))])
    print(i)

dir_path_stripped = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1B27/brs_A41C_A83C/stripped_structures/'   
protein_CM = []
params = CP.read_charmm_FF()
cmc = md.CenterOfMassCalculator(params)
for i in range(0,989):
    structure = struct.PDBParser(QUIET=True).get_structure('pdb', dir_path_stripped+str(i)+'.pdb') 
    rig = MRM.Molecular_Rigid_Manipulation(params)
    protein_CM.append(cmc.get_center_of_mass(structure))

np_delta = np.zeros((diffusion.shape[0],diffusion.shape[1]-3))
np_delta_CM = np.zeros((diffusion.shape[0],diffusion.shape[1]-3))
for i in range(0,diffusion.shape[1]-3):
    # This only works when the center of the box is at (0,0,0) for CHARMM not for Gromacs
    x_q = np.abs(protein_XYZ_length[i][0][1])+np.abs(protein_XYZ_length[i][0][0])
    max_diff_x = (x_q/2)
    y_q = np.abs(protein_XYZ_length[i][1][1])+np.abs(protein_XYZ_length[i][1][0])
    max_diff_y = (y_q/2)
    z_q = np.abs(protein_XYZ_length[i][2][1])+np.abs(protein_XYZ_length[i][2][0])
    max_diff_z = (z_q/2)
    # anything outside the two inner halves need to be corrected for image transform. In other words:
    # if a water molecule difuses more than a quarter of the box, it gets corrected for image convertion.
    # Warming: there might be temperatures and pressures where this diffusion is common, and not due to image
    # conversion, and could cause errors. to check for that, plot diffusion[a]-diffusion[b] histograms for X,Y,Z
    # to check the range of displacements in each dimension.
    a = diffusion.columns[i+1]
    b = diffusion.columns[i]
    delta_xyz = []
    j_count = 0
    for j in list(diffusion[a]-diffusion[b]):
        displacement = abs(j)
        if displacement[0] > max_diff_x:
            displacement[0] = x_q - displacement[0]
        if displacement[1] > max_diff_y:
            displacement[1] = y_q - displacement[1]
        if displacement[2] > max_diff_z:
            displacement[2] = z_q - displacement[2]
        np_delta[j_count,i] = np.linalg.norm(displacement)
        np_delta_CM[j_count,i] = np.linalg.norm(j-np.array(protein_CM[i]))
        j_count += 1
    print(i)

fig = plt.gcf()
fig.set_size_inches(7.297, 0.990)
plt.imshow(np_delta, aspect='auto')
plt.show()
plt.imshow(np_delta_CM, aspect='auto')
plt.show()

np.savetxt(dir_path+'np_delta.txt', np_delta)
np.savetxt(dir_path+'np_delta_CM.txt', np_delta_CM)
#b = np.loadtxt(dir_path+'np_delta.txt', dtype=float)

a = np_delta.flat
plt.hist(a, bins=[i for i in range(int(np.ceil(np.max(a))))])
b = np_delta_CM.flat
plt.hist(b, bins=[i for i in range(int(np.ceil(np.max(b))))])

total_displacement = []
mean_displacement = []
for i in range(np_delta.shape[0]):
    suma = np.sum(np_delta[i][:])
    total_displacement.append(suma)
    mean_displacement.append(suma/np_delta.shape[1])
plt.plot(total_displacement)
plt.show()
plt.plot(mean_displacement)
plt.hist(mean_displacement)
plt.hist(total_displacement)
# Todo: search for waters that displace less than 1A at least once and count them
#       search for waters that displace less than 1A at least twice consecutively and count them
#       search for waters that displace less than 1A at least three times consecutively and count them and so on... 
####################################################
# ok up to here
array_width = len(np_delta[0])
counts_per_frame = np.zeros((np_delta.shape[0],11))
for i in range(np_delta.shape[0]):
    for j in range(array_width):
        if np_delta[i][j] <= 0.5:
            counts_per_frame[i][0] += 1
        elif (np_delta[i][j] > 0.5) and (np_delta[i][j] <= 1.0):
            counts_per_frame[i][1] += 1
        elif (np_delta[i][j] > 1.0) and (np_delta[i][j] <= 1.5):
            counts_per_frame[i][2] += 1
        elif (np_delta[i][j] > 1.5) and (np_delta[i][j] <= 2.0):
            counts_per_frame[i][3] += 1
        elif (np_delta[i][j] > 2.0) and (np_delta[i][j] <= 2.5):
            counts_per_frame[i][4] += 1
        elif (np_delta[i][j] > 2.5) and (np_delta[i][j] <= 3.0):
            counts_per_frame[i][5] += 1
        elif (np_delta[i][j] > 3.0) and (np_delta[i][j] <= 3.5):
            counts_per_frame[i][6] += 1
        elif (np_delta[i][j] > 3.5) and (np_delta[i][j] <= 4.0):
            counts_per_frame[i][7] += 1
        elif (np_delta[i][j] > 4.0) and (np_delta[i][j] <= 4.5):
            counts_per_frame[i][8] += 1
        elif (np_delta[i][j] > 4.5) and (np_delta[i][j] <= 5.0):
            counts_per_frame[i][9] += 1
        elif (np_delta[i][j] > 5.0):
            counts_per_frame[i][10] += 1

fig = plt.gcf()
fig.set_size_inches(7.297, 0.990)
plt.imshow(counts_per_frame, aspect='auto')

fixed_waters = {}
for i in diffusion.index:
    fixed_waters_index = []
    if diffusion.loc[i,'mean'] < 8:
        if not i[0] in fixed_waters:
            fixed_waters[i[0]] = []
            fixed_waters[i[0]].append((i[1:],diffusion.ix[i]['mean']))
        else:
            fixed_waters[i[0]].append((i[1:],diffusion.ix[i]['mean']))

list_of_non_movable_W = []
for i in fixed_waters['W']:
    if i[1] < 2:
        list_of_non_movable_W.append(i)
        print(i[0])
# delta numpy plot
count = 0
for i in (diffusion['0']-diffusion['1']):
    if i[2] >20:
        print(i)
        count += 1
np_delta = np.zeros((diffusion.shape[0],diffusion.shape[1]-2))
i_count = 0
for i in range(np_delta.shape[1]-1):
    a = diffusion.columns[i+1]
    b = diffusion.columns[i]
    resultA = diffusion[a].apply(lambda x:np.fromstring(x.replace('[','').replace(']','').replace('  ',''), sep=''))
    resultB = diffusion[b].apply(lambda x:np.fromstring(x.replace('[','').replace(']','').replace('  ',''), sep=''))
    resultC = resultA - resultB
    j_count = 0
    for j in resultC:
        np_delta[j_count,i_count] = np.linalg.norm(j)
        j_count += 1
    print(i_count)
    i_count += 1
     