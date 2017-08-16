#!/usr/bin/python

"""
Created on Fri Jun 24 16:49:07 2016
@author: noel
For a description of the program, type:
python flower.py
"""
import os
import sys
import copy
import optparse
import string
import numpy as np
import em.describe.Molecular_Descriptors as md
import em.manipulate.Molecular_Rigid_Manipulations as MRM
import em.tools.CHARMM_Parser as CP
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Vector import *

def main():
    usage = "usage: %prog [options]"
    d = "It takes two pdbs and aligns the structures in different \
    orientations. A center-pdb is placed at (0,0,0),and a rotate-pdb \
    is place around at given angle intervals and distances. \
    The program assumes that the rotated protein's N-terminal will \
    become the N-Terminal when joined to the center protein. \
    WARNING: The programs needs to be extended to allow the reading \
    of the dipole vectors of the two proteins, and align them in any \
    way the user wants. As is, it centers the center-protein in whatever \
    orientation it is, and it alings the center of charge vector with \
    vector <0, 1, 0>. From that starting point, it rotates rotate-protein \
    by the angle given as input."
    option_parser = optparse.OptionParser(usage,description=d)
    option_parser.add_option("-c","--center", type="str", \
                             help="Path to pdb structure to place at (0,0,0) \
                             according to its center of mass.")
    option_parser.add_option("-r","--rotate", type="str", \
                             help="Path to pdb structure to rotated around \
                             center structure at X degree intervals. Angle \
                             intervals must be multiple of 360.")
    option_parser.add_option("-a","--angle", type="int", \
                             help="Angle intervals used to place rotated \
                             structured around center one.")
    option_parser.add_option("-d","--distance", type="float",\
                             help="Center of mass distance between rotated and \
                            centered structures.")
    option_parser.add_option('-i', '--id', type="str",action = "store", 
                             default = 's', help = "A prefix id to identify the \
                             output structures. Default is s.")    
    option_parser.add_option('-m', '--map', type="str",action = "store", 
                             default = 'no', help = "Optional arg if not present it \
                             is 'no', and it will add rotate with different chain \
                             identifiers. if it is 'yes' then it expects -l or --link \
                             to be something different from 'X' and in the right format.")
    option_parser.add_option('-l', '--link', type="str",action = "store", 
                             default = "X", help = "Optional arg if not present it \
                             is 'X', and it will add \'rotate\' structure with different chain \
                             identifiers. If it is 'yes', it will add the \'rotate\' structure \
                             with the same identifiers as \'center\' structure in the following \
                             format: A:A,B:B or chain A in \'rotate\' will link with A in \'center\' \
                             and B to B respectively. When both proteins have two chains, the only \
                             other option is A:B,B:A. There can be as many chains linked separated \
                             by commas but be careful because the program does not checks the \
                             physicality of this connections.")
    options, args = option_parser.parse_args()     
    if not os.path.exists(options.center):
        print "Error: File path for molecule to be centered does not exist."
        print("Type -h or --help for description and options.")
        sys.exit(1)        
    if not os.path.exists(options.rotate):
        print "Error: File path for molecule to be rotated does not exist."
        print("Type -h or --help for description and options.")
        sys.exit(1)
    if options.map.lower() == 'yes' and options.link.lower() == 'x':
        print("Error: when option map is equal to yes, option link must be")
        print("       a series of chains identifiers for association as described")
        print("       in the link paramater help.")
        sys.exit(1)
    pdb_parser = PDBParser()
    Angle = options.angle
    distance = options.distance
    #directory = os.path.dirname(options.center)
    filepath1 = options.center
    filepath2 = options.rotate
    map_o = options.map
    lnk_o = options.link
    ###########################################################################
    # Uncomment to test from spyder IDE
    #pdb_parser = PDBParser()
    #Angle = 45
    #distance = 45
    #file_name = os.path.basename(options.out).split('.')[0]
    #directory = '/home/noel/Projects/Protein_design/Constructs/Insulin_n_Barnase_barstar/'
    #filepath1 = directory+'2hiu_1rr_2.pdb'
    #filepath2 = directory+'1brs_comp_1rr.pdb'
    #map_o = "yes"
    #lnk_o = "A:C,B:F"
    ###########################################################################
    lnk_o = lnk_o.split(',')
    lnk_o = [i.split(':') for i in lnk_o]
    lnk_label = ''
    for i in lnk_o:
        lnk_label += '_'
        for j in i:
            lnk_label += j
    params = CP.read_charmm_FF(param_path)
    cmc = md.CenterOfMassCalculator(params)
    rig = MRM.Molecular_Rigid_Manipulation(param_path)
    s1 = pdb_parser.get_structure('Centered', filepath1)
    countS1 = 0
    modelS1 = -1
    for i in s1.get_models():
        countS1 += 1
        modelS1 = i.id
    if countS1 != 1:
        print("ERROR: Number of models cannot be different from 1. Models found:"+str(countS1))
        print("       Make sure Centered PDBs have only one model.")
        sys.exit(1)
    rig.translate_molecule(s1,modelS1,rig.center_molecule(cmc.get_center_of_mass(s1)))
    s2 = pdb_parser.get_structure('Rotated', filepath2) 
    countS2 = 0
    modelS2 = -1
    for i in s2.get_models():
        countS2 += 1
        modelS2 = i.id
    if countS2 != 1:
        print("ERROR: Number of models cannot be different from 1. Models found:"+str(countS2))
        print("       Make sure Rotated PDBs have only one model.")
        sys.exit(1)       
    rig.translate_molecule(s2,modelS2,rig.center_molecule(cmc.get_center_of_mass(s2)))
    #TODO: This works only for angles between 0 and 90 not including 0, and 90 and will generate angles in all 8 
    # quadrants of the cartesian coordinate system (I do not see why using an angle other than 45 for now.)
    locations = []
    angles = []
    for h in range(0,3):
        for i in range(0,360,Angle):
            for j in range(0,90/Angle-1):
                if h == 0:
                    z = 0
                    angles.append(str(i)+"_"+str(0))
                elif h == 1:
                    z = np.cos(45*np.pi/180)
                    angles.append(str(i)+"_"+str(45))
                elif h == 2:
                    z = -1*np.cos(45*np.pi/180)
                    angles.append(str(i)+"_"+str(315))
                locations.append([np.cos(i*np.pi/180),np.sin(i*np.pi/180),z])
    locations.append([0,0,1])
    angles.append(str(0)+"_"+str(90))
    locations.append([0,0,-1])
    angles.append(str(0)+"_"+str(270))
    for i in range(0,len(locations)):
        locations[i] = list(locations[i]/np.linalg.norm(locations[i]))
    ###############################################################################
    if map_o.lower() == 'no' and options.link.lower() == 'x':
        ids = {}
        for i in string.ascii_uppercase:
            ids[i] = False
        # First We used model 0 of structure 1 and turn ids for chains to True
        # identifier to True.
        for i in s1[0]:
            ids[i.id] = True
        # Now we go through structure 2 and if there are any chains with the same id
        # as those found in structure 1, we will change the chain ids to something else
        for i in s2[0]:
            if ids[i.id]:
                id_found = False
                for j in string.ascii_uppercase:
                    if not ids[j]:
                        i.id = j
                        ids[j] = True
                        # ID founds means id is unused so far and will be reserved and changed on S2
                        id_found = True
                        break
            # If the number of chains and identifiers exceeds letters in the alphabet,
            # It is necessary to modify the code. Until then, let's check this wont happen.
                if not id_found:
                    print("ERROR: Number of chains in both structures exceeds\n\
                           letters in the alphabet. No ID identifiers available.\n\
                           Program will exit without output. Fix the code.\n\
                           June 24, 2016")
                    sys.exit(1)
        s3 = copy.deepcopy(s1)
        s3.id = 'Ensamble'
        for i in s2.get_chains():
            s3[0].add(i)
    else:
        chain_info_s1 = {}
        for i in s1.get_chains():
            chain_info_s1[i.get_id()] = {}
            min_res = -1
            max_res = -1
            for j in i.get_residues():
                if min_res == -1:
                    min_res = j.get_id()[1]
                if j.get_id()[1] > max_res:
                    max_res = j.get_id()[1]
            chain_info_s1[i.get_id()]['min_res'] = min_res
            chain_info_s1[i.get_id()]['max_res'] = max_res
            # No need for min and max atoms.
            #min_atm = -1
            #max_atm = -1
            #for j in i.get_atoms():
            #    if min_atm == -1:
            #        min_atm = j.get_full_id()[3][1]
            #    if j.get_full_id()[3][1] > max_atm:
            #        max_atm = j.get_full_id()[3][1]
            #chain_info_s1[i.get_id()]['min_atm'] = min_atm
            #chain_info_s1[i.get_id()]['max_atm'] = max_atm  
        chain_info_s2 = {}
        for i in s2.get_chains():
            chain_info_s2[i.get_id()] = {}
            min_res = -1
            max_res = -1
            for j in i.get_residues():
                if min_res == -1:
                    min_res = j.get_id()[1]
                if j.get_id()[1] > max_res:
                    max_res = j.get_id()[1]
            chain_info_s2[i.get_id()]['min_res'] = min_res
            chain_info_s2[i.get_id()]['max_res'] = max_res
    ###############################################################################
    structure_id = 0
    for i in locations:
        if structure_id == 0:
            ccc = md.ChargeCalculator(params)
            cm = ccc.calculate_center_of_charge(s2)
            cm = cm/np.linalg.norm(cm)# ix orientation of LZ here only the first time.
            RM = rig.alignVectors(i,cm)
            for j in s2.get_atoms():
                v2 = [j.get_coord()[0],j.get_coord()[1],j.get_coord()[2]]
                jj = np.dot(v2,RM)
                j.set_coord(jj)
            # After aligning along centerofcharge/dipolemoment, the structure
            # is flipped to have the cterm closest to insulin. This only works with
            # LZ because it is a homodimer with both helices aligned in parallel.
            # For any other structure, this might not work. 
            # TODO Make this optional
            m = rotaxis(np.pi, Vector(0, 1, 0))
            for j2 in s2.get_atoms():
                v2 = Vector([j2.get_coord()[0],j2.get_coord()[1],j2.get_coord()[2]]) 
                v3 = v2.left_multiply(m)
                j2.set_coord(v3.get_array())
            last_direction = i
        else:
            # FIX: Something is wrong with rig.alignVectors([0,0,1],[0,0,-1])
            # the following inside the if statement patch just avoids the problem 
            # but does not solve it.
            if i == [0.0,0.0,-1.0] and last_direction == [0.0,0.0,1.0]:
                RM = rig.alignVectors([0.0,1.0,1.0],[0.0,0.0,1.0])
                for j in s2.get_atoms():
                    v2 = [j.get_coord()[0],j.get_coord()[1],j.get_coord()[2]]
                    jj = np.dot(v2,RM)
                    j.set_coord(jj)
                last_direction = [0.0,1.0,1.0]
            RM = rig.alignVectors(i,last_direction)
            for j in s2.get_atoms():
                v2 = [j.get_coord()[0],j.get_coord()[1],j.get_coord()[2]]
                jj = np.dot(v2,RM)
                j.set_coord(jj)
            last_direction = i
        # This 45 changes to explore different distaces between insulin and LZ
        ii = [k*distance for k in i]
        rig.translate_molecule(s2,modelS2,ii)
        # We want to join the two structures into one structure, with one model
        # and the chains of structure 1 and 2.
        if map_o.lower() == 'no' and lnk_o.lower() == 'x':
            pass
        else:
            s3 = s2.copy()
            # TODO: There is a bug. This algorithm cross conects two chains when both structures have the same
            # chain identifiers. When the chain identifiers are different, it fails. In such case. S3 remians identical
            # to S2.
            '''
            for h in lnk_o:
                for ii in s3.get_chains():
                    if ii.get_id() == h[0]:
                        for jj in s1.get_chains():
                            if jj.get_id() == h[1]:
                                count = chain_info_s2[ii.get_id()]['max_res'] + 1
                                for k in jj.get_residues():
                                    # first just do it for when they have same id
                                    k.id = (' ', count, ' ')
                                    ii.add(k)
                                    count += 1
            '''
            for h in lnk_o:
                change_rot_chain_id_to = h[0]
                chain_id_to_be_changed = h[1]
                # Go on the rotating chain and find the chain to be changed
                for iii in s3.get_chains():
                    if iii.get_id() == chain_id_to_be_changed: #C-C
                        count = chain_info_s2[iii.get_id()]['max_res'] + 1
                        for jjj in s1.get_chains():
                            if jjj.get_id() == change_rot_chain_id_to: #A-A
                                for kkk in jjj.get_residues():
                                    kkk.id = (' ', count, ' ')
                                    iii.add(kkk)
                                    count += 1
        io = PDBIO()
        io.set_structure(s3)
        io.save('s'+'_'+angles[structure_id]+"_"+str(structure_id)+lnk_label+'.pdb')
        #io.save(directory+options.id+'_'+angles[structure_id]+"_"+str(structure_id)+lnk_label+'.pdb')
        structure_id = structure_id + 1
        rig.translate_molecule(s2,modelS2,rig.center_molecule(cmc.get_center_of_mass(s2)))
if __name__ == '__main__':
    main()

