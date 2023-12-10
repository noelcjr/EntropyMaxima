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
import em.code.Molecular_Descriptors as md
import em.code.Molecular_Rigid_Manipulations as MRM
import em.code.CHARMM_Parser as CP
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Vector import *
import Bio.PDB as struct

def join_chains(center_struct,rotate_struct,joint_struct,fusion_order,linker):
    """ Join chains in to different structures and adds them to a third structure for output. This can only join two
        chains. TO join more chains, just add more chains to already joined chains by calling this function multiple times. All
        amino acids in the chain will be joined. Chains are joined by sequence number, but the chains might be sepa-
        rated, and they might require minimization or translation to join.
        first_struct:       Structure where the joininig chains will be together. New Chain.
        begining_struct:    This one will be the one rotated.
        end_struct:         This one will be at the centers.
        linker              linker. Describes chains to connect and the connencting terminals in any order.
    """
    #join_chains(s3,s2,s1,lnk_o[0],join_range1)
    linker = link_o.split(':')
    linker = [i.split(',') for i in linker]
    lnk_label = ''
    for i in linker:
        lnk_label += '_'
        for j in i:
            lnk_label += j.lower()
    chain_info_s1 = get_chains_info(s1)#center_struct)
    chain_info_s2 = get_chains_info(s2)#rotate_struct)
    chains_range = []
    chains_range.append([(chain_info_s1[linker[0][0][0]]['min_res'],chain_info_s1[linker[0][0][0]]['max_res']),
                        (chain_info_s2[linker[0][1][0]]['min_res'],chain_info_s2[linker[0][1][0]]['max_res'])])
    chains_range.append([(chain_info_s1[linker[1][0][0]]['min_res'],chain_info_s1[linker[1][0][0]]['max_res']),
                        (chain_info_s2[linker[1][1][0]]['min_res'],chain_info_s2[linker[1][1][0]]['max_res'])])
    if len(linker) != 2:
        print("ERROR: There can only be two chains in the chains_to_connect list. Program exits with an ERROR.")
        sys.exit(1)
    cntr = 0
    for i in linker:
        add_center_chain_res = []
        center_chain, center_direction = i[0].split('.')
        for ii in center_struct[0].get_chains():
            if ii.get_id() == center_chain:
                for jj in ii.get_residues(): 
                    add_center_chain_res.append(jj.copy())
        if center_direction == 'r':
            add_center_chain_res = add_center_chain_res[::-1]
        elif center_direction == 'f':
            pass
        else:
            print("Error: Wrong option for chain's amino acid order entered. Only 'f' for forward, or 'r' for reverse.")
            sys.exit(1)
        add_rotate_chain_res = []       
        rotate_chain, rotate_direction = i[1].split('.')
        for ii in rotate_struct[0].get_chains():
            if ii.get_id() == rotate_chain:
                for jj in ii.get_residues(): 
                    add_rotate_chain_res.append(jj.copy())
        if rotate_direction == 'r':
            add_rotate_chain_res = add_rotate_chain_res[::-1]
        elif rotate_direction == 'f':
            pass
        else:
            print("Error: Wrong option for chain's amino acid order. Only 'f' for forward, or 'r' for reverse.")
            sys.exit(1)
        newChain = struct.Chain.Chain(linker[cntr][0][0])
        count = 1
        if fusion_order == 'CR': 
            add_res = add_center_chain_res + add_rotate_chain_res
        elif fusion_order == 'RC':
            add_res = add_rotate_chain_res + add_center_chain_res
        else:
            print("Error: Wrong option for fusion order. Only 'CR' or 'RC', for Center to Rotate or viceversa")
        for j in add_res:
            #if j.get_id() != (' ',count,' '):
            j.id = (' ',count,' ')
            newChain.add(j)
            count += 1
        for ii in newChain.get_residues():
            if ii.get_parent().get_id() != rotate_chain:
                ii.get_parent().id = rotate_chain
        joint_struct[0].add(newChain)
        cntr += 1
    return ''.join(lnk_label.split('.'))


def get_chains_info(struct):
    """For now it only gets the maximum and minimum amino acid for each chain, but we could give out more info later.
       This function has been moved to structure.py and enhanced there to include multiple models.
    """
    dictionary_info = {}
    for i in struct.get_chains():
        dictionary_info[i.get_id()] = {}
        min_res = -1
        max_res = -1
        for j in i.get_residues():
            # We assume that residue numbers in chains monotonically increase so min_res is assigned the first residue
            if min_res == -1:
                min_res = j.get_id()[1]
            if j.get_id()[1] > max_res:
                max_res = j.get_id()[1]
        dictionary_info[i.get_id()]['min_res'] = min_res
        dictionary_info[i.get_id()]['max_res'] = max_res 
    return dictionary_info
def main():
    usage = "usage: %prog [options]"
    d = "It takes two pdbs and aligns the structures in different orientations. A center-pdb is placed at (0,0,0),and \
    a rotate-pdb is place around center at given angle intervals and distances. You can choose to join chains to make \
    a fusion peptide, or, if no, you can just align proteins for searching protein-protein interactions."
    option_parser = optparse.OptionParser(usage,description=d)
    option_parser.add_option("--center", type="str", \
                             help="Path to pdb structure to place at (0,0,0) according to its center of mass.")
    option_parser.add_option("--rotate", type="str", \
                             help="Path to pdb structure to rotated around center structure at X degree intervals. \
                             Angle intervals must be multiple of 360.")
    option_parser.add_option("--angle", type="int", \
                             help="Angle intervals used to place rotated structured around center one.")
    option_parser.add_option("--distance", type="float",\
                             help="Center of mass distance between rotated and centered structures.")
    option_parser.add_option('--id', type="str",action = "store", 
                             default = 's', help = "A prefix id to identify the output structures. Default is s.")    
    option_parser.add_option('--link', type="str",action = "store", default = 'x', \
                             help = "Option only works with a 'yes' value from the --join option. It will join chains \
                             in the following format, for example: 'A.r,B.f:B.f,A.f'. This --link option can be read as\
                             the centered structure's chain 'A' amino acid sequence numbering is reversed and joind to \
                             the rotated structure's chain B with amino acid sequenced numbers not reversed, or kept \
                             forward. The values after the column are read similarly. Default if missing: 'x'")
    option_parser.add_option('--join', type="str",action = "store", default = 'no', \
                             help = "Only two possible options 'CR' and 'RC'. In both options C stands for Center, and \
                             R for rotate. CR means that the numbering starts at the first amino acid of --center \
                             structure up to its end, and it is followed by the first amino acid of --rotate up to its \
                             end. This option might overide the directions (r or f) of amino acid numbering for each\
                             chain in --link option. This gives total manipulation flexibility for joining chains.")
    options, args = option_parser.parse_args()     
    if not os.path.exists(options.center):
        print("Error: File path for molecule to be centered does not exist.")
        print("Type -h or --help for description and options.")
        sys.exit(1)        
    if not os.path.exists(options.rotate):
        print("Error: File path for molecule to be rotated does not exist.")
        print("Type -h or --help for description and options.")
        sys.exit(1)
    if options.map.lower() == 'yes' and options.link.lower() == 'x':
        print("Error: when option map is equal to yes, option link must be")
        print("       a series of chains identifiers for association as described")
        print("       in the link paramater help.")
        sys.exit(1)
    pdb_parser = PDBParser(QUIET = True)
    Angle = options.angle
    distance = options.distance
    #directory = os.path.dirname(options.center)
    filepath1 = options.center
    filepath2 = options.rotate
    join_o = options.join
    link_o = options.link
    ###########################################################################
    # Uncomment to test from spyder IDE
    pdb_parser = PDBParser(QUIET = True)
    Angle = 45
    distance = 45
    # file_name = os.path.basename(options.out).split('.')[0]
    directory1 = "/home/noel/Projects/Protein_design/Insulin/struct_prep/2hiu/init_setup/"
    filepath1 = directory1+'2hiu_1rr.pdb'
    directory2 = "/home/noel/Projects/Protein_design/Leucine_zipper/Struct_Prep/init_setup/"
    filepath2 = directory2+'2zta_1rr.pdb'
    # param_path = "/home/noel/Projects/Protein_design/EntropyMaxima/params/charmm27.ff/"
    #join_o = "RC"
    #link_o = "A.f,A.f:B.f,B.f"
    ###########################################################################
    # Process strig that the determines how the centered and rotated structures
    # will be connected.
    params = CP.read_charmm_FF()
    cmc = md.CenterOfMassCalculator(params)
    rig = MRM.Molecular_Rigid_Manipulation(params)
    ###########################################################################
    # Check that the structures only have one model, and Place the structures' 
    # center of mass at (0,0,0) to give an idea of their location in the 
    # cartesian coordinate system 
    s1 = pdb_parser.get_structure('Centered', filepath1)
    countS1 = 0
    modelS1 = -1
    for i in s1.get_models():
        countS1 += 1
        modelS1 = i.id
    if countS1 != 1:
        print("ERROR: Number of models cannot be different from 1.")
        print("       Models found:"+str(countS1))
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
        print("ERROR: Number of models cannot be different from 1.")
        print("       Models found:"+str(countS2))
        print("       Make sure Rotated PDBs have only one model.")
        sys.exit(1)
    rig.translate_molecule(s2,modelS2,rig.center_molecule(cmc.get_center_of_mass(s2)))
    ###########################################################################
    # TODO: This works only for angles between 0 and 90 not including 0, and 90
    # and will generate angles in all 8 quadrants of the cartesian coordinate 
    # system (I do not see why using an angle other than 45 for now.)
    # TODO: quaternions might work better.
    # The location list has a list of normalized vectors releative to (0,0,0)
    # that will be use to place s2's in the right orientation relative to s1
    locations = []
    angles = []
    count = 0
    for h in range(0,3):
        for i in range(0,360,Angle):
            for j in range(0,int(90/Angle-1)):
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
                print(count,h,i,j,angles[count],locations[count], np.linalg.norm(locations[count]))
                count += 1
    locations.append([0,0,1])
    angles.append(str(0)+"_"+str(90))
    locations.append([0,0,-1])
    angles.append(str(0)+"_"+str(270))
    for i in range(0,len(locations)):
        locations[i] = list(locations[i]/np.linalg.norm(locations[i]))
        print(np.linalg.norm(locations[i]))
    ###########################################################################
    # This works when you do not want to link two proteins but instead you want
    # to place one around the other and check
    # protein-protein interaction's binding affinity.
    if join_o.lower() == 'no' and link_o.lower() == 'x':
        ids = {}
        for i in string.ascii_uppercase:
            ids[i] = False
        # First We used model 0 of structure 1 and turn ids for chains to True 
        # identifier to True.
        for i in s1[0]:
            ids[i.id] = True
        # Now we go through structure 2 and if there are any chains with the 
        # same id as those found in structure 1, we will change the chain ids 
        # to something else becase if chain ids are repeated in the same 
        # structure it will be considered one chain when they are actually 
        # separated.
        for i in s2[0]:
            if ids[i.id]:
                id_found = False
                for j in string.ascii_uppercase:
                    if not ids[j]:
                        i.id = j
                        ids[j] = True
                        # ID founds means id is unused so far and will be 
                        # reserved and changed on S2
                        id_found = True
                        break
                # If the number of chains and identifiers exceeds letters in 
                # the alphabet, It is necessary to modify the code. Until then,
                # let's check this won't happen by exiting.
                if not id_found:
                    print("ERROR: Number of chains in both structures exceeds \
                           letters in the alphabet. No ID identifiers availabl\
                           e. Program will exit without output. Fix the code. \
                           June 24, 2016")
                    sys.exit(1)
        s3 = copy.deepcopy(s1)
        s3.id = 'Ensamble'
        for i in s2.get_chains():
            s3[0].add(i)
    ###########################################################################
    structure_id = 0
    for i in locations:
        if structure_id == 0:
            ccc = md.ChargeCalculator(params)
            cm = ccc.calculate_center_of_charge(s2)
            cm = cm/np.linalg.norm(cm)
            # ix orientation of LZ here only the first time.
            RM = rig.alignVectors(i,cm)
            for j in s2.get_atoms():
                v2 = [j.get_coord()[0],j.get_coord()[1],j.get_coord()[2]]
                jj = np.dot(v2,RM)
                j.set_coord(jj)
            last_direction = i
            # After aligning along centerofcharge/dipolemoment, the structure
            # is flipped to have the cterm closest to insulin. This only works 
            # with LZ because it is a homodimer with both helices aligned in 
            # parallel. For any other structure, this might not work. 
            # TODO This Vector needs to be picked by the user from information 
            #      from the pdb_cif.py --summary
            # 
            # For some reason after modifications to flower the next comented 
            # out code in this if step does not seem necessery
            #m = rotaxis(np.pi, Vector(0, 1, 0))
            #for j2 in s2.get_atoms():
            #    v2 = Vector([j2.get_coord()[0],j2.get_coord()[1],j2.get_coord()[2]]) 
            #    v3 = v2.left_multiply(m)
            #    j2.set_coord(v3.get_array())
        else:
            # FIX: Something is wrong with rig.alignVectors([0,0,1],[0,0,-1]) 
            # it could just be a trigonometry case that gives some sort of 
            # singularity. The following code inside the if statement patch 
            # just avoids the problem but does not explain it. Find out.
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
        ###
        # We want to join the two structures into one structure, with one model
        # and the chains of structure 1 and 2 joined and named consistently.
        if join_o.lower() == 'no' and link_o.lower() == 'x':
            pass
        else:
            s3 = struct.StructureBuilder.Structure('newrot')
            s3.add(struct.Model.Model(0))
            lnk_label = join_chains(s1,s2,s3,join_o,link_o)
        io = PDBIO()
        io.set_structure(s3)
        io.save('s'+'_'+angles[structure_id]+"_"+str(structure_id)+lnk_label+'.pdb')
        #io.save(directory+options.id+'_'+angles[structure_id]+"_"+str(structure_id)+lnk_label+'.pdb')
        structure_id = structure_id + 1
        rig.translate_molecule(s2,modelS2,rig.center_molecule(cmc.get_center_of_mass(s2)))
if __name__ == '__main__':
    main()

