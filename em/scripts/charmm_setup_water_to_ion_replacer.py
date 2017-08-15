#!/usr/bin/python
import sys
import os
import optparse
import numpy as np
from Bio.PDB.PDBParser import PDBParser

usage = "usage: %prog [options] -c arg1 -f arg2 -d arg3\n \
This file generates information to add salt ,in the form of Na and Cl ions,\n \
to a a box of water with a protein in it, and Box values for NAMD's PME parameters.\n"
option_parser = optparse.OptionParser(usage=usage)

option_parser.add_option('-c',
                        type='float',
                        help='Enter the charge of the system. (Get it from CHARMM\'s output to waterbox.inp.')

option_parser.add_option('-f',
                        type='str',
                        help='Enter the path to the pdb files with water molecules to be deleted.')
                        
option_parser.add_option('-d',
                        type='float',
                        help='Enter the salt concentration in miliMolars. Physiological = 0.145 mM.')
                        
options, args = option_parser.parse_args()
pdb_parser = PDBParser()

def gen_random_water_number():
    pass

def main():
    #if not args:
    #    print "Error: No input options provided. Type -h or --help."
    #    sys.exit(1)
    if not os.path.exists(options.f):
        print "Error: File path for molecule to be rotated does not exist."
        sys.exit(1)
    #pdb_file = pdb_parser.get_structure('waterbox', pdb_path)
    count_atoms = 0
    count_waters = 0
    with open(options.f) as f:
        content = f.readlines() 
    
    for i in content:
        print i
        j = i.split()
        if j[0] != '*':
            if len(j) == 1 or len(j) == 2:
                 print('Number of attoms',j[0])
            elif len(j) == 11:
                if j[2] == 'TIP3':
                    count_waters += 1
                if count_atoms == 0:
                    max_x = 0
                    max_y = 0
                    max_z = 0
                    min_x = 0
                    min_y = 0
                    min_z = 0
                X = float(j[4])
                Y = float(j[5])
                Z = float(j[6])
                if X > max_x:
                    max_x = X
                if X < min_x:
                    min_x = X
                if Y > max_y:
                    max_y = Y
                if Y < min_y:
                    min_y = Y
                if Z > max_z:
                    max_z = Z
                if Z < min_z:
                    min_z = Z
                count_atoms += 1
    count_waters = count_waters/3

    M = options.d*1e-27
    Num_molecA = M*6.022140857e23
    VOLUME = (max_x-min_x) * (max_y-min_y) * (max_z-min_z)
    center_x = min_x + (max_x - min_x)/2.0
    length_x = max_x - min_x    
    center_y = min_y + (max_y - min_y)/2.0
    length_y = max_y - min_y
    center_z = min_x + (max_z - min_z)/2.0
    length_z = max_z - min_z
    print('The following parameters are important for PME grid values in NAMD.')
    print('Volume:          ',VOLUME)
    print('Number of atoms in the solute:',count_atoms-(count_waters*3))
    print('Number of atoms in the solvent:',(count_waters*3))
    print('Total Number of atoms:',count_atoms)
    print('Numer of waters:',count_waters)
    print('Center  (x,y,z): ',center_x,center_y,center_z)
    print('Lengths (x,y,z):',length_x,length_y,length_z)

    number_of_NA = int(VOLUME*Num_molecA)
    number_of_CL = number_of_NA
    charge = options.c
    if charge > 0:
        number_of_NA = number_of_NA + charge
    elif charge < 0:
        number_of_CL = number_of_CL + charge
    list_wat_delete = np.random.rand(number_of_NA+number_of_CL)*count_waters
    resid_int = []
    for i in list_wat_delete:
        number = int(np.random.rand(1)*count_waters)
        if not number in resid_int:
            resid_int.append(number)
        else:
            is_in = True
            number = int(np.random.rand(1)*count_waters)
            while is_in:
                if not number in resid_int:
                    resid_int.append(number)
                    is_in = False
            
    if len(set(resid_int)) == len(resid_int):
        if len(resid_int) == (number_of_NA+number_of_CL):
            print('CHECKED: unique waters\' resid to delte.')
            stri = '('
            count = 1
            numbers_per_line = 6
            print('For NA atom removal')
            for i in range(0,number_of_NA):
                if i == (number_of_NA-1):
                    stri = stri + 'resid '+str(resid_int[i])
                else:
                    stri = stri + 'resid '+str(resid_int[i])+' .or. '
                if count % numbers_per_line == 0:
                    print(stri+' -')
                    stri = ''
                    count += 1
                else:
                    count += 1
            print(stri+') end')
            
            stri = '('
            count = 1
            numbers_per_line = 6
            print('For CL atom removal')
            for i in range(number_of_NA,number_of_NA+number_of_CL):
                if i == (number_of_NA+number_of_CL-1):
                    stri = stri + 'resid '+str(resid_int[i])
                else:
                    stri = stri + 'resid '+str(resid_int[i])+' .or. '
                if count % numbers_per_line == 0:
                    print(stri+' -')
                    stri = ''
                    count += 1
                else:
                    count += 1
            print(stri+') end')
        
if __name__ == '__main__':
    main()
