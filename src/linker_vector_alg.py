# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 00:16:21 2016

@author: noel
"""

from Bio.PDB import *

parser = PDBParser()

structure = parser.get_structure('Insulin', '/home/noel/Projects/Protein_design/Insuline/2hiu.pdb')

for model in structure:
    for chain in model:
        print chain
        for residue in chain:
            print residue
            for atom in residue:
                print atom

# Iterate over all atoms in a structure
counter = 0
for atom in structure.get_atoms():
    print(counter,atom)
    counter = counter + 1
 
# Iterate over all residues in a model
counter = 0
for residue in model.get_residues():
    print(counter,residue)
    counter = counter + 1
    
model = structure[0]
chain = model['A']
residue = chain[1]
atom = residue['CA']

atom = structure[0]['A'][10]['CA']

atom.get_name()           # atom name (spaces stripped, e.g. 'CA')
atom.get_id()             # id (equals atom name)
atom.get_coord()          # atomic coordinates
atom.get_vector()         # atomic coordinates as Vector object
atom.get_bfactor()        # isotropic B factor
atom.get_occupancy()      # occupancy
atom.get_altloc()         # alternative location specifier
atom.get_sigatm()         # std. dev. of atomic parameters
atom.get_siguij()         # std. dev. of anisotropic B factor
atom.get_anisou()         # anisotropic B factor
atom.get_fullname()       # atom name (with spaces, e.g. '.CA.')

# Get some atoms and distance between
ca1 = residue1['CA']
ca2 = residue2['CA']
# Simply subtract the atoms to get their distance
distance = ca1-ca2

# Get angles
vector1 = atom1.get_vector()
vector2 = atom2.get_vector()
vector3 = atom3.get_vector()
angle = calc_angle(vector1, vector2, vector3)

# Get torsion angles
vector1 = atom1.get_vector()
vector2 = atom2.get_vector()
vector3 = atom3.get_vector()
vector4 = atom4.get_vector()
angle = calc_dihedral(vector1, vector2, vector3, vector4)

#eAngx = 0.0;      eAngy = 0.0;       eAngz = 0.0;
#a1 = 0.5 * eAngy;
#a2 = 0.5 * (eAngx - eAngz);
#a3 = 0.5 * (eAngx + eAngz);
#System.out.println("   a123 = "+a1+" "+a2+" "+a3);
#TM.q_u1[0] = Math.sin(a1) * Math.cos(a2);
#TM.q_u2[0] = Math.sin(a1) * Math.sin(a2);
#TM.q_u3[0] = Math.cos(a1) * Math.sin(a3);
#TM.q_u4[0] = Math.cos(a1) * Math.cos(a3);
#System.out.println("   q1234 = "+TM.q_u1[0]+" "+TM.q_u2[0]+" "+TM.q_u3[0]+" "+TM.q_u4[0]);
#tq[0] = TM.q_u1[0];  tq[1] = TM.q_u2[0];  tq[2] = TM.q_u3[0];    tq[3] = TM.q_u4[0];
#for(k = 0, k2 = 0; k2 < 4; k2++){
#    for(k1 = k2; k1 < 4; k1++, k++){
#        p[k] = 2.0*tq[k1]*tq[k2];
#    }
#}
#TM.rMatT[0] = p[0] + p[9] - 1;  TM.rMatT[4] = p[4] + p[9] - 1;   TM.rMatT[8] = p[7] + p[9] - 1;
#s = 1.0;    //Transpose = 1
#TM.rMatT[1] = p[1] + s * p[8];   TM.rMatT[3] = p[1] - s * p[8];   TM.rMatT[2] = p[2] - s * p[6];
#TM.rMatT[6] = p[2] + s * p[6];   TM.rMatT[5] = p[5] + s * p[3];   TM.rMatT[7] = p[5] - s * p[3];
