# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 12:57:28 2016

@author: noel
"""
import os
import sys
import CHARMM_Parser as CP
#import RTPParser as rp

class CenterOfMassCalculator(object):
    def __init__(self,charmm_params):
        self.mx_total = 0
        self.my_total = 0
        self.mz_total = 0
        self.total_mass = 0
        self.charmm = charmm_params
# TODO: Atom mass should be obtained from CHARMM parameters, BIO.PDB
    def add(self, atom):
        coords = atom.coord.tolist()
        self.total_mass += atom.mass
        self.mx_total += coords[0] * atom.mass
        self.my_total += coords[1] * atom.mass
        self.mz_total += coords[2] * atom.mass
    
    def get_center_of_mass(self, structure):
        self.mx_total = 0
        self.my_total = 0
        self.mz_total = 0
        self.total_mass = 0  
        for atom in structure.get_atoms():
            self.add(atom)
        return self.center_of_mass()

    def center_of_mass(self):
        if self.total_mass == 0:
            raise Exception("total mass is zero")
        return [self.mx_total/self.total_mass, self.my_total/self.total_mass, self.mz_total/self.total_mass]

def _get_residue_data(charmmdir):
    aminoacids_path = os.path.join(charmmdir, 'aminoacids.rtp')
    if not os.path.exists(aminoacids_path):
        print "Error: Failed to find aminoacids.rtp file."
        sys.exit(1)
    return rp.RTP(aminoacids_path)

def _do_histidine_magic(res, rtp):
    return rtp.residues.get('HSD')

def _do_atom_magic(rtp_residue, atom):
    name = atom.name
    replacement_atom = None
    if rtp_residue.name == 'ACE':
        if name == 'H1':
            replacement_atom = rtp_residue.atoms.get('HH31')
        elif name == 'H2':
            replacement_atom = rtp_residue.atoms.get('HH32')
        elif name == 'H3':
            replacement_atom = rtp_residue.atoms.get('HH33')
    elif name == 'H':
        replacement_atom = rtp_residue.atoms.get('HN')
    elif name == 'HA3':
        replacement_atom = rtp_residue.atoms.get('HA1')
    elif name == 'HB3':
        replacement_atom = rtp_residue.atoms.get('HB1')
    elif name == 'HD3':
        replacement_atom = rtp_residue.atoms.get('HD1')
    elif name == 'HE3':
        replacement_atom = rtp_residue.atoms.get('HE1')
    elif name == 'HG3':
        replacement_atom = rtp_residue.atoms.get('HG1')
    elif rtp_residue.name == 'SER' and name == 'HG':
        replacement_atom = rtp_residue.atoms.get('HG1')
    return replacement_atom

class ChargeCalculator(object):
    def __init__(self, charmm_params):
        self.qx_total = 0
        self.qy_total = 0
        self.qz_total = 0
        self.q_total = 0
        self.charmm = charmm_params
    
    def calculate_center_of_charge(self, structure):
        self.qx_total = 0
        self.qy_total = 0
        self.qz_total = 0
        self.q_total = 0

        for atom in structure.get_atoms():
            if atom.get_parent().get_resname() == 'HOH':
                continue      # To next item in the loop
            #self.add(atom)
            name = atom.get_name()
            # This is a fix because terminals are part of a residue but has atoms
            # that are not found in the parameter for those atoms.
            if name in ['OT1','OT2']:
                res = 'CTER' 
                if(name == 'OT1'):
                    self.q_total -= 0.17
                    self.qx_total -= coords[0] * 0.17
                    self.qy_total -= coords[1] * 0.17
                    self.qz_total -= coords[2] * 0.17
            elif name in ['CAY', 'HY1', 'HY2', 'HY3','CY', 'OY']:
                res = 'ACE'
            else:
                res = atom.get_parent().resname
                
            charge = self.charmm.AA[res].atom_chrg[name]
            coords = atom.coord.tolist()
    
            self.q_total += charge
            self.qx_total += coords[0] * charge
            self.qy_total += coords[1] * charge
            self.qz_total += coords[2] * charge
        qr_total = self.qr_total()       
        return [qr_total[0], qr_total[1], qr_total[2]]

    def qr_total(self):
        return [self.qx_total, self.qy_total, self.qz_total]

    def total_charge(self):
        return self.q_total

""" 
    def add(self, atom):
        #rtp_residue = self.rtp.residues.get(res.resname)
        #if not rtp_residue:
        #    rtp_residue = _do_histidine_magic(res, self.rtp)
        #
        #rtp_atom = rtp_residue.atoms.get(name)
        #if not rtp_atom:
        #    rtp_atom = _do_atom_magic(rtp_residue, atom)
        name = atom.get_name()
        if name in ['OT1','OT2']:
            res = 'CTER' 
            if(name == 'OT1'):
        
        else:
            res = atom.get_parent()
            
        charge = self.charmm.AA[res.resname].atom_chrg[name]
        coords = atom.coord.tolist()

        self.q_total += charge
        self.qx_total += coords[0] * charge
        self.qy_total += coords[1] * charge
        self.qz_total += coords[2] * charge
"""