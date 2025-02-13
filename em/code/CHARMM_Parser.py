# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 20:50:51 2016
@author: noel
"""
import re, pkg_resources
    
class amino_acid():  # To diferentiate from _Residue in RTPParser
    def __init__(self,typeRES):
        self.type_res = typeRES
        self.atom_type = {}
        self.atom_chrg = {}
        self.atoms = []
        self.bonds = []
        self.double = []
        self.impr = []
        self.cmap = []
        self.donor = []
        self.acceptor = []
        self.ic = []
        self.bild = []
        self.thet = []
        self.dihe = []
        self.delete_type = []
        self.delete_object = []
        self.current_groups = -1
    def add_group(self):
        self.current_groups += 1
        self.atoms.append([])
    def get_number_of_groups(self):
        return self.current_groups
    def add_atom(self,a):
        self.atoms[self.current_groups].append(a[0])
        self.atom_type[a[0]] = a[1]
        self.atom_chrg[a[0]] = float(a[2])
    def add_bond(self,b):
        self.bonds.append(b)
    def add_double(self,b):
        self.double.append(b)
    def add_impr(self,imp):
        self.impr.append(imp)
    def add_cmap(self,cm):
        self.cmap.append(cm)
    def add_donor(self,d):
        self.donor.append(d)
    def add_acceptor(self,ac):
        self.acceptor.append(ac)
    def add_ic(self,ic):
        self.ic.append(ic)
    def add_bild(self,bild):
        self.bild.append(bild)
    def add_thet(self,bild):
        self.thet.append(bild)
    def add_dihe(self,bild):
        self.dihe.append(bild)        
    def add_delete(self,a):
        self.delete_type.append(a[0])
        self.delete_object.append(a[1])
    
class atom_mass():
    def __init__(self):
        self.MASS = {}
        
    def add_atom(self, ln):
        s = ln.split()
        d = ln[ln.find('!')+2:-1]
        if s[2] in self.MASS:
            print("Atom already in dictionary")
        else:
            self.MASS[s[2]] = [s[1],s[3],s[4],d]

class read_charmm_FF():
    def __init__(self, verbose=False):
        self.toppath = pkg_resources.resource_filename('em', 'params/charmm27.ff/' + 'top_all27_prot_na.rtf')
        self.parpath = pkg_resources.resource_filename('em', 'params/charmm27.ff/' + 'par_all27_prot_na.prm')
        self.v = verbose
        # FIXME: Do not load DNA because of formating differences with amino acids.
        self.exceptions = ['GUA','ADE','CYT','THY','URA','5MET','5PHO','5POM','3PHO','3POM','3PO3','DELB','CY35','LKNA']
        self.AA = {}
        self.am = atom_mass()
        self.current_aa = ''
        self.load_topology()
        self.BONDS = {}
        self.ANGLES = {}
        self.DIHEDRALS = {}
        self.IMPROPER = {}
        self.NONBONDED = {}
        self.load_parameters()
        
    def load_parameters(self):
        lineNum = 0
        bonds = angles = dihedrals = improper = cmap = nonbonded = False
        with open(self.parpath, 'r') as f:
            for line in f:
                s = line[0:line.find('!')].split()
                lineNum += 1
                # CHECKS that I am parsing from the right section of the file
                if len(s) >= 1:
                    if s[0] == "BONDS":
                        bonds = True
                    elif s[0] == "ANGLES":
                        bonds = False
                        angles = True
                    elif s[0] == "DIHEDRALS":
                        angles = False
                        dihedrals = True
                    elif s[0] == "IMPROPER":
                        dihedrals = False
                        improper = True
                    elif s[0] == "CMAP":
                        improper = False
                        cmap = True
                    elif s[0] == "NONBONDED":
                        cmap = False
                        nonbonded = True
                # THROWS a message when not parsing from a section that does
                # not have the expected number of columns
                if bonds:
                    if len(s) == 4:
                        if ' '.join(s[:2]) not in self.BONDS:
                            self.BONDS[' '.join(s[:2])] = s[2:]
                        else:
                            if self.v:
                                print("ERROR: BOND "+' '.join(s[:2])+" entry already in dictionary.")
                if angles:
                    if len(s) == 5:
                        if ' '.join(s[:3]) not in self.ANGLES:
                            self.ANGLES[' '.join(s[:3])] = s[3:]
                        else:
                            if self.v:
                                print("ERROR: ANGLES "+' '.join(s[:3])+" entry already in dictionary.")
                if dihedrals:
                    if len(s) == 7:
                        if ' '.join(s[:4]) not in self.DIHEDRALS:
                            self.DIHEDRALS[' '.join(s[:4])] = s[4:]
                        else:
                            if self.v:
                                print("ERROR: DIHEDRAL "+' '.join(s[:4])+" entry already in dictionary.")
                if improper:
                    if len(s) == 7:
                        if ' '.join(s[:4]) not in self.IMPROPER:
                            self.IMPROPER[' '.join(s[:4])] = s[4:]
                        else:
                            if self.v:
                                print("ERROR: IMPROPER "+' '.join(s[:4])+" entry already in dictionary.")
                if cmap:
                    pass
                if nonbonded:                    
                    if ' '.join(s[:1]) not in self.NONBONDED:
                        self.NONBONDED[' '.join(s[:1])] = s[1:]
                    else:
                        if self.v:
                            print("ERROR: NONBONDED "+' '.join(s[:1])+" entry already in dictionary.")
            if not self.v:
                print("WARNING: in CHARMM_parser(), load_parameters. Some parameters are missing. Read \n parmaters with verbose set to True for details.")

    def load_topology(self):
        nuc_add_first_group = {'GUA':True,'ADE':True,'CYT':True,'THY':True,'URA':True}
        with open(self.toppath, 'r') as f:
            for line in f:
                s = line[1:line.find('!')].split()
                if re.compile(r'^(MASS)').match(line):
                    self.am.add_atom(line)            
                if re.compile(r'^(RESI)').match(line):
                    self.current_aa = s[1]
                    if self.current_aa not in self.exceptions:
                        self.AA[self.current_aa] = amino_acid('RESI')
                    if self.current_aa in self.exceptions[0:5]:
                        self.AA[self.current_aa] = amino_acid('RESI')
                if re.compile(r'^(PRES)').match(line):
                    self.current_aa = s[1]
                    if self.current_aa not in self.exceptions:
                        self.AA[self.current_aa] = amino_acid('PRES')
                    if self.current_aa in self.exceptions[0:5]:
                        self.AA[self.current_aa] = amino_acid('PRES')
                if self.current_aa not in self.exceptions:
                    if re.compile(r'^(GROUP)').match(line):
                        self.AA[self.current_aa].add_group()
                    if re.compile(r'^(ATOM)').match(line):
                        self.AA[self.current_aa].add_atom([s[1],s[2],s[3]])
                    if re.compile(r'^(DELETE)').match(line):
                        self.AA[self.current_aa].add_delete([s[1],s[2]])
                    if re.compile(r'^(BOND)').match(line):
                        self.AA[self.current_aa].add_bond([[s[i],s[i+1]] for i in range(1,len(s),2)])
                    if re.compile(r'^(DOUBLE)').match(line):
                        self.AA[self.current_aa].add_double([[s[i],s[i+1]] for i in range(1,len(s),2)])
                    if re.compile(r'^(IMPR)').match(line):
                        self.AA[self.current_aa].add_impr([[s[i],s[i+1],s[i+2],s[i+3]] for i in range(1,len(s),4)])
                    if re.compile(r'^(CMAP)').match(line):
                        self.AA[self.current_aa].add_cmap([[s[i],s[i+1],s[i+2],s[i+3]] for i in range(1,len(s),4)])
                    if re.compile(r'^(DONO)').match(line):
                        if len(s) == 2:
                            self.AA[self.current_aa].add_donor([[s[i],'&'] for i in range(1,len(s),1)])
                        elif len(s) == 3:
                            self.AA[self.current_aa].add_donor([[s[i],s[i+1]] for i in range(1,len(s),2)])
                        else:
                            print("ERROR. A DONOR line with more than two atoms is not allowed.")
                    if re.compile(r'^(ACCEPTOR)').match(line):
                        if len(s) == 2:
                            self.AA[self.current_aa].add_acceptor([[s[i],'&'] for i in range(1,len(s),1)])
                        elif len(s) == 3:
                            self.AA[self.current_aa].add_acceptor([[s[i],s[i+1]] for i in range(1,len(s),2)])
                        else:
                            print("ERROR. An ACCEPTOR line with more than two atoms is not allowed.")
                    if re.compile(r'^(IC)').match(line):
                        self.AA[self.current_aa].add_ic([s[i] for i in range(1,len(s))])
                    if re.compile(r'^(THET)').match(line):
                        self.AA[self.current_aa].add_thet([s[i] for i in range(1,len(s))])
                    if re.compile(r'^(BILD)').match(line):
                        self.AA[self.current_aa].add_bild([s[i] for i in range(1,len(s))])
                    if re.compile(r'^(DIHE)').match(line):
                        self.AA[self.current_aa].add_dihe([s[i] for i in range(1,len(s))])
                if self.current_aa in self.exceptions[0:5]:
                    if nuc_add_first_group[self.current_aa]:
                        self.AA[self.current_aa].add_group()
                        nuc_add_first_group[self.current_aa] = False
                    else:
                        if re.compile(r'^(GROUP)').match(line):
                            self.AA[self.current_aa].add_group()
                        if re.compile(r'^(ATOM)').match(line):
                            self.AA[self.current_aa].add_atom([s[1],s[2],s[3]])
                        if re.compile(r'^(DELETE)').match(line):
                            self.AA[self.current_aa].add_delete([s[1],s[2]])
                        if re.compile(r'^(BOND)').match(line):
                            self.AA[self.current_aa].add_bond([[s[i],s[i+1]] for i in range(1,len(s),2)])
                        if re.compile(r'^(DOUBLE)').match(line):
                            self.AA[self.current_aa].add_double([[s[i],s[i+1]] for i in range(1,len(s),2)])
                        if re.compile(r'^(IMPR)').match(line):
                            self.AA[self.current_aa].add_impr([[s[i],s[i+1],s[i+2],s[i+3]] for i in range(1,len(s),4)])
                        if re.compile(r'^(CMAP)').match(line):
                            self.AA[self.current_aa].add_cmap([[s[i],s[i+1],s[i+2],s[i+3]] for i in range(1,len(s),4)])
                        if re.compile(r'^(DONO)').match(line):
                            if len(s) == 2:
                                self.AA[self.current_aa].add_donor([[s[i],'&'] for i in range(1,len(s),1)])
                            elif len(s) == 3:
                                self.AA[self.current_aa].add_donor([[s[i],s[i+1]] for i in range(1,len(s),2)])
                            else:
                                print("ERROR. A DONOR line with more than two atoms is not allowed.")
                        if re.compile(r'^(ACCE)').match(line):
                            if len(s) == 2:
                                self.AA[self.current_aa].add_acceptor([[s[i],'&'] for i in range(1,len(s),1)])
                            elif len(s) == 3:
                                self.AA[self.current_aa].add_acceptor([[s[i],s[i+1]] for i in range(1,len(s),2)])
                            else:
                                print("ERROR. An ACCEPTOR line with more than two atoms is not allowed.")
                        if re.compile(r'^(IC)').match(line):
                            self.AA[self.current_aa].add_ic([s[i] for i in range(1,len(s))])
                        if re.compile(r'^(BILD)').match(line):
                            self.AA[self.current_aa].add_bild([s[i] for i in range(1,len(s))])
                #Unccoment for debuging
                #if self.current_aa in self.exceptions[5:9]:
                #    print(count,self.current_aa,line)
                #    count += 1
