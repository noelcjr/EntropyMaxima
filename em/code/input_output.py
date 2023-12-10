# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 23:26:00 2016
@author: noel
"""
import sys
import os
import pandas as pd
import Bio.PDB as struct
import em.code.utilities as utilities
import em.code.Molecular_Rigid_Manipulations as MRM
#from em.charmm.gen import *
#from em.charmm import run

"""
def CIF_summary(path):
    raise NotImplementedError("CIF summary needs to be coded.")

# Deprecated for __init__ in structure.py 
def _read_structure(path, pdb_id='pdb', cif_id='cif' ):
    file_name = os.path.basename(path).split('.')[0]
    file_sufix = os.path.basename(path).split('.')[1]
    dir_path = os.path.dirname(path)
    if file_sufix == 'pdb':
        parser = struct.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, path)
    elif file_sufix == 'cif':
        parser = struct.MMCIFParser()
        structure = parser.get_structure(cif_id, path)
    else:
        print("ERROR: Unreognized file type " + file_sufix + " in " + file_name)
        sys.exit(1)
    return structure, dir_path, file_name

def _save_structure(structure, path):
    io = struct.PDBIO()
    io.set_structure(structure)
    io.save(path)

Moved to structure.py
def PDBs_from_CIF(inFile, optionsmodels, optionschains, optionsgroups):
    structure, dir_path, file_name = _read_structure(inFile)
    if optionsmodels:
        for i in structure.get_models():
            _save_structure(i, dir_path + file_name + "_" + str(i.get_id()) + ".pdb")
    if optionschains:
        if not optionsgroups:
            print("ERROR: When outputing PDBs by chain identifier, and option groups must be included.")
            sys.exit(1)
        a = optionsgroups
        b = a.split(',')
        for h in b:
            for i in structure.get_models():
                s2 = struct.Structure.Structure(h)
                m2 = struct.Model.Model(i.get_id())
                s2.add(m2)
                for j in i.get_chains():
                    if j.get_id() in h:
                        m2.add(j)
                _save_structure(s2, dir_path + file_name + "_" + str(i.get_id()) + "_" + h + ".pdb")

def align_pdbs(referencePath, fitPath, optionsrefatomsA, optionsfitatomsA, optionsoutA,optionsaddatomsA=""):
    if not os.path.exists(referencePath):
        print "Error: File path for reference PDB or CIF file does not exist."
        print("Type -h or --help for description and options.")
        sys.exit(1)
    if not os.path.exists(fitPath):
        print "Error: File path for PDB or CIF file to fit to reference does not exist."
        print("Type -h or --help for description and options.")
        sys.exit(1)
    ref_al = []
    for i in optionsrefatomsA.split(':'):
        ref_al.append(tuple(i.split(',')))
    fit_al = []
    for i in optionsfitatomsA.split(':'):
        fit_al.append(tuple(i.split(',')))

    ref_structure, dir_path, _ = _read_structure(referencePath, 'reference', 'reference')
    fit_structure, dir_path, _ = _read_structure(fitPath, 'fit', 'fit')

    # Use the first model in the pdb-files for alignment
    # Change the number 0 if you want to align to another structure
    ref_model = ref_structure[0]
    fit_model = fit_structure[0]
    # Make a list of the atoms (in the structures) you wish to align.
    # In this case we use CA atoms whose index is in the specified range
    ref_atoms = []
    fit_atoms = []
    # Iterate of all chains in the model in order to find all residues
    for ref_chain in ref_model:
        for r_a in ref_al:
            if ref_chain.get_id() == r_a[1]:
                for ref_res in ref_chain:
                    # Check if residue number ( .get_id() ) is in the list
                    if ref_res.get_id()[1] in range(int(r_a[2]), int(r_a[3]) + 1):
                        # Append CA atom to list
                        ref_atoms.append(ref_res[r_a[0]])
    # Do the same for the sample structure
    for fit_chain in fit_model:
        for r_a in fit_al:
            if fit_chain.get_id() == r_a[1]:
                for fit_res in fit_chain:
                    if fit_res.get_id()[1] in range(int(r_a[2]), int(r_a[3]) + 1):
                        fit_atoms.append(fit_res[r_a[0]])
    # Now we initiate the superimposer:
    super_imposer = struct.Superimposer()
    super_imposer.set_atoms(ref_atoms, fit_atoms)
    super_imposer.apply(fit_model.get_atoms())
    if optionsaddatomsA != "":
        # A region is defined by a chain, an initial residue number, a followup residue number
        # to start of region, and another residue number fotr the end of the region B,3,10
        # both begining and end are included.
        add_region = []
        for i in optionsaddatomsA.split(':'):
            add_region.append(tuple(i.split(',')))
        if len(add_region) != 2:
            print("ERROR: Only two entries in the addatom option. One for reference, and one for fit.")
            sys.exit(1)
        print("Adding residues from (" + referencePath + "," + add_region[0][0] + "," + add_region[0][1] + "," +
              add_region[0][2] + ") ")
        print(
            " to (" + fitPath + "," + add_region[1][0] + "," + add_region[1][1] + "," + add_region[1][2] + ")")
        add_res = []
        # Add residues before missing segment from incomplete chain (fit)
        for i in fit_model:
            if i.get_id() == add_region[1][0]:
                for j in i.get_residues():
                    if j.get_id()[1] < int(add_region[0][1]):
                        add_res.append(j)
        # add residues in missing segment from complete chain (reference)
        for i in ref_model:
            if i.get_id() == add_region[0][0]:
                for j in i.get_residues():
                    if j.get_id()[1] in range(int(add_region[0][1]), int(add_region[0][2]) + 1):
                        add_res.append(j)
        # add residues after missing segment.from incomplete (fit)
        for i in fit_model:
            if i.get_id() == add_region[1][0]:
                for j in i.get_residues():
                    if j.get_id()[1] > int(add_region[0][2]):
                        add_res.append(j)
        newChain = struct.Chain.Chain(add_region[1][0])
        for i in add_res:
            newChain.add(i)
        # put chains in a list
        chain_order = []
        for i in fit_model:
            chain_order.append(i)
        # delete all chains from the fit model
        for i in chain_order:
            fit_model.detach_child(i.get_id())
        # Add chains back in fit_model making sure that newChain replaces the incoplete chain
        for i in chain_order:
            if i.get_id() == add_region[1][0]:
                fit_model.add(newChain)
            else:
                fit_model.add(i)
                ###########   TEST  ###############
                # for i in fit_model:
                #    print(i,i.get_id())
                #    for j in i.get_residues():
                #        print(j.get_full_id(),j.get_resname())
                # sys.exit(1)
                # TODO it is possible that a section with different resids is added to a fitted section. In this case residd
                #      fit the chains that were just added need to be renumbered approprietly.
                # Last check to make sure residue numbers in completed chain are monotonically increased.
                # if add_res_range_from_to[0][1:2] == add_res_range_from_to[1][1:2]:
                #    pass
                # else:
                #    get_fit_first_res_id = True
                #    fit_first_res_id = -1
                #    for i in fit_model:
                #        if i.get_id() == add_res_range_from_to[1][0]:
                #            for j in i.get_residues():
                #                if get_fit_first_res_id:
                #                    get_fit_first_res_id = False
                #                    fit_first_res_id = j.get_id()[1]
                #    get_ref_first_res_id = True
                #    ref_first_res_id = -1
                #    for i in ref_model:
                #        if i.get_id() == add_res_range_from_to[0][0]:
                #            for j in i.get_residues():
                #                if get_ref_first_res_id:
                #                    get_ref_first_res_id = False
                #                    ref_first_res_id = j.get_id()[1]
    # Print RMSD:
    print("Fiting " + fitPath + " by " + optionsrefatomsA + " to " + referencePath + " by " + optionsfitatomsA + "\
 RMSD=" + str(super_imposer.rms))
    # Save the aligned version
    _save_structure(fit_model, dir_path + optionsoutA)

def write_pdb_from_crd(lines):
    line = '{:6}{:>5} {:4}{:1}{:3} {:1}{:>4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:2}'
    pdb_converted_lines = []
    for i in lines:
        t = i.split()
        pdb_converted_lines.append(line.format('ATOM',t[0],t[3],'',t[2],t[7],int(t[1]),'',float(t[4]),float(t[5]),\
                                               float(t[6]),float('0'),float('0'),t[3][0],''))
    return pdb_converted_lines
"""
def write_pdb(SS, basedir, filename, output='original'):
    ''' TODO: write pdbdirectly with all chains the right way.
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    -------------------------------------------------------------------------------------
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 20        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.
    line = '{:6}{:>5} {:4}{:1}{:3} {:1}{:>4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:2}'
    '''
    line = '{:6}{:>5} {:4}{:1}{:3} {:1}{:>4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:2}'
    # TODO: label_aaid shoould move to a utilities class
    # TODO: check that label id still works for entities that have multiple chains
    # TODO: I changed back to pronting the chain directly without labell_aaid. Check this is right.
    # label_aaid = {1:'A',2:'B',3:'C',4:'D'}
    for g in SS.models:
        lines = []
        count = 1
        # self.get_missing_aa_DataFrame(g)
        # origianl atomic structures
        if output == 'original':
            for h in range(SS.Full_Structure.shape[0]):
                if SS.Full_Structure.loc[h, 'atminfo' + g]:
                    lines.append(line.format('ATOM', str(count), SS.Full_Structure.loc[h, 'atmtyp1'], \
                                             '', SS.Full_Structure.loc[h, 'aa'], SS.Full_Structure.loc[h, 'chain'], \
                                             int(SS.Full_Structure.loc[h, 'aaid']), '',
                                             SS.Full_Structure.loc[h, 'x' + g], \
                                             SS.Full_Structure.loc[h, 'y' + g], SS.Full_Structure.loc[h, 'z' + g], \
                                             SS.Full_Structure.loc[h, 'occupancy'],
                                             SS.Full_Structure.loc[h, 'B_factor'], \
                                             SS.Full_Structure.loc[h, 'atmtyp3'], ''))
                    count += 1
        elif output == 'minimum':
            '''Only heavy atoms for missing amino acids plus all atoms already
            in the x-structure'''
            for h in SS.Full_Structure.index:
                if SS.Full_Structure.loc[h, 'aainfo' + g]:
                    if SS.Full_Structure.loc[h, 'atminfo' + g]:
                        lines.append(line.format('ATOM', str(count), SS.Full_Structure.loc[h, 'atmtyp1'], \
                                                 '', SS.Full_Structure.loc[h, 'aa'],
                                                 SS.Full_Structure.loc[h, 'chain'], \
                                                 int(SS.Full_Structure.loc[h, 'aaid']), '',
                                                 SS.Full_Structure.loc[h, 'x' + g], \
                                                 SS.Full_Structure.loc[h, 'y' + g],
                                                 SS.Full_Structure.loc[h, 'z' + g], \
                                                 SS.Full_Structure.loc[h, 'occupancy'],
                                                 SS.Full_Structure.loc[h, 'B_factor'], \
                                                 SS.Full_Structure.loc[h, 'atmtyp3'], ''))
                        count += 1
                    else:
                        if SS.Full_Structure.loc[h, 'atmtyp3'] != 'H':
                            lines.append(line.format('ATOM', str(count), SS.Full_Structure.loc[h, 'atmtyp1'], \
                                                     '', SS.Full_Structure.loc[h, 'aa'],
                                                     SS.Full_Structure.loc[h, 'chain'], \
                                                     int(SS.Full_Structure.loc[h, 'aaid']), '',
                                                     SS.Full_Structure.loc[h, 'x' + g], \
                                                     SS.Full_Structure.loc[h, 'y' + g],
                                                     SS.Full_Structure.loc[h, 'z' + g], \
                                                     SS.Full_Structure.loc[h, 'occupancy'],
                                                     SS.Full_Structure.loc[h, 'B_factor'], \
                                                     SS.Full_Structure.loc[h, 'atmtyp3'], ''))
                            count += 1
                else:
                    if SS.Full_Structure.loc[h, 'atmtyp3'] != 'H':
                        lines.append(line.format('ATOM', str(count), SS.Full_Structure.loc[h, 'atmtyp1'], \
                                                 '', SS.Full_Structure.loc[h, 'aa'],
                                                 SS.Full_Structure.loc[h, 'chain'], \
                                                 int(SS.Full_Structure.loc[h, 'aaid']), '',
                                                 SS.Full_Structure.loc[h, 'x' + g], \
                                                 SS.Full_Structure.loc[h, 'y' + g],
                                                 SS.Full_Structure.loc[h, 'z' + g], \
                                                 SS.Full_Structure.loc[h, 'occupancy'],
                                                 SS.Full_Structure.loc[h, 'B_factor'], \
                                                 SS.Full_Structure.loc[h, 'atmtyp3'], ''))
                        count += 1
        elif output == 'all':
            for h in SS.Full_Structure.index:
                x = pd.notnull(SS.Full_Structure.loc[h, 'x' + g])
                y = pd.notnull(SS.Full_Structure.loc[h, 'y' + g])
                z = pd.notnull(SS.Full_Structure.loc[h, 'z' + g])
                if x & y & z:
                    lines.append(line.format('ATOM', str(count), SS.Full_Structure.loc[h, 'atmtyp1'], \
                                             '', SS.Full_Structure.loc[h, 'aa'], SS.Full_Structure.loc[h, 'chain'], \
                                             int(SS.Full_Structure.loc[h, 'aaid']), '',
                                             SS.Full_Structure.loc[h, 'x' + g], \
                                             SS.Full_Structure.loc[h, 'y' + g], SS.Full_Structure.loc[h, 'z' + g], \
                                             SS.Full_Structure.loc[h, 'occupancy'],
                                             SS.Full_Structure.loc[h, 'B_factor'], \
                                             SS.Full_Structure.loc[h, 'atmtyp3'], ''))
                    count += 1
        else:
            print("ERROR: output type is correct. Choose only 'original','minimum' and 'all'.")
            print("       program will exit.")
            sys.exit(1)
        if basedir == '':
            basedir = '.'
        if basedir[-1] != "/":
            basedir += "/"
        outFile = open(basedir + filename + '_' + g + '.pdb', 'w')
        for i in lines:
            outFile.write(i + '\n')
        outFile.close()
"""
def fix_pdb_from_CHARMM(PDB_file, a=21, b=72):
    ''' This fixes the problem of CHARMM not generating a chain id.
    But needs to be tested in more systems '''
    if not os.path.exists(PDB_file):
        print "Error: No --pdbin option entered. Type -h or --help."
        sys.exit(1)
    inFile = open(PDB_file, 'r')
    contents = inFile.read()
    contents_list = contents.split('\n')
    for i in range(len(contents_list)):
        if contents_list[i][0:4] == 'ATOM':
            temp = list(contents_list[i])
            # The following two indexes are off by one becasue list index start at 0, and lines in files at 1.
            temp[a] = temp[b]
            contents_list[i] = ''.join(temp)
    outFile = open(PDB_file, 'w')
    for i in contents_list:
        outFile.write(i + '\n')
    outFile.close()


def _print_template(outFile,template):
    for i in template:
        outFile.write(i)

def prepare_pdb_for_charmm(optionspdbin,prot_ends, reduceopt="-HIS -FLIP -OH -ROTEXOH -BUILD -OCC0.0 -H2OOCC0.0 -H2OB1000"):
    if not os.path.exists(optionspdbin):
        print "Error: No path --input option."
        sys.exit(1)
    input_name = os.path.basename(optionspdbin).split('.')[0]
    output_name = input_name.lower()
    os.system("reduce {} {} > {}".format(reduceopt,optionspdbin,output_name+"r.pdb"))
    if not os.path.exists(output_name+"r.pdb"):
        print "Error: System call to run reduce did not produce an output file."
        sys.exit(1)
    pdb_parser = struct.PDBParser(QUIET=True)
    sp = pdb_parser.get_structure(input_name, optionspdbin)
    # This fixes the problem of CHARMM not generating a chain id.
    # But needs to be tested in more systems
    print("For pdb=", optionspdbin, input_name)
    # TODO the problem with c_10 not being process seems to be here.
    for i in sp.get_chains():
        #print("chain fix", i, i.id, i.get_level())
        if i.id == ' ':
            i.id = i.get_level()
    chain_resnum = []
    for i in sp.get_chains():
        min_max = [1000000, -1000000]
        for j in i.get_residues():
            if j.get_full_id()[3][1] < min_max[0]:
                min_max[0] = j.get_full_id()[3][1]
            if j.get_full_id()[3][1] > min_max[1]:
                min_max[1] = j.get_full_id()[3][1]
        chain_resnum.append(min_max)
    print("chain_resnum=", chain_resnum)
    # TODO: the following 2 loops also found in flower.py. They need to be moved to a class to avoid repetitive code
    ids = {}
    for i in string.ascii_uppercase:
        ids[i] = False
    for i in sp[0]:
        ids[i.id] = True
    # TODO: FOR some reason ids are not lased consecutively but functions works
    # def get_next_avail_id():
    #    id_val = -1
    #    for i in ids:
    #        if not ids[i]:
    #            ids[i] = True
    #            id_val = i
    #            break
    #    return id_val
    num_residues = {}
    count_ions = 0
    HC = MRM.Histidin_Correction()
    for i in sp.get_models():
        for j in i.get_chains():
            count = 0
            for k in j.get_residues():
                count += 1
                if k.resname == 'HIS':
                    newname = HC.rename_HIS(k)
                    print('change HIS', k.get_id()[1], 'to', newname)
                    if not newname == 'HIS':
                        for l in k.get_atom():
                            k.resname = newname
                elif k.resname == ' NA':
                    k.resname = 'SOD'
                    count_ions += 1
                    # k.id = (' ',count_ions,' ')
                    for l in k.get_atom():
                        l.id = 'SOD'
                        l.fullname = ' SOD'
                elif k.resname == '  A':
                    k.resname = 'ADE'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == ' DA':
                    k.resname = 'ADE'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == '  G':
                    k.resname = 'GUA'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == ' DG':
                    k.resname = 'GUA'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == '  C':
                    k.resname = 'CYT'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == ' DC':
                    k.resname = 'CYT'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == '  T':
                    k.resname = 'THY'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == ' DT':
                    k.resname = 'THY'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == '  U':
                    k.resname = 'URA'
                    for l in k.get_atom():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
            num_residues[j.id] = count
    #This fixes when there are multiple chains and the residue number gets
    #reset at the beginning of each chain. With this fix, residue numbers
    #will be renumbered
    # f = '/home/noel/Projects/Protein_design/Insulin/Struct_Prep/pdbs/init_setup_bridge/cen/2hiu_1.pdb'
    # cen = pdb_parser.get_structure('test', f)
    line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
    atom_count = 0
    resi_count = 0
    lines = []
    for a in sp.get_models():
        if a.id == 0:
            for b in a.get_chains():
                for c in b.get_residues():
                    resi_count += 1
                    for d in c.get_atom():
                        atom_count += 1
                        lines.append(line.format(atom_count, resi_count, \
                                                 c.get_resname(), d.get_name(), d.get_coord()[0], \
                                                 d.get_coord()[1], d.get_coord()[2], b.id, \
                                                 c.get_full_id()[3][1], 0))

    outFile = open(output_name+".crd", 'w')
    outFile.write('* CRD generated with fixed HS \n')
    outFile.write('* Name: NN\n')
    outFile.write('* Title: New Molecule Generated by pdb_preparation.py\n')
    outFile.write('*\n')
    outFile.write('{:>5}'.format(len(lines)) + '\n')
    for i in lines:
        outFile.write(i + '\n')
    outFile.close()

    # The Following loops work for fixers that do not require an offset
    # Code needs to be checked for other possible ways to do this renaming
    #a_pair = {
    #    'id': 'A',
    #    'seq': "A.SEQ",
    #    'first': 'none',
    #    'last': 'none',
    #    'inp': "A_FIXRES.INP"
    #}
    # TODO we need a way to check that prot_ends match number of chains in structure.
    chain_terminals = prot_ends.split(':')
    temp = {}
    for i in chain_terminals:
        temp2 = i.split(',')
        temp[temp2[0]] = {}
        temp[temp2[0]]['first'] = temp2[1]
        temp[temp2[0]]['last'] = temp2[2]
    chains_info = {}
    for i in sp.get_chains():
        file_name = i.id + '.SEQ'
        chains_info[i.id] = {}
        chains_info[i.id]['id'] = i.id
        chains_info[i.id]['seq'] = file_name
        chains_info[i.id]['first'] = temp[i.id]['first']
        chains_info[i.id]['last'] = temp[i.id]['last']
        outFile = open('./'+file_name, 'w')
        outFile.write('* Sequence of Chain ' + i.id + '\n')
        outFile.write('* Name: Molecule\n')
        outFile.write('* Title: New Molecule Generated by pdb_tp_crd.py\n')
        outFile.write('*\n')
        x = 5 - len(str(num_residues[i.id]))
        outFile.write(' ' * x + str(num_residues[i.id]) + '\n')
        for j in i.get_residues():
            outFile.write(j.resname.strip(' ') + '\n')
        outFile.close()
    offset = 0
    residue_count = 1
    for i in sp.get_chains():
        file_name = i.id + '_FIXRES.INP'
        chains_info[i.id]['inp'] = file_name
        outFile = open('./'+file_name, 'w')
        outFile.write("! CHARMM Script to renumber RESID's after reading from SEQ file\n")
        outFile.write('! Generated by pdb_to_crd.py\n')
        outFile.write('\n')
        outFile.write('! First pass .... assigning temporary RESID\n')
        count = 1
        for j in sp.get_residues():
            if j.get_parent().id == i.id:
                if count == 1:
                    offset = j.get_id()[1] - 1
                outFile.write('! ' + str(offset + count) + ' , ' + str(count) + '\n')
                outFile.write('rename resid x' + str(count) + ' sele segid ' + i.id + ' .and. resid ' + str(
                    count) + ' end\n')
                count += 1
                residue_count += 1
        outFile.write('! First pass complete.\n')
        outFile.write('\n')
        outFile.write('! Second pass .... re-assigning correct RESID\n')
        count = 1
        for j in sp.get_residues():
            if j.get_parent().id == i.id:
                outFile.write(
                    'rename resid ' + str(offset + count) + ' sele segid ' + i.id + ' .and. resid x' + str(
                        count) + ' end\n')
                count += 1
        outFile.write('! Second pass complete.\n')
        outFile.write('\n')
        outFile.write('! End of generated script\n')
        outFile.close()
    outFile = open("setup_"+input_name+".inp", 'w')
    _print_template(outFile,gen_parameters().split('/n'))
    for i in sp.get_chains():
        _print_template(outFile,fix_sequence(chains_info[i.id]))
    _print_template(outFile,wmain('charge'))
    _print_template(outFile,out_psf(input_name.lower()+'r'))
    _print_template(outFile,out_ic(input_name.lower()+'r'))
    _print_template(outFile,in_crd(input_name.lower()))
    _print_template(outFile,build_heavy_atoms())
    _print_template(outFile,build_hydrogens())
    _print_template(outFile,wmain('charge'))
    _print_template(outFile,out_crd(output_name+'r'))
    _print_template(outFile,out_pdb(output_name+'r'))
    _print_template(outFile,out_psf_xplor(output_name+'r'))
    _print_template(outFile,stop())
    outFile.close()
    run.run("setup_"+input_name+".inp","setup_"+input_name+".out")
    # TODO WARNING: Check output if Ions are present in the original pdb structure.
    #               Check INF and SEQ files for systems with missing amino acids.')
    #               Modifications to the code might be required.')

def _get_atom_coordinates(sp):
    atoms = []
    for i in sp.get_models():
        for j in i.get_chains():
            for k in j.get_atoms():
                atoms.append(k.get_coord())
    return atoms


def min_max(path):
    sp, _, _2 = _read_structure(path)
    atom_coordinates = _get_atom_coordinates(sp)
    x,y,z = zip(*atom_coordinates)

    Xmin,Xmax = min(x), max(x)
    Ymin,Ymax = min(y), max(y)
    Zmin,Zmax = min(z), max(z)

    print('Results from ', path)
    print('Number of Atoms ', len(atom_coordinates))
    print('Xmin - Xmax', Xmin, ' - ', Xmax)
    print('Ymin - Ymax', Ymin, ' - ', Ymax)
    print('Zmin - Zmax', Zmin, ' - ', Zmax)
    print('--------------------------------')
    print('Celbasis X, Y, Z:', (Xmax - Xmin), (Ymax - Ymin), (Zmax - Zmin))

class crd(object):
    '''
     This fixes when there are multiple chains and the residue number gets
     reset at the beginning of each chain. With this fix, residue numbers
     will be renumbered
      1 - 5 Integer Atom no. sequential
     6 - 10 Integer ires Residue position from file beginning
    11 - 11 Space
    12 - 15 Achar resname Residue name
    16 - 16 Space
    17 - 20 Achar type Atom type, IUPAC name of atom left justified
    21 - 30 Real(10.5) x Orthogonal coordinates for X, Angstroms
    31 - 40 Real(10.5) y Orthogonal coordinates for Y, Angstroms
    41 - 50 Real(10.5) z Orthogonal coordinates for Z, Angstroms
    51 - 51 Space
    52 - 55 Achar segid Segment identifier
    56 - 56 Space
    57 - 60 Achar resid Residue number within the segment
    61 - 70 Real(10.5) Weighting array value
    line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
    '''

    # https://www.charmm.org/charmm/documentation/by-version/c40b1/params/doc/io/#IOFORM
    def __init__(self, crd_file):
        self.CRD_file = crd_file
        self.crd_lines = self.read_crd_file()

    def read_crd_file(self):
        inFile = open(self.CRD_file, 'r')
        contents = inFile.read()
        crd_contents = []
        #  the next loop could probably be simpler, but I made sure it can handle atoms with 
        atom_number = 0
        number_atoms = 0
        for i in contents.split('\n'):
            ii = i.strip()
            if ii.isdigit():
                number_atoms = int(ii)
            else:
                if atom_number >= number_atoms:
                    pass
                else:
                    if len(ii) > 0 and ii[0] != '*':
                        crd_contents.append(ii)
                        atom_number += 1
        return crd_contents

class psf(object):
    # http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win-html/node24.html
    def __init__(self, psf_file):
        self.PSF_file = psf_file
        self.psf_lines = self.read_psf_file()

    def read_psf_file(self):
        inFile = open(self.PSF_file, 'r')
        contents = inFile.read()
        psf_contents = []
        atom_number = 0
        number_atoms = 0
        is_natom = False
        for i in contents.split('\n'):
            ii = i.strip().split()
            if '!NATOM' in ii:
                number_atoms = int(ii[0])
                is_natom = True
                continue
            if is_natom:
                psf_contents.append(ii)
                atom_number += 1
                if atom_number == number_atoms:
                    is_natom = False
        return psf_contents
class crd(object):
    '''
     This fixes when there are multiple chains and the residue number gets
     reset at the beginning of each chain. With this fix, residue numbers
     will be renumbered
      1 - 5 Integer Atom no. sequential
     6 - 10 Integer ires Residue position from file beginning
    11 - 11 Space
    12 - 15 Achar resname Residue name
    16 - 16 Space
    17 - 20 Achar type Atom type, IUPAC name of atom left justified
    21 - 30 Real(10.5) x Orthogonal coordinates for X, Angstroms
    31 - 40 Real(10.5) y Orthogonal coordinates for Y, Angstroms
    41 - 50 Real(10.5) z Orthogonal coordinates for Z, Angstroms
    51 - 51 Space
    52 - 55 Achar segid Segment identifier
    56 - 56 Space
    57 - 60 Achar resid Residue number within the segment
    61 - 70 Real(10.5) Weighting array value
    line = '{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
    '''

    # https://www.charmm.org/charmm/documentation/by-version/c40b1/params/doc/io/#IOFORM
    def __init__(self, crd_file):
        self.CRD_file = crd_file
        self.crd_lines = self.read_crd_file()

    def read_crd_file(self):
        inFile = open(self.CRD_file, 'r')
        contents = inFile.read()
        crd_contents = []
        #  the next loop could probably be simpler, but I made sure it can handle atoms with 
        atom_number = 0
        number_atoms = 0
        for i in contents.split('\n'):
            ii = i.strip()
            if ii.isdigit():
                number_atoms = int(ii)
            else:
                if atom_number >= number_atoms:
                    pass
                else:
                    if len(ii) > 0 and ii[0] != '*':
                        crd_contents.append(ii)
                        atom_number += 1
        return crd_contents

class psf(object):
    # http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win-html/node24.html
    def __init__(self, psf_file):
        self.PSF_file = psf_file
        self.psf_lines = self.read_psf_file()

    def read_psf_file(self):
        inFile = open(self.PSF_file, 'r')
        contents = inFile.read()
        psf_contents = []
        atom_number = 0
        number_atoms = 0
        is_natom = False
        for i in contents.split('\n'):
            ii = i.strip().split()
            if '!NATOM' in ii:
                number_atoms = int(ii[0])
                is_natom = True
                continue
            if is_natom:
                psf_contents.append(ii)
                atom_number += 1
                if atom_number == number_atoms:
                    is_natom = False
        return psf_contents
"""
