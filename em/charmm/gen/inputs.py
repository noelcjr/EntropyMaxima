# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 23:26:00 2016
@author: noel
"""
import sys
import os
import pandas as pd
import Bio.PDB as struct
#import em.describe.utilities as utilities
import string
import em.manipulate.Molecular_Rigid_Manipulations as MRM
from em.charmm.gen import *
from em.charmm import run

def fix_pdb_from_CHARMM(PDB_file, a=21, b=72):
    """ This fixes the problem of CHARMM not generating a chain id.
    But needs to be tested in more systems """
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
        print("chain fix", i, i.id, i.get_level())
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
    '''
     This fixes when there are multiple chains and the residue number gets
     reset at the beginning of each chain. With this fix, residue numbers
     will be renumbered
    '''
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
