# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 23:26:00 2016
@author: noel
"""
import sys
import os
import Bio.PDB as struct
#import em.code.utilities as utilities
import string
import em.code.Molecular_Rigid_Manipulations as MRM
from em.charmm.gen import *
from . import *
from em.charmm import run
# import re

def _print_template(outFile,template):
    for i in template:
        outFile.write(i)

def prepare_pdb_for_charmm(pdbin,
                           prot_ends,
                           dir_path = ''):
    if not os.path.exists(pdbin):
        print("Error: Invalid path in --structure option.")
        sys.exit(1)
    if (os.path.basename(pdbin).split('.')[1] != 'pdb'):
        print("""ERROR: A PDB file with lower case 'pdb' suffix is required.
              Program will exit without results.""")
        sys.exit(1)
    input_name = os.path.basename(pdbin).split('.')[0]
    if dir_path == '':
        output_name = input_name.lower()
    else:
        if dir_path[-1] == '/':
            pass
        else:
            dir_path += '/'
        output_name = dir_path+input_name.lower()
    reduceopt="-HIS -FLIP -OH -ROTEXOH -BUILD -OCC0.0 -H2OOCC0.0 -H2OB1000"
    os.system("reduce {} {} > {} &".format(reduceopt,pdbin,"r.pdb"))
    if not os.path.exists(output_name+"r.pdb"):
        print("Error: System call to run reduce didn't output a PDB file.")
        sys.exit(1)
    pdb_parser = struct.PDBParser(QUIET=True)
    sp = pdb_parser.get_structure(input_name, pdbin)
    for i in self.structure.get_models():
        count_models += 1
    if count_models != 1:
        print("ERROR: This function works on a PDB file with only one model.")
        print("       Extract PDBs from CIF if needed. Simulations can only be done")
        print("       from a single structure and not multiple structure models.")
        sys.exit(1)
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
                    print('changed HIS', k.get_id()[1], 'to', newname)
                    if not newname == 'HIS':
                        for l in k.get_atoms():
                            k.resname = newname
                elif k.resname == ' NA':
                    k.resname = 'SOD'
                    count_ions += 1
                    # k.id = (' ',count_ions,' ')
                    for l in k.get_atoms():
                        l.id = 'SOD'
                        l.fullname = ' SOD'
                elif k.resname == '  A':
                    k.resname = 'ADE'
                    for l in k.get_atoms():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == ' DA':
                    k.resname = 'ADE'
                    for l in k.get_atoms():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == '  G':
                    k.resname = 'GUA'
                    for l in k.get_atoms():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == ' DG':
                    k.resname = 'GUA'
                    for l in k.get_atoms():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == '  C':
                    k.resname = 'CYT'
                    for l in k.get_atoms():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == ' DC':
                    k.resname = 'CYT'
                    for l in k.get_atoms():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == '  T':
                    k.resname = 'THY'
                    for l in k.get_atoms():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == ' DT':
                    k.resname = 'THY'
                    for l in k.get_atoms():
                        if l.id == 'OP1':
                            l.id = 'O1P'
                            l.fullname = ' O1P'
                        elif l.id == 'OP2':
                            l.id = 'O2P'
                            l.fullname = ' O2P'
                elif k.resname == '  U':
                    k.resname = 'URA'
                    for l in k.get_atoms():
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
    # f ='/home/noel/Projects/Protein_design/Insulin/struct_prep/2hiu/2hiu.pdb'
    # cen = pdb_parser.get_structure('test', f)
    line='{:>5}{:>5}{:>4}  {:4}{:>10.5f}{:>10.5f}{:>10.5f} {:4} {:<4}{:>12.5f}'
    atom_count = 0
    resi_count = 0
    lines = []
    for a in sp.get_models():
        if a.id == 0:
            for b in a.get_chains():
                for c in b.get_residues():
                    resi_count += 1
                    for d in c.get_atoms():
                        atom_count += 1
                        lines.append(line.format(atom_count, resi_count, \
                                                 c.get_resname(), \
                                                 d.get_name(), \
                                                 d.get_coord()[0], \
                                                 d.get_coord()[1], \
                                                 d.get_coord()[2], b.id, \
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

    # TODO we need way to check prot_ends match number of chains in structure
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
        outFile = open(dir_path+file_name, 'w')
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
        outFile = open(dir_path+file_name, 'w')
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
    print(dir_path+"setup_"+input_name+".inp")
    outFile = open(dir_path+"setup_"+input_name+".inp", 'w')
    _print_template(outFile,gen_parameters().split('/n'))
    for i in sp.get_chains():
        _print_template(outFile,fix_sequence(chains_info[i.id]))
    _print_template(outFile,wmain('charge'))
    _print_template(outFile,write_psf(input_name.lower()+'r'))
    _print_template(outFile,write_ic(input_name.lower()+'r'))
    _print_template(outFile,read_crd(input_name.lower()))
    _print_template(outFile,build_heavy_atoms())
    _print_template(outFile,build_hydrogens())
    _print_template(outFile,wmain('charge'))
    _print_template(outFile,write_crd(output_name+'r'))
    _print_template(outFile,write_pdb(output_name+'r'))
    _print_template(outFile,write_psf_xplor(output_name+'r'))
    _print_template(outFile,stop())
    outFile.close()
    run.run(dir_path+"setup_"+input_name+".inp",dir_path+"setup_"+input_name+".out")
    
def simple_minimization():
    pass

def minimization_1(input_file,input_info):
    input_name = os.path.basename(input_file).split('.')[0]
    dir_name = os.path.dirname(input_file)
    output_name = input_name.lower()
    outFile = open("min1_"+input_name+".inp", 'w')
    _print_template(outFile,gen_parameters().split('/n'))
    _print_template(outFile,stream_nina_radii())
    if len(dir_name) == 0:
        _print_template(outFile,read_psf(input_name))
        _print_template(outFile,read_crd2(input_name))
    else:
        _print_template(outFile,read_psf(dir_name+"/"+input_name))
        _print_template(outFile,read_crd2(dir_name+"/"+input_name))
    input_list = input_info.split(':')
    input_list = [i.replace('-', ',') for i in input_list]
    input_list = [i.split(',') for i in input_list]
    _print_template(outFile,minimization1(output_name,input_list))
    _print_template(outFile,stop())
    outFile.close()
    run.run("min1_"+input_name+".inp","min1_"+input_name+".out")
