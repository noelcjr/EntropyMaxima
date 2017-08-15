#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Tue May  2 19:44:35 2017

@author: noel
"""
import os
import sys
# The paths below allows spyder to find .pyc, and the program to find param 
# files, It is now set for my environment only but for distribution it needs
# to be setup for each user installation or Errors will occur.
folder = 'EntropyMaxima/'
path = '/home/noel/Projects/Protein_design/'+folder
sys_path = path+'/src'
if not os.path.isdir(sys_path):
    print("Error: The parameter directory "+sys_path+" is not found.")
    sys.exit(1)
param_path = path+'params/charmm27.ff/'
if not os.path.isdir(param_path):
    print("Error: The parameter directory "+param_path+" is not found.")
    sys.exit(1)
is_in_path = False
for i in sys.path:
    if i == sys_path:
        is_in_path = True
if not is_in_path:
    sys.path.append(sys_path)
import Bio.PDB as struct
import input_output as IO
import optparse

def main():
    usage = "usage: %prog [options] arg"
    d = "This script is used to work with CIF files."
    opt_parser = optparse.OptionParser(usage,description=d)
    group = optparse.OptionGroup(opt_parser, "To get a CIF's file Summary")
    group.add_option("--summary", action="store_true", help="Extracts Basic General information from CIF files. (Not Pr\
ogrammed yet)")
    group.add_option("--inp1",type="str",help="Path to input cif file.")
    opt_parser.add_option_group(group)

    group = optparse.OptionGroup(opt_parser, "To extract PDBs from a CIF file")
    group.add_option("--extract", action="store_true", help="It outputs PDB files from CIF files by Models or specific \
Chain groups, and it gets rid of unwanted HETEROATOMS.")
    group.add_option("--chains", action="store_true", help="Extracts PDBs by groups of chains using identifiers. Must s\
pecify groups. Read help for '--groups' option for method to specify PDB outputs by chain.")
    group.add_option("--models", action="store_true", help="Extracts PDBs by models. Outputs all models separately.")
    group.add_option("--inp2",type="str",help="Path to input cif file.")
    group.add_option("--groups", type="str",help="Chain groups to be extracted together as separate PDBs. User must spe\
cify the chains to be output together to each PDB. There is no checking for completeness. Ex: \"ABC,DEF\" will output t\
wo pdbs from a CIF with six chains, or \"AB,CD,EF\" will output three pdbs from a CIF with six chains.")
    opt_parser.add_option_group(group)
    options, args = opt_parser.parse_args()
############################################ Check Options Entered #####################################################
    if options.summary:
        if options.extract:
            opt_parser.error("Two options can't be selected at the same time. Enter -h for help.")
            sys.exit(1)
    inFile = ""
    if options.summary:
        inFile = options.inp1
    elif options.extract:
        inFile = options.inp2
    if not os.path.exists(inFile):
        print("Error: CIF input file does not exist.")
        print("Type -h or --help for description and options.")
        sys.exit(1)
    cif = IO.cif()
    if options.summary:
        cif.CIF_summary()
    elif options.extract:
        cif.PDBs_from_CIF(inFile,options.models,options.chains,options.groups)
if __name__ == '__main__':
    main()
