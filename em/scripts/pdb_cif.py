#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 22:53:59 2017

@author: noel
"""
import os
import sys
import em.tools.input_output as IO

import optparse

def main():
    usage = "usage: %prog --inputfile [input file] [options]"
    d = "This script is used to work with PDB or CIF files as input."
    opt_parser = optparse.OptionParser(usage,description=d)
    opt_parser.add_option("--inputfile", metavar="FILE", help="PDB or CIF input FILE needed for all operations li\
sted below. Do not mix options from different operations. Each operation only works with its options enumerated below.")
    group = optparse.OptionGroup(opt_parser, "1.To get a CIF or PDB structure's gaps")
    group.add_option("--gaps", action="store_true", help=" Flag to get amino acid structural gaps in X-structures")
    opt_parser.add_option_group(group)
    
    group = optparse.OptionGroup(opt_parser,"2.To align regions according to specified reference/input and fit sections\
 of\n    the structures. The fit structure is moved over the input structure. If \n    the optional --addatoms is inclu\
ded, it adds atoms from the reference to\n    the fit structure's specified region. The format to specify --refatoms \
\n    and --fitatoms is ATOMTYPE,CHAIN,FROMAMINOACID,TOAMINOACID Ex: CA,B,2,6")
    group.add_option("--align", action="store_true", help="Flag to aligns two structures.")
    group.add_option("--refatoms", type="str",help="String with atoms for alignment in reference structure.")
    group.add_option("--fit", metavar="FILE",type="str",help="File path of structure to fit to reference structure.")
    group.add_option("--fitatoms", type="str",help="String with atoms for alignment in fit structure.")
    group.add_option("--out", metavar="FILE",type="str",help="File output name for the translated fit structure.")
    group.add_option("--addatoms",default = "",type="str",help="Atoms to be added from reference to fit. Default is no\
 atoms added.")
    opt_parser.add_option_group(group)
    
    group = optparse.OptionGroup(opt_parser,"3.To prepare a pdb file that has been run through the reduce program to a\
dd\n    hydrogen atoms to the Histidines in the most likely configurations. Based on\n    that configuration, the progr\
am ranames HIS to HSE, HSD or HSP. It also\n    generates *.SEQ and *_FIXERS.INP for proper structure preparation by CH\
ARMM")
    group.add_option("--prepare", action="store_true", help="Flag to prepare structure files for CHARMM setup.")
    group.add_option('--crdout', metavar="FILE", type='str', help="Output CRD file. Coordinates are the same as PDB's.")
    group.add_option('--seqfix', type='str', help='If followed by yes, it generates SEQ and _FiXRES files. if followed\
by no or absent, nothing happens.')
    opt_parser.add_option_group(group)
    
    group = optparse.OptionGroup(opt_parser,"4.To fix a pdb file that was output by CHARMM. The chain identifier is pl\
aced\n    by CHARMM in a column that is different from what Biopython and this script\n    are programmed to handle. Th\
is option will onlly work on PDB files, and not\n    on CIF files")
    group.add_option("--fixpdb", action="store_true", help="Flag to fix PDB files output by CHARMM")
    group.add_option('--frm', type="int", action = "store", default = 72, help = 'Optional arg if not present it is col\
umn 72 in the pdb file. It should be at the column where the chain identifier is located.')
    group.add_option('--to', type="int",action = "store", default = 21, help = 'Optional arg if not present it is colum\
n 21 in the pdb file. It should be at the column where the chain identifier has to be placed')
    opt_parser.add_option_group(group)
    
    group = optparse.OptionGroup(opt_parser,"5.To find the maximum and minimum coordinate values in each dimension")
    group.add_option("--minmax", action="store_true", help="Flag to get minimum and maximum XYZ values of a structure.")
    opt_parser.add_option_group(group)
    
    group = optparse.OptionGroup(opt_parser,"6.Extracts Basic General information from structural files.(Not Programmed yet)")
    group.add_option("--summary", action="store_true", help="FLag to output summary.")
    opt_parser.add_option_group(group)

    group = optparse.OptionGroup(opt_parser, "7.To extract models or chain groups in the structure to separate PDB \
files, and\n    it gets rid of unwanted HETEROATOMS.\n    The flags --chains and --models cannot be input at the same time")
    group.add_option("--extract", action="store_true", help="Flag to indicate extraxtion.")
    group.add_option("--chains", action="store_true", help="Flag to extract PDBs by groups of chains using \
identifiers. Must specify --groups options and a string describing the grouping of chains in each output PDB by a comma\
 separated string as follows: \"ABC,DEF\" will output two pdbs from a structure with six chains, or \"AB,CD,EF\" will \
output three pdbs from a CIF with six chains")
    group.add_option("--models", action="store_true", help="Flag to extract all models that are outputs separately.")
    group.add_option("--groups", type="str",help="Chain groups to be extracted together as separate PDBs. User must spe\
cify the chains to be output together to each PDB. There is no checking for completeness. Ex: \"ABC,DEF\" will output  t\
wo pdbs from a CIF with six chains, or \"AB,CD,EF\" will output three pdbs from a CIF with six chains.")
    opt_parser.add_option_group(group)
    options, args = opt_parser.parse_args()
    # The program now does a few check to make sure options above were entered correctly
    # First, make sure there is at least an input file in the arguments.
    if not options.inputfile:
        opt_parser.error("An --inputfile option is necessary.")
        sys.exit(1)
    # Second, a general check to make sure only one of the main options are selected.
    options_dictionary = {'gaps':options.gaps,'align':options.align,'prepare':options.prepare,'fixpdb':options.fixpdb,\
                          'minmax':options.minmax,'summary':options.summary,'extract':options.extract}
    countTrue = 0
    for i in options_dictionary:
        if options_dictionary[i]:
            countTrue += 1
    if countTrue > 1:
        print("Error: More than one option selected.")
        for i in options_dictionary:
            if options_dictionary[i]:
                countTrue += 1
                print("Opption "+i+" "+str(options_dictionary[i]))
        opt_parser.error("Only one option permited at the time.")
        sys.ecit(1)
    # Third: Check path to input file is correct
    if not os.path.exists(options.inputfile):
        print("Error: Path to the input file does not exist,")
        print("or input file name is incorrect.")
        sys.exit(1)
    # All checks passed, now do the work! 
    if options.gaps:
        pdb = IO.pdb()
        pdb.gap_report(options.inputfile) 
    elif options.align:
        pdb = IO.pdb()
        pdb.align_pdbs(options.inputfile,options.fit,options.refatoms,options.fitatoms,options.out,options.addatoms)
    elif options.prepare:
        pdb = IO.pdb()
        pdb.prepare_pdb_for_charmm(options.inputfile,options.crdout,options.seqfix)
    elif options.fixpdb:
        pdb = IO.pdb(options.inputfile)
        pdb.fix_pdb_from_CHARMM(options.to,options.frm)
    elif options.minmax:
        pdb = IO.pdb()
        pdb.min_max(options.inputfile)
    elif options.summary:
        cif = IO.cif()
        cif.CIF_summary()
    elif options.extract:
        cif = IO.cif()
        cif.PDBs_from_CIF(options.inputfile,options.models,options.chains,options.groups)
if __name__ == '__main__':
    main()
