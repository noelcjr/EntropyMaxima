#!/usr/bin/python

"""
Created on Mon Nov 21 22:50:36 2016

@author: noel
"""
import os
import sys
import em.tools.input_output as IO

import optparse

def main():
    usage = "usage: %prog [options] arg"
    d = "This script is used to work with PDB files."
    opt_parser = optparse.OptionParser(usage,description=d)
    group = optparse.OptionGroup(opt_parser, "1. To get a PDB's structural gaps")
    group.add_option("--gaps", action="store_true", help="Reports structural gaps in the PDB file from X-structures.")
    group.add_option("-i","--inp", type="str",help="Path to PDB file for gap report.")
    opt_parser.add_option_group(group)
    group = optparse.OptionGroup(opt_parser, "2. To align regions according to specified reference and fit sections of \
the PDB. If addatoms option is included, it adds a specified region from the reference structure to the fit structure. \
MAKE SURE atoms added equals the size of the gap in the fit structure.")
    group.add_option("--align", action="store_true", help="Aligns two PDBs.")
    group.add_option("--refA",type="str",help="Path to input PDB reference file.")
    group.add_option("--refatomsA", type="str",help="Ref atoms for alignment.")
    group.add_option("--fitA",type="str",help="Path to input PDB fit file.")
    group.add_option("--fitatomsA", type="str",help="Fit atoms for alignment.")
    group.add_option("--outA",type="str",help="Name of output file for the Fit structure.")
    group.add_option("--addatomsA",default = "",type="str",help="Atoms to be added from reference to fit. Default is no\
 atoms added.")
    opt_parser.add_option_group(group)
    group = optparse.OptionGroup(opt_parser,"3. To prepare a pdb file that has been run through the reduce program to a\
dd hydrogen atoms to the Histidines in the most likely configurations. Based on that configuration, the program ranmes \
HIS to HSE, HSD or HSP. It also generates *.SEQ and *_FIXERS.INP for proper structure preparation")
    group.add_option("--prepare", action="store_true", help="Prepares PDB files for CHARMM setup.")
    group.add_option('--pdbin1', type='str', help="PDB file to be prepared and modified.")
    group.add_option('--crdout', type='str', help="Output CRD file. Coordinates are identical to PDB's.")
    group.add_option('--seqfix', type='str', help='If followed by yes, it generates SEQ and _FiXRES files. if followed \
by no or absent, nothing happens.')
    opt_parser.add_option_group(group)
    group = optparse.OptionGroup(opt_parser,"4. To fix a pdb file that was output by CHARMM. The chain identifier is pl\
lace by CHARMM in a column that is different from what Biopython and Pipeprot are programmed to handle.\n")
    group.add_option("--fixpdb", action="store_true", help="Prepares PDB files for CHARMM setup")
    group.add_option('--pdbin2', type='str', help='Path to PDB file to be prepared and modified.')
    group.add_option('--frm', type="int", action = "store", default = 72, help = 'Optional arg if not present it is col\
umn 72 in the pdb file. It should be at the column where the chain identifier is located.')
    group.add_option('--to', type="int",action = "store", default = 21, help = 'Optional arg if not present it is colum\
n 21 in the pdb file. It should be at the column where the chain identifier has to be placed.')
    opt_parser.add_option_group(group)
    group = optparse.OptionGroup(opt_parser,"5. To find the maximum and minimum coordinate in each dimension of a pdb f\
ile")
    group.add_option("--minmax", action="store_true", help="Gets minimum and maximum XYZ values of a PDB file.")
    group.add_option("--pdbin3", type='str', help='Path to PDB file to be analyzed.')
    # TODO
    #group.add_option("--heavy", action="store_true", help="For minimum and maximum of only non-hydrogen/heavy atoms.")
    options, args = opt_parser.parse_args()
    if options.gaps:
        if options.align:
            if options.prepare:
                if options.fixpdb:
                    opt_parser.error("Two options can't be selected at the same time. Enter -h for help.")
                    sys.exit(1)
    if options.inp:
        pdb = IO.pdb()
        pdb.gap_report(options.inp)
    elif options.align:
        pdb = IO.pdb()
        pdb.align_pdbs(options.refA,options.fitA,options.refatomsA,options.fitatomsA,options.outA,options.addatomsA)
    elif options.prepare:
        pdb = IO.pdb()
        pdb.prepare_pdb_for_charmm(options.pdbin1,options.crdout,options.seqfix)
    elif options.fixpdb:
        pdb = IO.pdb(options.pdbin2)
        pdb.fix_pdb_from_CHARMM(options.to,options.frm)
    elif options.minmax:
        pdb = IO.pdb()
        pdb.min_max(options.pdbin3)

if __name__ == '__main__':
    main()
