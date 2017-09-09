#!/usr/bin/python

"""
Created on Fri Jun 24 16:49:07 2016 
@author: noel
For a description of the program, type:
python mmcif_to_pdb.py 
Future functionality: The following functions only make sense when 
working with the whole Protein Databank.
1. In the future this file could be used to read FASTA files and compare to the
   PDB sequence information for a check across the protein databank as a sanity
   check.
2. It could check other information in the header of the CIF file as a first
   check.
"""
import os
import sys
from Bio.PDB.PDBParser import PDBParser
import em.tools.CHARMM_Parser as CP
import em.tools.Super_Structures as SS
import em.tools.input_output as IO
import optparse

def main():
    usage = "usage: %prog [options] arg"
    d = "This program reads a CIF file and checks that all residues in the file\
         are found in the CHARMM top_27 parameters. Residues not found \
         are added to the structure. The full structure is outputed to a \
         CSV file where Charmm, CIF and additional information is stored. \
         Added residues are copied from a peptide structure with all amino acids \
         present in the local CHARMM parameters files with fixed dihedral angles. \
         Info in the CSV file is all there is to explore the conformational \
         space of added atoms."
    opt_parser = optparse.OptionParser(usage,description=d)
    
    group = optparse.OptionGroup(opt_parser, "Generates CSV and PDB files for each model from a CIF file.")
    group.add_option("--fromcif", action="store_true", help="Flag to generate a CSV frile from a CIF file.")
    group.add_option("-i","--cif",type="str",help="Path to input cif file.")
    group.add_option("-o","--out1", type="str",help="Path to output csv.")
    group.add_option("-p","--pep", type="str",help="Path to CHARMM peptide file.")
    opt_parser.add_option_group(group)
    
    group = optparse.OptionGroup(opt_parser, "Generates a CSV file from CRD and PSF files.")
    group.add_option("--frompsfcrd", action="store_true", help="Flag to generates a CSV frile from a CRD and PSF file.")
    group.add_option("-f","--psf",type="str",help="Path to input PSF file in XPLOR format.")
    group.add_option("-d","--crd", type="str",help="Path to input CRD file.")
    opt_parser.add_option_group(group)
    
    options, args = opt_parser.parse_args()
############################################  Options Entered ##########################################################
    if options.fromcif:
        if options.frompsfcrd:
            opt_parser.error("Two option flags can't be selected at the same time. Enter -h for help.")
########################################################################################################################
    if options.fromcif:
        if not os.path.exists(options.cif):
            print "Error: File path for input file does not exist."
            print("Type -h or --help for description and options.")
            sys.exit(1)
        if not os.path.exists(options.pep):
            print "Error: File with CHARMM's generated peptides does not exists."
            print("Type -h or --help for description and options.")
            sys.exit(1)
        params = CP.read_charmm_FF()
        parser2 = PDBParser(QUIET = True)
        p1 = parser2.get_structure('Peptides', options.pep)
        ###########################################################################
        # The peptide construct is build with charmm so corrections for some atom
        # names to PDB/Databank atom types is needed.
        # TODO: this might not be necessary as the correction and inv_correction dictionary in Super Structure takes care of it.
        # Check before removing the correction here.
        for i in p1.get_models():
            for j in i.get_chains():
                for k in j.get_residues():
                    for l in k.get_atom():
                        if k.get_resname() == 'ILE' and l.get_id() == 'CD':
                            l.name = 'CD1'
                            l.id = 'CD1'
        ###########################################################################
        # Create Super Structure
        myCIF = SS.Super_Structure(params, options.cif, 'setup')
        myCIF.build_pep_and_anchers(p1)
        myCIF.read_dict_into_dataframes()
        myCIF.check_models()
        myCIF.create_super_structure_df()
        ###########################################################################
        # Find missing residues to add to the Super Structure. Missing residues
        # are group in lists of contiguous residues and aded to another list.
        myCIF.build_missing_aa()
        file_name = os.path.basename(options.cif).split('.')[0]
        myCIF.write_csv('',file_name)
        #outPDB = IO.pdb()
        IO.write_pdb(myCIF,'',file_name,'all')
    if options.frompsfcrd:
        if not os.path.exists(options.psf):
            print "Error: File path for PSF file does not exist."
            print("Type -h or --help for description and options.")
            sys.exit(1)        
        if not os.path.exists(options.crd):
            print "Error: File path for CRD file does not exist."
            print("Type -h or --help for description and options.")
            sys.exit(1)
        directory, filename = os.path.split(options.crd)
        crd_file = IO.crd(options.crd)
        psf_file = IO.psf(options.psf)
        file_name = filename.split('.')[0]
        ################################################################################################################
        ###################### After reading files, Generate and Index a Super Structure  ##############################
        params = CP.read_charmm_FF()
        myCSV = SS.Super_Structure(params, directory,'charmm_input')
        # At this point, a XPLOR psf could only have been creted from a complete structure, so no worries of gaps.
        myCSV.create_super_structure_df_from_CRD_PSF(crd_file,psf_file)
        myCSV.write_csv(directory,file_name)
if __name__ == '__main__':
    main()
