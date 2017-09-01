#!/usr/bin/python
# -*- coding: utf-8 -*-
import em.tools.input_output as IO

import optparse

usage = "usage: %prog --inputfile [input file] [options] [args]"
d = "This script is used to modify structures from PDB or CIF files as input."
opt_parser = optparse.OptionParser(usage,description=d)
group = optparse.OptionGroup(opt_parser,"An input file is necessary for all command options.")
group.add_option("--input", metavar="FILE", help="Input File.")
opt_parser.add_option_group(group)

group = optparse.OptionGroup(opt_parser,"Command options that work on both PDB or CIF files types. Type 'more' after an\
 option for \n  instructions to run these commands")
group.add_option("--gaps", action="store_true", metavar="FILE", help="To find missing regions in structures.")
group.add_option("--align", action="store_true", metavar="FILE", help="To align two structures and add atoms to missing\
 regions from one structure to the other.")
group.add_option("--summary", action="store_true", metavar="FILE", help="For relevant information on structure files.")
group.add_option("--extract", action="store_true", help="Gets PDBs grouped by models or chains from PDB or CIF files.")
opt_parser.add_option_group(group)

group = optparse.OptionGroup(opt_parser,"Command options that work on PDB files only. Type 'more' after an option for i\
nstructions\n  to run these commands")
group.add_option("--fixpdb", action="store_true", metavar="FILE", help="Corrects biopython formating of PDB structures \
so that PDBs can be read by CHARMM.")
group.add_option("--prepare", action="store_true", metavar="FILE", help="Assigns right histidine type to PDBs and gener\
ates files for preparing structures in CHARMM.")
group.add_option("--maxmin", action="store_true", metavar="FILE", help="Outputs maximin and minimum values for XYZ coor\
dinates of a PDB file.")
opt_parser.add_option_group(group)


group = optparse.OptionGroup(opt_parser,"FOR MORE HELP: Type 'more' after any of the above options to get more detail\
ed\n  instructions about the arguments used by each command")
# None for gaps
# Five for align
group.add_option("--refatoms", type="str",help=optparse.SUPPRESS_HELP)
group.add_option("--fit", metavar="FILE",type="str",help=optparse.SUPPRESS_HELP)
group.add_option("--fitatoms", type="str",help=optparse.SUPPRESS_HELP)
group.add_option("--out", metavar="FILE",type="str",help=optparse.SUPPRESS_HELP)
group.add_option("--addatoms",default = "",type="str",help=optparse.SUPPRESS_HELP)
# None for summary
# Three for extract
group.add_option("--chains", action="store_true", help=optparse.SUPPRESS_HELP)
group.add_option("--models", action="store_true", help=optparse.SUPPRESS_HELP)
group.add_option("--groups", type="str",help=optparse.SUPPRESS_HELP)
# Three for fixpdb
group.add_option('--frm', type="int", action = "store", default = 72, help = optparse.SUPPRESS_HELP)
group.add_option('--to', type="int",action = "store", default = 21, help = optparse.SUPPRESS_HELP)
# Two for prepare
group.add_option('--crdout', metavar="FILE", type='str', help=optparse.SUPPRESS_HELP)
group.add_option('--seqfix', type='str', help=optparse.SUPPRESS_HELP)
opt_parser.add_option_group(group)
options, args = opt_parser.parse_args()
# None for maxmin
########################################################################################################################
# Just Check for the argument 'more' to output more help.
run_command = True
for i in args:
    if i == "more":
        run_command = False
########################################################################################################################
if options.gaps:
    if run_command:
        IO.gap_report(options.input)
    else:
        print("Description and usage of --gaps:")
        print("    This command detects gaps in the crystal structure of a protein. The search for gaps")
        print("    is based on the amino acid residue number. It will detect missing amino acids at the")
        print("    N-terminal of the protein if the structure's residue numbers begin with a number greater")
        print("    than one, and it will miss gaps in the C-terminal in PDB and CIF files because it assumes")
        print("    that the last residue number is actually the last one of the crystal structure.(TODO: use")
        print("    sequence iformation found in CIF files to detect gaps in the C-terminal end of the protein)")
        print("        --gaps                 Flag to signal the program to check for gaps.")
        print("        --input = FILE         Follow this option by the path to a file.")
        print("    Example: pdb_cif.py --gaps --input 1BRS.pdb")
elif options.align:
    if run_command:
        #group = optparse.OptionGroup(opt_parser,"")
        #group.add_option("--refatoms", type="str",help="")
        #group.add_option("--fit", metavar="FILE",type="str",help="")
        #group.add_option("--fitatoms", type="str",help="")
        #group.add_option("--out", metavar="FILE",type="str",help="")
        #group.add_option("--addatoms",default = "",type="str",help="")
        #opt_parser.add_option_group(group)
        #options, args = opt_parser.parse_args()
        IO.align_pdbs(options.input,options.fit,options.refatoms,options.fitatoms,options.out,options.addatoms)
    else:
        print("Description and usage of --align:")
        print("    This command aligns regions according to specific atoms and amino acid contigous sections of two")
        print("    crystal structures. The fit structure is moved over the input structure according to reference and")
        print("    fit atoms selections respectively. If the opion to addatoms is included, it adds atoms from the")
        print("    reference to the fit structure's specified region. The format to specify reference and fit atoms")
        print("    is ATOMTYPE,CHAIN,FROMAMINOACID,TOAMINOACID Ex: CA,B,2,6. Multiple regions for aligment can be")
        print("    defined by separating them with a column. Ex: CA,B,2,6:CA,B,10,26 ")
        print("        --align              Flag to signal the programm to align two structures.")
        print("        --input = FILE       Follow this option by the path to a file.")
        print("        --refatoms = String  String with atoms for alignment in reference structure.")
        print("        --fit = FILE         Path to file for the fit structure.")
        print("        --fitatoms = String  String with atoms for alignment in fit structure.")
        print("        --out = FILE         Path to output file of reference structure aligned over the fit struscture")
        print("        --addatoms = String  Optional option, default no atoms. Cuts and pastes all atoms corresponding")
        print("                             to amino acids defined by a string such as A,202,217:A,202,217. The region")
        print("                             before the column corresponds to the region to be cut from the reference")
        print("                             structure, and the region defined after the colum to are to be pasted to")
        print("                             the fit structure. Unlike the string that defines reference and fit")
        print("                             atoms, atom type is not needed, and all atoms in the amino acid are added")
        print("    Example: pdb_cif.py --align --input 1GIA.cif --refatoms CA,B,199,201:CA,B,218,220 --fit 1GDD.pdb --f\
itatoms CA,B,199,201:CA,B,218,220 --out 1GDD_aligned_completed_with_1GIA.pdb --addatoms A,202,218:A,202,218")
elif options.summary:
    if run_command:
        IO.prepare_pdb_for_charmm(options.inputfile,options.crdout,options.seqfix)
    else:
        print("Description and usage of --summary:")
        print("    This command gives a brief or summarized output of the structure as a guide for other commands.")
        print("    The output gives information that is relevant to biopython, and that options that will help decide")
        print("    command line inputs for other options for this program such as --align and --extract.")
        print("    The output gives information on number of models and chains.(TODO: This option is not coded yet)")
        print("        --summary              Flag to signal the program to do a summary.")
        print("        --input = FILE         Follow this option by the path to a file.\n")
        print("    Example: pdb_cif.py --summary --input 1HIU.pdb")
elif options.extract:
    run_command = True
    for i in args:
        if i == "more":
            run_command = False
    if run_command:
        #group = optparse.OptionGroup(opt_parser,"")
        #group.add_option("--chains", action="store_true", help="")
        #group.add_option("--models", action="store_true", help="")
        #group.add_option("--groups", type="str",help="")
        #opt_parser.add_option_group(group)
        IO.PDBs_from_CIF(options.input,options.models,options.chains,options.groups)
    else:
        print("Description and usage of --extract:")
        print("    It extracts models or chain groups in the structure to separate PDB files, and it also gets rid of")
        print("    unwanted HETEROATOMS if they have a different chain identifier. The flags --chains and --models")
        print("    cannot be input at the same time. Use --summary to know what models and chains are inside the structure.")
        print("        --exract             Flag to extract models or chains from a structure file.")
        print("        --input = FILE       Follow this option by the path to a file.")
        print("        --chains             Flag to extract PDBs by groups of chains using chain identifiers.")
        print("        --models             Flag to extract all models to separate outputs. No need to specify groups")
        print("        --groups             Chain groups to be extracted together as separate PDBs. User must specify")
        print("                             the chains to be output together to each PDB. There is no checking for")
        print("                             completeness. Ex: \"ABC,DEF\" will output two pdbs from a CIF with six")
        print("                             chains, or \"AB,CD,EF\" will output three pdbs from a CIF with six chains.")
        print("    Example: pdb_cif.py --extract --input 1BRS.pdb --chains --groups AB,CD,EF")
        print("    Example: pdb_cif.py --extract --input 1BRS.pdb --models")
elif options.fixpdb:
    if run_command:
        #group = optparse.OptionGroup(opt_parser,"")
        #group.add_option("--fixpdb", action="store_true", help="")
        #group.add_option('--frm', type="int", action = "store", default = 72, help = "")
        #group.add_option('--to', type="int",action = "store", default = 21, help = "")
        #opt_parser.add_option_group(group)
        IO.fix_pdb_from_CHARMM(options.to,options.frm)
    else:
        print("Description and usage of --fixpdb:")
        print("    It fixes a PDB file that was output by CHARMM. The chain identifier is placed by CHARMM in a column")
        print("    that is different from what Biopython and Entropy Maxima can handle. The chain identifier is swapped")
        print("    from column 72 to column 21 in the PDB file.")
        print("        --fixpdb             Flag to extract models or chains from a structure file.")
        print("        --input = FILE       Follow this option by the path to a file.")
        print("        --frm = INT          This option takes an integer value for initial column. Default is 72.")
        print("        --to = INT           This option takes an integer value for destination column. Default is 21.")
        print("    Example: pdb_cif.py --fixpdb --input 2GIA.pdb ")
        print("    Example: pdb_cif.py --fixpdb --input 2GIA.pdb --frm 72 --to 21")
elif options.prepare:
    if run_command:
        #group = optparse.OptionGroup(opt_parser,"")
        #group.add_option('--crdout', metavar="FILE", type='str', help="")
        #group.add_option('--seqfix', type='str', help="")
        #opt_parser.add_option_group(group)
        IO.prepare_pdb_for_charmm(options.input,options.crdout,options.seqfix)
    else:
        print("Description and usage of --prepare:")
        print("    This command takes a PDB file that has been run through the reduce program to determine missing")
        print("    hydrogen atoms in hystidines, and based on the added hydrogens by reduce, to identify HIS residues")
        print("    as any of the three histidine types in CHARMM, HSD, HSE or HSP. Runing reduce before this program")
        print("    is strongly recomended because hydrogen atoms added without running reduce are likely to be wrong.")
        print("    The input PDB file is outputed as a CRD structure file which is the type CHARMM takes as input.")
        print("    Prepare has an optional feature that outputs to files needed by CHARMM to finish setting up the")
        print("    structure. The files generated have the following suffixes .SEQ and _FIXERS.INP.")
        print("        --crdout = FILE      Path and name to output file with a 'crd' suffix. This file is needed for")
        print("                             CHARMM to prepare a structure for simulation.")
        print("        --seqfix = String    This option takes only a 'Yes' or 'no', default is no. This file is needed")
        print("                             for CHARMM to an amino acid sequence of a structure for simulations.")
        print("    Example: pdb_cif.py --prepare --input 2GIA.pdb --crdout 1GIA.crd --seqfix yes")
elif options.maxmin:
    if run_command:
        IO.min_max(options.input)
    else:
        print("Description and usage of --maxmin:")
        print("    This command gives the maximum and minim XYZ coordinate values of all atoms in the structure.")
        print("    This information is used to add a water box around the protein that extends a defined number of")
        print("    agnstroms from the protein. This values are used to determine the right box dimensions for ")
        print("    minimum-image conversion.")
        print("        --maxmin                 Flag to signal the program to ")
        print("        --inputfile = FILE       Follow this option by the path to a file.\n")
        print("    Example: pdb_cif.py --maxmin --input 1HIU.pdb")
else:
    print("Follow the options by --help for instructions. Ex: %prog --help")
    options, args = opt_parser.parse_args()
