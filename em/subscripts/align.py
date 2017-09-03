import argparse

import em.tools.input_output as IO


def register_parser(subparsers):
    parser = subparsers.add_parser('align', usage=usage(), description=description())
    add_arguments(parser)


def add_arguments(parser):
    parser.add_argument("--input", metavar="FILE", help="Input File.", required=True)
    parser.add_argument("--refatoms", type=str, required=True)
    parser.add_argument("--fit", metavar="FILE", type=str, required=True)
    parser.add_argument("--fitatoms", type=str, required=True)
    parser.add_argument("--out", metavar="FILE", type=str)
    parser.add_argument("--addatoms", default="", type=str)
    parser.set_defaults(func=run)


def run(options):
    IO.align_pdbs(options.input, options.fit, options.refatoms, options.fitatoms, options.out, options.addatoms)


def description():
    return '''Description and usage of align:
        This command aligns regions according to specific atoms and amino acid contigous sections of two
        crystal structures. The fit structure is moved over the input structure according to reference and
        fit atoms selections respectively. If the opion to addatoms is included, it adds atoms from the
        reference to the fit structure's specified region. The format to specify reference and fit atoms
        is ATOMTYPE,CHAIN,FROMAMINOACID,TOAMINOACID Ex: CA,B,2,6. Multiple regions for aligment can be
        defined by separating them with a column. Ex: CA,B,2,6:CA,B,10,26
            --refatoms = String  String with atoms for alignment in reference structure.
            --fit = FILE         Path to file for the fit structure.
            --fitatoms = String  String with atoms for alignment in fit structure.
            --out = FILE         Path to output file of reference structure aligned over the fit struscture
            --addatoms = String  Optional option, default no atoms. Cuts and pastes all atoms corresponding
                                 to amino acids defined by a string such as A,202,217:A,202,217. The region
                                 before the column corresponds to the region to be cut from the reference
                                 structure, and the region defined after the colum to are to be pasted to
                                 the fit structure. Unlike the string that defines reference and fit
                                 atoms, atom type is not needed, and all atoms in the amino acid are added'''


def usage():
    return 'pdb_cif.py align --input 1GIA.cif --refatoms CA,B,199,201:CA,B,218,220 --fit 1GDD.pdb --f\
itatoms CA,B,199,201:CA,B,218,220 --out 1GDD_aligned_completed_with_1GIA.pdb --addatoms A,202,218:A,202,218'


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
