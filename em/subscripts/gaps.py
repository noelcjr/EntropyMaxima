import argparse

from em.tools.input_output import gap_report


def register_parser(subparsers):
    parser = subparsers.add_parser('gaps', usage=usage(), description=description())
    add_arguments(parser)


def add_arguments(parser):
    parser.add_argument("--input", metavar="FILE", help="Input File.", required=True)
    parser.set_defaults(func=run)


def run(args):
    gap_report(args.input)


def description():
    return '''
    This command detects gaps in the crystal structure of a protein. The search for gaps
    is based on the amino acid residue number. It will detect missing amino acids at the
    N-terminal of the protein if the structure's residue numbers begin with a number greater
    than one, and it will miss gaps in the C-terminal in PDB and CIF files because it assumes
    that the last residue number is actually the last one of the crystal structure.(TODO: use
    sequence iformation found in CIF files to detect gaps in the C-terminal end of the protein)
    '''


def usage():
    return 'pdb_cif.py gaps --input 1BRS.pdb'

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)