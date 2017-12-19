import argparse

import em.tools.input_output as IO

def register_parser(subparsers):
    parser = subparsers.add_parser('maxmin', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--input", metavar="FILE", help="Input File.", required=True)
    parser.set_defaults(func=run)

def run(options):
    IO.min_max(options.input)

def description():
    return '''This command gives the maximum and minim XYZ coordinate values of all atoms in the structure.
    This information is used to add a water box around the protein that extends a defined number of
    agnstroms from the protein. This values are used to determine the right box dimensions for
    minimum-image conversion.'''

def usage():
    return '\npdb_cif.py maxmin --input 1BRS.pdb\n'

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
