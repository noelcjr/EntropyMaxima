import argparse

import em.tools.input_output as IO


def register_parser(subparsers):
    parser = subparsers.add_parser('fixpdb', usage=usage(), description=description())
    add_arguments(parser)


def add_arguments(parser):
    parser.add_argument("--input", metavar="FILE", help="Input File.", required=True)
    parser.add_argument('--frm', type=int, action="store", default=72, help="This option takes an integer value for initial column. Default is 72.")
    parser.add_argument('--to', type=int, action="store", default=21, help="This option takes an integer value for destination column. Default is 21.")
    parser.set_defaults(func=run)


def run(options):
    IO.fix_pdb_from_CHARMM(options.input,options.to, options.frm)


def description():
        return '''It fixes a PDB file that was output by CHARMM. The chain identifier is placed by CHARMM in a column
        that is different from what Biopython and Entropy Maxima can handle. The chain identifier is swapped
        from column 72 to column 21 in the PDB file.
        '''

def usage():
    return '\npdb_cif.py fixpdb --input 2GIA.pdb \n' \
           'pdb_cif.py fixpdb --input 2GIA.pdb --frm 72 --to 21'


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
