import em.tools.input_output as IO


def register_parser(subparsers):
    parser = subparsers.add_parser('fixpdb', usage=usage(), description=description())
    add_arguments(parser)


def add_arguments(parser):
    parser.add_argument("--input", metavar="FILE", help="Input File.", required=True)
    parser.add_argument('--frm', type=int, action="store", default=72, help="")
    parser.add_argument('--to', type=int, action="store", default=21, help="")
    parser.set_defaults(func=run)


def run(options):
    IO.fix_pdb_from_CHARMM(options.to, options.frm)


def description():
        return '''It fixes a PDB file that was output by CHARMM. The chain identifier is placed by CHARMM in a column
        that is different from what Biopython and Entropy Maxima can handle. The chain identifier is swapped
        from column 72 to column 21 in the PDB file.
            --fixpdb             Flag to extract models or chains from a structure file.
            --input = FILE       Follow this option by the path to a file.
            --frm = INT          This option takes an integer value for initial column. Default is 72.
            --to = INT           This option takes an integer value for destination column. Default is 21.'''


def usage():
    return 'pdb_cif.py --fixpdb --input 2GIA.pdb \n' \
           'pdb_cif.py --fixpdb --input 2GIA.pdb --frm 72 --to 21'


if __name__ == '__main__':
    pass
