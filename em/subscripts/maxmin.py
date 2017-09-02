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
    return '''This command gives the maximum and minim XYZ coordinate values of all atoms in the structure.")
    This information is used to add a water box around the protein that extends a defined number of")
    agnstroms from the protein. This values are used to determine the right box dimensions for ")
    minimum-image conversion.")'''


def usage():
    return 'pdb_cif.py --extract --input 1BRS.pdb --chains --groups AB,CD,EF\n' \
           'pdb_cif.py --extract --input 1BRS.pdb --models'


if __name__ == '__main__':
    pass