import argparse

import em.charmm.gen.inputs as IO
#import em.tools.input_output as IO

def register_parser(subparsers):
    parser = subparsers.add_parser('prepare', usage=usage(), description=description())
    add_arguments(parser)


def add_arguments(parser):
    parser.add_argument("--input", metavar="FILE", help="Input File.", required=True)
    parser.add_argument('--terminals', type=str, help=
    "Format: A,none,CTER:B,ACE,none (e.i. for chain A, not ACE and a CTER. For chain B, an Ace, and no CTER.)")
    parser.set_defaults(func=run)


def run(options):
    IO.prepare_pdb_for_charmm(options.input,options.terminals)


def description():
    return '''This command takes a PDB file that has been output from gen_csv.py, and it prepares it so
        hydrogen and heavy atoms are added when missing. It first runs reduce to identify HIS residues
        as any of the three histidine types in CHARMM, HSD, HSE or HSP, depending on their particular 
        protonation state. Then a CHARMM script is automatically generated and run to check for other
        missing atoms in the protein, and to add ACE or CTER terminasl to the protein ends.'''

def usage():
    return '\npdb_cif.py prepare --input 1BRS.pdb\n'


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
