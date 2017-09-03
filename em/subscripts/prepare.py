import argparse

import em.tools.input_output as IO


def register_parser(subparsers):
    parser = subparsers.add_parser('prepare', usage=usage(), description=description())
    add_arguments(parser)


def add_arguments(parser):
    parser.add_argument("--input", metavar="FILE", help="Input File.", required=True)
    parser.add_argument('--crdout', metavar="FILE", type=str, help=
    "Path and name to output file with a 'crd' suffix. This file is needed for CHARMM to prepare a structure for simulation.")
    parser.add_argument('--seqfix', type=str, help=
    "This option takes only a 'Yes' or 'no', default is no. This file is needed for CHARMM to an amino acid sequence of a structure for simulations.")
    parser.set_defaults(func=run)


def run(options):
    IO.prepare_pdb_for_charmm(options.input, options.crdout, options.seqfix)


def description():
    return '''This command takes a PDB file that has been run through the reduce program to determine missing
        hydrogen atoms in hystidines, and based on the added hydrogens by reduce, to identify HIS residues
        as any of the three histidine types in CHARMM, HSD, HSE or HSP. Runing reduce before this program
        is strongly recomended because hydrogen atoms added without running reduce are likely to be wrong.
        The input PDB file is outputed as a CRD structure file which is the type CHARMM takes as input.
        Prepare has an optional feature that outputs to files needed by CHARMM to finish setting up the
        structure. The files generated have the following suffixes .SEQ and _FIXERS.INP.
            --crdout = FILE      Path and name to output file with a 'crd' suffix. This file is needed for
                                 CHARMM to prepare a structure for simulation.
            --seqfix = String    This option takes only a 'Yes' or 'no', default is no. This file is needed
                                 for CHARMM to an amino acid sequence of a structure for simulations.'''


def usage():
    return 'pdb_cif.py prepare --input 1BRS.pdb --crdout 1BRS.crd --seqfix yes\n'


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
