#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import optparse

from em.subscripts import align, gaps, extract, fixpdb, maxmin, prepare, summary
subscripts = [align, gaps, extract, fixpdb, maxmin, prepare, summary]
usage = "usage: %prog --inputfile [input file] [options] [args]"
d = "This script is used to modify structures from PDB or CIF files as input."
opt_parser = optparse.OptionParser(usage, description=d)

arg_parser = argparse.ArgumentParser(description=d)
subparsers = arg_parser.add_subparsers()

for s in subscripts:
    s.register_parser(subparsers)

args = arg_parser.parse_args()
args.func(args)
