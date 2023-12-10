#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 21:51:31 2018

@author: noel
"""
import argparse
import optparse

from em.subscripts import dbinding, dhydration
#, dcabinding, dcahydration, rawsum
subscripts = [dbinding, dhydration]
#, dcabinding, dcahydration, rawsum]
usage = "usage: %prog --inputfile [input file] [options] [args]"
d = "This script is used to modify structures from PDB or CIF files as input."
opt_parser = optparse.OptionParser(usage, description=d)

arg_parser = argparse.ArgumentParser(description=d)
subparsers = arg_parser.add_subparsers()

for s in subscripts:
    s.register_parser(subparsers)

args = arg_parser.parse_args()
args.func(args)