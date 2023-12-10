#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 21:26:58 2018

@author: noel
"""
import argparse
import em.code.MMGBSA_CA_L as MCL
import os

def register_parser(subparsers):
    parser = subparsers.add_parser('dbinding', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):  
    parser.add_argument("--self1", metavar="FILE", help="Path to MMGBSA self-electrostatics file.", required=True)
    parser.add_argument("--pair1", metavar="FILE", help="Path to MMGBSA pair electrostatics file.", required=True)
    parser.set_defaults(func=run)

def run(options):
    analysis = MCL.mmgbsa_ca_analysis(os.getcwd())
    if not os.path.exists(options.self1):
        print("Error: MMGBSA SELF file not found.\n Type --help for description and options.")
    if not os.path.exists(options.pair1):
        print("Error: MMGBSA PAIR file not found.\n Type --help for description and options.")
    analysis.mmgbsa_binding(options)

def description():
    return '''Generates de delta_MMGBSA between the bound and unbound state.
              Generates a file with delta (bound -unbound) MMGBSA CA of binding.
           '''
def usage():
    return '''\npdb_cif.py '''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
