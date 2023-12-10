#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 21:44:52 2018

@author: noel
"""
import argparse
import em.code.MMGBSA_CA_L as MCL
import em.code.structure as strctr
import os

def register_parser(subparsers):
    parser = subparsers.add_parser('dhydration', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):  
    parser.add_argument("--crd", metavar="FILE", help="Coordinate file.", required=True)
    parser.add_argument("--gb0", metavar="FILE", help="gb file at z=0.", required=True)
    parser.add_argument("--gb500", metavar="FILE", help="gb file at z=500.", required=True)
    parser.set_defaults(func=run)

def run(options):
    analysis = MCL.hydration_shell_from_gb(os.getcwd()+'/'+options.crd,os.getcwd()+'/'+options.gb0,os.getcwd()+'/'+options.gb500)
    if not os.path.exists(options.crd):
        print("Error: Coordinate file not found.\n Type --help for description and options.")
    if not os.path.exists(options.gb0):
        print("Error: GB file of molecule in solution not found.\n Type --help for description and options.")
    if not os.path.exists(options.gb500):
        print("Error: GB file of molecule in solution not found.\n Type --help for description and options.")
    analysis.get_conserved_waters()

def description():
    return '''Generates de delta_MMGBSA between the gas and condensed state.
              Generates a file with delta (gas - hydr) MMGBSA CA of binding.
           '''
def usage():
    return '''\nmmgbsa.py dhydration'''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
    
dir_path = '/home/noel/Projects/Protein_design/Barnase_barstar/Structures/1B27/brs_A41C_A83C/structures/'
crd_file = strctr.crd(dir_path+'0.crd')
