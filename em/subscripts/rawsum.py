#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 22:09:20 2018

@author: noel
"""
import argparse
import em.code.MMGBSA_CA_L as MCL
import os

def register_parser(subparsers):
    parser = subparsers.add_parser('rawsum', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):  
    parser.add_argument("--self2", metavar="FILE", help="Path to MMGBSA self-electrostatics file.", required=True)
    parser.add_argument("--pair2", metavar="FILE", help="Path to MMGBSA pair electrostatics file.", required=True)
    parser.add_argument("--cut",type="float",help="A cutoff value in kCal/mols.")
    parser.add_argument("--range",type="str",default="(-10,10,200):(-5,5,300)",help="RANge is a string with TWO tuples sep\
arated by a column ':'. The tuples have three elements (min,max,bins) for the min and max values for the histogram, and\
 the number of bins for the histogram. Ex: '(-10,10,200):(-5,5,300)' is the default. This gives flexibility to the user\
 to visualize the components' energy distributions for self and pair wise energies and select and appropriate cutoff. I\
f the labels cover the histogram, play with min and max until the histogram shows properly. Four histograms are output.\
 The first two are always output and cannot be modified. They correspond to distributions within the full range of ener\
gies, and distributions for the components' energies within the selected cutoff. The other two are '(-10,10,200):(-5,5,\
300)' by default or any specified values by the user. This option is not required.")
    parser.set_defaults(func=run)

def run(options):
    analysis = MCL.mmgbsa_ca_analysis(os.getcwd())
    if not os.path.exists(options.self2):
        print("Error: MMGBSA SELF file not found.\n Type --help for description and options.")
    if not os.path.exists(options.pair2):
        print("Error: MMGBSA PAIR file not found.\n Type --help for description and options.")
    if not options.cut:
        print("Error: Cutoff for energies not specified.\n Type --help for description and options.")
    analysis.mmgbsa_raw_summary(options)

def description():
    return '''Initial Raw Sumamry Option","This histograms self and pair MMGBSA energie\
s as an initial overview of the values. It histograms contributions from GB, SA, EE, VDW.The default number of bins for\
 the histograms is that given by the max and min energies found divided by the cutoff. A cutoff is needed because most \
components' contributions to the energy are close to zero and can be discarted for further analysis. This summary will \
help find an appropriate cutoff. The cutoff is entered as a positive number, but it applies to both positive and negati\
ve energies, keeping energies grater thant a positive cutoff and lower thatn a negative cutoff.It generates other basic\
statistics. It also generates a plot of number of components vs. cutoff to guide the search for an optimal cutoff. The \
SUMMMARY IS OUTPUT TO THE TERMINAL, and to a TXT file.
Generates a raw summary from raw MMGBSA CA of binding.
           '''
def usage():
    return '''\npdb_cif.py '''

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)