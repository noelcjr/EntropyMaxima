#!/usr/bin/python

"""
Created on Tue Oct 25 15:19:00 2016

@author: noel
"""
import os
import sys
import em.charmm.gen.inputs as IO
import optparse

# #bash setup_flowers_3.sh DDDDK_DDDDK 2hiu_a_ddddk_b_ddddk_1rr.pdb 1brs_comp_1rr_1rr.pdb "A:F,B:C" 45 5 45 5

def main():
    usage = "usage: %prog [options] arg"
    d = "This program does a minimization on an PDB file, and it takes the outputs and places them in a single CSV file"
    opt_parser = optparse.OptionParser(usage,description=d)

    opt_parser.add_option("--input", type="str", help="Inpput PDB file.")
    opt_parser.add_option("--csv", type="str", help="A path to a csv file where minimization results are stored.")
    opt_parser.add_option("--min", type="str", help="Enter a minimization template type.")

    options, args = opt_parser.parse_args()
############################################  Options Entered ##########################################################

if __name__ == '__main__':
    main()
