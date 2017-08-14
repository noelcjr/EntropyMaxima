7# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 12:10:40 2016

@author: noel

"""
import os
import sys
import py_compile
# The paths below allows spyder to find .pyc, and the program to find param 
# files, It is now set for my environment only but for distribution it needs
# to be setup for each user installation or Errors will occur.
folder = 'EntropyMaxima/'
path = '/home/noel/Projects/Protein_design/'+folder
sys_path = path+'/src'
if not os.path.isdir(sys_path):
    print("Error: The parameter directory "+sys_path+" is not found.")
    sys.exit(1)

py_compile.compile(sys_path+'/CHARMM_Parser.py')
py_compile.compile(sys_path+'/Molecular_Descriptors.py')
py_compile.compile(sys_path+'/Molecular_Rigid_Manipulations.py')
py_compile.compile(sys_path+'/input_output.py')
py_compile.compile(sys_path+'/Super_Structures.py')
py_compile.compile(sys_path+'/utilities.py')
py_compile.compile(sys_path+'/Energy_Functions.py')
py_compile.compile(sys_path+'/MMGBSA_CA_L.py')
