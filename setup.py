#!/usr/bin/env python

"""
Created on Fri Jun 24 12:10:40 2016

@author: noel

"""

# References:
# https://gist.github.com/blackfalcon/8428401
# https://docs.python.org/2/distutils/setupscript.html
# https://docs.python.org/2/distutils/introduction.html#a-simple-example

from distutils.core import setup

setup(name='Entropy Maxima',
      version='0.1.0',
      description='Python Distribution for Protein Modeling and Design',
      author='Noel Carrascal, PhD',
      author_email='noelcarrascal@gmail.com',
      url='https://github.com/noelcjr/EntropyMaxima',
      packages=['em','em.describe','em.energy','em.manipulate','em.tools'],
      license='GNU GENERAL PUBLIC LICENSE',
      long_description='A program for the classification and manipulation of protein structures.',
      platforms='Tested on a Ubuntu linux machine with python 2.7.12. Should work with Windows or Macs but not tested.',
      scripts=['em/scripts/add_residues.py',
               'em/scripts/charmm_setup_water_to_ion_replacer.py','em/scripts/cif.py','em/scripts/del_residue.py',
               'em/scripts/flower.py','em/scripts/gen_csv.py','em/scripts/MMGBSA_CA.py','em/scripts/pdb.py'],
     )
