#!/usr/bin/env python

"""
Created on Fri Jun 24 12:10:40 2016

@author: noel

"""

# References:
# https://docs.python.org/2/distutils/setupscript.html
# https://docs.python.org/2/distutils/introduction.html#a-simple-example
# https://codeyarns.com/2014/05/30/how-to-install-and-uninstall-python-package-from-source/

import setuptools

setuptools.setup(
      name='Entropy Maxima',
      version='0.1.0',
      description='Python Distribution for Protein Modeling and Design',
      author='Noel Carrascal, PhD',
      author_email='noelcarrascal@gmail.com',
      url='https://github.com/noelcjr/EntropyMaxima',
      packages=['em','em.code','em.subscripts','em.charmm','em.charmm.gen'],
      package_data={'em': ['params/*.str','params/*.pdb','params/charmm27.ff/*']},
      license='GNU GENERAL PUBLIC LICENSE',
      long_description='A program for the classification and manipulation of protein structures.',
      platforms='Tested on a Ubuntu biolinux 16.04 LTS machine with python 3.7.3. Should work with Windows or Macs but not tested.',
      scripts=['em/scripts/charmm_setup_water_to_ion_replacer.py',
               'em/scripts/MMGBSA_CA.py',
               'em/scripts/flower.py','em/scripts/gen_csv.py','em/scripts/mmgbsa.py',
               'em/scripts/minimization.py'],
      install_requires=[
        'biopython',
        'pandas',
        'jinja2'
      ]
)
