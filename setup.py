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
      version='0.1',
      description='Python Distribution for Protein Modeling and Design',
      author='Noel Carrascal, PhD',
      author_email='noelcarrascal@gmail.com',
      url='https://github.com/noelcjr/EntropyMaxima',
      packages=['em','em.describe','em.energy','em.manipulate','em.tools'],
     )
