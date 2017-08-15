#!/usr/bin/python

"""
Created on Mon Nov 21 08:53:06 2016

@author: noel
"""
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

file_dir = '/home/noel/Projects/Protein_design/Cas9/aln_pdb/cross_alignments/'
file_name = '4OO8_5CZZ.fasta.aln'
path = file_dir + file_name

#clustalw instead of clustalw2
cline = ClustalwCommandline("clustalw2", infile=path)
cline()

align = AlignIO.read(path, "clustal")

identical_count = 0
identical_aa = ''
for i in range(align.get_alignment_length()):
    if align.get_seq_by_num(0)[i] == align.get_seq_by_num(1)[i]:
        identical_aa += align.get_seq_by_num(0)[i]
        identical_count += 1

print identical_count

aa = 'ARNDCEQGHILKMFPSTWYV'
for char in aa:
  aa_count = identical_aa.count(char)
  if aa_count > 1:
    print char, aa_count

polar_charged_aa = 'RHKDESTNQ'
aa_count = 0
for char in polar_charged_aa:
    aa_count += identical_aa.count(char)
print('Polar or Charged AA:',aa_count)

non_polar_hydrophobic = 'CGPAVILMFYW'
aa_count = 0
for char in non_polar_hydrophobic:
    aa_count += identical_aa.count(char)
print('Non-Polar or Hydrophibic AA:',aa_count)
