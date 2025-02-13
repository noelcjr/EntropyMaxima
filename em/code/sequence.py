#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 10:36:56 2017

@author: noel
"""
import os
import sys
import Bio.AlignIO as AlignIO
import Bio.Seq as sq

class sequence(object):
    def __init__(self):
        self.chain_sequences = {}
        
    def add_Fasta_sequence(self,path):
        self.file_name = os.path.basename(path).split('.')[0]
        self.file_sufix = os.path.basename(path).split('.')[1]
        self.dir_path = os.path.dirname(path)
        self.chain_sequences = {}
        for alignment in AlignIO.parse(path, 'fasta', seq_count=1):
            for record in alignment:
                record_id = record.id.split()
                #record_id = record.id.split('|')
                #if ' '.join(record_id[1:]) == 'PDBID CHAIN SEQUENCE':
                self.chain_sequences[record_id[0]] = record.seq
                #else:
                #    print("ERROR: Formating of fasta files is not compatible with this method.")
                #    print("       It should be 'PDBID CHAIN SEQUENCE', and it is "+' '.join(record_id[1:]))
                #    sys.exit(1)
    
    def clear_sequences(self):
        self.chain_sequences = {}

    def add_sequence(self, key, sequence):
        self.chain_sequences[key] = sq.Seq(sequence)
