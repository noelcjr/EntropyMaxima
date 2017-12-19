#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 10:36:56 2017

@author: noel
"""
import os
import sys
import Bio.AlignIO as AlignIO

class sequence(object):
    def __init__(self,path, seq_id='', file_type='fasta'):
        if file_type == 'str':
            self.chain_sequences = {}
            self.chain_sequences[seq_id] = path
        elif file_type == 'fasta':
            self.file_name = os.path.basename(path).split('.')[0]
            self.file_sufix = os.path.basename(path).split('.')[1]
            self.dir_path = os.path.dirname(path)
            self.chain_sequences = {}
            for alignment in AlignIO.parse(path, 'fasta', seq_count=1):
                for record in alignment:
                    record_id = record.id.split('|')
                    if ' '.join(record_id[1:]) == 'PDBID CHAIN SEQUENCE':
                        self.chain_sequences[record_id[0]] = record.seq
                    else:
                        print("ERROR: Formating of fasta files is not compatible with this method.")
                        print("       It should be 'PDBID CHAIN SEQUENCE', and it is "+' '.join(record_id[1:]))
                        sys.exit(1)
        else:
            pass
    
    def add_sequence(self):
        pass