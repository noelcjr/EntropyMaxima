# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 22:54:08 2016

@author: noel
"""

class utilities(object):
    residueDict1_1 = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS',\
                      'E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE',\
                      'L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO',\
                      'S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
    residueDict1_2 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C', \
                      'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I', \
                      'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P', \
                      'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
    def __init__(self):
        pass   
    def check_gaps(self,data):
        ''' The following code identify missing indexs in an monotonically 
        increasing array. Full sequence is obtained from begining of sequence
        to the last amino acid in the sequence. Then we look for missing amino 
        acid indexes in the available structural information. 
        '''
        #data=[3,4,5,8,9]
        count = 0
        inserts = []
        for i in range(1,data[-1]+1):
            if i not in data:
                inserts.append((count,i))
            count += 1
        data2 = [i for i in data]
        for i in inserts:
            data2.insert(i[0],i[1])
        return inserts
    def gap_report(self,inserts):
        count = 0
        init = True
        first = 0
        last = 0
        report = []
        for i in inserts:
            if init:
                first = i[1]
                last = i[1]
                init = False
                count += 1
            if i[1] == last+1:
                last = i[1]
            else:
                if first != last:
                    report.append((count,first,last))
                    count += 1
                first = i[1]
                last = i[1]
        report.append((count,first,last))
        return report
    def count_gaps(self,present):
        count = 0
        for i in range(len(present)-1):
            if (present[i]+1) != present[i+1]:
                count += 1
        return count
