#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 09:49:19 2018

@author: noel
"""
import sys
import re

class structural_expression(object):
    # Check Model is present in structure or exit with error
    message = """
    ERROR: You are trying to extract from an entity that is not present in 
    the structure. The structural Expression must refer to existing models.
    Program will EXIT because results might not be what you expect. 
    The entity not found would have the following values:"""
    def __init__(self,SE):
        self.structural_expression = SE
        self.sel_dic = self._SE_read_selections(self.structural_expression)

    def atom_list(self, strct):
        #######################################################################
        # The code here checks the entities requested are present. if a 
        # structural expression calls for an entity, or model, chain, residue 
        # or atom that is not present in the structure, the operation will not 
        # proceed. This is to force the user to be aware of what is present and
        # missing in the structure.                    
        # for i in range(count_entries):
        all_models = False
        if '*' in self.sel_dic['models']:
            if len(self.sel_dic['models']) == 1:
                all_models = True
        all_chains = False
        if '*' in self.sel_dic['chains']:
            if len(self.sel_dic['chains']) == 1:
                all_chains = True
        all_residues = False
        if '*' in self.sel_dic['residues']:
            if len(self.sel_dic['residues']) == 1:
                all_residues = True
        all_atoms = False
        if '*' in self.sel_dic['atoms']:
            if len(self.sel_dic['atoms']) == 1:
                all_atoms = True
        # Check that structures has at least one model
        if len(strct) == 0:
            print("ERROR: No models in the structure. Exit.")
            sys.exit(1)
        # Check models in structure
        if all_models:
            pass
        else:
            self._chck_mdls_n_strct(strct)
        # Check chains in models
        if all_chains:
            if all_models:
                pass
        else:
            if all_models:
                # check chain(s) is in all models
                for jm in strct.get_models():
                    self._chck_chns_n_mdls(jm)
            else:   
                for jm in strct.get_models():
                    if jm.get_id() in self.sel_dic['models']:
                        self._chck_chns_n_mdls(jm)
        # Check residues in all chains and models
        if all_residues:
            if all_models:
                if all_chains:
                    pass
        else:
            if all_models and all_chains:
                for jm in strct.get_models():
                    for r in self.sel_dic['residues']:
                        for kc in jm.get_chains():
                            if r not in kc:
                                print(self.message)
                                print("        Model  :",jm.get_id())
                                print("        Chain  :",kc.get_id())
                                print("        Residue:",r)
                                sys.exit(1)
            elif all_models and not all_chains:
                for jm in strct.get_models():
                    #self._chck_res_n_chns_mdls(jm,i,False)
                    for r in self.sel_dic['residues']:
                        for kc in self.sel_dic['chains']:
                            if (' ',r,' ') not in jm[kc]:
                                print(self.message)
                                print("        Model  :",jm.get_id())
                                print("        Chain  :",kc)
                                print("        Residue:",r)
                                sys.exit(1)
            elif not all_models and all_chains:
                for jm in self.sel_dic['models']:
                    #self._chck_res_n_chns_mdls(strct[jm],i,True)
                    for r in self.sel_dic['residues']:
                        for kc in strct[jm].get_chains():
                            if r not in kc:
                                print(self.message)
                                print("        Model  :",jm)
                                print("        Chain  :",kc.get_id())
                                print("        Residue:",r)
                                sys.exit(1)
            elif not all_models and not all_chains:
                for jm in self.sel_dic['models']:
                    #self._chck_res_n_chns_mdls(strct[jm],i,False)
                    for r in self.sel_dic['residues']:
                        for kc in self.sel_dic['chains']:
                            if (' ',r,' ') not in strct[jm][kc]:
                                print(self.message)
                                print("        Model  :",jm)
                                print("        Chain  :",kc)
                                print("        Residue:",r)
                                sys.exit(1)
        if all_atoms:
            if all_models:
                if all_chains:
                    if all_residues:
                        pass
        else:
            if all_models and all_chains and all_residues:
                for r in self.sel_dic['atoms']:
                    for rm in strct.get_models():
                        for rc in rm.get_chains():
                            for rr in rc.get_residues():
                                if r not in rr:
                                    print(self.message)
                                    print("        Model  :",rm.get_id())
                                    print("        Chain  :",rc.get_id())
                                    print("        Residue:",rr.get_id())
                                    print("           Atom:",r)
                                    sys.exit(1)
            elif all_models and all_chains and not all_residues:
                for r in self.sel_dic['atoms']:
                    for rm in strct.get_models():
                        for rc in rm.get_chains():
                            for rr in rc.get_residues():
                                if rr.get_id()[1] in self.sel_dic['residues']:
                                    if r not in rr:
                                        print(self.message)
                                        print("        Model  :",rm.get_id())
                                        print("        Chain  :",rc.get_id())
                                        print("        Residue:",rr.get_id())
                                        print("           Atom:",r)
                                        sys.exit(1)
            elif all_models and not all_chains and all_residues:
                for r in self.sel_dic['atoms']:
                    for rm in strct.get_models():
                        for rc in rm.get_chains():
                            if rc.get_id() in self.sel_dic['chains']:
                                for rr in rc.get_residues():
                                    if r not in rr:
                                        print(self.message)
                                        print("        Model  :",rm.get_id())
                                        print("        Chain  :",rc.get_id())
                                        print("        Residue:",rr.get_id())
                                        print("           Atom:",r)
                                        sys.exit(1)
            elif not all_models and all_chains and all_residues:
                for r in self.sel_dic['atoms']:
                    for rm in strct.get_models():
                        if rm.get_id() in self.sel_dic['models']:
                            for rc in rm.get_chains():
                                for rr in rc.get_residues():
                                    if r not in rr:
                                        print(self.message)
                                        print("        Model  :",rm.get_id())
                                        print("        Chain  :",rc.get_id())
                                        print("        Residue:",rr.get_id())
                                        print("           Atom:",r)
                                        sys.exit(1)
            elif all_models and not all_chains and not all_residues:
                for r in self.sel_dic['atoms']:
                    for rm in strct.get_models():
                        for rc in rm.get_chains():                                
                            if rc.get_id() in self.sel_dic['chains']:
                                for rr in rc.get_residues():
                                    if rr.get_id()[1] in self.sel_dic['residues']:
                                        if r not in rr:
                                            print(self.message)
                                            print("        Model  :",rm.get_id())
                                            print("        Chain  :",rc.get_id())
                                            print("        Residue:",rr.get_id())
                                            print("           Atom:",r)
                                            sys.exit(1)
            elif not all_models and all_chains and not all_residues:
                for r in self.sel_dic['atoms']:
                    for rm in strct.get_models():
                        if rm.get_id() in self.sel_dic['models']:
                            for rc in rm.get_chains():
                                for rr in rc.get_residues():
                                    if rr.get_id()[1] in self.sel_dic['residues']:
                                        if r not in rr:
                                            print(self.message)
                                            print("        Model  :",rm.get_id())
                                            print("        Chain  :",rc.get_id())
                                            print("        Residue:",rr.get_id())
                                            print("           Atom:",r)
                                            sys.exit(1)
            elif not all_models and not all_chains and all_residues:
                for r in self.sel_dic['atoms']:
                    for rm in strct.get_models():
                        if rm.get_id() in self.sel_dic['models']:
                            for rc in rm.get_chains():
                                if rc.get_id() in self.sel_dic['chains']:
                                    for rr in rc.get_residues():
                                        if r not in rr:
                                            print(self.message)
                                            print("        Model  :",rm.get_id())
                                            print("        Chain  :",rc.get_id())
                                            print("        Residue:",rr.get_id())
                                            print("           Atom:",r)
                                            sys.exit(1)
            elif not all_models and not all_chains and not all_residues:
                for r in self.sel_dic['atoms']:
                    for rm in strct.get_models():
                        if rm.get_id() in self.sel_dic['models']:
                            for rc in rm.get_chains():
                                if rc.get_id() in self.sel_dic['chains']:
                                    for rr in rc.get_residues():
                                        if rr.get_id()[1] in self.sel_dic['residues']:
                                            if r not in rr:
                                                print(self.message)
                                                print("        Model  :",rm.get_id())
                                                print("        Chain  :",rc.get_id())
                                                print("        Residue:",rr.get_id())
                                                print("           Atom:",r)
                                                sys.exit(1)
        # Now that the structural expression has pass some syntax test, it is
        # time to put all atoms location (not atom objects) in a list. This
        # is a complex loop made simpler with helper functions.
        all_srct_exp = []
        #for i in range(count_entries):
        tpl_sel = []
        for m in strct.get_models():
            if '*' in self.sel_dic['models'] and len(self.sel_dic['models']) == 1:
                self._pick_chains(m, tpl_sel)
            elif m.get_id() in self.sel_dic['models']:
                self._pick_chains(m, tpl_sel)
        all_srct_exp.append(tpl_sel)
        #######################################################################
        # should tuple be sorted to keep the atom order which might be changed 
        # when structural expresions have random order of selction?
        # if so, do something like this:
        # enumerated_selections = sorted(set(enumerated_selections), 
        #                               key=operator.itemgetter(1, 2, 3, 4))
        return all_srct_exp
    
    def _pick_chains(self, m, tpl_sel):
        for c in m.get_chains():
            if '*' in self.sel_dic['chains'] and len(self.sel_dic['chains']) == 1:
                self._pick_residues(c, tpl_sel)
            elif c.get_id() in self.sel_dic['chains']:
                self._pick_residues(c, tpl_sel)
    
    def _pick_residues(self, c, tpl_sel):
        for r in c.get_residues():
            if '*' in self.sel_dic['residues'] and len(self.sel_dic['residues']) == 1:
                self._pick_atoms(r, tpl_sel)
            elif r.get_id()[1] in self.sel_dic['residues']:
                self._pick_atoms(r, tpl_sel)
    
    def _pick_atoms(self, r, tpl_sel):
        chain = r.get_parent().get_id()
        model = r.get_parent().get_parent().get_id()
        for a in r.get_atoms():
            if '*' in self.sel_dic['atoms'] and len(self.sel_dic['atoms']) == 1:
                tpl_sel.append((model,chain,r.get_id()[1],a))
            elif a.get_id() in self.sel_dic['atoms']:
                tpl_sel.append((model,chain,r.get_id()[1],a))
    
    def _chck_atms_n_res_chns_mdls(self, strct, i):
        for r in self.sel_dic['atoms'][i]:
            for rm in strct.get_models():
                for rc in rm.get_chains():
                    for rr in rc.get_residues():
                        if r not in rr:
                            print(self.message)
                            print("        Model  :",rm.get_id())
                            print("        Chain  :",rc.get_id())
                            print("        Residue:",rr.get_id())
                            print("           Atom:",r)
                            sys.exit(1) 
    
    def _chck_res_n_chns_mdls(self,jm,i,all_chns):
        for r in self.sel_dic['residues'][i]:
            if all_chns:
                for kc in jm.get_chains():
                    if r not in kc:
                        print(self.message)
                        print("        Model  :",jm.get_id())
                        print("        Chain  :",kc.get_id())
                        print("        Residue:",r)
                        sys.exit(1)
            else:
                for kc in self.sel_dic['chains'][i]:
                    if r not in kc:
                        print(self.message)
                        print("        Model  :",jm.get_id())
                        print("        Chain  :",kc)
                        print("        Residue:",r)
                        sys.exit(1)

    def _chck_chns_n_mdls(self,jm):
        for c in self.sel_dic['chains']: 
            if c not in jm:
                print(self.message)
                print("        Model  :",jm.get_id())
                print("        Chain  :",c)
                sys.exit(0)
    
    def _chck_mdls_n_strct(self, strct):
        for m in self.sel_dic['models']:
            if m not in strct:
                print(self.message)
                print("        Model:",m)
                sys.exit(1)

    def _SE_read_selections(self, selection):
        section = {'models':[],'chains':[],'residues':[],'atoms':[]}
        is_section_present = {'models':False,'chains':False,
                              'residues':False,'atoms':False}
        b = re.split('(m\[|c\[|r\[|a\[|\])',selection)
        # To removes empty strings from list generated by re.split:
        c = list(filter(None, b))
        # Now, we want to restrict only one m,c,r and a for each structural
        # expression. More than one or zero entries will exit the program
        for d in ['m[','c[','r[','a[']:
            count = 0
            for j in c:
                if j == d:
                    count += 1
            if count == 0:
                print("ERROR: "+d+" is missing "+str(count))
                print("Program will exit now.")
                sys.exit(1)
            if count > 1:
                print("ERROR: "+d+" is present more than once "+str(count))
                print("Program will exit now.")
                sys.exit(1)
        # The arguments for the structural expression must be in a strict 
        # format: m[0]c[A,B]r[10:13]a[CA] which when split must give a list
        # with 12 elements. Missign a m[, c[, r[, a[ or ] would cause an 
        # error. It does not matter what is inside as long they are not
        # characters used by split and are digits or integers alone or 
        # separated by ',' or ':'. The format of this list is rigid to 
        # avoid confusions, and it expects to find certain elements in 
        # certain list positions or it will give an error.
        if c[0].lower() == 'm[':
            if c[1] == ']' or c[1] == 'c[':
                print("ERROR: Structural expresion format is wrong. \
                      Model (m) is missing a numerical entry")
                sys.exit(1)
            section['models'] = self._SE_expand_selections('model',c[1],'int')
            is_section_present['models'] = True
        else:
            print("ERROR: Wrong structural expression format. \
                  Model (m) us missing a numerical entry.")
            sys.exit(1)
        if c[3].lower() == 'c[':
            if c[4] == ']' or c[4] == 'r[':
                print("ERROR: Structural expresion format is wrong. Chain \
                      (c) is missing an alphabetical entry.")
                sys.exit(1)
            section['chains'] = self._SE_expand_selections('chain',c[4],'str')
            is_section_present['chains'] = True
        else:
            print("ERROR: Wrong structural expression format. \
                  Chain (c) needs a valid entry.")
            sys.exit(1)
        if c[6].lower() == 'r[':
            if c[7] == ']' or c[7] == 'a[':
                print("ERROR: Structural expresion format is wrong. \
                      Residue (r) is missing a numerical entry.")
                sys.exit(1)
            section['residues'] = self._SE_expand_selections('residue',c[7],'int')
            is_section_present['residues'] = True
        else:
            print("ERROR: Wrong structural expression format.\
                   Residue (r) has an invalid entry.")
            sys.exit(1)
        if c[9].lower() == 'a[':
            if c[10] == ']' or c[10] == 'a[':
                print("ERROR: Structural expresion format is wrong. \
                      Atom (a) is missing an alphabetical entry.")
                sys.exit(1)
            section['atoms'] = self._SE_expand_selections('atom',c[10],'str')
            is_section_present['atoms'] = True
        else:
            print("ERROR: Wrong structural expression format at 9.\
                  Atom (a) is mising an alphabetical entry.")
            sys.exit(1)
        # The next loop makes sure that a structural expression is not
        # missing m,c,r or a sections            
        for k in is_section_present:
            if not is_section_present[k]:
                print("ERROR: Estructural expressions need to specify \
                      models (m), chains (c),residues (r) and atoms (a).\
                      If any of these is missing the expression is wrong. \
                      Example: m[0]c[A,B]r[10:13]a[CA]")
                sys.exit(1)
        return section
    
    def _SE_expand_selections(self, label, value, value_type):
        """ Takes an argument inside a structural expression bracket and
            expands it taking into consideration if it is int or string and
            columns and commas used in the structural expression synthax
        """
        results = []
        if value_type == 'int':
            for i in value.split(','):
                columns = i.split(":")
                if len(columns) == 1:
                    if self.RepresentsInt(columns[0]):
                        results.append(int(columns[0]))
                    else:
                        if columns[0] == '*':
                            results.append('*')
                        else:
                            print("Error: wrong character "+columns[0]+
                                  " in "+label+" only integers allowed. Exit.")
                            sys.exit(1)
                elif len(columns) == 2:
                    for k in columns:
                        if self.RepresentsInt(k):
                            pass
                        else:
                            print('ERROR: A '+label+' requires only integers \
                                  inside the brakets and/or separated by \
                                  columns. (2)\
                                  Program will exit without running.')
                            sys.exit(1)
                    if int(columns[0]) >= int(columns[1]):
                        print("ERROR: A range of numbers defined by two va-\
                              riables separated by a column (i.e. 2:5) requi-\
                              res the value before the column to be less than \
                              the one after the column. You entered:\n    "+
                              str(columns[0])+":"+str(columns[1])+"\n"+
                              "Program will exit without running.")
                        sys.exit(1)
                    if isinstance(int(columns[0]), int) and isinstance(
                            int(columns[1]), int):
                        for j in range(int(columns[0]),int(columns[1])+1):
                            results.append(j)
                    else:
                        print("Error: only integer values allowed in ranges.\n\
                              Ex: 2:4")
                        sys.exit(1)
                else:
                    print("Error: Wrong Integer value entered in structural \
                          expression. See instructions for valid types.")
                    sys.exit(1)
        elif value_type == 'str':
            for i in value.split(','):
                if i == '*':
                    results.append(i)
                else:
                    results.append(i)
        return results
    
    def RepresentsInt(self, s):
        try: 
            int(s)
            return True
        except ValueError:
            return False
