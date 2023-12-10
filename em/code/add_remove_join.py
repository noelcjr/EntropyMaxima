#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 21:05:41 2019

@author: noel
"""
import copy
import sys
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
import Bio.PDB as struct
import em.code.utilities as ut
import Bio.PDB.Superimposer as SuperImposer
import em.code.structural_expression as SE
#import pkg_resources

class add_residues(object):
    def __init__(self, strct_path):
        """ Modifications are done on a single chain for all models of a 
        structure. Modifyin different models in a structure separately would
        mean that a PDB file would have models of different composition, and it
        should not be in a single PDB structure file but in separate ones.
        Models should vary in conformation only, not composition.
        """
        # TODO uncomment the following line when installing as a package
        #      also uncomment import pkg_resources above
        #self.pep_f = pkg_resources.resource_filename('em',
        #                                             'params/peptides.pdb')
        self.LS = []
        self.mdls = {}
        self.pdb_parser = PDBParser(QUIET = True)
        self.strct = self.pdb_parser.get_structure('main', strct_path)
        self.min_max_resid_name()
        # TODO this file must be set at installation because the path is system dependent
        self.pep_path = "/home/noel/Code/EntropyMaxima/em/params/peptides.pdb"
        self.s_pep = self.pdb_parser.get_structure('peptides', self.pep_path)

    def min_max_resid_name(self):
        """ This generates min and max information for all chains in all models.
        The min/max info is used to renumber residues to chains that have been 
        extended at either of the C or N terminal. It gives the answer for all
        models even though the anser should be the same for all models.
        """
        for g in self.strct.get_models():
            chns = {}
            for h in g.get_chains():
                first_res_id = True
                resd = {}
                for i in h.get_residues():
                    if first_res_id:
                        resd['min_id'] = i.get_id()
                        resd['min_resname'] = i.resname
                        resd['max_id'] = i.get_id()
                        resd['max_resname'] = i.resname
                        first_res_id = False
                    if i.get_id()[1] > resd['max_id'][1]:
                        resd['max_id'] = i.get_id()
                        resd['max_resname'] = i.resname
                    if i.get_id()[1] < resd['min_id'][1]:
                        resd['min_id'] = i.get_id()
                        resd['min_resname'] = i.resname
                chns[h.get_id()] = resd
            self.mdls[g.get_id()] = chns
            
    def init_new_structure(self, chn_id):
        """ Creates a new structure where new residues will be placed to a
        chain at the time.
        """
        self.new_s = struct.StructureBuilder.Structure("AddLink")
        model_number = -1
        for i in self.strct.get_models():
            model_number = i.get_id()
            new_m = struct.Model.Model(model_number)
            for h in i.get_chains():
                if h.get_id() == chn_id:
                    new_c = struct.Chain.Chain(chn_id)
                    new_m.add(new_c)
                else:
                    new_m.add(h)
            self.new_s.add(new_m)
            
    def gen_linker_strct_list(self, chn_id, linker):
        """ For every residue added create a structure and place it in a list
        this is because alignment of atoms works with structures. The previous 
        and next residue from peptide.pdb is added to each structure in order 
        to have atoms for proper geometrical alignment. peptide.pdb has three 
        GLY, two at the ends, to have atoms for alignment for the residues next
        to them.
        """
        lnk_c = 0
        res_num = 1
        for h in linker:
            tmp_s = struct.StructureBuilder.Structure("TempStct "+str(lnk_c))
            for i in self.strct.get_models():
                tmp_m = struct.Model.Model(i.get_id())
                for j in i.get_chains():
                    if j.get_id() == chn_id:
                        tmp_c = struct.Chain.Chain(chn_id)
                        first_res = True
                        second_gly = False
                        next_res = False
                        for k in self.s_pep.get_residues():
                            if first_res:
                                prev_r = copy.deepcopy(k)
                                prev_r.detach_parent()
                                prev_r.id = (' ', 1, ' ')
                                prev_r.resname = k.resname
                                prev_r.segid = k.segid
                                first_res = False
                            elif next_res:
                                nxt_r = copy.deepcopy(k)
                                nxt_r.detach_parent()
                                nxt_r.id = (' ', 3, ' ')
                                nxt_r.resname = k.resname
                                nxt_r.segid = k.segid
                                tmp_c.add(nxt_r)
                                next_res = False
                            else:
                                if k.resname == ut.residueDict1_1[h]:
                                    # This bypasses third GLY in peptides.pdb
                                    if k.resname == "GLY" and second_gly:
                                        pass
                                    else:
                                        temp = copy.deepcopy(k)
                                        temp.detach_parent()
                                        temp.id = (' ', 2, ' ')
                                        temp.resname = k.resname
                                        temp.segid = k.segid
                                        tmp_c.add(prev_r)
                                        tmp_c.add(temp)
                                        next_res = True
                                        if k.resname == "GLY":
                                            second_gly = True
                                else:
                                    prev_r = copy.deepcopy(k)
                                    prev_r.detach_parent()
                                    prev_r.id = (' ', 1, ' ')
                                    prev_r.resname = k.resname
                                    prev_r.segid = k.segid
                        tmp_m.add(tmp_c)
                tmp_s.add(tmp_m)
            self.LS.append(tmp_s)
            lnk_c += 1
            res_num += 1
    
    # TODO delete after testing use the check_for_position
    def check_for_gaps_CTER(self, chn_id, position):
        position_exist = False
        get_next_res = False
        next_res_num = -1
        for i in self.strct.get_models():
            for j in i.get_chains():
                if j.get_id() == chn_id:
                    for k in j.get_residues():
                        if get_next_res:
                            next_res_num = k.get_id()[1]
                            get_next_res = False
                        if k.get_id()[1] == position:
                            position_exist = True
                            get_next_res = True
        return position_exist, next_res_num             

    def add_CTER_dir(self, chn_id, linker, position):
        self.init_new_structure(chn_id)
        self.gen_linker_strct_list(chn_id, linker)
        len_LS = len(self.LS)
        # Check that position exists, and if so, return gap number 0 -> infinity
        gap_in_pos = self.check_for_gaps_CTER(chn_id, position)
        # if position does not exists, exit. 
        if not gap_in_pos[0]:
            print("ERROR: Trying to attach new residues to a residue that ")
            print("       is not present is impossible. Exit with error.")
            sys.exit(1)
        else:
            fixed = "m[*]c["+chn_id+"]r["+str(position)+"]a[N,CA,C]"
            SE_F_M = SE.structural_expression(fixed)
            fixed_atoms = [i[3] for i in SE_F_M.atom_list(self.strct)[0]]
            
            moved = "m[*]c["+chn_id+"]r[1]a[N,CA,C]"
            SE_M_A = SE.structural_expression(moved)
            moved_atoms = [i[3] for i in SE_M_A.atom_list(self.LS[0])[0]]
            super_imposer = SuperImposer()
            super_imposer.set_atoms(fixed_atoms, moved_atoms)
            super_imposer.apply(self.LS[0].get_atoms())
            
            fixed = "m[*]c["+chn_id+"]r[2]a[N,CA,C]"
            for i in range(len_LS-1):
                SE_F = SE.structural_expression(fixed)
                fixed_atoms = [j[3] for j in SE_F.atom_list(self.LS[i])[0]]
                SE_M = SE.structural_expression(moved)
                moved_atoms = [j[3] for j in SE_M.atom_list(self.LS[i+1])[0]]
                super_imposer.set_atoms(fixed_atoms, moved_atoms)
                super_imposer.apply(self.LS[i+1].get_atoms())

            for i in self.new_s.get_models():
                i_mod = i.get_id()
                count = 1
                for j in i.get_chains():
                    if j.get_id() == chn_id:
                        for x in self.strct.get_models():
                            residues_per_model = []
                            residues_per_model_inx = []
                            for y in x.get_chains():
                                if y.get_id() == chn_id:
                                    id_count = self.mdls[x.get_id()][y.get_id()]['max_id'][1] + len(linker) + gap_in_pos[1] + 1
                                    gap_size = gap_in_pos[1] - position - 1
                                    count2 = id_count
                                    track_resids = []
                                    for z in y.get_residues():
                                        if z.get_id()[1] <= position:
                                            temp = copy.deepcopy(z)
                                            temp.detach_parent()
                                            self.new_s[i_mod][chn_id].add(temp)
                                            count = temp.get_id()[1]
                                            track_resids.append(count)
                                    count += 1       
                                    for z in range(len_LS):
                                        add_res = self.LS[z][i_mod][chn_id][2]
                                        add_res.detach_parent()
                                        add_res.id = (' ',count,' ')
                                        self.new_s[i_mod][chn_id].add(add_res)
                                        track_resids.append(count)
                                        count += 1 
                                    for z in y.get_residues():
                                        if z.get_id()[1] > position:
                                            temp = copy.deepcopy(z)
                                            temp.detach_parent()
                                            residues_per_model_inx.append(temp.get_id()[1])
                                            temp.id = (' ', count2, ' ')
                                            residues_per_model.append(temp)
                                            count2 += 1

                                    count_indx = 0
                                    for z in residues_per_model:
                                        self.new_s[i_mod][chn_id].add(z)
                                        temp = z.get_id()
                                        # In different situations, displaced residue numbers get
                                        # different new_indexes or numbers. The following if 
                                        # statements deal only with situations I could imagine
                                        if gap_size < len(linker):
                                            new_index = residues_per_model_inx[count_indx]+(len(linker)-gap_size)
                                        elif gap_size > len(linker):
                                            new_index = residues_per_model_inx[count_indx]
                                        else:
                                            new_index = residues_per_model_inx[count_indx]
                                        self.new_s[i_mod][chn_id][temp[1]].id = (temp[0], 
                                                  new_index, temp[2])
                                        count_indx += 1

    def check_for_gaps_NTER(self, chn_id, position):
        position_exist = False
        prev_residue = 0
        for i in self.strct.get_models():
            for j in i.get_chains():
                if j.get_id() == chn_id:
                    for k in j.get_residues():
                        if k.get_id()[1] == position:
                            position_exist = True
                        elif k.get_id()[1] < position:
                            prev_residue = k.get_id()[1]           
        return position_exist, prev_residue
    
    def add_NTER_dir(self, chn_id, linker, position):
        self.init_new_structure(chn_id)
        self.gen_linker_strct_list(chn_id, linker)
        len_LS = len(self.LS)
        # Check that position exists, and if so, return gap number 0 -> infinity
        gap_in_pos = self.check_for_gaps_NTER(chn_id, position)
        # if position does not exists, exit. 
        if not gap_in_pos[0]:
            print("ERROR: Trying to attach new residues to a residue that ")
            print("       is not present is impossible. Exit with error.")
            sys.exit(1)
        else:
            fixed = "m[*]c["+chn_id+"]r["+str(position)+"]a[N,CA,C]"
            SE_F_M = SE.structural_expression(fixed)
            fixed_atoms = [i[3] for i in SE_F_M.atom_list(self.strct)[0]]
            
            moved = "m[*]c["+chn_id+"]r[3]a[N,CA,C]"
            SE_M_A = SE.structural_expression(moved)
            moved_atoms = [i[3] for i in SE_M_A.atom_list(self.LS[-1])[0]]
            super_imposer = SuperImposer()
            super_imposer.set_atoms(fixed_atoms, moved_atoms)
            super_imposer.apply(self.LS[-1].get_atoms())
            
            fixed = "m[*]c["+chn_id+"]r[2]a[N,CA,C]"
            for i in range(len_LS-1,0,-1):
                SE_F = SE.structural_expression(fixed)
                fixed_atoms = [j[3] for j in SE_F.atom_list(self.LS[i])[0]]
                SE_M = SE.structural_expression(moved)
                moved_atoms = [j[3] for j in SE_M.atom_list(self.LS[i-1])[0]]
                super_imposer.set_atoms(fixed_atoms, moved_atoms)
                super_imposer.apply(self.LS[i-1].get_atoms())
            # This relabels residues in the chain that is going to be extended
            # so that biopython allows insertions in the Nterm
            # min_max = self.min_max_resid_name()
            gap = position - gap_in_pos[1] - 1
            if gap == len(linker):
                print("gap == len(linker)")
                offset = 0
                linker_index = gap_in_pos[1] + 1
            elif gap < len(linker):
                print("gap < len(linker)")
                offset = 0
                linker_index = gap_in_pos[1] + 1
            elif gap > len(linker):
                print("gap > len(linker)")
                offset = len(linker)- gap_in_pos[1]
                linker_index = position - len(linker)
            print(gap,position,gap_in_pos[0],gap_in_pos[1],offset,linker_index)  
            for j in self.strct.get_models():
                temp_idx = self.mdls[j.get_id()][chn_id]['max_id'][1] + len(linker) + 1
                residues_per_model_before_p = []
                residues_per_model_after_p = []
                residues_per_model_before_p_idx = []
                residues_per_model_after_p_idx = []
                for k in j.get_chains():
                    if k.get_id() == chn_id:
                        for l in k.get_residues():
                            if l.get_id()[1] < position:
                                temp = copy.deepcopy(l)
                                temp.detach_parent()
                                if gap == len(linker):
                                    residues_per_model_before_p.append(temp)
                                    residues_per_model_before_p_idx.append(temp.get_id()[1])
                                elif gap < len(linker):
                                    residues_per_model_before_p.append(temp)
                                    residues_per_model_before_p_idx.append(temp.get_id()[1])
                                elif gap > len(linker):
                                    residues_per_model_before_p.append(temp)
                                    residues_per_model_before_p_idx.append(temp.get_id()[1])
                                else:
                                    pass
                            else:
                                temp = copy.deepcopy(l)
                                temp.detach_parent()
                                if gap == len(linker):
                                    residues_per_model_after_p.append(temp)
                                    residues_per_model_after_p_idx.append(temp.get_id()[1])
                                elif gap < len(linker):
                                    residues_per_model_after_p.append(temp)
                                    residues_per_model_after_p_idx.append(temp.get_id()[1]+temp_idx)
                                elif gap > len(linker):
                                    residues_per_model_after_p.append(temp)
                                    residues_per_model_after_p_idx.append(temp.get_id()[1])
                                else:
                                    pass
                                
            for i in self.new_s.get_models():
                i_mod = i.get_id()
                for j in i.get_chains():
                    if j.get_id() == chn_id:        
                        count = self.mdls[i.get_id()][chn_id]['min_id']
                        count = 0
                        for k in residues_per_model_before_p:
                            if gap == len(linker):
                                self.new_s[i.get_id()][chn_id].add(k)
                            elif gap < len(linker):
                                self.new_s[i.get_id()][chn_id].add(k)
                            elif gap > len(linker):
                                self.new_s[i.get_id()][chn_id].add(k)
                            else:
                                self.new_s[i.get_id()][chn_id].add(k)
                                res_num = residues_per_model_before_p_idx[count]
                                temp = self.new_s[i.get_id()][chn_id][res_num].get_id()[1] + offset
                                self.new_s[i_mod][chn_id][count].id = (temp[0], count, temp[2])
                            count += 1
                        for k in range(len_LS):
                            add_res = self.LS[k][i_mod][chn_id][2]
                            add_res.detach_parent()
                            add_res.id = (' ',linker_index,' ')
                            self.new_s[i_mod][chn_id].add(add_res)
                            linker_index += 1
                        count = 0
                        for k in residues_per_model_after_p:
                            if gap == len(linker):
                                self.new_s[i.get_id()][chn_id].add(k)
                            elif gap < len(linker):
                                res_num = residues_per_model_after_p_idx[count]-temp_idx+len(linker)-gap
                                k.id = (' ',res_num,' ')
                                self.new_s[i.get_id()][chn_id].add(k)
                            elif gap > len(linker):
                                self.new_s[i.get_id()][chn_id].add(k)
                            else:
                                self.new_s[i.get_id()][chn_id].add(k)
                                res_num = residues_per_model_after_p_idx[count]
                                temp = self.new_s[i.get_id()][chn_id][res_num].get_id()[1]
                                self.new_s[i_mod][chn_id][count].id = (temp[0], temp, temp[2])
                            count += 1
            
    def view_chn_modified(self):
        for i in self.new_s.get_models():
            for j in i.get_chains():
                for k in j.get_residues():
                    print(i,j,k,k.get_id())
                        
    def output_structure(self,out_filepath):
        io = PDBIO()
        io.set_structure(self.new_s)
        io.save(out_filepath)
        
    def remove_from_se(self, remove):
        rSE = SE.structural_expression(remove)
        remove_atoms = rSE.atom_list(self.strct)
        
        keys_rm = {}
        for h in remove_atoms[0]:
            if (h[0],h[1],h[2]) in keys_rm:
                keys_rm[(h[0],h[1],h[2])].append(h[3].get_id())
            else:
                keys_rm[(h[0],h[1],h[2])] = []
                keys_rm[(h[0],h[1],h[2])].append(h[3].get_id())
            
        for h in keys_rm:
            for i in self.strct.get_models():
                if i.get_id() == h[0]:
                    for j in i.get_chains():
                        if j.get_id() == h[1]:
                            for k in j.get_residues():
                                if k.get_id()[1] == h[2]:
                                    for l in keys_rm[(h[0], h[1], h[2])]:
                                        k.detach_child(l)
        # Update mdls
        self.min_max_resid_name()

