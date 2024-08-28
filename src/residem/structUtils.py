import iotbx.pdb
import mmtbx
import mmtbx.model
import pandas as pd

""" pdb_hierarch.select(sel)
sel=pdb_hierarchy.atom_selection_cache().selection("protein")

"""


class ModelSelection:
    def __init__(self, pdb):
        self.pdb = pdb
        self.pdb_in = iotbx.pdb.input(pdb)
        self.model = mmtbx.model.manager(model_input=self.pdb_in)
        self.pdb_hierarchy = self.model.get_hierarchy()
        self.pdb_noh = self.non_h()
        self.atom_list = self.get_selection(self.pdb_noh)
        self.atom_df = pd.DataFrame(self.atom_list, columns=['chainid','resnum', 'resname', 'atomname', 'X', 'Y', 'Z'])

    def protein(self, copy=True):
        """ Select all the protein for pdb pdb_hierarchy"""
        cache = self.pdb_hierarchy.atom_selection_cache()
        sel = cache.selection('(not element H) and pepnames')
        return self.pdb_hierarchy.select(sel, copy_atoms=copy)

    def not_protein(self, copy=True):
        """ Select all the non protein for pdb pdb_hierarchy"""
        cache = self.pdb_hierarchy.atom_selection_cache()
        sel = cache.selection('(not element H) and (not pepnames)')
        return self.pdb_hierarchy.select(sel, copy_atoms=copy)

    def calphas(self, copy=True):
        """ Select all the C-alpha for pdb pdb_hierarchy"""
        cache = self.pdb_hierarchy.atom_selection_cache()
        sel = cache.selection('(not element H) and pepnames and (name CA)')
        return self.pdb_hierarchy.select(sel, copy_atoms=copy)

    def non_h(self, copy=True):
        """ Select all the non hydrogen for pdb_hierarchy"""
        cache = self.pdb_hierarchy.atom_selection_cache()
        sel = cache.selection('(not element H)')
        return self.pdb_hierarchy.select(sel, copy_atoms=copy)

    def backbone(self, copy=True):
        """ Select all the backbone for pdb_hierarchy"""
        cache = self.pdb_hierarchy.atom_selection_cache()
        sel = cache.selection('(not element H) and pepnames and (name C or name CA or name N or name O)')
        return self.pdb_hierarchy.select(sel, copy_atoms=copy)

    def sidechains(self, copy=True):
        """ Select all the side chains for pdb_hierarchy"""
        cache = self.pdb_hierarchy.atom_selection_cache()
        sel = cache.selection('(not element H) and pepnames and not (name C or name CA or name N or name O)')
        return self.pdb_hierarchy.select(sel, copy_atoms=copy)

    def non_water(self, copy=True):
        """ Select all the non water for pdb_hierarchy"""
        cache = self.pdb_hierarchy.atom_selection_cache()
        sel = cache.selection('(not resname HOH)')
        return self.pdb_hierarchy.select(sel, copy_atoms=copy)

    def water(self, copy=True):
        """ Select all the non water for pdb_hierarchy"""
        cache = self.pdb_hierarchy.atom_selection_cache()
        sel = cache.selection('(resname HOH)')
        return self.pdb_hierarchy.select(sel, copy_atoms=copy)

    def get_selection(self, pdb_hierarchy):
        """ pdb hierarchy selection"""
        atom_list = []
        item = pdb_hierarchy
        for m in item.models():
            for chain in m.chains():
                for residue_group in chain.residue_groups():
                    for conformers in residue_group.conformers():
                        for residue in conformers.residues():
                            for atoms in residue.atoms():
                                atom_list.append(
                                    ([chain.id,int(residue.resseq), residue.resname.strip(), atoms.name.strip(),
                                      atoms.xyz[0], atoms.xyz[1], atoms.xyz[2]]))
        return atom_list
