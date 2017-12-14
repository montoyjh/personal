"""
Script that generates the necessary calculations for collaboration
with John Gregoire.  Specifically, multiphase orderings of the
Mn-Sb chemical system.
"""

from pymatgen import MPRester
from pymatgen.transformations.advanced_transformations import EnumerateStructureTransformation 
from pymatgen.analysis.structure_matcher import StructureMatcher
import tqdm
import itertools
import numpy as np

mpr = MPRester()

def get_orderings(template, compositions=None):
    """
    Gets ordered structures from a template and a substitution

    Args:
        Composition
    """
    if not compositions: 
        red = template.composition.reduced_composition.get_el_amt_dict()
        total = red['Mn'] + red['Sb']
        step = 1 / (8 - 8 % total)
        compositions = [{"Mn": x, "Sb": 1-x} for x in  np.arange(0, 1.01, step)]

    structures = []
    for composition in compositions:
        structures.extend(order_structure(template, composition))
        """
        structure = template.copy()
        structure.replace_species({"Mn": composition, "Sb": composition})
        if structure.is_ordered:
            structures.append(structure)
        else:
            # import pdb; pdb.set_trace()
            est = EnumerateStructureTransformation(refine_structure=True, max_cell_size=8)
            structures.extend(est.apply_transformation(structure))
        """
    return structures

def get_all_structs_by_material_id(mp_ids):
    return {mpid: get_orderings(mpr.get_structure_by_material_id(mpid))
            for mpid in tqdm.tqdm(mp_ids)}


def order_structure(structure, composition):
    """
    Enumlib doesn't seem to work, so maybe just
    do it in a standard way?
    """
    # find minimum supercell for at least 8 sites
    this_structure = structure.copy()
    elamt = structure.composition.get_el_amt_dict()
    total = elamt['Mn'] + elamt['Sb']
    if total < 8:
        this_structure.make_supercell([2]*3)
        total *= 8
    indices = [structure.index(site) for site in structure 
               if site.species_string in ['Mn', 'Sb']]
    combos = list(itertools.combinations(indices, int(len(indices) * composition['Mn'])))
    unique_structures = []
    structure_matcher = StructureMatcher()
    name = '{} - {}'.format(structure.formula, composition)
    unique_structs = []
    # import pdb; pdb.set_trace()
    for combo in tqdm.tqdm(combos, desc='enumeration of {}'.format(name)):
        new_struct = this_structure.copy()
        mn_sites = combo
        sb_sites = list(set(indices) - set(combo)) 
        for si in mn_sites:
            new_struct[si] = 'Mn'
        for si in sb_sites:
            new_struct[si] = 'Sb'
        # Check if in unique list or not
        if not any([structure_matcher.fit(new_struct, struct) for struct in unique_structs]):
            unique_structs.append(new_struct)
    return unique_structs

if __name__=="__main__":
    structures = get_all_structs_by_material_id(["mp-19231", 
                                                 "mp-25043", 
                                                 "mp-565203",
                                                 "mp-18759",
                                                 "mp-19006",
                                                 "mp-19395",
                                                 "mp-2136",
                                                 "mp-1705",
                                                 "mp-230"])
    print("\n".join("{}: {}".format(k, len(v)) for k, v in structures.items()))
    import pdb; pdb.set_trace()
