"""
Script that generates the necessary calculations for collaboration
with John Gregoire.  Specifically, multiphase orderings of the
Mn-Sb chemical system.
"""

from pymatgen import MPRester
from pymatgen.transformations.advanced_transformations import EnumerateStructureTransformation 
import tqdm
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
        structure = template.copy()
        structure.replace_species({"Mn": composition, "Sb": composition})
        est = EnumerateStructureTransformation(max_cell_size=8)
        structures.append(est.apply_transformation(structure))
    return structures

def get_all_structs_by_material_id(mp_ids):
    return {mpid: get_orderings(mpr.get_structure_by_material_id(mpid))
            for mpid in tqdm.tqdm(mp_ids)}

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
