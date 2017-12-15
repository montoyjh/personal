"""
Script that generates the necessary calculations for collaboration
with John Gregoire.  Specifically, multiphase orderings of the
Mn-Sb chemical system.

Notes:
    - originally tried to make this work for every intermediate composition
        which got crazy fast, so limited it to substitutions of single atoms
        for neighboring compositions (up to 2)
"""

from pymatgen import MPRester
from pymatgen.analysis.structure_matcher import StructureMatcher
import tqdm
import itertools
from atomate.vasp.workflows.presets.core import wf_structure_optimization,\
    wf_static
from atomate.vasp.powerups import add_tags, add_modify_incar
from fireworks import LaunchPad
import sys

mpr = MPRester()

def get_structures_by_template(template, perturbations=1):
    """
    This takes a template, finds all structures with a composition
    with up to n fewer Mn atoms (replaced by Sb) and up to n fewer Sb atoms
    (replaced by Mn).

    Args:
        template (Structure): template structure
        perturbations (int): number of atoms to make substitutions
            for relative to the template composition
    """
    comp = template.composition.reduced_composition.get_el_amt_dict()
    all_structures = []
    for elt, sub in ['Mn', 'Sb'], ['Sb', 'Mn']:
        # I'm sure there's a more elegant way to do this, sorry to future self
        if elt in comp:
            indices = [template.index(site) for site in template
                       if site.species_string==elt]
            for n in range(perturbations+1):
                combos = itertools.combinations(indices, n)
                for combo in list(combos):
                    new_structure = template.copy()
                    for idx in combo:
                        new_structure[idx] = sub
                        all_structures.append(new_structure)
    unique_structs = []
    structure_matcher = StructureMatcher()
    for new_struct in tqdm.tqdm(all_structures, desc='reducing structures'):
        if not any([structure_matcher.fit(new_struct, struct) for struct in unique_structs]):
            unique_structs.append(new_struct)
    assert len(set([s.composition for s in unique_structs])) < 2 * n + 1
    return unique_structs

def get_all_structs_by_material_id(mp_ids):
    return {mpid[0]: get_structures_by_template(
        mpr.get_structure_by_material_id(mpid[0]), mpid[1])
            for mpid in tqdm.tqdm(mp_ids)}

if __name__=="__main__":
    structures = get_all_structs_by_material_id([("mp-19231", 2),
                                                 ("mp-25043", 2),
                                                 ("mp-565203",1),
                                                 ("mp-18759", 2),
                                                 ("mp-19006", 2),
                                                 ("mp-19395", 2),
                                                 ("mp-2136",  2),
                                                 ("mp-1705",  2),
                                                 ("mp-230", 2)])
    # all_structures = itertools.chain.from_iterable(structures.values())
    lpad = LaunchPad.from_file(sys.argv[1])
    for mpid, structures in structures.items():
        for structure in structures:
            wf = wf_structure_optimization(structure)
            wf.append_wf(wf_static(structure), wf.leaf_fw_ids)
            wf = add_modify_incar(wf)
            wf = add_tags(wf, ["mn_sb_calcs", mpid[0]])
            lpad.add_wf(wf)
    # print("\n".join("{}: {}".format(k, len(v)) for k, v in structures.items()))
    # import pdb; pdb.set_trace()
