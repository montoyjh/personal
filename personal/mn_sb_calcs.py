"""
Script that generates the necessary calculations for collaboration
with John Gregoire.  Specifically, multiphase orderings of the
Mn-Sb chemical system.

Notes:
    - originally tried to make this work for every intermediate composition
        which got crazy fast, so limited it to substitutions of single atoms
        for neighboring compositions (up to 2)
    - an early iteration of this workflow had a bug that didn't properly
        generate static workflows to copy vasp outputs
"""

from pymatgen import MPRester, Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.vasp.sets import MPStaticSet
import tqdm
import itertools
from atomate.utils.utils import get_wf_from_spec_dict
from atomate.vasp.powerups import add_tags, add_modify_incar
from atomate.vasp.workflows.base.core import get_wf as get_atomate_wf

from fireworks import LaunchPad
from pymatgen.command_line.bader_caller import *
from maggma.stores import MongoStore
from monty.serialization import loadfn

import sys

mpr = MPRester()
module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

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
    all_structures = [template]
    subs_template = template.copy()
    if subs_template.num_sites < 7:
        subs_template.make_supercell([2, 2, 2])
    for elt, sub in ['Mn', 'Sb'], ['Sb', 'Mn']:
        # I'm sure there's a more elegant way to do this, sorry to future self
        if elt in comp:
            indices = [subs_template.index(site) for site in subs_template
                       if site.species_string==elt]
            for n in range(perturbations+1):
                combos = itertools.combinations(indices, n)
                for combo in list(combos):
                    new_structure = subs_template.copy()
                    for idx in combo:
                        new_structure[idx] = sub
                        all_structures.append(new_structure)
    unique_structs = []
    structure_matcher = StructureMatcher()
    for new_struct in tqdm.tqdm(all_structures, desc='reducing structures'):
        if not any([structure_matcher.fit(new_struct, struct) for struct in unique_structs]):
            unique_structs.append(new_struct)
    assert len(set([s.composition for s in unique_structs])) <= 2 * n + 1

    # ensure total mn+sb/o stoichiometry preserved
    ratio = (sum(comp.values()) - comp['O']) / comp['O']
    for struct in unique_structs:
        unique_comp = struct.composition.get_el_amt_dict()
        unique_ratio = (sum(unique_comp.values()) - unique_comp['O']) / unique_comp['O']
        assert ratio == unique_ratio

    return unique_structs

def get_all_structs_by_material_id(mp_ids):
    return {mpid[0]: get_structures_by_template(
        mpr.get_structure_by_material_id(mpid[0]), mpid[1])
            for mpid in tqdm.tqdm(mp_ids)}

def get_opt_static_wf(structure, tags=None):
    wf_spec = loadfn(os.path.join(module_dir, 'opt_static.yaml'))
    wf = get_wf_from_spec_dict(structure, wf_spec)
    wf = add_modify_incar(wf)
    if tags:
        wf = add_tags(wf, tags)
    return wf

def add_bader():
    test_store = MongoStore.from_db_file(os.path.join(module_path, "tasks.yaml"))
    test_store.connect()
    docs = test_store.query(['dir_name', 'task_id'], {"tags": "mn_sb_calcs_4",
                                                      "task_label": "static", 
                                                      "bader": {"$exists": False}})
    for doc in tqdm.tqdm(docs, total=docs.count()):
        run_dir = doc['dir_name'].split(':')[-1]
        bader = BaderAnalysis.from_path(run_dir)
        summary = bader.summary
        oxi_structure = bader.get_oxidation_state_decorated_structure()
        summary.update({"oxi_structure": oxi_structure.as_dict()})
        test_store.collection.update({"task_id": doc['task_id']}, {"$set": {"bader": summary}})
        # import pdb; pdb.set_trace()

def get_high_fft_grid_wfs():
    test_store = MongoStore.from_db_file("tasks.yaml")
    test_store.connect()
    docs = test_store.query(['dir_name', 'output.structure', 'task_id', 'tags'], 
                            {"tags": "mn_sb_calcs_4", "task_label": "static"})
    wfs = []
    for doc in docs:
        structure = Structure.from_dict(doc['output']['structure'])
        wf = get_atomate_wf(structure, "static_only.yaml", vis=MPStaticSet(structure,
                            reciprocal_density=500),
                            common_params={"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"})
        # wf = wf_static(structure)
        wf = add_modify_incar(wf)
        wf = add_modify_incar(wf, {"incar_update": {"ENAUG": 5000, "PREC": "High"}})
        wf = add_tags(wf, doc['tags'] + ['dense_grid'])
        wfs.append(wf)
    return wfs

mode = 'add_3_trirutile'
launch = True

if __name__=="__main__":
    if mode == 'submit':
        structures = get_all_structs_by_material_id([("mp-19231", 2),
                                                     ("mp-25043", 2),
                                                     ("mp-565203",1),
                                                     ("mp-18759", 2),
                                                     ("mp-19006", 2),
                                                     ("mp-19395", 2),
                                                     ("mp-2136",  2),
                                                     ("mp-1705",  2),
                                                     ("mp-230", 2)])

        # Also get TiO2 as template
        rutile = mpr.get_structure_by_material_id("mp-2657")
        ti_indices = [rutile.index(s) for s in rutile if s.species_string=='Ti']
        rutile[ti_indices[0]] = 'Mn'
        rutile[ti_indices[1]] = 'Sb'
        structures["mp-2657"] = get_structures_by_template(rutile, 2)

        # Also get tri rutile template
        tri_rutile = mpr.get_structure_by_material_id("mp-24845")
        tri_rutile.replace_species({"Co": "Mn"})
        structures['mp-24845'] = get_structures_by_template(tri_rutile, 2)

        lpad = LaunchPad.from_file(sys.argv[1])
        wfs = []
        for mpid, structures in structures.items():
            for structure in structures:
                wfs.append(get_opt_static_wf(structure, tags=["mn_sb_calcs_4", mpid]))

        if launch:
            for wf in wfs:
                lpad.add_wf(wf)

    # Lazy way of adding larger perturbation of Tri-rutile template
    elif mode == 'add_3_trirutile':
        tri_rutile = mpr.get_structure_by_material_id("mp-24845")
        tri_rutile.replace_species({"Co": "Mn"})
        new_structs = get_structures_by_template(tri_rutile, 3)
        import nose; nose.tools.set_trace()
        new_structs = [s for s in new_structs if s.reduced_composition == Composition()]

        lpad = LaunchPad.from_file(sys.argv[1])
        wfs = []
        for mpid, structures in structures.items():
            for structure in structures:
                wfs.append(get_opt_static_wf(structure, tags=["mn_sb_calcs_4", mpid]))

        if launch:
            for wf in wfs:
                lpad.add_wf(wf)

    elif mode == 'do_bader':
        add_bader()
    elif mode == 'high_fft':
        wfs = get_high_fft_grid_wfs()
        lpad = LaunchPad.from_file(sys.argv[1])
        if launch:
            for wf in wfs:
                lpad.add_wf(wf)
    else:
        raise ValueError("Mode {} not supported".format(mode))
