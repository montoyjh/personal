"""
Script used to generate the workflows for collaboration with Seoin Back,
Samira Siaohrostami of SUNCAT.  Gets the structures from Ivano's old
database and submits atomate workflows to the launchpad
"""

from ase import db
import os
import pandas as pd
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.analysis.structure_analyzer import oxide_type
import tqdm
from monty.json import jsanitize
from fireworks import LaunchPad
from personal.functions import get_opt_static_wf
from maggma.stores import MongoStore
from atomate.utils.utils import get_mongolike

aaa = AseAtomsAdaptor()

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

def generate_structures():
    """
    Short function to generate all the necessary structures from the
    csv file provided and the cubic_perovskites.db file
    """
    con = db.connect(os.path.join(module_dir, "cubic_perovskites.db"))
    dataframe = pd.read_csv(os.path.join(module_dir, "screened.csv"),
                            index_col="Formula")
    data = con.select(combination='ABO3')
    all_atoms_objects = {"{}{}O3".format(row.A_ion, row.B_ion): row.toatoms()
                         for row in data}
    structures = []
    #for formula in dataframe.index:
    #    structures.append(aaa.get_structure(all_atoms_objects[formula]))
    # just do everything
    for formula, atoms in all_atoms_objects.items():
        if not formula in dataframe.index:
            structures.append(aaa.get_structure(atoms))
    return structures

def get_simplified_docs(store, criteria):
    """
    
    Args:
        store (Store): store of task documents
        criteria (dict): filter criteria

    Returns:
        list of simplified documents
        
    returns a list of documents to be inserted into external database,
    primarily for sharing with SUNCAT collaborators
    """
    store.connect()
    results = store.query(["dir_name", "output", "task_id", "input"],
                          criteria=criteria)
    docs = []
    for result in tqdm.tqdm(results, total=results.count()):
        result.pop('_id')
        output = result.pop('output')
        structure = Structure.from_dict(output['structure'])
        result['formula'] = '{}{}O3'.format(*[s.species_string for s in structure[:2]])
        params = ['input.pseudo_potential', 'input.is_hubbard',
                  'input.hubbards']
        params = {k.split('.')[-1]: get_mongolike(result, k) for k in params}
        params['potcar_symbols'] = [
            "%s %s" % (params["pseudo_potential"]["functional"], l)
            for l in params["pseudo_potential"]["labels"]]
        params['oxide_type'] = oxide_type(structure)
        result['entry'] = ComputedStructureEntry(structure, output['energy'],
                                                 parameters=params)
        result.pop('input')
        result = jsanitize(result, strict=True)
        docs.append(result)
    return docs

launchpad_filename = "/Users/josephmontoya/fw_config/config_test/my_launchpad.yaml"
test_db_file = os.path.join(module_dir, "tasks.yaml")

mode = "make_db"

if __name__ == "__main__":
    if mode == "create_workflows":
        structures = generate_structures()
        # This uses the input from the script running to add info
        lpad = LaunchPad.from_file(launchpad_filename)
        for structure in structures:
            workflow = get_opt_static_wf(structure, tags=['peroxide_catalysts_2'])
            lpad.add_wf(workflow)
    elif mode == "make_db":
        test_store = MongoStore.from_db_file(test_db_file)
        docs = get_simplified_docs(test_store,
                                   criteria={"tags": "peroxide_catalysts_2"})

        target = MongoStore.from_db_file(os.path.join(module_dir, "mangomiracle.yaml"))
        target.collection.drop()
        target.connect()
        target.update(docs, key='task_id')
    else:
        raise ValueError("Mode {} not supported".format(mode))