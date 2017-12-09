"""
This module is intended to be a temporary placeholder for my work on transitioning the elastic workflow

Main level actions:
    - Add elastic fireworks for all materials
    - Prioritize fireworks by symmetrically distinct deformations
    - Defuse fws with entries already in materials collection
        
Todo:
    - Put my and Wei Chen's tasks into a single database (translating to atomate schema)
    - Build elasticity collection
    - Build elasticity into materials collection
    - Hook into propjockey
"""

from maggma.stores import MongoStore
from atomate.vasp.workflows.presets.core import wf_elastic_constant
from atomate.vasp.powerups import add_tags
from atomate.utils.utils import get_fws_and_tasks
from pymatgen.analysis.elasticity.tensors import Tensor, SquareTensor,\
        get_tkd_value, symmetry_reduce
from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.core.operations import SymmOp
from pymatgen import Structure
from fireworks import LaunchPad
import numpy as np
import tqdm
import logging
from multiprocessing import Pool

logger = logging.getLogger(__name__)

def get_structures(materials_store, lu_filter=None, criteria=None, limit=None):
    """
    Simple function to get all structures 

    Args:
        materials_store (Store): store of materials
        lu_filter (datetime?): date past which to query for materials
        criteria (dict) criteria to filter documents accessed from materials
        limit (int) limit for number of documents to get, primarily for testing
    """
    criteria = criteria or {}
    if lu_filter:
        criteria.update(lu_filter)
    # data = materials_store.collection.aggregate(pipeline).next()
    data = materials_store.query(["structure", "task_id"], criteria)
    count = limit or data.count()
    if limit:
        data = data.limit(limit)
    result = {d["task_id"]: Structure.from_dict(d['structure'])
              for d in tqdm.tqdm(data, total=count)}
    return result

# Just for multiprocessing.pool.imap
def pgenerate_workflow(input_args):
    structure, tags = input_args
    return generate_workflow(structure, tags)

def generate_workflow(structure, tags=[]):
    """
    Generates a standard production workflow.

    Notes:
        Uses a primitive structure transformed into
        the conventional basis (for equivalent deformations).

        Adds the "minimal" category to the minimal portion
        of the workflow necessary to generate the elastic tensor,
        and the "minimal_full_stencil" category to the portion that
        includes all of the strain stencil, but is symmetrically complete
    """
    # transform the structure
    ieee_rot = Tensor.get_ieee_rotation(structure)
    assert SquareTensor(ieee_rot).is_rotation(tol=0.005)
    symm_op = SymmOp.from_rotation_and_translation(ieee_rot)
    ieee_structure = structure.copy()
    ieee_structure.apply_operation(symm_op)

    # construct workflow
    wf = wf_elastic_constant(ieee_structure)
    
    # Set categories, starting with optimization
    opt_fws = get_fws_and_tasks(wf, fw_name_constraint="optimization")
    wf.fws[opt_fws[0][0]].spec['elastic_category'] = "minimal"

    # find minimal set of fireworks using symmetry reduction
    fws_by_strain = {Strain(fw.tasks[-1]['pass_dict']['strain']): n
                     for n, fw in enumerate(wf.fws) if 'deformation' in fw.name}
    unique_tensors = symmetry_reduce(fws_by_strain.keys(), ieee_structure)
    for unique_tensor in unique_tensors:
        fw_index = get_tkd_value(fws_by_strain, unique_tensor)
        if np.isclose(unique_tensor, 0.005).any():
            wf.fws[fw_index].spec['elastic_category'] = "minimal"
        else:
            wf.fws[fw_index].spec['elastic_category'] = "minimal_full_stencil"
            
    # Add tags
    if tags:
        wf = add_tags(wf, tags)

    return wf

def defuse_wflows_with_elasticity_data(materials_store, lpad=None):
    """
    function to defuse workflows with existing data perhaps to
    be reignited later
    
    Args:
        materials_store (Store): store of materials documents
        lpad (LaunchPad): launchpad object to which to add fws,
            if None, defaults to auto-loaded launchpad
    """
    lpad = lpad or LaunchPad.auto_load()
    q = {'elasticity': {"$ne": None}}
    mpids_to_defuse = materials_store.distinct("task_id", criteria=q)
    tags = lpad.fireworks.distinct("spec.tags")
    mpids_to_defuse = list(set(tags).intersection(set(mpids_to_defuse)))
    for mpid in mpids_to_defuse:
        wf_id = lpad.fireworks.find_one({"spec.tags": mpid}, ['fw_id'])['fw_id']
        lpad.defuse_wf(wf_id)
        logger.info("Defusing wf with wf_id {} and mp_id {}".format(wf_id, mpid))

if __name__ == "__main__":
    materials_store = MongoStore.from_db_file("tests/materials.yaml")
    materials_store.connect()
    lpad = LaunchPad.from_file("tests/my_launchpad.yaml")
    lpad.reset(password='', require_password=False, max_reset_wo_password=101)
    # still testing so limit to 100 docs
    structures_by_mpid = get_structures(materials_store, limit=100)
    parallel = True
    if parallel:
        gwf_args = [(structure, [mpid, 'production_elastic']) 
                    for mpid, structure in structures_by_mpid.items()]
        pool = Pool(7)
        wfs = list(tqdm.tqdm(pool.imap_unordered(pgenerate_workflow, gwf_args), 
                                                 total=len(gwf_args)))
        pool.close()
    else:
        wfs = [generate_workflow(structure, tags=[mpid, 'production_elastic'])
               for mpid, structure in tqdm.tqdm(structures_by_mpid.items())]

    import nose; nose.tools.set_trace()
    #for wf in wfs:
    #    lpad.add_wf(wf)
    defuse_wflows_with_elasticity_data(materials_store, lpad)
