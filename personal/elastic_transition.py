"""
This module is intended to be a temporary placeholder for my work on transitioning the elastic workflow
"""

from maggma.stores import MongoStore
from atomate.vasp.workflows.presets.core import wf_elastic_constant
from atomate.utils.utils import get_fws_and_tasks
from pymatgen.analysis.elasticity.tensors import Tensor, SquareTensor,\
        get_tkd_value, symmetry_reduce
from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.core.operations import SymmOp
import numpy as np
import tqdm
import logging

logger = logging.get_logger(__module__)

def get_structures(materials_store, lu_filter=None, criteria=None, limit=None):
    """
    Simple function to get all structures 

    Args:
        materials_store (Store): store of materials
        lu_filter (datetime?): date past which to query for materials
    """
    criteria = criteria or {}
    if lu_filter:
        criteria.update(lu_filter)
    # data = materials_store.collection.aggregate(pipeline).next()
    data = materials_store.query(["structure", "task_id"], criteria)
    count = limit or data.count()
    result = {d["task_id"]: d['structure'] for d in tqdm.tqdm(data, total=count)}
    return result

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
    assert SquareTensor(ieee_rot).is_rotation()
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
    """
    lpad = lpad or LaunchPad.auto_load()
    q = {'elasticity': {"$ne": None}}
    mpids_to_defuse = materials_store.distinct("material_id", criteria=q)
    for mpid in mpids_to_defuse:
        wf_id = lpad.fireworks.find_one({"spec.tags": mpid}, ['fw_id'])['fw_id']
        lpad.defuse_wflows(wf_id)
        logger.info("Defusing wf with wf_id {} and mp_id {}".format(wf_id, mpid))
