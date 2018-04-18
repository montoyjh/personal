"""
This module is intended to be a temporary placeholder for my work on
transitioning the elastic workflow

Main level actions:
    - Add elastic fireworks for all materials
    - Prioritize fireworks by symmetrically distinct deformations
    - Defuse fws with entries already in materials collection

Todo:
    - Build elasticity into materials collection
    - Hook into propjockey
    - update web code for getting elastic properties
"""

from maggma.stores import MongoStore
from maggma.runner import Runner
from emmet.vasp.elastic import ElasticBuilder
from emmet.vasp.mpworks import MPWorksCompatibilityBuilder
from emmet.materials.property_workflows import PropertyWorkflowBuilder,\
    get_elastic_wf_builder
from fireworks import LaunchPad
import logging
from monty.serialization import dumpfn
import os

logger = logging.getLogger(__name__)




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

def set_priority(lpad):
    lpad.fireworks.update_many({"spec.elastic_category": "minimal"}, {"$inc": {"spec._priority": 2000}})
    lpad.fireworks.update_many({"spec.elastic_category": "minimal_full_stencil"},
                               {"$inc": {"spec._priority": 1000}})
    lpad.fireworks.update_many({"spec.elastic_category": "minimal", "name": {"$regex": "deformation"}},
                               {"$inc": {"spec._priority": 1}})

def chunk_iterable(iterable, chunk_size):
    """
    Yield successive n-sized chunks from l. Stolen from stackoverflow.

    Args:
        iterable (iterable)
        chunk_size (int)
    """
    for i in range(0, len(iterable), chunk_size):
        yield list(iterable)[i:i + chunk_size]

# Global params, should set these up with an argparser maybe
material_filter = None
"""
material_filter = {"nsites": {"$lt": 5}, 
                   "elasticity.warnings": None, 
                   "elasticity": {"$ne": None}}
                   """
limit = None
reset = True
parallel = True
defuse_existing = False
chunk_size = 300

# Files
mat_file = os.path.join(os.path.dirname(__file__), 'tests', 'materials.yaml')
fw_file = os.path.join(os.path.dirname(__file__), 'tests', 'my_launchpad.yaml')

# Testing
test = False
test_mat = os.path.join(os.path.dirname(__file__), 'tests', 'test_materials.yaml')

# Dump runner to file
if __name__ == "__main__":
    materials = MongoStore.from_db_file("materials.yaml")
    elasticity = MongoStore.from_db_file("elasticity.yaml")
    jhm_mpworks_tasks = MongoStore.from_db_file("jhm_mpworks_tasks.yaml")
    wc_mpworks_tasks = MongoStore.from_db_file("wc_mpworks_tasks.yaml")

    atomate_tasks = MongoStore.from_db_file("atomate_tasks.yaml")
    lpad = LaunchPad.from_file(
        "/Users/josephmontoya/fw_config/config_kpoints/my_launchpad.yaml")
    ewf_builder = get_elastic_wf_builder(elasticity, materials, lpad)

    wc_mpworks_builder = MPWorksCompatibilityBuilder(wc_mpworks_tasks, atomate_tasks)
    jhm_mpworks_builder = MPWorksCompatibilityBuilder(jhm_mpworks_tasks, atomate_tasks,
                                                      incremental=True)
    elastic_builder = ElasticBuilder(atomate_tasks, elasticity, materials)
    runner_wc = Runner([wc_mpworks_builder])
    runner = Runner([jhm_mpworks_builder, elastic_builder, ewf_builder])
    dumpfn(runner_wc, "wc_runner.json", indent=4)
    dumpfn(runner, "elastic_wf_builder.json", indent=4)
    dumpfn(Runner([elastic_builder]), "elastic_builder.json", indent=4)
