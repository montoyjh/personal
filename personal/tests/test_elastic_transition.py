import unittest

from maggma.stores import MongoStore
from personal.elastic_transition import get_structures, generate_workflow
from pymatgen.util.testing import PymatgenTest

# Note: you have to supply your own materials collection file
db_filename = 'materials.yaml'

class ElasticTransitionTest(PymatgenTest):
    def setUp(self):
        self.materials_store = MongoStore.from_db_file(db_filename)
        self.materials_store.connect()

    def test_get_structures(self):
        structures = get_structures(self.materials_store, limit=100)
        blargh

    def test_generate_workflow(self):
        si = self.get_structure("Si")
        structure = generate_workflow(si)
        blargh
