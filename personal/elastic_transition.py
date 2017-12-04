"""
This module is intended to be a temporary placeholder for my work on transitioning the elastic workflow
"""

from maggma.stores import MongoStore

def get_structures(materials_collection_filename, lu_date):
    """
    Simple function to get all structures 
    """
    materialsstore = MongoStore.from_db_file(materials_collection_filename)
    pipeline = [{"$match": {"last_updated": {"$gt": lu_date}}}] if lu_date else []
    pipeline += [{"$group": {"_id": "$material_id", "structure": "$output.structure"}}]
    data = materialsstore.collection.aggregate(pipeline)
    return data
