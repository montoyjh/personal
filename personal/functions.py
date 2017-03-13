import os, sys, traceback, pdb
from atomate.vasp.database import MMVaspDb
from monty.serialization import loadfn
from pymatgen import vis

def get_db(fworker_filename=None):
    """
    This function gets a db_file based on the db_file key in a given fireworker
    """
    if fworker_filename is None:
        fw_config = os.environ.get("FW_CONFIG_FILE", None)
        if fw_config is None:
            raise ValueError("No FW_CONFIG_FILE is set.")
        fworker = loadfn(os.path.join(os.path.dirname(fw_config), "my_fworker.yaml"))
    else:
        fworker = loadfn(fworker_filename)
    db_file = MMVaspDb.from_db_file(fworker["env"]["db_file"])
    return db_file.db


def get_colors(key="Jmol"):
    """
    Gets dictionary of matplotlib colors correponding to elements
    from pymatgen color library
    """
    colors = loadfn(os.path.join(os.path.dirname(vis.__file__), "ElementColorSchemes.yaml"))
    color_dict = {el:[j / 256. for j in colors[key][el]] for el in colors[key].keys()}
    return color_dict

if __name__=="__main__":
    try:
        db = get_db()
        colors = get_colors()
    except:
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
