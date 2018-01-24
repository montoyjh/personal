import os, sys, traceback, pdb
from atomate.vasp.database import VaspCalcDb
from atomate.utils.utils import get_wf_from_spec_dict
from atomate.vasp.powerups import add_modify_incar, add_tags
from monty.serialization import loadfn
from pymatgen import vis

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

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
    db_file = VaspCalcDb.from_db_file(fworker["env"]["db_file"])
    return db_file.db


def get_colors(key="Jmol"):
    """
    Gets dictionary of matplotlib colors correponding to elements
    from pymatgen color library
    """
    colors = loadfn(os.path.join(os.path.dirname(vis.__file__), "ElementColorSchemes.yaml"))
    color_dict = {el:[j / 256. for j in colors[key][el]] for el in colors[key].keys()}
    return color_dict

def pdb_function(function, *args, **kwargs):
    """
    This is a function that will attempt to execute another function
    and invoke the debugger if that function's execution produces
    an error
    """
    try:
        function(*args, **kwargs)
    except:
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
    print("Successful")

def get_opt_static_wf(structure, tags=None):
    wf_spec = loadfn(os.path.join(module_dir, 'opt_static.yaml'))
    wf = get_wf_from_spec_dict(structure, wf_spec)
    wf = add_modify_incar(wf)
    if tags:
        wf = add_tags(wf, tags)
    return wf
