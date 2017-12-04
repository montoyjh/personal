"""
A function dump for stuff I find useful, often
to be reorganized into external modules later
if prudent
"""

import os, sys, traceback, pdb
from atomate.vasp.database import VaspCalcDb
from monty.serialization import loadfn
from pymatgen import vis
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view

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
    print "Successful"

def pymatview(struct):
    """
    Helper function to use ase-gui with pymatgen structure
    """
    aaa = AseAtomsAdaptor()
    if isinstance(struct, list):
        view([aaa.get_atoms(s) for s in struct])
    else:
        view(aaa.get_atoms(struct))

if __name__=="__main__":
    pdb_function(get_db)
