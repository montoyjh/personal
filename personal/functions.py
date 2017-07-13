import os, sys, traceback, pdb

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
