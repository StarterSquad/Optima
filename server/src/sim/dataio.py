"""
DATAIO

Data input/output. Uses JSON format.

Version: 2014nov17 by cliffk
"""

from printv import printv
import os


DATADIR="/tmp/uploads"
TEMPLATEDIR = "/tmp/templates"
PROJECTDIR = "/tmp/projects"



def fullpath(filename, datadir=DATADIR):
    """
    "Normalizes" filename:  if it is full path, leaves it alone. Otherwise, prepends it with datadir.
    """
    
    result = filename

    # get user dir path
    datadir = upload_dir_user(datadir)

    if not(os.path.exists(datadir)):
        os.makedirs(datadir)
    if os.path.dirname(filename)=='' and not os.path.exists(filename):
        result = os.path.join(datadir, filename)

    return result

def templatepath(filename):
    return fullpath(filename, TEMPLATEDIR)

def projectpath(filename):
    return fullpath(filename, PROJECTDIR)




def savedata(filename, data, update=True, verbose=2, path=None):
    """
    Saves the pickled data into the file (either updates it or just overwrites).
    """
    printv('Saving data...', 1, verbose)
    from cPickle import dump, load
    
    filename = projectpath(filename)

    try: # First try loading the file and updating it
        rfid = open(filename,'rb') # "Read file ID" -- This will fail if the file doesn't exist
        origdata = load(rfid)
        if update: origdata.update(data)
        else: origdata = data
        wfid = open(filename,'wb')
        dump(data, wfid)
        printv('..updated file', 3, verbose)
    except: # If that fails, save a new file
        wfid = open(filename,'wb')
        dump(data, wfid)
        printv('..created new file', 3, verbose)
    printv(' ...done saving data at %s.' % filename, 2, verbose)
    return filename




def loaddata(filename, verbose=2):
    """
    Loads the file and unpickles data from it.
    """
    from cPickle import load
    printv('Loading data...', 1, verbose)
    if not os.path.exists(filename):
        filename = projectpath(filename)
    rfid = open(filename,'rb')
    data = load(rfid)

    printv('...done loading data.', 2, verbose)
    return data



def upload_dir_user(dirpath):
    
    try:
        from flask.ext.login import current_user

        # get current user 
        if current_user.is_anonymous() == False:
    
            # user_path
            user_path = os.path.join(dirpath, str(current_user.id))
    
            # if dir does not exist
            if not(os.path.exists(dirpath)):
                os.makedirs(dirpath)
    
            # if dir with user id does not exist
            if not(os.path.exists(user_path)):
                os.makedirs(user_path)
            
            return user_path
    except:
        return dirpath
    
    return dirpath
