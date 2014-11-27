"""
TEST_OPTIMA

While optima.py is a demonstration of everything Optima can do, this is used to
test specific features.

Version: 2014nov27 by cliffk
"""


print('WELCOME TO OPTIMA')

## Set parameters
projectname = 'example'
verbose = 4
show_wait = False

print('\n\n\n1. Making project...')
from makeproject import makeproject
D = makeproject(projectname=projectname, pops=['']*6, progs = ['']*7, datastart=2000, dataend=2015, verbose=verbose)

print('\n\n\n2. Updating data...')
from updatedata import updatedata
D = updatedata(D, verbose=verbose)

print('\n\n\n3. Running simulation...')
from runsimulation import runsimulation
D = runsimulation(D, startyear=2000, endyear=2015, verbose=verbose)

print('\n\n\n4. Viewing results...')
from viewresults import viewepiresults
viewepiresults(D.plot.E, whichgraphs={'prev':[1,1], 'inci':[0,1], 'daly':[0,1], 'death':[0,1], 'dx':[0,1], 'tx1':[0,1], 'tx2':[0,1]}, startyear=2000, endyear=2015, onefig=True, verbose=verbose, show_wait=show_wait)

print('\n\n\nDONE.')