"""
TEST_OPTIMIZATION

This function tests that the optimization is working.

Version: 2015jan19 by cliffk
"""

dotimevarying = True # True False


print('WELCOME TO OPTIMA')

## Set parameters
projectname = 'example'
verbose = 2
ntimepm = 2 # AS: Just use 1 or 2 parameters... using 3 or 4 can cause problems that I'm yet to investigate
timelimit = 500

print('\n\n\n1. Making project...')
from makeproject import makeproject
D = makeproject(projectname=projectname, pops=['']*6, progs = ['']*7, datastart=2000, dataend=2015, verbose=verbose)

print('\n\n\n2. Updating data...')
from updatedata import updatedata
D = updatedata(D, verbose=verbose, savetofile=False)

print('\n\n\n3. Running simulation...')
from runsimulation import runsimulation
D = runsimulation(D, startyear=2000, endyear=2030, verbose=verbose, dosave=False)

print('\n\n\n4. Running optimization...')
from optimize import optimize
optimize(D, objectives={"year":{"start":2000,"end":2030}}, ntimepm=2, timelimit=timelimit, verbose=verbose)

    
print('\n\n\n5. Viewing optimization...')
from viewresults import viewmultiresults, viewallocpies
viewmultiresults(D.plot.OM)
viewallocpies(D.plot.OA)

print('\n\n\nDONE.')