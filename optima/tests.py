"""
Test script to see if Optima works.

To use: comment out lines in the definition of 'tests' to not run those tests.

NOTE: for best results, run in interactive mode, e.g.

python -i tests.py

Version: 2015nov01 by cliffk
"""



## Define tests to run here!!!
tests = [
'makespreadsheet',
'makeproject',
'saveload',
'loadspreadsheet',
'runsim',
'makeprograms']
#'gui'
#]

numericalassertions = True # Whether or not to actually run things and test their values
doplot = True # Whether or not to show diagnostic plots

runalltests=True

##############################################################################
## Initialization
##############################################################################

from utils import tic, toc, blank, pd # analysis:ignore

def done(t=0):
    print('Done.')
    toc(t)
    blank()
    





blank()
print('Running tests:')
for i,test in enumerate(tests): print(('%i.  '+test) % (i+1))
blank()



##############################################################################
## The tests
##############################################################################

T = tic()


## Spreadsheet creation test
if 'makespreadsheet' in tests:
    t = tic()
    print('Running make spreadsheet test...')
    from makespreadsheet import makespreadsheet
    makespreadsheet()
    done(t)



## Project creation test
if 'makeproject' in tests:
    t = tic()
    print('Running make project test...')
    from project import Project
    P = Project()
    print(P)
    done(t)




## Project save/load test
if 'saveload' in tests:
    t = tic()
    print('Running save/load test...')
    
    from project import Project, load
    filename = 'testproject.prj'
    
    print('  Checking saving...')
    P = Project()
    P.save(filename)
    
    print('  Checking loading...')
    Q = load(filename)
    Q.save()
    Q.loadfromfile()
    
    print('  Checking defaults...')
    Z = Project()
    Z.save()
    
    done(t)




## Load spreadsheet test
if 'loadspreadsheet' in tests:
    t = tic()
    print('Running loadspreadsheet test...')
    from project import Project
    
    print('  Create a project from a spreadsheet')
    P = Project(spreadsheet='test.xlsx')
    
    print('  Load a project, then load a spreadsheet')
    Q = Project()
    Q.loadspreadsheet('test.xlsx')
    
    if numericalassertions:
        assert Q.data['const']['effcondom'][0]==0.05, 'Condom efficacy not 95% or not being read in properly'
    
    done(t)




## Run simulation test
if 'runsim' or 'gui' in tests:
    t = tic()
    print('Running runsim test...')
    
    from project import Project
    P = Project(spreadsheet='test.xlsx')
    results = P.runsim('default')
    
    done(t)


## Project creation test
if 'makeprograms' in tests:
    t = tic()

    print('Running make programs test...')
    from programs import Program, Programset

    # First set up some programs. Programs need to be initialized with a name. Often they will also be initialized with targetpars and targetpops
    HTC = Program(name='HTC', targetpars=[{'param': 'hivtest', 'pop': 'MSM'},{'param': 'hivtest', 'pop': 'FSW'}],targetpops=['FSW','MSM'])

    # Run additional tests if asked
    if runalltests:
        FSW = Program(name='FSW programs', targetpars=[{'param': 'hivtest', 'pop': 'FSW'},{'param': 'condoms', 'pop': 'FSW'}], targetpops=['FSW'])
        MGT = Program('MGT')
        ART = Program(name='ART', targetpars=[{'param': 'numtx', 'pop': 'Total'}],targetpops=['Total'])
    
        # Testing methods of program class
        # 1. Adding a target parameter to a program
        HTC.addtargetpar({'param': 'hivtest', 'pop': 'Males 15-49'})
            
        # 2. Removing a target parameter from a program
        HTC.rmtargetpar({'param': 'hivtest', 'pop': 'Males 15-49'})
    
        # 3. Add historical cost-coverage data point
        HTC.addcostcovdatum({'t':2013,'cost':1e6,'coverage':3e5})
    
        # 4. Overwrite historical cost-coverage data point
        HTC.addcostcovdatum({'t':2013,'cost':2e6,'coverage':3e5}, overwrite=True)
    
        # 5. Remove historical cost-coverage data point - specify year only
        HTC.rmcostcovdatum(2013)
    
        # 6. Add parameters for defining cost-coverage function.
        HTC.costcovfn.addccopar({'saturation': 0.8, 't': 2013.0, 'unitcost': 30})
        HTC.costcovfn.addccopar({'t': 2016.0, 'unitcost': 30})
        HTC.costcovfn.addccopar({'t': 2017.0, 'unitcost': 30})
    
        # 7. Overwrite parameters for defining cost-coverage function.
        HTC.costcovfn.addccopar({'t': 2016.0, 'unitcost': 25},overwrite=True)
    
        # 8. Remove parameters for defining cost-coverage function.
        HTC.costcovfn.rmccopar(2017)
    
        # 9. Get parameters for defining cost-coverage function for any given year (even if not explicitly entered).
        HTC.costcovfn.getccopar(2014)
    
        # 10. Evaluate cost-coverage function to get coverage for a given year, spending amount and population size
        HTC.costcovfn.evaluate(1e6,1e5,2015)
        HTC.getcoverage(1e6,1e5,2015) # Two equivalent ways to do this, probably redundant

    print('Running make programs set test...')
    R = Programset(programs={'HTC':HTC,'FSW':FSW,'MGT':MGT})

    # Run additional tests if asked
    if runalltests:
        # Testing methods of program class
        # 1. Adding a program
        R.addprog({'ART':ART})
    
        # 2. Removing a program
        R.rmprog('ART')
        
        # 3. See which programs are optimizable
        R.optimizable()
    
        # 4. Produce a dictionary whose keys are populations targeted by some 
        #    program, and values are the programs that target them
        R.progs_by_targetpop()
    
        # 5. Produce a dictionary whose keys are paramter types targeted by some 
        #    program, and values are the programs that target them
        R.progs_by_targetpartype()
    
        # 6. Produce a dictionary whose keys are paramter types targeted by some 
        #    program, and values are dictionaries  whose keys are populations 
        #    targeted by some program, and values are the programs that target them
        R.progs_by_targetpar()
    
        # 7. Get a vector of coverage levels corresponding to a vector of program allocations
        R.getcoverage(tvec=None,budget=None)
    
        # 8. Get a set of parameter values corresponding to a vector of program allocations
        R.getoutcomes(tvec=None,budget=None)

    done(t)


## Run the GUI
if 'gui' in tests:
    t = tic()
    print('Running GUI test...')
    
    from gui import gui
    gui(results)
    
    done(t)




print('\n\n\nDONE: ran %i tests' % len(tests))
toc(T)