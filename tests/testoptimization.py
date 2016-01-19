"""
Test optimization

To use: comment out lines in the definition of 'tests' to not run those tests.
NOTE: for best results, run in interactive mode, e.g.
python -i tests.py

Version: 2016jan05 by cliffk
"""

## Define tests to run here!!!
tests = [
'setup',
'minimizeoutcomes',
]


##############################################################################
## Initialization -- same for every test script
##############################################################################

from optima import tic, toc, blank, pd # analysis:ignore

if 'doplot' not in locals(): doplot = True

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






## Set up project etc.
if 'setup' in tests:
    t = tic()

    print('Running standard scenarios test...')
    from optima import Project, Program, Programset
    
    P = Project(spreadsheet='test7pops.xlsx')
    pops = P.data['pops']['short']
    caspships = P.data['pships']['cas']
    compships = P.data['pships']['com']
    
    condprog = Program(name='Condoms',
                  targetpars=[{'param': 'condcas', 'pop': caspship} for caspship in caspships],
                  targetpops=pops,
                  category='Prevention',
                  short='Condoms',
                  criteria = {'hivstatus': 'allstates', 'pregnant': False})    
                  
    fswprog = Program(name='FSW_programs',
                  targetpars=[{'param': 'condcom', 'pop': compship} for compship in [x for x in compships if 'FSW' in x]] + [{'param': 'condcas', 'pop': caspship} for caspship in [x for x in caspships if 'FSW' in x]] + [{'param': 'hivtest', 'pop': 'FSW'}],
                  targetpops=['FSW'],
                  category='Prevention',
                  short='FSW programs',
                  criteria = {'hivstatus': 'allstates', 'pregnant': False})
    
    condprog.costcovfn.addccopar({'saturation': (0.75,0.75), 't': 2016.0, 'unitcost': (30,40)})
    fswprog.costcovfn.addccopar({'saturation': (0.9,0.9), 't': 2016.0, 'unitcost': (50,80)})
    
    R = Programset(programs=[condprog, fswprog]) 
    
    R.covout['condcas'][('Clients', 'FSW')].addccopar({'intercept': (0.3,0.35), 't': 2016.0, 'Condoms':(0.45,0.55), 'FSW_programs':(0.55,0.65)})
    R.covout['condcas'][('Clients', 'F 15+')].addccopar({'intercept': (0.2,0.3), 't': 2016.0, 'Condoms':(0.35,0.45)})
    R.covout['condcas'][('MSM', 'MSM')].addccopar({'intercept': (0.5,0.55), 't': 2016.0, 'Condoms':(0.55,0.65)})
    R.covout['condcas'][('M 15+', 'FSW')].addccopar({'intercept': (0.3,0.35), 't': 2016.0, 'Condoms':(0.45,0.55), 'FSW_programs':(0.55,0.65)})
    R.covout['condcas'][('M 15+', 'F 15+')].addccopar({'intercept': (0.2,0.3), 't': 2016.0, 'Condoms':(0.35,0.45)})
    R.covout['condcom'][('Clients', 'FSW')].addccopar({'intercept': (0.6,0.65), 't': 2016.0, 'FSW_programs':(0.9,0.95)})
    R.covout['hivtest']['FSW'].addccopar({'intercept': (0.35,0.45), 't': 2016.0, 'FSW_programs':(0.6,0.65)})
    P.addprogset(name='default', progset=R)
    
    done(t)






## Minimize money test
if 'minimizeoutcomes' in tests:
    t = tic()

    print('Running minimize outcomes test...')
    from optima import defaultobjectives
    objectives = defaultobjectives()
    objectives['budget'] = 5e6 # Change default budget to optimize
    P.minoutcomes(name='optim', parset='default', progset='default', objectives=objectives, method='asd', maxiters=5)
    
    print('Original allocation: '),
    print(P.results[-1].budget['orig'])
    print('Optimal allocation: '),
    print(P.optims[-1].getresults().budget['optim']) # Showing that results are "stored" in the optimization -- same object as before
    if doplot: 
        from optima import plotmismatch, plotallocs, plotepi
        plotmismatch(P.results[-1])
        plotallocs(P.results[-1])
        plotepi(P.results[-1])
    
    done(t)



print('\n\n\nDONE: ran %i tests' % len(tests))
toc(T)