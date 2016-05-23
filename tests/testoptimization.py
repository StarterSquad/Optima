"""
Test optimization

To use: comment out lines in the definition of 'tests' to not run those tests.
NOTE: for best results, run in interactive mode, e.g.
python -i tests.py

Version: 2016feb07
"""

## Define tests to run here!!!
tests = [
'minimizeoutcomes',
#'minimizemoney',
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

from optima import defaults

T = tic()








## Minimize outcomes test
if 'minimizeoutcomes' in tests:
    t = tic()

    print('Running minimize outcomes test...')
    from optima import defaultobjectives, defaultconstraints
    
    P = defaults.defaultproject('best') 
    
    objectives = defaultobjectives(P.progsets[0]) # This or P
    objectives['budget'] *= 0.5
    constraints = defaultconstraints(P) # This or P.progsets[0]
    P.optimize(name='minoutcome', parsetname='default', progsetname='default', objectives=objectives, method='asd', maxtime=10, ccsample='best')
    
    print('Original allocation: '),
    print(P.results[-1].budget[0])
    print('Optimal allocation: '),
    print(P.optims[-1].getresults().budget[1]) # Showing that results are "stored" in the optimization -- same object as before
    if doplot: 
        from optima import pygui
        pygui(P.results[-1], toplot=['budget', 'improvement', 'prev-tot', 'prev-per', 'numinci-tot'])
    
    done(t)






## Minimize money test
if 'minimizemoney' in tests:
    t = tic()

    print('Running minimize money test...')
    from optima import defaultobjectives, defaultconstraints
    
    P = defaults.defaultproject('best')
    
    objectives = defaultobjectives(which='money', progset=P.progsets[0])
    objectives['deathfrac'] = 0.1 # Yes, this means an increase in deaths
    objectives['incifrac'] = 0.2
    constraints = defaultconstraints(P.progsets[0])
    P.optimize(name='minmoney', parsetname='default', progsetname='default', objectives=objectives, constraints=constraints, maxtime=10)
    
    print('Original allocation: ($%g)' % sum(P.results[-1].budget[0][:]))
    print(P.results[-1].budget[0])
    print('Optimal allocation: ($%g)' % sum(P.optims[-1].getresults().budget[1][:]))
    print(P.optims[-1].getresults().budget[1]) # Showing that results are "stored" in the optimization -- same object as before
    if doplot: 
        from optima import pygui
        pygui(P.results[-1], toplot=['budget', 'improvement', 'prev-tot', 'prev-per', 'numinci-tot'])
    
    done(t)





print('\n\n\nDONE: ran %i tests' % len(tests))
toc(T)