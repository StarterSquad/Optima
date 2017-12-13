"""
Test calibration

To use: comment out lines in the definition of 'tests' to not run those tests.
NOTE: for best results, run in interactive mode, e.g.
python -i tests.py

Version: 2016feb06 by cliffk
"""

## Define tests to run here!!!
tests = [
#'attributes',
#'sensitivity',
'manualfit',
#'autofit',
#'autofitmulti',
# 'longfit',
# 'debugautofit',
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




## Attributes test
if 'attributes' in tests:
    t = tic()

    print('Running attributes test...')
    from optima import Project
    P = Project(spreadsheet='simple.xlsx')
    P.parsets[0].listattributes()

    done(t)







## Sensitivity test
if 'sensitivity' in tests:
    t = tic()

    print('Running sensitivity test...')
    from optima import demo
    
    P = demo(doplot=False)
    P.parset().updateprior() # Make sure it's up to date with the calibration
    results = P.sensitivity(orig='default', name='sensitivity', tosample=None, n=5, randseed=15000) # tosample=['initprev','force']
    
    if doplot:
        from optima import pygui
        pygui(results)
    
    done(t)






## Manual calibration test
if 'manualfit' in tests and doplot:
    t = tic()

    print('Running manual calibration test...')
    from optima import Project, manualfit
    
    P = Project(spreadsheet='generalized.xlsx')
    manualfit(project=P, orig='default')
    
    done(t)






## Autofit test
if 'autofit' in tests:
    t = tic()

    print('Running autofit test...')
    from optima import Project
    
    P = Project(spreadsheet='generalized.xlsx')
    P.autofit(name='autofit', orig='default', fitwhat=['force'], maxtime=None, maxiters=5) # Run automatic fitting
    
    if doplot:
        from optima import plotresults
        results1 = P.runsim('default', end=2015) # Generate results
        results2 = P.runsim('autofit', end=2015)
        plotresults(P.parsets['default'].getresults(), toplot=['prev-total', 'prev-population', 'numinci-population'])
        plotresults(P.parsets['autofit'].getresults(), toplot=['improvement', 'prev-total', 'prev-population', 'numinci-population'])
    
    done(t)







## Autofit test
if 'autofitmulti' in tests:
    t = tic()

    print('Running autofitmulti test...')
    from optima import Project
    
    P = Project(spreadsheet='generalized.xlsx')
    P.sensitivity(orig='default', name='sensitivity', n=5, span=0.5) # Create MC initialization
    P.runsim('sensitivity', end=2015) # Generate results
    P.autofit(name='autofit', orig='sensitivity', fitwhat=['force'], maxtime=None, maxiters=30, inds=None) # Run automatic fitting
    
    
    if doplot:
        from optima import plotresults
        plotresults(P.parsets['sensitivity'].getresults(), toplot=['prev-tot', 'numinci-sta'])
        plotresults(P.parsets['autofit'].getresults(), toplot=['improvement', 'prev-tot', 'numinci-sta'])
    
    done(t)












## Autofit test
if 'longfit' in tests:
    t = tic()

    print('Running long autofit test...')
    from optima import Project
    
    P = Project(spreadsheet='generalized.xlsx')
    P.autofit(name='autofit', orig='default', fitwhat=['init','popsize','force','const'], maxiters=1000, inds=None, verbose=2) # Run automatic fitting
    results1 = P.runsim('default', end=2015) # Generate results
    results2 = P.runsim('autofit', end=2015)
    
    if doplot:
        from optima import plotresults
        plotresults(results1, toplot=['prev-tot', 'prev-pops', 'numinci-pops'])
        plotresults(results2, toplot=['improvement', 'prev-tot', 'prev-pops', 'numinci-pops'])
    
    done(t)






## Debug autofit test
if 'debugautofit' in tests:
    t = tic()

    print('Running autofit debugging test...')
    from optima import Project, plotresults
    P = Project(spreadsheet='concentrated.xlsx')
    
    # Run automatic fitting
    P.autofit(name='autofit', orig='default', fitwhat='force', fitto='prev',  method='mse', maxiters=50, doplot=doplot, verbose=2) # Set doplot=True and verbose=4 to see full debugging information
    if doplot:
        plotresults(P.results['parset-default'])
        plotresults(P.results['parset-autofit'])
    
    done(t)



print('\n\n\nDONE: ran %i tests' % len(tests))
toc(T)