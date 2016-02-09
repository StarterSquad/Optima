"""
Test script to see if Optima works.

To use: comment out lines in the definition of 'tests' to not run those tests.

NOTE: for best results, run in interactive mode, e.g.

python -i tests.py

Version: 2016feb06
"""



## Define tests to run here!!!
tests = [
#'makeprograms',
#'compareoutcomes',
'reconcilepars',
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

T = tic()

##############################################################################
## The tests
##############################################################################






## Programs creation test
if 'makeprograms' in tests:
    t = tic()

    print('Running make programs test...')
    from optima import Project, Program, Programset, odict
    
    P = Project(spreadsheet='concentrated.xlsx')
    pops = P.data['pops']['short']

    # First set up some programs. Programs need to be initialized with a name. Often they will also be initialized with targetpars and targetpops
    HTC = Program(short='HTC',
                  name='HIV testing and counseling',
                  targetpars=[{'param': 'hivtest', 'pop': 'F 15+'},
                              {'param': 'hivtest', 'pop': 'M 15+'},
                              {'param': 'hivtest', 'pop': 'FSW'},
                              {'param': 'hivtest', 'pop': 'MSM'},
                              {'param': 'hivtest', 'pop': 'Clients'}],
                  targetpops=['F 15+','FSW', 'Clients', 'MSM'],
                  costcovdata = {'t':[2013],
                                 'cost':[1e6],
                                 'coverage':[3e5]})

    SBCC = Program(short='SBCC',
                   name='Social and behaviour change communication',
                   targetpars=[{'param': 'condcas', 'pop': ('F 15+','M 15+')},
                               {'param': 'hivtest', 'pop': 'F 15+'}],
                   targetpops=['F 15+']) 

    MGT = Program(short='MGT')

    ART = Program(short='ART',
                  targetpars=[{'param': 'numtx', 'pop': 'tot'}],
                  targetpops=pops,
                  criteria={'hivstatus': ['lt50', 'gt50', 'gt200', 'gt350'], 'pregnant': False})

    Adherence = Program(short='Adherence',
                  targetpars=[{'param': 'stoprate', 'pop': pop} for pop in ['F 15+', 'M 15+', 'FSW', 'Clients', 'MSM']],
                  targetpops=['F 15+', 'M 15+', 'FSW', 'Clients', 'MSM'],
                  criteria={'hivstatus': ['lt50', 'gt50', 'gt200', 'gt350'], 'pregnant': False})

    PMTCT = Program(short='PMTCT',
                  targetpars=[{'param': 'numtx', 'pop': 'tot'}, {'param': 'numpmtct', 'pop': 'tot'}],
                  targetpops=['tot'],
                  category='Care and treatment',
                  name='Prevention of mother-to-child transmission',
                  criteria = {'hivstatus': 'allstates', 'pregnant': True})
                  
    VMMC = Program(short='VMMC',
                  targetpars=[{'param': 'numcirc', 'pop': 'M 15+'},
                              {'param': 'numcirc', 'pop': 'MSM'},
                              {'param': 'numcirc', 'pop': 'Clients'},
                              {'param': 'numcirc', 'pop': 'PWID'}],
                  targetpops=['M 15+', 'MSM', 'Clients', 'PWID'],
                  category='Prevention',
                  name='Voluntary medical male circumcision',
                  criteria = {'hivstatus': 'allstates', 'pregnant': False})              
                  
    # Testing methods of program class
    # 1. Adding a target parameter to a program
    HTC.addtargetpar({'param': 'hivtest', 'pop': 'M 15+'})
    
    ## NOTE that adding a targeted parameter does NOT automatically add a targeted population! Do this separately, e.g.
    HTC.targetpops.append('M 15+')
        
    # 2. Removing a target parameter from a program (then readding it)
    HTC.rmtargetpar({'param': 'hivtest', 'pop': 'F 15+'})
    HTC.addtargetpar({'param': 'hivtest', 'pop': 'F 15+'})

    # 3. Add historical cost-coverage data point                         
    SBCC.addcostcovdatum({'t':2011,
                         'cost':2e7,
                         'coverage':8e5})
    SBCC.addcostcovdatum({'t':2014,
                         'cost':4e7,
                         'coverage':10e5})
    HTC.addcostcovdatum({'t':2015,
                         'cost':1e7,
                         'coverage':4e5})
    MGT.addcostcovdatum({'t':2015,
                         'cost':1e6})
    ART.addcostcovdatum({'t':2015,
                         'cost':1e7})
    PMTCT.addcostcovdatum({'t':2015,
                         'cost':1e7})
    VMMC.addcostcovdatum({'t':2015,
                         'cost':1e7})     
    Adherence.addcostcovdatum({'t':2015,
                         'cost':2e6})
                         
    # 4. Overwrite historical cost-coverage data point
    HTC.addcostcovdatum({'t':2013,
                         'cost':2e6,
                         'coverage':3e5},
                         overwrite=True)

    # 5. Remove historical cost-coverage data point - specify year only
    HTC.rmcostcovdatum(2013)

    # 6. Add parameters for defining cost-coverage function.
    HTC.costcovfn.addccopar({'saturation': (0.75,0.85),
                             't': 2013.0,
                             'unitcost': (30,40)})
                             
    HTC.costcovfn.addccopar({'t': 2016.0,
                             'unitcost': (25,35)})
                             
    HTC.costcovfn.addccopar({'t': 2017.0,
                             'unitcost': (30,35)})
                             
    SBCC.costcovfn.addccopar({'saturation': (0.4,0.5),
                              't': 2013.0,
                              'unitcost': (8,12)})

    ART.costcovfn.addccopar({'saturation': (0.8,0.9),
                              't': 2013.0,
                              'unitcost': (100,200)})

    PMTCT.costcovfn.addccopar({'saturation': (0.9,0.9),
                             't': 2016.0,
                             'unitcost': (100,200)})

    VMMC.costcovfn.addccopar({'saturation': (.5,.6),
                             't': 2016.0,
                             'unitcost': (15,25)})
                             
    Adherence.costcovfn.addccopar({'saturation': (.2,.3),
                                   't': 2016.0,
                                   'unitcost': (300,500)})
                             
    # 7. Overwrite parameters for defining cost-coverage function.
    HTC.costcovfn.addccopar({'t': 2016.0,
                             'unitcost': (20,30)},
                             overwrite=True)

    # 8. Remove parameters for defining cost-coverage function.
    HTC.costcovfn.rmccopar(2017)

    # 9. Get parameters for defining cost-coverage function for any given year (even if not explicitly entered).
    HTC.costcovfn.getccopar(2014)

    # 10. Get target population size
    HTC.gettargetpopsize(t=[2013,2015],parset=P.parsets['default'])

    # 11. Evaluate cost-coverage function to get coverage for a given year, spending amount and population size
    from numpy import linspace, array
    HTC.getcoverage(x=linspace(0,1e6,3),t=[2013,2015,2017],parset=P.parsets['default'],total=False,bounds=None)
    HTC.targetcomposition = {'Clients': array([ 0.01]),
                       'F 15+': array([ 0.3]),
                       'FSW': array([ 0.24]),
                       'M 15+': array([ 0.3]),
                       'MSM': [ 0.15]}
    
    HTC.getcoverage(x=[2e7],t=[2016],parset=P.parsets['default'],total=False)
    HTC.getbudget(x=linspace(0,1e6,3),t=[2013,2015,2017],parset=P.parsets['default'],proportion=False)

    # NB, if you want to evaluate it for a particular population size, can also do...
    HTC.costcovfn.evaluate(x=[1e6],popsize=[1e5],t=[2015],toplot=False)

    # 12. Plot cost-coverage function
    caption = 'Spending data includes all HTC spending. Global Fund spending on HTC in 2012 was $40,245. '\
                  'In the reporting period, a total of 676 MARPs received HTC services, which makes a cumulative '\
                  'total of 1,102 MARPs that received HTC including provision of results. Due to changes in '\
                  'the definition and focus of the indicator, PWID that received HTC in DST Centers and prisoners '\
                  'are included, both of them previously omitted in the reports.'
    plotoptions = {}
    plotoptions['caption'] = caption
    plotoptions['xupperlim'] = 2e9
    plotoptions['perperson'] = False

    if doplot:
        ART.plotcoverage(t=[2014,2015],parset=P.parsets['default'],doplot=doplot)
        HTC.plotcoverage(t=[2014,2015],parset=P.parsets['default'],plotoptions=plotoptions,doplot=doplot)
        PMTCT.plotcoverage(t=[2014,2015],parset=P.parsets['default'],doplot=doplot)
        VMMC.plotcoverage(t=[2014,2015],parset=P.parsets['default'],doplot=doplot)

    print('Running make programs set test...')

    # Initialise with or without programs
    R = Programset()
    R = Programset(programs=[HTC,SBCC,MGT,ART,PMTCT,VMMC,Adherence])

    # Testing methods of programset class
    # 1. Adding a program
#    R.addprograms(ART)

    # 2. Removing a program
#    R.rmprogram(ART) # Alternative syntax: R.rmprogram('ART')
    
    # 3. See which programs are optimizable
    R.optimizable() # True/False lists
    R.optimizableprograms() # odict of optimizable programs


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

    # 7. Get a dictionary of coverage levels corresponding to a dictionary of program allocations
    budget=odict({'SBCC':array([1e6,1.2e6,1.5e6]),
                  'HTC':array([1e7,1.2e7,1.5e7]),
                  'ART':array([1e7,1.2e7,1.5e7]),
                  'PMTCT':array([1e7,1.2e7,1.5e7]),
                  'VMMC':array([1e7,1.2e7,1.5e7]),
                  'Adherence':array([1e6,1.2e6,1.5e6]),
                  'MGT':array([2e5,3e5,3e5])})
            
    coverage=odict({'HTC': array([ 368122.94593941, 467584.47194668, 581136.7363055 ]),
              'MGT': None,
              'ART':array([1e5,1.2e5,1.5e5]),
              'PMTCT':array([1e3,1.2e3,1.5e3]),
              'VMMC':array([1e5,1.2e5,1.5e5]),
              'Adherence':array([1e6,1.2e6,1.5e6]),
              'SBCC': array([ 97615.90198599, 116119.80759447, 143846.76414342])})
              
    budget = budget.sort([p.short for p in R.programs.values()])
    coverage = coverage.sort([p.short for p in R.programs.values()])

    defaultbudget = R.getdefaultbudget()
    defaultcoverage = R.getdefaultcoverage(t=2015, parset=P.parsets['default'])

    R.getprogcoverage(budget=budget,
                      t=[2015,2016,2020],
                      parset=P.parsets['default'])
                        
    R.getprogbudget(coverage=coverage,
                      t=[2015,2016,2020],
                      parset=P.parsets['default'])
                        
    R.getpopcoverage(budget=budget,
                     t=[2015,2016,2020],
                     parset=P.parsets['default'])

    # 8. Add parameters for defining coverage-outcome function.
    R.covout['hivtest']['F 15+'].addccopar({'intercept': (0.25,0.35),
                                                    't': 2013.0,
                                                  'HTC': (0.85,0.95),
                                                 'SBCC': (0.35,0.45)})
                                                                                                    
    R.covout['hivtest']['M 15+'].addccopar({'intercept': (0.25,0.35),
                                                  't': 2016.0,
                                                  'HTC': (0.9,1.)})
                                                  
    R.covout['hivtest']['F 15+'].addccopar({'intercept': (0.25,0.35),
                                                    't': 2015.0,
                                                    'HTC': (0.75,0.85),
                                                    'SBCC':(0.4,0.5)})
                                                    
    R.covout['hivtest']['F 15+'].addccopar({'intercept': (0.35,0.45),
                                                    't': 2017.0,
                                                    'HTC': (0.8,0.85),
                                                    'SBCC':(0.6,0.65)})

    R.covout['hivtest']['Clients'].addccopar({'intercept': (0.35,0.45),
                                                    't': 2017.0,
                                                    'HTC': (0.8,0.85)})
                                                    
    R.covout['hivtest']['MSM'].addccopar({'intercept': (0.35,0.45),
                                                    't': 2017.0,
                                                    'HTC': (0.8,0.85)})
                                                    
    R.covout['hivtest']['FSW'].addccopar({'intercept': (0.35,0.45),
                                                    't': 2017.0,
                                                    'HTC': (0.8,0.85)})
                                                    
    R.covout['condcas'][('F 15+', 'M 15+')].addccopar({'intercept': (0.3,0.35), 
                                                    't': 2015.0,
                                                    'SBCC':(0.45,0.55)})
                                                    
    R.covout['numtx']['tot'].addccopar({'intercept': (100.0,150.0), 't': 2016.0})
    R.covout['numpmtct']['tot'].addccopar({'intercept': (100.0,150.0), 't': 2016.0})

    R.covout['numcirc']['MSM'].addccopar({'intercept': (0,0), 't': 2016.0})
    R.covout['numcirc']['Clients'].addccopar({'intercept': (0,0), 't': 2016.0})
    R.covout['numcirc']['PWID'].addccopar({'intercept': (0,0), 't': 2016.0})
    R.covout['numcirc']['M 15+'].addccopar({'intercept': (0,0), 't': 2016.0})

#    R.covout['stoprate']['MSM'].addccopar({'intercept': (0,0), 'Adherence': (.5,.6), 't': 2016.0})
#    R.covout['stoprate']['Clients'].addccopar({'intercept': (0,0), 'Adherence': (.5,.6), 't': 2016.0})
#    R.covout['stoprate']['FSW'].addccopar({'intercept': (0,0), 'Adherence': (.5,.6), 't': 2016.0})
#    R.covout['stoprate']['M 15+'].addccopar({'intercept': (0,0), 'Adherence': (.5,.6), 't': 2016.0})
#    R.covout['stoprate']['F 15+'].addccopar({'intercept': (0,0), 'Adherence': (.5,.6), 't': 2016.0})

    # 9. Overwrite parameters for defining coverage-outcome function.
    R.covout['hivtest']['F 15+'].addccopar({'intercept': (0.35,0.45),
                                                    't': 2015.0,
                                                    'HTC': (0.85,0.95),
                                                    'SBCC':(0.55,0.65)},
                                                    overwrite=True)

    # 10. Remove parameters for defining coverage-outcome function.
    R.covout['hivtest']['F 15+'].rmccopar(2017)
    
    # 11. Get parameters for defining cost-coverage function for any given year (even if not explicitly entered).
    R.covout['hivtest']['F 15+'].getccopar(2014)

    # 12. Get a dictionary of only the program-affected parameters corresponding to a dictionary of program allocations or coverage levels
    outcomes = R.getoutcomes(coverage=coverage,
                                t=[2015,2016,2020],
                                parset=P.parsets['default'])
    
    R.getoutcomes(defaultcoverage, t=2015, parset=P.parsets['default'])
    R.getoutcomes(t=2015, parset=P.parsets['default'])
            

    # 13. Get an odict of the ALL parameter values corresponding to a vector of program allocations
    P.addprogset(name='default', progset=R)
    P.runbudget(budget=budget, budgetyears=[2015,2016,2020], progsetname='default', parsetname='default')
    
    

    done(t)
    




## Project creation test
if 'compareoutcomes' in tests:
    comparison = P.progsets[0].compareoutcomes(parset=P.parsets[0], year=2016, doprint=True)
    done(t)



# Reconciliation test
if 'reconcilepars' in tests:
    import optima as op
    P = op.defaults.defaultproject('best')
    P.progsets[0].reconcile(parset=P.parsets[0], year=2016)
    



print('\n\n\nDONE: ran %i tests' % len(tests))
toc(T)