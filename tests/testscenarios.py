"""
Test scenarios

Version: 2016feb07
"""


## Define tests to run here!!!
tests = [
'standardscen',
#'maxcoverage',
#'maxbudget',
#'90-90-90'
#'VMMC'
]

##############################################################################
## Initialization -- same for every test script
##############################################################################

from optima import tic, toc, blank, pd # analysis:ignore

if 'doplot' not in locals(): doplot = True
showstats = True

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


## Standard scenario test
if 'standardscen' in tests:
    t = tic()

    print('Running standard scenarios test...')
    from optima import Parscen, Budgetscen, Coveragescen, defaults
    from numpy import array
    
    P = defaults.defaultproject('concentrated')
    pops = P.data['pops']['short']
    malelist = [i for i in range(len(pops)) if P.data['pops']['male'][i]]
    
    caspships = P.parsets['default'].pars[0]['condcas'].y.keys()
    
    
    ## Create a scenario
    thisname = 'Get lots of people on treatment'
    thisparset = 0
    thispar = 'numtx'
    thisfor = 0
    thisstartyear = 2016.
    thisstartval = P.parsets[thisparset].pars[0][thispar].interp(thisstartyear)[thisfor]
    thisendyear = 2020.
    thisendval = 100000.
    
    scenlist = [
        Parscen(name='Get lots of people on treatment',
                parsetname=0,
                pars=[{
                 'name': 'numtx',
                 'for': 'tot',
                 'startyear': 2016.,
                 'endyear': thisendyear,
                 'endval': thisendval
                 }]),

        Parscen(name='Imagine that no-one gets circumcised',
             parsetname='default',
             pars=[{
                 'name': 'propcirc',
                 'for': malelist,
                 'startyear': 2015,
                 'endyear': 2020,
                 'endval': 0.,
                 }]),

        Parscen(name='Increase numpmtct',
             parsetname='default',
             pars=[{
                 'name': 'numpmtct',
                 'for': 'tot',
                 'startyear': 2015.,
                 'endyear': 2020,
                 'endval': 0.9,
                 }]),

        Parscen(name='Full casual condom use',
             parsetname='default',
             pars=[{
                 'name': 'condcas',
                 'for': caspships,
                 'startyear': 2005,
                 'endyear': 2015,
                 'endval': 1.,
                 }]),

         Parscen(name='More casual acts',
              parsetname='default',
              pars=[{
                  'name': 'actscas',
                  'for': caspships,
                  'startyear': 2005,
                  'endyear': 2015,
                  'endval': 2.,
                  }]),

         Parscen(name='100% testing',
              parsetname='default',
              pars=[{
                  'name': 'hivtest',
                  'for': ['FSW', 'Clients', 'MSM', 'M 15+', 'F 15+'],
                  'startyear': 2000.,
                  'endyear': 2020,
                  'endval': 1.,
                  }]),

         Parscen(name='Increased STI prevalence in FSW',
              parsetname='default',
              pars=[{
                  'name': 'stiprev',
                  'for': 0,
                  'startyear': 2005.,
                  'endyear': 2015,
                  'endval': 0.8,
                  }]),

         Parscen(name='Get 50K people on OST',
              parsetname='default',
              pars=[{
                  'name': 'numost',
                  'for': 0,
                  'startyear': 2005.,
                  'endyear': 2015,
                  'endval': 50000,
                  }]),

         Budgetscen(name='Keep current investment in condom program',
              parsetname='default',
              progsetname='default',
              t=2016,
              budget={'Condoms': 1e7,
                           'FSW programs': 1e6,
                           'HTC':2e7,
                           'ART':1e6}),

         Budgetscen(name='Double investment in condom program',
              parsetname='default',
              progsetname='default',
              t=[2016,2020],
              budget={'Condoms': array([1e7,2e7]),
                           'FSW programs':array([1e6,1e6]),
                           'HTC':array([2e7,2e7]),
                           'ART':array([1e6,1e6])}),

         Coveragescen(name='A million people covered by the condom program',
              parsetname='default',
              progsetname='default',
              t=[2016,2020],
              coverage={'Condoms': array([285706.,1e6]),
                           'FSW programs':array([15352.,15352.]),
                           'HTC':array([1332862.,1332862.]),
                           'ART':array([3324.,3324.])}),

         Budgetscen(name='Double investment in ART, HTC and OST',
              parsetname='default',
              progsetname='default',
              t=[2016,2018,2020],
              budget={'Condoms': array([1e7,1e7,1e7]),
                           'FSW programs':array([1e6,1e6,1e6]),
                           'HTC':array([2e7,3e7,4e7]),
                           'ART':array([1e6,1.5e6,2e6])}),

         Budgetscen(name='Test some progs only',
              parsetname='default',
              progsetname='default',
              t=2016,
              budget={'Condoms': 1e7,
                           'ART':1e6})

        ]
    
    # Store these in the project
    P.addscenlist(scenlist)

    # Run the scenarios
    P.runscenarios() 
     
    if doplot:
        from optima import pygui
        pygui(P.results[-1], toplot='default')

    if showstats:
        from optima import Settings, findinds
        from numpy import arange
        settings = Settings()
        tvec = arange(settings.start,settings.end+settings.dt,settings.dt)
        yr = 2020
        blank()
        for scenno, scen in enumerate([scen for scen in P.scens.values() if scen.active]):
            output = '===================================\n'
            output += scen.name
            output += '\n'           
            output += 'PLHIV: %s\n' % (P.results[-1].raw[scenno][0]['people'][settings.allplhiv,:,findinds(tvec,yr)].sum(axis=(0,1)))
            output += 'Prop aware: %s\n' % (P.results[-1].raw[scenno][0]['people'][settings.alldx,:,findinds(tvec,yr)].sum(axis=(0,1))/P.results[-1].raw[scenno][0]['people'][settings.allplhiv,:,findinds(tvec,yr)].sum(axis=(0,1)))
            output += 'Number treated: %s\n' % (P.results[-1].raw[scenno][0]['people'][settings.alltx,:,findinds(tvec,yr)].sum(axis=(0,1)))
            output += 'Prop treated: %s\n' % (P.results[-1].raw[scenno][0]['people'][settings.alltx,:,findinds(tvec,yr)].sum(axis=(0,1))/P.results[-1].raw[scenno][0]['people'][settings.allplhiv,:,findinds(tvec,yr)].sum(axis=(0,1)))
            print output


    done(t)



## 90-90-90 scenario test
if '90-90-90' in tests:
    t = tic()

    print('Running standard scenarios test...')
    from optima import Parscen, defaults, pygui, findinds, plotpeople
    from numpy import nan
    
    P = defaults.defaultproject('best')
    P.runsim(debug=True)
    
    pops = P.data['pops']['short']
    
    startyear = 2014.
    endyear = 2020.
    res_startind = findinds(P.results[-1].tvec, startyear)
    res_endind = findinds(P.results[-1].tvec, endyear)
    
    start_propdx = P.results[-1].main['numdiag'].tot[0][res_startind]/P.results[-1].main['numplhiv'].tot[0][res_startind]
    start_propincare = P.results[-1].main['numincare'].tot[0][res_startind]/P.results[-1].main['numdiag'].tot[0][res_startind]
    start_proptx = P.results[-1].main['numtreat'].tot[0][res_startind]/P.results[-1].main['numincare'].tot[0][res_startind]
    start_propsupp = P.results[-1].main['numsuppressed'].tot[0][res_startind]/P.results[-1].main['numtreat'].tot[0][res_startind]
    

    ## Define scenarios
    scenlist = [
        Parscen(name='Current conditions',
                parsetname='default',
                pars=[]),

         Parscen(name='90-90-90',
              parsetname='default',
              pars=[
              {'name': 'propdx',
              'for': ['tot'],
              'startyear': startyear,
              'endyear': endyear,
              'startval': start_propdx,
              'endval': .9,
              },
              
              {'name': 'propcare',
              'for': ['tot'],
              'startyear': startyear,
              'endyear': endyear,
              'startval': start_propincare,
              'endval': .9,
              },
              
              {'name': 'proptx',
              'for': ['tot'],
              'startyear': startyear,
              'endyear': endyear,
              'startval': start_proptx,
              'endval': .9,
              },
              
              {'name': 'propsupp',
              'for': ['tot'],
              'startyear': startyear,
              'endyear': endyear,
              'startval': start_propsupp,
              'endval': .9,
              },
                ]),
                
         Parscen(name='Increase numtx',
              parsetname='default',
              pars=[
              {'name': 'numtx',
              'for': ['tot'],
              'startyear': startyear,
              'endyear': 2030.,
              'startval': 48100.,
              'endval': 68000.,
              }]),
                
         Parscen(name='Constant numtx',
              parsetname='default',
              pars=[
              {'name': 'proptx',
              'for': ['tot'],
              'startyear': startyear,
              'startval': nan,
#              'endyear': 2030.,
#              'endval': 48100.,
              }]),
                
        ]

    # Store these in the project
    P.addscenlist(scenlist)
    P.scens[0].active = True # Turn off 90-90-90 scenario

    # Run the scenarios
    P.runscenarios(debug=True) 

    if showstats:
        blank()
        for scind, scen in enumerate([scen for scen in P.scens.values() if scen.active]):
            end_plhiv = P.results[-1].main['numplhiv'].tot[scind][res_endind]
            end_propdx = P.results[-1].main['numdiag'].tot[scind][res_endind]/P.results[-1].main['numplhiv'].tot[scind][res_endind]
            end_propincare = P.results[-1].main['numincare'].tot[scind][res_endind]/P.results[-1].main['numdiag'].tot[scind][res_endind]
            end_proptx = P.results[-1].main['numtreat'].tot[scind][res_endind]/P.results[-1].main['numincare'].tot[scind][res_endind]
            end_propsupp = P.results[-1].main['numsuppressed'].tot[scind][res_endind]/P.results[-1].main['numtreat'].tot[scind][res_endind]
            output = '===================================\n'
            output += scen.name
            output += '\nOutcomes in Year %i\n' % (endyear)
            output += 'PLHIV: %s\n' % (end_plhiv)
            output += 'Prop aware: %s\n' % (end_propdx)
            output += 'Prop in care: %s\n' % (end_propincare)
            output += 'Prop treated: %s\n' % (end_proptx)
            output += 'Prop suppressed: %s\n' % (end_propsupp)
            print output
   
     
    if doplot:
#        ppl = P.results[-1].raw['90-90-90'][0]['people']
#        plotpeople(P, ppl)
        pygui(P.results[-1], toplot='cascade')

    done(t)





#################################################################################################################
## Coverage
#################################################################################################################

if 'maxcoverage' in tests:
    t = tic()

    print('Running maximum coverage scenario test...')
    from optima import Coveragescen, Parscen, defaults, dcp
    from numpy import array
    
    ## Set up default project
    P = defaults.defaultproject('generalized')
    
    ## Define scenarios
    defaultbudget = P.progsets['default'].getdefaultbudget()
    maxcoverage = dcp(defaultbudget) # It's just an odict, though I know this looks awful
    for key in maxcoverage: maxcoverage[key] = array([maxcoverage[key]+1e9])
    scenlist = [
        Parscen(name='Current conditions', parsetname='default', pars=[]),
        Coveragescen(name='Full coverage', parsetname='default', progsetname='default', t=[2016], coverage=maxcoverage),
        ]
    
    # Run the scenarios
    P.addscenlist(scenlist)
    P.runscenarios() 
     
    if doplot:
        from optima import pygui
        pygui(P.results[-1], toplot='default')






## Set up project etc.
if 'maxbudget' in tests:
    t = tic()

    print('Running maximum budget scenario test...')
    from optima import Budgetscen, defaults, dcp
    from numpy import array
    
    ## Set up default project
    P = defaults.defaultproject('generalized')
    
    ## Define scenarios
    defaultbudget = P.progsets['default'].getdefaultbudget()
    maxbudget = dcp(defaultbudget)
    for key in maxbudget: maxbudget[key] += 1e14
    zerobudget = dcp(defaultbudget)
    for key in zerobudget: zerobudget[key] = array([0.]) # Alternate way of setting to zero   
    scenlist = [
        Budgetscen(name='Current conditions', parsetname='default', progsetname='default', t=[2016], budget=defaultbudget),
        Budgetscen(name='Unlimited spending', parsetname='default', progsetname='default', t=[2016], budget=maxbudget),
        Budgetscen(name='Zero spending', parsetname='default', progsetname='default', t=[2016], budget=zerobudget),
        ]
    
    # Run the scenarios
    P.addscenlist(scenlist)
    P.runscenarios() 
     
    if doplot:
        from optima import pygui, plotpars
        pygui(P.results[-1], toplot='default')
        apd = plotpars([scen.scenparset.pars[0] for scen in P.scens.values()])




## Set up project etc.
if 'VMMC' in tests:
    t = tic()

    print('Running VMMC scenario test...')
    from optima import Parscen, Budgetscen, findinds, defaults
    
    P = defaults.defaultproject('generalized')
    pops = P.data['pops']['short']

    malelist = findinds(P.data['pops']['male'])
    caspships = P.parsets['default'].pars[0]['condcas'].y.keys()
    
    ## Define scenarios
    scenlist = [
        Parscen(name='Current conditions',
                parsetname='default',
                pars=[]),

        Parscen(name='Imagine that no-one gets circumcised',
             parsetname='default',
             pars=[{'endval': 0.2,
                'endyear': 2020,
                'name': 'propcirc',
                'for': malelist,
                'startval': .85,
                'startyear': 2015.2}]),
        
        Budgetscen(name='Default budget',
              parsetname='default',
              progsetname='default',
              t=2015,
              budget=P.progsets['default'].getdefaultbudget()),

         Budgetscen(name='Scale up VMMC program',
              parsetname='default',
              progsetname='default',
              t=2016,
              budget={'VMMC': 1e8}),

        ]
    
    # Store these in the project
    P.addscenlist(scenlist)
    
    # Run the scenarios
    P.runscenarios()
     
    if doplot:
        from optima import pygui, plotpeople, plotpars
        ppl1 = P.results[-1].raw['Scale up VMMC program'][0]['people']
        ppl2 = P.results[-1].raw['Imagine that no-one gets circumcised'][0]['people']
        plotpeople(P, ppl1, start=0, end=None, pops=[-2], animate=False)
        apd = plotpars([scen.scenparset.pars[0] for scen in P.scens.values()])
        pygui(P.results[-1], toplot='default')
        

    done(t)
