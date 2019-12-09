'''
Define classes and functions for handling scenarios and ICERs.

Version: 2017jun03
'''

## Imports
from numpy import append, array, inf
from optima import OptimaException, Link, Multiresultset # Core classes/functions
from optima import dcp, today, odict, printv, findinds, defaultrepr, getresults, vec2obj, isnumber, uuid, promotetoarray # Utilities

class Scen(object):
    ''' The scenario base class -- not to be used directly, instead use Parscen or Progscen '''
    def __init__(self, name=None, parsetname=-1, progsetname=-1, t=None, active=True):
        self.uid = uuid()
        self.name = name
        self.parsetname  = parsetname
        self.progsetname = progsetname
        self.t = t
        self.active = active
        self.resultsref = None
        self.scenparset = None # Store the actual parset generated
    
    def __repr__(self):
        ''' Print out useful information when called'''
        output = defaultrepr(self)
        return output

    def getresults(self):
        ''' Returns the results '''
        if self.resultsref is not None and self.projectref() is not None:
            results = getresults(project=self.projectref(), pointer=self.resultsref)
            return results
        else:
            print('WARNING, no results associated with this scenario')
    
    

class Parscen(Scen):
    ''' An object for storing a single parameter scenario '''
    def __init__(self, pars=None, **defaultargs):
        Scen.__init__(self, **defaultargs)
        if pars is None: pars = []
        self.pars = pars



class Progscen(Scen):
    ''' The program scenario base class -- not to be used directly, instead use Budgetscen or Coveragescen '''
    def __init__(self, progsetname=-1, **defaultargs):
        Scen.__init__(self, **defaultargs)
        self.progsetname = progsetname # Programset


class Budgetscen(Progscen):
    ''' Stores a single budget scenario. Initialised with a budget. Coverage added during makescenarios()'''
    def __init__(self, budget=None, **defaultargs):
        Progscen.__init__(self, **defaultargs)
        self.budget = budget


class Coveragescen(Progscen):
    ''' Stores a single coverage scenario. Initialised with a coverage. Budget added during makescenarios()'''
    def __init__(self, coverage=None, **defaultargs):
        Progscen.__init__(self, **defaultargs)
        self.coverage = coverage


def runscenarios(project=None, verbose=2, name=None, defaultparset=-1, debug=False, nruns=None, base=None, ccsample=None, randseed=None, separateoutput=None, **kwargs):
    """
    Run all the scenarios.
    Version: 2017aug15
    """
    
    printv('Running scenarios...', 1, verbose)
    
    # Make sure scenarios exist
    if project is None: raise OptimaException('First argument to runscenarios() must be a project')
    if len(project.scens)==0:  # Create scenario list if not existing
        baselinescen = baselinescenario(project.parsets[defaultparset], verbose=verbose)
        project.addscens(baselinescen)
    scenlist = [scen for scen in project.scens.values() if scen.active==True]
    nscens = len(scenlist)
    
    # Handle inputs
    if ccsample is None: ccsample = 'best' 
    if nruns is None:    nruns = 1
    if base is None:     base = 0
    
    # Convert the list of scenarios to the actual parameters to use in the model
    scenparsets = makescenarios(project=project, scenlist=scenlist, ccsample=ccsample, randseed=randseed, verbose=verbose)        

    # Run scenarios
    allresults = []
    for scenno, scen in enumerate(scenparsets):
        printv('Running scenario "%s" (%i/%i)...' % (scen, scenno+1, nscens), 2, verbose)
        scenparset = scenparsets[scen]
        project.scens[scenno].scenparset = scenparset # Copy into scenarios objects
        

        # Items specific to program (budget or coverage) scenarios
        budget = scenlist[scenno].budget if isinstance(scenlist[scenno], Progscen) else None
        coverage = scenlist[scenno].coverage if isinstance(scenlist[scenno], Progscen) else None
        budgetyears = scenlist[scenno].t if isinstance(scenlist[scenno], Progscen) else None
        progsetname = scenlist[scenno].progsetname if isinstance(scenlist[scenno], Progscen) else None

        # Run model and add results
        result = project.runsim(pars=scenparset.pars, name=scenlist[scenno].parsetname, progsetname=progsetname, budget=budget, coverage=coverage, budgetyears=budgetyears, verbose=0, debug=debug, resultname=project.name+'-scenarios', addresult=False, n=nruns, **kwargs)

        result.name = scenlist[scenno].name # Give a name to these results so can be accessed for the plot legend
        allresults.append(result) 
        printv('... completed scenario: %i/%i' % (scenno+1, nscens), 3, verbose)
    
    if name is None: name='scenarios'

    multires = Multiresultset(resultsetlist=allresults, name=name)
    for scen in scenlist: scen.resultsref = multires.uid # Copy results into each scenario that's been run
    scenres = odict()
    scenres[name] = multires
    if separateoutput: return allresults # Return the individual parameter sets
    else:              return scenres





def makescenarios(project=None, scenlist=None, verbose=2, ccsample=False, randseed=None):
    """ Convert dictionary of scenario parameters into parset to model parameters """

    scenparsets = odict()
    for scenno, scen in enumerate(scenlist):
        
        try: 
            thisparset = dcp(project.parsets[scen.parsetname])
            thisparset.projectref = Link(project) # Replace copy of project with pointer -- TODO: improve logic
        except: raise OptimaException('Failed to extract parset "%s" from this project:\n%s' % (scen.parsetname, project))
        npops = len(thisparset.popkeys)

        if isinstance(scen,Parscen):
            if scenlist[scenno].pars is None: scenlist[scenno].pars = [] # Turn into empty list instead of None
            for scenpar in scenlist[scenno].pars: # Loop over all parameters being changed

                # Get the parameter object
                thispar = thisparset.pars[scenpar['name']]

                # Parse inputs to figure out which population(s) are affected
                if type(scenpar['for'])==tuple: # If it's a partnership...
                    if not scenpar['for'] in thispar.keys():
                        errormsg = 'Partnership %s not associated with parameter %s' % (scenpar['for'],thispar.short)
                        raise OptimaException(errormsg)
                    else: pops = [scenpar['for']] 

                elif isnumber(scenpar['for']): #... if its asingle  population.
                    pops = [range(npops)] if scenpar['for'] > npops else [scenpar['for']]

                elif type(scenpar['for']) in [list, type(array([]))]: #... if its list of population.
                    pops = scenpar['for']

                elif type(scenpar['for'])==str: 
                    if not scenpar['for'] in thispar.keys():
                        errormsg = 'Population %s not associated with parameter %s' % (scenpar['for'],thispar.short)
                        raise OptimaException(errormsg)
                    else: pops = [scenpar['for']] 

                else: 
                    errormsg = 'Unrecognized population or partnership type: %s' % scenpar['for']
                    raise OptimaException(errormsg)

                # Find last good value 
                last_t = scenpar['startyear'] - project.settings.dt # Last timestep before the scenario starts
                last_y = thispar.interp(tvec=last_t, dt=project.settings.dt, asarray=False, sample=False) # Find what the model would get for this value

                # Loop over populations
                for pop in pops:

                    # Get the index of the population
                    if isnumber(pop): popind = pop
                    else: popind = thispar.keys().index(pop)
                    
                    # Find or set new value 
                    if scenpar.get('startval') is not None:
                        this_y = promotetoarray(scenpar['startval']) # Use supplied starting value if there is one
                    else:
                        if int(thispar.fromdata): # If it's a regular parameter made from data, we get the default start value from the data
                            this_y = thispar.interp(tvec=scenpar['startyear'], sample=False)[popind] # Find what the model would get for this value
                        else:
                            this_y = inf # Another special value, indicating this should be filled in to the maximum                 
                            fixproppar = thisparset.pars['fix'+scenpar['name']] # Pull out e.g. fixpropdx
                            fixproppar.t = min(fixproppar.t, scenpar['startyear']) # Reset start year to the lower of these
                    
                    # Remove years after the last good year
                    if last_t < max(thispar.t[popind]):
                        thispar.t[popind] = thispar.t[popind][findinds(thispar.t[popind] <= last_t)]
                        thispar.y[popind] = thispar.y[popind][findinds(thispar.t[popind] <= last_t)]
                    
                    # Append the last good year, and then the new years
                    thispar.t[popind] = append(thispar.t[popind], last_t)
                    thispar.y[popind] = append(thispar.y[popind], last_y[popind]) 
                    thispar.t[popind] = append(thispar.t[popind], scenpar['startyear'])
                    thispar.y[popind] = append(thispar.y[popind], this_y) 
                    
                    # Add end year values if supplied
                    if scenpar.get('endyear') is not None: 
                        thispar.t[popind] = append(thispar.t[popind], scenpar['endyear'])
                        thispar.y[popind] = append(thispar.y[popind], scenpar['endval'])
                    
                    if len(thispar.t[popind])!=len(thispar.y[popind]):
                        raise OptimaException('Parameter lengths must match (t=%i, y=%i)' % (len(thispar.t), len(thispar.y)))
                        
        elif isinstance(scen,Progscen):

            try: thisprogset = dcp(project.progsets[scen.progsetname])
            except: raise OptimaException('Failed to extract progset "%s" from this project:\n%s' % (scen.progset, project))
            
            try: results = project.parsets[scen.parsetname].getresults() # See if there are results already associated with this parset
            except: results = None

            scen.t = promotetoarray(scen.t)
            
            if isinstance(scen, Budgetscen):
                
                # If the budget has been passed in as a vector, convert it to an odict & sort by program names
                tmpbudget = dcp(thisprogset.getdefaultbudget())
                if isinstance(scen.budget, list) or isinstance(scen.budget,type(array([]))):
                    scen.budget = vec2obj(orig=tmpbudget, newvec=scen.budget) # It seems to be a vector: convert to odict
                if not isinstance(scen.budget,dict): raise OptimaException('Currently only accepting budgets as dictionaries.')
                if not isinstance(scen.budget,odict): scen.budget = odict(scen.budget)

                # Update, ensuring a consistent number of programs, using defaults where not provided 
                tmpbudget.update(scen.budget)
                scen.budget = tmpbudget
                
                # Ensure budget values are lists
                for budgetkey, budgetentry in scen.budget.iteritems():
                    if isnumber(budgetentry):
                        scen.budget[budgetkey] = [budgetentry]
                
                # Figure out coverage
                scen.coverage = thisprogset.getprogcoverage(budget=scen.budget, t=scen.t, parset=thisparset, results=results, sample=ccsample)

            elif isinstance(scen, Coveragescen):
                
                # If the coverage levels have been passed in as a vector, convert it to an odict & sort by program names
                tmpbudget = dcp(thisprogset.getdefaultbudget())
                tmpcoverage = thisprogset.getprogcoverage(budget=tmpbudget, t=2000, parset=thisparset) # NOTE, IT DOES NOT MATTER THE VALUE OF t YOU USE HERE!!!

                if isinstance(scen.coverage, list) or isinstance(scen.coverage, type(array([]))):
                    scen.coverage = vec2obj(scen.progset.getdefaultbuget(), newvec=scen.coverage) # It seems to be a vector: convert to odict -- WARNING, super dangerous!!
                if not isinstance(scen.coverage,dict): raise OptimaException('Currently only accepting coverage as dictionaries.')
                if not isinstance(scen.coverage,odict): scen.coverage = odict(scen.coverage)

                # Update, ensuring a consistent number of programs, using defaults where not provided 
                tmpcoverage.update(scen.coverage)
                scen.coverage = tmpcoverage

                # Ensure coverage level values are lists
                for covkey, coventry in scen.coverage.iteritems():
                    if isnumber(coventry):
                        scen.coverage[covkey] = [coventry]

                # Figure out coverage
                scen.budget = thisprogset.getprogbudget(coverage=scen.coverage, t=scen.t, parset=thisparset, results=results)

            # Create parameter dictionary
            thisparsdict = thisprogset.getpars(coverage=scen.coverage, t=scen.t, parset=thisparset, results=results, sample=ccsample)
            scen.pars = thisparsdict
            thisparset.pars = thisparsdict
            
        else: 
            errormsg = 'Unrecognized program scenario type.'
            raise OptimaException(errormsg)
            
        thisparset.modified = today()
        thisparset.name = scen.name
        scenparsets[scen.name] = thisparset
        
    return scenparsets





def baselinescenario(parset=None, verbose=2):
    """ Define the baseline scenario -- "Baseline" by default """
    if parset is None: raise OptimaException('You need to supply a parset to generate default scenarios')
    
    scenlist = [Parscen()]
    
    ## Current conditions
    scenlist[0].name = 'Baseline'
    scenlist[0].parset = parset
    scenlist[0].pars = [] # No changes
    
    return scenlist



def setparscenvalues(parset=None, parname=None, forwhom=None, startyear=None, verbose=2):
    """ Set the values of a parameter scenario """
    if parset is None: raise OptimaException('You need to supply a parset to generate default scenarios')
    
    if parname is None: raise OptimaException('Please supply a parameter')
    
    ## Generate dictionary
    if parset.pars[parname].fromdata: # If it's a regular parameter made from data, we get the default start value from the data
        if startyear is None: startyear = parset.pars[parname].t[forwhom][-1]
        startval = parset.pars[parname].interp(startyear,asarray=False)[forwhom][0]
    else: # Otherwise, give up -- we can't predict a proportion until the model is run
        if startyear is None: startyear = parset.projectref().settings.now
        startval = None

    
    return {'startval':startval,'startyear':startyear}



def defaultscenarios(project=None, which=None, startyear=2020, endyear=2025, parset=-1, progset=-1, dorun=True, doplot=True, **kwargs):
    '''
    Add default scenarios to a project...examples include min-max budgets and 90-90-90.
    Keyword arguments are passed to runscenarios().
    
    '''
    
    if which is None: which = 'budgets'
    
    if which=='budgets':
        parsetname = 'default-scenarios'
        project.copyparset(orig=parset, new=parsetname)
        project.parsets['default-scenarios'].fixprops(False) # Ensure they're not fixed
        defaultbudget = project.progsets[progset].getdefaultbudget()
        maxbudget = dcp(defaultbudget)
        nobudget = dcp(defaultbudget)
        for key in maxbudget: maxbudget[key] += project.settings.infmoney
        for key in nobudget: nobudget[key] *= 1e-6
        scenlist = [
            Parscen(   name='Baseline',         parsetname=parsetname, pars=[]),
            Budgetscen(name='Zero budget',      parsetname=parsetname, progsetname=0, t=[startyear], budget=nobudget),
            Budgetscen(name='Baseline budget',  parsetname=parsetname, progsetname=0, t=[startyear], budget=defaultbudget),
            Budgetscen(name='Unlimited budget', parsetname=parsetname, progsetname=0, t=[startyear], budget=maxbudget),
            ]
    
    # WARNING, this may not entirely work
    elif which=='90-90-90' or which=='909090':
        scenlist = [
            Parscen(name='Baseline', parsetname=0, pars=[]),
            Parscen(name='90-90-90',
                  parsetname=0,
                  pars=[
                  {'name': 'propdx',
                  'for': ['tot'],
                  'startyear': startyear,
                  'endyear': endyear,
                  'startval': None,
                  'endval': 0.9,
                  },
                  
                  {'name': 'proptx',
                  'for': ['tot'],
                  'startyear': startyear,
                  'endyear': endyear,
                  'startval': None,
                  'endval': 0.9,
                  },
                  
                  {'name': 'propsupp',
                  'for': ['tot'],
                  'startyear': startyear,
                  'endyear': endyear,
                  'startval': None,
                  'endval': .9,
                  },
                  ]),]
    else:
        errormsg = 'Default scenario options are "budgets" and "90-90-90", not %s' % which
        raise OptimaException(errormsg)

    
    # Run the scenarios
    project.addscens(scenlist)
    if dorun: project.runscenarios(**kwargs)
    if doplot: 
        from optima import pygui
        pygui(project)
    return None # Can get it from project.scens