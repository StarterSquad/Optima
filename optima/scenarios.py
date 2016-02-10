## Imports
from numpy import append, array
from optima import OptimaException, dcp, today, odict, printv, findinds, runmodel, Multiresultset, defaultrepr, getresults, vec2obj, isnumber


class Scen(object):
    ''' The scenario base class -- not to be used directly, instead use Parscen or Progscen '''
    def __init__(self, name=None, parsetname=None, t=None, active=True):
        self.name = name
        self.parsetname = parsetname
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
        if self.resultsref is not None and self.project is not None:
            results = getresults(project=self.project, pointer=self.resultsref)
            return results
        else:
            print('WARNING, no results associated with this scenario')
    
    

class Parscen(Scen):
    ''' An object for storing a single parameter scenario '''
    def __init__(self, pars=None, **defaultargs):
        Scen.__init__(self, **defaultargs)
        self.pars = pars



class Progscen(Scen):
    ''' The program scenario base class -- not to be used directly, instead use Budgetscen or Coveragescen '''
    def __init__(self, progsetname=None, **defaultargs):
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


def runscenarios(project=None, verbose=2, defaultparset=0):
    """
    Run all the scenarios.
    Version: 2016jan22 by cliffk
    """
    
    printv('Running scenarios...', 1, verbose)
    
    # Make sure scenarios exist
    if project is None: raise OptimaException('First argument to runscenarios() must be a project')
    if len(project.scens)==0:  # Create scenario list if not existing
        defaultscens = defaultscenarios(project.parsets[defaultparset], verbose=verbose)
        project.addscenlist(defaultscens)
    scenlist = [scen for scen in project.scens.values() if scen.active==True]
    nscens = len(scenlist)
    
    # Convert the list of scenarios to the actual parameters to use in the model
    scenparsets = makescenarios(project=project, scenlist=scenlist, verbose=verbose)

    # Run scenarios
    allresults = []
    for scenno, scen in enumerate(scenparsets):
        scenparset = scenparsets[scen]
        project.scens[scenno].scenparset = scenparset # Copy into scenarios objects

        # Items specific to program (budget or coverage) scenarios
        budget = scenlist[scenno].budget if isinstance(scenlist[scenno], Progscen) else None
        coverage = scenlist[scenno].coverage if isinstance(scenlist[scenno], Progscen) else None
        budgetyears = scenlist[scenno].t if isinstance(scenlist[scenno], Progscen) else None
        progset = project.progsets[scenlist[scenno].progsetname] if isinstance(scenlist[scenno], Progscen) else None

        # Run model and add results
        result = runmodel(pars=scenparset.pars[0], parset=scenparset, progset=progset, project=project, budget=budget, coverage=coverage, budgetyears=budgetyears, verbose=0)
        result.name = scenlist[scenno].name # Give a name to these results so can be accessed for the plot legend
        allresults.append(result) 
        printv('Scenario: %i/%i' % (scenno+1, nscens), 2, verbose)
    
    multires = Multiresultset(resultsetlist=allresults, name='scenarios')
    for scen in scenlist: scen.resultsref = multires.uid # Copy results into each scenario that's been run
    
    return multires





def makescenarios(project=None, scenlist=None, verbose=2):
    """ Convert dictionary of scenario parameters into parset to model parameters """

    scenparsets = odict()
    for scenno, scen in enumerate(scenlist):
        
        try: thisparset = dcp(project.parsets[scen.parsetname])
        except: raise OptimaException('Failed to extract parset "%s" from this project:\n%s' % (scen.parsetname, project))
        thisparset.modified = today()
        thisparset.name = scen.name
        npops = len(thisparset.popkeys)

        if isinstance(scen,Parscen):
            for pardictno in range(len(thisparset.pars)): # Loop over all parameter sets
                if scenlist[scenno].pars is None: scenlist[scenno].pars = [] # Turn into empty list instead of None
                for scenpar in scenlist[scenno].pars: # Loop over all parameters being changed
                    thispar = thisparset.pars[pardictno][scenpar['name']]
                    if type(scenpar['for'])==tuple: # If it's a partnership...
                        par2 = (scenpar['for'][1],scenpar['for'][0])
                        pops = [scenpar['for'], par2] # This is confusing - for partnership parameters, pops is a list of the two different partnership orderings.
                    elif type(scenpar['for'])==int: #... if its a population.
                        pops = range(npops) if scenpar['for'] > npops else [scenpar['for']]
                    elif type(scenpar['for']) in [list, type(array([]))]: #... if its a population.
                        pops = scenpar['for']
                    elif scenpar['for']=='tot': #... if its a population.
                        pops = [scenpar['for']]
                    else: 
                        errormsg = 'Unrecognized population or partnership type: %s' % scenpar['for']
                        raise OptimaException(errormsg)

                    # Find last good value # WARNING, FIX!!!!
                    last_t = scenpar['startyear'] - project.settings.dt # Last timestep before the scenario starts
                    last_y = thispar.interp(tvec=last_t, dt=project.settings.dt, asarray=False) # Find what the model would get for this value

                    # Loop over populations
                    for pop in pops:

                        # Remove years after the last good year
                        if isnumber(pop): pass
                        else: pop = tuple([str(p) for p in pop]) if isinstance(pop, (list, tuple)) else str(pop)
                        if last_t < max(thispar.t[pop]):
                            thispar.t[pop] = thispar.t[pop][findinds(thispar.t[pop] <= last_t)]
                            thispar.y[pop] = thispar.y[pop][findinds(thispar.t[pop] <= last_t)]
                        
                        # Append the last good year, and then the new years
                        thispar.t[pop] = append(thispar.t[pop], last_t)
                        thispar.y[pop] = append(thispar.y[pop], last_y[pop]) 
                        thispar.t[pop] = append(thispar.t[pop], scenpar['startyear'])
                        thispar.y[pop] = append(thispar.y[pop], scenpar['startval']) 
                        
                        # Add end year values if supplied
                        if scenpar['endyear']: 
                            thispar.t[pop] = append(thispar.t[pop], scenpar['endyear'])
                            thispar.y[pop] = append(thispar.y[pop], scenpar['endval'])
                        
                        if len(thispar.t[pop])!=len(thispar.y[pop]):
                            raise OptimaException('Parameter lengths must match (t=%i, y=%i)' % (len(thispar.t), len(thispar.y)))
                        
                    thisparset.pars[pardictno][scenpar['name']] = thispar # WARNING, not sure if this is needed???
    
        elif isinstance(scen,Progscen):

            try: thisprogset = dcp(project.progsets[scen.progsetname])
            except: raise OptimaException('Failed to extract progset "%s" from this project:\n%s' % (scen.progset, project))
            
            try: results = project.parsets[scen.parsetname].getresults() # See if there are results already associated with this parset
            except:
                results = None

            if isnumber(scen.t): scen.t = [scen.t]
            
            if isinstance(scen, Budgetscen):
                
                # If the budget has been passed in as a vector, convert it to an odict & sort by program names
                tmpbudget = dcp(thisprogset.getdefaultbudget())
                if isinstance(scen.budget, list) or isinstance(scen.budget,type(array([]))):
                    scen.budget = vec2obj(orig=tmpbudget, newvec=scen.budget) # It seems to be a vector: convert to odict
                if not isinstance(scen.budget,dict): raise OptimaException('Currently only accepting budgets as dictionaries.')
                if not isinstance(scen.budget,odict): scen.budget = odict(scen.budget)

                # Update, ensuring a consistent number of programs, using defaults where not provided -- WARNING, ugly
                tmpbudget.update(scen.budget)
                scen.budget = tmpbudget
                
                # Ensure budget values are lists
                for budgetkey, budgetentry in scen.budget.iteritems():
                    if isnumber(budgetentry):
                        scen.budget[budgetkey] = [budgetentry]
                
                # Figure out coverage
                scen.coverage = thisprogset.getprogcoverage(budget=scen.budget, t=scen.t, parset=thisparset, results=results)

            elif isinstance(scen, Coveragescen):
                
                # If the coverage levels have been passed in as a vector, convert it to an odict & sort by program names
                tmpbudget = dcp(thisprogset.getdefaultbudget())
                tmpcoverage = thisprogset.getprogcoverage(budget=tmpbudget, t=2000, parset=thisparset) # IT DOES NOT MATTER THE VALUE OF t YOU USE HERE!!!

                if isinstance(scen.coverage, list) or isinstance(scen.coverage, type(array([]))):
                    scen.coverage = vec2obj(scen.progset.getdefaultbuget(), newvec=scen.coverage) # It seems to be a vector: convert to odict -- WARNING, super dangerous!!
                if not isinstance(scen.coverage,dict): raise OptimaException('Currently only accepting coverage as dictionaries.')
                if not isinstance(scen.coverage,odict): scen.coverage = odict(scen.coverage)

                # Update, ensuring a consistent number of programs, using defaults where not provided -- WARNING, ugly
                tmpcoverage.update(scen.coverage)
                scen.coverage = tmpcoverage

                # Ensure coverage level values are lists
                for covkey, coventry in scen.coverage.iteritems():
                    if isnumber(coventry):
                        scen.coverage[covkey] = [coventry]

                # Figure out coverage
                scen.budget = thisprogset.getprogbudget(coverage=scen.coverage, t=scen.t, parset=thisparset, results=results)

            # TODO @Robyn have a look
            thisparsdict = thisprogset.getpars(coverage=scen.coverage, t=scen.t, parset=thisparset, results=results)
            scen.pars = thisparsdict
            for pardictno in range(len(thisparset.pars)): # Loop over all parameter dictionaries
                thisparset.pars[pardictno] = thisparsdict
            
        else: 
            errormsg = 'Unrecognized program scenario type.'
            raise OptimaException(errormsg)
            

        scenparsets[scen.name] = thisparset
        
    return scenparsets





def defaultscenarios(parset=None, verbose=2):
    """ Define a list of default scenarios -- only "Current conditions" by default """
    if parset is None: raise OptimaException('You need to supply a parset to generate default scenarios')
    
    scenlist = [Parscen()]
    
    ## Current conditions
    scenlist[0].name = 'Current conditions'
    scenlist[0].parset = parset
    scenlist[0].pars = [] # No changes
    
    return scenlist



def getparvalues(parset, scenpar):
    """
    Return the default parameter values from simpars for a given scenario parameter.
    
    defaultvals = getparvalues(parset, scenariolist[1]['pars'][2])
    
    Version: 2016jan30
    """
    npops = len(parset.pars[0]['popkeys'])
    simpars = parset.interp(start=scenpar['startyear'], end=scenpar['endyear'])

    original = simpars[scenpar['names'][0]]
    
    if scenpar['pops'] < npops: # It's for a specific population, get the value
        original = original[scenpar['pops'],:]
    else:
        original = original[:,:].mean(axis=0)
    initialindex = findinds(simpars['tvec'], scenpar['startyear'])
    finalindex = findinds(simpars['tvec'], scenpar['endyear'])

    startval = original[initialindex][0]
    endval = original[finalindex][0]
    return [startval, endval]
        



