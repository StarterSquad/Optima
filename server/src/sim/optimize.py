from printv import printv
from bunch import Bunch as struct
from copy import deepcopy

default_startyear = 2000
default_endyear = 2030

def optimize(D, objectives=None, constraints=None, timelimit=60, verbose=2):
    """
    Allocation optimization code:
        D is the project data structure
        objectives is a dictionary defining the objectives of the optimization
        constraints is a dictionary defining the constraints on the optimization
        timelimit is the maximum time in seconds to run optimization for
        verbose determines how much information to print.
        
    Version: 2014dec01 by cliffk
    """
    
    from model import model
    from copy import deepcopy
    from ballsd import ballsd
    from getcurrentbudget import getcurrentbudget
    from makemodelpars import makemodelpars
    from numpy import array
    printv('Running optimization...', 1, verbose)
    
    # Set options to update year range
    from setoptions import setoptions
    startyear = objectives.get("year").get("start") or default_startyear
    endyear = objectives.get("year").get("end") or default_endyear
    D.opt = setoptions(D.opt, startyear=startyear, endyear=endyear)
    # Make sure objectives and constraints exist
    if not isinstance(objectives, struct): objectives = defaultobjectives(D, verbose=verbose)
    if not isinstance(constraints, struct): constraints = defaultconstraints(D, verbose=verbose)

    objectives = deepcopy(objectives)
    constraints = deepcopy(constraints)

    # Convert weightings from percentage to number
    if objectives.outcome.inci: objectives.outcome.inciweight = float( objectives.outcome.inciweight ) / 100.0
    if objectives.outcome.daly: objectives.outcome.dalyweight = float( objectives.outcome.dalyweight ) / 100.0
    if objectives.outcome.death: objectives.outcome.deathweight = float( objectives.outcome.deathweight ) / 100.0
    if objectives.outcome.cost: objectives.outcome.costweight = float( objectives.outcome.costweight ) / 100.0

    for ob in objectives.money.objectives.keys():
        if objectives.money.objectives[ob].use: objectives.money.objectives[ob].by = float(objectives.money.objectives[ob].by) / 100.0

    for prog in objectives.money.costs.keys():
        objectives.money.costs[prog] = float(objectives.money.costs[prog]) / 100.0

    for prog in constraints.decrease.keys():
        if constraints.decrease[prog].use: constraints.decrease[prog].by = float(constraints.decrease[prog].by) / 100.0

    # Run optimization # TODO -- actually implement :)
    nallocs = 1 # WARNING, will want to do this better
    D.A = deepcopy([D.A[0]])
    for alloc in range(nallocs): D.A.append(deepcopy(D.A[0])) # Just copy for now
    D.A[0].label = 'Original'
    D.A[1].label = 'Optimal'
    origalloc = deepcopy(array(D.A[1].alloc))
    D.A = D.A[:2] # TODO WARNING KLUDGY
    

    
    def objectivecalc(alloc):
        """ Calculate the objective function """
        alloc /= sum(alloc)/sum(origalloc)
        newD = deepcopy(D)
        newD, newcov, newnonhivdalysaverted = getcurrentbudget(newD, alloc)
        newD.M = makemodelpars(newD.P, newD.opt, withwhat='c', verbose=0)
        S = model(newD.G, newD.M, newD.F[0], newD.opt, verbose=0)
        objective = S.death.sum() # TEMP
        
        return objective
                
        
    # Run the optimization algorithm
    optalloc, fval, exitflag, output = ballsd(objectivecalc, origalloc, xmin=0*array(origalloc), timelimit=timelimit, verbose=verbose)
    
    # Update the model
    for i,alloc in enumerate([origalloc,optalloc]):
        D, D.A[i].coverage, D.A[i].nonhivdalysaverted = getcurrentbudget(D, alloc)
        D.M = makemodelpars(D.P, D.opt, withwhat='c', verbose=2)
        D.A[i].S = model(D.G, D.M, D.F[0], D.opt, verbose=verbose)
        D.A[i].alloc = alloc # Now that it's run, store total program costs
    
    # Calculate results
    from makeresults import makeresults
    for alloc in range(len(D.A)):
        D.A[alloc].R = makeresults(D, [D.A[alloc].S], D.opt.quantiles, verbose=verbose)
    
    # Gather plot data
    from gatherplotdata import gatheroptimdata, gathermultidata
    D.plot.OA = gatheroptimdata(D, D.A, verbose=verbose)
    D.plot.OM = gathermultidata(D, D.A, verbose=verbose)
    
    printv('...done optimizing programs.', 2, verbose)
    return D

def saveoptimization(D, name, objectives, constraints, result = None, verbose=2):
    #save the optimization parameters
    new_optimization = struct()
    new_optimization.name = name
    new_optimization.objectives = objectives
    new_optimization.constraints = constraints
    if result: new_optimization.result = result

    if not "optimizations" in D:
        D.optimizations = [new_optimization]
    else:
        try:
            index = [item.name for item in D.optimizations].index(name)
            D.optimizations[index] = deepcopy(new_optimization)
        except:
            D.optimizations.append(new_optimization)
    return D

def removeoptimization(D, name):
    if "optimizations" in D:
        try:
            index = [item.name for item in D.optimizations].index(name)
            D.optimizations.pop(index)
        except:
            pass
    return D

def defaultobjectives(D, verbose=2):
    """
    Define default objectives for the optimization.
    """

    printv('Defining default objectives...', 3, verbose=verbose)

    ob = struct() # Dictionary of all objectives
    ob.year = struct() # Time periods for objectives
    ob.year.start = 2015 # "Year to begin optimization from"
    ob.year.numyears = 5 # "Number of years to optimize funding for"
    ob.year.end = 2030 # "Year to end optimization"
    ob.year.until = 2030 # "Year to project outcomes to"
    ob.what = 'outcome' # Alternative is "money"
    
    ob.outcome = struct()
    ob.outcome.fixed = 1e6 # "With a fixed amount of money available"
    ob.outcome.inci = True # "Minimize cumulative HIV incidence"
    ob.outcome.inciweight = 100 # "Incidence weighting"
    ob.outcome.daly = False # "Minimize cumulative DALYs"
    ob.outcome.dalyweight = 100 # "DALY weighting"
    ob.outcome.death = False # "Minimize cumulative AIDS-related deaths"
    ob.outcome.deathweight = 100 # "Death weighting"
    ob.outcome.cost = False # "Minimize cumulative DALYs"
    ob.outcome.costweight = 100 # "Cost weighting"
    ob.funding = "constant" #that's how it works on FE atm
    
    ob.money = struct()
    ob.money.objectives = struct()
    for objective in ['inci', 'incisex', 'inciinj', 'mtct', 'mtctbreast', 'mtctnonbreast', 'deaths', 'dalys']:
        ob.money.objectives[objective] = struct()
        ob.money.objectives[objective].use = False # TIck box: by default don't use
        ob.money.objectives[objective].by = 50 # "By" text entry box: 0.5 = 50% reduction
        ob.money.objectives[objective].to = 0 # "To" text entry box: don't use if set to 0
    ob.money.objectives.inci.use = True # Set incidence to be on by default
    
    ob.money.costs = struct()
    for prog in D.programs.keys():
        ob.money.costs[prog] = 100 # By default, use a weighting of 100%
    
    return ob

def defaultconstraints(D, verbose=2):
    """
    Define default constraints for the optimization.
    """
    
    printv('Defining default constraints...', 3, verbose=verbose)
    
    con = struct()
    con.txelig = 4 # 4 = "All people diagnosed with HIV"
    con.dontstopart = True # "No one who initiates treatment is to stop receiving ART"
    con.decrease = struct()
    for prog in D.programs.keys(): # Loop over all defined programs
        con.decrease[prog] = struct()
        con.decrease[prog].use = False # Tick box: by default don't use
        con.decrease[prog].by = 50 # Text entry box: 0.5 = 50% per year
    
    con.coverage = struct()
    for prog in D.programs.keys(): # Loop over all defined programs
        con.coverage[prog] = struct()
        con.coverage[prog].use = False # Tick box: by default don't use
        con.coverage[prog].level = 0 # First text entry box: default no limit
        con.coverage[prog].year = 2030 # Year to reach coverage level by
        
    return con


def defaultoptimizations(D, verbose=2):
    """ Define a list of default optimizations (one default optimization) """
    
    # Start at the very beginning, a very good place to start :)
    optimizations = [struct()]
    
    ## Current conditions
    optimizations[0].name = 'Default'
    optimizations[0].constraints = defaultconstraints(D, verbose)
    optimizations[0].objectives = defaultobjectives(D, verbose)
    return optimizations


