"""
Allocation optimization code:
    D is the project data structure
    objectives is a dictionary defining the objectives of the optimization
    constraints is a dictionary defining the constraints on the optimization
    timelimit is the maximum time in seconds to run optimization for
    verbose determines how much information to print.
    
Version: 2015feb03 by cliffk
"""

from printv import printv
from bunch import Bunch as struct
from copy import deepcopy
from numpy import ones, zeros, concatenate, arange, inf, hstack, argmin, array, ndim
from utils import findinds
from makeresults import makeresults
from timevarying import timevarying, multiyear
from getcurrentbudget import getcurrentbudget
from model import model
from makemodelpars import makemodelpars
from quantile import quantile
from ballsd import ballsd


def runmodelalloc(D, thisalloc, parindices, randseed, financial=True, verbose=2):
    """ Little function to do calculation since it appears so many times """
    newD = deepcopy(D)
    newD, newcov, newnonhivdalysaverted = getcurrentbudget(newD, thisalloc, randseed=randseed) # Get cost-outcome curves with uncertainty
    newM = makemodelpars(newD.P, newD.opt, withwhat='c', verbose=verbose)
    newD.M = partialupdateM(D.M, newM, parindices)
    S = model(newD.G, newD.M, newD.F[0], newD.opt, verbose=verbose)
    R = makeresults(D, allsims=[S], financial=financial, verbose=0)
    return R



def objectivecalc(optimparams, options):
    """ Calculate the objective function """
    if 'ntimepm' in options.keys():
        thisalloc = timevarying(optimparams, ntimepm=options.ntimepm, nprogs=options.nprogs, tvec=options.D.opt.partvec, totalspend=options.totalspend, fundingchanges=options.fundingchanges) 
    elif 'years' in options.keys():
        thisalloc = multiyear(optimparams, years=options.years, totalspends=options.totalspends, nprogs=options.nprogs, tvec=options.D.opt.partvec) 
    else:
        raise Exception('Cannot figure out what kind of allocation this is since neither options.ntimepm nor options.years is defined')
    
    financial=True if options.weights['costann'] else False
    R = runmodelalloc(options.D, thisalloc, options.parindices, options.randseed, financial=financial) # Actually run
    
    outcome = 0 # Preallocate objective value 
    for key in options.outcomekeys:
        if options.weights[key]>0: # Don't bother unless it's actually used
            if key!='costann': thisoutcome = R[key].tot[0][options.outindices].sum()
            else: thisoutcome = R[key].total.total[0][options.outindices].sum() # Special case for costann
            outcome += thisoutcome * options.weights[key] / float(options.normalizations[key]) * options.D.opt.dt # Calculate objective
        
    return outcome
    
    
    
def optimize(D, objectives=None, constraints=None, maxiters=1000, timelimit=None, verbose=5, name='Default', stoppingfunc = None):
    """ Perform the actual optimization """
    from time import sleep
    
    printv('Running optimization...', 1, verbose)
    
    # Set up parameter vector for time-varying optimisation...
    stepsize = 100000
    growsize = 0.01
    verbose = 5
    print('TEMP')

    origR = deepcopy(D.R)
    origalloc = D.data.origalloc
    
    # Make sure objectives and constraints exist, and overwrite using saved ones if available
    if objectives is None: objectives = defaultobjectives(D, verbose=verbose)
    if constraints is None: constraints = defaultconstraints(D, verbose=verbose)

    # Do this so if e.g. /100 won't have problems
    objectives = deepcopy(objectives)
    constraints = deepcopy(constraints)
    ntimepm=1 + int(objectives.timevarying)*int(objectives.funding=='constant') # Either 1 or 2, but only if funding==constant

    nprogs = len(D.programs)
    totalspend = objectives.outcome.fixed # For fixed budgets
    
    # Define constraints on funding -- per year and total
    fundingchanges = struct()
    keys1 = ['year','total']
    keys2 = ['dec','inc']
    abslims = {'dec':0, 'inc':1e9}
    rellims = {'dec':-1e9, 'inc':1e9}
    for key1 in keys1:
        fundingchanges[key1] = struct()
        for key2 in keys2:
            fundingchanges[key1][key2] = []
            for p in range(nprogs):
                fullkey = key1+key2+'rease'
                this = constraints[fullkey][p] # Shorten name
                if key1=='total':
                    if this.use and objectives.funding != 'variable': # Don't constrain variable-year-spend optimizations
                        newlim = this.by/100.*origalloc
                        fundingchanges[key1][key2].append(newlim)
                    else: 
                        fundingchanges[key1][key2].append(abslims[key2])
                elif key1=='year':
                    if this.use and objectives.funding != 'variable': # Don't constrain variable-year-spend optimizations
                        newlim = this.by/100.-1 # Relative change in funding
                        fundingchanges[key1][key2].append(newlim)
                    else: 
                        fundingchanges[key1][key2].append(rellims[key2])
                
    
    ## Define indices, weights, and normalization factors
    initialindex = findinds(D.opt.partvec, objectives.year.start)
    finalparindex = findinds(D.opt.partvec, objectives.year.end)
    finaloutindex = findinds(D.opt.partvec, objectives.year.until)
    parindices = arange(initialindex,finalparindex)
    outindices = arange(initialindex,finaloutindex)
    weights = dict()
    normalizations = dict()
    outcomekeys = ['inci', 'death', 'daly', 'costann']
    if sum([objectives.outcome[key] for key in outcomekeys])>1: # Only normalize if multiple objectives, since otherwise doesn't make a lot of sense
        for key in outcomekeys:
            thisweight = objectives.outcome[key+'weight'] * objectives.outcome[key] / 100.
            weights.update({key:thisweight}) # Get weight, and multiply by "True" or "False" and normalize from percentage
            if key!='costann': thisnormalization = origR[key].tot[0][outindices].sum()
            else: thisnormalization = origR[key].total.total[0][outindices].sum() # Special case for costann
            normalizations.update({key:thisnormalization})
    else:
        for key in outcomekeys:
            weights.update({key:int(objectives.outcome[key])}) # Weight of 1
            normalizations.update({key:1}) # Normalizatoin of 1
        
    # Initiate probabilities of parameters being selected
    stepsizes = zeros(nprogs * ntimepm)
    
    # Easy access initial allocation indices and turn stepsizes into array
    ai = range(nprogs)
    gi = range(nprogs,   nprogs*2) if ntimepm >= 2 else []
    si = range(nprogs*2, nprogs*3) if ntimepm >= 3 else []
    ii = range(nprogs*3, nprogs*4) if ntimepm >= 4 else []
    
    # Turn stepsizes into array
    stepsizes[ai] = stepsize
    stepsizes[gi] = growsize if ntimepm > 1 else 0
    stepsizes[si] = stepsize
    stepsizes[ii] = growsize # Not sure that growsize is an appropriate starting point
    
    # Initial values of time-varying parameters
    growthrate = zeros(nprogs)   if ntimepm >= 2 else []
    saturation = origalloc       if ntimepm >= 3 else []
    inflection = ones(nprogs)*.5 if ntimepm >= 4 else []
    
    # Concatenate parameters to be optimised
    optimparams = concatenate((origalloc, growthrate, saturation, inflection)) # WARNING, not used for multi-year optimizations
        
    
    
    
    ###########################################################################
    ## Constant budget optimization
    ###########################################################################
    if objectives.funding == 'constant' and objectives.timevarying == False:
        
        ## Define options structure
        options = struct()
        options.ntimepm = ntimepm # Number of time-varying parameters
        options.nprogs = nprogs # Number of programs
        options.D = deepcopy(D) # Main data structure
        options.outcomekeys = outcomekeys # Names of outcomes, e.g. 'inci'
        options.weights = weights # Weights for each parameter
        options.outindices = outindices # Indices for the outcome to be evaluated over
        options.parindices = parindices # Indices for the parameters to be updated on
        options.normalizations = normalizations # Whether to normalize a parameter
        options.totalspend = totalspend # Total budget
        options.fundingchanges = fundingchanges # Constraints-based funding changes
        
        
        ## Run with uncertainties
        allocarr = []
        fvalarr = []
        for s in range(len(D.F)): # Loop over all available meta parameters
            print('========== Running uncertainty optimization %s of %s... ==========' % (s+1, len(D.F)))
            options.D.F = [D.F[s]] # Loop over fitted parameters
            options.randseed = s
            optparams, fval, exitflag, output = ballsd(objectivecalc, optimparams, options=options, xmin=fundingchanges.total.dec, xmax=fundingchanges.total.inc, absinitial=stepsizes, MaxIter=maxiters, timelimit=timelimit, fulloutput=True, stoppingfunc=stoppingfunc, verbose=verbose)
            optparams = optparams / optparams.sum() * options.totalspend # Make sure it's normalized -- WARNING KLUDGY
            allocarr.append(optparams)
            fvalarr.append(output.fval)
        
        ## Find which optimization was best
        bestallocind = -1
        bestallocval = inf
        for s in range(len(fvalarr)):
            if fvalarr[s][-1]<bestallocval:
                bestallocval = fvalarr[s][-1]
                bestallocind = s
        if bestallocind == -1: print('WARNING, best allocation value seems to be infinity!')
        
        # Update the model and store the results
        result = struct()
        result.kind = 'constant'
        result.fval = fvalarr[bestallocind] # Append the best value noe
        result.allocarr = [] # List of allocations
        result.allocarr.append(quantile([origalloc])) # Kludgy -- run fake quantile on duplicated origalloc just so it matches
        result.allocarr.append(quantile(allocarr)) # Calculate allocation arrays 
        labels = ['Original','Optimal']
        result.Rarr = []
        for params in [origalloc, allocarr[bestallocind]]: # CK: loop over original and (the best) optimal allocations
            sleep(0.1)
            alloc = timevarying(params, ntimepm=len(params)/nprogs, nprogs=nprogs, tvec=D.opt.partvec, totalspend=totalspend, fundingchanges=fundingchanges)   
            R = runmodelalloc(options.D, alloc, options.parindices, options.randseed, verbose=verbose) # Actually run
            result.Rarr.append(struct()) # Append a structure
            result.Rarr[-1].R = deepcopy(R) # Store the R structure (results)
            result.Rarr[-1].label = labels.pop(0) # Store labels, one at a time
        
        
        
    
    
    ###########################################################################
    ## Time-varying budget optimization
    ###########################################################################
    if objectives.funding == 'constant' and objectives.timevarying == True:
        
        ## Define options structure
        options = struct()
        options.ntimepm = ntimepm # Number of time-varying parameters
        options.nprogs = nprogs # Number of programs
        options.D = deepcopy(D) # Main data structure
        options.outcomekeys = outcomekeys # Names of outcomes, e.g. 'inci'
        options.weights = weights # Weights for each parameter
        options.outindices = outindices # Indices for the outcome to be evaluated over
        options.parindices = parindices # Indices for the parameters to be updated on
        options.normalizations = normalizations # Whether to normalize a parameter
        options.totalspend = totalspend # Total budget
        options.fundingchanges = fundingchanges # Constraints-based funding changes
        parammin = concatenate((fundingchanges.total.dec, ones(nprogs)*-1e9))  
        parammax = concatenate((fundingchanges.total.inc, ones(nprogs)*1e9))  
        options.randseed = None
        
        
        
        ## Run time-varying optimization
        print('========== Running time-varying optimization ==========')
        optparams, fval, exitflag, output = ballsd(objectivecalc, optimparams, options=options, xmin=parammin, xmax=parammax, absinitial=stepsizes, MaxIter=maxiters, timelimit=timelimit, fulloutput=True, stoppingfunc=stoppingfunc, verbose=verbose)
        optparams = optparams / optparams.sum() * options.totalspend # Make sure it's normalized -- WARNING KLUDGY
        
        # Update the model and store the results
        result = struct()
        result.kind = 'timevarying'
        result.fval = output.fval # Append the objective sequence
        result.Rarr = []
        labels = ['Original','Optimal']
        for params in [origalloc, optparams]: # CK: loop over original and (the best) optimal allocations
            sleep(0.1)
            alloc = timevarying(params, ntimepm=len(params)/nprogs, nprogs=nprogs, tvec=D.opt.partvec, totalspend=totalspend, fundingchanges=fundingchanges) #Regenerate allocation
            R = runmodelalloc(options.D, alloc, options.parindices, options.randseed, verbose=verbose) # Actually run
            result.Rarr.append(struct()) # Append a structure
            result.Rarr[-1].R = deepcopy(R) # Store the R structure (results)
            result.Rarr[-1].label = labels.pop(0) # Store labels, one at a time
        result.xdata = R.tvec # Store time data
        result.alloc = alloc[:,0:len(R.tvec)] # Store allocation data, and cut to be same length as time data
        
    
    
    
        
    ###########################################################################
    ## Multiple-year budget optimization
    ###########################################################################
    if objectives.funding == 'variable':
        
        ## Define options structure
        options = struct()
        
        options.nprogs = nprogs # Number of programs
        options.D = deepcopy(D) # Main data structure
        options.outcomekeys = outcomekeys # Names of outcomes, e.g. 'inci'
        options.weights = weights # Weights for each parameter
        options.outindices = outindices # Indices for the outcome to be evaluated over
        options.parindices = parindices # Indices for the parameters to be updated on
        options.normalizations = normalizations # Whether to normalize a parameter
        
        options.randseed = None # Death is enough randomness on its own
        options.fundingchanges = fundingchanges # Constraints-based funding changes
        
        options.years = []
        options.totalspends = []
        yearkeys = objectives.outcome.variable.keys()
        yearkeys.sort() # God damn I hate in-place methods
        for key in yearkeys: # Stored as a list of years:
            options.years.append(float(key)) # Convert from string to number
            options.totalspends.append(objectives.outcome.variable[key]) # Append this year
        
        
        
        ## Define optimization parameters
        nyears = len(options.years)
        optimparams = array(origalloc.tolist()*nyears).flatten() # Duplicate parameters
        parammin = zeros(len(optimparams))
        stepsizes = stepsize + zeros(len(optimparams))
        keys1 = ['year','total']
        keys2 = ['dec','inc']
        abslims = {'dec':0, 'inc':1e9}
        rellims = {'dec':-1e9, 'inc':1e9}
        for key1 in keys1:
            for key2 in keys2:
                options.fundingchanges[key1][key2] *= nyears # I know this just points to the list rather than copies, but should be fine. I hope
        
        ## Run time-varying optimization
        print('========== Running multiple-year optimization ==========')
        optparams, fval, exitflag, output = ballsd(objectivecalc, optimparams, options=options, xmin=fundingchanges.total.dec, xmax=fundingchanges.total.inc, MaxIter=maxiters, timelimit=timelimit, fulloutput=True, stoppingfunc=stoppingfunc, verbose=verbose)
        
        # Normalize
        proginds = arange(nprogs)
        optparams = array(optparams)
        for y in range(nyears):
            theseinds = proginds+y*nprogs
            optparams[theseinds] *= options.totalspends[y] / float(sum(optparams[theseinds]))
        optparams = optparams.tolist()
        
        # Update the model and store the results
        result = struct()
        result.kind = 'multiyear'
        result.fval = output.fval # Append the objective sequence
        result.Rarr = []
        labels = ['Original','Optimal']
        for params in [origalloc, optparams]: # CK: loop over original and (the best) optimal allocations
            sleep(0.1)
            alloc = multiyear(optimparams, years=options.years, totalspends=options.totalspends, nprogs=options.nprogs, tvec=options.D.opt.partvec) 
            R = runmodelalloc(options.D, alloc, options.parindices, options.randseed, verbose=verbose) # Actually run
            result.Rarr.append(struct()) # Append a structure
            result.Rarr[-1].R = deepcopy(R) # Store the R structure (results)
            result.Rarr[-1].label = labels.pop(0) # Store labels, one at a time
        result.xdata = R.tvec # Store time data
        result.alloc = alloc[:,0:len(R.tvec)] # Store allocation data, and cut to be same length as time data
        
    
    
    
    
        
        
        
    ###########################################################################
    ## Multiple budgets optimization
    ###########################################################################
    if objectives.funding == 'range':
        
        ## Define options structure
        options = struct()
        options.ntimepm = 1 # Number of time-varying parameters -- always 1 in this case
        options.nprogs = nprogs # Number of programs
        options.D = deepcopy(D) # Main data structure
        options.outcomekeys = outcomekeys # Names of outcomes, e.g. 'inci'
        options.weights = weights # Weights for each parameter
        options.outindices = outindices # Indices for the outcome to be evaluated over
        options.parindices = parindices # Indices for the parameters to be updated on
        options.normalizations = normalizations # Whether to normalize a parameter
        options.fundingchanges = fundingchanges # Constraints-based funding changes
        options.totalspend = totalspend # Total budget
        options.randseed = None
        
        ## Run multiple budgets
        budgets = arange(objectives.outcome.budgetrange.minval, objectives.outcome.budgetrange.maxval+objectives.outcome.budgetrange.step, objectives.outcome.budgetrange.step)
        closesttocurrent = argmin(abs(budgets-1)) + 1 # Find the index of the budget closest to current and add 1 since prepend current budget
        nbudgets = len(budgets)
        budgets = hstack([1,budgets]) # Include current budget
        allocarr = [origalloc] # Original allocation
        fvalarr = [objectivecalc(optimparams, options=options)] # Outcome for original allocation
        for b in range(nbudgets):
            print('========== Running budget optimization %s of %s... ==========' % (b+1, nbudgets))
            options.totalspend = totalspend*budgets[b+1] # Total budget, skipping first
            optparams, fval, exitflag, output = ballsd(objectivecalc, optimparams, options=options, xmin=fundingchanges.total.dec, xmax=fundingchanges.total.inc, absinitial=stepsizes, MaxIter=maxiters, timelimit=timelimit, fulloutput=True, stoppingfunc=stoppingfunc, verbose=verbose)
            optparams = optparams / optparams.sum() * options.totalspend # Make sure it's normalized -- WARNING KLUDGY
            allocarr.append(optparams)
            fvalarr.append(fval) # Only need last value
        
        # Update the model and store the results
        result = struct()
        result.kind = objectives.funding
        result.budgets = budgets
        result.budgetlabels = ['Original budget']
        for b in range(nbudgets): result.budgetlabels.append('%i%% budget' % (budgets[b+1]*100./float(budgets[0])))
            
        result.fval = fvalarr # Append the best value
        result.allocarr = allocarr # List of allocations
        labels = ['Original','Optimal']
        result.Rarr = []
        for params in [origalloc, allocarr[closesttocurrent]]: # CK: loop over original and (the best) optimal allocations
            sleep(0.1)
            alloc = timevarying(params, ntimepm=len(params)/nprogs, nprogs=nprogs, tvec=D.opt.partvec, totalspend=totalspend, fundingchanges=fundingchanges)   
            R = runmodelalloc(options.D, alloc, options.parindices, options.randseed, verbose=verbose) # Actually run
            result.Rarr.append(struct()) # Append a structure
            result.Rarr[-1].R = deepcopy(R) # Store the R structure (results)
            result.Rarr[-1].label = labels.pop(0) # Store labels, one at a time        
        
        
    
    ## Gather plot data
    from gatherplotdata import gatheroptimdata
    optim = gatheroptimdata(D, result, verbose=verbose)
    if 'optim' not in D.plot: D.plot.optim = [] # Initialize list if required
#    D.plot.optim.append(optim) # In any case, append
    D.plot.optim=[optim]
    
    ## Save optimization to D
    saveoptimization(D, name, objectives, constraints, result, verbose=2)
    
    printv('...done optimizing programs.', 2, verbose)
    return D











def saveoptimization(D, name, objectives, constraints, result, verbose=2):
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
    ob.year.end = 2020 # "Year to end optimization"
    ob.year.until = 2030 # "Year to project outcomes to"
    ob.what = 'outcome' # Alternative is "money"
    
    ob.outcome = struct()
    ob.outcome.fixed = sum(D.data.origalloc) # "With a fixed amount of money available"
    ob.outcome.inci = True # "Minimize cumulative HIV incidence"
    ob.outcome.inciweight = 100 # "Incidence weighting"
    ob.outcome.daly = False # "Minimize cumulative DALYs"
    ob.outcome.dalyweight = 100 # "DALY weighting"
    ob.outcome.death = False # "Minimize cumulative AIDS-related deaths"
    ob.outcome.deathweight = 100 # "Death weighting"
    ob.outcome.costann = False # "Minimize cumulative DALYs"
    ob.outcome.costannweight = 100 # "Cost weighting"
    ob.outcome.variable = [] # No variable budgets by default
    ob.outcome.budgetrange = struct() # For running multiple budgets
    ob.outcome.budgetrange.minval = None
    ob.outcome.budgetrange.maxval = None
    ob.outcome.budgetrange.step = None
    ob.funding = "constant" #that's how it works on FE atm
    
    # Other settings
    ob.timevarying = False # Do not use time-varying parameters
    ob.artcontinue = 1 # No one currently on ART stops
    ob.otherprograms = "remain" # Other programs remain constant after optimization ends
    
    ob.money = struct()
    ob.money.objectives = struct()
    for objective in ['inci', 'incisex', 'inciinj', 'mtct', 'mtctbreast', 'mtctnonbreast', 'deaths', 'dalys']:
        ob.money.objectives[objective] = struct()
        ob.money.objectives[objective].use = False # TIck box: by default don't use
        ob.money.objectives[objective].by = 50 # "By" text entry box: 0.5 = 50% reduction
        ob.money.objectives[objective].to = 0 # "To" text entry box: don't use if set to 0
    ob.money.objectives.inci.use = True # Set incidence to be on by default
    
    ob.money.costs = []
    for p in range(D.G.nprogs):
        ob.money.costs.append(100) # By default, use a weighting of 100%
    
    return ob

def defaultconstraints(D, verbose=2):
    """
    Define default constraints for the optimization.
    """
    
    printv('Defining default constraints...', 3, verbose=verbose)
    
    con = struct()
    con.txelig = 4 # 4 = "All people diagnosed with HIV"
    con.dontstopart = True # "No one who initiates treatment is to stop receiving ART"
    con.yeardecrease = []
    con.yearincrease = []
    for p in range(D.G.nprogs): # Loop over all defined programs
        con.yeardecrease.append(struct())
        con.yeardecrease[p].use = False # Tick box: by default don't use
        con.yeardecrease[p].by = 80 # Text entry box: 0.5 = 50% per year
        con.yearincrease.append(struct())
        con.yearincrease[p].use = False # Tick box: by default don't use
        con.yearincrease[p].by = 120 # Text entry box: 0.5 = 50% per year
    con.totaldecrease = []
    con.totalincrease = []
    for p in range(D.G.nprogs): # Loop over all defined programs
        con.totaldecrease.append(struct())
        con.totaldecrease[p].use = False # Tick box: by default don't use
        con.totaldecrease[p].by = 50 # Text entry box: 0.5 = 50% per total
        con.totalincrease.append(struct())
        con.totalincrease[p].use = False # Tick box: by default don't use
        con.totalincrease[p].by = 200 # Text entry box: 0.5 = 50% total
    
    con.coverage = []
    for p in range(D.G.nprogs): # Loop over all defined programs
        con.coverage.append(struct())
        con.coverage[p].use = False # Tick box: by default don't use
        con.coverage[p].level = 0 # First text entry box: default no limit
        con.coverage[p].year = 2030 # Year to reach coverage level by
        
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




def partialupdateM(oldM, newM, indices, setbefore=False, setafter=True):
    """ 
    Update M, but only for certain indices. If setbefore is true, reset all values before indices to new value; similarly for setafter. 
    WARNING: super ugly code!!
    """
    from makemodelpars import totalacts
    output = deepcopy(oldM)
    for key in output.keys():
        if key not in ['transit', 'pships', 'const', 'tvec', 'hivprev', 'totalacts']: # Exclude certain keys that won't be updated
            if hasattr(output[key],'keys'): # It's a dict or a bunch, loop again
                for key2 in output[key].keys():
                    try:
                        if ndim(output[key][key2])==1:
                            output[key][key2][indices] = newM[key][key2][indices]
                            if setbefore: output[key][key2][:min(indices)] = newM[key][key2][:min(indices)]
                            if setafter: output[key][key2][max(indices):] = newM[key][key2][max(indices):]
                        elif ndim(output[key][key2])==2:
                            output[key][key2][:,indices] = newM[key][key2][:,indices]
                            if setbefore: output[key][key2][:,:min(indices)] = newM[key][key2][:,:min(indices)]
                            if setafter: output[key][key2][:,max(indices):] = newM[key][key2][:,max(indices):]
                        else:
                            raise Exception('%i dimensions for parameter M.%s.%s' % (ndim(output[key][key2][indices]), key, key2))
                    except:
                        print('Could not set indices for parameter M.%s.%s, indices %i-%i' % (key, key2, min(indices), max(indices)))
                        import traceback; traceback.print_exc(); import pdb; pdb.set_trace()
            else:
                try:
                    if ndim(output[key])==1:
                        output[key][indices] = newM[key][indices]
                        if setbefore: output[key][:min(indices)] = newM[key][:min(indices)]
                        if setafter: output[key][max(indices):] = newM[key][max(indices):]
                    elif ndim(output[key])==2:
                        output[key][:,indices] = newM[key][:,indices]
                        if setbefore: output[key][:,:min(indices)] = newM[key][:,:min(indices)]
                        if setafter: output[key][:,max(indices):] = newM[key][:,max(indices):]
                    else:
                        raise Exception('%i dimensions for parameter M.%s' % (ndim(output[key][indices]), key, key2))
                except:
                    print('Could not set indices for parameter M.%s, indices %i-%i' % (key, min(indices), max(indices)))
                    import traceback; traceback.print_exc(); import pdb; pdb.set_trace()
    
    output.totalacts = totalacts(output, len(output.tvec)) # Update total acts
    return output
