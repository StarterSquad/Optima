## Imports
from math import pow as mpow
from numpy import zeros, exp, maximum, minimum, hstack, inf, array, isnan
from optima import OptimaException, printv, dcp, odict, findinds, makesimpars, Resultset

def model(simpars=None, settings=None, verbose=None, die=False, debug=False):
    """
    Runs Optima's epidemiological model.
    
    Version: 1.4 (2016mar04)
    """
    
    ##################################################################################################################
    ### Setup
    ##################################################################################################################

    # Hard-coded parameters that hopefully don't matter too much
    cd4transnorm = 1.2 # Was 3.3 -- estimated overestimate of infectiousness by splitting transmissibility multiple ways -- see commit 57057b2486accd494ef9ce1379c87a6abfababbd for calculations
    
    # Initialize basic quantities
    if simpars is None: raise OptimaException('model() requires simpars as an input')
    if settings is None: raise OptimaException('model() requires settings as an input')
    popkeys      = simpars['popkeys']
    npops        = len(popkeys)
    simpars      = dcp(simpars)
    tvec         = simpars['tvec']
    dt           = simpars['dt']      # Shorten dt
    npts         = len(tvec) # Number of time points
    ncd4         = settings.ncd4      # Shorten number of CD4 states
    nstates      = settings.nstates   # Shorten number of health states
    people       = zeros((nstates, npops, npts)) # Matrix to hold everything
    allpeople    = zeros((npops, npts)) # Population sizes
    effhivprev   = zeros((npops, 1))    # HIV effective prevalence (prevalence times infectiousness)
    inhomo       = zeros(npops)    # Inhomogeneity calculations
    usecascade   = settings.usecascade # Whether or not the full treatment cascade should be used
    safetymargin = settings.safetymargin # Maximum fraction of people to move on a single timestep
    eps          = settings.eps # Define another small number to avoid divide-by-zero errors
    if verbose is None: verbose = settings.verbose # Verbosity of output
    
    # Would be at the top of the script, but need to figure out verbose first
    printv('Running model...', 1, verbose)
    
    # Initialize arrays
    raw_inci       = zeros((npops, npts)) # Total incidence
    raw_mtct       = zeros((npops, npts)) # Number of mother-to-child transmissions to each population
    raw_diag       = zeros((npops, npts)) # Number diagnosed per timestep
    raw_newtreat   = zeros((npops, npts)) # Number initiating ART1 per timestep
    raw_death      = zeros((npops, npts)) # Number of deaths per timestep
    raw_otherdeath = zeros((npops, npts)) # Number of other deaths per timestep
    
    # Biological and failure parameters -- death etc
    prog       = array([simpars['progacute'], simpars['proggt500'], simpars['proggt350'], simpars['proggt200'], simpars['proggt50']]) # Ugly, but fast
    recov      = array([simpars['recovgt500'], simpars['recovgt350'], simpars['recovgt200'], simpars['recovgt50']])
    death      = array([simpars['deathacute'], simpars['deathgt500'], simpars['deathgt350'], simpars['deathgt200'], simpars['deathgt50'], simpars['deathlt50']])
    deathtx    = simpars['deathtreat']   # Death rate whilst on treatment
    cd4trans   = array([simpars['cd4transacute'], simpars['cd4transgt500'], simpars['cd4transgt350'], simpars['cd4transgt200'], simpars['cd4transgt50'], simpars['cd4translt50']])


    # Defined for total (not by populations) and time dependent [npts]
    treatvs     = simpars['treatvs']     # viral suppression - ART initiators (P)
    if usecascade:
        biofailure    = simpars['biofailure']  # biological treatment failure rate (P/T)
        freqvlmon     = simpars['freqvlmon']     # Viral load monitoring frequency (N/T)
        restarttreat  = simpars['restarttreat']  # Rate of ART re-inititation (P/T)
        progusvl      = simpars['progusvl']      # Proportion of people who progress when on unsuppressive ART
        recovusvl     = simpars['recovusvl']     # Proportion of people who recover when on unsuppressive ART
        stoppropcare  = simpars['stoppropcare']  # Proportion of people lost-to-follow-up who are actually still in care (transferred)
        # Behavioural transitions between stages [npop,npts]
        immediatecare = simpars['immediatecare'] # Linkage to care from diagnosis within 1 month (%) (P)
        linktocare    = simpars['linktocare']    # rate of linkage to care (P/T)
        stoprate      = simpars['stoprate']      # Percentage of people who receive ART in year who stop taking ART (%/year) (P/T)
        leavecare     = simpars['leavecare']     # Proportion of people in care then lost to follow-up per year (P/T)
        
        
    
    # Calculate other things outside the loop
    transinj = simpars['transinj']          # Injecting
    cd4trans /= cd4transnorm # Normalize CD4 transmission
    dxfactor = (1.0-simpars['effdx']) # Include diagnosis efficacy
    if usecascade:
        efftxunsupp = (1-simpars['efftxunsupp']) * dxfactor # (~30%) reduction in transmission probability for usVL
        efftxsupp  = (1-simpars['efftxsupp'])  * dxfactor # (~96%) reduction in transmission probability for sVL
    else:
        txfactor = dxfactor * ((1-simpars['efftxsupp'])*treatvs + (1-simpars['efftxunsupp'])*(1-treatvs)) # Roughly calculate treatment efficacy based on ART success rate; should be 92%*90% = 80%, close to 70% we had been using

    # Disease state indices
    susreg   = settings.susreg      # Susceptible, regular
    progcirc = settings.progcirc     # Susceptible, programmatically circumcised
    sus      = settings.sus      # Susceptible, both circumcised and uncircumcised
    undx     = settings.undx     # Undiagnosed
    dx       = settings.dx       # Diagnosed
    alldx    = settings.alldx    # All diagnosed
    alltx    = settings.alltx    # All on treatment
    allplhiv = settings.allplhiv # All PLHIV
    aidsind  = settings.aidsind  # Index for when people have AIDS (used for assinging AIDS testing rate)
    if usecascade:
        care    = settings.care    # in care
        usvl    = settings.usvl    # On treatment - Unsuppressed Viral Load
        svl     = settings.svl     # On treatment - Suppressed Viral Load
        lost    = settings.lost    # Not on ART (anymore) and lost to follow-up
        off     = settings.off     # off ART but still in care
        allcare = settings.allcare # All people in care
    else:
        tx   = settings.tx  # Treatment -- equal to settings.svl, but this is clearer
    if debug and len(sus)!=2:
        errormsg = 'Definition of susceptibles has changed: expecting regular circumcised + VMMC, but actually length %i' % len(sus)
        raise OptimaException(errormsg)
    
    # Proportion aware and treated (for 90/90/90)
    propdx = simpars['propdx']
    if usecascade: 
        propcare = simpars['propcare']
        propsupp = simpars['propsupp']
    proptx = simpars['proptx']

    # Population sizes
    popsize = dcp(simpars['popsize'])
    
    # Population characteristics
    male    = simpars['male']          # Boolean array, true for males
    female  = simpars['female']      # Boolean array, true for females
    injects = simpars['injects']    # Boolean array, true for PWID

    # Intervention uptake (P=proportion, N=number)
    sharing   = simpars['sharing']   # Sharing injecting equiptment (P)
    numtx     = simpars['numtx']     # 1st line treatement (N) -- tx already used for index of people on treatment [npts]
    hivtest   = simpars['hivtest']   # HIV testing (P) [npop,npts]
    aidstest  = simpars['aidstest']  # HIV testing in AIDS stage (P) [npts]
    numcirc   = simpars['numcirc']   # Number of programmatic circumcisions performed (N)
    numpmtct  = simpars['numpmtct']  # Number of people receiving PMTCT (N)
    
    # Uptake of OST
    numost = simpars['numost']                  # Number of people on OST (N)
    if any(injects):
        numpwid = popsize[injects,:].sum(axis=0)  # Total number of PWID
        try: 
            ostprev = numost/numpwid # Proportion of PWID on OST (P)
            ostprev = minimum(ostprev, 1.0) # Don't let more than 100% of PWID be on OST :)
        except: 
            errormsg = 'Cannot divide by the number of PWID (numost=%f, numpwid=5f' % (numost, numpwid)
            if die: raise OptimaException(errormsg)
            else: 
                printv(errormsg, 1, verbose)
                ostprev = zeros(npts) # Reset to zero
    else: # No one injects
        if sum(numost): 
            errormsg = 'You have entered non-zero value for the number of PWID on OST, but you have not specified any populations who inject'
            if die: raise OptimaException(errormsg)
            else: 
                printv(errormsg, 1, verbose)
                ostprev = zeros(npts)
        else: # No one on OST
            ostprev = zeros(npts)
    
    # Other interventions
    effcondom = simpars['effcondom']         # Condom effect
    circconst = 1 - simpars['effcirc'] # Actual efficacy 
    circeff   = 1 - simpars['propcirc']*simpars['effcirc']  # Actual efficacy 
    prepeff   = 1 - simpars['effprep']*simpars['prep']  # PrEP effect
    osteff    = 1 - simpars['effost']*ostprev  # OST effect
    stieff    = 1 + simpars['effsti']*simpars['stiprev'] # STI effect
    effmtct   = simpars['mtctbreast']*simpars['breast'] + simpars['mtctnobreast']*(1-simpars['breast']) # Effective MTCT transmission
    pmtcteff  = (1 - simpars['effpmtct']) * effmtct # Effective MTCT transmission whilst on PMTCT

    # Calculate these things outside of the time loop
    prepsti = prepeff*stieff
    prepsticirceff = prepsti*circeff
    prepsticircconst = prepsti*circconst
    


    # Force of infection metaparameter
    force = simpars['force']
    inhomopar = simpars['inhomo'] # WARNING, name is not consistent -- should be "inhomo"

    # More parameters...should maybe be moved somewhere else?
    birth = simpars['birth']
    agetransit = simpars['agetransit']*dt # Multiply by dt here so don't have to later
    risktransit = simpars['risktransit']*dt # Multiply by dt here so don't have to later
    birthtransit = simpars['birthtransit']*dt # Multiply by dt here so don't have to later
    
    # Shorten to lists of key tuples so don't have to iterate over every population twice for every timestep!
    agetransitlist = []
    risktransitlist = []
    for p1 in range(npops):
            for p2 in range(npops):
                if agetransit[p1,p2]: agetransitlist.append((p1,p2))
                if risktransit[p1,p2]: risktransitlist.append((p1,p2))
    
    # Figure out which populations have age inflows -- don't force population
    ageinflows   = agetransit.sum(axis=0) # Find populations with age inflows
    birthinflows = birthtransit.sum(axis=0) # Find populations with age inflows
    noinflows = findinds(ageinflows+birthinflows==0) # Find populations with no inflows
    
    
    
    
    
    
    #################################################################################################################
    ### Set initial epidemic conditions 
    #################################################################################################################
    
    # NB, to debug, use: for h in range(len(settings.statelabels)): print(settings.statelabels[h], sum(initpeople[h,:]))
    
    # Set parameters
    durationpreaids = 8.0 # Assumed duration of undiagnosed HIV pre-AIDS...used for calculating ratio of diagnosed to undiagnosed. WARNING, KLUDGY
    efftreatmentrate = 0.1 # Inverse of average duration of treatment in years...I think
    initpropcare = 0.8 # roughly estimating equilibrium proportion of diagnosed people in care
    initproplost = 0.3 # roughly estimating equilibrium proportion of people on treatment who are lost to follow-up

    # Shorten key variables
    initpeople = zeros((nstates, npops)) 
    allinfected = simpars['popsize'][:,0] * simpars['initprev'][:] # Set initial infected population
    
    # Can calculate equilibrium for each population separately
    for p in range(npops):
        # Set up basic calculations
        popinfected = allinfected[p]
        uninfected = simpars['popsize'][p,0] - popinfected # Set initial susceptible population -- easy peasy! -- should this have F['popsize'] involved?
        
        # Treatment & treatment failure
        fractotal =  popinfected / sum(allinfected) # Fractional total of infected people in this population
        treatment = simpars['numtx'][0] * fractotal # Number of people on 1st-line treatment
        if debug and treatment>popinfected: # More people on treatment than ever infected, uh oh!
            errormsg = 'More people on treatment (%f) than infected (%f)!' % (treatment, popinfected)
            if die: raise OptimaException(errormsg)
            else:
                printv(errormsg, 1, verbose)
                treatment = popinfected
        
        # Diagnosed & undiagnosed
        nevertreated = popinfected - treatment
        fracundiagnosed = exp(-durationpreaids*simpars['hivtest'][p,0])
        undiagnosed = nevertreated * fracundiagnosed     
        diagnosed = nevertreated * (1-fracundiagnosed)
        
        # Set rates within
        progratios = hstack([prog, simpars['deathlt50']]) # For last rate, use CD4<50 death as dominant rate
        progratios = (1/progratios)  / sum(1/progratios) # Normalize
        recovratios = hstack([inf, recov, efftreatmentrate]) # Not sure if this is right...inf since no progression to acute, treatmentrate since main entry here -- check
        recovratios = (1/recovratios)  / sum(1/recovratios) # Normalize
        
        # Final calculations
        undiagnosed *= progratios
        diagnosed *= progratios
        treatment *= recovratios
        
        # Populated equilibrated array
        initpeople[susreg, p]   = uninfected
        initpeople[progcirc, p] = 0.0 # This is just to make it explicit that the circ compartment only keeps track of people who are programmatically circumcised while the model is running
        initpeople[undx, p]     = undiagnosed
        if usecascade:
            initpeople[dx,   p] = diagnosed*(1.-initpropcare)
            initpeople[care, p] = diagnosed*initpropcare
            initpeople[usvl, p] = treatment * (1.-treatvs[0]) * (1.-initproplost)
            initpeople[svl,  p] = treatment * treatvs[0]      * (1.-initproplost)
            initpeople[off,  p] = treatment * initproplost * stoppropcare
            initpeople[lost, p] = treatment * initproplost * (1.-stoppropcare)
        else:
            initpeople[dx, p]   = diagnosed
            initpeople[tx, p]   = treatment
    
        if debug and not((initpeople>=0).all()): # If not every element is a real number >0, throw an error
            errormsg = 'Non-positive people found during epidemic initialization! Here are the people:\n%s' % initpeople
            if die: raise OptimaException(errormsg)
            else:
                printv(errormsg, 1, verbose)
                initpeople[initpeople<0] = 0.0
            
    people[:,:,0] = initpeople # No it hasn't, so run equilibration
    
    
    
    ##################################################################################################################
    ### Compute the effective numbers of acts outside the time loop
    ##################################################################################################################
    sexactslist = []
    injactslist = []
    
    # Sex
    for act in ['reg','cas','com']:
        for key in simpars['acts'+act]:
            this = odict()
            this['acts'] = simpars['acts'+act][key]
            if simpars['cond'+act].get(key) is not None:
                condkey = simpars['cond'+act][key]
            elif simpars['cond'+act].get((key[1],key[0])) is not None:
                condkey = simpars['cond'+act][(key[1],key[0])]
            else:
                errormsg = 'Cannot find condom use between "%s" and "%s", assuming there is none.' % (key[0], key[1]) # NB, this might not be the most reasonable assumption
                if die: raise OptimaException(errormsg)
                else: 
                    printv(errormsg, 1, verbose)
                    condkey = 0.0
                
            this['cond'] = 1.0 - condkey*effcondom
            this['pop1'] = popkeys.index(key[0])
            this['pop2'] = popkeys.index(key[1])
            if     male[this['pop1']] and   male[this['pop2']]: this['trans'] = (simpars['transmmi'] + simpars['transmmr'])/2.0 # Note: this looks horrible and stupid but it's correct! Ask Kedz
            elif   male[this['pop1']] and female[this['pop2']]: this['trans'] = simpars['transmfi']  
            elif female[this['pop1']] and   male[this['pop2']]: this['trans'] = simpars['transmfr']
            else: raise OptimaException('Not able to figure out the sex of "%s" and "%s"' % (key[0], key[1]))
            sexactslist.append(this)
            
            # Error checking
            for key in ['acts', 'cond']:
                if debug and not(all(this[key]>=0)):
                    errormsg = 'Invalid sexual behavior parameter "%s": values are:\n%s' % (key, this[key])
                    if die: raise OptimaException(errormsg)
                    else: 
                        printv(errormsg, 1, verbose)
                        this[key][this[key]<0] = 0.0 # Reset values
    
    # Injection
    for key in simpars['actsinj']:
        this = odict()
        this['acts'] = simpars['actsinj'][key]
        this['pop1'] = popkeys.index(key[0])
        this['pop2'] = popkeys.index(key[1])
        injactslist.append(this)
    
    # Convert from dicts to tuples to be faster
    for i,this in enumerate(sexactslist): sexactslist[i] = tuple([this['pop1'],this['pop2'],this['acts'],this['cond'],this['trans']])
    for i,this in enumerate(injactslist): injactslist[i] = tuple([this['pop1'],this['pop2'],this['acts']])
    
    
    ## Births precalculation
    birthslist = []
    for p1 in range(npops): # WARNING, should only loop over child populations
        alleligbirthrate = zeros(npts)
        for t in range(npts): # WARNING, could be made more efficient, it's just that the matrix multiplications get complicated -- npops x npts...
            alleligbirthrate[t] = sum(birthtransit[p1, :] * birth[p1, t]) # Births to diagnosed mothers eligible for PMTCT
        for p2 in range(npops): # WARNING, should only loop over female populations
            birthrates = birthtransit[p1, p2] * birth[p1, :] # WARNING, vector multiplication!!!! Need to create array
            if birthrates.any():
                birthslist.append(tuple([p1,p2,birthrates,alleligbirthrate]))
                
                
                
                
                
                
                
                
                
                
    ##################################################################################################################
    ### Run the model -- numerically integrate over time
    ##################################################################################################################

    for t in range(npts): # Loop over time
        printv('Timestep %i of %i' % (t+1, npts), 4, verbose)
        
        ## Calculate "effective" HIV prevalence -- taking diagnosis and treatment into account
        for pop in range(npops): # Loop over each population group
            allpeople[pop,t] = sum(people[:,pop,t]) # All people in this population group at this time point
            if debug and not(allpeople[pop,t]>0): 
                errormsg = 'No people in population %i at timestep %i (time %0.1f)' % (pop, t, tvec[t])
                if die: raise OptimaException(errormsg)
                else: printv(errormsg, 1, verbose)
            effundx = sum(cd4trans * people[undx,pop,t]); # Effective number of infecious undiagnosed people
            effdx   = sum(dxfactor * cd4trans * people[dx,pop,t]) # ...and diagnosed/failed -- WARNING, reinstating cd4trans because array multiplication gets ugly...but this should be fixed
            if usecascade:
                effcare = sum(dxfactor * cd4trans * people[care,pop,t]) # the diagnosis efficacy also applies to those in care??
                efftxus = sum(dxfactor * cd4trans * efftxunsupp * people[usvl,pop,t]) # ...and treated
                efftxs  = sum(dxfactor * cd4trans * efftxsupp  * people[svl,pop,t]) # ...and suppressed viral load
                efflost = sum(dxfactor * cd4trans * people[lost,pop,t]) # the diagnosis efficacy also applies to those lost to follow-up??
                effoff  = sum(dxfactor * cd4trans * people[off,pop,t])  # the diagnosis efficacy also applies to those off-ART but in care??
                # Calculate HIV "prevalence", scaled for infectiousness based on CD4 count; assume that treatment failure infectiousness is same as corresponding CD4 count
                effhivprev[pop] = (effundx+effdx+effcare+efftxus+efftxs+efflost+effoff) / allpeople[pop,t]
            else:
                efftx   = sum(dxfactor * cd4trans * txfactor[t] * people[tx,pop,t]) # ...and treated
                effhivprev[pop] = (effundx+effdx+efftx) / allpeople[pop,t] # Calculate HIV "prevalence", scaled for infectiousness based on CD4 count; assume that treatment failure infectiousness is same as corresponding CD4 count

            if debug and not(effhivprev[pop]>=0): 
                errormsg = 'HIV prevalence invalid in population %s! (=%f)' % (pop, effhivprev[pop])
                if die: raise OptimaException(errormsg)
                else:
                    printv(errormsg, 1, verbose)
                    effhivprev[pop] = 0.0
        
        ## Calculate inhomogeneity in the force-of-infection based on prevalence
        for pop in range(npops):
            c = inhomopar[pop]
            thisprev = sum(people[allplhiv,pop,t]) / allpeople[pop,t] 
            inhomo[pop] = (c+eps) / (exp(c+eps)-1) * exp(c*(1-thisprev)) # Don't shift the mean, but make it maybe nonlinear based on prevalence

        
        
        
        ###############################################################################
        ## Calculate force-of-infection (forceinf)
        ###############################################################################
        
        # Reset force-of-infection vector for each population group, handling circs and uncircs separately
        forceinfvec = zeros((len(sus), npops))
        thisforceinfsex = zeros(2)
        
        # Loop over all acts (partnership pairs) -- force-of-infection in pop1 due to pop2
        for pop1,pop2,acts,cond,thistrans in sexactslist:
            dtcondacts = dt*cond[t]*acts[t] # Make it so this only has to be calculated once
            
            if male[pop1]: # Separate FOI calcs for circs vs uncircs -- WARNING, could be shortened with a loop but maybe not simplified
                thisforceinfsex[0]     = 1 - mpow((1-thistrans*prepsticirceff[pop1,t]),   (dtcondacts*effhivprev[pop2]))
                thisforceinfsex[1]     = 1 - mpow((1-thistrans*prepsticircconst[pop1,t]), (dtcondacts*effhivprev[pop2]))
                forceinfvec[:,pop1] = 1 - (1-forceinfvec[:,pop1])   * (1-thisforceinfsex)
            else: # Only have uncircs for females
                thisforceinfsex[0] = 1 - mpow((1-thistrans*prepsti[pop1,t]), (dtcondacts*effhivprev[pop2]))
                forceinfvec[susreg,pop1] = 1 - (1-forceinfvec[susreg,pop1]) * (1-thisforceinfsex[0])
                
            if debug and not all(forceinfvec[:,pop1]>=0):
                errormsg = 'Sexual force-of-infection is invalid in population %s, time %0.1f, FOI:\n%s)' % (popkeys[pop1], tvec[t], forceinfvec)
                for var in ['thistrans', 'circeff[pop1,t]', 'prepeff[pop1,t]', 'stieff[pop1,t]', 'cond', 'acts', 'effhivprev[pop2]']:
                    errormsg += '\n%20s = %f' % (var, eval(var)) # Print out extra debugging information
                raise OptimaException(errormsg)
            
        # Injection-related infections -- force-of-infection in pop1 due to pop2
        for pop1,pop2,effinj in injactslist:
            
            thisforceinfinj = 1 - mpow((1-transinj), (dt*sharing[pop1,t]*effinj[t]*osteff[t]*effhivprev[pop2]))
            for index in sus: # Assign the same injecting FOI to circs and uncircs, as it doesn't matter
                forceinfvec[index,pop1] = 1 - (1-forceinfvec[index,pop1]) * (1-thisforceinfinj)
            
            if debug and not all(forceinfvec[:,pop1]>=0):
                errormsg = 'Injecting force-of-infection is invalid in population %s, time %0.1f, FOI:\n%s)' % (popkeys[pop1], tvec[t], forceinfvec)
                for var in ['transinj', 'sharing[pop1,t]', 'effinj', 'osteff[t]', 'effhivprev[pop2]']:
                    errormsg += '\n%20s = %f' % (var, eval(var)) # Print out extra debugging information
                raise OptimaException(errormsg)
        
        



        ##############################################################################################################
        ### The ODEs
        ##############################################################################################################
    
        ## Set up
    
        # New infections -- through pre-calculated force of infection
        newinfections = zeros((len(sus), npops)) 
        for index in sus:
            newinfections[index,:] = forceinfvec[index,:] * force * inhomo * people[index,:,t] 
    
        # Initalise / reset arrays
        dU = []; dD = []
        if usecascade: dC = []; dUSVL = []; dSVL = []; dL = []; dO = []; # Reset differences for cascade compartments
        else: dT = []; # Reset differences for simple compartments
        testingrate  = [0] * ncd4
        newdiagnoses = [0] * ncd4
        if usecascade:
            newtreat        = [0] * ncd4
            restarters      = [0] * ncd4
            newlinkcaredx   = [0] * ncd4
            newlinkcarelost = [0] * ncd4
            leavecareCD     = [0] * ncd4
            leavecareOL     = [0] * ncd4
            virallysupp     = [0] * ncd4
            failing         = [0] * ncd4
            stopUSlost      = [0] * ncd4
            stopSVLlost     = [0] * ncd4
            stopUSincare    = [0] * ncd4
            stopSVLincare   = [0] * ncd4
        else:
            newtreat        = [0] * ncd4

        background   = simpars['death'][:, t] # make OST effect this death rates
        
        ## Susceptibles
        otherdeaths = zeros((len(sus), npops)) 
        for index in sus:
            otherdeaths[index] = dt * people[sus[index],:,t] * background
            raw_otherdeath[:,t] += otherdeaths[index]/dt    # Save annual other deaths 
        dS = -newinfections - otherdeaths # Change in number of susceptibles -- death rate already taken into account in pm.totalpop and dt
        raw_inci[:,t] = (newinfections.sum(axis=0) + raw_mtct[:,t])/float(dt)  # Store new infections AND new MTCT births

        ## Undiagnosed
        if not(isnan(propdx[t])):
            currplhiv = people[allplhiv,:,t].sum(axis=0)
            currdx = people[alldx,:,t].sum(axis=0)
            currundx = currplhiv[:] - currdx[:]
            fractiontodx = maximum(0, (propdx[t]*currplhiv[:] - currdx[:])/(currundx[:] + eps)) # Don't allow to go negative -- note, this equation is right, I just checked it!

        for cd4 in range(ncd4):
            if cd4>0: 
                progin = dt*prog[cd4-1]*people[undx[cd4-1],:,t]
            else: 
                progin = 0 # Cannot progress into acute stage
            if cd4<ncd4-1: 
                progout = dt*prog[cd4]*people[undx[cd4],:,t]
                testingrate[cd4] = hivtest[:,t] # Population specific testing rates
                if cd4>=aidsind:
                    testingrate[cd4] = maximum(hivtest[:,t], aidstest[t]) # Testing rate in the AIDS stage (if larger!)
            else: 
                progout = 0  # Cannot progress out of AIDS stage
                testingrate[cd4] = maximum(hivtest[:,t], aidstest[t]) # Testing rate in the AIDS stage (if larger!)
            if not(isnan(propdx[t])):
                newdiagnoses[cd4] = fractiontodx * people[undx[cd4],:,t]
            else:
                newdiagnoses[cd4] =  testingrate[cd4] * dt * people[undx[cd4],:,t]
            hivdeaths   = dt * people[undx[cd4],:,t] * death[cd4]
            otherdeaths = dt * people[undx[cd4],:,t] * background
            inflows = progin  # Add in new infections after loop
            outflows = progout + newdiagnoses[cd4] + hivdeaths + otherdeaths
            dU.append(inflows - outflows)
            raw_diag[:,t]    += newdiagnoses[cd4]/dt # Save annual diagnoses 
            raw_death[:,t] += hivdeaths/dt    # Save annual HIV deaths 
            raw_otherdeath[:,t] += otherdeaths/dt    # Save annual other deaths 

        dU[0] = dU[0] + newinfections.sum(axis=0) # Now add newly infected people
        



        ############################################################################################################
        ## Here, split and decide whether or not to use the cascade for the rest of the ODEs to solve
        ############################################################################################################
        if usecascade:

            ## Diagnosed
            if not(isnan(propcare[t])):
                curralldx = people[alldx,:,t].sum(axis=0)
                currcare  = people[allcare,:,t].sum(axis=0)
                curruncare = curralldx[:] - currcare[:]
                fractiontocare = (propcare[t]*curralldx[:] - currcare[:])/(curruncare[:] + eps)
                fractiontocare = maximum(0, fractiontocare) # Don't allow to go negative -- note, this equation is right, I just checked it!
                fractiontocare = minimum(safetymargin, fractiontocare) # Cap at safetymargin rate
    
            for cd4 in range(ncd4):
                if cd4>0: 
                    progin = dt*prog[cd4-1]*people[dx[cd4-1],:,t]
                else: 
                    progin = 0 # Cannot progress into acute stage
                if cd4<ncd4-1: 
                    progout = dt*prog[cd4]*people[dx[cd4],:,t]
                else: 
                    progout = 0 # Cannot progress out of AIDS stage
                hivdeaths   = dt * people[dx[cd4],:,t] * death[cd4]
                otherdeaths = dt * people[dx[cd4],:,t] * background
                if not(isnan(propcare[t])):
                    newlinkcaredx[cd4]   = fractiontocare * people[dx[cd4],:,t] # diagnosed moving into care
                    newlinkcarelost[cd4] = fractiontocare * people[lost[cd4],:,t] # lost moving into care
                else:
                    newlinkcaredx[cd4]   = linktocare[:,t] * dt * people[dx[cd4],:,t] # diagnosed moving into care
                    newlinkcarelost[cd4] = linktocare[:,t] * dt * people[lost[cd4],:,t] # lost moving into care
                inflows = progin + newdiagnoses[cd4]*(1.-immediatecare[:,t]) # some go immediately into care after testing
                outflows = progout + hivdeaths + otherdeaths + newlinkcaredx[cd4] # NB, only newlinkcaredx flows out from here!
                dD.append(inflows - outflows)
                raw_death[:,t]  += hivdeaths/dt # Save annual HIV deaths 
                raw_otherdeath[:,t] += otherdeaths/dt    # Save annual other deaths 
            

            ## In care
            currentincare = people[care,:,t] # how many people currently in care (by population)

            if not(isnan(proptx[t])): # WARNING, newtreat should remove people not just from 'care' but also from 'off'
                currcare = people[allcare,:,t].sum(axis=0) # This assumed proptx referes to the proportion of diagnosed who are to be on treatment 
                currtx = people[alltx,:,t].sum(axis=0)
                totnewtreat =  max(0,(proptx[t]*currcare - currtx).sum()) # this is not meant to be split by population -- WARNING, not sure about max
            else:
                totnewtreat = max(0,numtx[t] - people[alltx,:,t].sum()) # Calculate difference between current people on treatment and people needed
                
            for cd4 in reversed(range(ncd4)):  # Going backwards so that lower CD4 counts move onto treatment first
                if cd4>0: 
                    progin = dt*prog[cd4-1]*people[care[cd4-1],:,t]
                else: 
                    progin = 0 # Cannot progress into acute stage
                if cd4<ncd4-1: 
                    progout = dt*prog[cd4]*people[care[cd4],:,t]
                else: 
                    progout = 0 # Cannot progress out of AIDS stage

                hivdeaths   = dt * people[care[cd4],:,t] * death[cd4]
                otherdeaths = dt * people[care[cd4],:,t] * background
                leavecareCD[cd4] = dt * people[care[cd4],:,t] * leavecare[:,t]
                inflows = progin + newdiagnoses[cd4]*immediatecare[:,t] + newlinkcaredx[cd4] + newlinkcarelost[cd4] # People move in from both diagnosed and lost states
                outflows = progout + hivdeaths + otherdeaths + leavecareCD[cd4]

                if totnewtreat: # Move people onto treatment if there are spots available
                    thisnewtreat = min(totnewtreat, sum(currentincare[cd4,:])) # Figure out how many spots are available
                    newtreat[cd4] = thisnewtreat * (currentincare[cd4,:]) / (eps+sum(currentincare[cd4,:])) # Pull out evenly from each population
                    newtreat[cd4] = minimum(newtreat[cd4], safetymargin*(currentincare[cd4,:]+inflows-outflows)) # RS: I think it would be much nicer to do this with rates
                    totnewtreat -= thisnewtreat # Adjust the number of available treatment spots
                    totnewtreat = max(totnewtreat,0.) # Prevent it going negative

                dC.insert(0, inflows - outflows - newtreat[cd4])
                dD[cd4] += leavecareCD[cd4]
                raw_newtreat[:,t] += newtreat[cd4]/dt # Save annual treatment initiation
                raw_death[:,t]  += hivdeaths/dt # Save annual HIV deaths 
                raw_otherdeath[:,t] += otherdeaths/dt    # Save annual other deaths 
            

            ## Unsuppressed/Detectable Viral Load (having begun treatment)
            currentusupp = people[usvl,:,t] # how many with suppressed viral load
            currentsupp  = people[svl,:,t]
            if not(isnan(propsupp[t])): # WARNING this will replace consequence of viral monitoring programs
                currsupp  = currentsupp.sum(axis=0)
                currusupp = currentusupp.sum(axis=0)
                newsupptot = (propsupp[t]*currusupp - currsupp).sum()
            # 40% progress, 40% recover, 20% don't change cd4 count
            for cd4 in range(ncd4):
                if cd4>0: 
                    progin = dt*prog[cd4-1]*people[usvl[cd4-1],:,t]*progusvl
                else: 
                    progin = 0 # Cannot progress into acute stage
                if cd4<ncd4-1: 
                    progout = dt*prog[cd4]*people[usvl[cd4],:,t]*progusvl
                else: 
                    progout = 0 # Cannot progress out of AIDS stage
                if (cd4>0 and cd4<ncd4-1): # CD4>0 stops people from moving back into acute
                    recovin = dt*recov[cd4-1]*people[usvl[cd4+1],:,t]*recovusvl
                else: 
                    recovin = 0 # Cannot recover in to acute or AIDS stage
                if cd4>1: # CD4>1 stops people from moving back into acute
                    recovout = dt*recov[cd4-2]*people[usvl[cd4],:,t]*recovusvl
                else: 
                    recovout = 0 # Cannot recover out of gt500 stage (or acute stage)
                hivdeaths         = dt * people[usvl[cd4],:,t] * death[cd4] * deathtx # Use death by CD4 state if lower than death on treatment
                otherdeaths       = dt * people[usvl[cd4],:,t] * background
                if not(isnan(propsupp[t])): # WARNING this will replace consequence of viral monitoring programs
                    virallysupp[cd4] = newsupptot * currentusupp[cd4,:] / (eps+currentusupp.sum()) # pull out evenly among usupp
                else:
                    virallysupp[cd4]  = dt * people[usvl[cd4],:,t] * freqvlmon[t]
                propdead          = dt * (death[cd4]*deathtx + background)
                stopUSincare[cd4] = dt * people[usvl[cd4],:,t] * stoprate[:,t] * stoppropcare  # People stopping ART but still in care
                stopUSlost[cd4]   = dt * people[usvl[cd4],:,t] * stoprate[:,t] * (1.-stoppropcare-propdead)  # People stopping ART and lost to followup
                inflows  = progin  + recovin  + newtreat[cd4]*(1.-treatvs[t]) # NB, treatvs will take care of the last 90... 
                outflows = progout + recovout + hivdeaths + otherdeaths + stopUSincare[cd4] + stopUSlost[cd4] + virallysupp[cd4]
                dUSVL.append(inflows - outflows)
                raw_death[:,t] += hivdeaths/dt # Save annual HIV deaths 
                raw_otherdeath[:,t] += otherdeaths/dt    # Save annual other deaths 
            

            ## Suppressed Viral Load
            currentsupp      = people[svl,:,t]
            for cd4 in range(ncd4):
                if (cd4>0 and cd4<ncd4-1): # CD4>0 stops people from moving back into acute
                    recovin = dt*recov[cd4-1]*people[svl[cd4+1],:,t]
                else: 
                    recovin = 0 # Cannot recover in to acute or AIDS stage
                if cd4>1: # CD4>1 stops people from moving back into acute
                    recovout = dt*recov[cd4-2]*people[svl[cd4],:,t]
                else: 
                    recovout = 0 # Cannot recover out of gt500 stage (or acute stage)
                hivdeaths          = dt * currentsupp[cd4,:] * death[cd4]
                otherdeaths        = dt * currentsupp[cd4,:] * background
                failing[cd4]       = dt * currentsupp[cd4,:] * biofailure[t]
                propdead           = dt * (death[cd4] + background)
                stopSVLincare[cd4] = dt * currentsupp[cd4,:] * stoprate[:,t] * stoppropcare  # People stopping ART but still in care
                stopSVLlost[cd4]   = dt * currentsupp[cd4,:] * stoprate[:,t] * (1.-stoppropcare-propdead) # People stopping ART and lost to followup
                inflows = recovin + virallysupp[cd4] + newtreat[cd4]*treatvs[t]
                outflows = recovout + hivdeaths + otherdeaths + failing[cd4] + stopSVLincare[cd4] + stopSVLlost[cd4]
                dSVL.append(inflows - outflows)
                dUSVL[cd4] += failing[cd4]
                raw_death[:,t]  += hivdeaths/dt # Save annual HIV deaths 
                raw_otherdeath[:,t] += otherdeaths/dt    # Save annual other deaths 


            ## Lost to follow-up (and not in care)
            for cd4 in range(ncd4):
                if cd4>0: 
                    progin = dt*prog[cd4-1]*people[lost[cd4-1],:,t]
                else: 
                    progin = 0 # Cannot progress into acute stage
                if cd4<ncd4-1: 
                    progout = dt*prog[cd4]*people[lost[cd4],:,t]
                else: 
                    progout = 0 # Cannot progress out of AIDS stage
                hivdeaths   = dt * people[lost[cd4],:,t] * death[cd4]
                otherdeaths = dt * people[lost[cd4],:,t] * background
                inflows  = progin + stopSVLlost[cd4] + stopUSlost[cd4]
                outflows = progout + hivdeaths + otherdeaths + newlinkcarelost[cd4] # These people move back into care
                dL.append(inflows - outflows) 
                raw_death[:,t]  += hivdeaths/dt # Save annual HIV deaths 
                raw_otherdeath[:,t] += otherdeaths/dt    # Save annual other deaths 


            ## Off ART but in care
            for cd4 in range(ncd4):
                if cd4>0: 
                    progin = dt*prog[cd4-1]*people[off[cd4-1],:,t]
                else: 
                    progin = 0 # Cannot progress into acute stage
                if cd4<ncd4-1: 
                    progout = dt*prog[cd4]*people[off[cd4],:,t]
                else: 
                    progout = 0 # Cannot progress out of AIDS stage
                hivdeaths   = dt * people[off[cd4],:,t] * death[cd4]
                otherdeaths = dt * people[off[cd4],:,t] * background
                leavecareOL[cd4] = dt * people[off[cd4],:,t] * leavecare[:,t]
                inflows  = progin + stopSVLincare[cd4] + stopUSincare[cd4]
                outflows = progout + hivdeaths + otherdeaths + leavecareOL[cd4]
                restarters[cd4] = dt * people[off[cd4],:,t] * restarttreat[t]
                restarters[cd4] = minimum(restarters[cd4], safetymargin*(people[off[cd4],:,t]+inflows-outflows)) # Allow it to go negative
                restarters[cd4] = maximum(restarters[cd4], -safetymargin*people[usvl[cd4],:,t]/(eps+1.-treatvs[t])) # Make sure it doesn't remove everyone from the usvl treatment compartment
                restarters[cd4] = maximum(restarters[cd4], -safetymargin*people[svl[cd4],:,t]/(eps+treatvs[t])) # Make sure it doesn't remove everyone from the svl treatment compartment
                dO.append(inflows - outflows - restarters[cd4])
                dL[cd4] += leavecareOL[cd4] 
                dUSVL[cd4] += restarters[cd4]*(1.-treatvs[t])
                dSVL[cd4]  += restarters[cd4]*treatvs[t]
                raw_death[:,t]  += hivdeaths/dt # Save annual HIV deaths 
                raw_otherdeath[:,t] += otherdeaths/dt    # Save annual other deaths 




        # Or, do not use the cascade
        else: 
            
            # WARNING, copied from above!!
            if not(isnan(proptx[t])):
                currdx = people[alldx,:,t].sum() # This assumed proptx referes to the proportion of diagnosed who are to be on treatment 
                currtx = people[alltx,:,t].sum()
                totnewtreat =  max(0,proptx[t] * currdx - currtx)
            else:
                totnewtreat = max(0, numtx[t] - people[alltx,:,t].sum()) # Calculate difference between current people on treatment and people needed
            tmpnewtreat = totnewtreat # Copy for modification later

            ## Diagnosed
            currentdiagnosed = people[dx,:,t] # Find how many people are diagnosed
            for cd4 in reversed(range(ncd4)): # Going backwards so that lower CD4 counts move onto treatment first
                if cd4>0: 
                    progin = dt*prog[cd4-1]*people[dx[cd4-1],:,t]
                else: 
                    progin = 0 # Cannot progress into acute stage
                if cd4<ncd4-1: 
                    progout = dt*prog[cd4]*people[dx[cd4],:,t]
                else: 
                    progout = 0 # Cannot progress out of AIDS stage

                hivdeaths   = dt * currentdiagnosed[cd4,:] * death[cd4] 
                otherdeaths = dt * currentdiagnosed[cd4,:] * background
                inflows = progin + newdiagnoses[cd4]
                outflows = progout + hivdeaths + otherdeaths

                if tmpnewtreat: # Move people onto treatment if there are spots available
                    thisnewtreat = min(tmpnewtreat, sum(currentdiagnosed[cd4,:])) # Figure out how many spots are available
                    newtreat[cd4] = thisnewtreat * (currentdiagnosed[cd4,:]) / (eps+sum(currentdiagnosed[cd4,:])) # Pull out evenly from each population
                    newtreat[cd4] = minimum(newtreat[cd4], safetymargin*(currentdiagnosed[cd4,:]+inflows-outflows)) # RS: I think it would be much nicer to do this with rates
                    tmpnewtreat -= thisnewtreat # Adjust the number of available treatment spots
                    tmpnewtreat = max(tmpnewtreat,0.) # Prevent it going negative

                dD.insert(0, inflows - outflows - newtreat[cd4])
                raw_newtreat[:,t] += newtreat[cd4]/dt # Save annual treatment initiation
                raw_death[:,t]  += hivdeaths/dt # Save annual HIV deaths 
                
            
            ## 1st-line treatment
            for cd4 in range(ncd4):
                if (cd4>0 and cd4<ncd4-1): # CD4>0 stops people from moving back into acute
                    recovin = dt*recov[cd4-1]*people[tx[cd4+1],:,t]
                else: 
                    recovin = 0 # Cannot recover in to acute or AIDS stage
                if cd4>1: # CD4>1 stops people from moving back into acute
                    recovout = dt*recov[cd4-2]*people[tx[cd4],:,t]
                else: 
                    recovout = 0 # Cannot recover out of gt500 stage (or acute stage)
                hivdeaths   = dt * people[tx[cd4],:,t] * death[cd4] * deathtx # Use death by CD4 state if lower than death on treatment
                otherdeaths = dt * people[tx[cd4],:,t] * background
                dT.append(recovin - recovout + newtreat[cd4] - hivdeaths - otherdeaths)
                if debug and not((people[tx[cd4],:,t]+dT[cd4] >= 0).all()):
                    errormsg = 'WARNING, Non-positive people found for treatment!\npeople[%s, :, %i] = people[%s, :, %s] = %s' % (tx[cd4], t, settings.statelabels[tx[cd4]], tvec[t], people[tx[cd4],:,t]+dT[cd4])
                    if die: raise OptimaException(errormsg)
                    else: printv(errormsg, 1, verbose=verbose)
                    
                raw_death[:,t] += hivdeaths/dt # Save annual HIV deaths 



        ##############################################################################################################
        ### Update next time point and check for errors
        ##############################################################################################################

        # Ignore the last time point, we don't want to update further
        if t<npts-1:
            change = zeros((nstates, npops))
            change[sus,:] = dS 
            for cd4 in range(ncd4): # this could be made much more efficient
                change[undx[cd4],:] = dU[cd4]
                change[dx[cd4],:]   = dD[cd4]
                if usecascade:
                    change[care[cd4],:] = dC[cd4]
                    change[usvl[cd4],:] = dUSVL[cd4]
                    change[svl[cd4],:]  = dSVL[cd4]
                    change[lost[cd4],:] = dL[cd4] 
                    change[off[cd4],:]  = dO[cd4]
                else:
                    change[tx[cd4],:]  = dT[cd4]
            people[:,:,t+1] = people[:,:,t] + change # Update people array
            
            
            
            ###############################################################################
            ## Calculate births, age transitions and mother-to-child-transmission
            ###############################################################################
            
            ## Calculate births
            for p1,p2,birthrates,alleligbirthrate in birthslist:
                thisbirthrate = birthrates[t]
                peopledx = people[alldx, p1, t].sum() # Assign to a variable since used twice
                popbirths      = thisbirthrate * people[:, p1, t].sum()
                mtctundx       = thisbirthrate * people[undx, p1, t].sum() * effmtct[t] # Births to undiagnosed mothers
                mtcttx         = thisbirthrate * people[alltx, p1, t].sum()  * pmtcteff[t] # Births to mothers on treatment
                thiseligbirths = thisbirthrate * peopledx # Births to diagnosed mothers eligible for PMTCT
            
                receivepmtct = min(numpmtct[t]*float(thiseligbirths)/(alleligbirthrate[t]*peopledx), thiseligbirths) # Births protected by PMTCT -- constrained by number eligible 
                
                mtctdx = (thiseligbirths - receivepmtct) * effmtct[t] # MTCT from those diagnosed not receiving PMTCT
                mtctpmtct = receivepmtct * pmtcteff[t] # MTCT from those receiving PMTCT
                popmtct = mtctundx + mtctdx + mtcttx + mtctpmtct # Total MTCT, adding up all components         
                
                raw_mtct[p2, t] += popmtct
                
                people[undx[0], p2, t+1] += popmtct # HIV+ babies assigned to undiagnosed compartment
                people[susreg, p2, t+1] += popbirths - popmtct  # HIV- babies assigned to uncircumcised compartment

            
            ## Age-related transitions
            for p1,p2 in agetransitlist:
                peopleaving = people[:, p1, t] * agetransit[p1,p2]
                peopleaving = minimum(peopleaving, safetymargin*people[:, p1, t]) # Ensure positive                     
                people[:, p1, t+1] -= peopleaving # Take away from pop1...
                people[:, p2, t+1] += peopleaving # ... then add to pop2
                
            
            ## Risk-related transitions
            for p1,p2 in risktransitlist:
                peoplemoving1 = people[:, p1, t] * risktransit[p1,p2]  # Number of other people who are moving pop1 -> pop2
                peoplemoving2 = people[:, p2, t] * risktransit[p1,p2] * (sum(people[:, p1, t])/sum(people[:, p2, t])) # Number of people who moving pop2 -> pop1, correcting for population size
                peoplemoving1 = minimum(peoplemoving1, safetymargin*people[:, p1, t]) # Ensure positive
                # Symmetric flow in totality, but the state distribution will ideally change.                
                people[:, p1, t+1] += peoplemoving2 - peoplemoving1
                people[:, p2, t+1] += peoplemoving1 - peoplemoving2
            
            
            
            
            
            
            ###############################################################################
            ## Reconcile things
            ###############################################################################
            
            # Reconcile population sizes for populations with no inflows
            thissusreg = people[susreg,noinflows,t+1] # WARNING, will break if susreg is not a scalar index!
            thisprogcirc = people[progcirc,noinflows,t+1]
            allsus = thissusreg+thisprogcirc
            newpeople = popsize[noinflows,t+1] - people[:,:,t+1][:,noinflows].sum(axis=0) # Number of people to add according to simpars['popsize'] (can be negative)
            people[susreg,noinflows,t+1]   += newpeople*thissusreg/allsus # Add new people
            people[progcirc,noinflows,t+1] += newpeople*thisprogcirc/allsus # Add new people
            
            # Handle circumcision
            circppl = maximum(0, minimum(numcirc[noinflows,t], safetymargin*people[susreg,noinflows,t+1])) # Don't circumcise more people than are available
            people[susreg,noinflows,t+1]   -= circppl
            people[progcirc,noinflows,t+1] += circppl # And add these people into the circumcised compartment
            
            # Check population sizes are correct
            actualpeople = people[:,:,t+1][:,noinflows].sum()
            wantedpeople = popsize[noinflows,t+1].sum()
            if debug and abs(actualpeople-wantedpeople)>1.0: # Nearest person is fiiiiine
                errormsg = 'model(): Population size inconsistent at time t=%f: %f vs. %f' % (tvec[t+1], actualpeople, wantedpeople)
                raise OptimaException(errormsg)
            
            # Check no negative people
            if debug and not((people[:,:,t+1]>=0).all()): # If not every element is a real number >0, throw an error
                for errstate in range(nstates): # Loop over all heath states
                    for errpop in range(npops): # Loop over all populations
                        if not(people[errstate,errpop,t+1]>=0):
                            errormsg = 'WARNING, Non-positive people found!\npeople[%i, %i, %i] = people[%s, %s, %s] = %s' % (errstate, errpop, t+1, settings.statelabels[errstate], popkeys[errpop], tvec[t+1], people[errstate,errpop,t+1])
                            if die: raise OptimaException(errormsg)
                            else:
                                printv(errormsg, 1, verbose=verbose)
                                people[errstate,errpop,t+1] = 0.0 # Reset
                
    # Append final people array to sim output
    if not (people>=0).all(): raise OptimaException('Non-positive people found!')
    
    raw               = odict()    # Sim output structure
    raw['tvec']       = tvec
    raw['popkeys']    = popkeys
    raw['people']     = people
    raw['inci']       = raw_inci
    raw['mtct']       = raw_mtct
    raw['diag']       = raw_diag
    raw['newtreat']   = raw_newtreat
    raw['death']      = raw_death
    raw['otherdeath'] = raw_otherdeath
    
    return raw # Return raw results





def runmodel(project=None, simpars=None, pars=None, parset=None, progset=None, budget=None, coverage=None, budgetyears=None, settings=None, start=2000, end=2030, dt=0.2, tvec=None, name=None, uid=None, data=None, debug=False, verbose=2):
    ''' 
    Convenience function for running the model. Requires input of either "simpars" or "pars"; and for including the data,
    requires input of either "project" or "data". All other inputs are optional.
    
    Version: 2016jan23 by cliffk    
    '''
    if simpars is None:
        if pars is None: raise OptimaException('runmodel() requires either simpars or pars input; neither was provided')
        simpars = makesimpars(pars, start=start, end=end, dt=dt, tvec=tvec, name=name, uid=uid)
    if settings is None:
        try: settings = project.settings 
        except: raise OptimaException('Could not get settings from project "%s" supplied to runmodel()' % project)
    try:
        raw = model(simpars=simpars, settings=settings, debug=debug, verbose=verbose) # RUN OPTIMA!!
    except: 
        printv('Running model failed; running again with debugging...', 1, verbose)
        raw = model(simpars=simpars, settings=settings, debug=True, verbose=verbose) # If it failed, run again, with tests
    results = Resultset(project=project, raw=raw, parset=parset, progset=progset, budget=budget, coverage=coverage, budgetyears=budgetyears, simpars=simpars, data=data, domake=True) # Create structure for storing results
    return results
