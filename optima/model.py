## Imports
from numpy import zeros, exp, maximum, minimum, inf, array, isnan, einsum, floor, ones, power as npow, concatenate as cat, interp, nan, squeeze, isinf, isfinite
from optima import OptimaException, printv, dcp, odict, findinds

def model(simpars=None, settings=None, initpeople=None, verbose=None, die=False, debug=False, label=None, startind=None):
    """
    Runs Optima's epidemiological model.
    
    Version: 1.8 (2017mar03)
    """
    
    ##################################################################################################################
    ### Setup 
    ##################################################################################################################

    # Initialize basic quantities
    
    if verbose is None:  verbose = settings.verbose # Verbosity of output
    if label is None:    label = ''
    else:                label += ': '# An optional label to add to error messages
    if startind is None: startind = 0 # Point to start from -- used with non-empty initpeople
    if simpars is None:  raise OptimaException(label+'model() requires simpars as an input')
    if settings is None: raise OptimaException(label+'model() requires settings as an input')
    printv('Running model...', 1, verbose)
    
    if initpeople is not None:
        print('WARNING, due to parameter interpolation, results are unreliable if initpeople is not None!')
    
    # Extract key items
    popkeys         = simpars['popkeys']
    npops           = len(popkeys)
    simpars         = dcp(simpars)
    tvec            = simpars['tvec']
    dt              = float(simpars['dt'])          # Shorten dt and make absolutely sure it's a float
    npts            = len(tvec)                     # Number of time points
    ncd4            = settings.ncd4                 # Shorten number of CD4 states
    nstates         = settings.nstates              # Shorten number of health states
    eps             = settings.eps                  # Define another small number to avoid divide-by-zero errors
    forcepopsize    = settings.forcepopsize         # Whether or not to force the population size to match the parameters
    treatbycd4      = settings.treatbycd4           # Whether or not to preferentially put people on treatment from lower CD4 counts
    fromto          = simpars['fromto']             # States to and from
    transmatrix     = simpars['transmatrix']        # Raw transitions matrix

    # Initialize people array
    people          = zeros((nstates, npops, npts)) # Matrix to hold everything

    # Initialize other arrays used for internatl calculations
    allpeople       = zeros((npops, npts))          # Population sizes
    effallprev      = zeros((nstates, npops))       # HIV effective prevalence (prevalence times infectiousness), by health state
    inhomo          = zeros(npops)                  # Inhomogeneity calculations

    # Initialize raw arrays -- reporting annual quantities (so need to divide by dt!)
    raw_inci        = zeros((npops, npts))          # Total incidence acquired by each population
    raw_incibypop   = zeros((npops, npts))          # Total incidence caused by each population
    raw_births      = zeros((npops, npts))          # Total number of births to each population
    raw_mtct        = zeros((npops, npts))          # Number of mother-to-child transmissions to each population
    raw_mtctfrom    = zeros((npops, npts))          # Number of mother-to-child transmissions from each population
    raw_hivbirths   = zeros((npops, npts))          # Number of births to HIV+ pregnant women
    raw_receivepmtct= zeros((npops, npts))          # Initialise a place to store the number of people in each population receiving PMTCT
    raw_diag        = zeros((npops, npts))          # Number diagnosed per timestep
    raw_newcare     = zeros((npops, npts))          # Number newly in care per timestep
    raw_newtreat    = zeros((npops, npts))          # Number initiating ART per timestep
    raw_newsupp     = zeros((npops, npts))          # Number newly suppressed per timestep
    raw_death       = zeros((nstates, npops, npts)) # Number of deaths per timestep
    raw_otherdeath  = zeros((npops, npts))          # Number of other deaths per timestep
    
    # Biological and failure parameters
    prog            = maximum(eps,1-exp(-dt/array([simpars['progacute'], simpars['proggt500'], simpars['proggt350'], simpars['proggt200'], simpars['proggt50'], 1./simpars['deathlt50']]) ))
    svlrecov        = maximum(eps,1-exp(-dt/array([inf,inf,simpars['svlrecovgt350'], simpars['svlrecovgt200'], simpars['svlrecovgt50'], simpars['svlrecovlt50']])))
    deathhiv        = array([simpars['deathacute'],simpars['deathgt500'],simpars['deathgt350'],simpars['deathgt200'],simpars['deathgt50'],simpars['deathlt50']])
    deathsvl        = simpars['deathsvl']           # Death rate whilst on suppressive ART
    deathusvl       = simpars['deathusvl']          # Death rate whilst on unsuppressive ART
    cd4trans        = array([simpars['cd4transacute'], simpars['cd4transgt500'], simpars['cd4transgt350'], simpars['cd4transgt200'], simpars['cd4transgt50'], simpars['cd4translt50']])
    background      = simpars['death']*dt           # Background death rates
    relhivdeath     = simpars['hivdeath']           # Relative HIV-related death rates
    deathprob       = zeros((nstates,npops))        # Initialise death probability array

    # Cascade-related parameters
    requiredvl      = simpars['requiredvl']                               # Number of VL tests required per year
    treatvs         = 1.-exp(-dt/(maximum(eps,simpars['treatvs'])))       # Probability of becoming virally suppressed after 1 time step
    treatfail       = simpars['treatfail']*dt                             # Probability of treatment failure in 1 time step
    linktocare      = 1.-exp(-dt/(maximum(eps,simpars['linktocare'])))    # Probability of being linked to care in 1 time step
    aidslinktocare  = 1.-exp(-dt/(maximum(eps,simpars['aidslinktocare'])))# Probability of being linked to care in 1 time step for people with AIDS
    leavecare       = simpars['leavecare']*dt                             # Proportion of people lost to follow-up per year
    aidsleavecare   = simpars['aidsleavecare']*dt                         # Proportion of people with AIDS being lost to follow-up per year
    returntocare    = 1.-exp(-dt/(maximum(eps,simpars['returntocare'])))  # Probability of being returned to care after loss to follow-up in 1 time step
    regainvs         = simpars['regainvs']                                # Proportion of people who switch regimens when found to be failing
        
    # Disease state indices
    susreg          = settings.susreg               # Susceptible, regular
    progcirc        = settings.progcirc             # Susceptible, programmatically circumcised
    sus             = settings.sus                  # Susceptible, both circumcised and uncircumcised
    nsus            = settings.nsus                 # Number of susceptible states
    undx            = settings.undx                 # Undiagnosed
    dx              = settings.dx                   # Diagnosed
    alldx           = settings.alldx                # All diagnosed
    allcare         = settings.allcare              # All in care
    alltx           = settings.alltx                # All on treatment
    allplhiv        = settings.allplhiv             # All PLHIV
    notonart        = settings.notonart             # All PLHIV who are not on ART    
    dxnotincare     = settings.dxnotincare          # Diagnosed + lost people not in care
    care            = settings.care                 # In care
    usvl            = settings.usvl                 # On treatment - Unsuppressed Viral Load
    svl             = settings.svl                  # On treatment - Suppressed Viral Load
    lost            = settings.lost                 # Not on ART (anymore) and lost to follow-up
    acute           = settings.acute                # Acute
    gt500           = settings.gt500                # >500
    gt350           = settings.gt350                # >350
    gt200           = settings.gt200                # >200
    gt50            = settings.gt50                 # >50
    lt50            = settings.lt50                 # <50
    aidsind         = settings.aidsind              # AIDS
    allcd4          = [acute,gt500,gt350,gt200,gt50,lt50]
    
    if debug and len(sus)!=2:
        errormsg = label + 'Definition of susceptibles has changed: expecting regular circumcised + VMMC, but actually length %i' % len(sus)
        raise OptimaException(errormsg)

    # Births, deaths and transitions
    birth           = simpars['birth']*dt           # Multiply birth rates by dt
    agetransit      = simpars['agetransit']         # Multiply age transition rates by dt
    risktransit     = simpars['risktransit']        # Don't multiply risk trnsitions by dt! These are stored as the mean number of years before transitioning, and we incorporate dt later
    birthtransit    = simpars['birthtransit']       # Don't multiply the birth transitions by dt as have already multiplied birth rates by dt
    
    # Shorten to lists of key tuples so don't have to iterate over every population twice for every timestep
    risktransitlist,agetransitlist = [],[]
    for p1 in range(npops):
        for p2 in range(npops):
            if agetransit[p1,p2]:   agetransitlist.append((p1,p2, (1.-exp(-dt/agetransit[p1,p2]))))
            if risktransit[p1,p2]:  risktransitlist.append((p1,p2, (1.-exp(-dt/risktransit[p1,p2]))))
    
    # Figure out which populations have age inflows -- don't force population
    ageinflows   = agetransit.sum(axis=0)               # Find populations with age inflows
    birthinflows = birthtransit.sum(axis=0)             # Find populations with birth inflows
    noinflows = findinds(ageinflows+birthinflows==0)    # Find populations with no inflows

    # Begin calculation of transmission probabilities -- continued in time loop
    cd4trans *= simpars['transnorm']                    # Normalize CD4 transmission
    dxfactor = (1.-simpars['effdx'])                    # Include diagnosis efficacy
    efftxunsupp = (1.-simpars['efftxunsupp'])*dxfactor  # Reduction in transmission probability for usVL
    efftxsupp = (1.-simpars['efftxsupp'])*dxfactor      # Reduction in transmission probability for sVL
    alltrans = zeros(nstates)
    alltrans[undx] = cd4trans
    alltrans[dx] = cd4trans*dxfactor
    alltrans[care] = cd4trans*dxfactor
    alltrans[usvl] = cd4trans*dxfactor*efftxunsupp
    alltrans[svl] = cd4trans*dxfactor*efftxsupp
    alltrans[lost] = cd4trans*dxfactor
    
    # Proportion aware and treated (for 90/90/90) -- these are the only arrays that get used directly, so have to be dcp'd
    propdx      = dcp(simpars['propdx'])
    propcare    = dcp(simpars['propcare'])
    proptx      = dcp(simpars['proptx'])
    propsupp    = dcp(simpars['propsupp'])
    proppmtct   = dcp(simpars['proppmtct'])
    
    #  Years to fix proportions
    def findfixind(fixyearname):
        fixyearpar = simpars[fixyearname]
        tvec = simpars['tvec']
        if isnan(fixyearpar) or     fixyearpar is None: fixind = nan # It's not defined, skip
        elif fixyearpar > tvec[-1]: fixind = nan # It's after the end, skip
        elif fixyearpar < tvec[0]:  fixind = 0 # It's before the beginning, set to beginning
        else:                       fixind = findinds(tvec>=fixyearpar)[0] # Main usage case
        return fixind
        
    fixpropdx      = findfixind('fixpropdx')
    fixpropcare    = findfixind('fixpropcare')
    fixproptx      = findfixind('fixproptx')
    fixpropsupp    = findfixind('fixpropsupp')
    
#    # These all have the same format, so we put them in tuples of (proptype, data structure for storing output, state below, state in question, states above (including state in question), numerator, denominator, data structure for storing new movers)
#    #                    name,       prop,    lower, to,    num,     denom,   raw_new,        fixyear
    propstruct = odict([('propdx',   [propdx,   undx, dx,    alldx,   allplhiv, raw_diag,       fixpropdx]),
                        ('propcare', [propcare, dx,   care,  allcare, alldx,    raw_newcare,    fixpropcare]),
                        ('proptx',   [proptx,   care, alltx, alltx,   allcare,  raw_newtreat,   fixproptx]),
                        ('propsupp', [propsupp, usvl, svl,   svl,     alltx,    raw_newsupp,    fixpropsupp])])
    
    # Population sizes
    popsize = dcp(simpars['popsize'])
    
    # Population characteristics
    male    = simpars['male']       # Boolean array, true for males
    female  = simpars['female']     # Boolean array, true for females
    injects = simpars['injects']    # Boolean array, true for PWID

    # Intervention uptake (P=proportion, N=number)
    sharing   = simpars['sharing']      # Sharing injecting equiptment (P)
    numtx     = simpars['numtx']        # 1st line treatement (N) -- tx already used for index of people on treatment [npts]
    numvlmon  = simpars['numvlmon']     # Number of viral load tests done per year (N)
    hivtest   = simpars['hivtest']*dt   # HIV testing (P) [npop,npts]
    aidstest  = simpars['aidstest']*dt  # HIV testing in AIDS stage (P) [npts]
    numcirc   = (simpars['numcirc']*dt)/(1-simpars['propcirc'])   # Number of programmatic circumcisions performed (N)
        # divided by (1-propcirc) to account for propcirc proportion of people in sus already being circumcised and not needing programmatic circumcision
        # this results in an increase in circumcised people matching the desired number.
    numpmtct  = simpars['numpmtct']     # Number of people receiving PMTCT (N)
    
    # Uptake of OST
    numost    = simpars['numost']                  # Number of people on OST (N)
    if any(injects):
        numpwid = popsize[injects,:].sum(axis=0)  # Total number of PWID
        try: 
            ostprev = numost/numpwid # Proportion of PWID on OST (P)
            ostprev = minimum(ostprev, 1.0) # Don't let more than 100% of PWID be on OST :)
        except: 
            errormsg = label + 'Cannot divide by the number of PWID (numost=%f, numpwid=5f' % (numost, numpwid)
            if die: raise OptimaException(errormsg)
            else:   printv(errormsg, 1, verbose)
            ostprev = zeros(npts) # Reset to zero
    else: # No one injects
        if numost.sum(): 
            errormsg = label + 'You have entered non-zero value for the number of PWID on OST, but you have not specified any populations who inject'
            if die: raise OptimaException(errormsg)
            else:   printv(errormsg, 1, verbose)
            ostprev = zeros(npts)
        else: # No one on OST
            ostprev = zeros(npts)
    
    # Other interventions
    transinj  = simpars['transinj']                         # Injecting transition probability
    effcondom = simpars['effcondom']                        # Condom effect
    circeff   = 1 - simpars['propcirc']*simpars['effcirc']  # Circumcision efficacy in group where a certain proportion of people are circumcised (susreg group)
    circconst = 1 - simpars['effcirc']                      # Circumcision efficacy in group where everyone is circumcised (progcirc group)
    prepeff   = 1 - simpars['effprep']*simpars['prep'] - simpars['effpep']*minimum(simpars['pep'], 1 - simpars['prep'])      # PEP + PrEP effect, assuming no overlap ->  if total of PEP + PrEP > 1 then reduce PEP
    osteff    = 1 - simpars['effost']*ostprev               # OST effect
    stieff    = 1 + simpars['effsti']*simpars['stiprev']    # STI effect
    effmtct   = simpars['mtctbreast']*simpars['breast'] + simpars['mtctnobreast']*(1-simpars['breast']) # Effective MTCT transmission
    pmtcteff  = (1 - simpars['effpmtct']) * effmtct         # Effective MTCT transmission whilst on PMTCT
    allcirceff = einsum('i,j',[1,circconst],male)+einsum('i,j',[1,1],female)
    alleff = einsum('ab,ab,ab,ca->abc',prepeff,stieff,circeff,allcirceff)

    # Force of infection metaparameter
    force = simpars['force']
    inhomopar = simpars['inhomo']



    ##################################################################################################################
    ### Make time-constant, dt-dependent transitions 
    ##################################################################################################################

    ## Progression and deaths for people not on ART
    for fromstate in notonart:
        fromhealthstate = [(fromstate in j) for j in allcd4].index(True) # CD4 count of fromstate
        for tostate in fromto[fromstate]: # Iterate over the states you could be going to  
            if fromstate not in lt50: # Cannot progress from this state
                if any([(tostate in j) and (fromstate in j) for j in allcd4]):
                    transmatrix[fromstate,tostate,:] *= 1.-prog[fromhealthstate]
                else:
                    transmatrix[fromstate,tostate,:] *= prog[fromhealthstate]
    
            # Death probabilities
            transmatrix[fromstate,tostate,:] *= 1.-deathhiv[fromhealthstate]*relhivdeath*dt 
            deathprob[fromstate,:] = deathhiv[fromhealthstate]*relhivdeath*dt
            
    ## Recovery and deaths for people on suppressive ART
    for fromstate in svl:
        fromhealthstate = [(fromstate in j) for j in allcd4].index(True) # CD4 count of fromstate
        for tostate in fromto[fromstate]: # Iterate over the states you could be going to  
            if (fromstate not in acute) and (fromstate not in gt500): # You don't recover from these states
                if any([(tostate in j) and (fromstate in j) for j in allcd4]):
                    transmatrix[fromstate,tostate,:] = 1.-svlrecov[fromhealthstate]
                else:
                    transmatrix[fromstate,tostate,:] = svlrecov[fromhealthstate]
            if fromstate in acute: # You can progress from acute
                if tostate in acute:
                    transmatrix[fromstate,tostate,:] = 1.-prog[0]
                elif tostate in gt500:
                    transmatrix[fromstate,tostate,:] = prog[0]
    
            # Death probabilities
            transmatrix[fromstate,tostate,:] *= (1.-deathhiv[fromhealthstate]*relhivdeath*deathsvl*dt)    
            deathprob[fromstate,:] = deathhiv[fromhealthstate]*relhivdeath*deathsvl*dt
            

    # Recovery and progression and deaths for people on unsuppressive ART
    for fromstate in usvl:
        fromhealthstate = [(fromstate in j) for j in allcd4].index(True) # CD4 count of fromstate
    
        # Iterate over the states you could be going to  
        for tostate in fromto[fromstate]:
            if fromstate in acute: # You can progress from acute
                if tostate in acute:
                    transmatrix[fromstate,tostate,:] = 1.-prog[0]
                elif tostate in gt500:
                    transmatrix[fromstate,tostate,:] = prog[0]
            elif fromstate in gt500: 
                if tostate in gt500:
                    transmatrix[fromstate,tostate,:] = 1.-simpars['usvlproggt500']*dt
                elif tostate in gt350:
                    transmatrix[fromstate,tostate,:] = simpars['usvlproggt500']*dt
            elif fromstate in gt350:
                if tostate in gt500:
                    transmatrix[fromstate,tostate,:] = simpars['usvlrecovgt350']*dt
                elif tostate in gt350:
                    transmatrix[fromstate,tostate,:] = 1.-simpars['usvlrecovgt350']*dt-simpars['usvlproggt350']*dt
                elif tostate in gt200:
                    transmatrix[fromstate,tostate,:] = simpars['usvlproggt350']*dt
            elif fromstate in gt200:
                if tostate in gt350:
                    transmatrix[fromstate,tostate,:] = simpars['usvlrecovgt200']*dt
                elif tostate in gt200:
                    transmatrix[fromstate,tostate,:] = 1.-simpars['usvlrecovgt200']*dt-simpars['usvlproggt200']*dt
                elif tostate in gt50:
                    transmatrix[fromstate,tostate,:] = simpars['usvlproggt200']*dt
            elif fromstate in gt50:
                if tostate in gt200:
                    transmatrix[fromstate,tostate,:] = simpars['usvlrecovgt50']*dt
                elif tostate in gt50:
                    transmatrix[fromstate,tostate,:] = 1.-simpars['usvlrecovgt50']*dt-simpars['usvlproggt50']*dt
                elif tostate in lt50:
                    transmatrix[fromstate,tostate,:] = simpars['usvlproggt50']*dt
            elif fromstate in lt50:
                if tostate in gt50:
                    transmatrix[fromstate,tostate,:] = simpars['usvlrecovlt50']*dt
                elif tostate in lt50:
                    transmatrix[fromstate,tostate,:] = 1.-simpars['usvlrecovlt50']*dt
                                    
            # Death probabilities
            transmatrix[fromstate,tostate,:] *= 1.-deathhiv[fromhealthstate]*relhivdeath*deathusvl*dt
            deathprob[fromstate,:] = deathhiv[fromhealthstate]*relhivdeath*deathusvl*dt
  

  
    #################################################################################################################
    ### Set initial epidemic conditions 
    #################################################################################################################
    
    # WARNING: Set parameters, remove hard-coding
    averagedurationinfected = 10.0/2.0   # Assumed duration of undiagnosed HIV pre-AIDS...used for calculating ratio of diagnosed to undiagnosed
    averagedurationdiagnosed = 1.   # Assumed duration of diagnosed HIV pre-treatment...used for calculating ratio of lost to in care
    averagedurationincare = 3.   # Assumed duration of diagnosed HIV pre-treatment...used for calculating ratio of lost to in care

    # Check wither the initial distribution was specified
    if initpeople is not None:
        initpeople = squeeze(initpeople)
        if debug and initpeople.shape != (nstates, npops):
            errormsg = label + 'Wrong shape of init distribution: should be (%i, %i) but is %s' % (nstates, npops, initpeople.shape)
            if die: raise OptimaException(errormsg)
            else:   printv(errormsg, 1, verbose)
            initpeople = None
                
    # If it wasn't specified, or if there's something wrong with it, determine what it should be here
    if initpeople is None:

        initpeople = zeros((nstates, npops)) # Initialise
        allinfected = simpars['popsize'][:,0] * simpars['initprev'][:] # Set initial infected population
        uninfected = simpars['popsize'][:,0] - allinfected
        initnumtx = minimum(simpars['numtx'][0], allinfected.sum()/(1+eps)) # Don't allow there to be more people on treatment than infected
        if sum(allinfected): fractotal = allinfected / sum(allinfected) # Fractional total of infected people in this population
        else:                fractotal = zeros(npops) # If there's no one infected, reset to 0
        treatment = initnumtx * fractotal # Number of people on 1st-line treatment
        if debug and any(treatment>allinfected): # More people on treatment than ever infected
            errormsg = label + 'More people on treatment (%f) than infected (%f)!' % (treatment, allinfected)
            if die: raise OptimaException(errormsg)
            else:   printv(errormsg, 1, verbose)
            treatment = maximum(allinfected, treatment)
                
        treatment = initnumtx * fractotal # Number of people on 1st-line treatment
        nevertreated = allinfected - treatment

        # Set initial distributions for cascade
        testingrates =  array([simpars['hivtest'][:,0]]*ncd4)
        linkagerates = array([1.-exp(-averagedurationdiagnosed/(maximum(eps,simpars['linktocare'][:,0])))]*ncd4)
        lossrates = array([simpars['leavecare'][:,0]]*ncd4)
        for cd4 in range(aidsind, ncd4):
            testingrates[cd4,:] = maximum(simpars['aidstest'][0],simpars['hivtest'][:,0])
            linkagerates[cd4,:] = maximum(linkagerates[cd4,:],1.-exp(-averagedurationinfected/(maximum(eps,simpars['aidslinktocare'][0]))))
            lossrates[cd4,:] = minimum(simpars['aidsleavecare'][0],simpars['leavecare'][:,0])
        dxfrac = 1.-exp(-averagedurationinfected*testingrates)
        linktocarefrac = linkagerates
        lostfrac = 1.-exp(-averagedurationincare*lossrates) # WARNING, this is not technically correct, but seems to work ok in practice 
        undxdist = 1.-dxfrac
        dxdist = dxfrac*(1.-linktocarefrac)
        incaredist = dxfrac*linktocarefrac*(1.-lostfrac)
        lostdist = dxfrac*linktocarefrac*lostfrac

        # Set initial distributions within treated & untreated 
        untxdist    = (1./prog) / sum(1./prog) # Normalize progression rates to get initial distribution
        txdist      = cat([[1.,1.], svlrecov[2:]]) # Use 1s for the first two entries so that the proportion of people on tx with acute infection is v small
        txdist      = (1./txdist)  / sum(1./txdist) # Normalize

        # Set initial distribution of PLHIV
        initundx    = einsum('ij,j,i->ij',undxdist,nevertreated,untxdist)
        initdx      = einsum('ij,j,i->ij',dxdist,nevertreated,untxdist)
        initcare    = einsum('ij,j,i->ij',incaredist,nevertreated,untxdist)
        initlost    = einsum('ij,j,i->ij',lostdist,nevertreated,untxdist)
        initusvl    = (1.-treatvs)*einsum('i,j->ji',treatment,txdist)
        initsvl     = treatvs*einsum('i,j->ji',treatment,txdist)
 
        # Populated equilibrated array
        initpeople[susreg, :]      = uninfected
        initpeople[progcirc, :]    = zeros(npops) # This is just to make it explicit that the circ compartment only keeps track of people who are programmatically circumcised while the model is running
        initpeople[undx, :]        = initundx
        initpeople[dx, :]          = initdx
        initpeople[care, :]        = initcare
        initpeople[usvl, :]        = initusvl
        initpeople[svl, :]         = initsvl
        initpeople[lost, :]        = initlost

    if debug and not((initpeople>=0).all()): # If not every element is a real number >0, throw an error
        errormsg = label + 'Non-positive people found during epidemic initialization! Here are the people:\n%s' % initpeople
        if die: raise OptimaException(errormsg)
        else:   printv(errormsg, 1, verbose)
        initpeople[initpeople<0] = 0.0
            
    people[:,:,startind] = initpeople

    
    ##################################################################################################################
    ### Compute the effective numbers of acts outside the time loop
    ##################################################################################################################
    sexactslist = []
    injactslist = []
    
    # Sex
    for act in ['reg','cas','com']:
        for key in simpars['acts'+act]:
            this = odict()
            this['wholeacts'] = floor(dt*simpars['acts'+act][key])
            this['fracacts'] = dt*simpars['acts'+act][key] - this['wholeacts'] # Probability of an additional act

            if simpars['cond'+act].get(key) is not None:
                condkey = simpars['cond'+act][key]
            elif simpars['cond'+act].get((key[1],key[0])) is not None:
                condkey = simpars['cond'+act][(key[1],key[0])]
            else:
                errormsg = label + 'Cannot find condom use between "%s" and "%s", assuming there is none.' % (key[0], key[1]) # NB, this might not be the most reasonable assumption
                if die: raise OptimaException(errormsg)
                else:   printv(errormsg, 1, verbose)
                condkey = 0.0
                    
            this['cond'] = 1.0 - condkey*effcondom
            this['pop1'] = popkeys.index(key[0])
            this['pop2'] = popkeys.index(key[1])
            if     male[this['pop1']] and   male[this['pop2']]: this['trans'] = (simpars['transmmi'] + simpars['transmmr'])/2.0 
            elif   male[this['pop1']] and female[this['pop2']]: this['trans'] = simpars['transmfi']  
            elif female[this['pop1']] and   male[this['pop2']]: this['trans'] = simpars['transmfr']
            else:
                errormsg = label + 'Not able to figure out the sex of "%s" and "%s"' % (key[0], key[1])
                printv(errormsg, 3, verbose)
                this['trans'] = (simpars['transmmi'] + simpars['transmmr'] + simpars['transmfi'] + simpars['transmfr'])/4.0 # May as well just assume all transmissions apply equally - will undersestimate if pop is predominantly biologically male and oversestimate if pop is predominantly biologically female                     
                    
            sexactslist.append(this)
            
            # Error checking
            for key in ['wholeacts', 'fracacts', 'cond']:
                if debug and not(all(this[key]>=0)):
                    errormsg = label + 'Invalid sexual behavior parameter "%s": values are:\n%s' % (key, this[key])
                    if die: raise OptimaException(errormsg)
                    else:   printv(errormsg, 1, verbose)
                    this[key][this[key]<0] = 0.0 # Reset values
                        
    # Injection
    for key in simpars['actsinj']:
        this = odict()
        this['wholeacts'] = floor(dt*simpars['actsinj'][key])
        this['fracacts'] = dt*simpars['actsinj'][key] - this['wholeacts']
        
        this['pop1'] = popkeys.index(key[0])
        this['pop2'] = popkeys.index(key[1])
        injactslist.append(this)
    
    # Convert from dicts to tuples to be faster
    for i,this in enumerate(sexactslist): sexactslist[i] = tuple([this['pop1'],this['pop2'],this['wholeacts'],this['fracacts'],this['cond'],this['trans']])
    for i,this in enumerate(injactslist): injactslist[i] = tuple([this['pop1'],this['pop2'],this['wholeacts'],this['fracacts']])
    
    
    ## Births precalculation
    birthslist = []
    for p1 in range(npops): 
        alleligbirthrate = einsum('i,j->j',birthtransit[p1, :],birth[p1, :])
        for p2 in range(npops):
            birthrates = birthtransit[p1, p2] * birth[p1, :]
            if birthrates.any():
                birthslist.append(tuple([p1,p2,birthrates,alleligbirthrate]))
    motherpops = set([thisbirth[0] for thisbirth in birthslist]) # Get the list of all populations who are mothers
    
    
    ##################################################################################################################
    ### Define error checking
    ##################################################################################################################
    
    def checkfornegativepeople(people, tind=None):
        tvec = simpars['tvec']
        if tind is None: tind = Ellipsis
        if not((people[:,:,tind]>=0).all()): # If not every element is a real number >0, throw an error
            for t in range(len(tvec)):
                for errstate in range(nstates): # Loop over all heath states
                    for errpop in range(npops): # Loop over all populations
                        if not(people[errstate,errpop,t]>=0):
                            errormsg = label + 'WARNING, Non-positive people found!\npeople[%i, %i, %i] = people[%s, %s, %s] = %s' % (errstate, errpop, t, settings.statelabels[errstate], popkeys[errpop], simpars['tvec'][t], people[errstate,errpop,t])
                            if die: raise OptimaException(errormsg)
                            else:   printv(errormsg, 1, verbose=verbose)
                            people[errstate,errpop,t] = 0.0 # Reset
                

    ##################################################################################################################
    ### Run the model 
    ##################################################################################################################
    
    # Preallocate here -- supposedly the most computationally efficient way to do this
    alltransmatrices = zeros((npts, transmatrix.shape[0], transmatrix.shape[1], transmatrix.shape[2]))
    alltransmatrices[:] = transmatrix
    
    for t in range(startind, npts): # Loop over time
        printv('Timestep %i of %i' % (t+1, npts), 4, verbose)
        
        ###############################################################################
        ## Initial steps
        ###############################################################################

        ## Pull out the transitions for this timestep
        thistransit = alltransmatrices[t]
        
        ## Calculate "effective" HIV prevalence -- taking diagnosis and treatment into account
        allpeople[:,t] = people[:, :, t].sum(axis=0)
        if debug and not( all(allpeople[:,t]>0)):
            errormsg = label + 'No people in populations %s at timestep %i (time %0.1f)' % (findinds(allpeople[:,t]<=0), t, tvec[t])
            if die: raise OptimaException(errormsg)
            else: printv(errormsg, 1, verbose)
                
        effallprev = einsum('i,ij->ij',alltrans,people[:,:,t]) / allpeople[:,t]                            
        if debug and not((effallprev[:,:]>=0).all()): 
            errormsg = label + 'HIV prevalence invalid at time %s' % (t)
            if die: raise OptimaException(errormsg)
            else:   printv(errormsg, 1, verbose)
            effallprev = minimum(effallprev,eps)
                
        ## Calculate inhomogeneity in the force-of-infection based on prevalence
        thisprev = people[allplhiv,:,t].sum(axis=0) / allpeople[:,t] 
        inhomo = (inhomopar+eps) / (exp(inhomopar+eps)-1) * exp(inhomopar*(1-thisprev)) # Don't shift the mean, but make it maybe nonlinear based on prevalence
        

        ###############################################################################
        ## Calculate probability of getting infected
        ###############################################################################
        
        # Probability of getting infected. In the first stage of construction, we actually store this as the probability of NOT getting infected
        # First dimension: infection acquired by (circumcision status). Second dimension:  infection acquired by (pop). Third dimension: infection caused by (pop). Fourth dimension: infection caused by (health/treatment state)
        forceinffull = ones((len(sus), npops, nstates, npops))

        # Loop over all acts (partnership pairs) -- probability of pop1 getting infected by pop2
        for pop1,pop2,wholeacts,fracacts,cond,thistrans in sexactslist:

            thisforceinfsex = (1-fracacts[t]*thistrans*cond[t]*einsum('a,b',alleff[pop1,t,:],effallprev[:,pop2]))
            if wholeacts[t]: thisforceinfsex  *= npow((1-thistrans*cond[t]*einsum('a,b',alleff[pop1,t,:],effallprev[:,pop2])), int(wholeacts[t]))
            forceinffull[:,pop1,:,pop2] *= thisforceinfsex 
            
            if debug and not((forceinffull[:,pop1,:,pop2]>=0).all()):
                errormsg = label + 'Sexual force-of-infection is invalid between populations %s and %s, time %0.1f, FOI:\n%s)' % (popkeys[pop1], popkeys[pop2], tvec[t], forceinffull[:,pop1,:,pop2])
                for var in ['thistrans', 'circeff[pop1,t]', 'prepeff[pop1,t]', 'stieff[pop1,t]', 'cond', 'wholeacts', 'fracacts', 'effallprev[:,pop2]']:
                    errormsg += '\n%20s = %f' % (var, eval(var)) # Print out extra debugging information
                raise OptimaException(errormsg)

        # Injection-related infections -- probability of pop1 getting infected by pop2
        for pop1,pop2,wholeacts,fracacts in injactslist:
            
            thisforceinfinj = 1-transinj*sharing[pop1,t]*osteff[t]*prepeff[pop1,t]*fracacts[t]*effallprev[:,pop2]
            if wholeacts[t]: thisforceinfinj *= npow((1-transinj*sharing[pop1,t]*osteff[t]*prepeff[pop1,t]*effallprev[:,pop2]), int(wholeacts[t]))

            for index in sus: # Assign the same probability of getting infected by injection to both circs and uncircs, as it doesn't matter
                forceinffull[index,pop1,:,pop2] *= thisforceinfinj
            
            if debug and not((forceinffull[:,pop1,:,pop2]>=0).all()):
                errormsg = label + 'Injecting force-of-infection is invalid between populations %s and %s, time %0.1f, FOI:\n%s)' % (popkeys[pop1], popkeys[pop2], tvec[t], forceinffull[:,pop1,:,pop2])
                for var in ['transinj', 'sharing[pop1,t]', 'wholeacts', 'fracacts', 'osteff[t]', 'effallprev[:,pop2]']:
                    errormsg += '\n%20s = %f' % (var, eval(var)) # Print out extra debugging information
                raise OptimaException(errormsg)
        
        # Probability of getting infected is one minus forceinffull times any scaling factors
        forceinffull  = einsum('ijkl,j,j,j->ijkl', 1.-forceinffull, force, inhomo,(1.-background[:,t]))
        infections_to = forceinffull.sum(axis=(2,3)) # Infections acquired through sex and injecting - by population who gets infected
        infections_to = minimum(infections_to, 1.0-eps-background[:,t].max()) # Make sure it never exceeds the limit

        # Add these transition probabilities to the main array
        si = susreg[0] # susreg is a single element, but needs an index since can't index a list with an array
        pi = progcirc[0] # as above 
        ui = undx[0]
        thistransit[si,si,:] *= (1.-background[:,t]) - infections_to[si] # Index for moving from sus to sus
        thistransit[si,ui,:] *= infections_to[si] # Index for moving from sus to infection
        thistransit[pi,pi,:] *= (1.-background[:,t]) - infections_to[pi] # Index for moving from circ to circ
        thistransit[pi,ui,:] *= infections_to[pi] # Index for moving from circ to infection

        # Calculate infections acquired and transmitted
        raw_inci[:,t]       = einsum('ij,ijkl->j', people[sus,:,t], forceinffull)/dt
        raw_incibypop[:,t]  = einsum('ij,ijkl->l', people[sus,:,t], forceinffull)/dt

            
        ##############################################################################################################
        ### Calculate deaths
        ##############################################################################################################

        # Adjust transition rates
        thistransit[nsus:,:,:] *= (1.-background[:,t]) 

        # Store deaths
        raw_death[:,:,t]    = einsum('ij,ij->ij', people[:,:,t], deathprob)/dt
        raw_otherdeath[:,t] = einsum('ij,j->j',   people[:,:,t], background[:,t])/dt


        ##############################################################################################################
        ### Calculate probabilities of shifting along cascade (if programmatically determined)
        ##############################################################################################################

        def userate(prop,t):
            if   t==0:             return True # Never force a proportion on the first timestep
            elif isnan(prop[t-1]): return True # The previous timestep had a rate, so keep the rate here
            else:                  return False # Neither: use a proportion instead of a rate

        # Undiagnosed to diagnosed
        if userate(propdx,t): # Need to project forward one year to avoid mismatch
            dxprob = [hivtest[:,t]]*ncd4
            for cd4 in range(aidsind, ncd4): dxprob[cd4] = maximum(aidstest[t],hivtest[:,t])
        else: dxprob = zeros(ncd4)
        
        for cd4, fromstate in enumerate(undx):
            for tostate in fromto[fromstate]:
                if tostate in undx: # Probability of not being tested
                    thistransit[fromstate,tostate,:] *= (1.-dxprob[cd4])
                else: # Probability of being tested
                    thistransit[fromstate,tostate,:] *= dxprob[cd4]
                    raw_diag[:,t] += people[fromstate,:,t]*thistransit[fromstate,tostate,:]/dt

        # Diagnosed/lost to care
        if userate(propcare,t):
            careprob = [linktocare[:,t]]*ncd4
            returnprob = [returntocare[:,t]]*ncd4
            for cd4 in range(aidsind, ncd4):
                careprob[cd4]   = maximum(aidslinktocare[t],linktocare[:,t])
                returnprob[cd4] = maximum(aidslinktocare[t],returntocare[:,t])
        else:
            careprob   = zeros(ncd4)
            returnprob = zeros(ncd4)
        for cd4ind, fromstate in enumerate(dx):  # 1 categories x 6 states per category = 6 states
            cd4 = cd4ind%ncd4 # Convert from state index to actual CD4 index
            for tostate in fromto[fromstate]:
                if tostate in dx: # Probability of not moving into care for diagnosed but not yet in care
                    thistransit[fromstate,tostate,:] *= (1.-careprob[cd4])
                else: # Probability of moving into care
                    thistransit[fromstate,tostate,:] *= careprob[cd4]
        for cd4ind, fromstate in enumerate(lost):  # 1 categories x 6 states per category = 6 states
            cd4 = cd4ind%ncd4 # Convert from state index to actual CD4 index
            for tostate in fromto[fromstate]:
                if tostate in lost: # Probability of not returning to care for diagnosed and lost to follow-up
                    thistransit[fromstate,tostate,:] *= (1.-returnprob[cd4])
                else: # Probability of returning to care
                    thistransit[fromstate,tostate,:] *= returnprob[cd4]


        # Care/USVL/SVL to lost
        if userate(propcare,t):
            lossprob = [leavecare[:,t]]*ncd4 
            for cd4 in range(aidsind, ncd4): lossprob[cd4] = minimum(aidsleavecare[t],leavecare[:,t])
        else: lossprob = zeros(ncd4)
        for cd4ind, fromstate in enumerate(allcare): # 3 categories x 6 states per category = 18 states
            cd4 = cd4ind%ncd4 # Convert from state index to actual CD4 index
            for tostate in fromto[fromstate]:
                if tostate in allcare: # Probability of not being lost and remaining in care
                    thistransit[fromstate,tostate,:] *= (1.-lossprob[cd4])
                else: # Probability of being lost
                    thistransit[fromstate,tostate,:] *= lossprob[cd4]
    
        # SVL to USVL
        usvlprob = treatfail[t] if userate(propsupp,t) else 0.
        for fromstate in svl:
            for tostate in fromto[fromstate]:
                if tostate in svl: # Probability of remaining suppressed
                    thistransit[fromstate,tostate,:] *= (1.-usvlprob)
                elif tostate in usvl: # Probability of becoming unsuppressed
                    thistransit[fromstate,tostate,:] *= usvlprob
        
        # USVL to SVL
        svlprob = min(regainvs[t]*numvlmon[t]/(eps+people[alltx,:,t].sum()*requiredvl),1) if userate(propsupp,t) else 0.
        for fromstate in usvl:
            for tostate in fromto[fromstate]:
                if tostate in usvl: # Probability of not receiving a VL test & thus remaining failed
                    thistransit[fromstate,tostate,:] *= (1.-svlprob)
                elif tostate in svl: # Probability of receiving a VL test, switching to a new regime & becoming suppressed
                    thistransit[fromstate,tostate,:] *= svlprob
        
        # Check that probabilities all sum to 1
        if debug:
            transtest = array([(abs(thistransit[j,:,:].sum(axis=0)/(1.-background[:,t])+deathprob[j]-ones(npops))>eps).any() for j in range(nstates)])
            if any(transtest):
                wrongstatesindices = findinds(transtest)
                wrongstates = [settings.statelabels[j] for j in wrongstatesindices]
                wrongprobs = array([thistransit[j,:,:].sum(axis=0)/(1.-background[:,t])+deathprob[j] for j in wrongstatesindices])
                errormsg = label + 'Transitions do not sum to 1 at time t=%f for states %s: sums are \n%s' % (tvec[t], wrongstates, wrongprobs)
                raise OptimaException(errormsg)
                
        # Check that no probabilities are less than 0
        if debug and any([(thistransit[k]<0).any() for k in range(nstates)]):
            wrongstatesindices = [k for k in range(nstates) if (thistransit[k]<0.).any()]
            wrongstates = [settings.statelabels[j] for j in wrongstatesindices]
            wrongprobs = array([thistransit[j][1] for j in wrongstatesindices])
            errormsg = label + 'Transitions are less than 0 at time t=%f for states %s: sums are \n%s' % (tvec[t], wrongstates, wrongprobs)
            raise OptimaException(errormsg)
            
        ## Shift people as required
        if t<npts-1:
            for fromstate,tostates in enumerate(fromto):
                people[tostates,:,t+1] += people[fromstate,:,t]*thistransit[fromstate,tostates,:]


        ##############################################################################################################
        ### Calculate births
        ##############################################################################################################

        # Precalculate proportion on PMTCT, whether numpmtct or proppmtct is used
        numhivpospregwomen = 0 # Initialize
        timestepsonpmtct = 1./dt # Specify the number of timesteps on which mothers are on PMTCT -- # WARNING: remove hard-coding
        fsums = dict() # Has to be a dict rather than an odict since the populations are numeric
        for p1 in motherpops: # Pull these out of the loop to speed computation
            fsumpop = people[:, p1, t] # Pull out the people who will be summed
            fsums[p1] = dict() # Use another dict, since only referring to by key
            fsums[p1]['all']       = fsumpop[:].sum()
            fsums[p1]['allplhiv']  = fsumpop[allplhiv].sum()
            fsums[p1]['undx']      = fsumpop[undx].sum()
            fsums[p1]['alldx']     = fsumpop[alldx].sum()
            fsums[p1]['alltx']     = fsumpop[alltx].sum()
        for p1,p2,birthrates,alleligbirthrate in birthslist: # p1 is mothers, p2 is children
            numhivpospregwomen += birthrates[t] * fsums[p1]['alldx'] * timestepsonpmtct # Divide by dt to get number of women
        if isnan(proppmtct[t]): calcproppmtct = numpmtct[t]/(eps+numhivpospregwomen) # Proportion on PMTCT is not specified: use number
        else:                   calcproppmtct = proppmtct[t] # Else, just use the proportion specified
        calcproppmtct = min(calcproppmtct, 1.)
        
        undxhivbirths = zeros(npops) # Store undiagnosed HIV+ births for this timestep
        dxhivbirths = zeros(npops) # Store diagnosed HIV+ births for this timestep
        
        # Calculate actual births, MTCT, and PMTCT
        for p1,p2,birthrates,alleligbirthrate in birthslist:
            thisbirthrate = birthrates[t]
            popbirths      = thisbirthrate * fsums[p1]['all']
            mtctundx       = thisbirthrate * fsums[p1]['undx'] * effmtct[t] # Births to undiagnosed mothers
            mtcttx         = thisbirthrate * fsums[p1]['alltx'] * pmtcteff[t] # Births to mothers on treatment
            thiseligbirths = thisbirthrate * fsums[p1]['alldx'] # Births to diagnosed mothers eligible for PMTCT

            mtctdx = (thiseligbirths * (1-calcproppmtct)) * effmtct[t] # MTCT from those diagnosed not receiving PMTCT
            mtctpmtct = (thiseligbirths * calcproppmtct) * pmtcteff[t] # MTCT from those receiving PMTCT
            thisreceivepmtct = thiseligbirths * calcproppmtct
            thispopmtct = mtctundx + mtctdx + mtcttx + mtctpmtct # Total MTCT, adding up all components         

            undxhivbirths[p2] += mtctundx                       # Births to add to undx  
            dxhivbirths[p2]   += (mtctdx + mtcttx + mtctpmtct)   # Births add to dx
            
            raw_receivepmtct[p1, t] += thisreceivepmtct * timestepsonpmtct
            raw_mtct[p2, t] += thispopmtct/dt
            raw_mtctfrom[p1, t] += thispopmtct/dt
            raw_births[p2, t] += popbirths/dt
            raw_hivbirths[p1, t] += thisbirthrate * fsums[p1]['allplhiv'] / dt
            
        raw_inci[:,t] += raw_mtct[:,t] # Update infections acquired based on PMTCT calculation
        raw_incibypop[:,t] += raw_mtctfrom[:,t] # Update infections caused based on PMTCT calculation

        if debug and abs(raw_inci[:,t].sum() - raw_incibypop[:,t].sum()) > eps:
            errormsg = label + 'Number of infections received (%f) is not equal to the number of infections caused (%f) at time %i' % (raw_inci[:,t].sum(), raw_incibypop[:,t].sum(), t)
            if die: raise OptimaException(errormsg)
            else: printv(errormsg, 1, verbose)

        ###############################################################################
        ## Shift numbers of people (circs, treatment, age transitions, risk transitions, prop scenarios)
        ###############################################################################
        if t<npts-1:
            
            ## Births 
            people[susreg, :, t+1] += raw_births[:,t]*dt - undxhivbirths - dxhivbirths # HIV- babies assigned to uncircumcised compartment
            people[undx[0], :, t+1] += undxhivbirths # HIV+ babies born to undiagnosed mothers assigned to undiagnosed compartment 
            people[dx[0], :, t+1]   += dxhivbirths   # HIV+ babies born to diagnosed mothers assigned to diagnosed compartment

            ## Circumcision 
            circppl = minimum(numcirc[:,t+1], people[susreg,:,t+1])
            people[susreg,:,t+1]   -= circppl
            people[progcirc,:,t+1] += circppl 


            ## Age-related transitions
            for p1,p2,thisagetransprob in agetransitlist:
                peopleleaving = people[:, p1, t+1] * thisagetransprob
                if debug and (peopleleaving > people[:, p1, t+1]).any():
                    errormsg = label + 'Age transitions between pops %s and %s at time %i are too high: the age transitions you specified say that %f%% of the population should age in a single time-step.' % (popkeys[p1], popkeys[p2], t+1, agetransit[p1, p2]*100.)
                    if die: raise OptimaException(errormsg)
                    else:   printv(errormsg, 1, verbose)
                    peopleleaving = minimum(peopleleaving, people[:, p1, t]) # Ensure positive  
                        
                                           
                people[:, p1, t+1] -= peopleleaving # Take away from pop1...
                people[:, p2, t+1] += peopleleaving # ... then add to pop2


            ## Risk-related transitions
            for p1,p2,thisrisktransprob in risktransitlist:
                peoplemoving1 = people[:, p1, t+1] * thisrisktransprob  # Number of other people who are moving pop1 -> pop2
                peoplemoving2 = people[:, p2, t+1] * thisrisktransprob * (sum(people[:, p1, t+1])/sum(people[:, p2, t+1])) # Number of people who moving pop2 -> pop1, correcting for population size
                # Symmetric flow in totality, but the state distribution will ideally change.                
                people[:, p1, t+1] += peoplemoving2 - peoplemoving1 # NOTE: this should not cause negative people; peoplemoving1 is guaranteed to be strictly greater than 0 and strictly less that people[:, p1, t+1]
                people[:, p2, t+1] += peoplemoving1 - peoplemoving2 # NOTE: this should not cause negative people; peoplemoving2 is guaranteed to be strictly greater than 0 and strictly less that people[:, p2, t+1]
            

            ###############################################################################
            ## Reconcile population sizes
            ###############################################################################
            
            # Reconcile population sizes for populations with no inflows
            thissusreg = people[susreg,noinflows,t+1] # WARNING, will break if susreg is not a scalar index!
            thisprogcirc = people[progcirc,noinflows,t+1]
            allsus = thissusreg+thisprogcirc
            if debug and not all(allsus>0): 
                errormsg = label + '100%% prevalence detected (t=%f, pop=%s)' % (t+1, array(popkeys)[findinds(allsus>0)][0])
                raise OptimaException(errormsg)
            newpeople = popsize[noinflows,t+1] - people[:,:,t+1][:,noinflows].sum(axis=0) # Number of people to add according to simpars['popsize'] (can be negative)
            people[susreg,noinflows,t+1]   += newpeople*thissusreg/allsus # Add new people
            people[progcirc,noinflows,t+1] += newpeople*thisprogcirc/allsus # Add new people
            
            # Check population sizes are correct
            actualpeople = people[:,:,t+1][:,noinflows].sum()
            wantedpeople = popsize[noinflows,t+1].sum()
            if debug and abs(actualpeople-wantedpeople)>1.0: # Nearest person is fiiiiine
                errormsg = label + 'Population size inconsistent at time t=%f: %f vs. %f' % (tvec[t+1], actualpeople, wantedpeople)
                raise OptimaException(errormsg)
            
            # If required, scale population sizes to exactly match the parameters
            if forcepopsize:
                relerr = 0.1 # Set relative error tolerance
                for p in range(npops):
                    susnotonart = cat([sus,notonart])
                    actualpeople = people[susnotonart,p,t+1].sum()
                    wantedpeople = popsize[p,t+1] - people[alltx,p,t+1].sum()
                    if actualpeople==0: raise Exception("ERROR: no people.")
                    ratio = wantedpeople/actualpeople
                    if abs(ratio-1)>relerr: # It's not OK
                        errormsg = label + 'Warning, expected population size is nowhere near calculated population size (t=%f, pop=%s, wanted=%f, actual=%f, ratio=%f)' % (t, popkeys[p], wantedpeople, actualpeople, ratio)
                        if die: raise OptimaException(errormsg)
                        else: printv(errormsg, 1, verbose=verbose)
                    people[susnotonart,p,t+1] *= ratio # It's OK, so scale to match
            

            #######################################################################################
            ## Proportions -- these happen after the Euler step, which is why it's t+1 instead of t
            #######################################################################################

            for name,proplist in propstruct.items():
                prop, lowerstate, tostate, numer, denom, raw_new, fixyear = proplist

                calcprop = people[numer,:,t].sum()/(eps+people[denom,:,t].sum()) # This is the value we fix it at
                if fixyear==t: # Fixing the proportion from this timepoint
                    naninds    = findinds(isnan(prop)) # Find the indices that are nan -- to be replaced by current values
                    infinds    = findinds(isinf(prop)) # Find indices that are infinite -- to be scaled up/down to a target value
                    finiteinds = findinds(isfinite(prop)) # Find indices that are defined
                    finiteind = npts-1 if not len(finiteinds) else finiteinds[0] # Get first finite index, or else just last point -- latter should not actually matter
                    naninds = naninds[naninds>t] # Trim ones that are less than the current point
                    infinds = infinds[infinds>t] # Trim ones that are less than the current point
                    ninterppts = len(infinds) # Number of points to interpolate over
                    if len(naninds): prop[naninds] = calcprop # Replace nans with current proportion
                    if len(infinds): prop[infinds] = interp(range(ninterppts), [0,ninterppts-1], [calcprop,prop[finiteind]]) # Replace infinities with scale-up/down
                
                # Figure out how many people we currently have...
                actual    = people[numer,:,t+1].sum() # ... in the higher cascade state
                available = people[denom,:,t+1].sum() # ... waiting to move up
                
                # Move the people who started treatment last timestep from usvl to svl
                if ~isfinite(prop[t+1]):
                    if name == 'proptx': wanted = numtx[t+1] # If proptx is nan, we use numtx
                    else:                wanted = None # If a proportion or number isn't specified, skip this
                else: # If the prop value is finite, we use it
                    wanted = prop[t+1]*available
                
                # Reconcile the differences between the number we have and the number we want
                if wanted is not None:
                    diff = wanted - actual # Wanted number minus actual number 
                    if diff>eps: # We need to move people forwards along the cascade 
                        ppltomoveup = people[lowerstate,:,t+1]
                        totalppltomoveup = ppltomoveup.sum()
                        if totalppltomoveup>eps:
                            diff = min(diff, totalppltomoveup-eps) # Make sure we don't move more people than are available
                            if name == 'proptx': # For treatment, we move people in lower CD4 states first
                                if treatbycd4:
                                    tmpdiff = diff
                                    newmovers = zeros((ncd4,npops))
                                    for cd4 in reversed(range(ncd4)): # Going backwards so that lower CD4 counts move up the cascade first
                                        if tmpdiff>eps: # Move people until you have the right proportions
                                            ppltomoveupcd4 = ppltomoveup[cd4,:]
                                            totalppltomoveupcd4 = ppltomoveupcd4.sum()
                                            if totalppltomoveupcd4>eps:
                                                tmpdiffcd4 = min(tmpdiff, totalppltomoveupcd4-eps)
                                                newmovers[cd4,:] = tmpdiffcd4*ppltomoveupcd4/totalppltomoveupcd4 # Pull out evenly from each population
                                                tmpdiff -= newmovers[cd4,:].sum() # Adjust the number of available spots
                                else:
                                    newmovers = diff*ppltomoveup/totalppltomoveup
                                # Need to handle USVL and SVL separately
                                people[care,:,t+1] -= newmovers # Shift people out of care
                                people[usvl,:,t+1]  += newmovers*(1.0-treatvs) # ... and onto treatment, according to existing proportions
                                people[svl,:,t+1]   += newmovers*treatvs # Likewise for SVL
                            else: # For everything else, we use a distribution based on the distribution of people waiting to move up the cascade
                                newmovers = diff*ppltomoveup/totalppltomoveup
                                people[lowerstate,:,t+1] -= newmovers # Shift people out of the less progressed state... 
                                people[tostate,:,t+1]    += newmovers # ... and into the more progressed state
                            raw_new[:,t+1]               += newmovers.sum(axis=0)/dt # Save new movers
                    elif diff<-eps: # We need to move people backwards along the cascade
                        ppltomovedown = people[tostate,:,t+1]
                        totalppltomovedown = ppltomovedown.sum()
                        if totalppltomovedown>eps: # To avoid having to add eps
                            diff = min(-diff, totalppltomovedown-eps) # Flip it around so we have positive people
                            newmovers = diff*ppltomovedown/totalppltomovedown
                            if name == 'proptx': # Handle SVL and USVL separately
                                newmoversusvl = newmovers[:ncd4,:] # First group of movers are from USVL
                                newmoverssvl  = newmovers[ncd4:,:] # Second group is SVL
                                people[usvl,:,t+1] -= newmoversusvl # Shift people out of USVL treatment
                                people[svl,:,t+1]  -= newmoverssvl  # Shift people out of SVL treatment
                                people[care,:,t+1] += newmoversusvl+newmoverssvl # Add both groups of movers into care
                            else:
                                people[tostate,:,t+1]    -= newmovers # Shift people out of the more progressed state... 
                                people[lowerstate,:,t+1] += newmovers # ... and into the less progressed state
                            raw_new[:,t+1]               -= newmovers.sum(axis=0)/dt # Save new movers, inverting again
            if debug: checkfornegativepeople(people, tind=t+1) # If ebugging, check for negative people on every timestep
        
    raw                 = odict()    # Sim output structure
    raw['tvec']         = tvec
    raw['popkeys']      = popkeys
    raw['people']       = people
    raw['inci']         = raw_inci
    raw['incibypop']    = raw_incibypop
    raw['mtct']         = raw_mtct
    raw['births']       = raw_births
    raw['hivbirths']    = raw_hivbirths
    raw['pmtct']        = raw_receivepmtct
    raw['diag']         = raw_diag
    raw['newtreat']     = raw_newtreat
    raw['death']        = raw_death
    raw['otherdeath']   = raw_otherdeath
    
    checkfornegativepeople(people) # Check only once for negative people, right before finishing
    
    return raw # Return raw results