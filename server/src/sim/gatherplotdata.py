"""
GATHERPLOTDATA

This file gathers all data that could be used for plotting and packs it into a
nice little convenient structure :)

Version: 2015jan06 by cliffk
"""

# Define labels
epititles = {'prev':'Prevalence', 'plhiv':'PLHIV', 'inci':'New infections', 'daly':'DALYs', 'death':'Deaths', 'dx':'Diagnoses', 'tx1':'First-line treatment', 'tx2':'Second-line treatment'}
epiylabels = {'prev':'HIV prevalence (%)', 'plhiv':'Number of PLHIV', 'inci':'New HIV infections per year', 'daly':'HIV-related DALYs per year', 'death':'AIDS-related deaths per year', 'dx':'New HIV diagnoses per year', 'tx1':'People on 1st-line treatment', 'tx2':'People on 2nd-line treatment'}

def gatheruncerdata(D, R, verbose=2):
    """ Gather standard results into a form suitable for plotting with uncertainties. """
    from numpy import zeros, nan, size, array, asarray
    from bunch import Bunch as struct
    from printv import printv
    printv('Gathering epidemiology results...', 3, verbose)
    
    uncer = struct()
    uncer.__doc__ = 'Output structure containing everything that might need to be plotted'
    uncer.tvec = R.tvec.tolist() # Copy time vector
    uncer.poplabels = D.G.meta.pops.short
    uncer.colorm = (0,0.3,1) # Model color
    uncer.colord = (0,0,0) # Data color
    uncer.legend = ('Model', 'Data')
    uncer.xdata = D.data.epiyears
    ndatayears = len(uncer.xdata)
    
    for key in ['prev', 'plhiv', 'inci', 'daly', 'death', 'dx', 'tx1', 'tx2']:
        percent = 100 if key=='prev' else 1 # Whether to multiple results by 100
        
        uncer[key] = struct()
        uncer[key].pops = [struct() for p in range(D.G.npops)]
        uncer[key].tot = struct()
        if key!='prev': # For stacked area plots -- an option for everything except prevalence
            uncer[key].popstacked = struct()
            uncer[key].popstacked.pops = []
            uncer[key].popstacked.legend = []
            uncer[key].popstacked.title = epititles[key]
            uncer[key].popstacked.ylabel = epiylabels[key]
        for p in range(D.G.npops):
            uncer[key].pops[p].best = (R[key].pops[0][p,:]*percent).tolist()
            uncer[key].pops[p].low = (R[key].pops[1][p,:]*percent).tolist()
            uncer[key].pops[p].high = (R[key].pops[2][p,:]*percent).tolist()
            uncer[key].pops[p].title = epititles[key] + ' - ' + D.G.meta.pops.short[p]
            uncer[key].pops[p].ylabel = epiylabels[key]
            if key!='prev':
                uncer[key].popstacked.pops.append(uncer[key].pops[p].best)
                uncer[key].popstacked.legend.append(D.G.meta.pops.short[p])
        uncer[key].tot.best = (R[key].tot[0]*percent).tolist()
        uncer[key].tot.low = (R[key].tot[1]*percent).tolist()
        uncer[key].tot.high = (R[key].tot[2]*percent).tolist()
        uncer[key].tot.title = epititles[key] + ' - Overall'
        uncer[key].tot.ylabel = epiylabels[key]
        uncer[key].xlabel = 'Years'
        
        if key=='prev':
            epidata = D.data.key.hivprev[0] # TODO: include uncertainties
            uncer.prev.ydata = zeros((D.G.npops,ndatayears)).tolist()
        if key=='plhiv':
            epidata = nan+zeros(ndatayears) # No data
            uncer.daly.ydata = zeros(ndatayears).tolist()
        if key=='inci':
            epidata = D.data.opt.numinfect[0]
            uncer.inci.ydata = zeros(ndatayears).tolist()
        if key=='death':
            epidata = D.data.opt.death[0]
            uncer.death.ydata = zeros(ndatayears).tolist()
        if key=='daly':
            epidata = nan+zeros(ndatayears) # No data
            uncer.daly.ydata = zeros(ndatayears).tolist()
        if key=='dx':
            epidata = D.data.opt.numdiag[0]
            uncer.dx.ydata = zeros(ndatayears).tolist()
        if key=='tx1':
            epidata = D.data.txrx.numfirstline[0]
            uncer.tx1.ydata = zeros(ndatayears).tolist()
        if key=='tx2':
            epidata = D.data.txrx.numsecondline[0]
            uncer.tx2.ydata = zeros(ndatayears).tolist()


        if size(epidata[0])==1: # TODO: make this less shitty, easier way of checking what shape the data is I'm sure
            uncer[key].ydata = (array(epidata)*percent).tolist()
        elif size(epidata)==D.G.npops:
            for p in range(D.G.npops):
                thispopdata = epidata[p]
                if len(thispopdata) == 1: 
                    thispopdata = nan+zeros(ndatayears) # If it's an assumption, just set with nans
                elif len(thispopdata) != ndatayears:
                    raise Exception('Expect data length of 1 or %i, actually %i' % (ndatayears, len(thispopdata)))
                uncer[key].ydata[p] = (asarray(thispopdata)*percent).tolist() # Stupid, but make sure it's an array, then make sure it's a list
        else:
            raise Exception("Can't figure out size of epidata; doesn't seem to be a vector or a matrix")

    
    # Financial outputs
    for key in ['costcur', 'costfut']:
        uncer[key] = struct()
        uncer[key].ann = struct()
        uncer[key].cum = struct()
        for ac in ['ann','cum']:
            if key=='costcur' and ac=='ann': origkey = 'annualhivcosts'
            if key=='costcur' and ac=='cum': origkey = 'cumulhivcosts'
            if key=='costfut' and ac=='ann': origkey = 'annualhivcostsfuture'
            if key=='costfut' and ac=='cum': origkey = 'cumulhivcostsfuture'
            uncer[key][ac].best = R[key][ac][0].tolist()
            uncer[key][ac].low = R[key][ac][1].tolist()
            uncer[key][ac].high = R[key][ac][2].tolist()
            uncer[key][ac].xdata = R['costshared'][origkey]['xlinedata'].tolist()
            uncer[key][ac].title = R['costshared'][origkey]['title']
            uncer[key][ac].xlabel = R['costshared'][origkey]['xlabel']
            uncer[key][ac].ylabel = R['costshared'][origkey]['ylabel']
            uncer[key][ac].legend = ['Model']
    
    
    printv('...done gathering uncertainty results.', 4, verbose)
    return uncer




def gathermultidata(D, Rarr, verbose=2):
    """ Gather multi-simulation results (scenarios and optimizations) into a form suitable for plotting. """
    from bunch import Bunch as struct
    from printv import printv
    printv('Gathering multi-simulation results...', 3, verbose)
    
    
    multi = struct()
    multi.__doc__ = 'Output structure containing everything that might need to be plotted'
    multi.nsims = len(Rarr) # Number of simulations
    multi.tvec = Rarr[0].R.tvec.tolist() # Copy time vector
    multi.poplabels = D.G.meta.pops.long
    
    for key in ['prev', 'plhiv', 'inci', 'daly', 'death', 'dx', 'tx1', 'tx2']:
        percent = 100 if key=='prev' else 1 # Whether to multiple results by 100
        multi[key] = struct()
        multi[key].pops = [struct() for p in range(D.G.npops)]
        for p in range(D.G.npops):
            multi[key].pops[p].data = []
            multi[key].pops[p].legend = []
            multi[key].pops[p].title = epititles[key] + ' - ' + D.G.meta.pops.short[p]
            multi[key].pops[p].ylabel = epiylabels[key]
            for sim in range(multi.nsims):
                thisdata = (Rarr[sim].R[key].pops[0][p,:]*percent).tolist()
                multi[key].pops[p].data.append(thisdata)
                multi[key].pops[p].legend.append(Rarr[sim].label)
        multi[key].tot = struct()
        multi[key].tot.data = []
        multi[key].tot.legend = []
        multi[key].tot.title = epititles[key] + ' - Overall'
        multi[key].tot.ylabel = epiylabels[key]
        multi[key].xlabel = 'Years'
        for sim in range(multi.nsims):
            thisdata =(Rarr[sim].R[key].tot[0]*percent).tolist()
            multi[key].tot.data.append(thisdata)
            multi[key].tot.legend.append(Rarr[sim].label) # Add legends
        
    
    # Financial outputs
    for key in ['costcur', 'costfut']:
        multi[key] = struct()
        for ac in ['ann','cum']:
            if key=='costcur' and ac=='ann': origkey = 'annualhivcosts'
            if key=='costcur' and ac=='cum': origkey = 'cumulhivcosts'
            if key=='costfut' and ac=='ann': origkey = 'annualhivcostsfuture'
            if key=='costfut' and ac=='cum': origkey = 'cumulhivcostsfuture'
            multi[key][ac] = struct()
            multi[key][ac].data = []
            multi[key][ac].legend = []
            for sim in range(multi.nsims):
                thisdata = Rarr[sim].R[key][ac][0].tolist()
                multi[key][ac].data.append(thisdata)
                multi[key][ac].legend.append(Rarr[sim].label) # Add legends
                multi[key][ac].xdata  = Rarr[sim].R['costshared'][origkey]['xlinedata'].tolist()
                multi[key][ac].title  = Rarr[sim].R['costshared'][origkey]['title']
                multi[key][ac].xlabel = Rarr[sim].R['costshared'][origkey]['xlabel']
                multi[key][ac].ylabel = Rarr[sim].R['costshared'][origkey]['ylabel']
        
    printv('...done gathering multi-simulation results.', 4, verbose)
    return multi


def gatheroptimdata(D, A, verbose=2):
    """ Return the data for plotting the two pie charts -- current allocation and optimal. """
    from bunch import Bunch as struct
    from printv import printv
    printv('Gathering optimization results...', 3, verbose)
    
    O = struct()
    O.legend = D.G.meta.progs.short
    
    O.pie1 = struct()
    O.pie1.name = 'Original'
    O.pie1.val = A[0].alloc.tolist()
    
    O.pie2 = struct()
    O.pie2.name = 'Optimal'
    O.pie2.val = A[1].alloc.tolist()
    
    printv('...done gathering optimization results.', 4, verbose)
    return O
