import optima,liboptima
from timevarying import timevarying

r = optima.Project.load_json('./projects/Dedza.json')
opt = optima.Optimization('Example',r,calibration=r.calibrations[0],programset=r.programsets[0],initial_alloc=r.data['origalloc']) # Make a sim, implicitly selecting a calibration and programset/ccocs
opt.optimize(timelimit = 5)
optima.plot([opt.initial_sim,opt.optimized_sim])
