from pylab import *
from numpy import *
import pymc as pm
from generic_model import *
import sys
import os

pr_type, continent, scale = sys.argv[1:]
scale = float(scale)
print pr_type, continent, scale
# AR low up PR
os.chdir('..')
from all_data import *
os.chdir('generic_analysis')

# pr_type='model_exp'
# pr_type = 'mixed'
if pr_type=='data':
    where_good = where(1-R.pr.mask)
    this_R = R[where_good]

if continent == 'Africa+':
    where_good = where(this_R['region']=='Africa+')

if continent == 'CSE Asia and Americas':
    where_good = where((this_R['region']=='America')+(this_R['region']=='CSE Asia'))

if continent == 'CSE Asia':
    where_good = where(this_R['region']=='CSE Asia')

if continent == 'Americas':
    where_good = where(this_R['region']=='America')

if continent == 'All':
    where_good=slice(None)

this_R = this_R[where_good].data

dbname = continent.replace(' ','_')+'_scale_'+str(scale)+'_'+pr_type

for i in xrange(50000):
    try:    
        M = pm.MCMC(make_model(this_R, curve_sub = model_salad.delta_one_spliner, curve_params=[model_salad.exp_model, 1, scale], 
                pr_type=pr_type), db='hdf5', name=dbname)
            
        nonpr_stochastics = []
        for s in M.stochastics:
            if s.__name__.find('pro_')==-1:
                nonpr_stochastics.append(s)
        M.use_step_method(pm.AdaptiveMetropolis,nonpr_stochastics,delay=200000)
        if pr_type=='unknown':
            M.use_step_method(pm.AdaptiveMetropolis, M.pr, scales={M.pr: .005*np.ones(M.pr.value.shape)}, delay=200000)
        break
    except pm.ZeroProbability:
        a,b,c = sys.exc_info()
        print b
        pass
if i==99:
    raise ValueError, 'Failed to create model'
M.isample(10000000,200000,1000)
# M.isample(10)
make_plots(M.db._h5file.root.chain0.PyMCsamples.cols, dbname, continent, this_R, pr_type)    

