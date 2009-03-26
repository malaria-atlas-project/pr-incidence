# AR low up PR
from generic_model import *
from pylab import *
from numpy import *
from all_data import *
import os

pr_type='model_exp'
# pr_type = 'mixed'
# pr_type = 'data_untrans'
scale = .6
# pr_type = 'data'

continent = 'Africa+'
# continent = 'CSE Asia and Americas'
# continent = 'CSE Asia'
# continent = 'America'
# continent = 'All'

this_R = filtered_data(pr_type, continent)

dbname = continent.replace(' ','_')+'_scale_'+str(scale)+'_'+pr_type

os.chdir('traces')
for i in xrange(500000):
    try:    
        # M = pm.MCMC(make_model(this_R, curve_sub = model_salad.polynomial_model, curve_params=[4]), db='hdf52', name='blah')
        # M = pm.MCMC(make_model(this_R, curve_sub = model_salad.nonparametric_model, curve_params=[model_salad.exp_model]), db='hdf52', name=continent+'_extra')
        # M = pm.MCMC(make_model(this_R, curve_sub = model_salad.delta_one_spliner, curve_params=[model_salad.exp_model, 1, scale], 
        #         pr_type=pr_type), db='hdf5', name=dbname)
        # M = pm.MCMC(make_model(this_R, curve_sub = model_salad.exp_plus_line_model, curve_params=[], pr_type=pr_type), db='hdf5', name=dbname)
        M = pm.MCMC(make_model(this_R, curve_sub = model_salad.delta_one_spliner, curve_params=[model_salad.exp_model, 1, scale], 
                pr_type=pr_type), db='hdf5', name=dbname)            
                
            
        nonpr_stochastics = []
        for s in M.stochastics:
            if s.__name__.find('pro_')==-1:
                nonpr_stochastics.append(s)
        M.use_step_method(pm.AdaptiveMetropolis,nonpr_stochastics,delay=200000,interval=1000, shrink_if_necessary=True)
        if pr_type=='unknown':
            M.use_step_method(pm.AdaptiveMetropolis, M.pr, scales={M.pr: .005*np.ones(M.pr.value.shape)}, delay=200000)
        # M.use_step_method(pm.AdaptiveMetropolis,list(M.stochastics - set([M.f_on_mesh])))
        # M.use_step_method(pm.AdaptiveMetropolis, [M.f_on_mesh], scales={M.f_on_mesh: .00000001*np.ones(len(M.mesh))})
        # M.use_step_method(pm.gp.GPParentMetropolis, metro_method=pm.AdaptiveMetropolis(scalar_stochastics))
        # M = pm.MCMC(make_model(this_R, curve_sub = model_salad.exp_model, check_inflec=False), db='hdf52', name=continent)
        break
    except pm.ZeroProbability:
        a,b,c = sys.exc_info()
        print b
        pass
os.chdir('..')
clf()
for i in xrange(10):
    OK = False
    while not OK:
        # M.m_params['asymp'].rand()
        M.m_params['rate'].rand()
        M.f_on_mesh.rand()
        try:
            M.check_trend.logp
            OK = True
        except pm.ZeroProbability:
            pass
    plot(xplot, M.fplot.value*1000)
    axis([0,.8,0,2500])
if i==99:
    raise ValueError, 'Failed to create model'
pl.close('all')
# M.isample(3000000,100000,100)
M.isample(10000000,200000,10000)
# # M.isample(50000,10000,100)
# make_plots(M, continent)

