from generic_model import make_model, model_salad, xplot
from pymc import *
from pylab import *
from numpy import *
from pylab import *
from numpy import *
from all_data import *

pr_type='model_exp'
# pr_type = 'mixed'
# pr_type = 'data_untrans'
scale = .6
# pr_type = 'data'
if pr_type=='data':
    where_good = where(1-np.isnan(R.pr))
    this_R = R[where_good]
elif pr_type=='data_untrans':
    where_good = where((1-np.isnan(R.pfpr))*(R.pfpr>0))
    this_R = R[where_good]

continent = 'Africa+'
where_good = where(this_R['region']=='Africa+')

# continent = 'CSE Asia and Americas'
# where_good = where((this_R['region']=='America')+(this_R['region']=='CSE Asia'))

# continent = 'CSE Asia'
# where_good = where(this_R['region']=='CSE Asia')

# continent = 'America'
# where_good = where(this_R['region']=='America')

# continent = 'All'
# where_good=slice(None)

this_R = this_R[where_good]

r_int_lims = [0,2.]
r_quad_lims = [0,.2]

for i in xrange(100):
    try:
        M = Model(make_model(this_R, curve_sub = model_salad.delta_one_spliner, curve_params=[model_salad.exp_model, 1, scale], pr_type=pr_type))
        break
    except:
        pass


        
def draw():
    while True:
        
        M.r_int.value = runiform(*r_int_lims)
        M.r_quad.value = runiform(*r_quad_lims)
        for g in M.generations:
            for s in g:
                if s._random is not None:
                    s.rand()
        try:
            for s in M.stochastics | M.observed_stochastics | M.potentials:
                s.logp
            break
        except:
            pass
    
    
    plot(xplot, M.fplot.value*1000)
    
clf()
for i in xrange(10):
    draw()
    
pl.xlabel(r'Expected prevalence $(Pf$PR$_{2-10})$')
pl.ylabel('Incidence (per 1000 p.a.)')

savefig('../figs/Figures/prior_draws.png')