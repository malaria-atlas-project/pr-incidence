import numpy as np
from scipy import interpolate as interp
import pylab as pl
import tables as tb
import mbgw
from generic_model import make_plots, xplot
# facs = mbgw.correction_factors.known_age_corr_factors(np.arange(1,6),10000)
PRs_before = np.array([10.7, 21.9, 26.4, 26.3, 25.2])/100
Ns_before = np.array([475, 425, 492, 400, 361])
PRs_after = np.array([3.6, 10.2, 11.2, 13.8, 12.5])/100
Ns_after = np.array([651, 628, 582, 627, 598])

pos_before = np.array(Ns_before * PRs_before, dtype=int)
pos_after = np.array(Ns_after * PRs_after, dtype=int)

P_mesh = np.linspace(.001,1,100)

def splreps_to_post_mean(l):
    all_p = np.sum([interp.splev(P_mesh, ll) for ll in l[1]],axis=0)
    all_p = np.exp(all_p-all_p.max())
    all_p /= (np.sum(all_p))
    return np.sum(all_p*P_mesh)

A = np.rec.fromarrays([np.arange(1,6), np.arange(1,6), pos_before, Ns_before], names=('LOW_AGE', 'UP_AGE', 'PF', 'EXAMINED'))
l = mbgw.correction_factors.age_corr_likelihoods(A, 10000, P_mesh, 'before')

PR_before = splreps_to_post_mean(l)

A = np.rec.fromarrays([np.arange(1,6), np.arange(1,6), pos_after, Ns_after], names=('LOW_AGE', 'UP_AGE', 'PF', 'EXAMINED'))
l = mbgw.correction_factors.age_corr_likelihoods(A, 10000, P_mesh, 'after')

PR_after = splreps_to_post_mean(l)

import os
os.chdir('..')
from all_data import *
os.chdir('generic_analysis')

# pr_type='model_exp'
pr_type = 'mixed'
# pr_type = 'data_untrans'
scale = .6
# pr_type = 'data'
if pr_type=='data':
    where_good = np.where(1-np.isnan(R.pr))
    this_R = R[where_good]
elif pr_type=='data_untrans':
    where_good = np.where((1-np.isnan(R.pfpr))*(R.pfpr>0))
    this_R = R[where_good]

continent = 'Africa+'
where_good = np.where(this_R['region']=='Africa+')

# continent = 'CSE Asia and Americas'
# where_good = where((this_R['region']=='America')+(this_R['region']=='CSE Asia'))

# continent = 'CSE Asia'
# where_good = where(this_R['region']=='CSE Asia')

# continent = 'America'
# where_good = where(this_R['region']=='America')

# continent = 'All'
# where_good=slice(None)

this_R = this_R[where_good]

dbname = continent.replace(' ','_')+'_scale_'+str(scale)+'_'+pr_type
hf = tb.openFile(dbname+'.hdf5')

envs=make_plots(hf.root.chain0.PyMCsamples.cols, dbname, continent, this_R, pr_type)

index_before=np.argmin(np.abs(xplot-PR_before))
index_after=np.argmin(np.abs(xplot-PR_after))

print 'Long-term average:'
print '\t2006: Estimated 2-10yo PR:', PR_before * 100
print '\t\t2.5%,25%,median,75%,97.5%:\n\t\t\t', envs[0][0].lo[index_before], envs[0][1].lo[index_before], envs[0][-1].value[index_before], envs[0][1].hi[index_before], envs[0][0].hi[index_before]
print '\t2008: Estimated 2-10yo PR:', PR_after * 100
print '\t\t2.5%,25%,median,75%,97.5%::\n\t\t\t', envs[0][0].lo[index_after], envs[0][1].lo[index_after], envs[0][-1].value[index_after], envs[0][1].hi[index_after], envs[0][0].hi[index_after]
print 'Any given year:'
print '\t2006: Estimated 2-10yo PR:', PR_before * 100
print '\t\t2.5%,25%,median,75%,97.5%::\n\t\t\t', envs[1][0].lo[index_before], envs[1][1].lo[index_before], envs[1][-1].value[index_before], envs[1][1].hi[index_before], envs[1][0].hi[index_before]
print '\t2008: Estimated 2-10yo PR:', PR_after * 100
print '\t\t2.5%,25%,median,75%,97.5%::\n\t\t\t', envs[1][0].lo[index_after], envs[1][1].lo[index_after], envs[1][-1].value[index_after], envs[1][1].hi[index_after], envs[1][0].hi[index_after]