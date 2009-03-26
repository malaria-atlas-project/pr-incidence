import numpy as np
from csv import reader
import pylab as pl

# =================
# = Emelda's data =
# =================
R = pl.csv2rec('datafiles/all_data.csv')
# R = pl.csv2rec('all_data.csv')
missing_fields = np.zeros(len(R))
for n in ['surv_int', 'pyor', 'cases', 'region', 'lat', 'lon']:
    missing_fields += R[n].mask
missing_fields = np.array(missing_fields, dtype=bool)
R_mis = R[np.where(missing_fields)]

R=R[np.where(1-missing_fields)].data

R.pr /= 100.
R.mbg_pr /= 100.
R.mix_pr /= 100.
R.pfpr /= 100

R_af = R[np.where(R.region=='Africa+')]
R_am = R[np.where(R.region=='America')]
R_as = R[np.where(R.region=='CSE Asia')]

def time_scaling(pcd, surv_int):
    out = np.ones(len(pcd))
    where_rescale = np.where((pcd!='Y') * (surv_int>7) + (surv_int<7))
    out[where_rescale] = surv_int[where_rescale]/7.
    return out

this_R = R
ar_data = this_R.cases/this_R.pyor
# where_rescale = np.where((R.pcd.data!='Y') * (R.surv_int.data>7) + (R.surv_int.data<7))
# ar_data[where_rescale]*=this_R.surv_int[where_rescale]/7.
ar_data*=time_scaling(R.pcd, R.surv_int)

def filtered_data(pr_type, continent, R=R):
    if pr_type=='data':
        where_good = np.where(1-np.isnan(R.pr))
        this_R = R[where_good]
    elif pr_type=='data_untrans':
        where_good = np.where((1-np.isnan(R.pfpr))*(R.pfpr>0))
        this_R = R[where_good]
    else:
        this_R = R
    
    if not continent=='All':
        where_good = np.where(this_R['region']==continent)        
    
    return this_R[where_good]