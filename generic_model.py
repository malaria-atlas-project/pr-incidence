from __future__ import division
import numpy as np
import pymc as pm
import pylab as pl
from IPython.Debugger import Pdb
import model_salad
import os, sys
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from tables import openFile

PR_to_column = {'data_untrans' : 'AB',
'data' : 'AC',
'model_exp' : 'AD',
'mixed' : 'AE'}

def time_scaling(pcd, surv_int):
    out = np.ones(len(pcd))
    where_rescale = np.where((pcd!='Y') * (surv_int>7) + (surv_int<7))
    out[where_rescale] = surv_int[where_rescale]/7.
    return out
    
xplot = pl.linspace(0.001,1,100)

def make_model(recs, curve_sub, curve_params=[], pr_type='mixed', pr_hists = None, pr_samps=None, check_inflec=True):
    input_dict=curve_sub(*curve_params)
    arfun = input_dict['arfun']
    fun_params = input_dict['fun_params']
    
    # if pr_type=='unknown':
    #     splreps = []
    #     for i in xrange(len(pr_hists)):
    #         where_ok = np.where(pr_hists[i][0]>0)
    #         pr_mesh = pr_hists[i][1][where_ok]
    #         lp_mesh = np.log(pr_hists[i][0][where_ok])
    #         splreps.append(UnivariateSpline(pr_mesh, lp_mesh, bbox=[0,1]))
    # 
    #     @pm.stochastic(dtype=float)
    #     def pr(value = pr_hists[:,1,10], splreps = splreps):
    #         out=0
    #         for i in xrange(len(value)):
    #             this_value = value[i]
    #             if this_value<0 or this_value>1:
    #                 return -np.inf
    #             else:
    #                 out += splreps[i](this_value)
    #         return out
    if pr_type=='model_exp': 
        pr = recs.mbg_pr
    elif pr_type=='data':
        pr=recs.pr
    elif pr_type=='mixed':
        pr=recs.mix_pr
    elif pr_type=='data_untrans':
        pr=recs.pfpr
    else:
        raise ValueError, 'PR type unknown'

    # # A deterministic that measures the change in attack rate given a certain change in PR.
    # delta_ar = pm.Lambda('delta_ar', lambda fp = fun_params: np.diff(arfun(diff_pts, *fp)))
    fboth = pm.Lambda('fboth', lambda fp=fun_params, pr=pr: arfun(np.hstack((pr, xplot)), *fp))
    
    # Evaluation of trend at PR values
    AR_trend = pm.Lambda('AR_trend', lambda fp = fun_params, pr=pr: arfun(pr, *fp))

    # The function evaluated on the display mesh
    fplot = pm.Lambda('fplot', lambda fp = fun_params: arfun(xplot, *fp))
    pl.clf()
    pl.plot(xplot,fplot.value)
    
    @pm.potential
    def check_trend(AR = AR_trend, f=fplot):
        if np.any(AR<=0) or np.any(f<=0):
            return -np.Inf
        if check_inflec:
            d2 = np.diff(f,2)
            d2=d2[np.where(np.abs(d2)>1e-6)]
            chgs = np.where(np.abs(np.diff(np.sign(d2)))>1)[0]
            if np.diff(f[-3:],2) >= 0 or len(chgs) > 1:
                return -np.Inf
        return 0

    # Negative-binomial parameters.
    r_int = pm.Exponential('r_int',.0001,value=.3)
    r_lin = pm.Uninformative('r_lin',value=1.)
    r_quad = pm.Uninformative('r_quad',value=.1)

    rplot = pm.Lambda('rplot', lambda r_int=r_int, r_lin=r_lin, r_quad=r_quad: r_int + r_lin*xplot + r_quad*xplot**2)

            
    @pm.potential
    def check_r(i=r_int, l=r_lin, q=r_quad):
        # if q>0:
        #     xhat = -l / 2 / q
        #     if i + l*xhat + q*xhat*xhat <= 0 and xhat>0:
        #         return -np.Inf
        if l <= 0 or l + 2.*q <= 0:
            return -np.Inf
        if i + l + q <= 0 or i < 0:
            return -np.Inf
        return 0

    # shape parameter of gamma process is multiplied by total survey time
    time_scale_fac = time_scaling(recs.pcd, recs.surv_int)
    tottime = (recs.yr_end - recs.yr_start + 1)
    scale_time = tottime / time_scale_fac
    pop = recs.pyor/tottime

    # Shape parameter of Poisson intensity is only multiplied by scaled survey time.
    r = pm.Lambda('r', lambda i=r_int, l=r_lin, q=r_quad, pr=pr: (i + l * pr + q * pr * pr)*scale_time)
    
    # scale parameter of Poisson intensity is multiplied by scaled survey time * number of people sampled.
    exp_rate = pm.Lambda('exp_rate', lambda t=AR_trend: scale_time*pop*t)
    
    # The data
    AR = pm.NegativeBinomial('AR', exp_rate, r, value=recs.cases, observed=True)
    
    @pm.deterministic(dtype=float)
    def AR_dev(AR=AR, mu=exp_rate, r=r):
        return np.array([pm.negative_binomial_like(AR[i], mu[i], r[i]) for i in xrange(len(AR))])

    out = locals()
    out.update(input_dict)
    return out

def make_plots_from_model(M, continent):
    recs = M.recs
    dbname = M.db.filename
    cols = M.db._h5file.root.chain0.PyMCsamples.cols
    pr_type = M.pr_type
    make_plots(cols, dbname, continent, recs, pr_type)

def make_plots(cols, dbname, continent, recs, pr_type, nyr = 1):
    
    samp_size=1000

    if continent.find('Africa') >= 0:
        lims = [.8,2.5]                
    elif continent.find('Asia')>=0:
        lims = [.5,1.5]
    elif continent.find('America')>=0:
        lims = [.2,1.]
    else:
        lims = [.8,2.5]

    model_id = dbname + '_' + PR_to_column[pr_type]
    time_scale_fac = time_scaling(recs.pcd, recs.surv_int)
    
    if pr_type=='model_exp': 
        pr = recs.mbg_pr
        
    elif pr_type=='data':
        pr=recs.pr
    elif pr_type=='mixed':
        pr=recs.mix_pr
    elif pr_type=='data_untrans':
        pr=recs.pfpr
    else:
        raise ValueError, 'PR type unknown'
    
    pl.clf()
    envs_post = pm.Matplot.func_envelopes(cols.fplot[:]*samp_size, [.25, .5, .9])
    for env in envs_post:
        env.display(xplot, .8, new=False) 
    pl.xlabel(r'Prevalence $(Pf$PR$_{2-10})$')
    pl.ylabel('Incidence (per 1000 p.a.)')
    # ar_data = recs.cases/recs.pyor/np.minimum(1,7./recs.surv_int)
    ar_data = recs.cases/recs.pyor*time_scale_fac
    # print ar_data.min()*samp_size, ar_data.max()*samp_size
    ar_in =  ar_data[np.where((pr<.25) * (pr > .10))]*samp_size
    print ar_in.min(), ar_in.max()
    pl.plot(pr, ar_data*samp_size, 'r.', label='data')
    # pl.title(continent)
    # pl.legend(loc=2)
    pl.axis([0,lims[0],0,2500])
    
    pl.savefig('../figs/Figures/%s_post.png'%model_id)
    
    pl.figure()
    Nsamps = len(cols.r)
    AR_pred = np.empty((Nsamps*100, len(xplot)))
    for i in xrange(Nsamps):
        this_mu = cols.fplot[i]
        this_r = (cols.r_int[i] + cols.r_lin[i] * xplot + cols.r_quad[i] * xplot**2)
        
        # Uncomment to make an actual sample
        # AR_pred[i*100:(i+1)*100,:] = pm.rnegative_binomial(r=this_r, mu=this_mu*samp_size, size=100)
        # Uncomment to multiply population-wide AR
        AR_pred[i*100:(i+1)*100,:] = pm.rgamma(beta=this_r*nyr/this_mu, alpha=this_r*nyr, size=100)*1000
    envs_pred = pm.Matplot.func_envelopes(AR_pred, [.25, .5, .9])
    for env in envs_pred:
        env.display(xplot, .8, new=False)
    thirty_index = np.argmin(np.abs(xplot-.3))
    print envs_pred[0].hi[thirty_index], envs_pred[0].lo[thirty_index]
    
    pl.xlabel(r'Prevalence $(Pf$PR$_{2-10})$')
    pl.ylabel('Incidence (per 1000 p.a.)')
    pl.plot(pr, ar_data*samp_size, 'r.', label='data')
    # pl.title(continent)
    # pl.legend(loc=2)
    
    pl.axis([0,lims[0],0,2500])
    
    pl.savefig('../figs/Figures/%s_pred.png'%model_id)
    
    # Pdb(color_scheme='Linux').set_trace()
    # if hasattr(recs.lat, 'mask'):
    # where_lonlat = np.where(1-recs.lat.mask)
    # # else:
    # #     where_lonlat = np.where(1-np.isnan(recs.lat))
    # lat = recs.lat[where_lonlat]
    # lon = recs.lon[where_lonlat]
    mean_dev = np.mean(cols.AR_dev[:], axis=0)#[where_lonlat]
    devs = np.rec.fromarrays([mean_dev, recs.lon, recs.lat], names=('mean_deviance','longitude','latitude'))
    pl.rec2csv(devs, '../figs/%s_deviance.csv'%model_id)
    # pl.close('all')
    return envs_post, envs_pred
    