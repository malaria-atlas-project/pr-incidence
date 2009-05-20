from __future__ import division
import numpy as np
import pymc as pm
import pylab as pl
from IPython.Debugger import Pdb   

# DONE: CSE_Asia_mixed, CSE_Asia_data, CSE_Asia_and_Americas_data, mixed, Africa+_data

def powerlaw_model():
    """
    Power-law model
    """
    power = pm.Uniform('power', 0, 10, 1)
    coef = pm.Uniform('coef',0,100,1)    
    fun_params = [power, coef]
    def arfun(x, power, coef):
        return coef*x**power
    return locals()
    
def polynomial_model(order):
    """
    Polynomial model of arbitrary order
    """
    coefs = []
    for i in xrange(order):
        coefs.append(pm.Uninformative('coef_%i'%i, 0))
    coefs[0].value = 1
    fun_params = coefs
    def arfun(x, *coefs):
        return np.sum([coefs[i]*x**(i+1) for i in xrange(len(coefs))], axis=0)
    return locals()

def smooth_cp_model():
    # # s1: Slope of segment 1.
    # s1 = pm.Uninformative('s1', 1)
    # 
    # 
    # # cp: Location of changepoint of decrease.
    # # Must be positive.
    # # Use informative prior to prevent it running off.
    # cp = pm.Uniform('cp',0,.4)
    # 
    # # s2: Slope of segment 2.   
    # s2_lb = pm.Lambda('s2_lb', lambda s1=s1,cp=cp: max(-6000,-s1*cp/(1-cp)))
    # s2 = pm.Uniform('s2', s2_lb, 6000)
    # 
    # # Smoothing factor for changepoint
    # smoove_facta = pm.Uniform('smoove_facta',0,10)    
    # fun_params = [s1, s2, cp, smoove_facta]
    # 
    # # The piecewise linear function
    # def arfun(x, s1, s2, cp, sf):
    #     o1 = x*s1
    #     o2 = (x-cp)*s2 + cp*s1
    #     w = 1./(x+sf)
    #     out = w*o1 + (1-w)*o2
    #     return out
    #     
    # # The function evaluated on the display mesh
    # return locals()

    # s1: Slope of segment 1.
    s1 = pm.Uninformative('s1', 1)
    
    # cp: Location of changepoint of decrease.
    # Must be positive.
    # Use informative prior to prevent it running off.
    cp = pm.Uniform('cp',.01,.4)
    
    s2 = pm.Uniform('s2', -6000, 6000,.5)
    sf = pm.Uniform('smoove_facta',0,1,.5)
    
    fun_params = [s1, s2, cp,sf]
    
    # The piecewise linear function
    def arfun(x, s1, s2, cp,sf):
        o1 = x*s1
        o2 = (x-cp)*s2 + cp*s1
        w = pm.invlogit((np.log(x)-np.log(cp))/sf)
        return w*o2  + (1-w)*o1
    return locals()
        
def cp_model():
    """
    Changepoint model
    """
    # s1: Slope of segment 1.
    s1 = pm.Uninformative('s1', 1)

    # cp: Location of changepoint of decrease.
    # Must be positive.
    # Use informative prior to prevent it running off.
    cp = pm.Uniform('cp',0,.4)

    s2 = pm.Uniform('s2', -6000, 6000, .5)
    
    fun_params = [s1, s2, cp]

    # The piecewise linear function
    def arfun(x, s1, s2, cp):
        out = np.empty(x.shape)
        where_less = np.where(x<cp)
        where_greater = np.where(x>=cp)
        out[where_less] = x[where_less] * s1
        out[where_greater] = (x[where_greater]-cp)*s2 + cp*s1
        return out
    return locals()
    
def bh_model():
    """
    Beverton-hold model
    """
    coef = pm.Uniform('coef',0,100,1)
    rate = pm.Uniform('rate',0,100,1)
    fun_params = [coef, rate]
    def arfun(x,coef, rate):
        return coef*x/(1+x/rate)
    return locals()

def exp_plus_line_model():
    """
    Exponential model
    """
    asymp = pm.Uniform('asymp',0,1,value=.2)
    rate = pm.Uniform('rate',.1,100,value=44)
    extra_slope = pm.Uniform('extra_slope',-1000,1000,value=1)
    fun_params = [asymp, rate, extra_slope]

    def arfun(x, asymp, rate, extra_slope):
        return (1 - np.exp(-x*rate))*asymp + extra_slope*x
    return locals()        

    
def exp_model():
    """
    Exponential model
    """
    asymp = pm.Uniform('asymp',0,100,value=.2)
    rate = pm.Uniform('rate',.1,100,value=44)
    fun_params = [asymp, rate]
    
    def arfun(x, asymp, rate):
        return (1 - np.exp(-x*rate))*asymp
    return locals()        
    
def nonparametric_model(mean_model=None):

    amp = pm.Uniform('amp',0,100,3)#,observed=True)
    scale = pm.Uniform('scale',0,100,.6)#,observed=True)
    if mean_model is None:
        m_fun = lambda x: np.zeros(len(x))
        m_params = {}
    else:
        m_dict=mean_model()
        m_fun=m_dict['arfun']
        m_params = dict([(fp.__name__, fp) for fp in m_dict['fun_params']])

    @pm.deterministic(trace=False)
    def M_and_C(amp=amp, scale=scale, mp=m_params):
        M = pm.gp.Mean(m_fun, **mp)
        C = pm.gp.NearlyFullRankCovariance(pm.gp.cov_funs.gaussian.euclidean, amp=amp, scale=scale)
        pm.gp.observe(M,C,[0],[0])
        return M,C
    M = pm.Lambda('M', lambda both=M_and_C: both[0], trace=False)
    C = pm.Lambda('C', lambda both=M_and_C: both[1], trace=False)
    f = pm.gp.GP('f',M,C,mesh=np.linspace(.1,.7,10),trace=False)

    fun_params = [f]

    def arfun(x,f):
        return f(x)
    return locals()
    
def delta_one_spliner(mean_model=None, amp=1, scale=.6):
    import scipy
    from scipy import interpolate
    
    print 'NOTICE here you have scale bounded above because the deviation just goes linear otherwise.'
    
    # amp = pm.Uniform('amp',0,100,amp,observed=True)
    # scale = pm.Uniform('scale',.1,2,scale,observed=True)
    if mean_model is None:
        m_fun = lambda x: np.zeros(len(x))
        m_params = {}
    else:
        m_dict=mean_model()
        m_fun=m_dict['arfun']
        m_params = dict([(fp.__name__, fp) for fp in m_dict['fun_params']])

    @pm.deterministic(trace=False)
    def M_and_C(amp=amp, scale=scale, mp=m_params):
        if amp<=0 or scale<=0:
            raise pm.ZeroProbability
        M = pm.gp.Mean(m_fun, **mp)
        C = pm.gp.NearlyFullRankCovariance(pm.gp.cov_funs.gaussian.euclidean, amp=amp, scale=scale)
        pm.gp.observe(M,C,[0],[0])
        return M,C

    mesh = np.linspace(.1,1,5)
    M = pm.Lambda('M', lambda both=M_and_C: both[0], trace=False)
    C = pm.Lambda('C', lambda both=M_and_C: both[1](mesh,mesh), trace=False)

    f_on_mesh = pm.MvNormalCov('f_on_mesh',0.*mesh,C)
    
    @pm.deterministic
    def f(f_on_mesh=f_on_mesh):
        return scipy.interpolate.interp1d(np.hstack(([0], mesh)),np.hstack(([0], f_on_mesh)), kind='cubic')

    fun_params = [f,M]

    def arfun(x,f,M):
        return f(x)+M(x)
    return locals()
