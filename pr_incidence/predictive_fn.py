from __future__ import division
import numpy as np
import pymc as pm
from scipy.interpolate import interp1d
from tables import openFile


xplot = np.linspace(0.001,1,100)
xplot_aug = np.concatenate(([0],xplot))


class BurdenPredictor(object):
    """
    Generates a callable object that converts PR to burden.
      - cols : cols attribute of a PyTables table containing the MCMC trace.
      - pop : A population surface, represented as a vector.
      - nyr : Integer, number of years to predict for.
    """
    
    
    def __init__(self, hf_name, nyr=1, burn=0):
        hf = openFile(hf_name)
        cols = hf.root.chain0.PyMCsamples.cols
        
        n = len(cols)
        self.nyr = nyr
        
        
        self.r_int = cols.r_int[burn::10]
        self.r_lin = cols.r_lin[burn::10]
        self.r_quad = cols.r_quad[burn::10]
        self.f = [interp1d(xplot_aug, np.concatenate(([0],cols.fplot[i])), 'linear') for i in xrange(burn,n,10)]
        self.n = len(self.r_int)
        
        hf.close()
        
        
    def __call__(self, pr, pop, pop_pr_res):
        """
        Expects a pr array. Should be of same shape as the pop array that was received as input.
        """
        if pr.shape != pop[::pop_pr_res,::pop_pr_res].shape:
            raise ValueError, 'PR input has shape %s, but the population input had shape %s.'%(pr.shape, pop.shape)

        out = np.zeros((pr.shape[0]*pop_pr_res, pr.shape[1]*pop_pr_res))        
        where_pos = np.where(pr > 0)
        if len(where_pos[0])==0:
            return out
        pr_where_pos = pr[where_pos]

        
        i = np.random.randint(self.n)
        mu = self.f[i](pr_where_pos)
        r = (self.r_int[i] + self.r_lin[i] * pr_where_pos + self.r_quad[i] * pr_where_pos**2)*self.nyr
        
        rate = pm.rgamma(beta=r/mu, alpha=r) * pop[where_pos]
        
        for j,k in zip(where_pos):
            out[k*pop_pr_res:(k+1)*pop_pr_res, j*pop_pr_res:(j+1)*pop_pr_res] \
                = np.random.poisson(rate[j,k],size=pop_pr_res*pop_pr_res).reshape((pop_pr_res,pop_pr_res))
        
        return out
        
        
if __name__ == '__main__':
    from tables import openFile
    from pylab import *

    N=10000
    pop=10000*np.ones(N)
    nyr = 10
    
    pop[::10] = 0

    p = BurdenPredictor('traces/Africa+_scale_0.6_model_exp.hdf5', pop, nyr, 0)
    pr_max = .6
    # p = BurdenPredictor('traces/CSE_Asia_and_Americas_scale_0.6_model_exp.hdf5', np.ones(N)*pop, nyr)
    # pr_max = .5
    
    # for i in xrange(10000):

    pr = pm.runiform(0,pr_max,size=N)
    
    # for i in xrange(p.n):
    #     clf()
    #     plot(xplot, p.cols.fplot[i],'r.',markersize=2)
    #     plot(pr, p.f[i](pr), 'b.',markersize=1)
    #     a=raw_input()
    
    b = p(pr)
    
    clf()
    plot(pr,b,'k.',markersize=1)