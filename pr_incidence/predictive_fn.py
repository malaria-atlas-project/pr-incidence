from __future__ import division
import numpy as np
import pymc as pm
from scipy.interpolate import UnivariateSpline, PiecewisePolynomial, interp1d
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
    
    
    def __init__(self, hf_name, pop, nyr=1, burn=0):
        hf = openFile(hf_name)
        cols = hf.root.chain0.PyMCsamples.cols
        
        n = len(cols)
        self.pop = pop
        self.nyr = nyr
        
        self.where_nonzero = np.where(self.pop > 0)
        self.pop_where_pos = self.pop[self.where_nonzero]
        
        self.r_int = cols.r_int[burn:]
        self.r_lin = cols.r_lin[burn:]
        self.r_quad = cols.r_quad[burn:]
        self.f = [interp1d(xplot_aug, np.concatenate(([0],cols.fplot[i])), 'linear') for i in xrange(burn,n)]
        self.n = len(self.r_int)
        
        hf.close()
        
        
    def __call__(self, pr):
        """
        Expects a pr array. Should be of same shape as the pop array that was received as input.
        """
        if pr.shape != self.pop.shape:
            raise ValueError, 'PR input has shape %s, but the population input had shape %s.'%(pr.shape, self.pop.shape)
            
        pr_where_pos = pr[self.where_nonzero]
        out = np.empty(pr.shape)
        
        i = np.random.randint(self.n)
        mu = self.f[i](pr_where_pos)
        r = (self.r_int[i] + self.r_lin[i] * pr_where_pos + self.r_quad[i] * pr_where_pos**2)*self.nyr
        
        rate = pm.rgamma(beta=r/mu, alpha=r) * self.pop_where_pos
        
        burden = pm.rpoisson(rate)
        out[self.where_nonzero] = burden
        
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