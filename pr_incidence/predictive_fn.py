from __future__ import division
import numpy as np
import pymc as pm
from scipy.interpolate import interp1d
from tables import openFile


xplot = np.linspace(0.001,1,100)
xplot_aug = np.concatenate(([0],xplot))

class MeanIncidencePredictor(object):
    """
    Generates a callable object that converts PR to burden.
      - cols : cols attribute of a PyTables table containing the MCMC trace.
      - pop : A population surface, represented as a vector.
      - nyr : Integer, number of years to predict for.
    """
    
    
    def __init__(self, hf_name, burn=0):
        hf = openFile(hf_name)
        cols = hf.root.chain0.PyMCsamples.cols
        
        n = len(cols)
                
        self.fs = np.array([np.concatenate(([0],cols.fplot[i])) for i in xrange(burn,n,10)])
        self.f = interp1d(xplot_aug, np.mean(self.fs, axis=0), 'linear')
        
        hf.close()
        
    def __call__(self, pr):
        """
        Expects a pr array. Should be of same shape as the pop array that was received as input.
        """
        out = pr.copy()
        out.data[np.where(True-pr.mask)] = self.f(pr.data[np.where(True-pr.mask)])
        return out

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
        
        
    def pr5km_pop1km(self, pr, pop, pop_pr_res):
        """
        Expects a pr array in 5km and pop array in 1km
        """

        #if len(pr.shape)>0:
        if len(pr.shape)>1:
            raise ValueError, 'PR is supposed to be 1d, dumbass.'

        #if pr.shape != (0,pop[::pop_pr_res].shape[1]):
        if pr.shape[0]!=np.shape(pop[:,::pop_pr_res])[1]:
            raise ValueError, 'PR input has shape %s, but the population input had shape %s.'%(pr.shape, pop.shape)

        # define blank 1km 2-d array to house burden
        burden_1km = np.zeros((pop_pr_res,pr.shape[0]*pop_pr_res))

        # extract vector of pr at 5km only wherre pr is non-zero - if all zero then return blank template        
        where_pr_pos_5km = np.where(pr > 0)
        if len(where_pr_pos_5km[0])==0:
            return burden_1km
        pr_where_pr_pos_5km = np.atleast_1d(pr[where_pr_pos_5km])
        
        # initialise 5km zero 1-d array for rate
        rate_5km = np.zeros(np.product(np.shape(pr))).reshape(np.shape(pr))
        
        # calculate rate for non-zero PR pixels
        i = np.random.randint(self.n)
        mu = self.f[i](pr_where_pr_pos_5km)
        r = (self.r_int[i] + self.r_lin[i] * pr_where_pr_pos_5km + self.r_quad[i] * pr_where_pr_pos_5km**2)*self.nyr
        rate_where_pr_pos_5km = np.atleast_1d(pm.rgamma(beta=r/mu, alpha=r))
        
        # re-map thse rate values onto full length 5km rate vector
        rate_5km[where_pr_pos_5km]=rate_where_pr_pos_5km
        
        # blow-up 5km rate vector to match 1km pop 2-d array
        temp1 = np.atleast_2d(np.arange(0,len(pr))) 
        temp2 = np.repeat(temp1,pop_pr_res,axis=0)
        temp3 = np.repeat(temp2,pop_pr_res,axis=1)
        rate_1km=rate_5km[temp3]

        if(np.shape(pop)!=np.shape(rate_1km)):
            raise ValueError, '1km rate array has shape %s, but the 1km population array has shape %s.'%(np.shape(rate_1km),np.shape(pop))
        
        # multiply 1km rate by 1km pop array 
        popRate = rate_1km*pop
        
        # extract non-zero pixels (now also excludes zero Pop as well as zero rate), and return all zeroes if no non-zero pixels
        where_popRate_pos = np.where(popRate > 0)
        if len(where_popRate_pos[0])==0:
            return burden_1km 
        popRate_where_popRate_pos = popRate[where_popRate_pos]
                    
        # carry out poisson draws to define burden in these non-zero pixels
        burden_where_popRate_pos = np.random.poisson(popRate_where_popRate_pos)

        # re-map burden values to full 1km 2-d burden array
        burden_1km[where_popRate_pos] = burden_where_popRate_pos
        
        #for l in xrange(0,len(where_pos[0])):
        #    j=where_pos[0][l]
        #    out[:,j*pop_pr_res:(j+1)*pop_pr_res] = np.random.poisson(rate[l]*pop[:,j*pop_pr_res:(j+1)*pop_pr_res],size=(pop_pr_res,pop_pr_res))
        
        return burden_1km

    def pr5km_pop5km(self, pr, pop):
        """
        Expects a pr array in 5km and pop array in 5km
        """

        #if pr.shape != (0,pop[::pop_pr_res].shape[1]):
        if pr.shape!=pop.shape:
            raise ValueError, 'PR input has shape %s, but the population input had shape %s.'%(pr.shape, pop.shape)

        # define blank 5km 2-d array to house burden
        burden_5km = np.zeros(np.product(pr.shape)).reshape(pr.shape)

        # extract vector of pr at 5km only where pr is non-zero - if all zero then return blank template        
        where_pr_pos_5km = np.where(pr > 0)
        if len(where_pr_pos_5km[0])==0:
            return burden_5km
        pr_where_pr_pos_5km = np.atleast_1d(pr[where_pr_pos_5km])
       
        # initialise 5km zero 1-d array for rate
        rate_5km = np.zeros(np.product(pr.shape)).reshape(pr.shape)
        
        # calculate rate for non-zero PR pixels
        i = np.random.randint(self.n)
        mu = self.f[i](pr_where_pr_pos_5km)
        r = (self.r_int[i] + self.r_lin[i] * pr_where_pr_pos_5km + self.r_quad[i] * pr_where_pr_pos_5km**2)*self.nyr
        rate_where_pr_pos_5km = np.atleast_1d(pm.rgamma(beta=r/mu, alpha=r))
        
        # re-map thse rate values onto full length 5km rate vector
        rate_5km[where_pr_pos_5km]=rate_where_pr_pos_5km
        
        if(np.shape(pop)!=np.shape(rate_5km)):
            raise ValueError, '1km rate array has shape %s, but the 1km population array has shape %s.'%(np.shape(rate_5km),np.shape(pop))
        
        # multiply 5km rate by 5km pop array 
        popRate = rate_5km*pop
        
        # extract non-zero pixels (now also excludes zero Pop as well as zero rate), and return all zeroes if no non-zero pixels
        where_popRate_pos = np.where(popRate > 0)
        if len(where_popRate_pos[0])==0:
            return burden_5km 
        popRate_where_popRate_pos = popRate[where_popRate_pos]
                    
        # carry out poisson draws to define burden in these non-zero pixels
        burden_where_popRate_pos = np.random.poisson(popRate_where_popRate_pos)

        # re-map burden values to full 1km 2-d burden array
        burden_5km[where_popRate_pos] = burden_where_popRate_pos
        
        #for l in xrange(0,len(where_pos[0])):
        #    j=where_pos[0][l]
        #    out[:,j*pop_pr_res:(j+1)*pop_pr_res] = np.random.poisson(rate[l]*pop[:,j*pop_pr_res:(j+1)*pop_pr_res],size=(pop_pr_res,pop_pr_res))
        
        return burden_5km
        
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