from __future__ import division
from pylab import *
from pymc import *
from numpy import *

T = 50
t = arange(T)+.5
shape_lo = .1
scale_lo = 50

shape_hi = 3
scale_hi = 200 

ts_lo = random.gamma(shape_lo, scale_lo/shape_lo, size=T)
ts_hi = random.gamma(shape_hi, scale_hi/shape_hi, size=T)

close('all')
# plot(t,ts_lo,'k.')
bar(t-.5,ts_lo,1.,linewidth=0,edgecolor='w',color='.8')
plot([0,T],[scale_lo,scale_lo],'k-.',label='Expected attack rate')
axis([0, T, 0, 1600])
xlabel('Year')
ylabel('Incidence (per 1000 p.a.)')
# legend(loc=0)
savefig('../figs/Figures/ts_low.png')


figure()
bar(t-.5,ts_hi,1.,linewidth=0,edgecolor='w',color='.8')
plot([0,T],[scale_hi,scale_hi],'k-.',label='Expected attack rate')
axis([0, T, 0, 1600]) 
xlabel('Year')
ylabel('Incidence (per 1000 p.a.)')
# legend(loc=0)
savefig('../figs/Figures/ts_hi.png')