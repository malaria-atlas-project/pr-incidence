from __future__ import division
from pylab import *
from pymc import *
from numpy import *

x=linspace(0,1,100)

allowed = [lambda x: 90*x/(1+100*x),
            # lambda x: 3*x/(1+x),
            lambda x: 1*x/(1+.1*x),
            lambda x: 3.5*(x - .7*x**2),
            lambda x: 1.5*(flib.invlogit((x-.2)*10))]

close('all')            
for f in allowed:
    plot(x,(f(x)-f(x).min())*1000,'k-')
xlabel(r'Prevalence ($Pf$PR$_{2-10}$)')
ylabel('Expected incidence (per 1000 p.a.)')
savefig('../figs/Figures/allowed.png')
axis([0,1,0,1600])


disallowed = [lambda x: 1.2*x**2,
                lambda x: 90*x/(1+100*x)*(.8+.2*x**2),
                lambda x: (sin(1.7*pi*x)**2 + x**2)*.7]
                
figure()
for f in disallowed:
    plot(x,(f(x)-f(x).min())*1000,'k-')
xlabel(r'Prevalence ($Pf$PR$_{2-10}$)')
ylabel('Expected incidence (per 1000 p.a.)')
savefig('../figs/Figures/disallowed.png')