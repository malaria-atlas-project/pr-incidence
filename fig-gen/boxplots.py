import os
from pylab import *
from numpy import *
import matplotlib
matplotlib.interactive(False)

os.chdir('..')
from all_data import R, ar_data

os.chdir('fig_gen')

# where_good=where(1-isnan(R.pr))
# 
# xax = R.pr[where_good]
# yax = ar_data[where_good]

xax_us = R.mix_pr
yax_us = ar_data * 1000

def nice_boxplot(xbx, *args, **kwargs):
    lines = boxplot(xbx, *args, **kwargs)

    for b in lines['boxes']:
        b.set_color('k')
        b.set_linewidth(1)
        fill(b.get_path().vertices[:,0], b.get_path().vertices[:,1], facecolor='.8')
    for w in lines['whiskers']:
        w.set_color('k')
        w.set_linestyle('-')
    for c in lines['caps']:
        c.set_color('k')
        c.set_linestyle('-')
    for m in lines['medians']:
        m.set_color('k')
        m.set_linewidth(3)
        m.set_solid_capstyle('round')
    for f in lines['fliers']:
        f.set_markerfacecolor('1')
        

def box(xax, yax, adj=True):
    f=figure(figsize=(9,9))
    a = f.gca()
    ybx = [yax[where(xax<=.1)], yax[where((xax>.1) * (xax<=.5))], yax[where(xax>.5)]]
    nice_boxplot(ybx, sym='o')
    labs = ['Hypo','Meso','Hyper-holo']
    labs = [labs[i] + '\nn=%i'%len(ybx[i]) for i in xrange(3)]
    a.set_xticklabels(labs)
    xlab= r'Prevalence ($Pf$PR)'
    if adj:
        xlab=xlab+r'$_{2-10}$'
    a.set_xlabel(xlab)
    a.set_yticks(arange(5)*500)
    a.set_ylabel('Incidence (per 1000 p.a.)')


box(xax_us, yax_us)
savefig('../figs/hay_etal_box.png')

old_data = csv2rec('snow_etal.csv')
old_data = old_data[where((1-isnan(old_data.pr))*(old_data.pr>0))]

box(old_data.pr/100., old_data.pfi, adj=False)
savefig('../figs/snow_etal_box.png')

os.system('open ../figs/*box.png')