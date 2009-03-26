# AR low up PR
from pylab import *
from numpy import *
from all_data import *

pr_type='model_exp'
# pr_type = 'mixed'
# pr_type = 'data_untrans'
scale = .6
# pr_type = 'data'

continent = 'Africa+'
# continent = 'CSE Asia and Americas'
# continent = 'CSE Asia'
# continent = 'America'
# continent = 'All'

this_R = filtered_data(pr_type, continent)

dbname = continent.replace(' ','_')+'_scale_'+str(scale)+'_'+pr_type
hf = openFile('traces/'+dbname+'.hdf5')

envs=make_plots(hf.root.chain0.PyMCsamples.cols, dbname, continent, this_R, pr_type, nyr=2)    

