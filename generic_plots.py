import model_salad as ms
from generic_model import *
from csv import reader
import numpy as np
import pymc as pm
from pylab import *
import os
os.chdir('..')
from all_data import *
os.chdir('generic_analysis')
from IPython.kernel import client
import pickle

all_continents=['Asia','Africa','All']
# conts = dict(zip(all_continents, range(len(all_continents))))
models = [('powerlaw_model', []),
            ('polynomial_model', range(1,7)),
            ('smooth_cp_model', []),
            ('cp_model', []),
            ('bh_model', []),
            ('exp_model', [])]
models_expanded = []
for model_pair in models:
    if len(model_pair[1]) > 0:
        for pval in model_pair[1]:
            models_expanded.append((model_pair[0], pval))
    else:
        models_expanded.append((model_pair[0], None))

mec = client.MultiEngineClient()
mec.reset()
# mec.block=True
mec.block=True
n_engines = len(mec.get_ids())
mec_ids = mec.get_ids()

job_queue = [[] for i in xrange(n_engines)]
ng = 0
for mod in models_expanded:
    for cont in all_continents:
        cont_clean = cont.replace(' ','_').replace('+','')
        # print (cont_clean,)+mod, ng%n_engines
        job_queue[mec_ids[ng%n_engines]].append((cont_clean,)+mod)
        ng += 1

mec.execute("import os")
# mec.execute("os.chdir('/Users/anand/renearch/malaria/Burden/code/generic_analysis')")
mec.execute("import matplotlib")
mec.execute("matplotlib.use('PDF')")
# mec.execute("matplotlib.interactive(False)")
mec.execute("import pymc as pm")
mec.execute("from generic_model import *")
mec.execute("import model_salad as ms")
mec.push({'R':R})
mec.execute('dics={}')
mcmcs_done = []
plots_done = []
dics_done = []
dics = {}
for task in ['setup', 'run']:
    if task=='setup':
        mec.block=True
    elif task=='run':
        mec.block=False
    for i in mec_ids:

        for continent, model, param in job_queue[i]:
            this_job_name = continent + '_' + model
            param_str = '()'
            params = ()
            if param is not None:
                this_job_name += '_'+str(param)
                param_str = '(%s,)'%str(param)
                params = (param,)
            if task=='setup':
                # mec.push({'continent':continent},i)            
                if continent == 'All':
                    mec.execute("this_R=R",i)
                    this_R = R
                elif continent == 'Africa':
                    mec.execute("this_R = R[np.where(R['region']=='Africa')]", i)
                    this_R = R[np.where(R['region']=='Africa')]
                else:
                    mec.execute("this_R = R[np.where((R['region']=='Americas')+(R['region']=='Asia'))]", i)
                    this_R = R[np.where((R['region']=='Americas')+(R['region']=='Asia'))]

                print 'Setting up %s on %i' %(this_job_name, i)
        
                for j in xrange(100):
                    # print j
                    try:
                        try:
                            os.remove('%s.hdf5'%this_job_name)
                        except:
                            pass
                        mec.execute("%s = pm.MCMC(make_model(this_R, curve_sub=ms.%s, curve_params=%s), name='%s', db='hdf52', complevel=9)"%(this_job_name, model, param_str, this_job_name),i)                
                        # M = pm.MCMC(make_model(this_R, curve_sub=getattr(ms,model), curve_params=params), name=this_job_name, db='hdf52', complevel=9)
                        break
                    except:
                        pass
                if j==99:
                    raise ValueError, 'Failed to create model'
            elif task=='run':
                print 'Running %s on %i' %(this_job_name, i)                
                # mcmcs_done.append(mec.execute("%s.sample(150, 20, 1, tune_interval=10)"%this_job_name,i))                
                mcmcs_done.append(mec.execute("%s.sample(500000, 200000, 100, tune_interval=1000)"%this_job_name,i))        
                plots_done.append(mec.execute("make_plots(%s,'%s')"%(this_job_name,continent),i))
                dics_done.append(mec.execute("dics['%s'] = %s.dic()"%(this_job_name,this_job_name),i))
                # M.sample(150, 20, 1, tune_interval=10)
                # make_plots(M, this_R, continent)

# Block until all jobs are done
for res in mcmcs_done + plots_done + dics_done:
    res.get_result()

mec.block=True        
[dics.update(mec.pull('dics',i)[0]) for i in mec_ids]
dic_file = file('dics.txt','w')
for cont in ['Africa', 'Asia', 'All']:
    dic_file.write('%s:\n'%cont.replace('_',' ').replace('Africa','Africa+').replace('Asia', 'CSE Asia and Americas'))
    these_dics = []
    for dcont, dic in dics.iteritems():
        if dcont.find(cont)==0:
            these_dics.append((dcont,dic))
    order = np.argsort([tup[1] for tup in these_dics])
    for item in order:
        dcont, dic = these_dics[item][0], these_dics[item][1]
        dic_file.write('\t-%s: DIC %s\n' % (dcont.replace(cont,'').replace('_',' '), dic))
dic_file.write("""Key:
	- cp = 'changepoint'
	- bh = 'Beverton-Holt', y=ax/(1+bx)
	- exp = 'exponential', y=a (1-exp(-bx))
	- smooth cp = 'smoothed changepoint'
	- powerlaw : y=ax**b""")
dic_file.close()

# import pickle
# pickle.dumps(dics,file('dics.pickle','w'))
# Remember to pickle DIC's !