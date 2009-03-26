from pylab import *
from numpy import *
import curses
from curses import ascii
import sys
import os
import subprocess

# pr_type, continent, scale = sys.argv[1:]
pr_types = ['mixed']
continents = ['Africa+', 'CSE Asia and Americas']
# continents = ['CSE Asia and Americas']
scales = [.4, .6, 1, 5]
pids = {}
for pr_type in pr_types:
    for continent in continents:
        for scale in scales:
            # os.execv('/usr/bin/screen',['ipython single_analysis_launcher.py %s %s %i'%(pr_type, continent, scale)])
            # os.spawnv(os.P_WAIT,'/usr/local/bin/ipython',['single_analysis_launcher.py', pr_type, continent, scale])
            # q=os.system('ipython single_analysis_launcher.py %s %s %f'%(pr_type, continent, scale))
            this_subprocess = subprocess.Popen(['screen', '/usr/local/bin/ipython','single_analysis_launcher.py', pr_type, continent, str(scale)])
            pids[pr_type + ', ' + continent + ', ' + str(scale)] = this_subprocess
            this_subprocess.communicate(curses.ascii.ctrl('a') + 'd')