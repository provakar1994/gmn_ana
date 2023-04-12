# This script will jcashe the 0th segment of a set of runs
# from the run DB.

import os
import pandas as pd

## ************* ##
##  User inputs  ##
## ************* ##
process = 0     #0=>jcache, 1=>grep
epicsPVs = ['MSUPERBIGBITEM', 'MBIGBITEM']
target = 'LH2'
replaypass = 1
outfile = 'epout/mag_redout_per_run_{}.csv'.format(target)
logfile = 'epout/grab_epics_PV_from_evio_{}.log'.format(target)
# ----- **  ---- #

# important directory paths
DATA_DIR = '/mss/halla/sbs/raw'
CACHE_DIR = '/cache/halla/sbs/raw'

def read_runDB(target, replaypass):
    '''Reading run DB'''
    runDB  = '../../DB/good_runList_GMn_nTPE_{}_pass_{}.csv'.format(target, replaypass)
    return pd.read_csv(runDB)

def jcache_0th_seg(runlist, logfile):
    '''jcache 0th element of all run in runlist'''
    with open(logfile, 'w') as f:
        for run in runlist:
            jcachecmd = './jcache_0th_seg.sh {}'.format(run)
            out = os.popen(jcachecmd).read()
            f.write(out + '\n')

# creating run list
runlist = read_runDB(target, replaypass)['runnum']

# jcache the runs
jcache_0th_seg(runlist, logfile)

# Let's jcashe the 0th segment of all runs in runlist
#jcachecmd = './jcache_0th_seg.sh {}'.format(runlist['runnum'][0])
#print(jcachecmd)
#for run in runlist['runnum']:
#    print(run)



#jcachecmd = './jcache_0th_seg.sh 12064'
#jcachecmd = 'grep -ai MSUPERBIGBITEM {}/e1209019_11593.evio.0.0'.format(CACHE_DIR)
#var = os.popen(jcachecmd).read()
#print(var)


