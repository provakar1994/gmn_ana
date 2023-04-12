# This script grabs the EPICS PV values directly from EVIO files.
# The idea is to first jcache the 0th segment of EVIO files for all
# GMn LH2 or LD2 runs by reading the run DB and then grab the user
# chosen EPICS PVs per run and write them in an output CSV file.
#
# Interdependency: jcache_0th_seg.sh is needed to do jcache
# -------
# P. Datta <pdbforce@jlab.org> Created 04-12-2023

import os
import pandas as pd

## ************* ##
##  User inputs  ##
## ************* ##
process = 1     #0=>jcache, 1=>grep
target = 'LH2'
replaypass = 1
logfile = 'epout/grab_epics_PV_from_evio_{}.log'.format(target)
# if process = 1
epicsPVs = ['MSUPERBIGBITEM', 'MBIGBITEM']
outfile = 'epout/mag_PVs_per_run_{}.csv'.format(target)
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
    print('Starting process 0..')
    with open(logfile, 'w') as f:
        for run in runlist:
            jcachecmd = './jcache_0th_seg.sh {}'.format(run)
            out = os.popen(jcachecmd).read()
            f.write(out + '\n')

def grab_epicsPV(runlist, epicsPVs, outfile):
    '''Grabbing EPICS PV values from evio file per run'''
    print('Starting process 1..')
    with open(outfile, 'w') as f:
        f.write('runnum')
        for pv in epicsPVs:
            f.write(',{}'.format(pv))
        itr = 0
        for run in runlist:
            itr += 1
            f.write('\n{}'.format(run))
            eviofile = '{}/e1209019_{}.evio.0.0'.format(CACHE_DIR, run)
            for pv in epicsPVs:
                grabpvcmd = 'grep -a -m 1 {} {}'.format(pv, eviofile)
                temp = os.popen(grabpvcmd).read()
                if (temp.find(pv) != -1):
                    val = temp.split(" ", 1)[1].strip()
                    f.write(',' + val)
                else:
                    val = '-9999'
                    f.write(',' + val)   
            print('Finished run: {}. Progress {}/{}'.format(run, itr, len(runlist)))

# creating run list
runlist = read_runDB(target, replaypass)['runnum']

if (process==0):
    # jcache the runs
    jcache_0th_seg(runlist, logfile)
elif (process==1):
    # grab values for EPICS PVs
    grab_epicsPV(runlist, epicsPVs, outfile)
else:
    print('Entered process ID is invalid!')


