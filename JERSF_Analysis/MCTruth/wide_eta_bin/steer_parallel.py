#!/usr/bin/env python

from datetime import datetime
import sys, time, subprocess
#from ROOT import *
from multiprocessing import Pool, Value
import random
import shutil
import glob
import os
import shutil

global year
n_jobs = 9
Ndone = 0 # keep track of finished jobs by each worker

global runpath

def ExecuteCommand(cmd, show=True):
	if show: print(cmd)
	os.system(cmd)

def RemoveTexFile(file):
        ExecuteCommand('rm '+file, show=False)

def ModifyFile(fOld, fNew, old, new, sample):
	with open(fOld, 'r') as f:
		filedata = f.read()
	for (o,n) in zip(old, new):
		filedata = filedata.replace(o, n)
        filedata = filedata.replace("Response_MCTruth_Parallel", "Response_MCTruth_"+sample)
	with open(fNew, 'w') as f:
		f.write(filedata)

def main():
    global n_jobs
    global runpath      # Needed to modify global copy of n_jobs

    pool = Pool(processes = n_jobs)
    USER = os.environ["USER"]

    common_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/MCTruth/wide_eta_bin"

    samples = [
        "QCDHT50to100",
        "QCDHT100to200",
        "QCDHT200to300",
        "QCDHT300to500",
        "QCDHT500to700",
        "QCDHT700to1000",
        "QCDHT1000to1500",
        "QCDHT1500to2000",
        # "QCDHT2000toInf",
    ]

    study = ["eta_common"]
    label = ["AK4Puppi_v15"]
    corr = ["Summer20UL18_V2"]

    old = ["-STUDY-", "-LABEL-", "-CORR-", "-SAMPLE-"]
    path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/MCTruth/wide_eta_bin/"
    file = "Response_MCTruth_Parallel.cxx"

    add = ''
    for s in study:
        for l in label:
            for c in corr:
                add = "files/parallel/"+s+"/"+c+"/"+l+"/"
                ExecuteCommand("mkdir -p "+add)
                ExecuteCommand("cp CombineFiles.py "+add+"CombineFiles.py")
                for sa in samples:
                    new = [s, l, c, sa]
                    fnew = file.replace("Parallel", sa)
                    ModifyFile(path+file, path+add+fnew, old, new, sa);

    runpath = path+add
    os.chdir(runpath)
    print os.getcwd()

    file_list = [runpath+"Response_MCTruth_"+sample+".cxx" for sample in samples]
    job_lists = []
    for i in range(n_jobs):
        job_lists.append( create_job_list(i, n_jobs, file_list) )

    print(job_lists)
    ndirs=len(samples)
    result = pool.map_async(run_job, job_lists, chunksize=1)
    pool.close()
    pool.join()

# ------------------------------------------------------------------------------
#
def run_job(dirlist):
    print dirlist
    for sample in dirlist:
        cmd = ['root','-l','-b','-q',sample]
        logfile = sample.replace(".cxx", ".txt")
        logfile = logfile.replace("Response_MCTruth", "log")
        subprocess.call(cmd, stdout=open(logfile,"w"), stderr=subprocess.STDOUT)
    return 1

# ------------------------------------------------------------------------------
# for every worker an array containing different dirs is created
# it is taken care of distributing the job as equally as possible
def create_job_list(i, nworkers, dirs):
    integer_div = len(dirs) // nworkers
    rest = len(dirs) - (integer_div * nworkers)
    new_list = []
    # the rest should not just be added to the last worker but distributed equally
    # therefore every worker with i+1 gets one additional job
    if i+1 <= rest:
        min_index = i*integer_div + i
        max_index = (i+1)*integer_div + 1 + i
    else:
        min_index = i*integer_div + rest
        max_index = (i+1)*integer_div + rest

    # loop over total list and only keep element if index is between min and max
    for index, val in enumerate(dirs):
        if index >= min_index and index < max_index:
            new_list.append( dirs[index] )

    return new_list

# ------------------------------------------------------------------------------
# needed to setup global counter
def init(args):
    global Ndone
    Ndone = args
    # global year

# ------------------------------------------------------------------------------
#
if __name__ == "__main__":
    main()
