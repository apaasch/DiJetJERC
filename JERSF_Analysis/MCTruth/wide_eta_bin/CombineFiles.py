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
import ROOT

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
    files = {}
    for sa in samples:
		files[sa] = ROOT.TFile(sa+'.root')


    # nhists = [key.GetName() for key in files[samples[0]].GetListOfKeys()]
    # print len(nhists)

    fout = ROOT.TFile('MCTruth_JER.root', 'recreate')
    hists = {}
    for object in files[samples[0]].GetListOfKeys():
		h = object.GetName()
		hist = files[samples[0]].Get(h)
		hist.Sumw2()
		values = []
		values.append(hist.GetEntries())
		for s in samples:
			if s==samples[0]: continue
			hist2 = files[s].Get(h)
			values.append(hist2.GetEntries())
			hist.Add(hist2, 1)
		values.append(hist.GetEntries())
		print h,values
		fout.cd()
		hist.Write(h)
    fout.Close()




# ------------------------------------------------------------------------------
#
if __name__ == "__main__":
    main()
