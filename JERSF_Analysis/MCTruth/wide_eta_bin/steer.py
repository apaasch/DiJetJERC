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

def ExecuteCommand(cmd, show=True):
	if show: print(cmd)
	os.system(cmd)

def RemoveTexFile(file):
        ExecuteCommand('rm '+file, show=False)

def ModifyFile(fOld, fNew, old, new):
	with open(fOld, 'r') as f:
		filedata = f.read()
	for (o,n) in zip(old, new):
		filedata = filedata.replace(o, n)
	with open(fNew, 'w') as f:
		f.write(filedata)

def main():

    # study = ["eta_common"]
    # label = ["AK4Puppi_v11"]
    # corr = ["Summer20UL18_V11"]

    study = ["eta_common"]
    label = ["AK4Puppi_v15"]
    corr = ["Summer20UL18_V2"]

    old = ["-STUDY-", "-LABEL-", "-CORR-"]

    path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/MCTruth/wide_eta_bin/"

    file = "Response_MCTruth.cxx"

    for s in study:
        for l in label:
            for c in corr:
                add = "files/"+s+"/"+c+"/"+l+"/"
                ExecuteCommand("mkdir -p "+add)
                new = [s, l, c]
                ModifyFile(path+file, path+add+file, old, new);

# ------------------------------------------------------------------------------
#
if __name__ == "__main__":
    main()
