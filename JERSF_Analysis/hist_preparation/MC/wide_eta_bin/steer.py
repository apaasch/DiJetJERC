#!/usr/bin/env python
import sys
import os
import time

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/conf/")
from utils import *

def main_program(path="", list_path="", out_path="", year="", JECVersions=[], JetLabels=[], systematics=[], samples=[]):
  list_path_=list_path
  out_path_=out_path
  for newJECVersion in JECVersions:
    for newJetLabel in JetLabels:
      for sys in set(systematics):
        if sys == "alpha":
          alpha_cut = 10
        else:
          alpha_cut = 15
        for dir in ["", "up", "down"]:
          if sys == "JER" and dir != "":
            continue
          if sys == "JER" and dir == "":
            dir = "nominal"
            print sys, dir
          if (sys == "" and dir != "") or (sys == "alpha" and dir != "") or ((sys != "" and sys != "alpha") and dir == ""):
            continue
          pattern = newJECVersion+"/"+newJetLabel+"/"+sys+"/"+dir+"/"
          if sys == "" or sys == "alpha":
            pattern = newJECVersion+"/"+newJetLabel+"/"
          source_path = path+pattern
          if not os.path.isdir(source_path):
            continue
          if sys == "alpha":
            pattern = newJECVersion+"/"+newJetLabel+"/"+sys+"/"
          if not os.path.isdir(list_path_+pattern):
            os.makedirs(list_path_+pattern)
          for sample in samples:
            run_list = list_path_+pattern+"file_QCD"+sample+".txt"
            with open(run_list, "w") as outputfile:
              for writable in sorted(os.listdir(source_path)):
                if not os.path.isfile(source_path+writable):
                  continue
                if sample in writable and ".root" in writable:
                  outputfile.write(source_path+writable+"\n")
            if not os.path.isfile(run_list):
              continue
            outdir = out_path_+pattern+"QCD"+sample+"_"+year+"/"
            if not os.path.isdir(outdir):
              os.makedirs(outdir)
            print "RUNNING ON ", run_list
            temp_time=time.time()
            cmd = "cp MySelector_full_DiJet.C MySelector.C"
            a = os.system(cmd)
            cmd = 'sed -i -e """s/jet_thr=15/jet_thr=%s/g" MySelector.C' % (alpha_cut)
            a = os.system(cmd)
            if study == "LowPtJets":
              cmd = 'sed -i -e """s/pt_bins_Si/pt_bins_MB/g" MySelector.C'
              a = os.system(cmd)
            cmd = "cp Makefile %s" % (outdir)
            a = os.system(cmd)
            cmd = "cp *.C %s" % (outdir)
            a = os.system(cmd)
            cmd = "cp *.h %s" % (outdir)
            a = os.system(cmd)
            cmd = "cp "+os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/include/"+"constants.h %s" % (outdir)
            a = os.system(cmd)
            os.chdir(outdir)
            time.sleep(3)
            cmd = "make clean"
            a = os.system(cmd)
            cmd = "make"
            a = os.system(cmd)
            logfilename = "log.txt"
            f = open(logfilename,'w')
            cmd = './Analysis.x %s >> log.txt &' % (run_list)
            command = [outdir+"Analysis.x", run_list, outdir, year]
            list_processes.append(command)
            list_logfiles.append(outdir+"log.txt")
            f.close()
            os.chdir(common_path+"wide_eta_bin/")
            print ("time needed: "+str((time.time()-temp_time))+" s")


USER = os.environ["USER"]

inputdir = "DiJetJERC_DiJetHLT"
#year = "2018"
year = "UL17"

common_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/hist_preparation/MC/"

samples = {}
samples["2018"] = ["HT"]
# samples["UL17"] = ["Pt"]
# samples["UL17"] = ["HT"]
samples["UL17"] = ["Pt", "HT"]

JECVersions = {}
JECVersions["UL17"] = ["Summer19UL17_V1_ComplexL1","Summer19UL17_V1_SimpleL1"]
JECVersions["UL17"] = ["Summer19UL17_V1_ComplexL1"]
JECVersions["2018"] = ["Autumn18_V19"]
# JetLabels = ["AK4CHS", "AK8Puppi", "AK4Puppi"]
JetLabels = ["AK4CHS"]
# systematics = ["", "alpha","PU", "JEC", "JER"]
systematics = ["", "alpha","PU", "JEC"]
# systematics = ["", "alpha","PU"]
# systematics = [""]

list_processes = []
list_logfiles = []

studies = []
studies.append("Standard")
# studies.append("L1L2Residual")
# studies.append("L1L2")
# studies.append("Simplified")
# studies.append("PuJetId")

for study in studies:
    list_path   = common_path+"lists/"+study+"/"+year+"/"
    out_path    = common_path+"wide_eta_bin/file/"+study+"/"+year+"/"
    os.chdir(common_path+"wide_eta_bin/")

    path = "/nfs/dust/cms/user/"+USER+"/sframe_all/"+inputdir+"/"+year+"/"+study+"/"

    main_program(path, list_path, out_path, year, JECVersions[year], JetLabels, systematics, samples[year])


for i in list_processes:
  print i

print len(list_processes)

parallelise(list_processes, 3, list_logfiles)
