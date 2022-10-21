#!/usr/bin/env python
import sys
import os
import time

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/conf/")
from utils import *

def main_program(path="", list_path="", out_path="", year="", study="", binning="", JECVersions=[], JetLabels=[], systematics=[], samples=[]):
  isRunII = year=="Legacy"
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
          if isRunII:
              source_path = source_path.replace(newJECVersion,"*").replace(year,"*")
          if not os.path.isdir(source_path) and not isRunII:
            continue
          if sys == "alpha":
            pattern = newJECVersion+"/"+newJetLabel+"/"+sys+"/"
          if not os.path.isdir(list_path_+pattern):
            os.makedirs(list_path_+pattern)
          for sample in samples:
            run_list = list_path_+pattern+"file_QCD"+sample+".txt"
            with open(run_list, "w") as outputfile:
              for file_ in sorted(glob.glob(source_path+"uhh2.AnalysisModuleRunner.MC.*root")):
                if not sample in file_ and not isRunII: continue
                outputfile.write(file_+"\n")
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
            command = [outdir+"Analysis.x", run_list, outdir, year, study, binning]
            list_processes.append(command)
            list_logfiles.append(outdir+"log.txt")
            f.close()
            os.chdir(common_path+"wide_eta_bin/")
            print ("time needed: "+str((time.time()-temp_time))+" s")


USER = os.environ["USER"]

inputdir = "DiJetJERC_DiJetHLT"

# year = "2018"
# year = "UL16preVFP_split"
# year = "UL16preVFP"
#year = "UL16postVFP"
#year = "UL17"
year = "UL18"
# year = "Legacy"
if len(sys.argv)<2:
    sys.exit("I need at least a year")
year = sys.argv[1]

common_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/hist_preparation/MC/"

samples = {}
# samples["2018"] = ["HT"]
# samples["2018"] = ["Pt"]
samples["UL16preVFP_split"] = ["HT"]
samples["UL16preVFP"] = ["HT"]
samples["UL16postVFP"] = ["HT"]
samples["UL17"] = ["HT"]
samples["UL18"] = ["HT"]
samples["Legacy"] = ["HT"]

JECVersions = {}
JECVersions["2018"] = ["Autumn18_V19"]
JECVersions["UL16preVFP_split"] = ["Summer19UL16APV_V3"]
JECVersions["UL16preVFP"] = ["Summer20UL16APV_V2"]
JECVersions["UL16postVFP"] = ["Summer20UL16_V2"]
JECVersions["UL17"] = ["Summer20UL17_V2"]
JECVersions["UL18"] = ["Summer20UL18_V2"]

JECVersions["Legacy"] = ["Summer19Legacy"]

# JetLabels = ["AK4CHS", "AK8Puppi", "AK4Puppi"]
# JetLabels = ["AK4Puppi", "AK8Puppi"]
JetLabels = ["AK4Puppi"]
# systematics = ["", "alpha","PU", "JEC", "JER", "Prefire"]
# systematics = ["alpha","PU", "JEC", "Prefire"]
# systematics = ["", "alpha","PU", "JEC", "Prefire"]
# systematics = ["", "alpha", "JEC", "JER"]
# systematics = ["PU", "JEC"]
# systematics = ["PU"]
# systematics = ["", "alpha","PU", "JEC"]
# systematics = ["", "PU", "JEC"]
# systematics = ["Prefire", "PU", "JEC"]
# systematics = ["Prefire"]
systematics = [""]

list_processes = []
list_logfiles = []

studies = []
# studies.append("Standard")
# studies.append("L1L2Residual")
# studies.append("L1L2")
# studies.append("Simplified")
# studies.append("PuJetId")
studies.append("eta_common")
# studies.append("eta_common")
# studies.append("eta_simple")

for study in studies:
    list_path   = common_path+"lists/"+study+"/"+year+"/"
    out = study
    binning = ''
    if len(sys.argv) == 3:
        binning = sys.argv[2]
        out += '_'+binning
    out_path    = common_path+"wide_eta_bin/file/"+out+"/"+year+"/"
    os.chdir(common_path+"wide_eta_bin/")

    path = "/nfs/dust/cms/user/"+USER+"/sframe_all/"+inputdir+"/"+year+"/"+study+"/"

    main_program(path, list_path, out_path, year, study, binning, JECVersions[year], JetLabels, systematics, samples[year])


for i in list_processes:
  print i

print len(list_processes)
# print "PRO - ", list_processes
# print ""
# print "LOG - ", list_logfiles

parallelise(list_processes, 8, list_logfiles)
