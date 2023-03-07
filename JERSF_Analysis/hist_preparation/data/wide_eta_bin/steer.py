#!/usr/bin/env python
import sys
import os
import time

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/conf/")
from utils import *

def main_program(path="", list_path="", out_path="", year="", study="", ptbins="", abins="", JECVersions=[], JetLabels=[], systematics=[], samples=[]):
  isRunII = year=="Legacy"
  list_path_=list_path
  out_path_=out_path
  dirs_sys = ["", "up", "down"]
  for newJECVersion in JECVersions:
    for newJetLabel in JetLabels:
      for sys in set(systematics):
        if sys == "alpha":
          alpha_cut = 10
        else:
          alpha_cut = 15
        dirs = dirs_PS if "PS" in sys else dirs_sys
        for dir in dirs:
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
            run_list = list_path_+pattern+"file_DataRun"+sample+"_"+year+".txt"
            with open(run_list, "w") as outputfile:
              for file_ in sorted(glob.glob(source_path+"uhh2.AnalysisModuleRunner.DATA.*root")):
                for el in sample:
                  if not "Run"+el in file_ and not isRunII: continue
                  outputfile.write(file_+"\n")
            if not os.path.isfile(run_list):
              continue
            outdir = out_path_+pattern+"Run"+sample+"_"+year+"/"
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
            command = [outdir+"Analysis.x", run_list, outdir, year, study, ptbins, abins]
            list_processes.append(command)
            list_logfiles.append(outdir+"log.txt")
            f.close()
            os.chdir(common_path+"wide_eta_bin/")
            print ("time needed: "+str((time.time()-temp_time))+" s")


USER = os.environ["USER"]

inputdir = "DiJetJERC_DiJetHLT"

#year = "2018"
# year = "UL16preVFP_split"
# year = "UL16preVFP"
# year = "UL16postVFP"
# year = "UL17"
year = "UL18"
if len(sys.argv)<2:
    sys.exit("I need at least a year")
year = sys.argv[1]
# year = "Legacy"

common_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/hist_preparation/data/"

samples = {}
# samples["2018"] = ["A", "B", "C", "D", "ABC", "ABCD"]
# samples["2018"] = ["ABC", "D", "ABCD"]
samples["UL16preVFP_split"] = ["BCD", "EF", "BCDEF"]
samples["UL16preVFP"] = ["BCD", "EF", "BCDEF"]
# samples["UL16preVFP"] = ["BCDEF"]
# samples["UL16preVFP"] = ["C"]
samples["UL16postVFP"] = ["F", "G", "H", "FG", "FGH"]
# samples["UL16postVFP"] = ["FGH"]
samples["UL17"] = ["B", "C", "D", "E", "F","BCDEF"]
# samples["UL17"] = ["BCDEF"]
# samples["UL18"] = ["A", "B", "C", "D", "ABC", "ABCD"]
# samples["UL18"] = ["ABCD", "A", "B", "C", "D"]
samples["UL18"] = ["ABCD"]
samples["Legacy"] = ["II"]

JECVersions = {}
JECVersions["2018"] = ["Autumn18_V19"]
JECVersions["UL16preVFP_split"] = ["Summer19UL16APV_V3"]
JECVersions["UL16preVFP"] = ["Summer20UL16APV_V2"]
JECVersions["UL16postVFP"] = ["Summer20UL16_V2"]
JECVersions["UL17"] = ["Summer20UL17_V2"]
JECVersions["UL18"] = ["Summer20UL18_V2"]

JECVersions["Legacy"] = ["Summer19Legacy"]

# JetLabels = ["AK4CHS", "AK8Puppi", "AK4Puppi"]
# JetLabels = ["AK4Puppi_v11"]
JetLabels = ["AK4Puppi"]
# systematics = ["", "alpha","PU", "JEC", "JER", "Prefire"]
# systematics = ["", "alpha","PU", "JEC", "Prefire"]
# systematics = ["", "alpha","PU", "JEC"]
# systematics = ["PU", "JEC", "Prefire"]
# systematics = ["JER"]
# systematics = ["", "alpha","PU"]
systematics = [""]

list_processes = []
list_logfiles = []

studies = []
# studies.append("Standard")
# studies.append("L1L2Residual")
# studies.append("L1L2")
# studies.append("Simplified")
# # studies.append("PuJetId")
# studies.append("eta_JER")
studies.append("eta_common")
# studies.append("eta_calo")
# studies.append("eta_simple")

global dirs_PS
dirs_PS = [p+d+'_'+f for p in ['FSR','ISR'] for d in ['up', 'down'] for f in ['4', '2']]
if 'PS' in systematics:
    print(dirs_PS)

for study in studies:
    # eta_bins = "_".join(study.split("_")[:2])
    # print("Using eta bins:",eta_bins)
    list_path   = common_path+"lists/"+study+"/"+year+"/"
    # ###############################################
    # eta bins in PreSelection only for reference jet
    # Don't need to run PreSel again
    if "eta_calo" in study:
      list_path   = common_path+"lists/eta_common/"+year+"/"
    out = study
    ptbins = ''
    abins = ''
    if len(sys.argv) >= 3:
        ptbins = sys.argv[2]
        out += '_'+ptbins
    if len(sys.argv) >= 4:
        if "alpha" in sys.argv[3]: 
          abins = sys.argv[3]
        out += '_'+sys.argv[3]
    if len(sys.argv) >= 5:
        out += '_'+sys.argv[4]
    print out
    out_path    = common_path+"wide_eta_bin/file/"+out+"/"+year+"/"
    os.chdir(common_path+"wide_eta_bin/")

    path = "/nfs/dust/cms/user/"+USER+"/sframe_all/"+inputdir+"/"+year+"/"+study+"/"
    if "eta_calo" in study:
        path = "/nfs/dust/cms/user/"+USER+"/sframe_all/"+inputdir+"/"+year+"/eta_common/"

    main_program(path, list_path, out_path, year, study, ptbins, abins, JECVersions[year], JetLabels, systematics, samples[year])


for i in list_processes:
  print i

print len(list_processes)

parallelise(list_processes, 10, list_logfiles)
