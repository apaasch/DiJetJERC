#!/usr/bin/env python
import sys
import os
import time
import argparse
from glob import glob

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/conf/")
from utils import *

# Add the options parser
def parse_args():
    parser = argparse.ArgumentParser(description="Your script description")

    parser.add_argument("-y", "--year", required=True, dest="year", help="Year (2022preEE, 2022postEE, 2023)")
    parser.add_argument("-p", "--pt", dest="ptbins", default="default", help="PT bins")
    parser.add_argument("-a", "--alpha", dest="alphabins", default="", help="Alpha bins")
    parser.add_argument("-add", "--addition", dest="addition", default="", help="Addition")
    parser.add_argument("--prescale", dest="prescale", action="store_true", help="Enable prescale")

    return parser.parse_args()

def main_program(path="", list_path="", out_path="", year="", study="", ptbins="", abins="", JECVersions=[], JetLabels=[], systematics=[], samples=[]):
  isRunII = year=="Legacy"
  isRun3 = "202" in year
  is23 = year=="2023"
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
                if not is23:
                  for el in sample:
                    if not "Run"+el in file_ and not isRunII: continue
                    outputfile.write(file_+"\n")
                else:
                    if not "Run"+sample in file_: continue
                    if sample=="C" and "C_v4" in file_: continue
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

args = parse_args()

# Update the variables with the command line arguments
year = args.year
ptbins = args.ptbins
abins = args.alphabins
addition = args.addition
prescale = args.prescale

common_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/hist_preparation/data/"

samples = {}
# samples["2018"] = ["A", "B", "C", "D", "ABC", "ABCD"]
# samples["2018"] = ["ABC", "D", "ABCD"]
samples["UL16preVFP_split"] = ["BCD", "EF", "BCDEF"]
# samples["UL16preVFP"] = ["BCD", "EF", "BCDEF"]
samples["UL16preVFP"] = ["BCDEF"]
# samples["UL16preVFP"] = ["C"]
# samples["UL16postVFP"] = ["F", "G", "H", "FG", "FGH"]
samples["UL16postVFP"] = ["FGH"]
# samples["UL17"] = ["B", "C", "D", "E", "F","BCDEF"]
samples["UL17"] = ["BCDEF"]
# samples["UL17"] = ["BCDEF"]
# samples["UL18"] = ["A", "B", "C", "D", "ABC", "ABCD"]
# samples["UL18"] = ["ABCD", "A", "B", "C", "D"]
samples["UL18"] = ["ABCD"]
samples["2022postEE"] = ["FG"]
# samples["2023"] = ["C0_v1", "C0_v2", "C0_v3", "C1_v1", "C1_v2", "C1_v3", "C0", "C1"]
samples["2023"] = ["C", "C_v4"]
samples["Legacy"] = ["II"]

JECVersions = {}
JECVersions["2018"] = ["Autumn18_V19"]
JECVersions["UL16preVFP_split"] = ["Summer19UL16APV_V3"]
JECVersions["UL16preVFP"] = ["Summer20UL16APV_V2"]
JECVersions["UL16postVFP"] = ["Summer20UL16_V2"]
JECVersions["UL17"] = ["Summer20UL17_V2"]
JECVersions["UL18"] = ["Summer20UL18_V2"]
JECVersions["2022postEE"] = ["Summer22EEPrompt22_V1"]
JECVersions["2023"] = ["Winter23Prompt23_V1"]

JECVersions["Legacy"] = ["Summer19Legacy"]

# JetLabels = ["AK4CHS", "AK8Puppi", "AK4Puppi"]
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
    if ptbins: out += '_' + ptbins
    if abins: out += '_' + abins
    if prescale: out += '_prescale'
    if addition: out += '_' + addition

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
