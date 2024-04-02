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
  is23 = "2023" in year
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
                if is23:
                  if "Cv" in sample:
                    if not "Run"+sample in file_ and not isRunII: continue
                    outputfile.write(file_+"\n")
                  else:
                    for el in sample:
                      if not "Run"+el in file_ and not isRunII: continue
                      outputfile.write(file_+"\n")
                else:
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

args = parse_args()

# Update the variables with the command line arguments
year = args.year
ptbins = args.ptbins
abins = args.alphabins
addition = args.addition
prescale = args.prescale

common_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/hist_preparation/data/"

samples = {}
samples["2022preEE"] = ["CD"] # CHT is considered with C
# samples["2022preEE"] = ["CD"] # CHT is considered with C
samples["2022postEE"] = ["EFG"]
samples["2023preBPix"] = ["C", "Cv123", "Cv4"] # Cv123, Cv4
samples["2023postBPix"] = ["D"]

JECVersions = {}
JECVersions["2022preEE"] = ["Summer22_19Dec2023_V2"]
JECVersions["2022postEE"] = ["Summer22EE_22Sep2023_V2"]
JECVersions["2023"] = ["Winter23Prompt23_V1"]
JECVersions["2023preBPix"] = ["Summer23Prompt23_V1"]
JECVersions["2023postBPix"] = ["Summer23BPixPrompt23_V1"]

JetLabels = ["AK4Puppi"]
# systematics = ["", "alpha","PU", "JEC", "JER"]
# systematics = ["alpha","PU", "JEC"]
systematics = [""]

list_processes = []
list_logfiles = []

studies = []
studies.append("eta_common")

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
