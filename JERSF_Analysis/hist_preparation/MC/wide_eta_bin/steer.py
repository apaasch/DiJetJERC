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
    parser.add_argument("-p", "--pt", dest="ptbins", default="default", help="PT bins (e.g. default, fine_v1,)")
    parser.add_argument("-a", "--alpha", dest="alphabins", default="", help="Alpha bins (e.g. finealpha)")
    parser.add_argument("-add", "--addition", dest="addition", default="", help="Addition to study name")
    parser.add_argument("-s", "--sample", dest="sample", default="", help="If only one QCD sample should be running")

    return parser.parse_args()

def main_program(path="", list_path="", out_path="", year="", study="", ptbins="", abins="", JECVersions=[], JetLabels=[], systematics=[], samples=[]):
  isRunII = year=="Legacy"
  isRun3 = "202" in year
  is23 = year=="2023"
  list_path_=list_path
  out_path_=out_path
  dirs = ["", "up", "down"]
  for newJECVersion in JECVersions:
    for newJetLabel in JetLabels:
      for sys in set(systematics):
        if sys == "alpha":
          alpha_cut = 10
        else:
          alpha_cut = 15
        for dir in dirs:
          # if sys == "JER" and dir != "":
          #   continue
          if sys == "JER" and dir == "":
            dir = "nominal"
            print sys, dir
          if (sys == "" and dir != "") or (sys == "alpha" and dir != "") or ((sys != "" and sys != "alpha") and dir == ""):
            continue
          pattern = newJECVersion+"/"+newJetLabel+"/"+sys+"/"+dir+"/"
          if sys == "" or sys == "alpha" or sys == "PS":
            pattern = newJECVersion+"/"+newJetLabel+"/"
          source_path = path+pattern
          if isRunII:
              source_path = source_path.replace(newJECVersion,"*").replace(year,"*")
          if not os.path.isdir(source_path) and not isRunII:
            continue
          if sys == "alpha":
            pattern = newJECVersion+"/"+newJetLabel+"/"+sys+"/"
          if sys == "PS":
            pattern = newJECVersion+"/"+newJetLabel+"/"+sys+"/"+dir+"/"
          if not os.path.isdir(list_path_+pattern):
            os.makedirs(list_path_+pattern)
          for sample in samples:
            run_list = list_path_+pattern+"file_QCD"+sample+".txt"
            with open(run_list, "w") as outputfile:
              fname="uhh2.AnalysisModuleRunner.MC.*root"
              for file_ in sorted(glob.glob(source_path+fname)):
                if not sample in file_ and not isRunII: continue
                if not binning == "":
                  if binning not in file_:
                    continue
                  else:
                    print "Run on ",file_, "with binning",binning
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
            command = [outdir+"Analysis.x", run_list, outdir, year, study, ptbins, abins, dir]
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

global binning
binning = args.sample

common_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/hist_preparation/MC/"

samples = {}
samples["2022preEE"] = ["HT"]
samples["2022postEE"] = ["HT"]
samples["2023"] = ["PT"]
samples["2023preBPix"] = ["HT"]
samples["2023postBPix"] = ["HT"]

JECVersions = {}
JECVersions["2022preEE"] = ["Summer22_22Sep2023_V2"]
JECVersions["2022postEE"] = ["Summer22EE_22Sep2023_V2"]
JECVersions["2023"] = ["Winter23Prompt23_V1"]
JECVersions["2023preBPix"] = ["Summer23Prompt23_V1"]
JECVersions["2023postBPix"] = ["Summer23BPixPrompt23_V1"]

# JetLabels = ["AK4CHS", "AK8Puppi", "AK4Puppi"]
JetLabels = ["AK4Puppi"]
# systematics = ["", "alpha","PU", "JEC", "JER", "Prefire", "PS"]
systematics = ["JER"]

list_processes = []
list_logfiles = []

studies = []
studies.append("eta_common")

for study in studies:
    list_path   = common_path+"lists/"+study+"/"+year+"/"
    # ###############################################
    # eta bins in PreSelection only for reference jet
    # Don't need to run PreSel again
    if "eta_calo" in study:
      list_path   = common_path+"lists/eta_common/"+year+"/"
    out = study
    if ptbins: out += '_' + ptbins
    if abins: out += '_' + abins
    if binning: out += '_' + binning
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

parallelise(list_processes, 12, list_logfiles)
