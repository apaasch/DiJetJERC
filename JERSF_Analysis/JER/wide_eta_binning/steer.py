 #!/usr/bin/env python
import sys
import os
import time
import numpy as np

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/conf/")
from utils import *

def getLabel(sample, year):
    LABEL_LUMI_INV_FB = "[MC 106X] "+year
    # if sample == "A":
    #     LABEL_LUMI_INV_FB += "Run2018A 14.00 fb^{-1}"
    # elif sample == "B":
    #     LABEL_LUMI_INV_FB += "Run2018B 7.10 fb^{-1}"
    # elif sample == "C":
    #     LABEL_LUMI_INV_FB += "Run2018C 6.94 fb^{-1}"
    # elif sample == "D":
    #     LABEL_LUMI_INV_FB += "Run2018D 31.93 fb^{-1}"
    # elif sample == "ABC":
    #     LABEL_LUMI_INV_FB += "Run2018 28.04 fb^{-1}"
    # elif sample == "ABCD":
    #     LABEL_LUMI_INV_FB += "Run2018 59.97 fb^{-1}"
    # else:
    #     LABEL_LUMI_INV_FB += "(2018)"
    return LABEL_LUMI_INV_FB


def main_function(gaustails=False, shiftForPLI="central", gaustail_num = 0.985):
    outdir = out_path+pattern+QCDsample+"/"+run+"/"
    shiftForPLI_num = 0.0
    ref_shift = 3
    if "barrel_check" in extraText:
        ref_shift = int(extraText[-2])
    if "eta_simple" in MC_file:
        ref_shift = 1
    if gaustails:
        outdir = out_path+year+"/"+newJECVersion+"/"+newJetLabel+"/gaustails_"+str(gaustail_num)+"/"+QCDsample+"/"+run+"/"
    if shiftForPLI=="up":
        outdir = out_path+year+"/"+newJECVersion+"/"+newJetLabel+"/PLI/up/"+QCDsample+"/"+run+"/"
        shiftForPLI_num = 0.25
    if shiftForPLI=="down":
        outdir = out_path+year+"/"+newJECVersion+"/"+newJetLabel+"/PLI/down/"+QCDsample+"/"+run+"/"
        shiftForPLI_num = -0.25
    print "outdir ", outdir
    if os.path.isdir(outdir):
        for el in sorted(os.listdir(outdir)):
            cmd = "rm -fr %s" % (outdir+el)
            a = os.system(cmd)
    else:
        os.makedirs(outdir)
    programm ="mainRun"
    programm ="mainRun_alpha"
    # if "AK8" in outdir: programm += "AK8"
    cmd = "cp %s.cxx %s" % (programm, outdir)
    a = os.system(cmd)
    cmd = "cp functions.C %s" % (outdir)
    a = os.system(cmd)
    cmd = "cp "+os.environ["CMSSW_BASE"]+"/src/UHH2/JERCProtoLab/macros/common_info/common_binning.hpp %s" % (outdir)
    a = os.system(cmd)
    cmd = "cp "+os.environ["CMSSW_BASE"]+"/src/UHH2/JERCProtoLab/macros/common_info/common_binning.cxx %s" % (outdir)
    a = os.system(cmd)
    cmd = "cp "+os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/include/tdrstyle_all.h %s" % (outdir)
    a = os.system(cmd)
    cmd = "cp "+os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/include/constants.h %s" % (outdir)
    a = os.system(cmd)
    os.chdir(outdir)
    os.makedirs("pdfy")
    os.makedirs("pdfy/MCTruth")
    os.makedirs("pdfy/SFs")
    os.makedirs("pdfy/kValues")
    os.makedirs("pdfy/NSC_SFs")
    os.makedirs("pdfy/JERs")
    os.makedirs("pdfy/widths")
    os.makedirs("pdfy/maps")
    os.makedirs("ClosureTest")
    os.makedirs("output")
    os.makedirs("output/asymmetries")
    # print pattern+run
    temp_time=time.time()
    # f = open("log.txt",'w')
    MC_type = '\\"MC\\"'
    data_type = '\\"Data\\"'
    trigger_type = '\\"'+study+'\\"'
    # cmd = 'root -l -b -q "%s%s.cxx(%s, false, %s, %s, %s, %s , %s, %s, %s, %s, %s, %s)" >> log.txt &' % (outdir, programm, year, MC_file, Data_file, LABEL_LUMI_INV_FB, MC_type, data_type, trigger_type, '\\"'+outdir+'\\"', gaustail_num, shiftForPLI_num, ref_shift)
    # cmd = 'root -l -b -q "%s%s.cxx(%s, false, %s, %s, %s, %s , %s, %s, %s, %s, %s, %s)" >> log.txt  ' % (outdir, programm, year, MC_file, Data_file, LABEL_LUMI_INV_FB, MC_type, data_type, trigger_type, '\\"'+outdir+'\\"', gaustail_num, shiftForPLI_num, ref_shift)
    cmd = 'root -l -b -q "%s%s.cxx(\\"%s\\", false, %s, %s, %s, %s , %s, %s, %s, %s, %s, %s)"' % (outdir, programm, year, MC_file, Data_file, LABEL_LUMI_INV_FB, MC_type, data_type, trigger_type, '\\"'+outdir+'\\"', gaustail_num, shiftForPLI_num, ref_shift)
    print "cmd", cmd
    a = os.system(cmd)
    print ("time needed: "+str((time.time()-temp_time))+" s")
    os.chdir(common_path)


source_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/hist_preparation/"
common_path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/"


# year = "2018"
year = "UL16preVFP_split"
# year = "UL16preVFP"
# year = "UL16postVFP"
# year = "UL17"
# year = "UL18"
year = "Legacy"

year = sys.argv[1]


samples = {}
# samples["2018"] = ["A", "B", "C", "D", "ABC", "ABCD"]
# samples["2018"] = ["D", "ABC", "ABCD"]
# samples["UL17"] = ["B", "C", "D", "E", "F","BCDEF"]
# samples["UL18"] = ["ABC", "ABCD", "A", "B", "C", "D"]

# samples["UL16preVFP"] = ["B", "C", "D", "E", "F", "BCD", "EF", "BCDEF"]
# samples["UL16postVFP"] = ["F", "G", "H", "FG", "FGH"]
# samples["UL16preVFP"] = ["BCD", "EF", "BCDEF"]
samples["UL16preVFP_split"] = ["BCDEF"]
samples["UL16preVFP"] = ["BCDEF"]
samples["UL16postVFP"] = ["FGH"]
samples["UL17"] = ["BCDEF"]
samples["UL18"] = ["ABCD"]
# samples["UL18"] = ["ABC"]


# samples["UL16preVFP"] = ["B", "C", "D", "E", "F"]
# samples["UL16postVFP"] = ["FG", "H"]
# samples["UL17"] = ["B", "C", "D", "E", "F"]
# samples["UL18"] = ["A", "B", "C", "D", "AB", "CD", "ABC"]



# samples["UL16preVFP"] = ["EF"]
# samples["UL17"] = ["EF"]



samples["Legacy"] = ["II"]



QCDSamples = {}
QCDSamples["2018"] = ["QCDHT"]
# QCDSamples["UL17"] = ["QCDPt"]
QCDSamples["UL16preVFP_split"] = ["QCDHT"]
QCDSamples["UL16preVFP"] = ["QCDHT"]
QCDSamples["UL16postVFP"] = ["QCDHT"]
QCDSamples["UL17"] = ["QCDHT"]
QCDSamples["UL18"] = ["QCDHT"]
QCDSamples["Legacy"] = ["QCDHT"]

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
JECVersions = {}

# Summer 19 campaign
JECVersions["2018"] = ["Autumn18_V19"]
JECVersions["UL16preVFP_split"] = ["Summer19UL16APV_V3"]
# JECVersions["UL16preVFP"] = ["Summer19UL16APV_V3"]
# JECVersions["UL16postVFP"] = ["Summer19UL16_V2"]
# JECVersions["UL18"] = ["Summer19UL18_V5"]
# JECVersions["UL17"] = ["Summer19UL17_V5"]
JECVersions["Legacy"] = ["Summer19Legacy"]

# Summer 20 campaign
JECVersions["UL16preVFP"] = ["Summer20UL16APV_V2"]
JECVersions["UL16postVFP"] = ["Summer20UL16_V2"]
JECVersions["UL17"] = ["Summer20UL17_V2"]
JECVersions["UL18"] = ["Summer20UL18_V2"]

# JetLabels=["AK4CHS", "AK8Puppi", "AK4Puppi"]
# JetLabels=["AK4CHS"]
JetLabels=["AK4Puppi"]

# systematics=["", "PU", "JEC", "alpha", "JER"]
#systematics=["", "PU", "JEC", "alpha", "Prefire"]
# systematics=["", "PLI", "gaus", "PU", "JEC", "alpha", "Prefire", "PS"]
# systematics=["PLI", "gaus", "PU", "JEC", "alpha", "Prefire"]
# systematics=["PU", "JEC", "alpha", "Prefire"]
# systematics=["PU", "JEC", "alpha", "JER"]
# systematics=["", "JER"]
# systematics=["Prefire"]
# systematics=["JEC", "PU", "Prefire"]
# systematics=["PU", "alpha"]
systematics=[""]
# systematics=["PLI", "gaus"]

dirs_sys = ["", "up", "down"]
dirs_PS = [p+d+'_'+f for p in ['ISR'] for d in ['up'] for f in ['2']]
# dirs_PS = [p+d+'_'+f for p in ['FSR','ISR'] for d in ['up', 'down'] for f in ['sqrt2']]
if 'PS' in systematics:
    print(dirs_PS)
dirs = dirs_sys

studies = []
# studies.append("Standard")
# studies.append("L1L2Residual")
# studies.append("L1L2")
# studies.append("eta_JER")
# studies.append("eta_JER_fine")
# studies.append("eta_JER_default")
#studies.append("eta_common_default")
# studies.append("eta_common_fine_v1")
# studies.append("eta_common_fine_v2")
# studies.append("eta_common_fine_v3")
# studies.append("eta_common_fine_v4")
# studies.append("eta_common_fine_v5_prescale")
studies.append("eta_common_fine_v5_highalpha")
# studies.append("eta_common_fine")
# studies.append("eta_common_central_fine_v2")
# studies.append("eta_common")
# studies.append("eta_simple")

for extraText in [""]:
    for study in studies:
        study_out = study
        if len(sys.argv)==3:
            study_out += '_'+sys.argv[2]
        out_path  = common_path+"file/"+study_out+"/"+extraText
        for newJECVersion in JECVersions[year]:
            for newJetLabel in JetLabels:
                for syst in systematics:
                    dirs = dirs_PS if "PS" in syst else dirs_sys
                    for dir in dirs:
                        if syst == "JER" and dir != "":
                          continue
                        if syst == "JER" and dir == "":
                          dir = "nominal"
                        if (syst in ["", "alpha", "PLI", "gaus"] and dir != "") or (not syst in ["", "alpha", "PLI", "gaus"] and dir == ""):
                            continue
                        pattern = year+"/"+newJECVersion+"/"+newJetLabel
                        if syst == "" or syst == "PLI" or syst == "gaus":
                            pattern += "/standard/"
                        elif syst == "alpha":
                            pattern += "/"+syst+"/"
                        else:
                            pattern += "/"+syst+"/"+dir+"/"
                        print "pattern ", pattern
                        for QCDsample in QCDSamples[year]:
                            for sample in samples[year]:
                                run = "Run"+sample
                                LABEL_LUMI_INV_FB=getLabel(sample, year)
                                LABEL_LUMI_INV_FB = '\\"'+LABEL_LUMI_INV_FB+'\\"'
                                MC_file   = '\\"'+source_path+"MC/wide_eta_bin/file/"+study+"/"+pattern.replace("/standard","")+QCDsample+"_"+year+extraText+"/histograms_mc_incl_full.root"+'\\"'
                                if 'prescale' in study:
                                    MC_file = MC_file.replace("_prescale", "")
                                    print 'Set back MC input file:\n',MC_file
                                Data_file = '\\"'+source_path+"data/wide_eta_bin/file/"+study+"/"+pattern.replace("/standard","")+run+"_"+year+extraText+"/histograms_data_incl_full.root"+'\\"'
                                if 'PS' in syst:
                                    Data_file = Data_file.replace("/"+syst+"/"+dir, "")
                                    print Data_file
                                print MC_file, Data_file
                                if not os.path.isfile(str(MC_file.replace("\\","").strip("\""))) or not os.path.isfile(str(Data_file.replace("\\","").strip("\""))):
                                    continue
                                # print MC_file, Data_file
                                if not syst in ["PLI", "gaus"]:
                                  main_function(gaustails=False)
                                if "PLI" == syst:
                                  main_function(gaustails=False, shiftForPLI="up")
                                  main_function(gaustails=False, shiftForPLI="down")
                                if "gaus" == syst:
                                  main_function(gaustails=True, shiftForPLI="central", gaustail_num = 0.95)
                                  main_function(gaustails=True, shiftForPLI="central", gaustail_num = 0.87)
                                    ## for gaustail_num in np.arange(0.8,1.0,0.005):
                                    ##    main_function(gaustails=True, shiftForPLI="central", gaustail_num=gaustail_num)
