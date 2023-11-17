#!/usr/bin/env python3

from createConfigFiles import *

@timeit
def condor_control(original_dir ="./SubmittedJobs/" , JECVersions_Data=["Autumn18_V4"], JetLabels=["AK4CHS"], systematics=["", "PU", "JEC", "JER"], internal_option="-l", processes=[], extratext=""):
    count = 0
    list_processes = []
    nProcess = 48
    time_ = 1
    for newJECVersion in JECVersions_Data:
        for newJetLabel in JetLabels:
            for sys in systematics:
                dirs = ["", "up", "down"]
                for dir in dirs:
                    if sys == "" and dir != "":
                        continue
                    if sys == "JER" and dir != "":
                        continue
                    if sys == "JER" and dir == "":
                        dir = "nominal"
                    path = original_dir+newJECVersion+"/"+newJetLabel+extratext+"/"+sys+"/"+dir+"/"
                    for sample in sorted(os.listdir(path)):
                        if not ".xml" in sample:
                            continue
                        if all(not control in sample for control in processes): continue
                        if internal_option:
                            command = ['sframe_batch.py', internal_option, path+sample]
                        else:
                            command = ['sframe_batch.py', path+sample]
                        command = [path]+command
                        list_processes.append(command)
                        if internal_option == "-f":
                            nProcess = 20
                        if internal_option == "":
                            time_ = 0.5
    print("Start parallelise! - ", len(list_processes))
    parallelise(list_processes, nProcess, cwd=True, time_=time_)


@timeit
def delete_workdir(original_dir ="./SubmittedJobs/" , JECVersions_Data=["Autumn18_V4", "Autumn18_V4"], JetLabels=["AK4CHS", "AK8Puppi"], systematics=["", "PU", "JEC", "JER"],extratext=""):
    add_name = original_dir[original_dir.find("SubmittedJobs")+len("SubmittedJobs"):-1]
    for sample in ["DATA", "QCD"]:
        for newJECVersion in JECVersions_Data:
            for newJetLabel in JetLabels:
                for sys in systematics:
                    dirs = ["", "up", "down"]
                    for dir in dirs:
                        if sys == "" and dir != "":
                            continue
                        if sys == "JER" and dir != "":
                            continue
                        if sys == "JER" and dir == "":
                       	    dir = "nominal"
                        if 'PS' in sys:
                            if 'DATA'==sample:
                                continue
                        path = userPathSframeOutput+"/"+newJECVersion+"/"+newJetLabel+extratext+"/"+sys+"/"+dir+"/"
                        if os.path.isdir(path):
                            for workdir in sorted(os.listdir(path)):
                                if "workdir" in workdir:
                                    cmd = "rm -fr %s" % (path+workdir)
                                    a = os.system(cmd)
                                    print(cmd)
                        path = original_dir+newJECVersion+"/"+newJetLabel+extratext+"/"+sys+"/"+dir+"/"
                        if os.path.isdir(path):
                            for workdir in sorted(os.listdir(path)):
                                if "workdir" in workdir:
                                    cmd = "rm -fr %s" % (path+workdir)
                                    a = os.system(cmd)




def main_program(option="", internal_option="", study="Standard", processes=[], others=[], JECVersions_Data=[], JECVersions_MC=[], JetLabels=[], systematics=[], original_dir="./SubmittedJobs/", original_file="JER2018.xml", year="2018", isMB=False, test_trigger=False, isThreshold=False, isLowPt=False, isECAL=False, extratext=""):
    if option == "new":
        createConfigFiles(study, processes, others, JECVersions_Data, JECVersions_MC, JetLabels, systematics, original_dir, original_file, outdir, year, isMB, test_trigger, isThreshold,isLowPt,isECAL,extratext)
    elif option == "remove" or option == "delete":
        delete_workdir(original_dir, JECVersions_Data, JetLabels, systematics, extratext)
    else:
        condor_control(original_dir, JECVersions_Data, JetLabels, systematics, internal_option, processes, extratext)



##################################################
#                                                #
#                   MAIN Program                 #
#                                                #
##################################################

USER = os.environ["USER"]

try:
    option = sys.argv[1]
except:
    option = ""

if option == "resubmit":
    internal_option = "-r"
elif option == "submit":
    internal_option = "-s"
elif option == "add" or option == "merge":
    internal_option = "-f"
elif option == "list":
    internal_option = "-l"
elif option == "new":
    internal_option = ""
elif option == "remove" or option == "delete":
    internal_option = ""
elif option == "split":
    internal_option = ""
else:
    internal_option = ""


QCD_process= []
Data_process= []

QCD_process.append("QCDHT40to70_2022preEE")
QCD_process.append("QCDHT70to100_2022preEE")
QCD_process.append("QCDHT100to200_2022preEE")
QCD_process.append("QCDHT200to400_2022preEE")
QCD_process.append("QCDHT400to600_2022preEE")
QCD_process.append("QCDHT600to800_2022preEE")
QCD_process.append("QCDHT800to1000_2022preEE")
QCD_process.append("QCDHT1000to1200_2022preEE")
QCD_process.append("QCDHT1200to1500_2022preEE")
QCD_process.append("QCDHT1500to2000_2022preEE")
QCD_process.append("QCDHT2000_2022preEE")

QCD_process.append("QCDPT50to80_2022preEE")
QCD_process.append("QCDPT80to120_2022preEE")
QCD_process.append("QCDPT120to170_2022preEE")
QCD_process.append("QCDPT170to300_2022preEE")
QCD_process.append("QCDPT300to470_2022preEE")
QCD_process.append("QCDPT470to600_2022preEE")
QCD_process.append("QCDPT600to800_2022preEE")
QCD_process.append("QCDPT800to1000_2022preEE")
QCD_process.append("QCDPT1000to1400_2022preEE")
QCD_process.append("QCDPT1400to1800_2022preEE")
QCD_process.append("QCDPT1800to2400_2022preEE")
QCD_process.append("QCDPT2400to3200_2022preEE")
QCD_process.append("QCDPT3200_2022preEE")
QCD_process.append("QCDPT50to80_2022postEE")
QCD_process.append("QCDPT80to120_2022postEE")
QCD_process.append("QCDPT120to170_2022postEE")
QCD_process.append("QCDPT170to300_2022postEE")
QCD_process.append("QCDPT300to470_2022postEE")
QCD_process.append("QCDPT470to600_2022postEE")
QCD_process.append("QCDPT600to800_2022postEE")
QCD_process.append("QCDPT800to1000_2022postEE")
QCD_process.append("QCDPT1000to1400_2022postEE")
QCD_process.append("QCDPT1400to1800_2022postEE")
QCD_process.append("QCDPT1800to2400_2022postEE")
QCD_process.append("QCDPT2400to3200_2022postEE")
QCD_process.append("QCDPT3200_2022postEE")
Data_process.append("DATA_RunB_JetHT_2022")
Data_process.append("DATA_RunC_JetHT_2022")
Data_process.append("DATA_RunC_2022preEE")
Data_process.append("DATA_RunD_2022preEE")
Data_process.append("DATA_RunE_2022preEE")
Data_process.append("DATA_RunF_2022postEE")
Data_process.append("DATA_RunG_2022postEE")

Data_process.append("DATA_RunC_2023")
Data_process.append("DATA_RunC_v4_2023")
# Data_process.append("DATA_RunC0_v1_2023")
# Data_process.append("DATA_RunC0_v2_2023")
# Data_process.append("DATA_RunC0_v3_2023")
# Data_process.append("DATA_RunC0_v4_2023")
# Data_process.append("DATA_RunC1_v1_2023")
# Data_process.append("DATA_RunC1_v2_2023")
# Data_process.append("DATA_RunC1_v3_2023")
# Data_process.append("DATA_RunC1_v4_2023")

# JECVersions_Data = ["Autumn18_V4"]
# JetLabels = ["AK4CHS", "AK8Puppi"]
# systematics = ["", "PU", "JEC", "JER"]

# year = "2023"
year = "2022preEE"
# year = "UL16preVFP"
# year = "UL16postVFP"
# year = "UL17"
# year = "UL18"


studies = []
# studies.append("Standard")
# studies.append("L1L2Residual")
# studies.append("L1L2")
# studies.append("eta_JER")
studies.append("eta_common")
# studies.append("eta_common")
# studies.append("eta_L2R")
# studies.append("eta_narrow")
#studies.append("eta_simple")

print("Running for: ", studies)
time.sleep(2)

# outdir = "DiJetJERC_DiJetHLT"
outdir = "DiJetJERC_DiJetHLT"
original_file = outdir+".xml"
original_dir_ = os.getcwd()


# QCDSamples = ["QCDPt","QCDHT", "DATA"]
# QCDSamples = ["QCDHT", "DATA"]
# QCDSamples = ["DATA"]
QCDSamples = ["QCDHT"]
processes = list(filter( lambda sample: year in sample and any(QCD in sample for QCD in QCDSamples) , QCD_process+Data_process))
others = list(set(QCD_process+Data_process)-set(processes))

JECVersions_Data = {}
JECVersions_MC = {}

JECVersions_Data["UL16preVFP"]  = ["Summer20UL16APV_V2"]
JECVersions_MC["UL16preVFP"]    = ["Summer20UL16APV_V2"]
JECVersions_Data["UL16postVFP"] = ["Summer20UL16_V2"]
JECVersions_MC["UL16postVFP"]   = ["Summer20UL16_V2"]
JECVersions_Data["UL17"]        = ["Summer20UL17_V2"]
JECVersions_MC["UL17"]          = ["Summer20UL17_V2"]
JECVersions_Data["UL18"]        = ["Summer20UL18_V2"]
JECVersions_MC["UL18"]          = ["Summer20UL18_V2"]

JECVersions_Data["2022"]        = ["Winter22Run3_V1"]
JECVersions_MC["2022"]          = ["Winter22Run3_V1"]
JECVersions_Data["2022preEE"]   = ["Winter22Run3_V1"]
JECVersions_MC["2022preEE"]     = ["Winter22Run3_V1"]
JECVersions_Data["2022postEE"]  = ["Summer22EEPrompt22_V1"] # V11 copy L2Res into L2L3Res due to nan values; V1 original
JECVersions_MC["2022postEE"]    = ["Summer22EEPrompt22_V1"] # V11 copy L2Res into L2L3Res due to nan values; V1 original

JECVersions_Data["2023"]        = ["Winter23Prompt23_V1"]
JECVersions_MC["2023"]          = ["Winter23Prompt23_V1"]

# JetLabels = ["AK4CHS", "AK4Puppi", "AK8CHS", "AK8Puppi"]
# JetLabels = ["AK4Puppi_v11"]
JetLabels = ["AK4Puppi"]
# JetLabels = ["AK4CHS", "AK4Puppi"]
# JetLabels = ["AK8Puppi", "AK4Puppi"]
# JetLabels = [ "AK8Puppi"]
# systematics = ["", "PU", "JEC", "JER"]
# systematics = ["", "PU", "JEC", "Prefire"]
# systematics = ["PU", "JEC", "Prefire", "PS"]
#systematics = ["PU", "JEC"]
# systematics = ["PS"]
systematics = [""]

print(year,studies, QCDSamples, JECVersions_Data[year], JetLabels, systematics)

for study in studies:

    userPathSframeOutput="/nfs/dust/cms/user/"+USER+"/sframe_all/"+outdir+"/"+year+"/"+study+"/"

    original_dir = original_dir_
    original_dir += "/SubmittedJobs/"+year+"/"+study+"/"

    main_program(option, internal_option, study, processes, others, JECVersions_Data[year], JECVersions_MC[year], JetLabels, systematics, original_dir, original_file, year)
