#!/usr/bin/env python3

from createConfigFiles import *

@timeit
def condor_control(original_dir ="./SubmittedJobs/" , JECVersions_Data=["Autumn18_V4"], JetLabels=["AK4CHS"], systematics=["", "PU", "JEC", "JER"], internal_option="-l", processes=[], extratext=""):
    count = 0
    list_processes = []
    nProcess = 48
    time_ = 1
    dirs_sys = ["", "up", "down"]
    dirs_PS = [p+d+'_'+f for p in ['FSR', 'ISR'] for d in ['up', 'down'] for f in ['sqrt2', '2', '4']]
    print(dirs_PS)
    for newJECVersion in JECVersions_Data:
        for newJetLabel in JetLabels:
            for sys in systematics:
                dirs = dirs_PS if "PS" in sys else dirs_sys
                for dir in dirs:
                    if sys == "" and dir != "":
                        continue
                    if sys == "JER" and dir != "":
                        continue
                    if sys == "JER" and dir == "":
                        dir = "nominal"
                    if 'PS' in sys:
                        if 'PS'!=sys and not dir in sys: # Is sys='PS' run all variations
                            continue
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
    print(len(list_processes))
    parallelise(list_processes, nProcess, cwd=True, time_=time_)


@timeit
def delete_workdir(original_dir ="./SubmittedJobs/" , JECVersions_Data=["Autumn18_V4", "Autumn18_V4"], JetLabels=["AK4CHS", "AK8Puppi"], systematics=["", "PU", "JEC", "JER"],extratext=""):
    add_name = original_dir[original_dir.find("SubmittedJobs")+len("SubmittedJobs"):-1]
    for sample in ["DATA", "QCD"]:
        for newJECVersion in JECVersions_Data:
            for newJetLabel in JetLabels:
                for sys in systematics:
                    for dir in ["", "up", "down"]:
                        if sys == "" and dir != "":
                            continue
                        if sys == "JER" and dir != "":
                            continue
                        if sys == "JER" and dir == "":
                       	    dir = "nominal"
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

QCD_process.append("QCDHT50to100_2018")
QCD_process.append("QCDHT100to200_2018")
QCD_process.append("QCDHT200to300_2018")
QCD_process.append("QCDHT300to500_2018")
QCD_process.append("QCDHT500to700_2018")
QCD_process.append("QCDHT700to1000_2018")
QCD_process.append("QCDHT1000to1500_2018")
QCD_process.append("QCDHT1500to2000_2018")
QCD_process.append("QCDHT2000toInf_2018")
Data_process.append("DATA_RunA_2018")
Data_process.append("DATA_RunB_2018")
Data_process.append("DATA_RunC_2018")
Data_process.append("DATA_RunD_2018")

QCD_process.append("QCDHT50to100_UL16preVFP")
QCD_process.append("QCDHT100to200_UL16preVFP")
QCD_process.append("QCDHT200to300_UL16preVFP")
QCD_process.append("QCDHT300to500_UL16preVFP")
QCD_process.append("QCDHT500to700_UL16preVFP")
QCD_process.append("QCDHT700to1000_UL16preVFP")
QCD_process.append("QCDHT1000to1500_UL16preVFP")
QCD_process.append("QCDHT1500to2000_UL16preVFP")
QCD_process.append("QCDHT2000toInf_UL16preVFP")
QCD_process.append("QCDHT50to100_UL16postVFP")
QCD_process.append("QCDHT100to200_UL16postVFP")
QCD_process.append("QCDHT200to300_UL16postVFP")
QCD_process.append("QCDHT300to500_UL16postVFP")
QCD_process.append("QCDHT500to700_UL16postVFP")
QCD_process.append("QCDHT700to1000_UL16postVFP")
QCD_process.append("QCDHT1000to1500_UL16postVFP")
QCD_process.append("QCDHT1500to2000_UL16postVFP")
QCD_process.append("QCDHT2000toInf_UL16postVFP")
Data_process.append("DATA_RunB_UL16preVFP")
Data_process.append("DATA_RunC_UL16preVFP")
Data_process.append("DATA_RunD_UL16preVFP")
Data_process.append("DATA_RunE_UL16preVFP")
Data_process.append("DATA_RunF_UL16preVFP")
Data_process.append("DATA_RunF_UL16postVFP")
Data_process.append("DATA_RunG_UL16postVFP")
Data_process.append("DATA_RunH_UL16postVFP")


QCD_process.append("QCDHT50to100_UL17")
QCD_process.append("QCDHT100to200_UL17")
QCD_process.append("QCDHT200to300_UL17")
QCD_process.append("QCDHT300to500_UL17")
QCD_process.append("QCDHT500to700_UL17")
QCD_process.append("QCDHT700to1000_UL17")
QCD_process.append("QCDHT1000to1500_UL17")
QCD_process.append("QCDHT1500to2000_UL17")
QCD_process.append("QCDHT2000toInf_UL17")
QCD_process.append("QCDPt15to30_UL17")
QCD_process.append("QCDPt30to50_UL17")
QCD_process.append("QCDPt50to80_UL17")
QCD_process.append("QCDPt80to120_UL17")
QCD_process.append("QCDPt120to170_UL17")
QCD_process.append("QCDPt170to300_UL17")
QCD_process.append("QCDPt300to470_UL17")
QCD_process.append("QCDPt470to600_UL17")
QCD_process.append("QCDPt600to800_UL17")
QCD_process.append("QCDPt800to1000_UL17")
QCD_process.append("QCDPt1000to1400_UL17")
QCD_process.append("QCDPt1400to1800_UL17")
QCD_process.append("QCDPt1800to2400_UL17")
QCD_process.append("QCDPt2400to3200_UL17")
QCD_process.append("QCDPt3200toInf_UL17")
Data_process.append("DATA_RunB_UL17")
Data_process.append("DATA_RunC_UL17")
Data_process.append("DATA_RunD_UL17")
Data_process.append("DATA_RunE_UL17")
Data_process.append("DATA_RunF_UL17")



QCD_process.append("QCDHT50to100_UL18")
QCD_process.append("QCDHT100to200_UL18")
QCD_process.append("QCDHT200to300_UL18")
QCD_process.append("QCDHT300to500_UL18")
QCD_process.append("QCDHT500to700_UL18")
QCD_process.append("QCDHT700to1000_UL18")
QCD_process.append("QCDHT1000to1500_UL18")
QCD_process.append("QCDHT1500to2000_UL18")
QCD_process.append("QCDHT2000toInf_UL18")
Data_process.append("DATA_RunA_UL18")
Data_process.append("DATA_RunB_UL18")
Data_process.append("DATA_RunC_UL18")
Data_process.append("DATA_RunD_UL18")


# JECVersions_Data = ["Autumn18_V4"]
# JetLabels = ["AK4CHS", "AK8Puppi"]
# systematics = ["", "PU", "JEC", "JER"]

# year = "2018"
#year = "UL16preVFP"
#year = "UL16postVFP"
# year = "UL17"
year = "UL18"


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
QCDSamples = ["QCDHT", "DATA"]
#QCDSamples = ["DATA"]
#QCDSamples = ["QCDHT"]
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

# JECVersions_Data["UL18"]        = ["Summer19UL18_V5"]
# JECVersions_MC["UL18"]          = ["Summer19UL18_V5"]

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
systematics = ["PS"]
# systematics = [""]

print(year,studies, QCDSamples, JECVersions_Data[year], JetLabels, systematics)

for study in studies:

    userPathSframeOutput="/nfs/dust/cms/user/"+USER+"/sframe_all/"+outdir+"/"+year+"/"+study+"/"

    original_dir = original_dir_
    original_dir += "/SubmittedJobs/"+year+"/"+study+"/"

    main_program(option, internal_option, study, processes, others, JECVersions_Data[year], JECVersions_MC[year], JetLabels, systematics, original_dir, original_file, year)
