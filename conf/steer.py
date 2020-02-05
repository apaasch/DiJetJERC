#!/usr/bin/env python

from createConfigFiles import *

def cont_event(paths ="./submittedJobs/" , JECVersions_Data=["Autumn18_V4"], JetLabels=["AK4CHS"], systematics=["", "PU", "JEC", "JER"],extratext=""):
    count = 0
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
                    path = paths+newJECVersion+"/"+newJetLabel+extratext+"/"+sys+"/"+dir+"/"
                    for sample in sorted(os.listdir(path)):
                        if not ".xml" in sample:
                            continue
                        count += 1
    return count

@timeit
def condor_control(original_dir ="./SubmittedJobs/" , JECVersions_Data=["Autumn18_V4"], JetLabels=["AK4CHS"], systematics=["", "PU", "JEC", "JER"], internal_option="-l", processes=[], extratext=""):
    count = 0
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
                    path = original_dir+newJECVersion+"/"+newJetLabel+extratext+"/"+sys+"/"+dir+"/"
                    for sample in sorted(os.listdir(path)):
                        if not ".xml" in sample:
                            continue
 			if all(not control in sample for control in processes): continue
                        count += 1
                        all_events = cont_event(original_dir, JECVersions_Data, JetLabels, systematics,extratext)
                        print "Already completed "+str(count)+" out of "+str(all_events)+" jobs --> "+str(float(count)/float(all_events)*100)+"%."
                        os.chdir(original_dir)
                        os.chdir(path)
                        if internal_option:
                            command = ['sframe_batch.py', internal_option, sample]
                        else:
                            command = ['sframe_batch.py', sample]
                        process = subprocess.Popen(command)
                        process.wait()
                        if internal_option == "-a":
                            time.sleep(5)
                        os.chdir(original_dir)


from createConfigFiles import *
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
                                    print cmd
                        path = original_dir+newJECVersion+"/"+newJetLabel+extratext+"/"+sys+"/"+dir+"/"
                        if os.path.isdir(path):
                            for workdir in sorted(os.listdir(path)):
                                if "workdir" in workdir:
                                    cmd = "rm -fr %s" % (path+workdir)
                                    a = os.system(cmd)




def main_program(option="", internal_option="", processes=[], others=[], JECVersions_Data=[], JECVersions_MC=[], JetLabels=[], systematics=[], original_dir="./SubmittedJobs/", original_file="JER2018.xml", year="2018", isMB=False, test_trigger=False, isThreshold=False, isLowPt=False, isL1Seed=False, isECAL=False, extratext=""):
    if option == "new":
        createConfigFiles(processes, others, JECVersions_Data, JECVersions_MC, JetLabels, systematics, original_dir, original_file, outdir, year, isMB, test_trigger, isThreshold,isLowPt,isL1Seed,isECAL,extratext)
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

QCD_process.append("QCDHT50to100_UL17")
QCD_process.append("QCDHT100to200_UL17")
QCD_process.append("QCDHT200to300_UL17")
QCD_process.append("QCDHT300to500_UL17")
QCD_process.append("QCDHT500to700_UL17")
QCD_process.append("QCDHT700to1000_UL17")
QCD_process.append("QCDHT1000to1500_UL17")
QCD_process.append("QCDHT1500to2000_UL17")
QCD_process.append("QCDHT2000toInf_UL17")

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


# JECVersions_Data = ["Autumn18_V4"]
# JetLabels = ["AK4CHS", "AK8Puppi"]
# systematics = ["", "PU", "JEC", "JER"]

#year = "UL17"
year = "2018"
outdir = "DiJetJERC_DiJetHLT"

userPathSframeOutput="/nfs/dust/cms/user/"+USER+"/sframe_all/"+outdir+"/"+year+"/"
original_file = outdir+".xml"
original_dir_ = os.getcwd()

original_dir = original_dir_
original_dir += "/SubmittedJobs/"+year+"/"

QCDSamples = ["QCDHT","DATA"]
processes = filter( lambda sample: year in sample and any(QCD in sample for QCD in QCDSamples) , QCD_process+Data_process)
others = list(set(QCD_process+Data_process)-set(processes))

JECVersions_Data = ["Fall17_17Nov2017_V32"]
JECVersions_MC   = ["Fall17_17Nov2017_V32"]
JECVersions_Data = ["Autumn18_V19"]
JECVersions_MC   = ["Autumn18_V19"]
# JetLabels = ["AK4CHS","AK8Puppi", "AK4Puppi"]
JetLabels = ["AK4CHS"]
systematics = ["", "PU", "JEC"]
#systematics = [""]

main_program(option, internal_option, processes, others, JECVersions_Data, JECVersions_MC, JetLabels, systematics, original_dir, original_file, year)
