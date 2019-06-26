#!/usr/bin/env python

from createConfigFiles import *

def cont_event(paths ="./submittedJobs/" , JECVersions_Data=["Autumn18_V4"], JetLabels=["AK4CHS"], systematics=["", "PU", "JEC", "JER"]):
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
                    path = paths+newJECVersion+"/"+newJetLabel+"/"+sys+"/"+dir+"/"
                    for sample in sorted(os.listdir(path)):
                        if not ".xml" in sample:
                            continue
                        count += 1
    return count

@timeit
def condor_control(original_dir ="./SubmittedJobs/" , JECVersions_Data=["Autumn18_V4"], JetLabels=["AK4CHS"], systematics=["", "PU", "JEC", "JER"], internal_option="-l", processes=[]):
#    print "path = ",original_dir
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
                    path = original_dir+newJECVersion+"/"+newJetLabel+"/"+sys+"/"+dir+"/"
                    #print "full path", path
                    for sample in sorted(os.listdir(path)):
                        if not ".xml" in sample:
                            continue
 			if all(not control in sample for control in processes): continue
                        count += 1
                        all_events = cont_event(original_dir, JECVersions_Data, JetLabels, systematics)
                        print "Already completed "+str(count)+" out of "+str(all_events)+" jobs --> "+str(float(count)/float(all_events)*100)+"%."
                        os.chdir(original_dir)
                        os.chdir(path)
                        if internal_option:
                            print "sample ",sample
                            print "internal_option", internal_option
                            command = ['sframe_batch.py', internal_option, sample]
                        else:
                            command = ['sframe_batch.py', sample]
                        #print "command",command
                        process = subprocess.Popen(command)
                        process.wait()
                        if internal_option == "-a":
                            time.sleep(5)
                        os.chdir(original_dir)


from createConfigFiles import *
@timeit
def delete_workdir(original_dir ="./SubmittedJobs/" , JECVersions_Data=["Autumn18_V4", "Autumn18_V4"], JetLabels=["AK4CHS", "AK8PUPPI"], systematics=["", "PU", "JEC", "JER"]):
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
                        path = "/nfs/dust/cms/user/karavdia/DiJetJERC_preselection/"+outdir+add_name+"_"+sample+"/"+newJECVersion+"/"+newJetLabel+"/"+sys+"/"+dir+"/"
                        if os.path.isdir(path):
                            for workdir in sorted(os.listdir(path)):
                                if "workdir" in workdir:
                                    cmd = "rm -fr %s" % (path+workdir)
                                    a = os.system(cmd)
                        path = original_dir+newJECVersion+"/"+newJetLabel+"/"+sys+"/"+dir+"/"
                        if os.path.isdir(path):
                            for workdir in sorted(os.listdir(path)):
                                if "workdir" in workdir:
                                    cmd = "rm -fr %s" % (path+workdir)
                                    a = os.system(cmd)




def main_program(option="", internal_option="", processes=[], JECVersions_Data=[], JECVersions_MC=[], JetLabels=[], systematics=[], original_dir="./SubmittedJobs/", original_file="JER2018.xml", isMB=False, test_trigger=False, isThreshold=False, isLowPt=False, isL1Seed=False, isECAL=False):
    if option == "new":
        createConfigFiles(processes, JECVersions_Data, JECVersions_MC, JetLabels, systematics, original_dir, original_file, outdir, isMB, test_trigger, isThreshold,isLowPt,isL1Seed,isECAL)
    elif option == "remove" or option == "delete":
        delete_workdir(original_dir, JECVersions_Data, JetLabels, systematics)
    else:
        condor_control(original_dir, JECVersions_Data, JetLabels, systematics, internal_option, processes)



##################################################
#                                                #
#                   MAIN Program                 #
#                                                #
##################################################

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
QCD_process.append("QCD_Flat2018")
QCD_process.append("QCD_Flat")
# # QCD_process.append("QCDPt15to30")
# # QCD_process.append("QCDPt30to50")
# # QCD_process.append("QCDPt50to80")
# # QCD_process.append("QCDPt80to120")
# # QCD_process.append("QCDPt120to170")
# # QCD_process.append("QCDPt170to300")
# # QCD_process.append("QCDPt300to470")
# # QCD_process.append("QCDPt470to600")
# # QCD_process.append("QCDPt600to800")
# # QCD_process.append("QCDPt800to1000")
# # QCD_process.append("QCDPt1000to1400")
# # QCD_process.append("QCDPt1400to1800")
# # QCD_process.append("QCDPt1800to2400")
# # QCD_process.append("QCDPt2400to3200")
# # QCD_process.append("QCDPt3200toInf")
#
QCD_process.append("QCDHT50to100")
QCD_process.append("QCDHT100to200")
QCD_process.append("QCDHT200to300")
QCD_process.append("QCDHT300to500")
QCD_process.append("QCDHT500to700")
QCD_process.append("QCDHT700to1000")
QCD_process.append("QCDHT1000to1500")
QCD_process.append("QCDHT1500to2000")
QCD_process.append("QCDHT2000toInf")


Data_process= []
Data_process.append("DATA_RunA")
Data_process.append("DATA_RunB")
Data_process.append("DATA_RunC")
Data_process.append("DATA_RunD")

processes = QCD_process+Data_process
#processes = Data_process

# JECVersions_Data = ["Autumn18_V4"]
# JetLabels = ["AK4CHS", "AK8PUPPI"]
# systematics = ["PU", "JEC", "JER"]

original_file = "DiJetJERC_DiJetHLT.xml"
#outdir = "JER2018_TEST2"
#outdir = "JER2018_woPUid"
#outdir = "JER2018_wPUid"
outdir = "DiJetJERC_JERSFtest"
original_dir_ = os.getcwd()

#JECVersions_Data = ["Autumn18_V10"]
#JECVersions_MC = ["Autumn18_V10"]

JECVersions_Data = ["Autumn18_V13h"]
JECVersions_MC = ["Autumn18_V13h"]

JetLabels = ["AK4CHS"]
# JetLabels = ["AK8PUPPI"]
#systematics = ["", "PU", "JEC", "JER"]
#systematics = ["PU", "JEC", "JER"]
#systematics = ["", "PU"]
systematics = [""]

isLowPt = False
isMB = False
test_trigger = False
isThreshold = False
isL1Seed = False
isECAL = False
original_dir = original_dir_
original_dir += "/SubmittedJobs/"
main_program(option, internal_option, processes, JECVersions_Data, JECVersions_MC, JetLabels, systematics, original_dir, original_file, isMB, test_trigger, isThreshold, isLowPt, isL1Seed)
