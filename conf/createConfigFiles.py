from utils import *

@timeit
def createConfigFiles(processes=["QCDPt15to30", "QCDPt15to30_MB", "DATA_RunF"], JECVersions_Data=["Autumn18_V4"], JECVersions_MC=["Autumn18_V4"], JetLabels=["AK4CHS"], systematics=["PU", "JEC", "JER"], original_dir = "./submittedJobs/", original_file = "JER2018.xml", outdir="JER2018", isMB = False, test_trigger=False, isThreshold=False, isLowPt=False, isL1Seed=False, isECAL=False,extratext=""):
    add_name = original_dir[original_dir.find("SubmittedJobs")+len("SubmittedJobs"):-1]
    print add_name
    check_dir = add_name!="" or isMB or test_trigger or isThreshold or isLowPt or isL1Seed
    for index_JEC, newJECVersion in enumerate(JECVersions_Data):
        for newJetLabel in JetLabels:
            add_path = newJECVersion+"/"+newJetLabel+extratext+"/"
            path = original_dir+add_path
            if not os.path.exists(path):
                os.makedirs(path)
            for process in processes:
                print process
                if (isMB and not "_MB" in process) or (not isMB and "_MB" in process):
                    continue
                filename = original_file[:len(original_file)-4]+"_"+process+".xml"
                cmd = "cp %s %s" % (original_file, path+filename)
                a = os.system(cmd)
                cmd = "cp %s %s" % ("JobConfig.dtd", path)
                a = os.system(cmd)
                controls = []
                for el in list(set(processes) - set([process])):
                    if "QCD" in el:
                        controls.append(["<InputData", "Type", "MC", '"'+el+'"'])
                    elif "DATA" in el:
                        controls.append(["<InputData", "Type", "DATA", '"'+el+'"'])
                comment_lines(path, filename, controls, remove=True)
                if isMB:
                    controls = []
                    controls.append(["<Item", "Name", "trigger", "HLT_PFJet"])
                    comment_lines(path, filename, controls, remove=True)
                    if test_trigger:
                        controls = []
                        controls.append(["<Item", "Name", "triggerMB", "HLT_ZeroBias_v"])
                        comment_lines(path, filename, controls, remove=True)
                    else:
                        controls = []
                        controls.append(["<Item", "Name", "triggerMB_part", "HLT_ZeroBias_part"])
                        comment_lines(path, filename, controls, remove=True)
                else:
                    controls = []
                    controls.append(["<Item", "Name", "triggerMB", "Value"])
                    comment_lines(path, filename, controls, remove=True)
                if "QCD" in process:
                    sample = "QCD"
                elif "DATA" in process:
                    sample = "DATA"
                controls = []
                if check_dir:
                    controls.append(["<ConfigSGE", "Workdir", "workdir_"+outdir, "workdir_"+outdir+add_name+"_"+process])
                    controls.append(["<!ENTITY", "OUTDIR", outdir , outdir+add_name+"_"+sample+"/"+add_path])
                else:
                    controls.append(["<ConfigSGE", "Workdir", "workdir_"+outdir, "workdir_"+outdir+"_"+process])
                    controls.append(["<!ENTITY", "OUTDIR", outdir , outdir+"_"+sample+"/"+add_path])
                if isThreshold:
                    controls.append(["<!ENTITY", "isThreshold", "false" , "true"])
                if isL1Seed:
                    controls.append(["<!ENTITY", "apply_L1seed_from_bx1_filter", "false" , "true"])
                if "DATA" in process:
                    if "ECAL" in process:
                        controls.append(["<!ENTITY", "JEC_VERSION", '"Autumn18_V4"', '"Fall17_09May2018_V1"'])
                    else:
                        controls.append(["<!ENTITY", "JEC_VERSION", '"Autumn18_V4"', '"'+newJECVersion+'"'])
                if "QCD" in process:
                    controls.append(["<!ENTITY", "JEC_VERSION", '"Autumn18_V4"', '"'+JECVersions_MC[index_JEC]+'"'])
                controls.append(["<!ENTITY", "JETLABEL", '"AK4CHS"', '"'+newJetLabel+'"'])
                if "AK8" in newJetLabel:
                    controls.append(["<Item", "JetCollection", '"jetAk4CHS"', '"jetsAk8Puppi"'])
                    controls.append(["<Item", "GenJetCollection", '"slimmedGenJets"', '"slimmedGenJetsAK8"'])
                #     if isLowPt:
                #         controls.append(["<Item", "JetCollection", '"jetAk4CHS"', '"jetsAk8Puppi"'])
                #         controls.append(["<Item", "GenJetCollection", '"slimmedGenJets"', '"ak8GenJets"'])
                # else:
                #     if isLowPt:
                #         controls.append(["<Item", "JetCollection", '"jetAk4CHS"', '"patJetsAK4PFCHS"'])
                #         controls.append(["<Item", "GenJetCollection", '"slimmedGenJets"', '"ak4GenJets"'])
                if "QCD" in process:
                    controls.append(["<!ENTITY", "PILEUP_DIRECTORY ", "MyMCPileupHistogram" , "MyMCPileupHistogram_"+process])
                change_lines(path, filename, [el[0:2] for el in controls ], [el[2:3] for el in controls ], [el[3:4] for el in controls ])

                for sys in systematics:
                    if sys == "":
                        continue
                    for dir in ["up", "down"]:
                        if "JER" in sys:
                            dir = "nominal"
                        add_path_sys = sys+"/"+dir+"/"
                        if not os.path.exists(path+add_path_sys):
                            os.makedirs(path+add_path_sys)
                        newfilename = original_file[:len(original_file)-4]+"_"+process+"_"+sys+dir+".xml"
                        cmd = "cp %s %s" % (path+filename, path+add_path_sys+newfilename)
                        a = os.system(cmd)
                        cmd = "cp %s %s" % ("JobConfig.dtd", path+add_path_sys)
                        a = os.system(cmd)
                        controls = []
                        if check_dir:
                            controls.append(["<ConfigSGE", "Workdir", "workdir_"+outdir+add_name+"_"+process, "workdir_"+outdir+add_name+"_"+process+"_"+sys+dir])
                            controls.append(["<!ENTITY", "OUTDIR", outdir+add_name+"_"+sample+"/"+add_path , outdir+add_name+"_"+sample+"/"+add_path+add_path_sys])
                        else:
                            controls.append(["<ConfigSGE", "Workdir", "workdir_"+outdir+"_"+process, "workdir_"+outdir+"_"+process+"_"+sys+dir])
                            controls.append(["<!ENTITY", "OUTDIR", outdir+"_"+sample+"/"+add_path , outdir+"_"+sample+"/"+add_path+add_path_sys])
                        if "JEC" in sys:
                            controls.append(["<!ENTITY", "JECSMEAR_DIRECTION", '"nominal"', '"'+dir+'"'])
                        elif "JER" in sys:
                            controls.append(["<!ENTITY", "DO_JERSMEAR", '"false"', '"true"'])
                        elif "PU" in sys:
                            controls.append(["<!ENTITY", "SYSTYPE_PU", '"central"', '"'+dir+'"'])
                        change_lines(path+add_path_sys, newfilename, [el[0:2] for el in controls ], [el[2:3] for el in controls ], [el[3:4] for el in controls ])
