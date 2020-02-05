from utils import *

newNumber = {
    "DATA_RunB_UL17":        100,
    "DATA_RunC_UL17":        100,
    "DATA_RunD_UL17":        100,
    "DATA_RunE_UL17":        100,
    "DATA_RunF_UL17":        100,
    "QCDPt15to30_UL17":      50,
    "QCDPt30to50_UL17":      60,
    "QCDPt50to80_UL17":      40,
    "QCDPt80to120_UL17":     30,
    "QCDPt120to170_UL17":    40,
    "QCDPt170to300_UL17":    40,
    "QCDPt300to470_UL17":    40,
    "QCDPt470to600_UL17":    50,
    "QCDPt600to800_UL17":    50,
    "QCDPt800to1000_UL17":   50,
    "QCDPt1000to1400_UL17":  50,
    "QCDPt1400to1800_UL17":  90,
    "QCDPt1800to2400_UL17":  100,
    "QCDPt2400to3200_UL17":  100,
    "QCDPt3200toInf_UL17":   150,

    "QCDHT50to100_2018":     125,
    "QCDHT100to200_2018":    105,
    "QCDHT200to300_2018":    90,
    "QCDHT300to500_2018":    90,
    "QCDHT500to700_2018":    80,
    "QCDHT700to1000_2018":   85,
    "QCDHT1000to1500_2018":  85,
    "QCDHT1500to2000_2018":  80,
    "QCDHT2000toInf_2018":   75,
    "DATA_RunA_2018":        155,
    "DATA_RunB_2018":        300,
    "DATA_RunC_2018":        270,
    "DATA_RunD_2018":        180,

}



lumi_file = {
    "2017": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.root",
    "2018": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2018/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.root",
    "UL17": os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.root",
}

@timeit
def createConfigFiles(processes=["QCDPt15to30", "QCDPt15to30_MB", "DATA_RunF"], others=[], JECVersions_Data=[""], JECVersions_MC=[""], JetLabels=["AK4CHS"], systematics=["PU", "JEC", "JER"], original_dir = "./submittedJobs/", original_file = "JER2018.xml", outdir="JER2018", year="2018", isMB = False, test_trigger=False, isThreshold=False, isLowPt=False, isL1Seed=False, isECAL=False,extratext=""):
    add_name = original_dir[original_dir.find("SubmittedJobs")+len("SubmittedJobs"):-1]
    time = int(filter(lambda x: "TIME" in x, open(original_dir[:original_dir.find("SubmittedJobs")]+original_file).readlines())[0].split("\"")[5])
    FileSplit = filter(lambda x: "FileSplit" in x, open(original_dir[:original_dir.find("SubmittedJobs")]+original_file).readlines())[0].split("\"")[3]
    print add_name
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
                comments = []
                for el in list(set(others+processes) - set([process])):
                    if "QCD" in el:
                        comments.append(["<InputData", "Type", "MC", '"'+el+'"'])
                    elif "DATA" in el:
                        comments.append(["<InputData", "Type", "DATA", '"'+el+'"'])
                comment_lines(path, filename, comments, remove=True)
                changes = []
                changes.append(["user", "amalara/WorkingArea", "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/", os.environ["CMSSW_BASE"]+"/"])
                change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
                changes = []
                changes.append(["user", "amalara", "amalara", os.environ["USER"]])
                change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
                changes = []
                changes.append(["<ConfigParse", 'FileSplit="'+FileSplit+'"', 'FileSplit="'+FileSplit+'"', 'FileSplit="'+str(int(newNumber[process]*time))+'"'])
                changes.append(["<!ENTITY", 'LUMI_FILE', 'lumifile.root', lumi_file[year]])
                changes.append(["<!ENTITY", "YEAR", "year", year])
                if "17" in year:
                    changes.append(["<Cycle", "TargetLumi", "158640", "41530"])
                    changes.append(["<!ENTITY", "PtBinsTrigger", '"DiJet"', '"SingleJet"'])
                if "UL17" == year:
                    changes.append(["<!ENTITY", "TRIGGER_FWD", '"true"', '"false"'])
                if "2018" == year:
                    changes.append(["<Cycle", "TargetLumi", "158640", "59740"])
                changes.append(["<!ENTITY", "OUTDIR", outdir , outdir+add_name+"/"+add_path])
                changes.append(["<ConfigSGE", "Workdir", "workdir_"+outdir, "workdir_"+outdir+"_"+process])
                if isThreshold:
                    changes.append(["<!ENTITY", "isThreshold", "false" , "true"])
                if isL1Seed:
                    changes.append(["<!ENTITY", "apply_L1seed_from_bx1_filter", "false" , "true"])
                if "DATA" in process:
                    changes.append(["<!ENTITY", "JEC_VERSION", '"defaultJEC"', '"'+newJECVersion+'"'])
                if "QCD" in process:
                    changes.append(["<!ENTITY", "JEC_VERSION", '"defaultJEC"', '"'+JECVersions_MC[index_JEC]+'"'])
                changes.append(["<!ENTITY", "JETLABEL", '"AK4CHS"', '"'+newJetLabel+'"'])
                if "AK4Puppi" in newJetLabel:
                    changes.append(["<Item", "JetCollection", '"jetsAk4CHS"', '"jetsAk4Puppi"'])
                if "AK8" in newJetLabel:
                    changes.append(["<!ENTITY", "PtBinsTrigger", '"DiJet"', '"SingleJet"'])
                    changes.append(["<Item", "JetCollection", '"jetsAk4CHS"', '"jetsAk8Puppi"'])
                    changes.append(["<Item", "GenJetCollection", '"slimmedGenJets"', '"slimmedGenJetsAK8"'])
                # if "QCD" in process:
                #     changes.append(["<!ENTITY", "PILEUP_DIRECTORY ", "MyMCPileupHistogram" , "MyMCPileupHistogram_"+process])
                change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])

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
                        changes = []
                        changes.append(["<!ENTITY", "OUTDIR", outdir+add_name+"/"+add_path , outdir+add_name+"/"+add_path+add_path_sys])
                        changes.append(["<ConfigSGE", "Workdir", "workdir_"+outdir+"_"+process, "workdir_"+outdir+"_"+process+"_"+sys+dir])
                        if "JEC" in sys:
                            changes.append(["<!ENTITY", "JECSMEAR_DIRECTION", '"nominal"', '"'+dir+'"'])
                        elif "JER" in sys:
                            changes.append(["<!ENTITY", "DO_JERSMEAR", '"false"', '"true"'])
                        elif "PU" in sys:
                            changes.append(["<!ENTITY", "SYSTYPE_PU", '"central"', '"'+dir+'"'])
                        change_lines(path+add_path_sys, newfilename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
