from utils import *

newNumber = {
    "DATA_RunB_UL17":        150,
    "DATA_RunC_UL17":        150,
    "DATA_RunD_UL17":        150,
    "DATA_RunE_UL17":        150,
    "DATA_RunF_UL17":        150,
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
    "QCDPt1800to2400_UL17":  110,
    "QCDPt2400to3200_UL17":  100,
    "QCDPt3200toInf_UL17":   50,
    "QCDHT50to100_UL17":     60,
    "QCDHT100to200_UL17":    50,
    "QCDHT200to300_UL17":    50,
    "QCDHT300to500_UL17":    50,
    "QCDHT500to700_UL17":    70,
    "QCDHT700to1000_UL17":   50,
    "QCDHT1000to1500_UL17":  65,
    "QCDHT1500to2000_UL17":  75,
    "QCDHT2000toInf_UL17":   70,

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

    "QCDHT50to100_UL18":     150,
    "QCDHT100to200_UL18":    130,
    "QCDHT200to300_UL18":    100,
    "QCDHT300to500_UL18":    100,
    "QCDHT500to700_UL18":    100,
    "QCDHT700to1000_UL18":   100,
    "QCDHT1000to1500_UL18":  100,
    "QCDHT1500to2000_UL18":  140,
    "QCDHT2000toInf_UL18":   100,
    "DATA_RunA_UL18":        220,
    "DATA_RunB_UL18":        200,
    "DATA_RunC_UL18":        200,
    "DATA_RunD_UL18":        70,

    "QCDHT50to100_UL16preVFP":      150,
    "QCDHT100to200_UL16preVFP":     130,
    "QCDHT200to300_UL16preVFP":     100,
    "QCDHT300to500_UL16preVFP":     100,
    "QCDHT500to700_UL16preVFP":     100,
    "QCDHT700to1000_UL16preVFP":    100,
    "QCDHT1000to1500_UL16preVFP":   100,
    "QCDHT1500to2000_UL16preVFP":   140,
    "QCDHT2000toInf_UL16preVFP":    100,
    "QCDHT50to100_UL16postVFP":     150,
    "QCDHT100to200_UL16postVFP":    130,
    "QCDHT200to300_UL16postVFP":    100,
    "QCDHT300to500_UL16postVFP":    100,
    "QCDHT500to700_UL16postVFP":    100,
    "QCDHT700to1000_UL16postVFP":   100,
    "QCDHT1000to1500_UL16postVFP":  100,
    "QCDHT1500to2000_UL16postVFP":  140,
    "QCDHT2000toInf_UL16postVFP":   100,
    "DATA_RunB_UL16preVFP":         150,
    "DATA_RunC_UL16preVFP":         150,
    "DATA_RunD_UL16preVFP":         150,
    "DATA_RunE_UL16preVFP":         150,
    "DATA_RunF_UL16preVFP":         150,
    "DATA_RunF_UL16postVFP":        150,
    "DATA_RunG_UL16postVFP":        150,
    "DATA_RunH_UL16postVFP":        150,
}


lumi_file = {
    "2017":         os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.root",
    "2018":         os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/2018/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.root",
    "UL16preVFP":   os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/UL16preVFP/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16preVFP_normtag.root",
    "UL16postVFP":  os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/UL16postVFP/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16postVFP_normtag.root",
    "UL17":         os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/UL17/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON_normtag.root",
    "UL18":         os.environ["CMSSW_BASE"]+"/src/UHH2/common/data/UL18/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_normtag.root",
}

TargetLumi = {
    "2016":         "35920",
    "2017":         "41530",
    "2018":         "59740",
    "UL16preVFP":   "35920",
    "UL16postVFP":  "35920",
    "UL17":         "41530",
    "UL18":         "59740",
    "RunII":       "137190",
}

@timeit
def createConfigFiles(study="Standard", processes=["QCDPt15to30", "QCDPt15to30_MB", "DATA_RunF"], others=[], JECVersions_Data=[""], JECVersions_MC=[""], JetLabels=["AK4CHS"], systematics=["PU", "JEC", "JER"], original_dir = "./submittedJobs/", original_file = "JER2018.xml", outdir="JER2018", year="2018", isMB = False, test_trigger=False, isThreshold=False, isLowPt=False, isL1Seed=False, isECAL=False,extratext=""):
    add_name = original_dir[original_dir.find("SubmittedJobs")+len("SubmittedJobs"):-1]
    try:
        time = int(filter(lambda x: "TIME" in x, open(original_dir[:original_dir.find("SubmittedJobs")]+original_file).readlines())[0].split("\"")[5])
    except:
        time = 3
    FileSplit = filter(lambda x: "FileSplit" in x, open(original_dir[:original_dir.find("SubmittedJobs")]+original_file).readlines())[0].split("\"")[3]
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
                changes.append(["Mail=", "USER@mail.desy.de", "USER@mail.desy.de", os.environ["USER"]+"@mail.desy.de"])
                changes.append(["<!ENTITY", "/nfs/dust/cms/user/USER", "USER", os.environ["USER"]])
                changes.append(["<!ENTITY", "CMSSW_BASE", "CMSSW_BASE", os.environ["CMSSW_BASE"]])
                change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
                changes = []
                changes.append(["user", "amalara", "amalara", os.environ["USER"]])
                change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
                changes = []
                changes.append(["<ConfigParse", 'FileSplit="'+FileSplit+'"', 'FileSplit="'+FileSplit+'"', 'FileSplit="'+str(int(newNumber[process]*0.9*time))+'"'])
                changes.append(["<!ENTITY", 'LUMI_FILE', 'lumifile.root', lumi_file[year]])
                changes.append(["<!ENTITY", "YEAR", "year", year])
                changes.append(["<Cycle", "TargetLumi", "defaultValue", TargetLumi[year]])
                changes.append(["<!ENTITY", "Study", "default", study])
                if study!= "Standard" and "L1" in study:
                    changes.append(["<!ENTITY", "JEC_LEVEL", "L1L2L3Residual", study])
                # if "17" in year:
                #     changes.append(["<!ENTITY", "PtBinsTrigger", '"DiJet"', '"SingleJet"'])
                # if "UL17" == year:
                #     changes.append(["<!ENTITY", "TRIGGER_FWD", '"true"', '"false"'])
                # if "18" in year:
                #     changes.append(["<!ENTITY", "APPLY_EtaPhi_HCAL", '"false"', '"true"'])
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
                            changes.append(["<!ENTITY", "jer_closuretest", '"false"', '"true"'])
                        elif "PU" in sys:
                            changes.append(["<!ENTITY", "SYSTYPE_PU", '"central"', '"'+dir+'"'])
                        change_lines(path+add_path_sys, newfilename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
