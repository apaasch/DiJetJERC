#!/usr/bin/env python3

from utils import *
import os, sys
sys.path.insert(1, os.environ["CMSSW_BASE"]+'/src/UHH2/common/UHH2-datasets')
from CrossSectionHelper import MCSampleValuesHelper

helper = MCSampleValuesHelper()

newNumber = {

    "QCD_Flat_2022preEE":                20,
    "QCDPt50to80_2022preEE":             10,
    "QCDPt80to120_2022preEE":            20,
    "QCDPt120to170_2022preEE":           20,
    "QCDPt170to300_2022preEE":           20,
    "QCDPt300to470_2022preEE":           15,
    "QCDPt470to600_2022preEE":           15,
    "QCDPt600to800_2022preEE":           20,
    "QCDPt800to1000_2022preEE":          15,
    "QCDPt1000to1400_2022preEE":         20,
    "QCDPt1400to1800_2022preEE":         20,
    "QCDPt1800to2400_2022preEE":         20,
    "QCDPt2400to3200_2022preEE":         20,
    "QCDPt3200toInf_2022preEE":          20,

    "QCDHT40to70_2022preEE":     150,
    "QCDHT70to100_2022preEE":    130,
    "QCDHT100to200_2022preEE":   100,
    "QCDHT200to400_2022preEE":   100,
    "QCDHT400to600_2022preEE":   100,
    "QCDHT600to800_2022preEE":   100,
    "QCDHT800to1000_2022preEE":  100,
    "QCDHT1000to1200_2022preEE": 140,
    "QCDHT1200to1500_2022preEE": 100,
    "QCDHT1500to2000_2022preEE": 100,
    "QCDHT2000_2022preEE":       100,

    # "QCDHT40to70_2022postEE":     150,
    "QCDHT70to100_2022postEE":    130,
    "QCDHT100to200_2022postEE":   100,
    "QCDHT200to400_2022postEE":   100,
    "QCDHT400to600_2022postEE":   100,
    "QCDHT600to800_2022postEE":   100,
    "QCDHT800to1000_2022postEE":  100,
    "QCDHT1000to1200_2022postEE": 140,
    "QCDHT1200to1500_2022postEE": 100,
    "QCDHT1500to2000_2022postEE": 100,
    "QCDHT2000_2022postEE":       100,

    # "DATA_RunC_2022preEE":         150,
    # "DATA_RunD_2022preEE":         150,

    "DATA_RunCHT_2022preEE":       150,
    "DATA_RunC_2022preEE":         150,
    "DATA_RunD_2022preEE":         150,
    "DATA_RunE_2022postEE":         150,
    "DATA_RunF_2022postEE":         150,
    "DATA_RunG_2022postEE":         150,

    "DATA_RunC_2023":         200,
    "DATA_RunC_v4_2023":       200,
    # "DATA_RunC0_v1_2023":         150,
    # "DATA_RunC0_v2_2023":         150,
    # "DATA_RunC0_v3_2023":         150,
    # "DATA_RunC1_v1_2023":         150,
    # "DATA_RunC1_v2_2023":         150,
    # "DATA_RunC1_v3_2023":         150,
}


lumi_file = {
    "2022":         os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2022/lumi_2022.root",
    "2022preEE":    os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2022/lumi_2022.root",
    "2022postEE":   os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2022/lumi_2022.root",
    "2023":         os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2023/lumi_runC.root",
}

TargetLumi = {
    "2022":         "35180", # TODO
    "2022preEE":    "35180", # TODO
    "2022postEE":   "35180", # TODO
    "2023":         "35180", # TODO
}

LumiWeight = {
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
}

@timeit
def createConfigFiles(study="Standard", processes=["QCDPt15to30", "QCDPt15to30_MB", "DATA_RunF"], others=[], JECVersions_Data=[""], JECVersions_MC=[""], JetLabels=["AK4CHS"], systematics=["PU", "JEC", "JER"], original_dir = "./submittedJobs/", original_file = "JER2018.xml", outdir="JER2018", year="2018", isMB = False, test_trigger=False, isThreshold=False, isLowPt=False, isECAL=False,extratext=""):
    add_name = original_dir[original_dir.find("SubmittedJobs")+len("SubmittedJobs"):-1]
    try:
        time = int(filter(lambda x: "TIME" in x, open(original_dir[:original_dir.find("SubmittedJobs")]+original_file).readlines())[0].split("\"")[5])
    except:
        time = 3
    FileSplit = list(filter(lambda x: "FileSplit" in x, open(original_dir[:original_dir.find("SubmittedJobs")]+original_file).readlines()))[0].split("\"")[3]
    isRun3 = "2022" in year or "2023" in year
    for index_JEC, newJECVersion in enumerate(JECVersions_Data):
        for newJetLabel in JetLabels:
            add_path = newJECVersion+"/"+newJetLabel+extratext+"/"
            path = original_dir+add_path
            if not os.path.exists(path):
                os.makedirs(path)
            for process in processes:
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
                if "QCD" in process:
                    proc = process.replace('_'+year, '').replace('HT','_HT').replace('postEE','').replace('preEE','').replace("QCDPt", "QCDPT")
                    lumi = helper.get_lumi(proc, '13TeV', year)
                    changes.append(["<InputData", "Type", "<CrossSection>", "{:.6f}".format(lumi)])
                    change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
                    changes = []
                # print(process)
                changes.append(["user", "amalara", "amalara", os.environ["USER"]])
                change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
                changes = []
                changes.append(["<ConfigParse", 'FileSplit="'+FileSplit+'"', 'FileSplit="'+FileSplit+'"', 'FileSplit="'+str(int(newNumber[process]*0.1*time))+'"']) # newNumber[process]*0.9*time
                changes.append(["<!ENTITY", 'LUMI_FILE', 'lumifile.root', lumi_file[year]])
                changes.append(["<!ENTITY", "YEAR", "year", year])
                changes.append(["<Cycle", "TargetLumi", "defaultValue", TargetLumi[year]])
                changes.append(["<!ENTITY", "Study", "default", study])
                if study!= "Standard" and "L1" in study:
                    changes.append(["<!ENTITY", "JEC_LEVEL", "L1L2L3Residual", study])
                # if 'DATA' in process and 'UL16' in year: # L2Residual not availble for UL16 in Summer20 yet -use from Summer19
                #     changes.append(["<!ENTITY", "JEC_LEVEL", "L1L2L3Residual", "L1L2"])
                # if 'DATA' in process and 'UL17' in year: # L2Residual not availble for UL17 in Summer20 yet -use from Summer19
                #     changes.append(["<!ENTITY", "JEC_LEVEL", "L1L2L3Residual", "L1L2"])
                # if 'DATA' in process and 'UL18' in year: # For first iteration
                #     if 'AK8' in newJetLabel:
                #         changes.append(["<!ENTITY", "JEC_LEVEL", "L1L2L3Residual", "L1L2"])
                #     else:
                #         changes.append(["<!ENTITY", "JEC_LEVEL", "L1L2L3Residual", "L1L2Residual"])
                # if "17" in year:
                #     changes.append(["<!ENTITY", "PtBinsTrigger", '"DiJet"', '"SingleJet"'])
                if "UL17" == year and ("RunB" in process or "RunC" in process):
                    changes.append(["<!ENTITY", "TRIGGER_FWD", '"true"', '"false"'])
                if "UL16" in year and "AK8" in newJetLabel:
                    changes.append(["<!ENTITY", "TRIGGER_FWD", '"true"', '"false"'])
                changes.append(["<!ENTITY", "OUTDIR", outdir , outdir+add_name+"/"+add_path])
                changes.append(["<ConfigSGE", "Workdir", "workdir_"+outdir, "workdir_"+outdir+"_"+process])
                if isThreshold:
                    changes.append(["<!ENTITY", "isThreshold", "false" , "true"])
                if "DATA" in process:
                    changes.append(["<!ENTITY", "JEC_VERSION", '"defaultJEC"', '"'+newJECVersion+'"'])
                if "QCD" in process:
                    changes.append(["<!ENTITY", "JEC_VERSION", '"defaultJEC"', '"'+JECVersions_MC[index_JEC]+'"'])
                if "AK4Puppi_" in newJetLabel:
                    changes.append(["<!ENTITY", "JETLABEL", '"AK4CHS"', '"AK4Puppi"'])
                else:
                    changes.append(["<!ENTITY", "JETLABEL", '"AK4CHS"', '"'+newJetLabel+'"'])
                if "Puppi" in newJetLabel or "AK8" in newJetLabel:
                    changes.append(["<!ENTITY", "APPLY_PUid_3rdjet", '"true"', '"false"'])
                if isRun3:
                    changes.append(["<Item", "ElectronCollection", '"slimmedElectronsUSER"', '"slimmedElectrons"'])
                    changes.append(["<Item", "MuonCollection", '"slimmedMuonsUSER"', '"slimmedMuons"'])
                    if "AK4Puppi" in newJetLabel:
                        changes.append(["<Item", "JetCollection", '"jetsAk4CHS"', '"slimmedJetsPuppi"'])
                    if "AK4CHS" in newJetLabel:
                        changes.append(["<Item", "JetCollection", '"jetsAk4CHS"', '"slimmedJets"'])
                else:
                    if "AK4Puppi" in newJetLabel:
                        if "v11" in newJetLabel:
                            changes.append(["<Item", "JetCollection", '"jetsAk4CHS"', '"patJetsAk4PuppiJetswithMultuplicity"'])
                        else:
                            changes.append(["<Item", "JetCollection", '"jetsAk4CHS"', '"jetsAk4Puppi"'])
                    if "AK8" in newJetLabel:
                        changes.append(["<!ENTITY", "PtBinsTrigger", '"DiJet"', '"SingleJet"'])
                        changes.append(["<Item", "GenJetCollection", '"slimmedGenJets"', '"slimmedGenJetsAK8"'])
                        if "Puppi" in newJetLabel: changes.append(["<Item", "JetCollection", '"jetsAk4CHS"', '"jetsAk8Puppi"'])
                        if "CHS" in newJetLabel: changes.append(["<Item", "JetCollection", '"jetsAk4CHS"', '"jetsAk8CHS"'])

                # if "QCD" in process:
                #     changes.append(["<!ENTITY", "PILEUP_DIRECTORY ", "MyMCPileupHistogram" , "MyMCPileupHistogram_"+process])
                change_lines(path, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])

                for sys in systematics:
                    if sys == "":
                        continue
                    dirs = ["", "up", "down"]
                    for dir in dirs:
                        if dir=="" and not "JER" in sys:
                            continue
                        if dir=="" and "JER" in sys:
                            dir = "nominal"
                        add_path_sys = sys+"/"+dir+"/"
                        if not os.path.exists(path+add_path_sys):
                            os.makedirs(path+add_path_sys)
                        print(path+add_path_sys)
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
                            changes.append(["<!ENTITY", "JERSMEAR_DIRECTION", '"nominal"', '"'+dir+'"'])
                            changes.append(["<!ENTITY", "jer_closuretest", '"false"', '"true"'])
                        elif "PU" in sys:
                            changes.append(["<!ENTITY", "SYSTYPE_PU", '"central"', '"'+dir+'"'])
                        elif "Prefire" in sys:
                            changes.append(["<!ENTITY", "PREFIRE_WEIGHT", '"nominal"', '"'+dir+'"'])

                        change_lines(path+add_path_sys, newfilename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
