import glob, os, sys, ROOT

# ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT.kError) + ";")
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT. kFatal) + ";")

class CheckFileNumbers():
    def __init__(self, year="", study="", jecversion = "", jetLabel="", systematic=""):
        USER = os.environ["USER"]
        # self.sframeDir = "/nfs/dust/cms/user/"+USER+"/sframe_all//DiJetJERC_DiJetHLT/"
        self.sframeDir = "/nfs/dust/cms/user/"+USER+"/sframe_all/DiJetJERC_DiJetHLT/"
        self.xmlDir = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/conf/SubmittedJobs/"
        self.systematic = systematic.replace("_","")
        if systematic == "":
            self.folder = year+"/"+study+"/"+jecversion+"/"+jetLabel+"/"
        else:
            self.folder = year+"/"+study+"/"+jecversion+"/"+jetLabel+"/"+systematic.split("_")[0]+"/"+systematic.split("_")[1]+"/"
        # self.module = "DiJetJERC_DiJetHLT"
        self.module = "DiJetJERC_DiJetHLT"
        self.PrefixrootFile = "uhh2.AnalysisModuleRunner."

    def Count(self):
        config_path  = self.xmlDir+self.folder
        storage_path = self.sframeDir+self.folder
        # print config_path
        tot_xml = 0
        tot_root = 0
        for workdir in glob.glob(config_path+"workdir_"+self.module+"_*"):
            file_dir = workdir.replace(config_path,storage_path)
            workdir_name = workdir.replace(config_path,"").replace("workdir_"+self.module+"_","").replace(self.systematic,"")
            # print "\t", file_dir
            # print "\t", workdir_name
            for xml in glob.glob(workdir+"/*xml"):
                xml_name = xml.replace(workdir,"")
                if "Result" in xml_name: continue
                if self.module in xml_name: continue
                tot_xml += 1
                number = xml_name.replace(workdir_name,"").replace(".xml","").replace("/","").replace("_","")
                # print xml_name, number
                number = int(number)
                root_file = file_dir+"/"+self.PrefixrootFile+("MC." if "QCD" in xml_name else "DATA.")+workdir_name+"_"+str(number-1)+".root"
                root_file = root_file.replace("__","_")
                # print number, xml, root_file
                if not os.path.isfile(root_file):
                    dummy=0
                    # print "FILE doesn't exist", root_file
                    # print "mkdir -p", file_dir.replace(self.xmlDir,self.sframeDir), "; sframe_main",xml
		else:
                    ntuple = ROOT.TFile(str(root_file))
                    if ntuple.IsZombie() or ntuple.ReadKeys()==0:
                        # print "sframe_main",xml
			dummy=1
                    else: tot_root += 1
        print self.folder, " "*(70-len(self.folder)) , "DONE. Counted: XML=", tot_xml, " "*(5-len(str(tot_xml))), "ROOT:", tot_root


if __name__ == '__main__':

    # years       = ["UL16postVFP", "UL17", "UL18"]
    # years       = ["UL17", "UL18"]
    # studies     = ["eta_common"]
    # JetLabels   = ["AK4CHS", "AK8CHS", "AK8Puppi"]

    years       = ["UL18"]
    studies     = ["eta_common"]
    # JetLabels   = ["AK4Puppi", "AK4CHS", "AK8CHS", "AK8Puppi"]
    JetLabels   = ["AK4Puppi"]

    # Systematics = ["", "PU_up", "PU_down", "JEC_up", "JEC_down", "JER_nominal"]
    # Systematics = ["", "PU_up", "PU_down", "JEC_up", "JEC_down"]
    Systematics = [""]

    JECVersions = {}
    JECVersions["UL16preVFP"]    = ["Summer20UL16APV_V2"]
    JECVersions["UL16postVFP"]   = ["Summer20UL16_V2"]
    JECVersions["UL17"]          = ["Summer20UL17_V2"]
    JECVersions["UL18"]          = ["Summer20UL18_V2"]

    for year in years:
        for study in studies:
            for jecversion in JECVersions[year]:
                for jetLabel in JetLabels:
                    for systematic in Systematics:
                        CFN = CheckFileNumbers(year=year, study=study, jecversion=jecversion, jetLabel=jetLabel, systematic=systematic)
                        CFN.Count()
