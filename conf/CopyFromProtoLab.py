import os, glob

newDir = "../../JECDatabase/textFiles/Summer19UL18"
newVersion = "V2"

for run in ["A","B","C","D"]:
    folder = "JERCProtoLab/Summer19UL18/L1Offset/V3/Run"+run
    dir = newDir+"_Run"+run+"_V2_DATA"
    os.system("mkdir -p "+dir)
    for f_ in glob.glob(folder+"/*"):
        os.system("cp "+f_+" "+f_.replace(folder,dir).replace("V3",newVersion))
    for f_ in glob.glob(newDir+"_V2_MC/*L2Relative*"):
        os.system("cp "+f_+" "+f_.replace("V2","Run"+run+"_"+newVersion).replace("MC","DATA"))
    for f_ in glob.glob(newDir+"_V2_MC/*L3Absolute*"):
        os.system("cp "+f_+" "+f_.replace("V2","Run"+run+"_"+newVersion).replace("MC","DATA"))
    folder = "JERCProtoLab/Summer19UL18/L2Residual/V1/Run"+run
    for f_ in glob.glob(folder+"/*"):
        if not "MPF_LOGLIN" in f_: continue
        os.system("cp "+f_+" "+f_.replace(folder,dir).replace("V1_MPF_LOGLIN","Run"+run+"_"+newVersion+"_DATA").replace("_pythia8",""))
        os.system("cp "+f_+" "+f_.replace(folder,dir).replace("V1_MPF_LOGLIN","Run"+run+"_"+newVersion+"_DATA").replace("_pythia8","").replace("L2Residual","L2L3Residual"))
    folder = newDir.replace("Summer19UL18","Summer19UL17_V5_MC")
    for f_ in glob.glob(folder+"/Summer19UL17_V5_MC_Uncertainty*"):
        os.system("cp "+f_+" "+f_.replace("Summer19UL17","Summer19UL18").replace("V5",newVersion))
        os.system("cp "+f_+" "+f_.replace("Summer19UL17","Summer19UL18").replace("V5_MC","Run"+run+"_"+newVersion+"_DATA"))

'''

- L1RC is different for data and MC.
- L1RC DataMcSF is different for each run.
- L1FastJet is different for data and MC.
- L2Relative is equal for data and MC. Derived from MC. It's MC Response. Include L2 and L3 steps.
- L3Absolute is dummy for data and MC.
- L2Residual is dummy for MC.
- L2L3Residual is dummy for MC.
- L2L3Residual is from L2Residual + L3Residual for pt dep.
- Uncertainty and UncertaintySources are equal for data and MC.

'''
