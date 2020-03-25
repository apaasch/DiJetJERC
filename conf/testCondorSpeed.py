from createConfigFiles import *
import glob, numpy as np

check = "*.o*"
rt = "Real time"
year = "UL17"
# year = "2018"
# version = "Autumn18_V19"
version = "Standard/Summer19UL17_V1_ComplexL1/"
JetLabels = "AK4CHS"

for Sample in newNumber:
    if not year in Sample: continue
    print Sample,
    path = os.environ["CMSSW_BASE"]+"/src/UHH2/DiJetJERC/conf/SubmittedJobs/"+year+"/"+version+"/"+JetLabels+"/workdir_DiJetJERC_DiJetHLT_*"+Sample+"*/**/"+check
    val = []
    nFilesMax = 0
    for el in glob.glob(path):
	nFiles = 0
        with open(el, "U") as file:
            lines = file.readlines()
            for line in lines:
		if "Input File" in line: nFiles +=1
                if not rt in line: continue
                sec = float(line.split()[10])
                val.append(sec)
	nFilesMax = nFilesMax if nFilesMax>nFiles else nFiles
    if len(val)==0 :
        print "\n"
        continue
    max = np.amax(np.array(val))
    min = np.amin(np.array(val))
    std = np.std(np.array(val))
    print " "*(30-len(Sample)), int(nFilesMax/round(max/3600,2)), nFilesMax, round(max/3600,2), "\t", round(min/3600,2), "\t", round(std*100/max,2), len(val), len(glob.glob(path))
