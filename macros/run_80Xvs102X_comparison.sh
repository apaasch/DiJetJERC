#!bash
# nevents=10000
nevents=1000000
# nevents=100000000
# root -l -b -q CompareTrees_80Xvs102X.C\(\"MET\",\"met_pt\",\"met_pt\",true,90,0,900,$nevents\)
# root -l -b -q CompareTrees_80Xvs102X.C\(\"AK4Jet1pt\",\"ljet_pt\",\"ak4jet1_pt\",true,90,0,900,$nevents\)
# root -l -b -q CompareTrees_80Xvs102X.C\(\"AK4Jet1eta\",\"ljet_eta\",\"ak4jet1_eta\",true,50,-2.5,2.5,$nevents\)
# root -l -b -q CompareTrees_80Xvs102X.C\(\"Muon1pt\",\"lep_pt\",\"lep1_pt\",true,90,0,900,$nevents\)
# root -l -b -q CompareTrees_80Xvs102X.C\(\"Muon1eta\",\"lep_eta\",\"lep1_eta\",true,50,-2.5,2.5,$nevents\)

# #root -l -b -q CompareTrees_80Xvs102X.C\(\"AK8Jet1pt\",\"ljet_pt\",\"ak8jet1_pt\",true,90,0,900,$nevents\)
# #root -l -b -q CompareTrees_80Xvs102X.C\(\"AK8Jet1eta\",\"ljet_eta\",\"ak8jet1_eta\",true,50,-2.5,2.5,$nevents\)

# root -l -b -q CompareTrees_80Xvs102X.C\(\"chi2\",\"rec_chi2\",\"rec_chi2\",true,200,0,100,$nevents\)
# root -l -b -q CompareTrees_80Xvs102X.C\(\"Mttbar\",\"Mttbar\",\"Mttbar\",true,300,0,6000,$nevents\)

# root -l -b -q CompareTrees_80Xvs102X.C\(\"MuonIDsf\",\"weight_sfmu_ID\",\"weight_sfmu_MuonID\",true,50,0,1.5,$nevents\)
# root -l -b -q CompareTrees_80Xvs102X.C\(\"MuonHLTsf\",\"weight_sfmu_HLT\",\"weight_sfmu_MuonTrigger\",true,50,0,1.5,$nevents\)


# root -l -b -q CompareTrees_80Xvs102X.C\(\"MuonHLTsf\",\"weight_sfmu_HLT\",\"weight_sfmu_MuonTrigger\",true,50,0,1.5,$nevents\)

# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunA\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunB\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunC\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunD\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &

# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunA\",\"pt_ave_asym_cut\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunA\",\"asymm\",\"asymmetry\",\"asymmetry\",true,100,-1,1,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunA\",\"asymm_cut\",\"asymmetry\",\"asymmetry\",true,100,-1,1,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunB\",\"asymm_cut\",\"asymmetry\",\"asymmetry\",true,100,-1,1,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunC\",\"asymm_cut\",\"asymmetry\",\"asymmetry\",true,100,-1,1,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunD\",\"asymm_cut\",\"asymmetry\",\"asymmetry\",true,100,-1,1,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunA\",\"asymm_abs\",\"asymmetry\",\"asymmetry\",true,100,0,1,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunB\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunC\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunD\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &


# root -l -b -q CompareTrees_80Xvs102X.C\(\"RunA\",\"asymm_cutHF\",\"asymmetry\",\"asymmetry\",true,100,-1,1,$nevents\) &
root -l -b -q CompareTrees_80Xvs102X.C\(\"RunA\",\"asymm_cutHF\",\"asymmetry\",\"asymmetry\",true,100,-1,1,-1\) &
root -l -b -q CompareTrees_80Xvs102X.C\(\"RunB\",\"asymm_cutHF\",\"asymmetry\",\"asymmetry\",true,100,-1,1,-1\) &
root -l -b -q CompareTrees_80Xvs102X.C\(\"RunC\",\"asymm_cutHF\",\"asymmetry\",\"asymmetry\",true,100,-1,1,-1\) &
root -l -b -q CompareTrees_80Xvs102X.C\(\"RunD\",\"asymm_cutHF\",\"asymmetry\",\"asymmetry\",true,100,-1,1,-1\) &

# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT100to200\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT1000to1500\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT1500to2000\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT2000toInf\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT200to300\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT300to500\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT500to700\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT50to100\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
# root -l -b -q CompareTrees_80Xvs102X.C\(\"QCDHT700to1000\",\"pt_ave\",\"pt_ave\",\"pt_ave\",true,100,0,1500,$nevents\) &
