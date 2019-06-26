For 2016 DATA
1a) hadd QCD Flat and QCD fwd to Full, e.g. uhh2.AnalysisModuleRunner.MC.QCDPt15to7000_pythia8_AK4CHS_Full_RunEFearly.root
1b) hadd QCD pT binned to uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_AK4CHS_RunH.root
2) hadd Data according to different periods (BCD, EFearly, FlateG, H)
3) Run main with "CalculateMCWeights()"
4) run BaconJetsTestModule_*.xml with "APPLY_WEIGHTS = true"
5a) Run main with "FullCycle_CorrectFormulae()"     (0 < |eta| < 5.2)
5b) Run main with "FullCycle_CorrectFormulae_eta()" (-5.2 < eta < 5.2)

How to compile the program:
source compile.sh
OR
source execute.sh (compile and execute)

$> g++ -std=c++0x -Wall -o main main.C CorrectionObject.cc ControlPlots.cc kFSR_CorrectFormulae.cc Pt_Extrapolation_Alternative_CorrectFormulae.cc L2ResOutput.cc L2ResAllRuns.cc InputForGlobalFit.cc FinalControlPlots_CorrectFormulae.cc useful_functions.cc CalculateMCWeights.cc tdrstyle_mod15.C `root-config  --cflags --evelibs`

$> g++ -std=c++0x -Wall -o main main.C CorrectionObject.cc ControlPlots.cc kFSR_CorrectFormulae.cc Pt_Extrapolation_Alternative_CorrectFormulae.cc L2ResOutput.cc L2ResAllRuns.cc InputForGlobalFit.cc FinalControlPlots_CorrectFormulae.cc useful_functions.cc CalculateMCWeights_TriggerThresholds.cc tdrstyle_mod15.C `root-config  --cflags --evelibs`
How to execute the program:
./main
