//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 13 16:05:52 2015 by ROOT version 5.34/32
// from TTree t/t
// found on file: /pnfs/desy.de/cms/tier2/store/user/mniedzie/Spring15_QCD_Pt-binned/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/JERmcSpring15/150810_142114/0000/output_mc_1.root
//////////////////////////////////////////////////////////

#ifndef MySelector_h
#define MySelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TDirectory.h>
#include <vector>
#include <algorithm>

#include "constants.h"


double Weight( TString filename );
bool JetInRange(double jet_eta, double min, double max);
bool JetInEtaBin(double jet_eta, std::vector<double> bins, int bin);

class MySelector : public TSelector {
public:

  TString outdir;
  std::string year, study, binning, abins;
  bool isAK8;

  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  // Declaration of leaf types
  float weight;
  float weight_pu;
  float weight_pu_down;
  float weight_pu_up;
  float nPU;
  int njet;
  bool is_JER_SM;
  float pt_ave;
  float barreljet_phi;
  float barreljet_eta;
  float barreljet_pt;
  float probejet_phi;
  float probejet_eta;
  float probejet_pt;
  float jet3_pt;
  float asymmetry;
  float rho;
  int HLT;
  int L1min;
  int L1max;
  float alpha_;
  float eventID;
  int run;
  int lumi_sec;

  TBranch *b_weight;
  TBranch *b_weight_pu;
  TBranch *b_weight_pu_down;
  TBranch *b_weight_pu_up;
  TBranch *b_nPU;
  TBranch *b_njet;
  TBranch *b_is_JER_SM;
  TBranch *b_pt_ave;
  TBranch *b_barreljet_phi;
  TBranch *b_barreljet_eta;
  TBranch *b_barreljet_pt;
  TBranch *b_probejet_phi;
  TBranch *b_probejet_eta;
  TBranch *b_probejet_pt;
  TBranch *b_jet3_pt;
  TBranch *b_asymmetry;
  TBranch *b_alpha;
  TBranch *b_rho;
  TBranch *b_prescale;
  TBranch *b_L1min;
  TBranch *b_L1max;

  TBranch *b_eventID;
  TBranch *b_run;
  TBranch *b_lumi_sec;


  MySelector(TString name, std::string year_, std::string study_, std::string binning_, std::string abins_, bool isAK8_, TTree * /*tree*/ =0) : fChain(0), outdir(name), year(year_), study(study_), binning(binning_), abins(abins_), isAK8(isAK8_) { }
  virtual ~MySelector() { }
  virtual int   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual bool  Notify();
  virtual bool  Process(Long64_t entry);
  virtual int   GetEntry(Long64_t entry, int getall = 0) {
    int treeentry = entry;
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  void BuildEvent();

  int TotalEvents;

  int EtaBins_SM, EtaBins_SM_control, EtaBins_FE_reference, EtaBins_FE_control, EtaBins_FE, EtaBins;
  int etaShift_SM, etaShift_SM_control, etaShift_FE_reference, etaShift_FE_control, etaShift_FE;
  int PtBins_Central, PtBins_HF, PtBins;
  int AlphaBins, AlphaBinsInc;

  std::vector<int> Pt_bins_Central;
  std::vector<int> Pt_bins_HF;
  std::vector<double> Eta_bins_SM;
  std::vector<double> Eta_bins_SM_control;
  std::vector<double> Eta_bins_FE_reference;
  std::vector<double> Eta_bins_FE_control;
  std::vector<double> Eta_bins_FE;
  std::vector<double> Alpha_bins, Alpha_bins_Inc;



  std::vector< std::vector< std::vector< double > > > nevents_central,nevents_HF;

  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_SM, 						asymmetries_pt_SM,						asymmetries_rho_SM,						asymmetries_pt3_SM;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_SM_control, 		asymmetries_pt_SM_control,		asymmetries_rho_SM_control,		asymmetries_pt3_SM_control;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_FE_reference, 	asymmetries_pt_FE_reference,	asymmetries_rho_FE_reference,	asymmetries_pt3_FE_reference;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_FE_control, 		asymmetries_pt_FE_control,		asymmetries_rho_FE_control,		asymmetries_pt3_FE_control;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_FE, 						asymmetries_pt_FE,						asymmetries_rho_FE,						asymmetries_pt3_FE;
  std::vector< std::vector< std::vector< TH1F* > > > pt_spectrum_SM, pt_spectrum_SM_control, pt_spectrum_FE_reference, pt_spectrum_FE_control, pt_spectrum_FE;

  std::vector< std::vector< TH1F* > > alpha_spectrum_SM, alpha_spectrum_SM_control, alpha_spectrum_FE_reference, alpha_spectrum_FE_control, alpha_spectrum_FE;

  std::vector< std::vector< TH1F* > > inclusive_rho_Central, inclusive_ptave_Central, inclusive_nevents_Central; // alpha, pt
  std::vector< std::vector< TH1F* > > inclusive_rho_HF, inclusive_ptave_HF, inclusive_nevents_HF; // alpha, pt
  std::vector< std::vector< TH1F* > > inclusive_rho_Central_HLT, inclusive_ptave_Central_HLT, inclusive_nevents_Central_HLT; // alpha, pt
  std::vector< std::vector< TH1F* > > inclusive_rho_HF_HLT, inclusive_ptave_HF_HLT, inclusive_nevents_HF_HLT; // alpha, pt
  std::vector< std::vector< TH1F* > > inclusive_rho_Central_L1HLT, inclusive_ptave_Central_L1HLT, inclusive_nevents_Central_L1HLT; // alpha, pt
  std::vector< std::vector< TH1F* > > inclusive_rho_HF_L1HLT, inclusive_ptave_HF_L1HLT, inclusive_nevents_HF_L1HLT; // alpha, pt

  TH1F *h_PU;
  TH1F *h_alpha_raw;
  TH1F *h_alpha_sel;
  TH1F *h_JetAvePt_SM;
  TH1F *h_Jet1Pt_SM;
  TH1F *h_Jet2Pt_SM;
  TH1F *h_Jet3Pt_SM;
  TH1F *h_JetAvePt_FE;
  TH1F *h_Jet1Pt_FE;
  TH1F *h_Jet2Pt_FE;
  TH1F *h_Jet3Pt_FE;

  TH2F *central_PS_pt, *central_PS_rho;
};

#endif

#ifdef MySelector_cxx
void MySelector::Init(TTree *tree){
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);

  TFile* currentFile = ((TChain*)fChain)->GetFile();

  fChain->SetBranchAddress("weight", &weight, &b_weight);
  fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
  fChain->SetBranchAddress("Njet", &njet, &b_njet);
  fChain->SetBranchAddress("is_JER_SM", &is_JER_SM, &b_is_JER_SM);
  fChain->SetBranchAddress("pt_ave", &pt_ave, &b_pt_ave);
  fChain->SetBranchAddress("barreljet_phi", &barreljet_phi, &b_barreljet_phi);
  fChain->SetBranchAddress("barreljet_eta", &barreljet_eta, &b_barreljet_eta);
  fChain->SetBranchAddress("barreljet_pt", &barreljet_pt, &b_barreljet_pt);
  fChain->SetBranchAddress("probejet_phi", &probejet_phi, &b_probejet_phi);
  fChain->SetBranchAddress("probejet_eta", &probejet_eta, &b_probejet_eta);
  fChain->SetBranchAddress("probejet_pt", &probejet_pt, &b_probejet_pt);
  fChain->SetBranchAddress("jet3_pt", &jet3_pt, &b_jet3_pt);
  fChain->SetBranchAddress("asymmetry", &asymmetry, &b_asymmetry);
  fChain->SetBranchAddress("alpha", &alpha_, &b_alpha);
  fChain->SetBranchAddress("rho", &rho, &b_rho);
  fChain->SetBranchAddress("prescale", &HLT, &b_prescale);
  fChain->SetBranchAddress("prescale_L1min", &L1min, &b_L1min);
  fChain->SetBranchAddress("prescale_L1max", &L1max, &b_L1max);

  fChain->SetBranchAddress("lumi_sec", &lumi_sec, &b_lumi_sec);
  fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
  fChain->SetBranchAddress("run", &run, &b_run);
}

bool MySelector::Notify(){
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  TFile *currentFile = fChain->GetCurrentFile();
  return kTRUE;
}

#endif // #ifdef MySelector_cxx
