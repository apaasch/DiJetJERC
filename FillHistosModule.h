#pragma once

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


#include <string>
#include <vector>
#include <memory>
#include <limits>


class MySelector : public TSelector {
public:

  TString outdir;
  TTree *fChain;   //!pointer to the analyzed TTree or TChain

  // Declaration of inputs

  std::vector<std::string> input_names = {"rho", "pt_ave", "asymmetry", "alpha", "barreljet_phi", "barreljet_eta", "barreljet_pt", "probejet_phi", "probejet_eta", "probejet_pt", "jet3_pt", "jet3_eta", "jet3_phi" };
  std::vector<std::string> input_names_int = {"Njet", "is_JER_SM", "hf_trigger", "bl_trigger"};

  std::map<std::string, float> inputs;
  std::map<std::string, int> inputs_int;
  std::map<std::string, TBranch*> branches;

  for (const std::string& name : input_names) {
    inputs[name]    = std::numeric_limits<float>::infinity();
    branches[name]  = NULL;
    fChain->SetBranchAddress(name, &(inputs[name]), &(branches[name]));
  }

  for (const std::string& name : input_names_int) {
    inputs_int[name]    = std::numeric_limits<float>::infinity();
    branches[name]  = NULL;
    fChain->SetBranchAddress(name, &(inputs_int[name]), &(branches[name]));
  }

  int TotalEvents;

  int EtaBins_SM, EtaBins_SM_control, EtaBins_FE_reference, EtaBins_FE_control, EtaBins_FE;
  int etaShift_SM, etaShift_SM_control, etaShift_FE_reference, etaShift_FE_control, etaShift_FE;
  int PtBins_Central, PtBins_HF;
  int AlphaBins;


  // std::vector< std::vector< std::vector< double > > > nevents_central,nevents_HF;
  // std::vector< double > nevents_central,nevents_HF;
  std::vector< std::vector< std::vector< double > > > nevents_central,nevents_HF;

  std::map<std::string, std::vector< std::vector< std::vector< TH1F* > > > > asyHistos;

  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_SM, 						asymmetries_pt_SM,						asymmetries_rho_SM,						asymmetries_pt3_SM;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_SM_control, 		asymmetries_pt_SM_control,		asymmetries_rho_SM_control,		asymmetries_pt3_SM_control;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_FE_reference, 	asymmetries_pt_FE_reference,	asymmetries_rho_FE_reference,	asymmetries_pt3_FE_reference;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_FE_control, 		asymmetries_pt_FE_control,		asymmetries_rho_FE_control,		asymmetries_pt3_FE_control;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_FE, 						asymmetries_pt_FE,						asymmetries_rho_FE,						asymmetries_pt3_FE;

  std::vector< std::vector< TH1F* > > alpha_spectrum_SM, alpha_spectrum_SM_control, alpha_spectrum_FE_reference, alpha_spectrum_FE_control, alpha_spectrum_FE;

  std::vector<TH1F*> histograms;

  TH1F *h_PU;
  TH1F *h_alpha_raw;
  TH1F *h_alpha_select;

  TH1F *h_JetAvePt_SM;
  TH1F *h_Jet1Pt_SM;
  TH1F *h_Jet2Pt_SM;
  TH1F *h_Jet3Pt_SM;

  TH1F *h_JetAvePt_FE;
  TH1F *h_Jet1Pt_FE;
  TH1F *h_Jet2Pt_FE;
  TH1F *h_Jet3Pt_FE;



  MySelector(TString name, TTree * /*tree*/ =0) : fChain(0), outdir(name) { }
  virtual ~MySelector() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) {
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
  fChain->SetBranchAddress("jet3_eta", &jet3_eta, &b_jet3_eta);
  fChain->SetBranchAddress("jet3_phi", &jet3_phi, &b_jet3_phi);
  fChain->SetBranchAddress("asymmetry", &asymmetry, &b_asymmetry);
  fChain->SetBranchAddress("alpha", &alpha_, &b_alpha);
  fChain->SetBranchAddress("hf_trigger", &pass_trigger_hf, &b_pass_trigger_hf);
  fChain->SetBranchAddress("bl_trigger", &pass_trigger_bl, &b_pass_trigger_bl);
  fChain->SetBranchAddress("rho", &rho, &b_rho);
}

Bool_t MySelector::Notify(){
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  TFile *currentFile = fChain->GetCurrentFile();
  return kTRUE;
}

#endif // #ifdef MySelector_cxx
