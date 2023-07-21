#define MySelector_cxx
// The class definition in MySelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("MySelector.C")
// Root > T->Process("MySelector.C","some options")
// Root > T->Process("MySelector.C+")
//

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include "MySelector.h"
#include "constants.h"

bool JetInRange(double jet_eta, double min, double max) {
  return (TMath::Abs(jet_eta) > min && TMath::Abs(jet_eta) < max);
}

bool JetInEtaBin(double jet_eta, std::vector<double> bins, int bin) {
  return JetInRange(jet_eta, bins[bin], bins[bin+1]);
}


#define FILL_HISTOS(region,method)                                                                          \
asymmetries_##region.at(r).at(k).at(m)->Fill( asy , weight);                                                \
asymmetries_pt_##region.at(r).at(k).at(m)->Fill( pt_ave, weight);                                           \
asymmetries_rho_##region.at(r).at(k).at(m)->Fill( rho, weight);                                             \
asymmetries_pt3_##region.at(r).at(k).at(m)->Fill( jet3_pt, weight);                                         \
pt_spectrum_##region.at(r).at(k).at(m)->Fill( pt_ave, weight);                                              \
if ( m == AlphaBins-1 ) {                                                                                   \
  h_JetAvePt_##method->Fill( pt_ave, weight);                                                               \
  h_Jet1Pt_##method->Fill( barreljet_pt, weight);                                                           \
  h_Jet2Pt_##method->Fill( probejet_pt, weight);                                                            \
  h_Jet3Pt_##method->Fill( jet3_pt, weight);                                                                \
}                                                                                                           \


#define SELECT_ETA_ALPHA_BIN(region,method,cond1,cond2)     \
if (cond1 || cond2) {                                       \
  h_alpha_sel->Fill(alpha, 1);                              \
  for ( int m = 0 ; m < AlphaBins ; m++ ) {                 \
    if ( alpha < Alpha_bins[m] ) {                          \
      if (dofill) nevents_central[k][r+shift][m] +=1;       \
      else nevents_HF[k][r+shift][m] +=1;                   \
      double asy = asymmetry;                               \
      FILL_HISTOS(region,method)                            \
      if ( excl_bin ) break;                                \
    }                                                       \
  }                                                         \
  alpha_spectrum_##region.at(r).at(k)->Fill(alpha, weight); \
}                                                           \

#define WRITE_HISTOS(region)                                    \
for( int m = 0; m < EtaBins_##region; m++ ) {                   \
  f->cd();                                                      \
  asymmetries_##region.at(m).at(p).at(r)->Write();              \
  asymmetries_pt_##region.at(m).at(p).at(r)->Write();           \
  pt_spectrum_##region.at(m).at(p).at(r)->Write();              \
  f1->cd();                                                     \
  asymmetries_rho_##region.at(m).at(p).at(r)->Write();          \
  asymmetries_pt3_##region.at(m).at(p).at(r)->Write();          \
  f_alpha->cd();                                                \
  alpha_spectrum_##region.at(m).at(p)->Write();                 \
}                                                               \

void MakeHistograms(std::vector< std::vector< std::vector< TH1F* > > > &asymmetries, std::vector< std::vector< std::vector< TH1F* > > > &asymmetries_pt, std::vector< std::vector< std::vector< TH1F* > > > &asymmetries_rho, std::vector< std::vector< std::vector< TH1F* > > > &asymmetries_pt3, std::vector< std::vector< std::vector< TH1F* > > > &pt_spectrum, std::vector< std::vector< TH1F* > > &alpha_spectrum, TString text, TString extraText, int etaBins, int ptBins, int AlphaBins, int etaShift, int ptShift, int alphaShift) {
  for( int m = etaShift; m < etaBins+etaShift; m++ ) {
    std::vector< std::vector< TH1F* > > temp2, temp2pt, temp2rho, temp2pt3, temp2pts;
    std::vector< TH1F* > alpha_temp2;


    for( int p = 0; p < ptBins; p++ ) {

      // =======================================================================
      // === Study alpha at low pt
      // int abins = AlphaBins;
      // if(p<3) abins = 13;
      // else if(p>5 && p<9) abins = 13;
      // else if(p>13 && p<16) abins = 13;
      // else if(p>21) abins = 13;

      std::vector< TH1F* > temp1, temp1pt, temp1rho, temp1pt3, temp1pts;
      TString name_alpha = "alpha"; name_alpha += extraText; name_alpha += "_eta"; name_alpha += m+1; name_alpha += "_pt"; name_alpha += p+1;
      TH1F *h1_alpha = new TH1F( name_alpha, name_alpha, 80, 0., 0.8); h1_alpha->SetXTitle("Alpha"); h1_alpha->SetYTitle("a.u.");  h1_alpha->Sumw2();  alpha_temp2.push_back(h1_alpha);
      for( int r = 0; r < AlphaBins; r++ ) {
      // for( int r = 0; r < abins; r++ ) {
        TString name     = text;        name     += extraText; name     += "_eta"; name     += m+1; name     += "_pt"; name     += p+1; name     += "_alpha"; name     += r+1;
        TString name_pt  = text+"pt";   name_pt  += extraText; name_pt  += "_eta"; name_pt  += m+1; name_pt  += "_pt"; name_pt  += p+1; name_pt  += "_alpha"; name_pt  += r+1;
        TString name_rho = text+"rho";  name_rho += extraText; name_rho += "_eta"; name_rho += m+1; name_rho += "_pt"; name_rho += p+1; name_rho += "_alpha"; name_rho += r+1;
        TString name_pt3 = text+"pt3";  name_pt3 += extraText; name_pt3 += "_eta"; name_pt3 += m+1; name_pt3 += "_pt"; name_pt3 += p+1; name_pt3 += "_alpha"; name_pt3 += r+1;
        TString name_pts = text+"ptspectrum";   name_pts += extraText; name_pts += "_eta"; name_pts += m+1; name_pts += "_pt"; name_pts += p+1; name_pts += "_alpha"; name_pts += r+1;
        TH1F *h1 = new TH1F(name,     name,     160,-0.8, 0.8);   h1->SetXTitle("Asymmetry"); h1->SetYTitle("a.u.");  h1->Sumw2();  temp1.push_back(h1);
        TH1F *h2 = new TH1F(name_pt,  name_pt,  100,  0,  7000);  h2->SetXTitle("Pt[GeV]");   h2->SetYTitle("a.u.");  h2->Sumw2();  temp1pt.push_back(h2);
        TH1F *h3 = new TH1F(name_rho, name_rho, 100,  0,  100);   h3->SetXTitle("rho");       h3->SetYTitle("a.u.");  h3->Sumw2();  temp1rho.push_back(h3);
        TH1F *h4 = new TH1F(name_pt3, name_pt3, 7000, 0,  7000);  h4->SetXTitle("Pt[GeV]");   h4->SetYTitle("a.u.");  h4->Sumw2();  temp1pt3.push_back(h4);
        TH1F *h5 = new TH1F(name_pts, name_pts, 7000, 0,  7000);  h5->SetXTitle("Pt[GeV]");   h5->SetYTitle("a.u.");  h5->Sumw2();  temp1pts.push_back(h5);
      }
      temp2.push_back(temp1); temp2pt.push_back(temp1pt); temp2rho.push_back(temp1rho);  temp2pt3.push_back(temp1pt3);  temp2pts.push_back(temp1pts);
    }
    asymmetries.push_back(temp2); asymmetries_pt.push_back(temp2pt); asymmetries_rho.push_back(temp2rho); asymmetries_pt3.push_back(temp2pt3); pt_spectrum.push_back(temp2pts);
    alpha_spectrum.push_back(alpha_temp2);
  }
}

void MakeEtaInclusiveHistograms(
  std::vector< std::vector< TH1F* > > &rho,
  std::vector< std::vector< TH1F* > > &pt,
  std::vector< std::vector< TH1F* > > &nevents,
  TString text, int ptBins, int AlphaBins
) {
  for( int p = 0; p < ptBins; p++ ) {
    std::vector< TH1F* > temp_rho, temp_pt, temp_ne;
    for( int r = 0; r < AlphaBins; r++ ) {
      TString name_pt  = "ptave";   name_pt  += text; name_pt  += "_pt"; name_pt  += p+1; name_pt  += "_alpha"; name_pt  += r+1;
      TString name_rho = "rho";     name_rho += text; name_rho += "_pt"; name_rho += p+1; name_rho += "_alpha"; name_rho += r+1;
      TString name_ne  = "nevents"; name_ne  += text; name_ne  += "_pt"; name_ne  += p+1; name_ne  += "_alpha"; name_ne  += r+1;

      TH1F *h1 = new TH1F(name_pt,  name_pt,  7000, 0, 7000);  h1->SetXTitle("Pt[GeV]");   h1->SetYTitle("Events");  h1->Sumw2();  temp_pt.push_back(h1);
      TH1F *h2 = new TH1F(name_rho, name_rho,   80, 0,   80);  h2->SetXTitle("rho");       h2->SetYTitle("Events");  h2->Sumw2();  temp_rho.push_back(h2);
      TH1F *h3 = new TH1F(name_ne,  name_ne,     1, 0,    1);  h3->SetXTitle("");          h3->SetYTitle("Events");  h3->Sumw2();  temp_ne.push_back(h3);
    }
    rho.push_back(temp_rho); pt.push_back(temp_pt); nevents.push_back(temp_ne);
  }
}

void MySelector::BuildEvent() {
}

void MySelector::Begin(TTree * /*tree*/) {
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  TotalEvents = 0;
}

void MySelector::SlaveBegin(TTree * /*tree*/) {

  TString option = GetOption();

  isQuick = (outdir.Contains("quick"))?true:false;
  std::cout << outdir << std::endl; 
  if(isQuick) std::cout << "Only run a few events " << std::endl; 

  std::vector<double> eta_bins;

  if (study.find("eta_narrow") != std::string::npos)      eta_bins = std::vector<double>(eta_bins_narrow, eta_bins_narrow + n_eta_bins_narrow);
  else if (study.find("eta_simple") != std::string::npos) eta_bins = std::vector<double>(eta_bins_simple, eta_bins_simple + n_eta_bins_simple);
  else if (study.find("eta_common") != std::string::npos) eta_bins = std::vector<double>(eta_bins_common, eta_bins_common + n_eta_bins_common);
  else if (study.find("eta_calo") != std::string::npos)   eta_bins = std::vector<double>(eta_bins_calo, eta_bins_calo + n_eta_bins_calo);
  else if (study.find("eta_L2R") != std::string::npos)    eta_bins = std::vector<double>(eta_bins_L2R, eta_bins_L2R + n_eta_bins_L2R);
  else                          eta_bins = std::vector<double>(eta_bins_JER, eta_bins_JER + n_eta_bins_JER);

  int n_eta_bins = eta_bins.size();
  EtaBins = n_eta_bins; // For SlaveTerminate

  EtaBins_SM            = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return i<eta_cut; });
  EtaBins_SM_control    = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return i>eta_cut; });
  EtaBins_FE_reference  = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return i<s_eta_barr;});
  EtaBins_FE_control    = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return (i>=s_eta_barr)&&(i<eta_cut);});
  EtaBins_FE            = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return i>eta_cut; });

  etaShift_SM           = 0;
  etaShift_SM_control   = EtaBins_SM;
  etaShift_FE_reference = 0;
  etaShift_FE_control   = EtaBins_FE_reference;
  etaShift_FE           = EtaBins_FE_reference + EtaBins_FE_control;

  for (size_t i = etaShift_SM;            i < etaShift_SM           + EtaBins_SM            + 1; i++)  Eta_bins_SM.push_back(eta_bins[i]);
  for (size_t i = etaShift_SM_control;    i < etaShift_SM_control   + EtaBins_SM_control    + 1; i++)  Eta_bins_SM_control.push_back(eta_bins[i]);
  for (size_t i = etaShift_FE_reference;  i < etaShift_FE_reference + EtaBins_FE_reference  + 1; i++)  Eta_bins_FE_reference.push_back(eta_bins[i]);
  for (size_t i = etaShift_FE_control;    i < etaShift_FE_control   + EtaBins_FE_control    + 1; i++)  Eta_bins_FE_control.push_back(eta_bins[i]);
  for (size_t i = etaShift_FE;            i < etaShift_FE           + EtaBins_FE            + 1; i++)  Eta_bins_FE.push_back(eta_bins[i]);

  std::string triggerName = isAK8? "SingleJet" : "DiJet";
  std::string name_pt_bin = triggerName+"_central_";
  if (isAK8) name_pt_bin += "AK8_";
  name_pt_bin += year+"_ptbins";
  if(binning != "") name_pt_bin += "_"+binning;
  PtBins_Central = pt_trigger_thr.at(name_pt_bin).size();
  for (auto &pt: pt_trigger_thr.at(name_pt_bin)) Pt_bins_Central.push_back(pt);
  name_pt_bin = triggerName+"_forward_";
  if (isAK8) name_pt_bin += "AK8_";
  name_pt_bin += year+"_ptbins";
  if(binning != "") name_pt_bin += "_"+binning;
  PtBins_HF = pt_trigger_thr.at(name_pt_bin).size();
  for (auto &pt: pt_trigger_thr.at(name_pt_bin)) Pt_bins_HF.push_back(pt);
  Pt_bins_Central.push_back(7000);
  Pt_bins_HF.push_back(7000);

  PtBins=(PtBins_Central>PtBins_HF)?PtBins_Central:PtBins_HF;

  isPrescale = (outdir.Contains("prescale"))?true:false;
  if(isPrescale) std::cout << "Apply prescale to data; set weight to combined HLT & L1 seed prescale weight." << std::endl;

  Alpha_bins = {0.05,0.10,0.15,0.20,0.25,0.30}; // default
  if(abins.find("finealpha") != std::string::npos) Alpha_bins = {0.05,0.075,0.10,0.125,0.15,0.175,0.20,0.225,0.25,0.275,0.30};
  if(abins.find("highalpha") != std::string::npos) Alpha_bins = {0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.};
  AlphaBins = Alpha_bins.size();

  Alpha_bins_Inc = {0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.};
  AlphaBinsInc = Alpha_bins_Inc.size();

  h_JetAvePt_SM = new TH1F("JetAvePt" ,   "Inclusive Jet Ave Pt",   50, 0., 2000);  h_JetAvePt_SM->SetXTitle("Pt_{ave}[GeV]");  h_JetAvePt_SM->SetYTitle("a.u."); h_JetAvePt_SM->Sumw2();
  h_Jet1Pt_SM   = new TH1F("Jet1Pt",      "Inclusive Jet 1 Pt",     50, 0., 2000);  h_Jet1Pt_SM->SetXTitle("Pt_1[GeV]");        h_Jet1Pt_SM->SetYTitle("a.u.");   h_Jet1Pt_SM->Sumw2();
  h_Jet2Pt_SM   = new TH1F("Jet2Pt",      "Inclusive Jet 2 Pt",     50, 0., 2000);  h_Jet2Pt_SM->SetXTitle("Pt_2[GeV]");        h_Jet2Pt_SM->SetYTitle("a.u.");   h_Jet2Pt_SM->Sumw2();
  h_Jet3Pt_SM   = new TH1F("Jet3Pt",      "Inclusive Jet 3 Pt",     50, 0., 2000);  h_Jet3Pt_SM->SetXTitle("Pt_3[GeV]");        h_Jet3Pt_SM->SetYTitle("a.u.");   h_Jet3Pt_SM->Sumw2();
  h_JetAvePt_FE = new TH1F("FEJetAvePt",  "Inclusive FEJet Ave Pt", 50, 0., 2000);  h_JetAvePt_FE->SetXTitle("Pt_{ave}[GeV]");  h_JetAvePt_FE->SetYTitle("a.u."); h_JetAvePt_FE->Sumw2();
  h_Jet1Pt_FE   = new TH1F("FEJet1Pt",    "Inclusive FEJet 1 Pt",   50, 0., 2000);  h_Jet1Pt_FE->SetXTitle("Pt_1[GeV]");        h_Jet1Pt_FE->SetYTitle("a.u.");   h_Jet1Pt_FE->Sumw2();
  h_Jet2Pt_FE   = new TH1F("FEJet2Pt",    "Inclusive FEJet 2 Pt",   50, 0., 2000);  h_Jet2Pt_FE->SetXTitle("Pt_2[GeV]");        h_Jet2Pt_FE->SetYTitle("a.u.");   h_Jet2Pt_FE->Sumw2();
  h_Jet3Pt_FE   = new TH1F("FEJet3Pt",    "Inclusive FEJet 3 Pt",   50, 0., 2000);  h_Jet3Pt_FE->SetXTitle("Pt_3[GeV]");        h_Jet3Pt_FE->SetYTitle("a.u.");   h_Jet3Pt_FE->Sumw2();
  h_PU          = new TH1F("PileUp",      "PU distribution",        60, 0., 60);    h_PU->SetXTitle("PU");                      h_PU->SetYTitle("a.u.");          h_PU->Sumw2();
  h_alpha_raw   = new TH1F("Alpha_raw",   "#alpha before selection",80, 0., 0.8);  h_alpha_raw->SetXTitle("#alpha_raw");       h_alpha_raw->SetYTitle("a.u.");   h_alpha_raw->Sumw2();
  h_alpha_sel   = new TH1F("Alpha",       "#alpha after selection", 80, 0., 0.8);   h_alpha_sel->SetXTitle("#alpha");           h_alpha_sel->SetYTitle("a.u.");   h_alpha_sel->Sumw2();

  if (true) {
    std::cout << "\nConstructor: " << PtBins_Central << " " << PtBins_HF;
    std::cout << "\nPt_bins_Central: ";         for (size_t i = 0; i < Pt_bins_Central.size(); i++) std::cout << " " << Pt_bins_Central[i];
    std::cout << "\nPt_bins_HF: ";              for (size_t i = 0; i < Pt_bins_HF.size(); i++) std::cout << " " << Pt_bins_HF[i];
    std::cout << "\nAlpha_bins: ";              for (size_t i = 0; i < Alpha_bins.size(); i++) std::cout << " " << Alpha_bins[i];
    std::cout << "\nEta_bins_SM : ";            for (size_t i = 0; i < Eta_bins_SM.size(); i++) std::cout << " " << Eta_bins_SM[i];
    std::cout << "\nEta_bins_SM_control : ";    for (size_t i = 0; i < Eta_bins_SM_control.size(); i++) std::cout << " " << Eta_bins_SM_control[i];
    std::cout << "\nEta_bins_FE_reference : ";  for (size_t i = 0; i < Eta_bins_FE_reference.size(); i++) std::cout << " " << Eta_bins_FE_reference[i];
    std::cout << "\nEta_bins_FE_control : ";    for (size_t i = 0; i < Eta_bins_FE_control.size(); i++) std::cout << " " << Eta_bins_FE_control[i];
    std::cout << "\nEta_bins_FE : ";            for (size_t i = 0; i < Eta_bins_FE.size(); i++) std::cout << " " << Eta_bins_FE[i];
    std::cout << "\n" << std::endl;
  }

  int neta = (int) eta_bins.size();
  for ( int k = 0 ; k < PtBins_Central ; k++ ) {
    std::vector< std::vector< double > >  temp;
    for (int r = 0; r < neta; r++) {
      std::vector< double >  temp2;
      for (int m = 0; m < AlphaBins; m++) temp2.push_back(0);
      temp.push_back(temp2);
    }
    nevents_central.push_back(temp);
  }

  for ( int k = 0 ; k < PtBins_HF ; k++ ) {
    std::vector< std::vector< double > >  temp;
    for (int r = 0; r < neta; r++) {
      std::vector< double >  temp2;
      for (int m = 0; m < AlphaBins; m++) temp2.push_back(0);
      temp.push_back(temp2);
    }
    nevents_HF.push_back(temp);
  }

  MakeHistograms(asymmetries_SM,            asymmetries_pt_SM,            asymmetries_rho_SM,            asymmetries_pt3_SM,            pt_spectrum_SM,           alpha_spectrum_SM,           "asymm", "_SM",            EtaBins_SM,            PtBins_Central, AlphaBins, etaShift_SM,            0, 0);
  MakeHistograms(asymmetries_SM_control,    asymmetries_pt_SM_control,    asymmetries_rho_SM_control,    asymmetries_pt3_SM_control,    pt_spectrum_SM_control,   alpha_spectrum_SM_control,   "asymm", "_SM_control",    EtaBins_SM_control,    PtBins_HF,      AlphaBins, etaShift_SM_control,    0, 0);
  MakeHistograms(asymmetries_FE_reference,  asymmetries_pt_FE_reference,  asymmetries_rho_FE_reference,  asymmetries_pt3_FE_reference,  pt_spectrum_FE_reference, alpha_spectrum_FE_reference, "asymm", "_FE_reference",  EtaBins_FE_reference,  PtBins_Central, AlphaBins, etaShift_FE_reference,  0, 0);
  MakeHistograms(asymmetries_FE_control,    asymmetries_pt_FE_control,    asymmetries_rho_FE_control,    asymmetries_pt3_FE_control,    pt_spectrum_FE_control,   alpha_spectrum_FE_control,   "asymm", "_FE_control",    EtaBins_FE_control,    PtBins_Central, AlphaBins, etaShift_FE_control,    0, 0);
  MakeHistograms(asymmetries_FE,            asymmetries_pt_FE,            asymmetries_rho_FE,            asymmetries_pt3_FE,            pt_spectrum_FE,           alpha_spectrum_FE,           "asymm", "_FE",            EtaBins_FE,            PtBins_HF,      AlphaBins, etaShift_FE,            0, 0);

  MakeEtaInclusiveHistograms(inclusive_rho_Central, inclusive_ptave_Central, inclusive_nevents_Central, "_central", PtBins, AlphaBinsInc);
  MakeEtaInclusiveHistograms(inclusive_rho_HF, inclusive_ptave_HF, inclusive_nevents_HF, "_HF", PtBins, AlphaBinsInc);

  // MakeEtaInclusiveHistograms(inclusive_rho_Central_HLT, inclusive_ptave_Central_HLT, inclusive_nevents_Central_HLT, "_central_HLT", PtBins, AlphaBinsInc);
  // MakeEtaInclusiveHistograms(inclusive_rho_HF_HLT, inclusive_ptave_HF_HLT, inclusive_nevents_HF_HLT, "_HF_HLT", PtBins, AlphaBinsInc);

  MakeEtaInclusiveHistograms(inclusive_rho_Central_L1HLT, inclusive_ptave_Central_L1HLT, inclusive_nevents_Central_L1HLT, "_central_L1HLT", PtBins, AlphaBinsInc);
  MakeEtaInclusiveHistograms(inclusive_rho_HF_L1HLT, inclusive_ptave_HF_L1HLT, inclusive_nevents_HF_L1HLT, "_HF_L1HLT", PtBins, AlphaBinsInc);

  central_PS_pt = new TH2F("central_PS_pt", "central_PS_pt", 3000, 0, 3000, 160, 0, 160); central_PS_pt->SetXTitle("Pt [GeV]"); central_PS_pt->SetYTitle("HLT Prescale");  central_PS_pt->Sumw2();
  central_PS_rho = new TH2F("central_PS_rho", "central_PS_rho", 80, 0, 80, 160, 0, 160); central_PS_rho->SetXTitle("rho [GeV]"); central_PS_rho->SetYTitle("HLT Prescale");  central_PS_rho->Sumw2();

}

bool MySelector::Process(Long64_t entry) {

  ++TotalEvents;
  if ( TotalEvents%1000000 == 0 ) {  std::cout << "\t\tAnalyzing event #" << TotalEvents << std::endl; }
  if(isQuick){ // skip many events
    if(probejet_pt>3000) return kTRUE;
  }

  GetEntry(entry);
  BuildEvent();

  bool cond1, cond2;

  // double jet_thr = jet_threshold_min;
  double jet_thr=15;
  double alpha_raw = alpha_;
  double alpha;

  // Below I choose what kind of asymmetries I want to study! excl_bin = true for exclusive bins
  bool excl_bin = false; // inclusive
  if (njet<2) return kTRUE;
  if (barreljet_pt < jet_thr && probejet_pt < jet_thr) return kTRUE;
  if ( TMath::Abs(TVector2::Phi_mpi_pi((probejet_phi - barreljet_phi))) < s_delta_phi ) { std::cout << "Jets are not back to back" << std::endl; return kTRUE;}

  if ( jet3_pt > jet_thr ) alpha = TMath::Abs(alpha_raw);
  else alpha = (njet <3) ? 0. : 1. ;

  h_PU->Fill(nPU, 1);
  h_alpha_raw->Fill(alpha_raw, 1);

  bool dofill; int shift;
  bool isHF = TMath::Abs(probejet_eta)>eta_cut? true : false;

  if(L1min!=L1max) std::cout << "Warning - More than one L1T Seed - L1min=" << std::setw(7) << L1min <<  " and L1max=" << std::setw(7) << L1max << std::setw(20) << probejet_pt << std::setw(20) << barreljet_pt << std::endl;
  double w_HLT = weight*HLT;
  double w_prescale = weight*HLT*L1max;
  if(isPrescale) weight = w_prescale;

  // DELETE LATER - Calculate is_JER_SM for eta_calo since not PreSel was run here
  if (study=="eta_calo"){
    is_JER_SM = false;
    for (unsigned int i = 0; i < EtaBins-1; i++) {
      if(fabs(barreljet_eta)>=eta_bins[i] && fabs(barreljet_eta)<eta_bins[i+1]) {
        if (fabs(probejet_eta)>=eta_bins[i] && fabs(probejet_eta)<eta_bins[i+1]) {
          is_JER_SM = true;
        }
      }
    }
  }
  // DELETE LATER - Calculate is_JER_SM for eta_calo since not PreSel was run here

  if (!isHF) {
    dofill=true;
    for ( int k = 0 ; k < PtBins_Central ; k++ ) {
      if ((pt_ave > Pt_bins_Central[k]) && (pt_ave < Pt_bins_Central[k+1]) ) {

        // if (!is_JER_SM){ // First only FE
        //   if( JetInRange(barreljet_eta, 0, s_eta_barr) || JetInRange(probejet_eta,  0, s_eta_barr) ) {
        //     for(int a = 0; a<AlphaBinsInc; a++){
        //       if ( alpha < Alpha_bins_Inc[a] ){
        //         inclusive_rho_Central[k][a]->Fill(rho, weight);
        //         inclusive_ptave_Central[k][a]->Fill(pt_ave, weight);
        //         inclusive_nevents_Central[k][a]->Fill(0.5, weight);
        //       }
        //     }
        //   }
        // }

        for(int a = 0; a<AlphaBinsInc; a++){
          if ( alpha < Alpha_bins_Inc[a] ){
            inclusive_rho_Central[k][a]->Fill(rho, weight);
            inclusive_ptave_Central[k][a]->Fill(pt_ave, weight);
            inclusive_nevents_Central[k][a]->Fill(0.5, weight);

            // inclusive_rho_Central_HLT[k][a]->Fill(rho, w_HLT);
            // inclusive_ptave_Central_HLT[k][a]->Fill(pt_ave, w_HLT);
            // inclusive_nevents_Central_HLT[k][a]->Fill(0.5, w_HLT);

            inclusive_rho_Central_L1HLT[k][a]->Fill(rho, w_prescale);
            inclusive_ptave_Central_L1HLT[k][a]->Fill(pt_ave, w_prescale);
            inclusive_nevents_Central_L1HLT[k][a]->Fill(0.5, w_prescale);

            central_PS_pt->Fill(pt_ave, HLT, weight);
            central_PS_rho->Fill(rho, HLT, weight);

          }
        }

        for ( int r = 0 ; r < EtaBins_SM ; r++ ) {
          if (!is_JER_SM) continue;
          cond1 = (JetInEtaBin(barreljet_eta, Eta_bins_SM, r) && JetInEtaBin(probejet_eta, Eta_bins_SM, r));
          cond2 = false;
          // std::cout << "Eta_bins_SM " << Eta_bins_SM[r] << "-" << Eta_bins_SM[r+1] << " " << cond1 << " " << cond2 << std::endl;
          shift = 0;
          SELECT_ETA_ALPHA_BIN(SM,SM,cond1,cond2)
        }
        for ( int r = 0 ; r < EtaBins_FE_reference ; r++ ) {
          if (is_JER_SM) continue;
          cond1 = (JetInRange(barreljet_eta, 0, s_eta_barr) && JetInEtaBin(probejet_eta,  Eta_bins_FE_reference, r));
          cond2 = (JetInRange(probejet_eta,  0, s_eta_barr) && JetInEtaBin(barreljet_eta, Eta_bins_FE_reference, r));
          shift = 0;
          // std::cout << "EtaBins_FE_reference " << Eta_bins_FE_reference[r] << "-" << Eta_bins_FE_reference[r+1] << " " << cond1 << " " << cond2 << std::endl;
          SELECT_ETA_ALPHA_BIN(FE_reference,FE,cond1,cond2)
        }
        for ( int r = 0 ; r < EtaBins_FE_control ; r++ ) {
          if (is_JER_SM) continue;
          cond1 = (JetInRange(barreljet_eta, 0, s_eta_barr) && JetInEtaBin(probejet_eta,  Eta_bins_FE_control, r));
          cond2 = (JetInRange(probejet_eta,  0, s_eta_barr) && JetInEtaBin(barreljet_eta, Eta_bins_FE_control, r));
          shift = etaShift_FE_control;
          // std::cout << "EtaBins_FE_control " << Eta_bins_FE_control[r] << "-" << Eta_bins_FE_control[r+1] << " " << cond1 << " " << cond2 << std::endl;
          SELECT_ETA_ALPHA_BIN(FE_control,FE,cond1,cond2)
        }
        break;
      }
    }
  } else {
    dofill=false;
    for ( int k = 0 ; k < PtBins_HF ; k++ ) {
      if ((pt_ave > Pt_bins_HF[k]) && (pt_ave < Pt_bins_HF[k+1]) ) {

        if (!is_JER_SM){ // First only FE
          if(JetInRange(barreljet_eta, 0, s_eta_barr) || JetInRange(probejet_eta,  0, s_eta_barr)){
            for(int a = 0; a<AlphaBinsInc; a++){
              if ( alpha < Alpha_bins_Inc[a] ){
                inclusive_rho_HF[k][a]->Fill(rho, weight);
                inclusive_ptave_HF[k][a]->Fill(pt_ave, weight);
                inclusive_nevents_HF[k][a]->Fill(0.5, weight);

                // inclusive_rho_HF_HLT[k][a]->Fill(rho, w_HLT);
                // inclusive_ptave_HF_HLT[k][a]->Fill(pt_ave, w_HLT);
                // inclusive_nevents_HF_HLT[k][a]->Fill(0.5, w_HLT);

                inclusive_rho_HF_L1HLT[k][a]->Fill(rho, w_prescale);
                inclusive_ptave_HF_L1HLT[k][a]->Fill(pt_ave, w_prescale);
                inclusive_nevents_HF_L1HLT[k][a]->Fill(0.5, w_prescale);
              }
            }
          }
        }

        for ( int r = 0 ; r < EtaBins_SM_control; r++ ) {
          if (!is_JER_SM) continue;
          cond1 = (JetInEtaBin(barreljet_eta, Eta_bins_SM_control, r) && JetInEtaBin(probejet_eta, Eta_bins_SM_control, r));
          cond2 = false;
          shift = etaShift_SM_control;
          // std::cout << "Eta_bins_SM_control " << Eta_bins_SM_control[r] << "-" << Eta_bins_SM_control[r+1] << " " << cond1 << " " << cond2 << std::endl;
          SELECT_ETA_ALPHA_BIN(SM_control,SM,cond1,cond2)
        }
        for ( int r = 0 ; r < EtaBins_FE ; r++ ) {
          if (is_JER_SM) continue;
          cond1 = (JetInRange(barreljet_eta, 0, s_eta_barr) && JetInEtaBin(probejet_eta, Eta_bins_FE, r));
          cond2 = (JetInRange(probejet_eta,  0, s_eta_barr) && JetInEtaBin(barreljet_eta, Eta_bins_FE, r));
          // std::cout << "Eta_bins_FE " << Eta_bins_FE[r] << "-" << Eta_bins_FE[r+1] << " " << cond1 << " " << cond2 << std::endl;
          shift = etaShift_FE;
          SELECT_ETA_ALPHA_BIN(FE,FE,cond1,cond2)
        }
        break;
      }
    }
  }
  return kTRUE;
}

void MySelector::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  std::cout <<"\t\tAnalyzed events #" <<  TotalEvents << std::endl;

  std::ofstream mytxtfile;
  mytxtfile.open (outdir+"counts.txt");
  mytxtfile << "Analyzed events #" <<  TotalEvents << "\n";
  mytxtfile.close();

  TFile *fpt = new TFile(outdir+"pt_data_incl_full.root","RECREATE"); ;
  fpt->cd();
  h_PU->Write();
  h_alpha_raw->Write();
  h_alpha_sel->Write();
  h_JetAvePt_FE->Write();
  h_Jet1Pt_FE->Write();
  h_Jet2Pt_FE->Write();
  h_Jet3Pt_FE->Write();
  h_JetAvePt_SM->Write();
  h_Jet1Pt_SM->Write();
  h_Jet2Pt_SM->Write();
  h_Jet3Pt_SM->Write();
  fpt->Close();

  TFile *f  = new TFile(outdir+"histograms_data_incl_full.root","RECREATE");
  TFile *f1 = new TFile(outdir+"histograms_data_incl_full_control.root","RECREATE");
  TFile *f2 = new TFile(outdir+"histograms_data_incl_full_AlphaCuts.root","RECREATE");
  // TFile *f3 = new TFile(outdir+"histograms_data_incl_full_AlphaCuts_HLT.root","RECREATE");
  TFile *f4 = new TFile(outdir+"histograms_data_incl_full_HLT_2D.root","RECREATE");
  TFile *f5 = new TFile(outdir+"histograms_data_incl_full_AlphaCuts_L1HLT.root","RECREATE");
  TFile *f_alpha = new TFile(outdir+"alpha_spectrum.root","RECREATE");

  f4->cd();
  central_PS_pt->Write();
  central_PS_rho->Write();

  for( int r = 0; r < AlphaBinsInc; r++ ) {
    for( int p = 0; p < PtBins_Central; p++ ) {
      f2->cd();
      inclusive_rho_Central.at(p).at(r)->Write();
      inclusive_ptave_Central.at(p).at(r)->Write();
      inclusive_nevents_Central.at(p).at(r)->Write();

      // f3->cd();
      // inclusive_rho_Central_HLT.at(p).at(r)->Write();
      // inclusive_ptave_Central_HLT.at(p).at(r)->Write();
      // inclusive_nevents_Central_HLT.at(p).at(r)->Write();

      f5->cd();
      inclusive_rho_Central_L1HLT.at(p).at(r)->Write();
      inclusive_ptave_Central_L1HLT.at(p).at(r)->Write();
      inclusive_nevents_Central_L1HLT.at(p).at(r)->Write();
    }
    for( int p = 0; p < PtBins_HF; p++ ) {
      f2->cd();
      inclusive_rho_HF.at(p).at(r)->Write();
      inclusive_ptave_HF.at(p).at(r)->Write();
      inclusive_nevents_HF.at(p).at(r)->Write();

      // f3->cd();
      // inclusive_rho_HF_HLT.at(p).at(r)->Write();
      // inclusive_ptave_HF_HLT.at(p).at(r)->Write();
      // inclusive_nevents_HF_HLT.at(p).at(r)->Write();

      f5->cd();
      inclusive_rho_HF_L1HLT.at(p).at(r)->Write();
      inclusive_ptave_HF_L1HLT.at(p).at(r)->Write();
      inclusive_nevents_HF_L1HLT.at(p).at(r)->Write();
    }
  }

  for( int r = 0; r < AlphaBins; r++ ) {
    for( int p = 0; p < PtBins_Central; p++ ) {
      WRITE_HISTOS(SM)
      WRITE_HISTOS(FE_reference)
      WRITE_HISTOS(FE_control)
    }
    for( int p = 0; p < PtBins_HF; p++ ) {
      WRITE_HISTOS(SM_control)
      WRITE_HISTOS(FE)
    }
  }

  f->Close();
  f1->Close();
  f2->Close();
  // f3->Close();
  f_alpha->Close();


  std::vector<TH2F*> h_nevents_central, h_nevents_HF;

  std::vector<double> Pt_bins_Central_D(Pt_bins_Central.begin(), Pt_bins_Central.end());
  std::vector<double> Pt_bins_HF_D(Pt_bins_HF.begin(), Pt_bins_HF.end());

  for (int m = 0; m < 6; m++){
    h_nevents_central.push_back(new TH2F(("central_"+std::to_string(m)).c_str(),("central_"+std::to_string(m)).c_str(),n_eta_bins_JER-1,&eta_bins_JER[0], Pt_bins_Central_D.size()-1,&Pt_bins_Central_D[0]));
    h_nevents_HF.push_back(new TH2F(("HF_"+std::to_string(m)).c_str(),("HF_"+std::to_string(m)).c_str(),n_eta_bins_JER-1,&eta_bins_JER[0], Pt_bins_HF_D.size()-1,&Pt_bins_HF_D[0]));
  }
  std::cout << "Pt_bins_Central: " << nevents_central.size() << std::endl;
  for (size_t i = 0; i < nevents_central.size(); i++) std::cout << "\t" << Pt_bins_Central_D[i];
  std::cout << std::endl;

  for (int r = 0; r < 14; r++) {
    for (int m = 0; m < 6; m++){
      std::cout << r << " " << m << " ";
      for ( int k = 0 ; k < nevents_central.size() ; k++ ) {
        std::cout << "\t" << nevents_central[k][r][m];
        h_nevents_central[m]->SetBinContent(h_nevents_central[m]->GetXaxis()->FindBin(eta_bins_JER[r]), h_nevents_central[m]->GetYaxis()->FindBin(Pt_bins_Central_D.at(k)), nevents_central[k][r][m]);
      }
      std::cout << std::endl;
    } std::cout << std::endl;
  }

  std::cout << "Pt_bins_HF: " << nevents_HF.size() << std::endl;
  for (size_t i = 0; i < nevents_HF.size(); i++) std::cout << "\t" << Pt_bins_HF_D[i];
  std::cout << std::endl;

  for (int r = 0; r < 14; r++) {
    for (int m = 0; m < 6; m++){
      std::cout << r << " " << m << " ";
      for ( int k = 0 ; k < nevents_HF.size() ; k++ ) {
        std::cout << "\t" << nevents_HF[k][r][m];
        h_nevents_HF[m]->SetBinContent(h_nevents_HF[m]->GetXaxis()->FindBin(eta_bins_JER[r]), h_nevents_HF[m]->GetYaxis()->FindBin(Pt_bins_HF_D.at(k)), nevents_HF[k][r][m]);
      }
      std::cout << std::endl;
    }std::cout << std::endl;
  }

  TFile *f_nevents  = new TFile(outdir+"histograms_nevents.root","RECREATE");
  for (int m = 0; m < 6; m++){
    h_nevents_central[m]->Write();
    h_nevents_HF[m]->Write();
  }
  f_nevents->Close();

}

void MySelector::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
