//#include "../include/parameters.h"
#include "../../include/constants.h"
#include "../include/tdrstyle_mod15.h"
#include "../include/CorrectionObject.h"

#include <iostream>
#include <cmath>
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "THStack.h"



using namespace std;


void CorrectionObject::FlavorCorrection_TTree(){
  cout << "--------------- FlavorCorrection_TTree() ---------------" << endl << endl;
  CorrectionObject::make_path(CorrectionObject::_outpath+"/plots/control/Flavor/");
  //  TString path="/nfs/dust/cms/user/karavdia/JERC/FlavorStudy/PFfractions_2017_MadGraph_AK4CHS/";
  //  TString _MCfile="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V6/AK4CHS/MC_NoReweighted_NewHF_Binning/uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_AK4CHS_RunBCDEFGH.root";
  TString jettag=CorrectionObject::_jettag;
  TString txttag=CorrectionObject::_generator_tag; 
  double al_cut=alpha_cut;

  TString JetDescrib;                                                                                                                            
  if (jettag=="AK4PFchs") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";                                                                             
  if (jettag=="AK4PFpuppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";                                                                         
  if (jettag=="AK8PFchs") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";                                                                             
  if (jettag=="AK8PFpuppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 

  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  TTreeReaderValue<Float_t> barreljet_eta(myReader_MC, "barreljet_eta");
  TTreeReaderValue<Float_t> probejet_eta(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> barreljet_pt(myReader_MC, "barreljet_pt");
  TTreeReaderValue<Float_t> probejet_pt(myReader_MC, "probejet_pt");
  TTreeReaderValue<Float_t> probejet_ptgen(myReader_MC, "probejet_ptgen");
  //  TTreeReaderValue<Int_t> flavorBarreljet(myReader_MC, "flavorBarreljet");
  TTreeReaderValue<Int_t> flavorProbejet(myReader_MC, "flavorProbejet");
  //  TTreeReaderValue<Float_t> leadingjet_response(myReader_MC, "leadingjet_response");
  TTreeReaderValue<Float_t> had_n_Efrac(myReader_MC, "probejet_neutHadEF");
  TTreeReaderValue<Float_t> had_ch_Efrac(myReader_MC, "probejet_chHadEF");
  TTreeReaderValue<Float_t> em_n_Efrac(myReader_MC, "probejet_neutEmEF");
  TTreeReaderValue<Float_t> em_ch_Efrac(myReader_MC, "probejet_chEmEF");
  TTreeReaderValue<Float_t> mu_Efrac(myReader_MC, "probejet_muonEF");
  TTreeReaderValue<Float_t> ph_Efrac(myReader_MC, "probejet_photonEF");

  TTreeReaderValue<Float_t> alpha(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> weight(myReader_MC, "weight");



  int n_leadingjet_u = 0, n_leadingjet_d = 0, n_leadingjet_s = 0, n_leadingjet_c = 0, n_leadingjet_b = 0, n_leadingjet_g = 0, n_leadingjet_unm = 0;

  TH1D* energy_fractions_u[5][n_eta-1]; TH1D* energy_fractions_d[5][n_eta-1]; TH1D* energy_fractions_s[5][n_eta-1]; 
  TH1D*  energy_fractions_c[5][n_eta-1]; TH1D* energy_fractions_b[5][n_eta-1]; TH1D* energy_fractions_g[5][n_eta-1]; TH1D* energy_fractions_unm[5][n_eta-1];

  // TH1D* energy_fractions_u_pt[5][n_pt-1][n_eta-1]; TH1D* energy_fractions_d_pt[5][n_pt-1][n_eta-1]; TH1D* energy_fractions_s_pt[5][n_pt-1][n_eta-1]; 
  // TH1D*  energy_fractions_c_pt[5][n_pt-1][n_eta-1]; TH1D* energy_fractions_b_pt[5][n_pt-1][n_eta-1]; TH1D* energy_fractions_g_pt[5][n_pt-1][n_eta-1]; TH1D* energy_fractions_unm_pt[5][n_pt-1][n_eta-1];
  for(int i=0; i<5; i++){
    for(int j=0; j<n_eta-1; j++){
      TString add = "_"+i;      add+="_";      add+=j;
      energy_fractions_u[i][j] = new TH1D("energy_fractions_u_"+add,"",100,0,1);
      energy_fractions_d[i][j] = new TH1D("energy_fractions_d_"+add,"",100,0,1);
      energy_fractions_s[i][j] = new TH1D("energy_fractions_s_"+add,"",100,0,1);
      energy_fractions_c[i][j] = new TH1D("energy_fractions_c_"+add,"",100,0,1);
      energy_fractions_b[i][j] = new TH1D("energy_fractions_b_"+add,"",100,0,1);
      energy_fractions_g[i][j] = new TH1D("energy_fractions_g_"+add,"",100,0,1);
      energy_fractions_unm[i][j] = new TH1D("energy_fractions_unm_"+add,"",100,0,1);
      // for(int k=0; k<n_pt-1; k++){
      // 	add +="_"; add +=k;
      // 	energy_fractions_u_pt[i][k][j] = new TH1D("energy_fractions_u_"+add,"",100,0,1);
      // 	energy_fractions_d_pt[i][k][j] = new TH1D("energy_fractions_d_"+add,"",100,0,1);
      // 	energy_fractions_s_pt[i][k][j] = new TH1D("energy_fractions_s_"+add,"",100,0,1);
      // 	energy_fractions_c_pt[i][k][j] = new TH1D("energy_fractions_c_"+add,"",100,0,1);
      // 	energy_fractions_b_pt[i][k][j] = new TH1D("energy_fractions_b_"+add,"",100,0,1);
      // 	energy_fractions_g_pt[i][k][j] = new TH1D("energy_fractions_g_"+add,"",100,0,1);
      // 	energy_fractions_unm_pt[i][k][j] = new TH1D("energy_fractions_unm_"+add,"",100,0,1);
      // }
    }
  }

  TH1D* response_u[n_eta]; TH1D* response_d[n_eta]; TH1D* response_s[n_eta]; //pt_rec/pt_gen
  TH1D* response_c[n_eta]; TH1D* response_b[n_eta]; TH1D* response_g[n_eta]; 
  TH1D* response_unm[n_eta];
  for(int j=0; j<n_eta-1; j++){
    TString add = "_"+j;
    response_u[j] = new TH1D("response_u_"+add,"",100,0,10);
    response_d[j] = new TH1D("response_d_"+add,"",100,0,10);
    response_s[j] = new TH1D("response_s_"+add,"",100,0,10);
    response_c[j] = new TH1D("response_c_"+add,"",100,0,10);
    response_b[j] = new TH1D("response_b_"+add,"",100,0,10);
    response_g[j] = new TH1D("response_g_"+add,"",100,0,10);
    response_unm[j] = new TH1D("response_unm_"+add,"",100,0,10);
  }
  //  int test_ev = 0;
  while (myReader_MC.Next()) {
    // while (test_ev<100) {
    // test_ev++;
    // myReader_MC.Next();

    if(*alpha>al_cut) continue;
        
    double flavor_leadingjet = 0;
    double eta_leadingjet = -999;


    //response plots
    //for each flavor, calculate mean response of the leading jet
    for(int j=0; j<n_eta-1; j++){
      //      if(0.5*((*probejet_pt)+(*barreljet_pt))*cosh(*probejet_eta)<3250) continue; //Select only bad events
      //      if((*probejet_pt)>50) continue;
      //      if((*probejet_pt)>150 || (*probejet_pt)<50) continue;
      //      if((*probejet_pt)<150) continue;
      //      if(((*barreljet_pt)+(*probejet_pt))<250) continue;
      //if((*barreljet_pt)>50) continue;
      //if((*barreljet_pt)>150 || (*barreljet_pt)<50) continue;
      if(fabs(*probejet_eta)>eta_bins[j+1] || fabs(*probejet_eta)<eta_bins[j]) continue;
      //      cout<<"probejet_pt = "<<(*probejet_pt)<<endl;
      if(*flavorProbejet == 1) {
	energy_fractions_u[0][j] ->Fill( *had_n_Efrac,*weight);
	energy_fractions_u[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_u[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_u[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_u[4][j] ->Fill(  *em_ch_Efrac,*weight);
	response_u[j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
      }
      else if(*flavorProbejet == 2){
	energy_fractions_d[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_d[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_d[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_d[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_d[4][j] ->Fill(  *em_ch_Efrac,*weight);
	response_d[j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
      }
      else if(*flavorProbejet == 3){
	energy_fractions_s[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_s[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_s[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_s[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_s[4][j] ->Fill(  *em_ch_Efrac,*weight);
	response_s[j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
      }
      else if(*flavorProbejet == 4){
	energy_fractions_c[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_c[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_c[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_c[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_c[4][j] ->Fill(  *em_ch_Efrac,*weight);
	response_c[j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
      }
      else if(*flavorProbejet == 5){
	energy_fractions_b[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_b[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_b[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_b[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_b[4][j] ->Fill(  *em_ch_Efrac,*weight);
	response_b[j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
      }
      else if(*flavorProbejet == 21){
	energy_fractions_g[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_g[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_g[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_g[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_g[4][j] ->Fill(  *em_ch_Efrac,*weight);
	response_g[j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
      }
      else if(*flavorProbejet == -1){
	energy_fractions_unm[0][j] ->Fill(  *had_n_Efrac,*weight);
	energy_fractions_unm[1][j] ->Fill(  *had_ch_Efrac,*weight);
	energy_fractions_unm[2][j] ->Fill(  *mu_Efrac,*weight);
	energy_fractions_unm[3][j] ->Fill(  *ph_Efrac,*weight);
	energy_fractions_unm[4][j] ->Fill(  *em_ch_Efrac,*weight);
	response_unm[j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
      }
      else {
	cout << "jetflavor: " << *flavorProbejet << endl;
	throw runtime_error("flavor of jet not known.");
      }
    }  
  } //event loop

  //scale Energy_fractions to stack to unity 
  for(int j=0; j<n_eta-1; j++){
    TH1D* h_neutral_hadron_Efraction = new TH1D("h_neutral_hadron_Efraction","Neutral Hadron energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_charged_hadron_Efraction = new TH1D("h_charged_hadron_Efraction","Charged Hadron energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_neutral_em_Efraction = new TH1D("h_neutral_em_Efraction","Neutral electromagnetic energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_charged_em_Efraction = new TH1D("h_charged_em_Efraction","Charged electromagnetic energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_photon_Efraction = new TH1D("h_photon_Efraction","Photon energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_muon_Efraction = new TH1D("h_muon_Efraction","Muon energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    double sum_fractions_u = 0, sum_fractions_d = 0, sum_fractions_s = 0, sum_fractions_c = 0, sum_fractions_b = 0, sum_fractions_g = 0, sum_fractions_unm = 0;

    TH1D* h_neutral_hadron_Efraction_relud = new TH1D("h_neutral_hadron_Efraction_relud","Neutral Hadron energy fraction (relative to u);Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_charged_hadron_Efraction_relud = new TH1D("h_charged_hadron_Efraction_relud","Charged Hadron energy fraction (relative to u);Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_neutral_em_Efraction_relud = new TH1D("h_neutral_em_Efraction_relud","Neutral electromagnetic energy fraction (relative to u);Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_charged_em_Efraction_relud = new TH1D("h_charged_em_Efraction_relud","Charged electromagnetic energy fraction (relative to u);Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_photon_Efraction_relud = new TH1D("h_photon_Efraction_relud","Photon energy fraction (relative to u);Flavor (Physics);Probe jet PF energy fraction",7,0,7);
    TH1D* h_muon_Efraction_relud = new TH1D("h_muon_Efraction_relud","Muon energy fraction (relative to u);Flavor (Physics);Probe jet PF energy fraction",7,0,7);

    TH1D* h_response = new TH1D("h_response","Response;Flavor (Physics);pt_{REC}/pt_{GEN}",7,0,7);
    h_response->GetXaxis()->SetBinLabel(1,"u");
    h_response->GetXaxis()->SetBinLabel(2,"d");
    h_response->GetXaxis()->SetBinLabel(3,"g");
    h_response->GetXaxis()->SetBinLabel(4,"c");
    h_response->GetXaxis()->SetBinLabel(5,"b");
    h_response->GetXaxis()->SetBinLabel(6,"s");
    h_response->GetXaxis()->SetBinLabel(7,"unmatch");

    for(int i=0; i<5; i++) {
      sum_fractions_u += energy_fractions_u[i][j]->GetMean();
      sum_fractions_d += energy_fractions_d[i][j]->GetMean();
      sum_fractions_s += energy_fractions_s[i][j]->GetMean();
      sum_fractions_c += energy_fractions_c[i][j]->GetMean();
      sum_fractions_b += energy_fractions_b[i][j]->GetMean();
      sum_fractions_g += energy_fractions_g[i][j]->GetMean();
      sum_fractions_unm += energy_fractions_unm[i][j]->GetMean();
    }

    h_response->Fill("u",response_u[j]->GetMean());
    h_response->Fill("d",response_d[j]->GetMean());
    h_response->Fill("g",response_g[j]->GetMean());
    h_response->Fill("c",response_c[j]->GetMean());
    h_response->Fill("b",response_b[j]->GetMean());
    h_response->Fill("s",response_s[j]->GetMean());
    h_response->Fill("unmatch",response_unm[j]->GetMean()); 
  
    h_response->SetBinError(1,response_u[j]->GetMeanError());
    h_response->SetBinError(2,response_d[j]->GetMeanError());
    h_response->SetBinError(3,response_g[j]->GetMeanError());
    h_response->SetBinError(4,response_c[j]->GetMeanError());
    h_response->SetBinError(5,response_b[j]->GetMeanError());
    h_response->SetBinError(6,response_s[j]->GetMeanError());
    h_response->SetBinError(7,response_unm[j]->GetMeanError()); 

    // for(int i=0; i<5; i++) { //TEST; normalisation to flavor response
    //   sum_fractions_u *= response_u[j]->GetMean();
    //   sum_fractions_d *= response_d[j]->GetMean();
    //   sum_fractions_s *= response_s[j]->GetMean();
    //   sum_fractions_c *= response_c[j]->GetMean();
    //   sum_fractions_b *= response_b[j]->GetMean();
    //   sum_fractions_g *= response_g[j]->GetMean();
    //   sum_fractions_unm *= response_unm[j]->GetMean();
    // }


    for(int i=0; i<5; i++) { //0: neutral had, 1: charged had, 2: muon, 3: photon, 4: charged em
      if(i==0){ 
	h_neutral_hadron_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
	h_neutral_hadron_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
	h_neutral_hadron_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
	h_neutral_hadron_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
	h_neutral_hadron_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
	h_neutral_hadron_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
	h_neutral_hadron_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 

	h_neutral_hadron_Efraction->SetBinError(1,(energy_fractions_u[i][j]->GetMeanError())/sum_fractions_u);
	h_neutral_hadron_Efraction->SetBinError(2,(energy_fractions_d[i][j]->GetMeanError())/sum_fractions_d);
	h_neutral_hadron_Efraction->SetBinError(3,(energy_fractions_g[i][j]->GetMeanError())/sum_fractions_g);
	h_neutral_hadron_Efraction->SetBinError(4,(energy_fractions_c[i][j]->GetMeanError())/sum_fractions_c);
	h_neutral_hadron_Efraction->SetBinError(5,(energy_fractions_b[i][j]->GetMeanError())/sum_fractions_b);
	h_neutral_hadron_Efraction->SetBinError(6,(energy_fractions_s[i][j]->GetMeanError())/sum_fractions_s);
	h_neutral_hadron_Efraction->SetBinError(7,(energy_fractions_unm[i][j]->GetMeanError())/sum_fractions_unm); 

	h_neutral_hadron_Efraction_relud = (TH1D*)h_neutral_hadron_Efraction->Clone();
	double value_ud = 0.5*(h_neutral_hadron_Efraction->GetBinContent(2)+h_neutral_hadron_Efraction->GetBinContent(1));
	double value_ud_err = TMath::Hypot(h_neutral_hadron_Efraction->GetBinError(2),h_neutral_hadron_Efraction->GetBinError(1));
	for(int k=1;k<8;k++){
	  double value = 0; double value_err = 0;
	  if(value_ud>0){
	    value = (h_neutral_hadron_Efraction->GetBinContent(k)-value_ud)/value_ud;
	    value_err = fabs((h_neutral_hadron_Efraction->GetBinContent(k)-value_ud)/value_ud)*(TMath::Hypot(h_neutral_hadron_Efraction->GetBinError(k)/h_neutral_hadron_Efraction->GetBinContent(k),value_ud_err/value_ud));
	  }
	  //	  cout<<"value = "<<value<<" h_neutral_hadron_Efraction->GetBinContent(k) = "<<h_neutral_hadron_Efraction->GetBinContent(k)<<" value_ud = "<<value_ud<<endl;
	  h_neutral_hadron_Efraction_relud->SetBinContent(k,value);
	  h_neutral_hadron_Efraction_relud->SetBinError(k,value_err);
	}
      }
      else if(i==1){ 
	h_charged_hadron_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
	h_charged_hadron_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
	h_charged_hadron_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
	h_charged_hadron_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
	h_charged_hadron_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
	h_charged_hadron_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
	h_charged_hadron_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 

	h_charged_hadron_Efraction->SetBinError(1,(energy_fractions_u[i][j]->GetMeanError())/sum_fractions_u);
	h_charged_hadron_Efraction->SetBinError(2,(energy_fractions_d[i][j]->GetMeanError())/sum_fractions_d);
	h_charged_hadron_Efraction->SetBinError(3,(energy_fractions_g[i][j]->GetMeanError())/sum_fractions_g);
	h_charged_hadron_Efraction->SetBinError(4,(energy_fractions_c[i][j]->GetMeanError())/sum_fractions_c);
	h_charged_hadron_Efraction->SetBinError(5,(energy_fractions_b[i][j]->GetMeanError())/sum_fractions_b);
	h_charged_hadron_Efraction->SetBinError(6,(energy_fractions_s[i][j]->GetMeanError())/sum_fractions_s);
	h_charged_hadron_Efraction->SetBinError(7,(energy_fractions_unm[i][j]->GetMeanError())/sum_fractions_unm); 

	h_charged_hadron_Efraction_relud = (TH1D*)h_charged_hadron_Efraction->Clone();
	double value_ud = 0.5*(h_charged_hadron_Efraction->GetBinContent(2)+h_charged_hadron_Efraction->GetBinContent(1));
	double value_ud_err = TMath::Hypot(h_charged_hadron_Efraction->GetBinError(2),h_charged_hadron_Efraction->GetBinError(1));
	for(int k=1;k<8;k++){
	  double value = 0; double value_err = 0;
	  if(value_ud>0){
	    value = (h_charged_hadron_Efraction->GetBinContent(k)-value_ud)/value_ud;
	    value_err = (fabs(h_charged_hadron_Efraction->GetBinContent(k)-value_ud)/value_ud)*(TMath::Hypot(h_charged_hadron_Efraction->GetBinError(k)/h_charged_hadron_Efraction->GetBinContent(k),value_ud_err/value_ud ));
	  }
	  //	  cout<<"value = "<<value<<" h_charged_hadron_Efraction->GetBinContent(k) = "<<h_charged_hadron_Efraction->GetBinContent(k)<<" value_ud = "<<value_ud<<endl;
	  h_charged_hadron_Efraction_relud->SetBinContent(k,value);
	  h_charged_hadron_Efraction_relud->SetBinError(k,value_err);
	}
      }
      else if(i==2){ 
	h_muon_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
	h_muon_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
	h_muon_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
	h_muon_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
	h_muon_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
	h_muon_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
	h_muon_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 

	h_muon_Efraction->SetBinError(1,(energy_fractions_u[i][j]->GetMeanError())/sum_fractions_u);
	h_muon_Efraction->SetBinError(2,(energy_fractions_d[i][j]->GetMeanError())/sum_fractions_d);
	h_muon_Efraction->SetBinError(3,(energy_fractions_g[i][j]->GetMeanError())/sum_fractions_g);
	h_muon_Efraction->SetBinError(4,(energy_fractions_c[i][j]->GetMeanError())/sum_fractions_c);
	h_muon_Efraction->SetBinError(5,(energy_fractions_b[i][j]->GetMeanError())/sum_fractions_b);
	h_muon_Efraction->SetBinError(6,(energy_fractions_s[i][j]->GetMeanError())/sum_fractions_s);
	h_muon_Efraction->SetBinError(7,(energy_fractions_unm[i][j]->GetMeanError())/sum_fractions_unm); 

	h_muon_Efraction_relud = (TH1D*)h_muon_Efraction->Clone();
	double value_ud = 0.5*(h_muon_Efraction->GetBinContent(2)+h_muon_Efraction->GetBinContent(1));
	double value_ud_err = TMath::Hypot(h_muon_Efraction->GetBinError(2),h_muon_Efraction->GetBinError(1));
	for(int k=1;k<8;k++){
	  double value = 0; double value_err = 0;
	  if(value_ud>0){
	    value = (h_muon_Efraction->GetBinContent(k)-value_ud)/value_ud;
	    value_err = fabs((h_muon_Efraction->GetBinContent(k)-value_ud)/value_ud)*(TMath::Hypot(h_muon_Efraction->GetBinError(k)/h_muon_Efraction->GetBinContent(k),value_ud_err/value_ud ));
	  }
	  //	  cout<<"value = "<<value<<" h_muon_Efraction->GetBinContent(k) = "<<h_muon_Efraction->GetBinContent(k)<<" value_ud = "<<value_ud<<endl;
	  h_muon_Efraction_relud->SetBinContent(k,value);
	  h_muon_Efraction_relud->SetBinError(k,value_err);
	}
      }
      else if(i==3){ 
	h_photon_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
	h_photon_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
	h_photon_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
	h_photon_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
	h_photon_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
	h_photon_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
	h_photon_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 

	h_photon_Efraction->SetBinError(1,(energy_fractions_u[i][j]->GetMeanError())/sum_fractions_u);
	h_photon_Efraction->SetBinError(2,(energy_fractions_d[i][j]->GetMeanError())/sum_fractions_d);
	h_photon_Efraction->SetBinError(3,(energy_fractions_g[i][j]->GetMeanError())/sum_fractions_g);
	h_photon_Efraction->SetBinError(4,(energy_fractions_c[i][j]->GetMeanError())/sum_fractions_c);
	h_photon_Efraction->SetBinError(5,(energy_fractions_b[i][j]->GetMeanError())/sum_fractions_b);
	h_photon_Efraction->SetBinError(6,(energy_fractions_s[i][j]->GetMeanError())/sum_fractions_s);
	h_photon_Efraction->SetBinError(7,(energy_fractions_unm[i][j]->GetMeanError())/sum_fractions_unm); 

	h_photon_Efraction_relud = (TH1D*)h_photon_Efraction->Clone();
	double value_ud = 0.5*(h_photon_Efraction->GetBinContent(2)+h_photon_Efraction->GetBinContent(1));
	double value_ud_err = TMath::Hypot(h_photon_Efraction->GetBinError(2),h_photon_Efraction->GetBinError(1));
	for(int k=1;k<8;k++){
	  double value = 0; double value_err = 0;
	  if(value_ud>0){
	    value = (h_photon_Efraction->GetBinContent(k)-value_ud)/value_ud;
	    value_err = fabs((h_photon_Efraction->GetBinContent(k)-value_ud)/value_ud)*(TMath::Hypot(h_photon_Efraction->GetBinError(k)/h_photon_Efraction->GetBinContent(k),value_ud_err/value_ud ));
	  }
	  //	  cout<<"value = "<<value<<" h_photon_Efraction->GetBinContent(k) = "<<h_photon_Efraction->GetBinContent(k)<<" value_ud = "<<value_ud<<endl;
	  h_photon_Efraction_relud->SetBinContent(k,value);
	  h_photon_Efraction_relud->SetBinError(k,value_err);
	}
      }
      else if(i==4){ 
	h_charged_em_Efraction->Fill("u",(energy_fractions_u[i][j]->GetMean())/sum_fractions_u);
	h_charged_em_Efraction->Fill("d",(energy_fractions_d[i][j]->GetMean())/sum_fractions_d);
	h_charged_em_Efraction->Fill("g",(energy_fractions_g[i][j]->GetMean())/sum_fractions_g);
	h_charged_em_Efraction->Fill("c",(energy_fractions_c[i][j]->GetMean())/sum_fractions_c);
	h_charged_em_Efraction->Fill("b",(energy_fractions_b[i][j]->GetMean())/sum_fractions_b);
	h_charged_em_Efraction->Fill("s",(energy_fractions_s[i][j]->GetMean())/sum_fractions_s);
	h_charged_em_Efraction->Fill("unmatch",(energy_fractions_unm[i][j]->GetMean())/sum_fractions_unm); 

	h_charged_em_Efraction->SetBinError(1,(energy_fractions_u[i][j]->GetMeanError())/sum_fractions_u);
	h_charged_em_Efraction->SetBinError(2,(energy_fractions_d[i][j]->GetMeanError())/sum_fractions_d);
	h_charged_em_Efraction->SetBinError(3,(energy_fractions_g[i][j]->GetMeanError())/sum_fractions_g);
	h_charged_em_Efraction->SetBinError(4,(energy_fractions_c[i][j]->GetMeanError())/sum_fractions_c);
	h_charged_em_Efraction->SetBinError(5,(energy_fractions_b[i][j]->GetMeanError())/sum_fractions_b);
	h_charged_em_Efraction->SetBinError(6,(energy_fractions_s[i][j]->GetMeanError())/sum_fractions_s);
	h_charged_em_Efraction->SetBinError(7,(energy_fractions_unm[i][j]->GetMeanError())/sum_fractions_unm); 

	h_charged_em_Efraction_relud = (TH1D*)h_charged_em_Efraction->Clone();
	double value_ud = 0.5*(h_charged_em_Efraction->GetBinContent(2)+h_charged_em_Efraction->GetBinContent(1));
	double value_ud_err = TMath::Hypot(h_charged_em_Efraction->GetBinError(2),h_charged_em_Efraction->GetBinError(1));
	for(int k=1;k<8;k++){
	  double value = 0; double value_err = 0;
	  if(value_ud>0){
	    value = (h_charged_em_Efraction->GetBinContent(k)-value_ud)/value_ud;
	    value_err = fabs((h_charged_em_Efraction->GetBinContent(k)-value_ud)/value_ud)*(TMath::Hypot(h_charged_em_Efraction->GetBinError(k)/h_charged_em_Efraction->GetBinContent(k),value_ud_err/value_ud ));
	  }
	  //	  cout<<"value = "<<value<<" h_charged_em_Efraction->GetBinContent(k) = "<<h_charged_em_Efraction->GetBinContent(k)<<" value_ud = "<<value_ud<<endl;
	  h_charged_em_Efraction_relud->SetBinContent(k,value);
	  h_charged_em_Efraction_relud->SetBinError(k,value_err);
	}
      }
    }

    // //main plot
    TString alVal;
    alVal.Form("%0.2f\n",al_cut);
    TString altitle = " #alpha<"+alVal;

    TString etaVal;
    etaVal.Form("%0.2f\n",eta_bins[n_etabarr]);
 
    TString Latextext1 = JetDescrib;
    TString Latextext2 = "QCD dijet, " + eta_range[j] + " < |#eta| < " + eta_range[j+1]  + altitle;
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045);  
 
 
    //TCanvas* c_stack =  tdrCanvas("c_stack", h_neutral_hadron_Efraction, 4, 10, kSquare);
    TCanvas* c_stack =  new TCanvas("c_stack","c_stack",1);
    // draw stack
    THStack* stack = new THStack("stack",";Flavor;Probe jet energy fraction");
    h_neutral_hadron_Efraction->SetLineColor(8);
    h_neutral_hadron_Efraction->SetFillColor(8);
    h_charged_hadron_Efraction->SetLineColor(7);
    h_charged_hadron_Efraction->SetFillColor(7);
    h_muon_Efraction->SetLineColor(46);
    h_muon_Efraction->SetFillColor(46);
    h_photon_Efraction->SetLineColor(41);
    h_photon_Efraction->SetFillColor(41);
    h_charged_em_Efraction->SetLineColor(15);
    h_charged_em_Efraction->SetFillColor(15);
    h_charged_hadron_Efraction->GetYaxis()->SetRangeUser(0,1.1);
    stack->Add(h_charged_hadron_Efraction);
    stack->Add(h_photon_Efraction);
    stack->Add(h_neutral_hadron_Efraction);
    stack->Add(h_charged_em_Efraction);
    stack->Add(h_muon_Efraction);

    stack->Draw("HIST");
    //  stack->Draw("E");

    TLegend* leg2 =  new TLegend(0.70,0.78,0.99,0.99);
    leg2 -> AddEntry(h_muon_Efraction, "#mu");
    leg2 -> AddEntry(h_photon_Efraction, "#gamma");
    leg2 -> AddEntry(h_charged_hadron_Efraction, "Ch. hadron");
    leg2 -> AddEntry(h_neutral_hadron_Efraction, "Neut. hadron");
    leg2 -> AddEntry(h_charged_em_Efraction, "Ch. EM");
    leg2 -> SetFillColor(10);
    leg2 -> SetBorderSize(0);
    leg2 -> SetTextSize(0.042);
    leg2 -> SetLineColor(1);
    leg2 -> SetTextFont(42);
    leg2 -> Draw();
    tex->DrawLatex(0.22,0.57,Latextext1);
    tex->DrawLatex(0.22,0.53,Latextext2);
 
    c_stack->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/PFEnergyFractions_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+".pdf");
    delete stack;
    delete c_stack;


    TCanvas* c_stack_relud =  new TCanvas("c_stack_relud","c_stack_relud",1);
    // draw stack_relud
    //  THStack* stack_relud = new THStack("stack_relud",";Flavor;Relative to ud PF energy fraction: (frac^{i} - frac^{ud})/frac_{ud}");
    THStack* stack_relud = new THStack("stack_relud",";Flavor;(PF frac^{i} - PF frac^{ud})/(PF frac_{ud})");
    h_neutral_hadron_Efraction_relud->SetLineColor(8);
    h_neutral_hadron_Efraction_relud->SetMarkerColor(8);
    h_neutral_hadron_Efraction_relud->SetMarkerStyle(20);
    h_neutral_hadron_Efraction_relud->SetFillColor(8);
    h_charged_hadron_Efraction_relud->SetLineColor(7);
    h_charged_hadron_Efraction_relud->SetMarkerColor(7);
    h_charged_hadron_Efraction_relud->SetMarkerColor(7);
    h_charged_hadron_Efraction_relud->SetMarkerStyle(20);
    h_charged_hadron_Efraction_relud->SetFillColor(7);
    h_muon_Efraction_relud->SetLineColor(46);
    h_muon_Efraction_relud->SetMarkerColor(46);
    h_muon_Efraction_relud->SetMarkerStyle(20);
    h_muon_Efraction_relud->SetFillColor(46);
    h_photon_Efraction_relud->SetLineColor(41);
    h_photon_Efraction_relud->SetMarkerColor(41);
    h_photon_Efraction_relud->SetMarkerStyle(20);
    h_photon_Efraction_relud->SetFillColor(41);
    h_charged_em_Efraction_relud->SetLineColor(15);
    h_charged_em_Efraction_relud->SetMarkerColor(15);
    h_charged_em_Efraction_relud->SetMarkerStyle(20);
    h_charged_em_Efraction_relud->SetFillColor(15);

    h_charged_hadron_Efraction_relud->GetYaxis()->SetRangeUser(-0.35,0.35);
    h_photon_Efraction_relud->GetYaxis()->SetRangeUser(-0.35,0.35);
    h_neutral_hadron_Efraction_relud->GetYaxis()->SetRangeUser(-0.35,0.35);
    stack_relud->Add(h_charged_hadron_Efraction_relud);
    stack_relud->Add(h_photon_Efraction_relud);
    stack_relud->Add(h_neutral_hadron_Efraction_relud);
    //stack_relud->Add(h_charged_em_Efraction_relud);
    //stack_relud->Add(h_muon_Efraction_relud);

    //    stack_relud->GetYaxis()->SetRangeUser(-0.5,0.5);

    stack_relud->Draw("nostack"); //as 'same'
    //  stack_relud->Draw("nostackb"); //draw next to each other
    //  stack->Draw("E");

    TLegend* leg3 =  new TLegend(0.70,0.78,0.99,0.99);
    leg3 -> AddEntry(h_muon_Efraction, "#mu");
    leg3 -> AddEntry(h_photon_Efraction, "#gamma");
    leg3 -> AddEntry(h_charged_hadron_Efraction, "Ch. hadron");
    leg3 -> AddEntry(h_neutral_hadron_Efraction, "Neut. hadron");
    leg3 -> AddEntry(h_charged_em_Efraction, "Ch. EM");
    //  leg3 -> AddEntry(h_neutral_em_Efraction, "Neut. EM");
    leg3 -> SetFillColor(10);
    leg3 -> SetBorderSize(0);
    leg3 -> SetTextSize(0.042);
    leg3 -> SetLineColor(1);
    leg3 -> SetTextFont(42);
    leg3 -> Draw();
    tex->DrawLatex(0.15,0.27,Latextext1);
    tex->DrawLatex(0.15,0.23,Latextext2);
 
    c_stack_relud->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/PFEnergyFractions_RelToUD_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+".pdf");
    delete stack_relud;
    delete c_stack_relud;


    //dummy for tdrCanvas
    TH1D *h = new TH1D("h","Response;Flavor (Physics);pt_{REC}/pt_{GEN}",7,0,7);
    // h->SetMaximum(1.2);
    // h->SetMinimum(0.8);
    h->GetXaxis()->SetBinLabel(1,"u");
    h->GetXaxis()->SetBinLabel(2,"d");
    h->GetXaxis()->SetBinLabel(3,"g");
    h->GetXaxis()->SetBinLabel(4,"c");
    h->GetXaxis()->SetBinLabel(5,"b");
    h->GetXaxis()->SetBinLabel(6,"s");
    h->GetXaxis()->SetBinLabel(7,"unmatch");

    TCanvas* c_response =  new TCanvas("c_response","c_response",1);
    tdrCanvas(c_response,"c_response",h,4,10,true,CorrectionObject::_lumitag);
    h->SetMarkerColor(2);
    h->SetLineColor(2);
    h->SetMarkerStyle(20);
    h->GetYaxis()->SetRangeUser(0.75,1.20);
    h->GetYaxis()->SetTitle("pt_{REC}/pt_{GEN}");
    h->GetXaxis()->SetTitle("Flavor");
    h_response->SetMarkerColor(2);
    h->GetYaxis()->SetRangeUser(0.75,1.20);
    h_response->SetLineColor(2);
    h_response->SetMarkerStyle(20);
    //    h_response->GetYaxis()->SetRangeUser(0.75,1.20);
    h_response->Draw("P SAME");

    // // draw response
    // //  THStack* response = new THStack("response",";Flavor;Relative to ud PF energy fraction: (frac^{i} - frac^{ud})/frac_{ud}");
    // THStack* response = new THStack("response",";Flavor;pt_{REC}/pt_{GEN}");
    // h_response->SetMarkerColor(2);
    // h_response->SetLineColor(2);
    // h_response->SetMarkerStyle(20);
    // h_response->GetYaxis()->SetRangeUser(0.75,1.20);
    // response->Add(h_response);
    // //    response->GetYaxis()->SetRangeUser(0.75,1.20);

    // response->Draw("nostack"); //as 'same'
    // response->GetYaxis()->SetRangeUser(0.75,1.20);
    // //    c_response->Update();
    // //  response->Draw("nostackb"); //draw next to each other
    // //  stack->Draw("E");

    TLegend* leg4 =  new TLegend(0.70,0.88,0.99,0.99);
    leg4 -> AddEntry(h_response, "response");
    leg4 -> SetFillColor(10);
    leg4 -> SetBorderSize(0);
    leg4 -> SetTextSize(0.042);
    leg4 -> SetLineColor(1);
    leg4 -> SetTextFont(42);
    leg4 -> Draw();
    tex->DrawLatex(0.22,0.27,Latextext1);
    tex->DrawLatex(0.22,0.23,Latextext2);
 
    c_response->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/ResponsesFlavor_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+".pdf");
    //    delete response;
    delete c_response;

  }

  cout << "Hello, I'm done" << endl;

}
