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


void CorrectionObject::Flavor_JetPFFractions(){
  cout << "--------------- Flavor_JetPFFractions() ---------------" << endl << endl;
  CorrectionObject::make_path(CorrectionObject::_outpath+"/plots/control/Flavor/");
  TString jettag=CorrectionObject::_jettag;
  TString txttag=CorrectionObject::_generator_tag; 
  double al_cut=alpha_cut;
  int n_pt_ = max(n_pt,n_pt_HF);

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
  bool eta_cut_bool = false;


  //  int n_leadingjet_u = 0, n_leadingjet_d = 0, n_leadingjet_s = 0, n_leadingjet_c = 0, n_leadingjet_b = 0, n_leadingjet_g = 0, n_leadingjet_unm = 0;
  TH1D* energy_fractions[7][5][n_eta-1]; //flavor, jet PF fraction, eta
  TH1D* energy_fractions_pt[7][5][n_pt_-1][n_eta-1];//flavor, jet PF fraction,pt,eta
  TH1D* response[7][n_eta-1]; //jet PF fraction, eta
  TH1D* response_pt[7][n_pt_-1][n_eta-1]; //jet PF fraction, pt, eta
  int pdgID[7]={1,2,3,4,5,21,-1};//u,d,s,c,b,g,unmatched
  TString pdgIDLabels[7]={"u","d","s","c","b","g","unmatch"};//u,d,s,c,b,g,unmatched

  for(int iflv=0; iflv<7; iflv++){
    for(int i=0; i<5; i++){
      for(int j=0; j<n_eta-1; j++){
	TString add = "_"+iflv; add+="_";  add+=i;    add+="_";      add+=j;
	energy_fractions[iflv][i][j] = new TH1D("energy_fractions_"+add,"",100,0,1);
	response[iflv][j] = new TH1D("response_"+add,"",100,0,10);
	for(int k=0; k<n_pt_-1; k++){
	  add +="_"; add +=k;
	  energy_fractions_pt[iflv][i][k][j] = new TH1D("energy_fractions_"+add,"",100,0,1);
	  response_pt[iflv][k][j] = new TH1D("response_"+add,"",100,0,10);
	}
      }
    }
  }


 
  while (myReader_MC.Next()) {
  // int test_ev=0;
  // while (test_ev<100) {
  //   test_ev++;
  //   myReader_MC.Next();
    double pt_binning_var = 0.5*(*barreljet_pt+*probejet_pt);//bin in pt_ave, RECO
    if(pt_binning_var*cosh(*probejet_eta)>3250) continue;
    if(*alpha>al_cut) continue;
        
    double flavor_leadingjet = 0;
    double eta_leadingjet = -999;


    //response plots
    //for each flavor, calculate mean response of the leading jet
    for(int j=0; j<n_eta-1; j++){
      if(fabs(*probejet_eta)>eta_bins[j+1] || fabs(*probejet_eta)<eta_bins[j]) continue;
      for(int ipdg=0;ipdg<7;ipdg++){
	if(*flavorProbejet == pdgID[ipdg]){
	  energy_fractions[ipdg][0][j] ->Fill( *had_n_Efrac,*weight);
	  energy_fractions[ipdg][1][j] ->Fill(  *had_ch_Efrac,*weight);
	  energy_fractions[ipdg][2][j] ->Fill(  *mu_Efrac,*weight);
	  //	  energy_fractions[ipdg][3][j] ->Fill(  *ph_Efrac,*weight);
	  energy_fractions[ipdg][3][j] ->Fill(  *em_n_Efrac,*weight);
	  energy_fractions[ipdg][4][j] ->Fill(  *em_ch_Efrac,*weight);
	  response[ipdg][j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
	  eta_cut_bool = fabs(eta_bins[j])>eta_cut;
	  //	  cout<<"eta_cut_bool = "<<eta_cut_bool<<" fabs(eta_bins[j]) = "<<fabs(eta_bins[j])<<" eta_cut="<<eta_cut<<endl;
	  for(int k= 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){
	    if(pt_binning_var < (eta_cut_bool?pt_bins_HF:pt_bins)[k] || pt_binning_var >= (eta_cut_bool?pt_bins_HF:pt_bins)[k+1]) continue;
	    //	    if(eta_cut_bool && k==3) cout<<" pt_binning_var = "<<pt_binning_var<<endl;
	    energy_fractions_pt[ipdg][0][k][j] ->Fill( *had_n_Efrac,*weight);
	    energy_fractions_pt[ipdg][1][k][j] ->Fill(  *had_ch_Efrac,*weight);
	    energy_fractions_pt[ipdg][2][k][j] ->Fill(  *mu_Efrac,*weight);
	    //	    energy_fractions_pt[ipdg][3][k][j] ->Fill(  *ph_Efrac,*weight);
	    energy_fractions_pt[ipdg][3][k][j] ->Fill(  *em_n_Efrac,*weight);
	    energy_fractions_pt[ipdg][4][k][j] ->Fill(  *em_ch_Efrac,*weight);
	    response_pt[ipdg][k][j] ->Fill((*probejet_pt)/(*probejet_ptgen),*weight);
	    //	    cout<<"pt = "<<pt_binning_var<<" pt_bin = "<<k<<", response = "<<(*probejet_pt)/(*probejet_ptgen)<<endl;
	  }
	}
      }

    }
  } //event loop



  //scale Energy_fractions to stack to unity 
  for(int j=0; j<n_eta-1; j++){
    eta_cut_bool = fabs(eta_bins[j])>eta_cut;
    TH1D* h_PF_Efraction[5]; TH1D* h_PF_Efraction_relud[5];
    double sum_fractions[7];
    TH1D* h_PF_Efraction_pt[5][n_pt_-1]; TH1D* h_PF_Efraction_relud_pt[5][n_pt_-1];
    double sum_fractions_pt[7][n_pt_-1];
    TH1D* h_response = new TH1D("h_response","Response;Flavor (Physics);pt_{REC}/pt_{GEN}",7,0,7);
    TH1D* h_response_pt[n_pt_-1];
    for(int i=0;i<7;i++){
      sum_fractions[i] = 0;
      TString name;// = "h_response_"; name +=i;
      for(int k = 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){ 
	sum_fractions_pt[i][k] = 0;
	name = "h_response_"; name +=i;  name+="_"; name+=j; 
	TString name2 = "_pT_";	name2 += (eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];
	TString name3 = name+name2;
	h_response_pt[k] = new TH1D(name3,"Response;Flavor (Physics);pt_{REC}/pt_{GEN}",7,0,7);
	name = "h_PF_EFraction_"; name+=i; name+="_";  name+=j; 
	name3 = name+name2;
	h_PF_Efraction_pt[i][k] = new TH1D(name3,"PF energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
	name = "h_PF_EFraction_relud_"; name+=i; name+="_"; name+=j; 
	name3 = name+name2;
	h_PF_Efraction_relud_pt[i][k] = new TH1D(name3,";Flavor (Physics);Probe jet PF energy fraction (relative to ud)",7,0,7);
	//	cout<<name3.Data()<<endl;
      }
      name = "h_PF_EFraction"; name+=i; name+="_"; name+=j; 
      h_PF_Efraction[i] = new TH1D(name,"PF energy fraction;Flavor (Physics);Probe jet PF energy fraction",7,0,7);
      name = "h_PF_EFraction_relud"; name+=i; name+="_"; name+=j; 
      h_PF_Efraction_relud[i] = new TH1D(name,";Flavor (Physics);Probe jet PF energy fraction (relative to ud)",7,0,7);
    }

    for(int i=0; i<5; i++) {
        for(int ipdg=0; ipdg<7; ipdg++) {
	h_response->GetXaxis()->SetBinLabel(ipdg+1,pdgIDLabels[ipdg]);
	h_PF_Efraction[i]->GetXaxis()->SetBinLabel(ipdg+1,pdgIDLabels[ipdg]);
	h_PF_Efraction_relud[i]->GetXaxis()->SetBinLabel(ipdg+1,pdgIDLabels[ipdg]);
	sum_fractions[ipdg] += energy_fractions[ipdg][i][j]->GetMean();
	h_response->SetBinContent(ipdg,response[ipdg][j]->GetMean());
	h_response->SetBinError(ipdg,response[ipdg][j]->GetMeanError());
	for(int k = 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){ 
	  h_response_pt[k]->GetXaxis()->SetBinLabel(ipdg+1,pdgIDLabels[ipdg]);
	  h_PF_Efraction_pt[i][k]->GetXaxis()->SetBinLabel(ipdg+1,pdgIDLabels[ipdg]);
	  h_PF_Efraction_relud_pt[i][k]->GetXaxis()->SetBinLabel(ipdg+1,pdgIDLabels[ipdg]);

	  sum_fractions_pt[ipdg][k] += energy_fractions_pt[ipdg][i][k][j]->GetMean();
	  h_response_pt[k]->SetBinContent(ipdg,response_pt[ipdg][k][j]->GetMean());
	  h_response_pt[k]->SetBinError(ipdg,response_pt[ipdg][k][j]->GetMeanError());
	}
      }

    }

    for(int i=0; i<5; i++) { //0: neutral had, 1: charged had, 2: muon, 3: photon, 4: charged em
      for(int ipdg=0;ipdg<7;ipdg++){
	h_PF_Efraction[i]->SetBinContent(ipdg,(energy_fractions[ipdg][i][j]->GetMean())/sum_fractions[ipdg]);
	h_PF_Efraction[i]->SetBinError(ipdg,(energy_fractions[ipdg][i][j]->GetMeanError())/sum_fractions[ipdg]);

	double value_ud = 0.5*(h_PF_Efraction[i]->GetBinContent(2)+h_PF_Efraction[i]->GetBinContent(1));
	double value_ud_err = TMath::Hypot(h_PF_Efraction[i]->GetBinError(2),h_PF_Efraction[i]->GetBinError(1));
	for(int kbin=1;kbin<8;kbin++){
	  double value = 0; double value_err = 0;
	  if(value_ud>0){
	    value = (h_PF_Efraction[i]->GetBinContent(kbin)-value_ud)/value_ud;
	    value_err = fabs((h_PF_Efraction[i]->GetBinContent(kbin)-value_ud)/value_ud)*(TMath::Hypot(h_PF_Efraction[i]->GetBinError(kbin)/h_PF_Efraction[i]->GetBinContent(kbin),value_ud_err/value_ud));
	  }
	  h_PF_Efraction_relud[i]->SetBinContent(kbin,value);
	  h_PF_Efraction_relud[i]->SetBinError(kbin,value_err);
	}

	for(int k = 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){ 
	  h_PF_Efraction_pt[i][k]->SetBinContent(ipdg,(energy_fractions_pt[ipdg][i][k][j]->GetMean())/sum_fractions_pt[ipdg][k]);
	  h_PF_Efraction_pt[i][k]->SetBinError(ipdg,(energy_fractions_pt[ipdg][i][k][j]->GetMeanError())/sum_fractions_pt[ipdg][k]);
	  value_ud = 0.5*(h_PF_Efraction_pt[i][k]->GetBinContent(2)+h_PF_Efraction_pt[i][k]->GetBinContent(1));
	  value_ud_err = TMath::Hypot(h_PF_Efraction_pt[i][k]->GetBinError(2),h_PF_Efraction_pt[i][k]->GetBinError(1));
	  for(int kbin=1;kbin<8;kbin++){
	    double value = 0; double value_err = 0;
	    if(value_ud>0){
	      value = (h_PF_Efraction_pt[i][k]->GetBinContent(kbin)-value_ud)/value_ud;
	      value_err = fabs((h_PF_Efraction_pt[i][k]->GetBinContent(kbin)-value_ud)/value_ud)*(TMath::Hypot(h_PF_Efraction_pt[i][k]->GetBinError(kbin)/h_PF_Efraction_pt[i][k]->GetBinContent(kbin),value_ud_err/value_ud));
	    }
	    h_PF_Efraction_relud_pt[i][k]->SetBinContent(kbin,value);
	    h_PF_Efraction_relud_pt[i][k]->SetBinError(kbin,value_err);
	  }
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
 
    int colorSchemePF[5]={8,7,46,41,15};//neutral_hadron,charged_hadron,muon,photon,charged_em
    //    TString LabeleSchemePF[5]={"Neut. hadron","Ch. hadron","#mu","#gamma","Ch. EM"};
    TString LabeleSchemePF[5]={"Neut. HAD","Ch. HAD","#mu","Neut. EM","Ch. EM"};

    //Integrated over pt ----------------
    TCanvas* c_stack =  new TCanvas("c_stack","c_stack",1);
    // draw stack
    THStack* stack = new THStack("stack",";Flavor;Probe jet energy fraction");
    THStack* stack_relud = new THStack("stack_relud",";Flavor;(PF frac^{i} - PF frac^{ud})/(PF frac_{ud})");

    TLegend* leg2 =  new TLegend(0.70,0.78,0.99,0.99);
    leg2 -> SetFillColor(10);
    leg2 -> SetBorderSize(0);
    leg2 -> SetTextSize(0.042);
    leg2 -> SetLineColor(1);
    leg2 -> SetTextFont(42);
    for(int iPF=0;iPF<5;iPF++){
      h_PF_Efraction[iPF]->SetLineColor(colorSchemePF[iPF]);
      h_PF_Efraction[iPF]->SetFillColor(colorSchemePF[iPF]);
      h_PF_Efraction[iPF]->SetMarkerColor(colorSchemePF[iPF]);
      h_PF_Efraction[iPF]->SetMarkerStyle(20);
      h_PF_Efraction_relud[iPF]->SetLineColor(colorSchemePF[iPF]);
      h_PF_Efraction_relud[iPF]->SetFillColor(colorSchemePF[iPF]);
      h_PF_Efraction_relud[iPF]->SetMarkerColor(colorSchemePF[iPF]);
      h_PF_Efraction_relud[iPF]->SetMarkerStyle(20);

      h_PF_Efraction_relud[iPF]->GetYaxis()->SetRangeUser(-0.35,0.35);
      h_PF_Efraction[iPF]->GetYaxis()->SetRangeUser(0,1.1);

      stack->Add(h_PF_Efraction[iPF]);
      stack_relud->Add(h_PF_Efraction_relud[iPF]);

      leg2 -> AddEntry(h_PF_Efraction[iPF],LabeleSchemePF[iPF]);
    }
   
    stack->Draw("HIST");
    leg2 -> Draw();
    tex->DrawLatex(0.22,0.57,Latextext1);
    tex->DrawLatex(0.22,0.53,Latextext2);
    tex->DrawLatex(0.15,0.95,txttag);

    c_stack->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/PFEnergyFractions_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+".pdf");
    delete stack;
    delete c_stack;


    TCanvas* c_stack_relud =  new TCanvas("c_stack_relud","c_stack_relud",1);
    stack_relud->Draw("nostack"); //as 'same'
    leg2 -> Draw();
    tex->DrawLatex(0.15,0.95,txttag);
    tex->DrawLatex(0.15,0.27,Latextext1);
    tex->DrawLatex(0.15,0.23,Latextext2);
 
    c_stack_relud->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/PFEnergyFractions_RelToUD_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+".pdf");
    delete stack_relud;
    delete c_stack_relud;


    //dummy for tdrCanvas
    TH1D *h = new TH1D("h","Response;Flavor (Physics);pt_{REC}/pt_{GEN}",7,0,7);
    for(int ipdg=0;ipdg<7;ipdg++){
      h->GetXaxis()->SetBinLabel(ipdg+1,pdgIDLabels[ipdg]);
    }
   

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
    h_response->Draw("P SAME");

    TLegend* leg4 =  new TLegend(0.60,0.68,0.70,0.75);
    leg4 -> AddEntry(h_response, "response");
    leg4 -> SetFillColor(10);
    leg4 -> SetBorderSize(0);
    leg4 -> SetTextSize(0.042);
    leg4 -> SetLineColor(1);
    leg4 -> SetTextFont(42);
    leg4 -> Draw();
    tex->DrawLatex(0.22,0.27,Latextext1);
    tex->DrawLatex(0.22,0.23,Latextext2);
    tex->DrawLatex(0.15,0.95,txttag);

    c_response->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/ResponsesFlavor_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+".pdf");
    delete c_response;
    //[END]Integrated over pt ----------------

    //in pt bins ----------------------------
    for(int k = 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){ 
      TString Latextext3 = " " + (eta_cut_bool?pt_range_HF:pt_range)[k]+"#leq p_{T}<"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];

      TCanvas* c_stack =  new TCanvas("c_stack","c_stack",1);
      // draw stack
      THStack* stack = new THStack("stack",";Flavor;Probe jet energy fraction");
      THStack* stack_relud = new THStack("stack_relud",";Flavor;(PF frac^{i} - PF frac^{ud})/(PF frac_{ud})");

      TLegend* leg2 =  new TLegend(0.70,0.78,0.99,0.99);
      leg2 -> SetFillColor(10);
      leg2 -> SetBorderSize(0);
      leg2 -> SetTextSize(0.042);
      leg2 -> SetLineColor(1);
      leg2 -> SetTextFont(42);
      for(int iPF=0;iPF<5;iPF++){
	h_PF_Efraction_pt[iPF][k]->SetLineColor(colorSchemePF[iPF]);
	h_PF_Efraction_pt[iPF][k]->SetFillColor(colorSchemePF[iPF]);
	h_PF_Efraction_pt[iPF][k]->SetMarkerColor(colorSchemePF[iPF]);
	h_PF_Efraction_pt[iPF][k]->SetMarkerStyle(20);
	h_PF_Efraction_relud_pt[iPF][k]->SetLineColor(colorSchemePF[iPF]);
	h_PF_Efraction_relud_pt[iPF][k]->SetFillColor(colorSchemePF[iPF]);
	h_PF_Efraction_relud_pt[iPF][k]->SetMarkerColor(colorSchemePF[iPF]);
	h_PF_Efraction_relud_pt[iPF][k]->SetMarkerStyle(20);
	
	h_PF_Efraction_relud_pt[iPF][k]->GetYaxis()->SetRangeUser(-0.35,0.35);
	h_PF_Efraction_pt[iPF][k]->GetYaxis()->SetRangeUser(0,1.1);

	stack->Add(h_PF_Efraction_pt[iPF][k]);
	stack_relud->Add(h_PF_Efraction_relud_pt[iPF][k]);
	
	leg2 -> AddEntry(h_PF_Efraction_pt[iPF][k],LabeleSchemePF[iPF]);
      }
      
      stack->Draw("HIST");
      leg2 -> Draw();
      tex->DrawLatex(0.22,0.57,Latextext1);
      tex->DrawLatex(0.22,0.53,Latextext2);
      tex->DrawLatex(0.22,0.49,Latextext3);
      tex->DrawLatex(0.15,0.95,txttag);

      //      c_stack->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/PFEnergyFractions_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+".pdf");
      c_stack->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/PFEnergyFractions_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]+".pdf");
      delete stack;
      delete c_stack;
      
      
      TCanvas* c_stack_relud =  new TCanvas("c_stack_relud","c_stack_relud",1);
      stack_relud->Draw("nostack"); //as 'same'
      leg2 -> Draw();
      tex->DrawLatex(0.15,0.27,Latextext1);
      tex->DrawLatex(0.15,0.23,Latextext2);
      tex->DrawLatex(0.15,0.19,Latextext3);
      tex->DrawLatex(0.15,0.95,txttag);
      c_stack_relud->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/PFEnergyFractions_RelToUD_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]+".pdf");
      delete stack_relud;
      delete c_stack_relud;
      
      
      //dummy for tdrCanvas
      TH1D *h = new TH1D("h","Response;Flavor (Physics);pt_{REC}/pt_{GEN}",7,0,7);
      for(int ipdg=0;ipdg<7;ipdg++){
	h->GetXaxis()->SetBinLabel(ipdg+1,pdgIDLabels[ipdg]);
      }
      
      
      TCanvas* c_response =  new TCanvas("c_response","c_response",1);
      tdrCanvas(c_response,"c_response",h,4,10,true,CorrectionObject::_lumitag);
      h->SetMarkerColor(2);
      h->SetLineColor(2);
      h->SetMarkerStyle(20);
      h->GetYaxis()->SetRangeUser(0.75,1.20);
      h->GetYaxis()->SetTitle("pt_{REC}/pt_{GEN}");
      h->GetXaxis()->SetTitle("Flavor");
      h_response_pt[k]->SetMarkerColor(2);
      h->GetYaxis()->SetRangeUser(0.75,1.20);
      h_response_pt[k]->SetLineColor(2);
      h_response_pt[k]->SetMarkerStyle(20);
      h_response_pt[k]->Draw("P SAME");
      
      TLegend* leg4 =  new TLegend(0.60,0.68,0.70,0.75);
      leg4 -> AddEntry(h_response, "response");
      leg4 -> SetFillColor(10);
      leg4 -> SetBorderSize(0);
      leg4 -> SetTextSize(0.042);
      leg4 -> SetLineColor(1);
      leg4 -> SetTextFont(42);
      leg4 -> Draw();
      tex->DrawLatex(0.22,0.27,Latextext1);
      tex->DrawLatex(0.22,0.23,Latextext2);
      tex->DrawLatex(0.22,0.19,Latextext3);
      tex->DrawLatex(0.15,0.95,txttag);

      c_response->SaveAs(CorrectionObject::_outpath+"/plots/control/Flavor/ResponsesFlavor_"+jettag+"_"+txttag+ "_eta_" + eta_range2[j] + "_" + eta_range2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]+".pdf");
      delete c_response;
    }
    //[END] in pt bins ----------------------
  }

  cout << "Hello, I'm done" << endl;

}
