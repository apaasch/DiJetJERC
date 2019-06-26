#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TProfile.h>



using namespace std;

void CorrectionObject::CalculateMCWeights_TriggerThresholds(bool CentralTriggers){
  cout << "--------------- Starting CalculateMCWeights() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  //Files
  TFile* f_mc, *f_data;
  if(CentralTriggers){
    f_mc   = new TFile(CorrectionObject::_MCpath_ForWeights_FLAT,"READ");
    f_data = new TFile(CorrectionObject::_DATApath_ForWeights_FLAT,"READ");
  }
  else{
    f_mc   = new TFile(CorrectionObject::_MCpath_ForWeights_FWD,"READ");
    f_data = new TFile(CorrectionObject::_DATApath_ForWeights_FWD,"READ");

  }

   int n_pt_bins;

  if(CentralTriggers){
    n_pt_bins = 10;
  }
  else{
    n_pt_bins = 7;
  }

  double bins[n_pt_bins];

  if(CentralTriggers){
    bins[0] = 51;
    bins[1] = 73;
    bins[2] = 95;
    bins[3] = 163;
    bins[4] = 230;
    bins[5] = 299;
    bins[6] = 365;
    bins[7] = 453;
    bins[8] = 566;
    bins[9] = 1000;
  }
  else{
    bins[0] = 100;
    bins[1] = 126;
    bins[2] = 152;
    bins[3] = 250;
    bins[4] = 316;
    bins[5] = 433;
    bins[6] = 1000;
  }
  TH1D* h_pt_ave_binned_mc = new TH1D("pt_ave_binned_mc","pt_ave binned mc;p_{T}^{ave}", n_pt_bins-1, bins);
  TH1D* h_pt_ave_binned_data = new TH1D("pt_ave_binned_data","pt_ave binned data;p_{T}^{ave}", n_pt_bins-1, bins);
  TH1D* h_pt_ave_binned_mc_scaled = new TH1D("pt_ave_binned_mc_scaled","pt_ave binned mc scaled;p_{T}^{ave}", n_pt_bins-1, bins);
  TH1D* h_pt_ave_mc_scaled = new TH1D("h1_pt_ave_mc_scaled", "pt_ave mc scaled;p_{T}^{ave};entries", 1000, 0, 5000); //cross-check 1d
  TH1D* h_pt_ave_data = new TH1D("h1_pt_ave_data", "pt_ave data;p_{T}^{ave};entries", 1000, 0, 5000);                //cross-check 1d
  


  //  std::cout<<"bins[0] = "<<bins[0]<<std::endl;
  //Fill histograms
  TTreeReader myReader_data("AnalysisTree", f_data);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_data, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_data, "probejet_eta");
  TTreeReaderValue<Float_t> weight_data(myReader_data, "weight");
  TTreeReaderValue<Float_t> alpha_data(myReader_data, "alpha");
  while (myReader_data.Next()){
    if(*alpha_data >= 0.3 || *pt_ave_data < bins[0]) continue;//|| *probejet_eta_data < eta_min || *probejet_eta_data > eta_max
    h_pt_ave_binned_data->Fill(*pt_ave_data, *weight_data);
    h_pt_ave_data->Fill(*pt_ave_data, *weight_data);
  }

  TTreeReader myReader_mc("AnalysisTree", f_mc);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_mc, "pt_ave");
  TTreeReaderValue<Float_t> weight_mc(myReader_mc, "weight");
  TTreeReaderValue<Float_t> alpha_mc(myReader_mc, "alpha");
  while (myReader_mc.Next()){
    if(*alpha_mc >= 0.3 || *pt_ave_mc<bins[0]) continue;
    h_pt_ave_binned_mc->Fill(*pt_ave_mc, *weight_mc);
  }

  TH1D* SF =  (TH1D*)h_pt_ave_binned_data->Clone();
  TH1D* for_SF_mc =  (TH1D*)h_pt_ave_binned_mc->Clone();
  
  
  for(int i=0; i<SF->GetNbinsX(); i++){
    double content = 0;
    double data = SF->GetBinContent(i+1);
    double mc = for_SF_mc->GetBinContent(i+1);
    if(mc > 0) content = data/mc;
    else content = 0;
    //    std::cout<<"bin #"<<i<<" SF = "<<content<<std::endl;
    SF->SetBinContent(i+1, content);
  }
  


  TTreeReader myReader_mc2("AnalysisTree", f_mc);
  TTreeReaderValue<Float_t> pt_ave_mc2(myReader_mc2, "pt_ave");
  TTreeReaderValue<Float_t> weight_mc2(myReader_mc2, "weight");
  TTreeReaderValue<Float_t> alpha_mc2(myReader_mc2, "alpha");

  int ind = 0;
  while (myReader_mc2.Next()){
    if(*alpha_mc2 >= 0.3) continue;
    double right_SF = 0;
    double pt = *pt_ave_mc2; 

    int idx = 0;
    while(pt > bins[idx]){ 
      idx++;
    }
    
    //Now i is the number of the bin, whose weight has to be taken for this event
    right_SF = SF->GetBinContent(idx);

    double right_weight = *weight_mc2 * right_SF;

    h_pt_ave_mc_scaled->Fill(pt, right_weight);
    h_pt_ave_binned_mc_scaled->Fill(pt, right_weight);

    ind++;
  } //looping over events for the 2nd time

  for(int i=0; i<n_pt_bins; i++){
    if(fabs(h_pt_ave_binned_mc_scaled->GetBinContent(i+1) - h_pt_ave_binned_data->GetBinContent(i+1)) > 0.05) cout << "pT_low: " << bins[i] << ", MC: " << h_pt_ave_binned_mc_scaled->GetBinContent(i+1) << ", DATA: " << h_pt_ave_binned_data->GetBinContent(i+1) << endl;
  }

  TCanvas* c1 = new TCanvas();
  h_pt_ave_mc_scaled->SetLineWidth(2);
  h_pt_ave_mc_scaled->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
  h_pt_ave_mc_scaled->GetXaxis()->SetRangeUser(20,1000);
  h_pt_ave_mc_scaled->Draw("HIST");
  h_pt_ave_data->SetLineColor(2);
  h_pt_ave_data->SetMarkerColor(2);
  h_pt_ave_data->SetMarkerSize(0.5);
  h_pt_ave_data->SetMarkerStyle(20);
  h_pt_ave_data->Draw("P SAME");

  TLegend* l1 = new TLegend(0.52,0.7,0.9,0.9);    
  l1->AddEntry(h_pt_ave_mc_scaled,"MC after reweighting","l");
  l1->AddEntry(h_pt_ave_data,"DATA","p");
  l1->Draw();

  if(CentralTriggers) c1->SaveAs(CorrectionObject::_weightpath_FLAT + "MC_scaled_Run"+CorrectionObject::_runnr+"_Central_TriggerThresholds_NoUnflat_Fine.pdf");
  else                c1->SaveAs(CorrectionObject::_weightpath_FWD + "MC_scaled_Run"+CorrectionObject::_runnr+"_Fwd_TriggerThresholds_NoUnflat_Fine.pdf");


  TCanvas* c2 = new TCanvas();
  h_pt_ave_binned_mc_scaled->SetLineWidth(2);
  h_pt_ave_binned_mc_scaled->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
  h_pt_ave_binned_mc_scaled->GetXaxis()->SetRangeUser(20,1000);
  h_pt_ave_binned_mc_scaled->Draw("HIST");
  h_pt_ave_binned_data->SetLineColor(2);
  h_pt_ave_binned_data->SetMarkerColor(2);
  h_pt_ave_binned_data->SetMarkerStyle(20);
  h_pt_ave_binned_data->Draw("P SAME");

  TLegend* l2 = new TLegend(0.52,0.2,0.90,0.4);    
  l2->AddEntry(h_pt_ave_binned_mc_scaled,"MC after reweighting","l");
  l2->AddEntry(h_pt_ave_binned_data,"DATA","p");
  l2->Draw();

  if(CentralTriggers) c2->SaveAs(CorrectionObject::_weightpath_FLAT + "MC_scaled_Run"+CorrectionObject::_runnr+"_Central_TriggerThresholds_NoUnflat_binned.pdf");
  else                c2->SaveAs(CorrectionObject::_weightpath_FWD + "MC_scaled_Run"+CorrectionObject::_runnr+"_Fwd_TriggerThresholds_NoUnflat_binned.pdf");


  TFile* out;
  if(CentralTriggers) out = new TFile(CorrectionObject::_weightpath_FLAT + "/MC_ReWeights_CENTRAL_Run" + CorrectionObject::_runnr  + ".root","RECREATE");
  else                out = new TFile(CorrectionObject::_weightpath_FWD + "/MC_ReWeights_FWD_Run" + CorrectionObject::_runnr  + ".root","RECREATE");
  SF->Write();
  out->Close();

  //delete out;
  delete l2;
  delete c2;
  delete l1;
  delete c1;
  delete for_SF_mc;
  delete SF;
  delete h_pt_ave_binned_mc_scaled;
  delete h_pt_ave_binned_data;
  delete h_pt_ave_binned_mc;
  delete h_pt_ave_mc_scaled;
  delete h_pt_ave_data;
  delete f_data;
  delete f_mc;

}
