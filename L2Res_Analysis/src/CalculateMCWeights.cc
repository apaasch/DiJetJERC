#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
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
#include <TProfile2D.h>



using namespace std;

void CorrectionObject::CalculateMCWeights(){

  //************************************************************************
  //
  // Calculate MC Weights to scale MC to data
  // Re-weighting is done in merged pT thresholds
  // and eta bins to take trigger differencies into account
  // Scale factors defined as data/mc
  // Save SF in 2D histogram
  //
  // Loop over events for second time to apply right weights
  // Fill pT_ave hists to Cross-check re-weighting procedure 
  //
  //************************************************************************


  cout << "--------------- Starting CalculateMCWeights() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  const int n_pt_bins = 100;


  //Check if Files are available
  cout<<"MC Path for Weights: "<<_MCpath_ForWeights<<endl;
   cout<<"Data Path for Weights: "<<_DATApath_ForWeights<<endl;

  TFile* f_mc   = new TFile(CorrectionObject::_MCpath_ForWeights,"READ");
  TFile* f_data = new TFile(CorrectionObject::_DATApath_ForWeights,"READ");
     if(f_mc->GetSize()==-1) throw runtime_error("In CorrectionObject.cc: File or Directory " + _MCpath_ForWeights+" does not exist!");
     if(f_data->GetSize()==-1) throw runtime_error("In CorrectionObject.cc: File or Directory " + _DATApath_ForWeights+" does not exist!");

  //pT and eta bins  
  // double bins[16] = {51, 73, 95, 100, 126,152, 163, 230, 250, 299, 316, 365, 433, 453, 566, 1000}; //take pt-bins from paramters.h
  double bins_eta[5] = {-5.2, -2.853, 0, 2.853, 5.2};

 //pT_ave-histograms for MC & DATA
  TH2D* h_pt_ave_mc = new TH2D("pt_ave_mc","pt_ave mc;p_{T}^{ave};|#eta|", n_pt_bins,0,5000,n_eta-1, eta_bins);
  TH2D* h_pt_ave_data = new TH2D("pt_ave_data","pt_ave data;p_{T}^{ave};|#eta|", n_pt_bins,0,5000,n_eta-1,eta_bins);
  TH2D* h_pt_ave_mc_scaled = new TH2D("pt_ave_mc_scaled","pt_ave mc scaled;p_{T}^{ave};|#eta|", n_pt_bins,0,5000,n_eta-1,eta_bins);
  TH1D* h1_pt_ave_mc_scaled = new TH1D("h1_pt_ave_mc_scaled", "pt_ave mc scaled;p_{T}^{ave};entries", n_pt_bins, 0, 5000); //cross-check 1d
  TH1D* h1_pt_ave_data = new TH1D("h1_pt_ave_data", "pt_ave data;p_{T}^{ave};entries", n_pt_bins, 0, 5000);                //cross-check 1d


  TH2Poly* h_pt_ave_binned_mc = new TH2Poly();
  h_pt_ave_binned_mc->SetName("pt_ave_binned_mc");  
  h_pt_ave_binned_mc->SetTitle("pt_ave binned mc");
  h_pt_ave_binned_mc->GetYaxis()->SetTitle("|#eta|");
  h_pt_ave_binned_mc->GetXaxis()->SetTitle("p_{T}^{ave}");
  bool eta_cut_bool;
  for(int i=0; i<5-1; i++){
    eta_cut_bool = fabs(bins_eta[i])>eta_cut && fabs(bins_eta[i+1])>eta_cut;
       for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
	 h_pt_ave_binned_mc->AddBin( (eta_cut_bool?pt_bins_HF:pt_bins)[j], bins_eta[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1], bins_eta[i+1] );
       }
  }
  
 
  TH2Poly* h_pt_ave_binned_data = new TH2Poly();
  h_pt_ave_binned_data->SetName("pt_ave_binned_data");  
  h_pt_ave_binned_data->SetTitle("pt_ave binned data");
  h_pt_ave_binned_data->GetYaxis()->SetTitle("#eta");
  h_pt_ave_binned_data->GetXaxis()->SetTitle("p_{T}"); 
  for(int i=0; i<5-1; i++){
    eta_cut_bool = fabs(bins_eta[i])>eta_cut && fabs(bins_eta[i+1])>eta_cut;
       for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
	 h_pt_ave_binned_data->AddBin( (eta_cut_bool?pt_bins_HF:pt_bins)[j], bins_eta[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1], bins_eta[i+1] );
       }
  }

  TH2Poly* h_pt_ave_binned_yield = new TH2Poly();
  h_pt_ave_binned_yield->SetName("pt_ave_binned_yield");  
  h_pt_ave_binned_yield->SetTitle("pt_ave binned yield");
  h_pt_ave_binned_yield->GetYaxis()->SetTitle("#eta");
  h_pt_ave_binned_yield->GetXaxis()->SetTitle("p_{T}"); 
  for(int i=0; i<5-1; i++){
    eta_cut_bool = fabs(bins_eta[i])>eta_cut && fabs(bins_eta[i+1])>eta_cut;
       for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
	 h_pt_ave_binned_yield->AddBin( (eta_cut_bool?pt_bins_HF:pt_bins)[j], bins_eta[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1], bins_eta[i+1] );
       }
  }  
  
  // TH2D* h_pt_ave_binned_yield = new TH2D("pt_ave_binned_yield","CENTRAL Triggers;p_{T}^{ave};|#eta|", 15, bins, n_eta-1, eta_bins);

  TH2Poly* h_pt_ave_binned_mc_scaled = new TH2Poly();
  h_pt_ave_binned_mc_scaled->SetName("pt_ave_binned_mc_scaled");  
  h_pt_ave_binned_mc_scaled->SetTitle("pt_ave binned mc scaled");
  h_pt_ave_binned_mc_scaled->GetYaxis()->SetTitle("|#eta|");
  h_pt_ave_binned_mc_scaled->GetXaxis()->SetTitle("p_{T}^{ave}"); 
  for(int i=0; i<5-1; i++){
    eta_cut_bool = fabs(bins_eta[i])>eta_cut && fabs(bins_eta[i+1])>eta_cut;
       for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
	 h_pt_ave_binned_mc_scaled->AddBin( (eta_cut_bool?pt_bins_HF:pt_bins)[j], bins_eta[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1], bins_eta[i+1] );
       }
  }  
  
  
  // TH1D* h1_pt_ave_binned_mc_scaled = new TH1D("h1_pt_ave_binned_mc_scaled", "pt_ave mc scaled;p_{T}^{ave};entries", 15, bins);   //cross-check 1d
  // TH1D* h1_pt_ave_binned_data = new TH1D("h1_pt_ave_binned_data", "pt_ave data;p_{T}^{ave};entries", 15, bins);                  //cross-check 1d 


  //Fill histograms
  //DATA:
  TTreeReader myReader_data("AnalysisTree", f_data);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_data, "pt_ave");
  TTreeReaderValue<Float_t> weight_data(myReader_data, "weight");
  TTreeReaderValue<Float_t> alpha_data(myReader_data, "alpha");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_data, "probejet_eta");
  while (myReader_data.Next()){
    if(*alpha_data >= 0.3 || *pt_ave_data < pt_bins[0]) continue;
    h_pt_ave_data->Fill(*pt_ave_data, *probejet_eta_data, *weight_data);
    h_pt_ave_binned_data->Fill(*pt_ave_data, *probejet_eta_data, *weight_data);
    h_pt_ave_binned_yield->Fill(*pt_ave_data, *probejet_eta_data, *weight_data);
    h1_pt_ave_data->Fill(*pt_ave_data, *weight_data);
    // h1_pt_ave_binned_data->Fill(*pt_ave_data, *weight_data);
  }

  h_pt_ave_binned_yield->SetMinimum(0);
  h_pt_ave_binned_yield->SetMaximum(1000000);

  //MC:
  TTreeReader myReader_mc("AnalysisTree", f_mc);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_mc, "pt_ave");
  TTreeReaderValue<Float_t> weight_mc(myReader_mc, "weight");
  TTreeReaderValue<Float_t> alpha_mc(myReader_mc, "alpha");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_mc, "probejet_eta");
  while (myReader_mc.Next()){
    if(*alpha_mc >= 0.3  || *pt_ave_mc<pt_bins[0]) continue;
    h_pt_ave_mc->Fill(*pt_ave_mc, *probejet_eta_mc, *weight_mc);
    h_pt_ave_binned_mc->Fill(*pt_ave_mc, *probejet_eta_mc, *weight_mc);
  }
  cout<<"\nDATA fine binned entries "<<h_pt_ave_data->GetEntries()<<endl;
  cout<<"DATA fine binned effective entries "<<h_pt_ave_data->GetEffectiveEntries()<<endl;
  cout<<"DATA fine binned integral "<<h_pt_ave_data->Integral()<<endl;
  cout<<"DATA pt binned entries "<<h_pt_ave_binned_data->GetEntries()<<endl;
  cout<<"DATA pt binned effective entries "<<h_pt_ave_binned_data->GetEffectiveEntries()<<endl;
  cout<<"DATA pt binned integral "<<h_pt_ave_binned_data->Integral()<<endl;
  cout<<"DATA pt binned mean "<<h_pt_ave_binned_data->GetMean()<<endl;  
  cout<<"DATA pt binned rms "<<h_pt_ave_binned_data->GetStdDev()<<endl;
  cout<<"\nMC fine binned entries "<<h_pt_ave_mc->GetEntries()<<endl;
  cout<<"MC fine binned effective entries "<<h_pt_ave_mc->GetEffectiveEntries()<<endl;
  cout<<"MC fine binned integral "<<h_pt_ave_mc->Integral()<<endl;
  cout<<"MC pt binned entries "<<h_pt_ave_binned_mc->GetEntries()<<endl;
  cout<<"MC pt binned effective entries "<<h_pt_ave_binned_mc->GetEffectiveEntries()<<endl;
  cout<<"MC pt binned integral "<<h_pt_ave_binned_mc->Integral()<<endl;
  cout<<"MC pt binned mean "<<h_pt_ave_binned_mc->GetMean()<<endl;  
  cout<<"MC pt binned rms "<<h_pt_ave_binned_mc->GetStdDev()<<endl;
   
 //Calculate scale factors
  TH2Poly* SF =  (TH2Poly*)h_pt_ave_binned_data->Clone();
  SF->SetMinimum(0);
  SF->SetMaximum(1);
  TH2Poly* for_SF_mc =  (TH2Poly*)h_pt_ave_binned_mc->Clone();
  
  int idx_runing=0;
  for(int j=0; j<4; j++){
    eta_cut_bool = fabs(bins_eta[j])>eta_cut && fabs(bins_eta[j+1])>eta_cut;
    for(int i= 0 ; i <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; i++ ){
      double content = 0;
      double data = SF->GetBinContent(idx_runing);
      double mc = for_SF_mc->GetBinContent(idx_runing);
      if(mc > 0) {content = data/mc; // cout<<"data/mc "<<data<<"/"<<mc<<endl;
      }
      else content = 0;
      SF->SetBinContent(idx_runing, content);
      //cout << "lower pt: " << 5+i*5 << ", eta-bin: " << j+1 << " DATA: " << h_pt_ave_data->GetBinContent(i+1,j+1) << ", MC: " << for_SF_mc->GetBinContent(i+1,j+1) << ", content: " << content << endl;
      idx_runing++;
    }
  }

  //do a dummy crosscheck first...
  for(int k=0; k<idx_runing; k++){
    if(fabs(h_pt_ave_binned_data->GetBinContent(k)-(h_pt_ave_binned_mc->GetBinContent(k)*SF->GetBinContent(k)))>0.5){
      cout<<"Bin "<<k<<" data: "<<h_pt_ave_binned_data->GetBinContent(k)<<" unscaled MC: "<<h_pt_ave_binned_mc->GetBinContent(k)<<" SF: "<<SF->GetBinContent(k)<<" scaled MC: "<<h_pt_ave_binned_mc->GetBinContent(k)*SF->GetBinContent(k)<<endl<<endl;
      if(int(h_pt_ave_binned_mc->GetBinContent(k))){
	throw runtime_error("did not pass dummy crosscheck");
      }
    }
  }
  cout<<"after dummy crosscheck"<<endl;

//looping over events for the 2nd time
//Create re-weighted pT_ave hists to cross-check procedure

  TTreeReader myReader_mc2("AnalysisTree", f_mc);
  TTreeReaderValue<Float_t> pt_ave_mc2(myReader_mc2, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc2(myReader_mc2, "probejet_eta");
  TTreeReaderValue<Float_t> weight_mc2(myReader_mc2, "weight");
  TTreeReaderValue<Float_t> alpha_mc2(myReader_mc2, "alpha");

  int ind = 0;
  while (myReader_mc2.Next()){
    if(*alpha_mc2 >= 0.3) continue;
    double right_SF = 0;
    double pt = *pt_ave_mc2;
    double eta = *probejet_eta_mc2;

    int idx_x = 0;
    int idx_y = 0;
    eta_cut_bool = fabs(eta)>eta_cut;
    idx_runing = 0;
    while(eta > bins_eta[idx_y]){
      idx_y++;
      idx_runing+= ( eta_cut_bool ?  n_pt_HF : n_pt );
    }
    while(pt > (eta_cut_bool?pt_bins_HF:pt_bins)[idx_x]){ 
      idx_x++;
      idx_runing++;
    }

    
    //Now i is the number of the bin, whose weight has to be taken for this event
    // cout<<"glob findBin: "<<SF->FindBin(pt,eta)<<" idx runing: "<<idx_runing<<endl;
    if(SF->FindBin(pt,eta)<0)
      if(SF->FindBin(2000-1e-7,eta)<0)
	cout<<"didnt found SF for pt,eta "<<pt<<","<<eta<<endl;
    right_SF = SF->GetBinContent(SF->FindBin((pt>2000?(2000-1e-7):pt),eta));

    double right_weight = *weight_mc2 * right_SF;

    h_pt_ave_mc_scaled->Fill(pt, eta, right_weight);
    h_pt_ave_binned_mc_scaled->Fill(pt, eta, right_weight);
    h1_pt_ave_mc_scaled->Fill(pt, right_weight);
    // h1_pt_ave_binned_mc_scaled->Fill(pt, right_weight);

    ind++;
  }

   cout<<"\nMC pt scaled entries "<<h_pt_ave_mc_scaled->GetEntries()<<endl;
   cout<<"MC pt scaled effective entries "<<h_pt_ave_mc_scaled->GetEffectiveEntries()<<endl;
   cout<<"MC pt scaled integral "<<h_pt_ave_mc_scaled->Integral()<<endl;
   cout<<"MC pt scaled mean "<<h_pt_ave_mc_scaled->GetMean()<<endl;  
   cout<<"MC pt scaled rms "<<h_pt_ave_mc_scaled->GetStdDev()<<endl;


   double count_xcheck_mc =0.;
   
    idx_runing = 0;
    for(int j=0; j<n_eta-1; j++){
      eta_cut_bool = fabs(bins_eta[j])>eta_cut && fabs(bins_eta[j+1])>eta_cut;
      for(int i= 0 ; i < 1000 ; i++ ){
	if(fabs(h_pt_ave_mc_scaled->GetBinContent(i,j) - h_pt_ave_data->GetBinContent(i,j))/h_pt_ave_mc_scaled->GetBinContent(i,j) > 0.05 && h_pt_ave_data->GetBinContent(i,j)>0){
	  cout << "\nData: " << h_pt_ave_data->GetBinContent(i,j) << ", MC: " << h_pt_ave_mc_scaled->GetBinContent(i,j) << endl;
	count_xcheck_mc+=h_pt_ave_mc_scaled->GetBinContent(i,j);
	cout << "Data: " << h_pt_ave_data->GetBinContent(i,j) << ", SF: " << SF->GetBinContent(SF->FindBin(i*5000./n_pt_bins,eta_bins[j])) << endl;
	  cout << "Eta_low: " << eta_bins[j] << ", pT: " << i*5000/n_pt_bins << ", MC: " << h_pt_ave_mc_scaled->GetBinContent(i,j) << ", DATA: " << h_pt_ave_data->GetBinContent(i,j) << endl;
	  cout<<"so MC is off by "<<fabs(h_pt_ave_mc_scaled->GetBinContent(i,j) - h_pt_ave_data->GetBinContent(i,j))/h_pt_ave_mc_scaled->GetBinContent(i,j)*100<<"%"<<endl;
	if(SF->GetBinContent(SF->FindBin(i*5000./n_pt_bins,eta_bins[j]))>0.){
	  cout<<"MC eta-pt binned count "<<for_SF_mc->GetBinContent(for_SF_mc->FindBin(i*5000./n_pt_bins,eta_bins[j]))<<endl;
	  // throw runtime_error("In CalculateMCWeights.cc: Scaled MC bin-content does not match the bin-content in data.");
	}
      }
      idx_runing++;
    }
  }

    cout<<"counted scaled fine binned MC to be "<< count_xcheck_mc<<endl;

  //Save re-weighted pT_ave hists
  //Setup Canvas Fine binnig 
  TCanvas* c1 = new TCanvas();
  h1_pt_ave_mc_scaled->SetLineWidth(2);
  h1_pt_ave_mc_scaled->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
  h1_pt_ave_mc_scaled->GetXaxis()->SetRangeUser(20,1000);
  h1_pt_ave_mc_scaled->Draw("HIST");
  h1_pt_ave_data->SetLineColor(2);
  h1_pt_ave_data->SetMarkerColor(2);
  h1_pt_ave_data->SetMarkerSize(0.5);
  h1_pt_ave_data->SetMarkerStyle(20);
  h1_pt_ave_data->Draw("P SAME");

  TLegend* l1 = new TLegend(0.52,0.7,0.9,0.9);    
  l1->AddEntry(h1_pt_ave_mc_scaled,"MC after reweighting","l");
  l1->AddEntry(h1_pt_ave_data,"DATA","p");
  l1->Draw();

  // //Save
  c1->SaveAs(CorrectionObject::_weight_path + "MC_scaled_PtEta_Fine.pdf");

  // //Setup Canvas binned 
  // TCanvas* c2 = new TCanvas();
  // h1_pt_ave_binned_mc_scaled->SetLineWidth(2);
  // h1_pt_ave_binned_mc_scaled->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
  // h1_pt_ave_binned_mc_scaled->GetXaxis()->SetRangeUser(20,1000);
  // h1_pt_ave_binned_mc_scaled->Draw("HIST");
  // h1_pt_ave_binned_data->SetLineColor(2);
  // h1_pt_ave_binned_data->SetMarkerColor(2);
  // h1_pt_ave_binned_data->SetMarkerStyle(20);
  // h1_pt_ave_binned_data->Draw("P SAME");

  // TLegend* l2 = new TLegend(0.52,0.2,0.90,0.4);    
  // l2->AddEntry(h1_pt_ave_binned_mc_scaled,"MC after reweighting","l");
  // l2->AddEntry(h1_pt_ave_binned_data,"DATA","p");
  // l2->Draw();

  // //Save
  // c2->SaveAs(CorrectionObject::_weight_path + "MC_scaled_PtEta_Fine_binned.pdf");

  if(!  CorrectionObject::make_path(CorrectionObject::_weight_path)){
    cout << "No new Directory was created" << endl;
  }

  //Write output
  TFile* out = new TFile(CorrectionObject::_weight_path + "/MC_ReWeights_Run" + CorrectionObject::_runnr  + CorrectionObject::_outpath_postfix +".root","RECREATE");
  SF->Write();
  h_pt_ave_binned_yield->Write();
  out->Close();

  //delete out;
  // delete l2;
  // delete c2;
  delete l1;
  delete c1;
  delete for_SF_mc;
  delete SF;
  delete h_pt_ave_binned_mc_scaled;
  delete h_pt_ave_binned_data;
  delete h_pt_ave_binned_mc;
  delete h_pt_ave_mc_scaled;
  delete h_pt_ave_data;
  delete h_pt_ave_mc;
  delete f_data;
  delete f_mc;

}
