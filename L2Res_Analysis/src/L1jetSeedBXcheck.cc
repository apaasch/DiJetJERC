#include "../include/parameters.h"
#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"
#include "../include/tdrstyle_mod15.h"

#include <TStyle.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TFile.h>
#include <TH2D.h>

using namespace std;

void CorrectionObject::L1jetSeedBXcheck(){
  cout << "--------------- Starting L1jetSeedBXcheck() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  bool eta_abs=true;
  int n_pt_ = max(n_pt,n_pt_HF);
  bool eta_cut_bool;
  int n_pt_cutted;
  
  // const int n_pt_ = (si_trg ? n_pt_Si : n_pt );

  // double pt_bins_[n_pt_]; = (si_trg ? pt_bins_Si : pt_bins);
  // TString pt_range_[n_pt_]; = (si_trg ? pt_range_Si : pt_range);

  // for(int i=0; i<n_pt_;i++)
  
  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/L1JetSeedBX/");

  
  //Table with number of events in each pT- and eta-bin
  TH1D *hdata_jet1[n_pt_-1][(eta_abs ? n_eta : n_eta_full)  -1];
  TH1D *hdata_jet2[n_pt_-1][(eta_abs ? n_eta : n_eta_full)  -1];
  TH1D *hdata_jet3[n_pt_-1][(eta_abs ? n_eta : n_eta_full)  -1];

  
  int count = 0;
  TString name1 = "hist_l1jetseedbx_jet1";
  TString name2 = "hist_l1jetseedbx_jet2";
  TString name3 = "hist_l1jetseedbx_jet3";

  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
      eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut;
      TString eta_name = "eta_"+(eta_abs ? eta_range2 : eta_range2_full)[j]+"_"+(eta_abs ? eta_range2 : eta_range2_full)[j+1];

    for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];
      

      TString name = name1 + eta_name + "_" + pt_name; 
      hdata_jet1[k][j] = new TH1D(name,"",21, -10, 11);
      name = name2 + eta_name + "_" + pt_name; 
      hdata_jet2[k][j] = new TH1D(name,"",21, -10, 11);
      name = name3 + eta_name + "_" + pt_name; 
      hdata_jet3[k][j] = new TH1D(name,"",21, -10, 11);

		 
      count++;
    }
  }


  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_data("AnalysisTree", CorrectionObject::_DATAFile);
  
  TTreeReaderValue<Float_t> bxjet1(myReader_data, "jet1_l1bx");
  TTreeReaderValue<Float_t> bxjet2(myReader_data, "jet1_l2bx");
  TTreeReaderValue<Float_t> bxjet3(myReader_data, "jet1_l3bx");
  TTreeReaderValue<Float_t> pt_ave_data(myReader_data, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_data, "probejet_eta");    

  int myCount = 0;
  int myCount_cut = 0;

  std::cout<<"\nstart reading data\n";
  
  myCount = 0;
  myCount_cut = 0;
  while (myReader_data.Next()) {

    if(myCount &&! myCount%10000)   std::cout<<myCount;
    
  
  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
	eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut;
    for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
	hdata_jet1[k][j]->Fill(*bxjet1);
	hdata_jet2[k][j]->Fill(*bxjet2);
	hdata_jet3[k][j]->Fill(*bxjet3);
    }
  }

  }
  

  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
    
  std::cout<<"\nSave\n";
  
  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
	eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut; 
  for(int k=0; k<n_pt-1; k++){
    TString eta_name = "eta_"+(eta_abs ? eta_range2 : eta_range2_full)[j]+"_"+(eta_abs ? eta_range2 : eta_range2_full)[j+1];
    TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];
    TString legname = "p_{T} #in [" +(eta_cut_bool?pt_range_HF:pt_range)[k] + "," +(eta_cut_bool?pt_range_HF:pt_range)[k+1]+ ",#eta #in [" +(eta_abs ? eta_range2 : eta_range2_full)[j] + "," +(eta_abs ? eta_range2 : eta_range2_full)[j+1]+"]+ jet2";
    TCanvas *c = new TCanvas();
    hdata_jet2[k][j]->GetXaxis()->SetTitle("L1 Jet Seed BX, +-10 are additional codes");
    hdata_jet2[k][j]->GetYaxis()->SetTitle("events per bin");
    hdata_jet2[k][j]->Draw();
    tex->DrawLatex(0.54,0.8,legname);
    c->SetLogy();
    c->SaveAs(CorrectionObject::_outpath+"plots/control/L1JetSeedBX/l1jetseedbx_jet1_"+pt_name+"_"+eta_name+".pdf","pdf");
  }
  }

   
  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
	eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut; 
  for(int k=0; k<n_pt-1; k++){
    TString eta_name = "eta_"+(eta_abs ? eta_range2 : eta_range2_full)[j]+"_"+(eta_abs ? eta_range2 : eta_range2_full)[j+1];
    TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];
    TString legname = "p_{T} #in [" +(eta_cut_bool?pt_range_HF:pt_range)[k] + "," +(eta_cut_bool?pt_range_HF:pt_range)[k+1]+ ",#eta #in [" +(eta_abs ? eta_range2 : eta_range2_full)[j] + "," +(eta_abs ? eta_range2 : eta_range2_full)[j+1]+"]+ jet2";
    TCanvas *c = new TCanvas();
    hdata_jet2[k][j]->GetXaxis()->SetTitle("L1 Jet Seed BX, +-10 are additional codes");
    hdata_jet2[k][j]->GetYaxis()->SetTitle("events per bin");
    hdata_jet2[k][j]->Draw();
    tex->DrawLatex(0.54,0.8,legname);
    c->SetLogy();
    c->SaveAs(CorrectionObject::_outpath+"plots/control/L1JetSeedBX/l1jetseedbx_jet2_"+pt_name+"_"+eta_name+".pdf","pdf");
  }
  }
 

  
  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
	eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut; 
  for(int k=0; k<n_pt-1; k++){
    TString eta_name = "eta_"+(eta_abs ? eta_range2 : eta_range2_full)[j]+"_"+(eta_abs ? eta_range2 : eta_range2_full)[j+1];
    TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];
    TString legname = "p_{T} #in [" +(eta_cut_bool?pt_range_HF:pt_range)[k] + "," +(eta_cut_bool?pt_range_HF:pt_range)[k+1]+ ",#eta #in [" +(eta_abs ? eta_range2 : eta_range2_full)[j] + "," +(eta_abs ? eta_range2 : eta_range2_full)[j+1]+"]+ jet3";
    TCanvas *c = new TCanvas();
    hdata_jet3[k][j]->GetXaxis()->SetTitle("L1 Jet Seed BX, +-10 are additional codes");
    hdata_jet3[k][j]->GetYaxis()->SetTitle("events per bin");
    hdata_jet3[k][j]->Draw();
    tex->DrawLatex(0.54,0.8,legname);
    c->SetLogy();
    c->SaveAs(CorrectionObject::_outpath+"plots/control/L1JetSeedBX/l1jetseedbx_jet3_"+pt_name+"_"+eta_name+".pdf","pdf");
  }
  }
  
}
    
    


