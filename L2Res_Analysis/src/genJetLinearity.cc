//#include "../include/parameters.h"
#include "../../include/constants.h"
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

void CorrectionObject::genJetLinearity(){
  cout << "--------------- StartingGenJetLinearity() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  bool eta_abs=true;
  int n_pt_ = max(n_pt,n_pt_HF);
  bool eta_cut_bool;
  int n_pt_cutted;
  
  // const int n_pt_ = (si_trg ? n_pt_Si : n_pt );

  // double pt_bins_[n_pt_]; = (si_trg ? pt_bins_Si : pt_bins);
  // TString pt_range_[n_pt_]; = (si_trg ? pt_range_Si : pt_range);

  // for(int i=0; i<n_pt_;i++)
  
  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/genJetlinearity/");

  int n_alpha = 50;
  
  //Table with number of events in each pT- and eta-bin
  TH1D *hmc_RelGenDiff[n_pt_-1][(eta_abs ? n_eta : n_eta_full)  -1][n_alpha];
  TH1D *hmc_alpha[n_pt_-1][(eta_abs ? n_eta : n_eta_full)  -1];

  TF1 *fitFuncs[n_pt_-1][(eta_abs ? n_eta : n_eta_full)  -1];

  TH1D *hmc_kfstrextrapolChi2Ndf[n_pt_-1];
  
  int count = 0;
  TString name1 = "hist_mc_RelGenDiff_";
  TString name2 = "hist_mc_alpha_";
  TString name3 = "hist_mc_ksfsrExtrapol_chi2PNDF";  

  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
      eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut;
      TString eta_name = "eta_"+(eta_abs ? eta_range2 : eta_range2_full)[j]+"_"+(eta_abs ? eta_range2 : eta_range2_full)[j+1];

    for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];

      if(j==0) hmc_kfstrextrapolChi2Ndf[k] = new TH1D();
      

      TString name = name2 + eta_name + "_" + pt_name; 
      hmc_alpha[k][j] = new TH1D(name,"",n_alpha, 0, 0.5);
      fitFuncs[k][j] = new TF1("","pol1",0.,0.5); 

      for(int i = 0; i<n_alpha; i++){
	name = name1 + eta_name + "_" + pt_name + "_" + Form("%d",i);
	hmc_RelGenDiff[k][j][i] = new TH1D(name,"",200, 0, 0.5);
      }
      
      name = name3  + pt_name;
	if(j==0) hmc_kfstrextrapolChi2Ndf[k] = new TH1D(name,"",n_eta, 0, 5.);
		 
      count++;
    }
  }


  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  
  TTreeReaderValue<Float_t> genjet1_pt(myReader_MC, "jet1_ptGen");
  TTreeReaderValue<Float_t> genjet2_pt(myReader_MC, "jet2_ptGen");
  TTreeReaderValue<Float_t> genjet3_pt(myReader_MC, "jet3_ptGen");
  TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");    

  int myCount = 0;
  int myCount_cut = 0;

  std::cout<<"\nstart reading MC\n";
  
  myCount = 0;
  myCount_cut = 0;
  while (myReader_MC.Next()) {

    if(myCount &&! myCount%10000)   std::cout<<myCount;
    
      for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
	eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut;
    for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      for(int i=0; i<n_alpha; i++){
	double pt_ave = (*genjet1_pt+*genjet2_pt)/2;
	double alpha = *genjet3_pt/pt_ave;
	if( alpha > i*0.01 ) continue;
      if(*pt_ave_mc<(eta_cut_bool?pt_bins_HF:pt_bins)[k] || *pt_ave_mc>(eta_cut_bool?pt_bins_HF:pt_bins)[k+1]) continue;
	if( (eta_abs ? fabs(*probejet_eta_mc) : *probejet_eta_mc) > (eta_abs ? eta_bins : eta_bins_full)[j+1] || (eta_abs ? fabs(*probejet_eta_mc) : *probejet_eta_mc)  < (eta_abs ? eta_bins : eta_bins_full)[j]) continue;
	else{
	  double genDiff = *genjet1_pt - *genjet2_pt;
	  hmc_RelGenDiff[k][j][i]->Fill(genDiff/pt_ave, *weight_mc);

	}
      }
    }
    }
	  myCount++;
  }
  
  std::cout<<"\nGet all the mean stuff\n";
  
  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
	eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut;
    for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      for(int i=0; i<n_alpha; i++){
	hmc_alpha[k][j]->SetBinContent(i+1,hmc_RelGenDiff[k][j][i]->GetMean());
	hmc_alpha[k][j]->SetBinError(i+1,hmc_RelGenDiff[k][j][i]->GetMeanError());
      }
      
    }
  }
  
  std::cout<<"\nFill the Chi2/NDF\n";
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
    
  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
	eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut;
    for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      if(hmc_alpha[k][j]){
	hmc_alpha[k][j]->Fit(fitFuncs[k][j],"","",0.14,0.36);
	hmc_kfstrextrapolChi2Ndf[k]->SetBinContent(j,fitFuncs[k][j]->GetChisquare()/fitFuncs[k][j]->GetNDF());
      }
      else hmc_kfstrextrapolChi2Ndf[k]->SetBinContent(j,0.);
    }
  }
  
  std::cout<<"\nSave the Chi2/NDF\n";
  
  for(int k=0; k<n_pt-1; k++){
    TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];
    TString legname = "p_{T} #in [" +pt_range[k] + "," +pt_range[k+1]+ "]";
    TCanvas *c = new TCanvas();
    hmc_kfstrextrapolChi2Ndf[k]->GetXaxis()->SetTitle("eta");
    hmc_kfstrextrapolChi2Ndf[k]->GetYaxis()->SetTitle("chi2 per ndf");
    hmc_kfstrextrapolChi2Ndf[k]->Draw();
    tex->DrawLatex(0.54,0.8,legname);	
    c->SaveAs(CorrectionObject::_outpath+"plots/control/genJetlinearity/Chi2PerNDFOfKFSRextrapol_"+pt_name+".pdf","pdf");
  }

  
    std::cout<<"\nLook for the Chi2/NDF<2\n";

    double up, down, chi2NDF;
    ofstream myfile;
  myfile.open (CorrectionObject::_outpath+"plots/control/genJetlinearity/fitrangeChi2NDFis2.txt",ios::trunc);
  myfile<<"\npt\t eta\t down\t up\n\n";
  for(int j=0; j<(eta_abs ? n_eta : n_eta_full)-1; j++){
	eta_cut_bool = fabs((eta_abs ? eta_bins : eta_bins_full)[j])>=eta_cut;
    for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      myfile<<pt_range[k];
      myfile<<"\t";
      myfile<<eta_range[j];
      myfile<<"\t";
      down = 0.14;
      up = 0.36;
      chi2NDF=fitFuncs[k][j]->GetChisquare()/fitFuncs[k][j]->GetNDF();
      if(chi2NDF<2) continue;
      while(chi2NDF<2){
	down-=0.05;
	up+=0.05;
	hmc_alpha[k][j]->Fit(fitFuncs[k][j],"","",down,up);
	chi2NDF=fitFuncs[k][j]->GetChisquare()/fitFuncs[k][j]->GetNDF();
      }
      myfile<<Form("%.2f\t%.2f\n",down,up);
      
    }
  }


 myfile.close();

}
    
    


