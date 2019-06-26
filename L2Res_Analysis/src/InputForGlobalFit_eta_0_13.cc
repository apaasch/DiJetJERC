// This script stores ratio between MC and DATA responses as function of pT in different eta bins and for different alpha values

#include "../include/parameters.h"
#include "../include/CorrectionObject.h"
#include "../include/useful_functions.h"

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


using namespace std;


void CorrectionObject::InputForGlobalFit_eta_0_13(){
  const int min_number_events = 20;//required number of events in response histogram
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(0);

  cout << "hello, it's mikkos macro" << endl;
  //calculate ratio in MC to DATA responses

  double mc_al_mpf[n_alpha_common][n_eta_common_2-1][n_pt-1]; 
  double err_mc_al_mpf[n_alpha_common][n_eta_common_2-1][n_pt-1]; 
  double mc_al_pTbal[n_alpha_common][n_eta_common_2-1][n_pt-1]; 
  double err_mc_al_pTbal[n_alpha_common][n_eta_common_2-1][n_pt-1]; 
  double data_al_mpf[n_alpha_common][n_eta_common_2-1][n_pt-1]; 
  double err_data_al_mpf[n_alpha_common][n_eta_common_2-1][n_pt-1]; 
  double data_al_pTbal[n_alpha_common][n_eta_common_2-1][n_pt-1]; 
  double err_data_al_pTbal[n_alpha_common][n_eta_common_2-1][n_pt-1]; 

  double ratio_al_mpf[n_alpha_common][n_eta_common_2-1][n_pt-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf[n_alpha_common][n_eta_common_2-1][n_pt-1]; //error of ratio at pt,eta,alpha bins
  double ratio_al_pTbal[n_alpha_common][n_eta_common_2-1][n_pt-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_pTbal[n_alpha_common][n_eta_common_2-1][n_pt-1]; //error of ratio at pt,eta,alpha bins

  TH1D *hdata_rel_r[n_pt-1][n_eta_common_2-1][n_alpha_common];// pT-balance responce for data
  TH1D *hdata_mpf_r[n_pt-1][n_eta_common_2-1][n_alpha_common];//MPF responce for data
  TH1D *hmc_rel_r[n_pt-1][n_eta_common_2-1][n_alpha_common];// pT-balance responce for MC
  TH1D *hmc_mpf_r[n_pt-1][n_eta_common_2-1][n_alpha_common];//MPF responce for MC
  int count = 0;
  TString name1 = "hist_data_rel_r_";
  TString name2 = "hist_data_mpf_r_";
  TString name3 = "hist_mc_rel_r_";
  TString name4 = "hist_mc_mpf_r_";

  //Initialize with 0 values
  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common_2-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_al_mpf[i][j][k] = 0;
	err_ratio_al_mpf[i][j][k] = 0;
	ratio_al_pTbal[i][j][k] = 0;
	err_ratio_al_pTbal[i][j][k] = 0;
	mc_al_mpf[i][j][k] = 0;
	err_mc_al_mpf[i][j][k] = 0; 
	mc_al_pTbal[i][j][k] = 0;  
	err_mc_al_pTbal[i][j][k] = 0;
	data_al_mpf[i][j][k] = 0; 
	err_data_al_mpf[i][j][k] = 0; 
	data_al_pTbal[i][j][k] = 0; 
	err_data_al_pTbal[i][j][k] = 0; 

	TString name = name1; name+=count;
	hdata_rel_r[k][j][i] = new TH1D(name,"",100, 0, 2.5);
	name = name2;name+=count;
	hdata_mpf_r[k][j][i] = new TH1D(name,"",100, 0, 2.5);
	name = name3; name+=count;
	hmc_rel_r[k][j][i] = new TH1D(name,"",100, 0, 2.5);
	name = name4; name+=count;
	hmc_mpf_r[k][j][i] = new TH1D(name,"",100, 0, 2.5);
	count++;

      }
    }
  }

 
 // Create the tree reader and its data containers
   TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
   TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
   TTreeReaderValue<Float_t> rel_r_data(myReader_DATA, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_data(myReader_DATA, "mpf_r");
   TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
   
   while (myReader_DATA.Next()) {
     for(int k=0; k<n_pt-1; k++){
   	   if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
	   for(int j=0; j<n_eta_common_2-1; j++){
	     if(fabs(*probejet_eta_data)>eta_common_bins_2[j+1] || fabs(*probejet_eta_data)<eta_common_bins_2[j]) continue;
	     for(int i=0; i<n_alpha_common; i++){
	       if(*alpha_data>alpha_bins_common[i]) continue;
	       else{
		 hdata_rel_r[k][j][i]->Fill(*rel_r_data,*weight_data);
		 hdata_mpf_r[k][j][i]->Fill(*mpf_r_data,*weight_data);
	       }
	   }
	 }
     }
   }

   TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
   TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
   TTreeReaderValue<Float_t> rel_r_mc(myReader_MC, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_mc(myReader_MC, "mpf_r");
   TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
   while (myReader_MC.Next()) {
     for(int k=0; k<n_pt-1; k++){
       if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
       for(int j=0; j<n_eta_common_2-1; j++){
   	 if(fabs(*probejet_eta_mc)>eta_common_bins_2[j+1] || fabs(*probejet_eta_mc)<eta_common_bins_2[j]) continue;
   	 for(int i=0; i<n_alpha_common; i++){
   	   if(*alpha_mc>alpha_bins_common[i]) continue;
   	   else{
   	     hmc_rel_r[k][j][i]->Fill(*rel_r_mc,*weight_mc);
   	     hmc_mpf_r[k][j][i]->Fill(*mpf_r_mc,*weight_mc);
   	   }
   	 }
       }
     }
   }




  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common_2-1; j++){
      for(int k=0; k<n_pt-1; k++){
	pair<double,double> res_mc_rel_r,res_data_rel_r;
	pair<double,double> res_mc_mpf_r,res_data_mpf_r;
	res_mc_rel_r = GetValueAndError(hmc_rel_r[k][j][i]);
	res_data_rel_r = GetValueAndError(hdata_rel_r[k][j][i]);
	res_mc_mpf_r = GetValueAndError(hmc_mpf_r[k][j][i]);
	res_data_mpf_r = GetValueAndError(hdata_mpf_r[k][j][i]);

	mc_al_mpf[i][j][k] = res_mc_mpf_r.first;
	err_mc_al_mpf[i][j][k] = res_mc_mpf_r.second;
	mc_al_pTbal[i][j][k] = res_mc_rel_r.first;
	err_mc_al_pTbal[i][j][k] = res_mc_rel_r.second;
	data_al_mpf[i][j][k] = res_data_mpf_r.first;
	err_data_al_mpf[i][j][k] = res_data_mpf_r.second;
	data_al_pTbal[i][j][k] = res_data_rel_r.first;
	err_data_al_pTbal[i][j][k] = res_data_rel_r.second;
	pair<double,double> ratio_res_mpf = Rmc_to_Rdata(res_mc_mpf_r,res_data_mpf_r);
	ratio_al_mpf[i][j][k] = ratio_res_mpf.first;
	err_ratio_al_mpf[i][j][k] = ratio_res_mpf.second;
	pair<double,double> ratio_res_pTbal = Rmc_to_Rdata(res_mc_rel_r,res_data_rel_r);
	ratio_al_pTbal[i][j][k] = ratio_res_pTbal.first;
	err_ratio_al_pTbal[i][j][k] = ratio_res_pTbal.second;
      }
    }
  }


  //Store results in TGraphErrors
  double res_xbin_tgraph[n_pt-1];// = {(pt_bins[0]+pt_bins[1])/2
  double res_zero[n_pt-1];
  for(int i=0;i<n_pt-1;i++){
    res_xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
    res_zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
  }

  TGraphErrors *mpf_data[n_alpha_common][n_eta_common_2-1];
  TGraphErrors *rrel_data[n_alpha_common][n_eta_common_2-1];
  TGraphErrors *mpf_mc[n_alpha_common][n_eta_common_2-1];
  TGraphErrors *rrel_mc[n_alpha_common][n_eta_common_2-1];
  TGraphErrors *mpf_ratio[n_alpha_common][n_eta_common_2-1];
  TGraphErrors *rrel_ratio[n_alpha_common][n_eta_common_2-1];

  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common_2-1; j++){
      mpf_data[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,data_al_mpf[i][j],res_zero,err_data_al_mpf[i][j]);
      mpf_mc[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,mc_al_mpf[i][j],res_zero,err_mc_al_mpf[i][j]);
      rrel_data[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,data_al_pTbal[i][j],res_zero,err_data_al_pTbal[i][j]);
      rrel_mc[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,mc_al_pTbal[i][j],res_zero,err_mc_al_pTbal[i][j]);
      mpf_ratio[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,ratio_al_mpf[i][j],res_zero,err_ratio_al_mpf[i][j]);
      rrel_ratio[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph,ratio_al_pTbal[i][j],res_zero,err_ratio_al_pTbal[i][j]);

    }
  }

  // Cleaning for empty points (with low statistic)
  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common_2-1; j++){
      mpf_data[i][j] = CleanEmptyPoints(mpf_data[i][j]);
      mpf_mc[i][j] = CleanEmptyPoints(mpf_mc[i][j]);
      rrel_data[i][j] = CleanEmptyPoints(rrel_data[i][j]);
      rrel_mc[i][j] = CleanEmptyPoints(rrel_mc[i][j]);
      mpf_ratio[i][j] = CleanEmptyPoints(mpf_ratio[i][j]);
      rrel_ratio[i][j] = CleanEmptyPoints(rrel_ratio[i][j]);

      mpf_data[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      mpf_mc[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      rrel_data[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      rrel_mc[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      mpf_ratio[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      rrel_ratio[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      }
    }

  CorrectionObject::make_path(std::string((_outpath + "output/").Data()));

  //Save results in root file
  TFile* outputfile = new TFile(CorrectionObject::_outpath+"output/JEC_L2_Dijet_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+"_eta_0_13.root","RECREATE");
  outputfile->Print();

  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common_2-1; j++){
      if(i==0){
	outputfile->mkdir("ratio/"+eta_output_2[j]);
	outputfile->mkdir("data/"+eta_output_2[j]);
	outputfile->mkdir("mc/"+eta_output_2[j]);
      }
      outputfile->cd("ratio/"+eta_output_2[j]);
      mpf_ratio[i][j]->Write();
      rrel_ratio[i][j]->Write();
    
      outputfile->cd("data/"+eta_output_2[j]);
      mpf_data[i][j]->Write();
      rrel_data[i][j]->Write();
      
      outputfile->cd("mc/"+eta_output_2[j]);
      mpf_mc[i][j]->Write();
      rrel_mc[i][j]->Write();
    }
 }

  cout<<"Draw result for alpha = "<<alpha_range_common[3]<<" eta = "<<eta_common_range_2[n_eta_common_2-2]<<" "<<eta_common_range_2[n_eta_common_2-1]<<endl;
  mpf_ratio[3][0]->Draw();

  outputfile->Write();
  outputfile->Close();








  //delete everything



  delete outputfile;

  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common_2-1; j++){
      delete mpf_data[i][j];
      delete mpf_mc[i][j];
      delete rrel_data[i][j];
      delete rrel_mc[i][j];
      delete mpf_ratio[i][j];
      delete rrel_ratio[i][j];

    }
  }




  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common_2-1; j++){
      for(int k=0; k<n_pt-1; k++){
	delete hdata_rel_r[k][j][i];
	delete hdata_mpf_r[k][j][i];
	delete hmc_rel_r[k][j][i];
	delete hmc_mpf_r[k][j][i];

      }
    }
  }



  delete m_gStyle;
}
