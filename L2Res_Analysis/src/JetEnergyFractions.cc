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


using namespace std;

// phi_binned not fully implemented yet
void CorrectionObject::JetEnergyFractions(double abs_asymmetry_cut, bool create_dir, bool phi_binned){
  cout << "--------------- Starting JetEnergyFractions() ---------------" << endl << endl;
  gStyle->SetOptStat(0);
  TString txttag=CorrectionObject::_generator_tag; 
  // TLatex *tex = new TLatex();
  // tex->SetNDC();
  // tex->SetTextSize(0.045);  
  int n_pt_ = max(n_pt,n_pt_HF);
  bool eta_cut_bool;
  int n_pt_cutted;

  if(create_dir) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/");
  if(create_dir) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/");
  if(create_dir) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/EnergyFractions/");
  if(create_dir) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/EnergyFractionsPhiBinned/");

  //Set up histos 
  TH1D *hdata_probejet_neutEmEF[n_pt_-1][n_eta-1]; //neutral EM energy fraction
  TH1D *hdata_probejet_neutHadEF[n_pt_-1][n_eta-1]; //neutral hadron energy fraction
  TH1D *hdata_probejet_chEmEF[n_pt_-1][n_eta-1]; //charged EM energy fraction
  TH1D *hdata_probejet_chHadEF[n_pt_-1][n_eta-1]; //charged hadron energy fraction
  TH1D *hdata_probejet_photonEF[n_pt_-1][n_eta-1]; //photon energy fraction
  TH1D *hdata_probejet_muonEF[n_pt_-1][n_eta-1]; //muon hadron energy fraction
  TH1D *hdata_probejet_phi[n_pt_-1][n_eta-1]; //phi

  TH1D *hmc_probejet_neutEmEF[n_pt_-1][n_eta-1]; //neutral EM energy fraction
  TH1D *hmc_probejet_neutHadEF[n_pt_-1][n_eta-1]; //neutral hadron energy fraction
  TH1D *hmc_probejet_chEmEF[n_pt_-1][n_eta-1]; //charged EM energy fraction
  TH1D *hmc_probejet_chHadEF[n_pt_-1][n_eta-1]; //charged hadron energy fraction
  TH1D *hmc_probejet_photonEF[n_pt_-1][n_eta-1]; //photon energy fraction
  TH1D *hmc_probejet_muonEF[n_pt_-1][n_eta-1]; //muon hadron energy fraction
  TH1D *hmc_probejet_phi[n_pt_-1][n_eta-1]; //phi

  TH1D *hdata_barreljet_neutEmEF[n_pt_-1][n_eta-1]; //neutral EM energy fraction
  TH1D *hdata_barreljet_neutHadEF[n_pt_-1][n_eta-1]; //neutral hadron energy fraction
  TH1D *hdata_barreljet_chEmEF[n_pt_-1][n_eta-1]; //charged EM energy fraction
  TH1D *hdata_barreljet_chHadEF[n_pt_-1][n_eta-1]; //charged hadron energy fraction
  TH1D *hdata_barreljet_photonEF[n_pt_-1][n_eta-1]; //photon energy fraction
  TH1D *hdata_barreljet_muonEF[n_pt_-1][n_eta-1]; //muon hadron energy fraction
  TH1D *hdata_barreljet_phi[n_pt_-1][n_eta-1]; //phi

  TH1D *hmc_barreljet_neutEmEF[n_pt_-1][n_eta-1]; //neutral EM energy fraction
  TH1D *hmc_barreljet_neutHadEF[n_pt_-1][n_eta-1]; //neutral hadron energy fraction
  TH1D *hmc_barreljet_chEmEF[n_pt_-1][n_eta-1]; //charged EM energy fraction
  TH1D *hmc_barreljet_chHadEF[n_pt_-1][n_eta-1]; //charged hadron energy fraction
  TH1D *hmc_barreljet_photonEF[n_pt_-1][n_eta-1]; //photon energy fraction
  TH1D *hmc_barreljet_muonEF[n_pt_-1][n_eta-1]; //muon hadron energy fraction
  TH1D *hmc_barreljet_phi[n_pt_-1][n_eta-1]; //phi


  //Set up histos 
  TH1D *hdata_probejet_neutEmEF_phi[n_pt_-1][n_eta-1][10]; //neutral EM energy fraction
  TH1D *hdata_probejet_neutHadEF_phi[n_pt_-1][n_eta-1][10]; //neutral hadron energy fraction
  TH1D *hdata_probejet_chEmEF_phi[n_pt_-1][n_eta-1][10]; //charged EM energy fraction
  TH1D *hdata_probejet_chHadEF_phi[n_pt_-1][n_eta-1][10]; //charged hadron energy fraction
  TH1D *hdata_probejet_photonEF_phi[n_pt_-1][n_eta-1][10]; //photon energy fraction
  TH1D *hdata_probejet_muonEF_phi[n_pt_-1][n_eta-1][10]; //muon hadron energy fraction
  TH1D *hdata_probejet_phi_phi[n_pt_-1][n_eta-1][10]; //phi

  TH1D *hmc_probejet_neutEmEF_phi[n_pt_-1][n_eta-1][10]; //neutral EM energy fraction
  TH1D *hmc_probejet_neutHadEF_phi[n_pt_-1][n_eta-1][10]; //neutral hadron energy fraction
  TH1D *hmc_probejet_chEmEF_phi[n_pt_-1][n_eta-1][10]; //charged EM energy fraction
  TH1D *hmc_probejet_chHadEF_phi[n_pt_-1][n_eta-1][10]; //charged hadron energy fraction
  TH1D *hmc_probejet_photonEF_phi[n_pt_-1][n_eta-1][10]; //photon energy fraction
  TH1D *hmc_probejet_muonEF_phi[n_pt_-1][n_eta-1][10]; //muon hadron energy fraction
  TH1D *hmc_probejet_phi_phi[n_pt_-1][n_eta-1][10]; //phi
  
  int count = 0;

  TString name9 = "hist_data_probejet_neutEmEF_";
  TString name10 = "hist_mc_probejet_neutEmEF_";
  TString name11 = "hist_data_probejet_neutHadEF_";
  TString name12 = "hist_mc_probejet_neutHadEF_";
  TString name13 = "hist_data_probejet_chEmEF_";
  TString name14 = "hist_mc_probejet_chEmEF_";
  TString name15 = "hist_data_probejet_chHadEF_";
  TString name16 = "hist_mc_probejet_chHadEF_";
  TString name17 = "hist_data_probejet_photonEF_";
  TString name18 = "hist_mc_probejet_photonEF_";
  TString name19 = "hist_data_probejet_muonEF_";
  TString name20 = "hist_mc_probejet_muonEF_";
  TString name21 = "hist_data_probejet_phi_";
  TString name22 = "hist_mc_probejet_phi_";

  TString name9b = "hist_data_barreljet_neutEmEF_";
  TString name10b = "hist_mc_barreljet_neutEmEF_";
  TString name11b = "hist_data_barreljet_neutHadEF_";
  TString name12b = "hist_mc_barreljet_neutHadEF_";
  TString name13b = "hist_data_barreljet_chEmEF_";
  TString name14b = "hist_mc_barreljet_chEmEF_";
  TString name15b = "hist_data_barreljet_chHadEF_";
  TString name16b = "hist_mc_barreljet_chHadEF_";
  TString name17b = "hist_data_barreljet_photonEF_";
  TString name18b = "hist_mc_barreljet_photonEF_";
  TString name19b = "hist_data_barreljet_muonEF_";
  TString name20b = "hist_mc_barreljet_muonEF_";
  TString name21b = "hist_data_barreljet_phi_";
  TString name22b = "hist_mc_barreljet_phi_";
 
  for(int j=0; j<n_eta-1; j++){
    TString eta_name = "eta_"+eta_range2[j]+"_"+eta_range2[j+1];
    eta_cut_bool = fabs(eta_bins[j])>eta_cut;
    n_pt_cutted = ( eta_cut_bool ?  n_pt_HF-2 : n_pt-2 );
    //    for(int k=0; k<n_pt_cutted; k++){
    for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      //      TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];
      TString name;
      name = name9  + eta_name + "_" + pt_name;
      hdata_probejet_neutEmEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name10  + eta_name + "_" + pt_name;
      hmc_probejet_neutEmEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name11 + eta_name + "_" + pt_name;
      hdata_probejet_neutHadEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name12  + eta_name + "_" + pt_name;
      hmc_probejet_neutHadEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name13 + eta_name + "_" + pt_name;
      hdata_probejet_chEmEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name14  + eta_name + "_" + pt_name;
      hmc_probejet_chEmEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name15  + eta_name + "_" + pt_name;
      hdata_probejet_chHadEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name16  + eta_name + "_" + pt_name;
      hmc_probejet_chHadEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name17  + eta_name + "_" + pt_name;
      hdata_probejet_photonEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name18  + eta_name + "_" + pt_name;
      hmc_probejet_photonEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name19  + eta_name + "_" + pt_name;
      hdata_probejet_muonEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name20  + eta_name + "_" + pt_name;
      hmc_probejet_muonEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name21  + eta_name + "_" + pt_name;
      hdata_probejet_phi [k][j] = new TH1D(name,"",30,-3.14,3.14);
      name = name22  + eta_name + "_" + pt_name;
      hmc_probejet_phi [k][j] = new TH1D(name,"",30,-3.14,3.14);

      name = name9b  + eta_name + "_" + pt_name;
      hdata_barreljet_neutEmEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name10b  + eta_name + "_" + pt_name;
      hmc_barreljet_neutEmEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name11b + eta_name + "_" + pt_name;
      hdata_barreljet_neutHadEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name12b  + eta_name + "_" + pt_name;
      hmc_barreljet_neutHadEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name13b + eta_name + "_" + pt_name;
      hdata_barreljet_chEmEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name14b  + eta_name + "_" + pt_name;
      hmc_barreljet_chEmEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name15b  + eta_name + "_" + pt_name;
      hdata_barreljet_chHadEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name16b  + eta_name + "_" + pt_name;
      hmc_barreljet_chHadEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name17b  + eta_name + "_" + pt_name;
      hdata_barreljet_photonEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name18b  + eta_name + "_" + pt_name;
      hmc_barreljet_photonEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name19b  + eta_name + "_" + pt_name;
      hdata_barreljet_muonEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name20b  + eta_name + "_" + pt_name;
      hmc_barreljet_muonEF [k][j] = new TH1D(name,"",30,0,1.1);
      name = name21b  + eta_name + "_" + pt_name;
      hdata_barreljet_phi [k][j] = new TH1D(name,"",30,-3.14,3.14);
      name = name22b  + eta_name + "_" + pt_name;
      hmc_barreljet_phi [k][j] = new TH1D(name,"",30,-3.14,3.14);

      

      if(phi_binned){
	for(int l = 0; l<10; l++){
	  double phi0 = -1.*M_PI;
	  double phid = 2.*M_PI/10.;
	  TString phi_name = "phi_";//+ string(l*phid+phi0) + "_" + string(l*phid+phid+phi0);
	  
	  name = name9  + eta_name + "_" + pt_name + "_" + phi_name;
	  hdata_probejet_neutEmEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name10  + eta_name + "_" + pt_name + "_" + phi_name;
	  hmc_probejet_neutEmEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name11 + eta_name + "_" + pt_name + "_" + phi_name;

	  hdata_probejet_neutHadEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name12  + eta_name + "_" + pt_name + "_" + phi_name;
	  hmc_probejet_neutHadEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name13 + eta_name + "_" + pt_name + "_" + phi_name;

	  hdata_probejet_chEmEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name14  + eta_name + "_" + pt_name + "_" + phi_name;
	  hmc_probejet_chEmEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name15  + eta_name + "_" + pt_name + "_" + phi_name;

	  hdata_probejet_chHadEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name16  + eta_name + "_" + pt_name + "_" + phi_name;
	  hmc_probejet_chHadEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name17  + eta_name + "_" + pt_name + "_" + phi_name;

	  hdata_probejet_photonEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name18  + eta_name + "_" + pt_name + "_" + phi_name;
	  hmc_probejet_photonEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name19  + eta_name + "_" + pt_name + "_" + phi_name;

	  hdata_probejet_muonEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name20  + eta_name + "_" + pt_name + "_" + phi_name;
	  hmc_probejet_muonEF_phi [k][j][l] = new TH1D(name,"",30,0,1.1);
	  name = name21  + eta_name + "_" + pt_name + "_" + phi_name;

	  hdata_probejet_phi_phi [k][j][l] = new TH1D(name,"",30,-3.14,3.14);
	  name = name22  + eta_name + "_" + pt_name + "_" + phi_name;
	  hmc_probejet_phi_phi [k][j][l] = new TH1D(name,"",30,-3.14,3.14);

	}
      }
      
      
      count++;
    }
  }


  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Float_t> probejet_phi_data(myReader_DATA, "probejet_phi");
  TTreeReaderValue<Float_t> probejet_pt_data(myReader_DATA, "probejet_pt");
  TTreeReaderValue<Float_t> barreljet_pt_data(myReader_DATA, "barreljet_pt");
  TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
  TTreeReaderValue<Float_t> asymmetry_data(myReader_DATA, "asymmetry");
  TTreeReaderValue<Float_t> B_data(myReader_DATA, "B");
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  TTreeReaderValue<Float_t> MET_data(myReader_DATA, "MET");
  TTreeReaderValue<Float_t> sum_jets_pt_data(myReader_DATA, "sum_jets_pt");
  TTreeReaderValue<Float_t> jet3_pt_data(myReader_DATA, "jet3_pt");

  
  TTreeReaderValue<Float_t> probejet_neutEmEF_data(myReader_DATA, "probejet_neutEmEF");
  TTreeReaderValue<Float_t> probejet_neutHadEF_data(myReader_DATA, "probejet_neutHadEF");
  TTreeReaderValue<Float_t> probejet_chEmEF_data(myReader_DATA, "probejet_chEmEF");
  TTreeReaderValue<Float_t> probejet_chHadEF_data(myReader_DATA, "probejet_chHadEF");
  TTreeReaderValue<Float_t> probejet_photonEF_data(myReader_DATA, "probejet_photonEF");
  TTreeReaderValue<Float_t> probejet_muonEF_data(myReader_DATA, "probejet_muonEF");
  // TTreeReaderValue<Float_t> probejet_phi_data(myReader_DATA, "probejet_phi");

  TTreeReaderValue<Float_t> barreljet_neutEmEF_data(myReader_DATA, "barreljet_neutEmEF");
  TTreeReaderValue<Float_t> barreljet_neutHadEF_data(myReader_DATA, "barreljet_neutHadEF");
  TTreeReaderValue<Float_t> barreljet_chEmEF_data(myReader_DATA, "barreljet_chEmEF");
  TTreeReaderValue<Float_t> barreljet_chHadEF_data(myReader_DATA, "barreljet_chHadEF");
  TTreeReaderValue<Float_t> barreljet_photonEF_data(myReader_DATA, "barreljet_photonEF");
  TTreeReaderValue<Float_t> barreljet_muonEF_data(myReader_DATA, "barreljet_muonEF");
  TTreeReaderValue<Float_t> barreljet_phi_data(myReader_DATA, "barreljet_phi");


  
  TTreeReaderValue<int> lumibin_data(myReader_DATA, "lumibin");
  
  int myCount = 0;
  int myCount_cut = 0;
  while (myReader_DATA.Next()) {
    if(*alpha_data>alpha_cut){
      myCount_cut++;
      continue;
    }
    if(abs_asymmetry_cut and abs_asymmetry_cut<fabs(*asymmetry_data)){
      myCount_cut++;
      continue;
    }

    //fill histos in bins of pt and eta
    //    for(int k=0; k<n_pt_cutted; k++){
    //    for(int k=n_pt-6; k<n_pt_cutted; k++){//only last 4
    for(int j=0; j<n_eta-1; j++){
      if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
      eta_cut_bool = fabs(eta_bins[j])>eta_cut;
      for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
	//	  if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
	if(*pt_ave_data<(eta_cut_bool?pt_bins_HF:pt_bins)[k] || *pt_ave_data>(eta_cut_bool?pt_bins_HF:pt_bins)[k+1]) continue;
	else{
	  hdata_probejet_neutEmEF[k][j]->Fill(*probejet_neutEmEF_data,*weight_data);
	  hdata_probejet_neutHadEF[k][j]->Fill(*probejet_neutHadEF_data,*weight_data);
	  hdata_probejet_chEmEF[k][j]->Fill(*probejet_chEmEF_data,*weight_data);
	  hdata_probejet_chHadEF[k][j]->Fill(*probejet_chHadEF_data,*weight_data);
	  hdata_probejet_photonEF[k][j]->Fill(*probejet_photonEF_data,*weight_data);
	  hdata_probejet_muonEF[k][j]->Fill(*probejet_muonEF_data,*weight_data);
	  hdata_probejet_phi[k][j]->Fill(*probejet_phi_data,*weight_data);

	  hdata_barreljet_neutEmEF[k][j]->Fill(*barreljet_neutEmEF_data,*weight_data);
	  hdata_barreljet_neutHadEF[k][j]->Fill(*barreljet_neutHadEF_data,*weight_data);
	  hdata_barreljet_chEmEF[k][j]->Fill(*barreljet_chEmEF_data,*weight_data);
	  hdata_barreljet_chHadEF[k][j]->Fill(*barreljet_chHadEF_data,*weight_data);
	  hdata_barreljet_photonEF[k][j]->Fill(*barreljet_photonEF_data,*weight_data);
	  hdata_barreljet_muonEF[k][j]->Fill(*barreljet_muonEF_data,*weight_data);
	  hdata_barreljet_phi[k][j]->Fill(*barreljet_phi_data,*weight_data);

	}
      }
    }
  }

  //DEBUG
  std::cout<<"\ncount data "<<myCount<<"  count cut data "<<myCount_cut<<std::endl<<std::endl;


  //Get relevant information from MC, loop over MC events 
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> probejet_pt_mc(myReader_MC, "probejet_pt");
  TTreeReaderValue<Float_t> barreljet_pt_mc(myReader_MC, "barreljet_pt");
  TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> asymmetry_mc(myReader_MC, "asymmetry");
  TTreeReaderValue<Float_t> B_mc(myReader_MC, "B");
  TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
  TTreeReaderValue<Float_t> MET_mc(myReader_MC, "MET");
  TTreeReaderValue<Float_t> sum_jets_pt_mc(myReader_MC, "sum_jets_pt");
  TTreeReaderValue<Float_t> jet3_pt_mc(myReader_MC, "jet3_pt");

  
  TTreeReaderValue<Float_t> probejet_neutEmEF_mc(myReader_MC, "probejet_neutEmEF");
  TTreeReaderValue<Float_t> probejet_neutHadEF_mc(myReader_MC, "probejet_neutHadEF");
  TTreeReaderValue<Float_t> probejet_chEmEF_mc(myReader_MC, "probejet_chEmEF");
  TTreeReaderValue<Float_t> probejet_chHadEF_mc(myReader_MC, "probejet_chHadEF");
  TTreeReaderValue<Float_t> probejet_photonEF_mc(myReader_MC, "probejet_photonEF");
  TTreeReaderValue<Float_t> probejet_muonEF_mc(myReader_MC, "probejet_muonEF");
  TTreeReaderValue<Float_t> probejet_phi_mc(myReader_MC, "probejet_phi");

  TTreeReaderValue<Float_t> barreljet_neutEmEF_mc(myReader_MC, "barreljet_neutEmEF");
  TTreeReaderValue<Float_t> barreljet_neutHadEF_mc(myReader_MC, "barreljet_neutHadEF");
  TTreeReaderValue<Float_t> barreljet_chEmEF_mc(myReader_MC, "barreljet_chEmEF");
  TTreeReaderValue<Float_t> barreljet_chHadEF_mc(myReader_MC, "barreljet_chHadEF");
  TTreeReaderValue<Float_t> barreljet_photonEF_mc(myReader_MC, "barreljet_photonEF");
  TTreeReaderValue<Float_t> barreljet_muonEF_mc(myReader_MC, "barreljet_muonEF");
  TTreeReaderValue<Float_t> barreljet_phi_mc(myReader_MC, "barreljet_phi");
  
  double myCount_mc = 0.;
  double  myCount_cut_mc = 0.;
  //  int ncount = 0;
  while (myReader_MC.Next()) {
  //  while (myReader_MC.Next() && ncount<1e5) { //TEST
    //    cout<<"ncount = "<<ncount<<endl;
    //    ncount++;
    if(*alpha_mc>alpha_cut) {
      myCount_cut_mc+=*weight_mc;
      continue;
    }
    if(abs_asymmetry_cut and (abs_asymmetry_cut<fabs(*asymmetry_mc))){
      myCount_cut_mc+=*weight_mc;
      continue;
    }


    //fill histos in bins of pt and eta
    //    for(int k=0; k<n_pt_cutted; k++){
    //    for(int k=n_pt-6; k<n_pt_cutted; k++){//only last 4
    for(int j=0; j<n_eta-1; j++){
      if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
      eta_cut_bool = fabs(eta_bins[j])>eta_cut;
      for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
	//      if(*pt_ave_mc<pt_bins[7]) continue; //plot only high pt staff
	//      if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
	if(*pt_ave_mc<(eta_cut_bool?pt_bins_HF:pt_bins)[k] || *pt_ave_mc>(eta_cut_bool?pt_bins_HF:pt_bins)[k+1]) continue;
	else{
	  hmc_probejet_neutEmEF[k][j]->Fill(*probejet_neutEmEF_mc,*weight_mc);
	  hmc_probejet_neutHadEF[k][j]->Fill(*probejet_neutHadEF_mc,*weight_mc);
	  hmc_probejet_chEmEF[k][j]->Fill(*probejet_chEmEF_mc,*weight_mc);
	  hmc_probejet_chHadEF[k][j]->Fill(*probejet_chHadEF_mc,*weight_mc);
	  hmc_probejet_photonEF[k][j]->Fill(*probejet_photonEF_mc,*weight_mc);
	  hmc_probejet_muonEF[k][j]->Fill(*probejet_muonEF_mc,*weight_mc);
	  hmc_probejet_phi[k][j]->Fill(*probejet_phi_mc,*weight_mc);

	  hmc_barreljet_neutEmEF[k][j]->Fill(*barreljet_neutEmEF_mc,*weight_mc);
	  hmc_barreljet_neutHadEF[k][j]->Fill(*barreljet_neutHadEF_mc,*weight_mc);
	  hmc_barreljet_chEmEF[k][j]->Fill(*barreljet_chEmEF_mc,*weight_mc);
	  hmc_barreljet_chHadEF[k][j]->Fill(*barreljet_chHadEF_mc,*weight_mc);
	  hmc_barreljet_photonEF[k][j]->Fill(*barreljet_photonEF_mc,*weight_mc);
	  hmc_barreljet_muonEF[k][j]->Fill(*barreljet_muonEF_mc,*weight_mc);
	  hmc_barreljet_phi[k][j]->Fill(*barreljet_phi_mc,*weight_mc);
	  
	}
      }
    }
  }
  cout<<"Filled MC hists"<<endl;

  // Dump 1-d distributions
  TFile* test_out_mc_B = new TFile(CorrectionObject::_outpath+"plots/control/EF_1d_mc"+ (abs_asymmetry_cut ? "_wAsymCut":"") +".root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    //    for(int k=0; k<n_pt_cutted; k++){     ///k=0 n_pt-1 
      eta_cut_bool = fabs(eta_bins[j])>eta_cut;
      n_pt_cutted = ( eta_cut_bool ?  n_pt_HF-2 : n_pt-2 );
      //      for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      for(int k=0; k<n_pt_cutted; k++){     ///k=0 n_pt-1 

      hmc_probejet_neutEmEF[k][j]->Write();
      hmc_probejet_neutHadEF[k][j]->Write();
      hmc_probejet_chEmEF[k][j]->Write();
      hmc_probejet_chHadEF[k][j]->Write();
      hmc_probejet_photonEF[k][j]->Write();
      hmc_probejet_muonEF[k][j]->Write();
      hmc_probejet_phi[k][j]->Write();

      hmc_barreljet_neutEmEF[k][j]->Write();
      hmc_barreljet_neutHadEF[k][j]->Write();
      hmc_barreljet_chEmEF[k][j]->Write();
      hmc_barreljet_chHadEF[k][j]->Write();
      hmc_barreljet_photonEF[k][j]->Write();
      hmc_barreljet_muonEF[k][j]->Write();
      hmc_barreljet_phi[k][j]->Write();
      
    }
  }
  test_out_mc_B->Close();
  delete test_out_mc_B;
  cout<<"Wrote hists to file"<<endl;
  TFile* test_out_data_B = new TFile(CorrectionObject::_outpath+"plots/control/EF_1d_data"+ (abs_asymmetry_cut ? "_wAsymCut":"") +".root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    eta_cut_bool = fabs(eta_bins[j])>eta_cut;
    n_pt_cutted = ( eta_cut_bool ?  n_pt_HF-2 : n_pt-2 );    
    for(int k=0; k<n_pt_cutted; k++){
     
      hdata_probejet_neutEmEF[k][j]->Write();
      hdata_probejet_neutHadEF[k][j]->Write();
      hdata_probejet_chEmEF[k][j]->Write();
      hdata_probejet_chHadEF[k][j]->Write();
      hdata_probejet_photonEF[k][j]->Write();
      hdata_probejet_muonEF[k][j]->Write();
      hdata_probejet_phi[k][j]->Write();

      hdata_barreljet_neutEmEF[k][j]->Write();
      hdata_barreljet_neutHadEF[k][j]->Write();
      hdata_barreljet_chEmEF[k][j]->Write();
      hdata_barreljet_chHadEF[k][j]->Write();
      hdata_barreljet_photonEF[k][j]->Write();
      hdata_barreljet_muonEF[k][j]->Write();
      hdata_barreljet_phi[k][j]->Write();
      
    }
  }
  test_out_data_B->Close();
  delete test_out_data_B;




  //dummy for tdrCanvas
  TH1D *h = new TH1D("h",";dummy;",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);

  TH1D *hEF = new TH1D("hEF",";dummy;",1000,0,5.191);

  TCanvas* c_0 = new TCanvas();
  tdrCanvas(c_0,"c_0",h,4,10,true,CorrectionObject::_lumitag);

  // for(int i=0; i<n_eta-1; i++){
  //   //Create and fill TGraphErrors
  //   double xbin_tgraph[n_pt-1];
  //   double zero[n_pt-1];
  //   for(int i=0;i<n_pt-1;i++){
  //     xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
  //     zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
  //   }
  // }


  //********************************************************************  Plot all Control Hists *******************************************************************
  
  //************************* Different energy fractions **************************************************************************************
    
  TFile* f_mpf_mc = new TFile(CorrectionObject::_outpath+"plots/control/EF_1d_mc"+ (abs_asymmetry_cut ? "_wAsymCut":"")+".root","READ");
  TFile* f_mpf_data = new TFile(CorrectionObject::_outpath+"plots/control/EF_1d_data"+ (abs_asymmetry_cut ? "_wAsymCut":"") + ".root","READ");
  for(int i=0; i<n_eta-1; i++){
    TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];    
    TLatex *tex = new TLatex();                                                                                                                                       
    tex->SetNDC();                    
    tex->SetTextSize(0.045);
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1]; 
    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 
    eta_cut_bool = fabs(eta_bins[i])>eta_cut;
    n_pt_cutted = ( eta_cut_bool ?  n_pt_HF-2 : n_pt-2 );
    TString PFcomp[6] = {"neutEmEF","neutHadEF","chEmEF","chHadEF","muonEF","phi"};
    for(int ipf = 0; ipf<6;ipf++){//loop over PF fractions

      TCanvas* c7 = new TCanvas();
      tdrCanvas(c7,"c7",hEF,4,10,kSquare,"MC");
      TLegend leg7 = tdrLeg(0.17,0.6,0.85,0.81);
      leg7.SetNColumns(2);

      TH1D* htemp_barreljet_neutEmEF_mc;
      TH1D* htemp_probejet_neutEmEF_mc;

      gPad->SetLogy();
      //      cout<<"n_pt_cutted = "<<n_pt_cutted<<endl;
      //      for(int j=n_pt_cutted-4; j<n_pt_cutted; j++){//only last 4
      for(int j=0; j<n_pt_cutted; j++){//all
	//      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
	TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
	TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
	TString name_probejet_neutEmEF_mc = "hist_mc_probejet_"+PFcomp[ipf]+"_"+eta_name+"_"+pt_name;
	htemp_probejet_neutEmEF_mc = (TH1D*)f_mpf_mc->Get(name_probejet_neutEmEF_mc);
	//      htemp_probejet_neutEmEF_mc->Print();
	int n_ev =  htemp_probejet_neutEmEF_mc->GetEntries();
	if(htemp_probejet_neutEmEF_mc->Integral() > 0)htemp_probejet_neutEmEF_mc->Scale(1/htemp_probejet_neutEmEF_mc->Integral());
	hEF->GetXaxis()->SetTitle("probejet "+PFcomp[ipf]);
	hEF->GetYaxis()->SetTitle("Norm. Entries");
	hEF->GetYaxis()->SetTitleOffset(1.5);
	// h->SetMaximum(0.3);
	hEF->GetXaxis()->SetLimits(0,1.1);
	//      h->GetYaxis()->SetLimits(0,0.8);
	hEF->SetMaximum(1.3);
	// hEF->SetMaximum(3);
	hEF->SetMinimum(0.001);
	if(j<9) htemp_probejet_neutEmEF_mc->SetLineColor(j+1);
	else    htemp_probejet_neutEmEF_mc->SetLineColor(j-8);
	htemp_probejet_neutEmEF_mc->SetLineWidth(2+j*0.2);
	if(n_ev>100) htemp_probejet_neutEmEF_mc->Draw("HIST SAME");
	leg7.AddEntry(htemp_probejet_neutEmEF_mc, legname,"l");
      }

      leg7.Draw();
      tex->DrawLatex(0.47,0.85,"MC, " + text);
      tex->DrawLatex(0.15,0.95,txttag);
      c7->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_"+PFcomp[ipf]+"_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");

      TCanvas* c7b = new TCanvas();
      tdrCanvas(c7b,"c7b",hEF,4,10,kSquare,"MC");
      TLegend leg7b = tdrLeg(0.17,0.6,0.85,0.81);
      leg7b.SetNColumns(2);
      gPad->SetLogy();
      //      for(int j=n_pt_cutted-4; j<n_pt_cutted; j++){//only last 4
      for(int j=0; j<n_pt_cutted; j++){//all
	TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
	TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
	TString name_barreljet_neutEmEF_mc = "hist_mc_barreljet_"+PFcomp[ipf]+"_"+eta_name+"_"+pt_name;
	htemp_barreljet_neutEmEF_mc = (TH1D*)f_mpf_mc->Get(name_barreljet_neutEmEF_mc);

	//      htemp_barreljet_neutEmEF_mc->Print();
	int n_ev =  htemp_barreljet_neutEmEF_mc->GetEntries();
	if(htemp_barreljet_neutEmEF_mc->Integral() > 0)htemp_barreljet_neutEmEF_mc->Scale(1/htemp_barreljet_neutEmEF_mc->Integral());
	hEF->GetXaxis()->SetTitle("barreljet "+PFcomp[ipf]);
	hEF->GetYaxis()->SetTitle("Norm. Entries");
	hEF->GetYaxis()->SetTitleOffset(1.5);
	hEF->GetXaxis()->SetLimits(0,1.1);
	hEF->SetMaximum(1.3);
	hEF->SetMinimum(0.001);
	if(j<9) htemp_barreljet_neutEmEF_mc->SetLineColor(j+1);
	else    htemp_barreljet_neutEmEF_mc->SetLineColor(j-8);
	htemp_barreljet_neutEmEF_mc->SetLineWidth(2+j*0.2);
	if(n_ev>100) htemp_barreljet_neutEmEF_mc->Draw("HIST SAME");
	leg7b.AddEntry(htemp_barreljet_neutEmEF_mc, legname,"l");
      }

      leg7b.Draw();
      tex->DrawLatex(0.47,0.85,"MC, " + text);
      tex->DrawLatex(0.15,0.95,txttag);
      c7b->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/barreljet_"+PFcomp[ipf]+"_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");

      TCanvas* c7r = new TCanvas();
      tdrCanvas(c7r,"c7r",hEF,4,10,kSquare,"MC");
      TLegend leg7r = tdrLeg(0.17,0.6,0.85,0.81);
      leg7r.SetNColumns(2);
      gPad->SetLogy();
      //      for(int j=n_pt_cutted-4; j<n_pt_cutted; j++){//only last 4
      for(int j=0; j<n_pt_cutted; j++){//all
	TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
	TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
	TString name_barreljet_neutEmEF_mc = "hist_mc_barreljet_"+PFcomp[ipf]+"_"+eta_name+"_"+pt_name;
	htemp_barreljet_neutEmEF_mc = (TH1D*)f_mpf_mc->Get(name_barreljet_neutEmEF_mc);

	TString name_probejet_neutEmEF_mc = "hist_mc_probejet_"+PFcomp[ipf]+"_"+eta_name+"_"+pt_name;
	htemp_probejet_neutEmEF_mc = (TH1D*)f_mpf_mc->Get(name_probejet_neutEmEF_mc);

	//      htemp_barreljet_neutEmEF_mc->Print();
	int n_ev =  htemp_barreljet_neutEmEF_mc->GetEntries();
	if(htemp_barreljet_neutEmEF_mc->Integral() > 0)htemp_barreljet_neutEmEF_mc->Scale(1/htemp_barreljet_neutEmEF_mc->Integral());
	if(htemp_probejet_neutEmEF_mc->Integral() > 0)htemp_probejet_neutEmEF_mc->Scale(1/htemp_probejet_neutEmEF_mc->Integral());
	if(htemp_barreljet_neutEmEF_mc->Integral() > 0)htemp_probejet_neutEmEF_mc->Divide(htemp_barreljet_neutEmEF_mc);

	hEF->GetXaxis()->SetTitle("probejet/barreljet "+PFcomp[ipf]);
	hEF->GetYaxis()->SetTitle("Norm. Entries");
	hEF->GetYaxis()->SetTitleOffset(1.5);
	hEF->GetXaxis()->SetLimits(0,1.1);
	hEF->SetMaximum(1000.3);
	hEF->SetMinimum(0.001);
	if(j<9) htemp_barreljet_neutEmEF_mc->SetLineColor(j+1);
	else    htemp_barreljet_neutEmEF_mc->SetLineColor(j-8);
	htemp_probejet_neutEmEF_mc->SetLineWidth(2+j*0.2);
	if(n_ev>100) htemp_probejet_neutEmEF_mc->Draw("HIST SAME");
	leg7r.AddEntry(htemp_probejet_neutEmEF_mc, legname,"l");
      }

      leg7r.Draw();
      tex->DrawLatex(0.47,0.85,"MC, " + text);
      tex->DrawLatex(0.15,0.95,txttag);
      c7r->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probebarrelratio_"+PFcomp[ipf]+"_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");


      TCanvas* c8 = new TCanvas();
      tdrCanvas(c8,"c8",hEF,4,10,kSquare,"DATA");
      //    TLegend leg8 = tdrLeg(0.62,0.46,0.85,0.81);
      TLegend leg8 = tdrLeg(0.17,0.6,0.85,0.81);
      leg8.SetNColumns(2);
      TH1D* htemp_probejet_neutEmEF_data;
      TH1D* htemp_barreljet_neutEmEF_data;
      gPad->SetLogy();
      //      for(int j=n_pt_cutted-4; j<n_pt_cutted; j++){//only last 4
      for(int j=0; j<n_pt_cutted; j++){//all
	TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
	TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
	TString name_probejet_neutEmEF_data = "hist_data_probejet_"+PFcomp[ipf]+"_"+eta_name+"_"+pt_name;
	htemp_probejet_neutEmEF_data = (TH1D*)f_mpf_data->Get(name_probejet_neutEmEF_data);
	//      htemp_probejet_neutEmEF_data->Print();
	int n_ev =  htemp_probejet_neutEmEF_data->GetEntries();
	if(htemp_probejet_neutEmEF_data->Integral() > 0)htemp_probejet_neutEmEF_data->Scale(1/htemp_probejet_neutEmEF_data->Integral());
	hEF->GetXaxis()->SetTitle("probejet "+PFcomp[ipf]);
	hEF->GetYaxis()->SetTitle("Norm. Entries");
	hEF->GetYaxis()->SetTitleOffset(1.5);
	hEF->GetXaxis()->SetLimits(0,1.1);
	hEF->SetMaximum(1.3);
	hEF->SetMinimum(0.001);
	if(j<9) htemp_probejet_neutEmEF_data->SetLineColor(j+1);
	else    htemp_probejet_neutEmEF_data->SetLineColor(j-8);      htemp_probejet_neutEmEF_data->SetLineWidth(2+j*0.2);
	if(n_ev>100) htemp_probejet_neutEmEF_data->Draw("HIST SAME");
	leg8.AddEntry(htemp_probejet_neutEmEF_data, legname);
      }

      leg8.Draw();
      tex->DrawLatex(0.47,0.85,"Data, " + text);
      c8->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_"+PFcomp[ipf]+"_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");

      TCanvas* c8b = new TCanvas();
      tdrCanvas(c8b,"c8",hEF,4,10,kSquare,"DATA");
      TLegend leg8b = tdrLeg(0.17,0.6,0.85,0.81);
      leg8b.SetNColumns(2);
      gPad->SetLogy();
      //      for(int j=n_pt_cutted-4; j<n_pt_cutted; j++){//only last 4
      for(int j=0; j<n_pt_cutted; j++){//all
	TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
	TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
	TString name_barreljet_neutEmEF_data = "hist_data_barreljet_"+PFcomp[ipf]+"_"+eta_name+"_"+pt_name;
	htemp_barreljet_neutEmEF_data = (TH1D*)f_mpf_data->Get(name_barreljet_neutEmEF_data);
	//      htemp_barreljet_neutEmEF_data->Print();
	int n_ev =  htemp_barreljet_neutEmEF_data->GetEntries();
	if(htemp_barreljet_neutEmEF_data->Integral() > 0)htemp_barreljet_neutEmEF_data->Scale(1/htemp_barreljet_neutEmEF_data->Integral());
	hEF->GetXaxis()->SetTitle("barreljet "+PFcomp[ipf]);
	hEF->GetYaxis()->SetTitle("Norm. Entries");
	hEF->GetYaxis()->SetTitleOffset(1.5);
	hEF->GetXaxis()->SetLimits(0,1.1);
	hEF->SetMaximum(1.3);
	hEF->SetMinimum(0.001);
	if(j<9) htemp_barreljet_neutEmEF_data->SetLineColor(j+1);
	else    htemp_barreljet_neutEmEF_data->SetLineColor(j-8);      htemp_barreljet_neutEmEF_data->SetLineWidth(2+j*0.2);
	if(n_ev>100) htemp_barreljet_neutEmEF_data->Draw("HIST SAME");
	leg8b.AddEntry(htemp_barreljet_neutEmEF_data, legname);
      }

      leg8b.Draw();
      tex->DrawLatex(0.47,0.85,"Data, " + text);
      c8b->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/barreljet_"+PFcomp[ipf]+"_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");
      TCanvas* c8r = new TCanvas();
      tdrCanvas(c8r,"c8r",hEF,4,10,kSquare,"DATA");
      TLegend leg8r = tdrLeg(0.17,0.6,0.85,0.81);
      leg8r.SetNColumns(2);
      gPad->SetLogy();
      //      for(int j=n_pt_cutted-4; j<n_pt_cutted; j++){//only last 4
      for(int j=0; j<n_pt_cutted; j++){//all
	TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
	//      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
	TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
	TString name_probejet_neutEmEF_data = "hist_data_probejet_"+PFcomp[ipf]+"_"+eta_name+"_"+pt_name;
	htemp_probejet_neutEmEF_data = (TH1D*)f_mpf_data->Get(name_probejet_neutEmEF_data);
	TString name_barreljet_neutEmEF_data = "hist_data_barreljet_"+PFcomp[ipf]+"_"+eta_name+"_"+pt_name;
	htemp_barreljet_neutEmEF_data = (TH1D*)f_mpf_data->Get(name_barreljet_neutEmEF_data);

	int n_ev =  htemp_probejet_neutEmEF_data->GetEntries();
	if(htemp_probejet_neutEmEF_data->Integral() > 0)htemp_probejet_neutEmEF_data->Scale(1/htemp_probejet_neutEmEF_data->Integral());
	if(htemp_barreljet_neutEmEF_data->Integral() > 0) htemp_barreljet_neutEmEF_data->Scale(1/htemp_barreljet_neutEmEF_data->Integral());
	if(htemp_barreljet_neutEmEF_data->Integral() > 0) htemp_probejet_neutEmEF_data->Divide(htemp_barreljet_neutEmEF_data);
	hEF->GetXaxis()->SetTitle("probejet/barreljet "+PFcomp[ipf]);
	hEF->GetYaxis()->SetTitle("Norm. Entries");
	hEF->GetYaxis()->SetTitleOffset(1.5);
	hEF->GetXaxis()->SetLimits(0,1.1);
	hEF->SetMaximum(1000.3);
	hEF->SetMinimum(0.001);
	if(j<9) htemp_probejet_neutEmEF_data->SetLineColor(j+1);
	else    htemp_probejet_neutEmEF_data->SetLineColor(j-8);      htemp_probejet_neutEmEF_data->SetLineWidth(2+j*0.2);
	if(n_ev>100) htemp_probejet_neutEmEF_data->Draw("HIST SAME");
	leg8r.AddEntry(htemp_probejet_neutEmEF_data, legname);
      }

      leg8r.Draw();
      tex->DrawLatex(0.47,0.85,"Data, " + text);
      c8r->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probebarrelratio_"+PFcomp[ipf]+"_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");

    }//loop over PF fractions


    //END Different energy fractions
    delete tex;
   
  }






}
