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


using namespace std;

// phi_binned not fully implemented yet
void CorrectionObject::JetEnergyFractions(double abs_asymmetry_cut, bool create_dir, bool phi_binned){
  cout << "--------------- Starting JetEnergyFractions() ---------------" << endl << endl;
  gStyle->SetOptStat(0);
  TString txttag=CorrectionObject::_generator_tag; 
  // TLatex *tex = new TLatex();
  // tex->SetNDC();
  // tex->SetTextSize(0.045);  

  if(create_dir) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/");
  if(create_dir) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/");
  if(create_dir) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/EnergyFractions/");
  if(create_dir) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/EnergyFractionsPhiBinned/");

  //Set up histos 
  TH1D *hdata_probejet_neutEmEF[n_pt-1][n_eta-1]; //neutral EM energy fraction
  TH1D *hdata_probejet_neutHadEF[n_pt-1][n_eta-1]; //neutral hadron energy fraction
  TH1D *hdata_probejet_chEmEF[n_pt-1][n_eta-1]; //charged EM energy fraction
  TH1D *hdata_probejet_chHadEF[n_pt-1][n_eta-1]; //charged hadron energy fraction
  TH1D *hdata_probejet_photonEF[n_pt-1][n_eta-1]; //photon energy fraction
  TH1D *hdata_probejet_muonEF[n_pt-1][n_eta-1]; //muon hadron energy fraction
  TH1D *hdata_probejet_phi[n_pt-1][n_eta-1]; //phi

  TH1D *hmc_probejet_neutEmEF[n_pt-1][n_eta-1]; //neutral EM energy fraction
  TH1D *hmc_probejet_neutHadEF[n_pt-1][n_eta-1]; //neutral hadron energy fraction
  TH1D *hmc_probejet_chEmEF[n_pt-1][n_eta-1]; //charged EM energy fraction
  TH1D *hmc_probejet_chHadEF[n_pt-1][n_eta-1]; //charged hadron energy fraction
  TH1D *hmc_probejet_photonEF[n_pt-1][n_eta-1]; //photon energy fraction
  TH1D *hmc_probejet_muonEF[n_pt-1][n_eta-1]; //muon hadron energy fraction
  TH1D *hmc_probejet_phi[n_pt-1][n_eta-1]; //phi


  //Set up histos 
  TH1D *hdata_probejet_neutEmEF_phi[n_pt-1][n_eta-1][10]; //neutral EM energy fraction
  TH1D *hdata_probejet_neutHadEF_phi[n_pt-1][n_eta-1][10]; //neutral hadron energy fraction
  TH1D *hdata_probejet_chEmEF_phi[n_pt-1][n_eta-1][10]; //charged EM energy fraction
  TH1D *hdata_probejet_chHadEF_phi[n_pt-1][n_eta-1][10]; //charged hadron energy fraction
  TH1D *hdata_probejet_photonEF_phi[n_pt-1][n_eta-1][10]; //photon energy fraction
  TH1D *hdata_probejet_muonEF_phi[n_pt-1][n_eta-1][10]; //muon hadron energy fraction
  TH1D *hdata_probejet_phi_phi[n_pt-1][n_eta-1][10]; //phi

  TH1D *hmc_probejet_neutEmEF_phi[n_pt-1][n_eta-1][10]; //neutral EM energy fraction
  TH1D *hmc_probejet_neutHadEF_phi[n_pt-1][n_eta-1][10]; //neutral hadron energy fraction
  TH1D *hmc_probejet_chEmEF_phi[n_pt-1][n_eta-1][10]; //charged EM energy fraction
  TH1D *hmc_probejet_chHadEF_phi[n_pt-1][n_eta-1][10]; //charged hadron energy fraction
  TH1D *hmc_probejet_photonEF_phi[n_pt-1][n_eta-1][10]; //photon energy fraction
  TH1D *hmc_probejet_muonEF_phi[n_pt-1][n_eta-1][10]; //muon hadron energy fraction
  TH1D *hmc_probejet_phi_phi[n_pt-1][n_eta-1][10]; //phi
  
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
 
  for(int j=0; j<n_eta-1; j++){
    TString eta_name = "eta_"+eta_range2[j]+"_"+eta_range2[j+1];
    for(int k=0; k<n_pt-1; k++){
      TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];
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
    //    for(int k=0; k<n_pt-1; k++){
    for(int k=n_pt-6; k<n_pt-1; k++){//only last 4
      //      if(*pt_ave_data<pt_bins[7]) continue; //plot only high pt staff
      if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
      for(int j=0; j<n_eta-1; j++){
	if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
	else{
	  hdata_probejet_neutEmEF[k][j]->Fill(*probejet_neutEmEF_data,*weight_data);
	  hdata_probejet_neutHadEF[k][j]->Fill(*probejet_neutHadEF_data,*weight_data);
	  hdata_probejet_chEmEF[k][j]->Fill(*probejet_chEmEF_data,*weight_data);
	  hdata_probejet_chHadEF[k][j]->Fill(*probejet_chHadEF_data,*weight_data);
	  hdata_probejet_photonEF[k][j]->Fill(*probejet_photonEF_data,*weight_data);
	  hdata_probejet_muonEF[k][j]->Fill(*probejet_muonEF_data,*weight_data);
	  hdata_probejet_phi[k][j]->Fill(*probejet_phi_data,*weight_data);
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
  
  double myCount_mc = 0.;
  double  myCount_cut_mc = 0.;
  while (myReader_MC.Next()) {
    if(*alpha_mc>alpha_cut) {
      myCount_cut_mc+=*weight_mc;
      continue;
    }
    if(abs_asymmetry_cut and (abs_asymmetry_cut<fabs(*asymmetry_mc))){
      myCount_cut_mc+=*weight_mc;
      continue;
    }


    //fill histos in bins of pt and eta
    //    for(int k=0; k<n_pt-1; k++){
    for(int k=n_pt-6; k<n_pt-1; k++){//only last 4
      //      if(*pt_ave_mc<pt_bins[7]) continue; //plot only high pt staff
      if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
      for(int j=0; j<n_eta-1; j++){
	if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
	else{	  
	  hmc_probejet_neutEmEF[k][j]->Fill(*probejet_neutEmEF_mc,*weight_mc);
	  hmc_probejet_neutHadEF[k][j]->Fill(*probejet_neutHadEF_mc,*weight_mc);
	  hmc_probejet_chEmEF[k][j]->Fill(*probejet_chEmEF_mc,*weight_mc);
	  hmc_probejet_chHadEF[k][j]->Fill(*probejet_chHadEF_mc,*weight_mc);
	  hmc_probejet_photonEF[k][j]->Fill(*probejet_photonEF_mc,*weight_mc);
	  hmc_probejet_muonEF[k][j]->Fill(*probejet_muonEF_mc,*weight_mc);
	  hmc_probejet_phi[k][j]->Fill(*probejet_phi_mc,*weight_mc);
	  
	}
      }
    }
  }
  

  // Dump 1-d distributions
  TFile* test_out_mc_B = new TFile(CorrectionObject::_outpath+"plots/control/EF_1d_mc"+ (abs_asymmetry_cut ? "_wAsymCut":"") +".root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){     ///k=0 n_pt-1 
      
      hmc_probejet_neutEmEF[k][j]->Write();
      hmc_probejet_neutHadEF[k][j]->Write();
      hmc_probejet_chEmEF[k][j]->Write();
      hmc_probejet_chHadEF[k][j]->Write();
      hmc_probejet_photonEF[k][j]->Write();
      hmc_probejet_muonEF[k][j]->Write();
      hmc_probejet_phi[k][j]->Write();
      
    }
  }
  test_out_mc_B->Close();
  delete test_out_mc_B;

  TFile* test_out_data_B = new TFile(CorrectionObject::_outpath+"plots/control/EF_1d_data"+ (abs_asymmetry_cut ? "_wAsymCut":"") +".root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      
      hdata_probejet_neutEmEF[k][j]->Write();
      hdata_probejet_neutHadEF[k][j]->Write();
      hdata_probejet_chEmEF[k][j]->Write();
      hdata_probejet_chHadEF[k][j]->Write();
      hdata_probejet_photonEF[k][j]->Write();
      hdata_probejet_muonEF[k][j]->Write();
      hdata_probejet_phi[k][j]->Write();
      
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
    TCanvas* c7 = new TCanvas();
    tdrCanvas(c7,"c7",hEF,4,10,kSquare,"MC");
    TLegend leg7 = tdrLeg(0.17,0.6,0.85,0.81);
    leg7.SetNColumns(2);

    TH1D* htemp_probejet_neutEmEF_mc;

    gPad->SetLogy();
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_neutEmEF_mc = "hist_mc_probejet_neutEmEF_"+eta_name+"_"+pt_name;
      htemp_probejet_neutEmEF_mc = (TH1D*)f_mpf_mc->Get(name_probejet_neutEmEF_mc);
      //      htemp_probejet_neutEmEF_mc->Print();
      int n_ev =  htemp_probejet_neutEmEF_mc->GetEntries();
      if(htemp_probejet_neutEmEF_mc->Integral() > 0)htemp_probejet_neutEmEF_mc->Scale(1/htemp_probejet_neutEmEF_mc->Integral());
      hEF->GetXaxis()->SetTitle("probejet neutralEmEF");
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
    c7->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_neutEmEF_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");



    TCanvas* c8 = new TCanvas();
    tdrCanvas(c8,"c8",hEF,4,10,kSquare,"DATA");
    //    TLegend leg8 = tdrLeg(0.62,0.46,0.85,0.81);
    TLegend leg8 = tdrLeg(0.17,0.6,0.85,0.81);
    leg8.SetNColumns(2);
    TH1D* htemp_probejet_neutEmEF_data;

    gPad->SetLogy();
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_neutEmEF_data = "hist_data_probejet_neutEmEF_"+eta_name+"_"+pt_name;
      htemp_probejet_neutEmEF_data = (TH1D*)f_mpf_data->Get(name_probejet_neutEmEF_data);
      //      htemp_probejet_neutEmEF_data->Print();
      int n_ev =  htemp_probejet_neutEmEF_data->GetEntries();
      if(htemp_probejet_neutEmEF_data->Integral() > 0)htemp_probejet_neutEmEF_data->Scale(1/htemp_probejet_neutEmEF_data->Integral());
      hEF->GetXaxis()->SetTitle("probejet neutralEmEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      hEF->SetMaximum(0.5);
      // hEF->SetMaximum(3);
      hEF->SetMinimum(0.001);
      if(j<9) htemp_probejet_neutEmEF_data->SetLineColor(j+1);
      else    htemp_probejet_neutEmEF_data->SetLineColor(j-8);      htemp_probejet_neutEmEF_data->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_neutEmEF_data->Draw("HIST SAME");
      leg8.AddEntry(htemp_probejet_neutEmEF_data, legname);
    }

    leg8.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    c8->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_neutEmEF_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");


    TCanvas* c9 = new TCanvas();
    tdrCanvas(c9,"c9",hEF,4,10,kSquare,"MC");
    TLegend leg9 = tdrLeg(0.17,0.6,0.85,0.81);
    leg9.SetNColumns(2);
    TH1D* htemp_probejet_neutHadEF_mc;

    gPad->SetLogy();
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_neutHadEF_mc = "hist_mc_probejet_neutHadEF_"+eta_name+"_"+pt_name;
      htemp_probejet_neutHadEF_mc = (TH1D*)f_mpf_mc->Get(name_probejet_neutHadEF_mc);
      //      htemp_probejet_neutHadEF_mc->Print();
      int n_ev =  htemp_probejet_neutHadEF_mc->GetEntries();
      if(htemp_probejet_neutHadEF_mc->Integral() > 0)htemp_probejet_neutHadEF_mc->Scale(1/htemp_probejet_neutHadEF_mc->Integral());
      hEF->GetXaxis()->SetTitle("probejet neutralHadEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMinimum(0.001);
      hEF->SetMaximum(1.3);
      if(j<9) htemp_probejet_neutHadEF_mc->SetLineColor(j+1);
      else    htemp_probejet_neutHadEF_mc->SetLineColor(j-8);
      htemp_probejet_neutHadEF_mc->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_neutHadEF_mc->Draw("HIST SAME");
      leg9.AddEntry(htemp_probejet_neutHadEF_mc, legname,"l");
    }

    leg9.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    c9->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_neutHadEF_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + (abs_asymmetry_cut ? "_wAsymCut":"") +".pdf");



    TCanvas* c10 = new TCanvas();
    tdrCanvas(c10,"c10",hEF,4,10,kSquare,"DATA");
    TLegend leg10 = tdrLeg(0.17,0.6,0.85,0.81);
    leg10.SetNColumns(2);
    TH1D* htemp_probejet_neutHadEF_data;

    gPad->SetLogy();
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_neutHadEF_data = "hist_data_probejet_neutHadEF_"+eta_name+"_"+pt_name;
      htemp_probejet_neutHadEF_data = (TH1D*)f_mpf_data->Get(name_probejet_neutHadEF_data);
      //      htemp_probejet_neutHadEF_data->Print();
      int n_ev =  htemp_probejet_neutHadEF_data->GetEntries();
      if(htemp_probejet_neutHadEF_data->Integral() > 0)htemp_probejet_neutHadEF_data->Scale(1/htemp_probejet_neutHadEF_data->Integral());
      hEF->GetXaxis()->SetTitle("probejet neutralHadEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMinimum(0.001);
      hEF->SetMaximum(1.3);
      if(j<9) htemp_probejet_neutHadEF_data->SetLineColor(j+1);
      else    htemp_probejet_neutHadEF_data->SetLineColor(j-8);
      htemp_probejet_neutHadEF_data->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_neutHadEF_data->Draw("HIST SAME");
      leg10.AddEntry(htemp_probejet_neutHadEF_data, legname);
    }

    leg10.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    c10->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_neutHadEF_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + ".pdf");



    TCanvas* c11 = new TCanvas();
    tdrCanvas(c11,"c11",hEF,4,10,kSquare,"MC");
    TLegend leg11 = tdrLeg(0.17,0.6,0.85,0.81);
    leg11.SetNColumns(2);

    TH1D* htemp_probejet_chEmEF_mc;

    gPad->SetLogy();
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_chEmEF_mc = "hist_mc_probejet_chEmEF_"+eta_name+"_"+pt_name;
      htemp_probejet_chEmEF_mc = (TH1D*)f_mpf_mc->Get(name_probejet_chEmEF_mc);
      //      htemp_probejet_chEmEF_mc->Print();
      int n_ev =  htemp_probejet_chEmEF_mc->GetEntries();
      if(htemp_probejet_chEmEF_mc->Integral() > 0)htemp_probejet_chEmEF_mc->Scale(1/htemp_probejet_chEmEF_mc->Integral());
      hEF->GetXaxis()->SetTitle("probejet chEmEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //  hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMaximum(0.5);
      hEF->SetMinimum(0.001);
      if(j<9) htemp_probejet_chEmEF_mc->SetLineColor(j+1);
      else    htemp_probejet_chEmEF_mc->SetLineColor(j-8);
      htemp_probejet_chEmEF_mc->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_chEmEF_mc->Draw("HIST SAME");
      leg11.AddEntry(htemp_probejet_chEmEF_mc, legname,"l");
    }

    leg11.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    c11->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_chEmEF_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + ".pdf");

    TCanvas* c12 = new TCanvas();
    tdrCanvas(c12,"c12",hEF,4,10,kSquare,"DATA");
    TLegend leg12 = tdrLeg(0.17,0.6,0.85,0.81);
    leg12.SetNColumns(2);

    TH1D* htemp_probejet_chEmEF_data;

    gPad->SetLogy();
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_chEmEF_data = "hist_data_probejet_chEmEF_"+eta_name+"_"+pt_name;
      htemp_probejet_chEmEF_data = (TH1D*)f_mpf_data->Get(name_probejet_chEmEF_data);
      //      htemp_probejet_chEmEF_data->Print();
      int n_ev =  htemp_probejet_chEmEF_data->GetEntries();
      if(htemp_probejet_chEmEF_data->Integral() > 0)htemp_probejet_chEmEF_data->Scale(1/htemp_probejet_chEmEF_data->Integral());
      hEF->GetXaxis()->SetTitle("probejet chEmEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMaximum(0.5);
      hEF->SetMinimum(0.001);
      if(j<9) htemp_probejet_chEmEF_data->SetLineColor(j+1);
      else    htemp_probejet_chEmEF_data->SetLineColor(j-8);
      htemp_probejet_chEmEF_data->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_chEmEF_data->Draw("HIST SAME");
      leg12.AddEntry(htemp_probejet_chEmEF_data, legname);
    }

    leg12.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    c12->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_chEmEF_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")  + ".pdf");


    TCanvas* c13 = new TCanvas();
    tdrCanvas(c13,"c13",hEF,4,10,kSquare,"MC");
    TLegend leg13 = tdrLeg(0.17,0.6,0.85,0.81);
    leg13.SetNColumns(2);
    TH1D* htemp_probejet_chHadEF_mc;

    gPad->SetLogy();
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_chHadEF_mc = "hist_mc_probejet_chHadEF_"+eta_name+"_"+pt_name;
      htemp_probejet_chHadEF_mc = (TH1D*)f_mpf_mc->Get(name_probejet_chHadEF_mc);
      //      htemp_probejet_chHadEF_mc->Print();
      int n_ev =  htemp_probejet_chHadEF_mc->GetEntries();
      if(htemp_probejet_chHadEF_mc->Integral() > 0)htemp_probejet_chHadEF_mc->Scale(1/htemp_probejet_chHadEF_mc->Integral());
      hEF->GetXaxis()->SetTitle("probejet chHadEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMaximum(0.5);
      hEF->SetMinimum(0.001);
      if(j<9) htemp_probejet_chHadEF_mc->SetLineColor(j+1);
      else    htemp_probejet_chHadEF_mc->SetLineColor(j-8);
      htemp_probejet_chHadEF_mc->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_chHadEF_mc->Draw("HIST SAME");
      leg13.AddEntry(htemp_probejet_chHadEF_mc, legname,"l");
    }

    leg13.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    c13->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_chHadEF_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + (abs_asymmetry_cut ? "_wAsymCut":"") +".pdf");

    TCanvas* c14 = new TCanvas();
    tdrCanvas(c14,"c14",hEF,4,10,kSquare,"DATA");
    TLegend leg14 = tdrLeg(0.17,0.6,0.85,0.81);
    leg14.SetNColumns(2);
    TH1D* htemp_probejet_chHadEF_data;

    gPad->SetLogy();
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_chHadEF_data = "hist_data_probejet_chHadEF_"+eta_name+"_"+pt_name;
      htemp_probejet_chHadEF_data = (TH1D*)f_mpf_data->Get(name_probejet_chHadEF_data);
      //      htemp_probejet_chHadEF_data->Print();
      int n_ev =  htemp_probejet_chHadEF_data->GetEntries();
      if(htemp_probejet_chHadEF_data->Integral() > 0)htemp_probejet_chHadEF_data->Scale(1/htemp_probejet_chHadEF_data->Integral());
      hEF->GetXaxis()->SetTitle("probejet chHadEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMaximum(0.5);
     hEF->SetMinimum(0.001);
      if(j<9) htemp_probejet_chHadEF_data->SetLineColor(j+1);
      else    htemp_probejet_chHadEF_data->SetLineColor(j-8);
      htemp_probejet_chHadEF_data->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_chHadEF_data->Draw("HIST SAME");
      leg14.AddEntry(htemp_probejet_chHadEF_data, legname);
    }

    leg14.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    c14->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_chHadEF_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + (abs_asymmetry_cut ? "_wAsymCut":"") +".pdf");


    TCanvas* c15 = new TCanvas();
    tdrCanvas(c15,"c15",hEF,4,10,kSquare,"MC");
    TLegend leg15 = tdrLeg(0.45,0.46,0.70,0.81);
    TH1D* htemp_probejet_photonEF_mc;
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_photonEF_mc = "hist_mc_probejet_photonEF_"+eta_name+"_"+pt_name;
      htemp_probejet_photonEF_mc = (TH1D*)f_mpf_mc->Get(name_probejet_photonEF_mc);
      //      htemp_probejet_photonEF_mc->Print();
      int n_ev =  htemp_probejet_photonEF_mc->GetEntries();
      if(htemp_probejet_photonEF_mc->Integral() > 0)htemp_probejet_photonEF_mc->Scale(1/htemp_probejet_photonEF_mc->Integral());
      hEF->GetXaxis()->SetTitle("probejet photonEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      // hEF->SetMaximum(3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
           hEF->SetMaximum(1.3);
      if(j<9) htemp_probejet_photonEF_mc->SetLineColor(j+1);
      else    htemp_probejet_photonEF_mc->SetLineColor(j-8);
      htemp_probejet_photonEF_mc->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_photonEF_mc->Draw("HIST SAME");
      leg15.AddEntry(htemp_probejet_photonEF_mc, legname,"l");
    }

    leg15.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    c15->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_photonEF_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + (abs_asymmetry_cut ? "_wAsymCut":"") +".pdf");

    TCanvas* c16 = new TCanvas();
    tdrCanvas(c16,"c16",hEF,4,10,kSquare,"DATA");
    TLegend leg16 = tdrLeg(0.45,0.46,0.70,0.81);
    TH1D* htemp_probejet_photonEF_data;
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_photonEF_data = "hist_data_probejet_photonEF_"+eta_name+"_"+pt_name;
      htemp_probejet_photonEF_data = (TH1D*)f_mpf_data->Get(name_probejet_photonEF_data);
      //      htemp_probejet_photonEF_data->Print();
      int n_ev =  htemp_probejet_photonEF_data->GetEntries();
      if(htemp_probejet_photonEF_data->Integral() > 0)htemp_probejet_photonEF_data->Scale(1/htemp_probejet_photonEF_data->Integral());
      hEF->GetXaxis()->SetTitle("probejet photonEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMaximum(1.3);
      if(j<9) htemp_probejet_photonEF_data->SetLineColor(j+1);
      else    htemp_probejet_photonEF_data->SetLineColor(j-8);
      htemp_probejet_photonEF_data->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_photonEF_data->Draw("HIST SAME");
      leg16.AddEntry(htemp_probejet_photonEF_data, legname);
    }

    leg16.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    c16->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_photonEF_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + ".pdf");


  TCanvas* c17 = new TCanvas();
    tdrCanvas(c17,"c17",hEF,4,10,kSquare,"MC");
    TLegend leg17 = tdrLeg(0.45,0.46,0.70,0.81);
    TH1D* htemp_probejet_muonEF_mc;
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_muonEF_mc = "hist_mc_probejet_muonEF_"+eta_name+"_"+pt_name;
      htemp_probejet_muonEF_mc = (TH1D*)f_mpf_mc->Get(name_probejet_muonEF_mc);
      //      htemp_probejet_muonEF_mc->Print();
      int n_ev =  htemp_probejet_muonEF_mc->GetEntries();
      if(htemp_probejet_muonEF_mc->Integral() > 0)htemp_probejet_muonEF_mc->Scale(1/htemp_probejet_muonEF_mc->Integral());
      hEF->GetXaxis()->SetTitle("probejet muonEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMaximum(1.3);
      if(j<9) htemp_probejet_muonEF_mc->SetLineColor(j+1);
      else    htemp_probejet_muonEF_mc->SetLineColor(j-8);
      htemp_probejet_muonEF_mc->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_muonEF_mc->Draw("HIST SAME");
      leg17.AddEntry(htemp_probejet_muonEF_mc, legname,"l");
    }

    leg17.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    c17->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_muonEF_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + (abs_asymmetry_cut ? "_wAsymCut":"") +".pdf");

    TCanvas* c18 = new TCanvas();
    tdrCanvas(c18,"c18",hEF,4,10,kSquare,"DATA");
    TLegend leg18 = tdrLeg(0.45,0.46,0.70,0.81);
    TH1D* htemp_probejet_muonEF_data;
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_muonEF_data = "hist_data_probejet_muonEF_"+eta_name+"_"+pt_name;
      htemp_probejet_muonEF_data = (TH1D*)f_mpf_data->Get(name_probejet_muonEF_data);
      htemp_probejet_muonEF_data->Print();
      int n_ev =  htemp_probejet_muonEF_data->GetEntries();
      if(htemp_probejet_muonEF_data->Integral() > 0)htemp_probejet_muonEF_data->Scale(1/htemp_probejet_muonEF_data->Integral());
      hEF->GetXaxis()->SetTitle("probejet muonEF");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      hEF->GetXaxis()->SetLimits(0,1.1);
      //hEF->GetYaxis()->SetLimits(0,0.1);
      // hEF->SetMaximum(3);
      hEF->SetMaximum(0.5);
      if(j<9) htemp_probejet_muonEF_data->SetLineColor(j+1);
      else    htemp_probejet_muonEF_data->SetLineColor(j-8);
      htemp_probejet_muonEF_data->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_muonEF_data->Draw("HIST SAME");
      leg18.AddEntry(htemp_probejet_muonEF_data, legname);
    }

    leg18.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    c18->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_muonEF_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + ".pdf");


    TCanvas* c19 = new TCanvas();
    tdrCanvas(c19,"c19",hEF,4,10,kSquare,"MC");
    TLegend leg19 = tdrLeg(0.45,0.46,0.70,0.81);
    TH1D* htemp_probejet_phi_mc;
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_phi_mc = "hist_mc_probejet_phi_"+eta_name+"_"+pt_name;
      htemp_probejet_phi_mc = (TH1D*)f_mpf_mc->Get(name_probejet_phi_mc);
      //      htemp_probejet_phi_mc->Print();
      int n_ev =  htemp_probejet_phi_mc->GetEntries();
      if(htemp_probejet_phi_mc->Integral() > 0)htemp_probejet_phi_mc->Scale(1/htemp_probejet_phi_mc->Integral());
      hEF->GetXaxis()->SetTitle("probejet phi");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      //      hEF->GetXaxis()->SetLimits(0,1.1);
      hEF->GetXaxis()->SetLimits(-3.15,3.15);
      hEF->GetYaxis()->SetLimits(0,0.1);
      //      hEF->SetMaximum(3);
      // hEF->SetMaximum(1.3);
      if(j<9) htemp_probejet_phi_mc->SetLineColor(j+1);
      else    htemp_probejet_phi_mc->SetLineColor(j-8);
      htemp_probejet_phi_mc->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_phi_mc->Draw("HIST SAME");
      leg19.AddEntry(htemp_probejet_phi_mc, legname,"l");
    }

    leg19.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    c19->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_phi_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + ".pdf");

    TCanvas* c20 = new TCanvas();
    tdrCanvas(c20,"c20",hEF,4,10,kSquare,"DATA");
    TLegend leg20 = tdrLeg(0.45,0.46,0.70,0.81);
    TH1D* htemp_probejet_phi_data;
    //    for(int j=0; j<n_pt-1; j++){
    for(int j=n_pt-6; j<n_pt-1; j++){//only last 3
    //    for(int j=0; j<5; j++){ //TEST
    //    for(int j=5; j<n_pt-1; j++){//TEST
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_probejet_phi_data = "hist_data_probejet_phi_"+eta_name+"_"+pt_name;
      htemp_probejet_phi_data = (TH1D*)f_mpf_data->Get(name_probejet_phi_data);
      //      htemp_probejet_phi_data->Print();
      int n_ev =  htemp_probejet_phi_data->GetEntries();
      if(htemp_probejet_phi_data->Integral() > 0)htemp_probejet_phi_data->Scale(1/htemp_probejet_phi_data->Integral());
      hEF->GetXaxis()->SetTitle("probejet phi");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      // hEF->GetXaxis()->SetLimits(0,1.1);
      hEF->GetXaxis()->SetLimits(-3.15,3.15);
      hEF->GetYaxis()->SetLimits(0,0.1);
      //      hEF->SetMaximum(3);
      //hEF->SetMaximum(1.3);
      if(j<9) htemp_probejet_phi_data->SetLineColor(j+1);
      else    htemp_probejet_phi_data->SetLineColor(j-8);
      htemp_probejet_phi_data->SetLineWidth(2+j*0.2);
      if(n_ev>100) htemp_probejet_phi_data->Draw("HIST SAME");
      leg20.AddEntry(htemp_probejet_phi_data, legname);
    }

    leg20.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    c20->SaveAs(CorrectionObject::_outpath+"plots/control/EnergyFractions/probejet_phi_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + ".pdf");

    
    //END Different energy fractions
    delete tex;
    // delete htemp_rel_data;
    // delete htemp_met_data;
    // delete c4;
    // delete htemp_rel_mc;
    // delete htemp_met_mc;
    // delete c3;
    // delete htemp_mpf_data;
    // delete c2;
    // delete htemp_mpf_mc;
    // delete c1;
}






}
