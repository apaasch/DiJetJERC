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

void CorrectionObject::FinalControlPlots_CorrectFormulae(double abs_asymmetry_cut, bool plot_full_A ,int lumi_bin){
  cout << "--------------- Starting FinalControlPlots_CorrectFormulae() ---------------" << endl << endl;
  gStyle->SetOptStat(0);
  TString txttag=CorrectionObject::_generator_tag; 
  int n_pt_ = max(n_pt,n_pt_HF);
  bool eta_cut_bool;
  int n_pt_cutted;
  
  if(plot_full_A) CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/fullAsym/");

  //Table with number of events in each pT- and eta-bin
  
  //Set up histos for ratios of responses
  TH1D *hdata_asymmetry[n_pt_-1][n_eta-1]; // A for data
  TH1D *hdata_B[n_pt_-1][n_eta-1];         // B for data
  TH1D *hdata_METoverJetsPt[n_pt_-1][n_eta-1];         // MET/sum_jets_pt for data
  TH1D *hdata_METoverSqrtJetsPt[n_pt_-1][n_eta-1];       //MET/Sqrt(sum_jets_pt)

  TH1D *hdata_probejet_neutEmEF[n_pt_-1][n_eta-1]; //neutral EM energy fraction
  TH1D *hdata_probejet_neutHadEF[n_pt_-1][n_eta-1]; //neutral hadron energy fraction
  TH1D *hdata_probejet_chEmEF[n_pt_-1][n_eta-1]; //charged EM energy fraction
  TH1D *hdata_probejet_chHadEF[n_pt_-1][n_eta-1]; //charged hadron energy fraction
  TH1D *hdata_probejet_photonEF[n_pt_-1][n_eta-1]; //photon energy fraction
  TH1D *hdata_probejet_muonEF[n_pt_-1][n_eta-1]; //muon hadron energy fraction
  TH1D *hdata_probejet_phi[n_pt_-1][n_eta-1]; //phi
  TH1D *hdata_MET[n_eta-1];//MET
  TH1D *hdata_alpha[n_eta-1];//alpha
  TH1D *hdata_jet3_pt[n_eta-1];//jet3_pt

  TH1D *hmc_asymmetry[n_pt_-1][n_eta-1];   // A for MC
  TH1D *hmc_B[n_pt_-1][n_eta-1];           // B for MC
  TH1D *hmc_METoverJetsPt[n_pt_-1][n_eta-1];         // MET/sum_jets_pt for MC
  TH1D *hmc_METoverSqrtJetsPt[n_pt_-1][n_eta-1];       //MET/Sqrt(sum_jets_pt)
  TH1F* hmc_pt_ave[n_eta-1];              // pt_ave for MC
  TH1F* hdata_pt_ave[n_eta-1];            // pt_ave for data

  TH1D *hmc_probejet_neutEmEF[n_pt_-1][n_eta-1]; //neutral EM energy fraction
  TH1D *hmc_probejet_neutHadEF[n_pt_-1][n_eta-1]; //neutral hadron energy fraction
  TH1D *hmc_probejet_chEmEF[n_pt_-1][n_eta-1]; //charged EM energy fraction
  TH1D *hmc_probejet_chHadEF[n_pt_-1][n_eta-1]; //charged hadron energy fraction
  TH1D *hmc_probejet_photonEF[n_pt_-1][n_eta-1]; //photon energy fraction
  TH1D *hmc_probejet_muonEF[n_pt_-1][n_eta-1]; //muon hadron energy fraction
  TH1D *hmc_probejet_phi[n_pt_-1][n_eta-1]; //phi
  TH1D *hmc_MET[n_eta-1];//MET
  TH1D *hmc_alpha[n_eta-1];//alpha
  TH1D *hmc_jet3_pt[n_eta-1];//jet3_pt
  

  int count = 0;
  TString name1 = "hist_data_A_";
  TString name2 = "hist_data_B_";
  TString name3 = "hist_mc_A_";
  TString name4 = "hist_mc_B_";
  TString name5 = "hist_data_METoverJetsPt_";
  TString name6 = "hist_mc_METoverJetsPt_";
  TString name7 = "hist_data_METoverSqrtJetsPt_";
  TString name8 = "hist_mc_METoverSqrtJetsPt_";

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
  TString name23 = "hist_mc_MET_";
  TString name24 = "hist_data_MET_";
 
  for(int j=0; j<n_eta-1; j++){
     eta_cut_bool = fabs(eta_bins[j])>eta_cut;
      TString eta_name = "eta_"+eta_range2[j]+"_"+eta_range2[j+1];
    for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];

      
      TString name = name1 + eta_name + "_" + pt_name; 
      hdata_asymmetry[k][j] = new TH1D(name,"",nResponseBins, Response_min, Response_max);
      name = name2 + eta_name + "_" + pt_name;
      hdata_B[k][j] = new TH1D(name,"",nResponseBins, Response_min, Response_max);
      name = name5 + eta_name + "_" + pt_name;
      hdata_METoverJetsPt[k][j] = new TH1D(name,"",60,0,1.2);
    
      name = name3 + eta_name + "_" + pt_name;
      hmc_asymmetry[k][j] = new TH1D(name,"",nResponseBins, Response_min, Response_max);
      name = name4 + eta_name + "_" + pt_name;
      hmc_B[k][j] = new TH1D(name,"",nResponseBins, Response_min, Response_max);
      name = name6 + eta_name + "_" + pt_name;
      hmc_METoverJetsPt[k][j] = new TH1D(name,"",50,0,1.2);

      name = name7 + eta_name + "_" + pt_name;
      hdata_METoverSqrtJetsPt[k][j] = new TH1D(name,"",50,0,5);
      name = name8 + eta_name + "_" + pt_name;
      hmc_METoverSqrtJetsPt[k][j] = new TH1D(name,"",50,0,5);

      name = name9  + eta_name + "_" + pt_name;
      hdata_probejet_neutEmEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name10  + eta_name + "_" + pt_name;
      hmc_probejet_neutEmEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name11 + eta_name + "_" + pt_name;

      hdata_probejet_neutHadEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name12  + eta_name + "_" + pt_name;
      hmc_probejet_neutHadEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name13 + eta_name + "_" + pt_name;

      hdata_probejet_chEmEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name14  + eta_name + "_" + pt_name;
      hmc_probejet_chEmEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name15  + eta_name + "_" + pt_name;

      hdata_probejet_chHadEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name16  + eta_name + "_" + pt_name;
      hmc_probejet_chHadEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name17  + eta_name + "_" + pt_name;

      hdata_probejet_photonEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name18  + eta_name + "_" + pt_name;
      hmc_probejet_photonEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name19  + eta_name + "_" + pt_name;

      hdata_probejet_muonEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name20  + eta_name + "_" + pt_name;
      hmc_probejet_muonEF [k][j] = new TH1D(name,"",60,0,1.5);
      name = name21  + eta_name + "_" + pt_name;

      hdata_probejet_phi [k][j] = new TH1D(name,"",60,-3.14,3.14);
      name = name22  + eta_name + "_" + pt_name;
      hmc_probejet_phi [k][j] = new TH1D(name,"",60,-3.14,3.14);

      count++;
    }

    TString name_jet3_pt = "hist_jet3_pt_";
    hdata_jet3_pt [j] = new TH1D(name_jet3_pt+"data"+eta_name,"",100,0,300);
    hmc_jet3_pt [j] = new TH1D(name_jet3_pt+"MC"+eta_name,"",100,0,300);

    TString name_alpha = "hist_alpha_";
    hdata_alpha [j] = new TH1D(name_alpha+"data"+eta_name,"",10,0,1);
    hmc_alpha [j] = new TH1D(name_alpha+"MC"+eta_name,"",10,0,1);


    TString name_MET = "hist_MET_";
    hdata_MET [j] = new TH1D(name_MET+"data"+eta_name,"",150,0,500);
    hmc_MET [j] = new TH1D(name_MET+"MC"+eta_name,"",150,0,500);

    //define pt_ave[eta] histos
    TString name_pt_ave = "hist_pt_ave_";
    hmc_pt_ave[j] = new TH1F(name_pt_ave+"MC_"+eta_name,"",1000,0,5000);
    hdata_pt_ave[j] = new TH1F(name_pt_ave+"data_"+eta_name,"",1000,0,5000);
  }


  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
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
  TTreeReaderValue<Float_t> probejet_phi_data(myReader_DATA, "probejet_phi");
  
  TTreeReaderValue<int> lumibin_data(myReader_DATA, "lumibin");
  int idx = 0;
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
    if(lumi_bin>=0 and *lumibin_data!=lumi_bin){
      myCount_cut++;
      continue;
    }
    //fill histos in bins of eta
    for(int i=0; i<n_eta-1; i++){
      if(fabs(*probejet_eta_data)<eta_bins[i+1] && fabs(*probejet_eta_data)>=eta_bins[i]){
	myCount++;
	hdata_pt_ave[i]->Fill(*pt_ave_data,*weight_data);
	hdata_MET[i]->Fill(*MET_data, *weight_data);
	hdata_alpha[i]->Fill(*alpha_data, *weight_data);
	hdata_jet3_pt[i]->Fill(*jet3_pt_data, *weight_data);
      }
    }

    //fill histos in bins of pt and eta
    for(int j=0; j<n_eta-1; j++){
    eta_cut_bool = fabs(eta_bins[j])>eta_cut; 
    for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      if(*pt_ave_data<(eta_cut_bool?pt_bins_HF:pt_bins)[k] || *pt_ave_data>(eta_cut_bool?pt_bins_HF:pt_bins)[k+1]) continue;
	if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
	  hdata_asymmetry[k][j]->Fill(*asymmetry_data,*weight_data);
	  hdata_B[k][j]->Fill(*B_data,*weight_data);
	  hdata_METoverJetsPt[k][j]->Fill((*MET_data)/(*sum_jets_pt_data+*probejet_pt_data+*barreljet_pt_data),*weight_data);
	  hdata_METoverSqrtJetsPt[k][j]->Fill((*MET_data)/(sqrt(*sum_jets_pt_data+*probejet_pt_data+*barreljet_pt_data)),*weight_data);
	  
	  hdata_probejet_neutEmEF[k][j]->Fill(*probejet_neutEmEF_data,*weight_data);
	  hdata_probejet_neutHadEF[k][j]->Fill(*probejet_neutHadEF_data,*weight_data);
	  hdata_probejet_chEmEF[k][j]->Fill(*probejet_chEmEF_data,*weight_data);
	  hdata_probejet_chHadEF[k][j]->Fill(*probejet_chHadEF_data,*weight_data);
	  hdata_probejet_photonEF[k][j]->Fill(*probejet_photonEF_data,*weight_data);
	  hdata_probejet_muonEF[k][j]->Fill(*probejet_muonEF_data,*weight_data);
	  hdata_probejet_phi[k][j]->Fill(*probejet_phi_data,*weight_data);
    }
    }
    idx++;
    if(idx%1000000==0) cout << "looping over data-TTree: Idx = " << idx << endl;
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
  idx=0;
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
    //fill histos in bins of eta
    
    for(int i=0; i<n_eta-1; i++){
      if(fabs(*probejet_eta_mc)<eta_bins[i+1] && fabs(*probejet_eta_mc)>=eta_bins[i]){
	myCount_mc+=*weight_mc;
	hmc_alpha[i]->Fill(*alpha_mc,*weight_mc);
	hmc_pt_ave[i]->Fill(*pt_ave_mc,*weight_mc);
	hmc_MET[i]->Fill(*MET_mc, *weight_mc);
	hmc_jet3_pt[i]->Fill(*jet3_pt_mc, *weight_mc);
      }
    }

    //fill histos in bins of pt and eta
    for(int j=0; j<n_eta-1; j++){
      eta_cut_bool = fabs(eta_bins[j])>eta_cut; 
    for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
	if(*pt_ave_mc<(eta_cut_bool?pt_bins_HF:pt_bins)[k] || *pt_ave_mc>(eta_cut_bool?pt_bins_HF:pt_bins)[k+1]) continue;
	if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
	else{
	  hmc_asymmetry[k][j]->Fill(*asymmetry_mc,*weight_mc);
	  hmc_B[k][j]->Fill(*B_mc,*weight_mc);
	  hmc_METoverJetsPt[k][j]->Fill((*MET_mc)/(*sum_jets_pt_mc+*probejet_pt_mc+*barreljet_pt_mc),*weight_mc);
	  hmc_METoverSqrtJetsPt[k][j]->Fill((*MET_mc)/(sqrt(*sum_jets_pt_mc+*probejet_pt_mc+*barreljet_pt_mc)),*weight_mc);
	  
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
    idx++;
    if(idx%1000000==0) cout << "looping over MC-TTree: Idx = " << idx << endl;
  }
  
  
  //DEBUG
  std::cout<<"\ncount mc "<<myCount_mc<<"  count cut mc "<<myCount_cut_mc<<std::endl<<std::endl;

  ofstream output;
  output.open(CorrectionObject::_outpath+"plots/control/Number_Events_Pt_Eta_bins_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")  + ".txt");

  output << "Number of events in each bin for MC" << endl;
  output << "|Eta|:          ";
  double n_tot_MC = 0;
  double n_tot_DATA = 0;
  for(int i=0; i<n_eta; i++) {
    if(i != n_eta-1) output << eta_range[i] << " -- ";
    else output << eta_range[i] << endl;
  }
  output << "central trg pt bins:\n";
  for(int i=0; i<n_pt-1; i++){
    if(i==0) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]  :    ";
    else if(i==1) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "] :    ";
    else output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]:    ";

    for(int j=0; j<n_eta-1; j++){
      if(fabs(eta_bins[j])>eta_cut){
	output << " ....   - ";
	continue;
      }
      if(j!=n_eta-2){
	if(hmc_B[i][j]->GetEntries()/1000 < 0.01)     output << hmc_B[i][j]->GetEntries() << "      - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 0.1) output << hmc_B[i][j]->GetEntries() << "     - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 1)   output << hmc_B[i][j]->GetEntries() << "    - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 <10)   output << hmc_B[i][j]->GetEntries() << "   - "; //<10000
	else if(hmc_B[i][j]->GetEntries()/1000 <100)  output << hmc_B[i][j]->GetEntries() << "  - ";
	else                                          output << hmc_B[i][j]->GetEntries() << " - ";
      }
      else output << hmc_B[i][j]->GetEntries();
      n_tot_MC+= hmc_B[i][j]->GetEntries();
    }
    output<<endl;

  }
  output << "HF trg pt bins:\n";
  for(int i=0; i<n_pt_HF-1; i++){
    if(i==0) output << "pT = ["  << fixed << setprecision(0) << pt_bins_HF[i] << "," << setprecision(0) << pt_bins_HF[i+1] << "]  :    ";
    else if(i==1) output << "pT = ["  << fixed << setprecision(0) << pt_bins_HF[i] << "," << setprecision(0) << pt_bins_HF[i+1] << "] :    ";
    else output << "pT = ["  << fixed << setprecision(0) << pt_bins_HF[i] << "," << setprecision(0) << pt_bins[i+1] << "]:    ";

    for(int j=0; j<n_eta-1; j++){
      if(fabs(eta_bins[j])<eta_cut){
	output << " ....   - ";
	continue;
      }
      if(j!=n_eta-2){
	if(hmc_B[i][j]->GetEntries()/1000 < 0.01)     output << hmc_B[i][j]->GetEntries() << "      - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 0.1) output << hmc_B[i][j]->GetEntries() << "     - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 1)   output << hmc_B[i][j]->GetEntries() << "    - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 <10)   output << hmc_B[i][j]->GetEntries() << "   - "; //<10000
	else if(hmc_B[i][j]->GetEntries()/1000 <100)  output << hmc_B[i][j]->GetEntries() << "  - ";
	else                                          output << hmc_B[i][j]->GetEntries() << " - ";
      }
      else output << hmc_B[i][j]->GetEntries();
      n_tot_MC+= hmc_B[i][j]->GetEntries();
    }
    output<<endl;

  }

  output << endl << endl << "Total number of events in MC: " << n_tot_MC << endl;

  output << endl << endl << endl << endl << "Number of events in each bin for DATA" << endl;
  output << "|Eta|:          ";
  for(int i=0; i<n_eta; i++) {
    if(i != n_eta-1) output << eta_range[i] << " -- ";
    else output << eta_range[i] << endl;
  }
  output << "central trg pt bins:\n";
  for(int i=0; i<n_pt-1; i++){
    output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]  :    ";
    for(int j=0; j<n_eta-1; j++){
      if(fabs(eta_bins[j])>eta_cut){
	output << " ....   - ";
	continue;
      }
      if(j!=n_eta-2){
	if(hdata_B[i][j]->GetEntries()/1000 < 0.01)     output << hdata_B[i][j]->GetEntries() << "      - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 < 0.1) output << hdata_B[i][j]->GetEntries() << "     - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 < 1)   output << hdata_B[i][j]->GetEntries() << "    - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 <10)   output << hdata_B[i][j]->GetEntries() << "   - "; //<10000
	else if(hdata_B[i][j]->GetEntries()/1000 <100)  output << hdata_B[i][j]->GetEntries() << "  - ";
	else                                                output << hdata_B[i][j]->GetEntries() << " - ";
      }
      else output << hdata_B[i][j]->GetEntries();
      n_tot_DATA += hdata_B[i][j]->GetEntries();
    }
    output<< endl;

  }
  output << "HF trg pt bins:\n";
  for(int i=0; i<n_pt_HF-1; i++){
    output << "pT = ["  << fixed << setprecision(0) << pt_bins_HF[i] << "," << setprecision(0) << pt_bins_HF[i+1] << "]  :    ";
    for(int j=0; j<n_eta-1; j++){
      if(fabs(eta_bins[j])<eta_cut){
	output << " ....   - ";
	continue;
      }
      if(j!=n_eta-2){
	if(hdata_B[i][j]->GetEntries()/1000 < 0.01)     output << hdata_B[i][j]->GetEntries() << "      - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 < 0.1) output << hdata_B[i][j]->GetEntries() << "     - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 < 1)   output << hdata_B[i][j]->GetEntries() << "    - "; //<1000
	else if(hdata_B[i][j]->GetEntries()/1000 <10)   output << hdata_B[i][j]->GetEntries() << "   - "; //<10000
	else if(hdata_B[i][j]->GetEntries()/1000 <100)  output << hdata_B[i][j]->GetEntries() << "  - ";
	else                                                output << hdata_B[i][j]->GetEntries() << " - ";
      }
      else output << hdata_B[i][j]->GetEntries();
      n_tot_DATA += hdata_B[i][j]->GetEntries();
    }
    output<< endl;

  }
  output << endl << endl << "Total number of events in DATA: " << n_tot_DATA << endl;



  // Dump 1-d distributions of A and B in bins of pT, eta

  TFile* test_out_mc_B = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_mc"+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+".root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    eta_cut_bool = fabs(eta_bins[j])>eta_cut; 
    for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      hmc_B[k][j]->Write();
      hmc_METoverJetsPt[k][j]->Write();
      hmc_METoverSqrtJetsPt[k][j]->Write();
      
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

  TFile* test_out_data_B = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_data"+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+".root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    eta_cut_bool = fabs(eta_bins[j])>eta_cut; 
    for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      hdata_B[k][j]->Write();
      hdata_METoverJetsPt[k][j]->Write();
      hdata_METoverSqrtJetsPt[k][j]->Write();
      
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

  TFile* test_out_mc_A = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_mc"+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+".root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    eta_cut_bool = fabs(eta_bins[j])>eta_cut; 
    for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      hmc_asymmetry[k][j]->Write();
    }
  }
  test_out_mc_A->Close();
  delete test_out_mc_A;

  TFile* test_out_data_A = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_data"+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+".root","RECREATE");
  for(int j=0; j<n_eta-1; j++){  
    eta_cut_bool = fabs(eta_bins[j])>eta_cut; 
    for(int k=0; k< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
      hdata_asymmetry[k][j]->Write();
    }
  }
  test_out_data_A->Close();
  delete test_out_data_A;

  //R_MC and R_DATA overlaid in the same plot as a function of pT, in bins of |eta|
  double val_rel_mc[n_eta-1][n_pt_-1]; //value at pt,eta
  double err_rel_mc[n_eta-1][n_pt_-1]; //error of ratio at pt,eta
  double val_mpf_mc[n_eta-1][n_pt_-1]; //ratio at pt,eta
  double err_mpf_mc[n_eta-1][n_pt_-1]; //error of ratio at pt,eta
  double val_rel_data[n_eta-1][n_pt_-1]; //value at pt,eta
  double err_rel_data[n_eta-1][n_pt_-1]; //error of ratio at pt,eta
  double val_mpf_data[n_eta-1][n_pt_-1]; //ratio at pt,eta
  double err_mpf_data[n_eta-1][n_pt_-1]; //error of ratio at pt,eta

  for(int i=0; i<n_eta-1; i++){
    eta_cut_bool = fabs(eta_bins[i])>eta_cut; 
    for(int j=0; j< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){

      //get <A> and error on <A>
      pair <double,double> A_mc = GetValueAndError(hmc_asymmetry[j][i]);
      pair <double,double> A_data = GetValueAndError(hdata_asymmetry[j][i]);
      pair <double,double> B_mc = GetValueAndError(hmc_B[j][i]);
      pair <double,double> B_data = GetValueAndError(hdata_B[j][i]);

      //build MPF and pt_bal and their errors
      pair<double,double> res_mc_rel_r,res_data_rel_r;
      pair<double,double> res_mc_mpf_r,res_data_mpf_r;
      res_mc_mpf_r.first = (1+B_mc.first)/(1-B_mc.first);
      res_mc_mpf_r.second = 2/(pow((1-B_mc.first),2)) * B_mc.second;
      res_data_mpf_r.first = (1+B_data.first)/(1-B_data.first);
      res_data_mpf_r.second = 2/(pow((1-B_data.first),2)) * B_data.second;
      res_mc_rel_r.first = (1+A_mc.first)/(1-A_mc.first);
      res_mc_rel_r.second = 2/(pow((1-A_mc.first),2)) * A_mc.second;
      res_data_rel_r.first = (1+A_data.first)/(1-A_data.first);
      res_data_rel_r.second = 2/(pow((1-A_data.first),2)) * A_data.second;

      val_rel_mc[i][j] = res_mc_rel_r.first;
      cout<<"---- HEY! eta="<<eta_bins[i]<<" rel_r="<<res_mc_rel_r.first<<endl;
      err_rel_mc[i][j] = res_mc_rel_r.second;
      val_mpf_mc[i][j] = res_mc_mpf_r.first;
      err_mpf_mc[i][j] = res_mc_mpf_r.second;
      val_rel_data[i][j] = res_data_rel_r.first;
      err_rel_data[i][j] = res_data_rel_r.second;
      val_mpf_data[i][j] = res_data_mpf_r.first;
      err_mpf_data[i][j] = res_data_mpf_r.second;
    }
  }

  //dummy for tdrCanvas
  TH1D *h = new TH1D("h",";dummy;",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);

  TH1D *hEF = new TH1D("hEF",";dummy;",1000,0,5.191);

  TCanvas* c_0 = new TCanvas();
  tdrCanvas(c_0,"c_0",h,4,10,true,CorrectionObject::_lumitag);

  for(int i=0; i<n_eta-1; i++){
    eta_cut_bool = fabs(eta_bins[i])>eta_cut;
    n_pt_cutted = ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 );
    //Create and fill TGraphErrors
    double xbin_tgraph[n_pt_cutted];
    double zero[n_pt_cutted];
    for(int i=0;i<n_pt_cutted ;i++){
      xbin_tgraph[i]=((eta_cut_bool?pt_bins_HF:pt_bins)[i]+(eta_cut_bool?pt_bins_HF:pt_bins)[i+1])/2;
      zero[i]=((eta_cut_bool?pt_bins_HF:pt_bins)[i+1]-(eta_cut_bool?pt_bins_HF:pt_bins)[i])/2 ;
    }

    TGraphErrors *graph_mpf_mc   = new TGraphErrors(n_pt_cutted, xbin_tgraph, val_mpf_mc[i], zero, err_mpf_mc[i]);
    TGraphErrors *graph_mpf_data = new TGraphErrors(n_pt_cutted, xbin_tgraph, val_mpf_data[i], zero, err_mpf_data[i]);
    TGraphErrors *graph_rel_mc   = new TGraphErrors(n_pt_cutted, xbin_tgraph, val_rel_mc[i], zero, err_rel_mc[i]);
    TGraphErrors *graph_rel_data = new TGraphErrors(n_pt_cutted, xbin_tgraph, val_rel_data[i], zero, err_rel_data[i]);
    graph_mpf_mc   = (TGraphErrors*)CleanEmptyPoints(graph_mpf_mc);
    graph_mpf_data = (TGraphErrors*)CleanEmptyPoints(graph_mpf_data);
    graph_rel_mc   = (TGraphErrors*)CleanEmptyPoints(graph_rel_mc);
    graph_rel_data = (TGraphErrors*)CleanEmptyPoints(graph_rel_data);

    graph_mpf_mc->SetTitle("");
    graph_mpf_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_mpf_mc->GetXaxis()->SetTitleSize(0.05);
    graph_mpf_mc->GetXaxis()->SetTitleOffset(0.80);
    graph_mpf_mc->GetXaxis()->SetLimits(30,(eta_cut_bool?pt_bins_HF:pt_bins)[n_pt_cutted]+100);
    graph_mpf_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_mpf_mc->SetMarkerColor(kRed);
    graph_mpf_mc->SetMarkerStyle(20);
    graph_mpf_mc->SetLineColor(kRed);

    graph_mpf_data->SetTitle("");
    graph_mpf_data->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_mpf_data->GetXaxis()->SetTitleSize(0.05);
    graph_mpf_data->GetXaxis()->SetTitleOffset(0.80);
    graph_mpf_data->GetXaxis()->SetLimits(30,(eta_cut_bool?pt_bins_HF:pt_bins)[n_pt_cutted]+100);
    graph_mpf_data->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_mpf_data->SetMarkerColor(kBlack);
    graph_mpf_data->SetMarkerStyle(20);
    graph_mpf_data->SetLineColor(kBlack);

    graph_rel_mc->SetTitle("");
    graph_rel_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_mc->GetXaxis()->SetTitleSize(0.05);
    graph_rel_mc->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_mc->GetXaxis()->SetLimits(30,(eta_cut_bool?pt_bins_HF:pt_bins)[n_pt_cutted]+100);
    graph_rel_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_mc->SetMarkerColor(kRed);
    graph_rel_mc->SetMarkerStyle(20);
    graph_rel_mc->SetLineColor(kRed);

    graph_rel_data->SetTitle("");
    graph_rel_data->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_data->GetXaxis()->SetTitleSize(0.05);
    graph_rel_data->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_data->GetXaxis()->SetLimits(30,(eta_cut_bool?pt_bins_HF:pt_bins)[n_pt_cutted]+100);
    graph_rel_data->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_data->SetMarkerColor(kBlack);
    graph_rel_data->SetMarkerStyle(20);
    graph_rel_data->SetLineColor(kBlack);

    TString alVal;
    alVal.Form("%0.2f\n",alpha_cut);
    TString altitle = "{#alpha<"+alVal+"}";
    TString axistitle_mc = "R^{MC}_";
    TString axistitle_data = "R^{DATA}_";
 
    axistitle_mc   += altitle;
    axistitle_data += altitle;

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 

    TCanvas* c_mpf = new TCanvas();
    tdrCanvas(c_mpf,"c_mpf",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(30,(eta_cut_bool?pt_bins_HF:pt_bins)[n_pt_cutted]+100);
    h->GetYaxis()->SetRangeUser(0.50,1.50);
    graph_mpf_mc->Draw("P SAME");
    graph_mpf_data->Draw("P SAME");
    gPad->SetLogx();

    TLegend *leg_mpf;
    leg_mpf = new TLegend(0.35,0.72,0.51,0.92,"","brNDC");//x+0.1
    leg_mpf->SetBorderSize(0);
    leg_mpf->SetTextSize(0.036);
    leg_mpf->SetFillColor(10);
    leg_mpf->SetFillStyle(0);
    leg_mpf->SetLineColor(1);
    leg_mpf->SetTextFont(42);
    leg_mpf->SetHeader("MPF response, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+", #alpha<"+s_alpha_cut);
    leg_mpf->AddEntry(graph_mpf_mc, axistitle_mc,"P");
    leg_mpf->AddEntry(graph_mpf_data, axistitle_data,"P");
    leg_mpf->Draw();

    c_mpf->SaveAs(CorrectionObject::_outpath+"plots/control/MPF_Response_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");

    TCanvas* c_rel = new TCanvas();
    tdrCanvas(c_rel,"c_rel",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(30,(eta_cut_bool?pt_bins_HF:pt_bins)[n_pt_cutted]+100);
    h->GetYaxis()->SetRangeUser(0.50,1.50);
    graph_rel_mc->Draw("P SAME");
    graph_rel_data->Draw("P SAME");
    gPad->SetLogx();

    TLegend *leg_rel;
    leg_rel = new TLegend(0.35,0.72,0.51,0.92,"","brNDC");//x+0.1
    leg_rel->SetBorderSize(0);
    leg_rel->SetTextSize(0.036);
    leg_rel->SetFillColor(10);
    leg_rel->SetFillStyle(0);
    leg_rel->SetLineColor(1);
    leg_rel->SetTextFont(42);
    leg_rel->SetHeader("p_{T}-balance response, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]); 
    leg_rel->AddEntry(graph_rel_mc, axistitle_mc,"P");
    leg_rel->AddEntry(graph_rel_data, axistitle_data,"P");
    leg_rel->Draw();

    c_rel->SaveAs(CorrectionObject::_outpath+"plots/control/Rel_Response_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");

    delete c_rel;
    delete c_mpf;
    delete tex;
    delete graph_rel_data;
    delete graph_rel_mc;
    delete graph_mpf_data;
    delete graph_mpf_mc;
  }



//*******************************************  Plot all Control Hists **************************

  //Plot 1d response distributions in a particular eta-bin for different pt-bins onto a single canvas

  //Get histo files
  TFile* f_mpf_mc = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_mc"+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+".root","READ");
  TFile* f_mpf_data = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_data"+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+".root","READ");
  TFile* f_rel_mc = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_mc"+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+".root","READ");
  TFile* f_rel_data = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_data"+ (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+".root","READ");
  for(int i=0; i<n_eta-1; i++){
    TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 
    

    TCanvas* c1 = new TCanvas();
    tdrCanvas(c1,"c1",h,4,10,kSquare,"MC");
    TLegend leg1 = tdrLeg(0.17,0.6,0.85,0.81);
    leg1.SetNColumns(2);
    TH1D* htemp_mpf_mc;
    gPad->SetLogy();

    eta_cut_bool = fabs(eta_bins[i])>eta_cut; 
    for(int j=0; j < ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
      TString name_mpf_mc = "hist_mc_B_"+eta_name+"_"+pt_name;
      htemp_mpf_mc = (TH1D*)f_mpf_mc->Get(name_mpf_mc);
      int n_ev = htemp_mpf_mc->GetEntries();
      if(htemp_mpf_mc->Integral() > 0)htemp_mpf_mc->Scale(1/htemp_mpf_mc->Integral());
      h->GetXaxis()->SetTitle("B");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-1.2,1.2);
      h->SetMinimum(0.001);
      h->SetMaximum(3);
      if(j<9) htemp_mpf_mc->SetLineColor(j+1);
      else    htemp_mpf_mc->SetLineColor(j+31);
      htemp_mpf_mc->SetLineWidth(3);
      if(n_ev>100) htemp_mpf_mc->Draw("HIST SAME");
      leg1.AddEntry(htemp_mpf_mc, legname, "l");
    }

    leg1.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    c1->SaveAs(CorrectionObject::_outpath+"plots/control/B_NormDistribution_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");


    TCanvas* c2 = new TCanvas();
    tdrCanvas(c2,"c2",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg2 = tdrLeg(0.17,0.6,0.85,0.79);
    leg2.SetNColumns(2);

    TH1D* htemp_mpf_data;
    gPad->SetLogy();
    for(int j=0; j<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){
   
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
      TString name_mpf_data = "hist_data_B_"+eta_name+"_"+pt_name;
      htemp_mpf_data = (TH1D*)f_mpf_data->Get(name_mpf_data);
      int n_ev = htemp_mpf_data->GetEntries();
      if(htemp_mpf_data->Integral() > 0)htemp_mpf_data->Scale(1/htemp_mpf_data->Integral());
      h->GetXaxis()->SetTitle("B");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-1.2,1.2);
      h->SetMinimum(0.001);
      h->SetMaximum(3);
      if(j<9) htemp_mpf_data->SetLineColor(j+1);
      else    htemp_mpf_data->SetLineColor(j+31);
      htemp_mpf_data->SetLineWidth(3);
      if(n_ev>100) htemp_mpf_data->Draw("HIST SAME");
      leg2.AddEntry(htemp_mpf_data, legname ,"l");
    }

    leg2.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    //tex_lumi->DrawLatex(0.50,0.91,CorrectionObject::_lumitag+"(13TeV)");
    c2->SaveAs(CorrectionObject::_outpath+"plots/control/B_NormDistribution_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");


    TCanvas* c3 = new TCanvas();
    tdrCanvas(c3,"c3",h,4,10,kSquare,"MC");
    TLegend leg3 = tdrLeg(0.17,0.6,0.85,0.79);
    leg3.SetNColumns(2);
    TH1D* htemp_rel_mc;
    gPad->SetLogy();

    for(int j=0; j<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
      TString name_rel_mc = "hist_mc_A_"+eta_name+"_"+pt_name;
      htemp_rel_mc = (TH1D*)f_rel_mc->Get(name_rel_mc);
      int n_ev =  htemp_rel_mc->GetEntries();
      if(htemp_rel_mc->Integral() > 0)htemp_rel_mc->Scale(1/htemp_rel_mc->Integral());
      h->GetXaxis()->SetTitle("A");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-1.2,1.2);
      h->SetMinimum(0.001);
      h->SetMaximum(3);
      if(j<9) htemp_rel_mc->SetLineColor(j+1);
      else    htemp_rel_mc->SetLineColor(j+31);
      htemp_rel_mc->SetLineWidth(3);
      if(n_ev>100) htemp_rel_mc->Draw("HIST SAME");
      leg3.AddEntry(htemp_rel_mc, legname);
    }

    leg3.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    //tex_lumi->DrawLatex(0.6,0.91,"MC");
    c3->SaveAs(CorrectionObject::_outpath+"plots/control/A_NormDistribution_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")  + ".pdf");



    if(plot_full_A){
      for(int j=0; j<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){
	        TCanvas* cFullA = new TCanvas();
		tdrCanvas(cFullA,"cFullA",h,4,10,kSquare,CorrectionObject::_lumitag);
		TH1D* htemp_rel_data;
		TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
		TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
		TString name_rel_data = "hist_data_A_"+eta_name+"_"+pt_name;
		htemp_rel_data = (TH1D*)f_rel_data->Get(name_rel_data);
		htemp_rel_data->Draw("E");
		htemp_rel_data->GetXaxis()->SetTitle("A");
		htemp_rel_data->GetYaxis()->SetTitle("Entries per Bin");
		htemp_rel_data->GetYaxis()->SetTitleOffset(1.5);
		htemp_rel_data->GetXaxis()->SetLimits(-1.2,1.2);
		htemp_rel_data->Draw("E");		
		tex->DrawLatex(0.47,0.85,"Data, " + text);
		tex->DrawLatex(0.54,0.8,legname);		
		cFullA->SaveAs(CorrectionObject::_outpath+"plots/control/fullAsym/A_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1]+"_pt_"+ (eta_cut_bool?pt_range_HF:pt_range)[j] + "_" + (eta_cut_bool?pt_range_HF:pt_range)[j+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");
		delete cFullA;
		delete htemp_rel_data;
	  }
    }
    
    TCanvas* c4 = new TCanvas();
    tdrCanvas(c4,"c4",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg4 = tdrLeg(0.17,0.6,0.85,0.79);
    leg4.SetNColumns(2);
    TH1D* htemp_rel_data;

   gPad->SetLogy();
    for(int j=0; j<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
      TString name_rel_data = "hist_data_A_"+eta_name+"_"+pt_name;
      htemp_rel_data = (TH1D*)f_rel_data->Get(name_rel_data);
      int n_ev = htemp_rel_data->GetEntries();
      if(htemp_rel_data->Integral() > 0)htemp_rel_data->Scale(1/htemp_rel_data->Integral());
      h->GetXaxis()->SetTitle("A");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-1.2,1.2);
      h->SetMinimum(0.001);
      h->SetMaximum(3);
      if(j<9) htemp_rel_data->SetLineColor(j+1);
      else    htemp_rel_data->SetLineColor(j+31);
      htemp_rel_data->SetLineWidth(3);
      if(n_ev>100) htemp_rel_data->Draw("HIST SAME");
      leg4.AddEntry(htemp_rel_data, legname);
    }

    leg4.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    //tex_lumi->DrawLatex(0.50,0.91,CorrectionObject::_lumitag+"(13TeV)");
    c4->SaveAs(CorrectionObject::_outpath+"plots/control/A_NormDistribution_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");    



    ///MET over sum pt
    TCanvas* c5 = new TCanvas();
    tdrCanvas(c5,"c5",h,4,10,kSquare,"MC");
    TLegend leg5 = tdrLeg(0.22,0.6,0.88,0.79);
    leg5.SetNColumns(2);
    TH1D* htemp_met_mc;
    for(int j=0; j<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
      TString name_met_mc = "hist_mc_METoverJetsPt_"+eta_name+"_"+pt_name;
      htemp_met_mc = (TH1D*)f_mpf_mc->Get(name_met_mc);
    
      int n_ev =  htemp_met_mc->GetEntries();
      if(htemp_met_mc->Integral() > 0)htemp_met_mc->Scale(1/htemp_met_mc->Integral());
      h->GetXaxis()->SetTitle("MET/#sum p_{T}");
      h->GetYaxis()->SetTitle("Norm. Entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      h->GetXaxis()->SetLimits(0,1.2);
      //      h->GetYaxis()->SetLimits(0,0.8);
      h->SetMaximum(0.3);
      if(j<9) htemp_met_mc->SetLineColor(j+1);
      else    htemp_met_mc->SetLineColor(j+31);
      htemp_met_mc->SetLineWidth(3);
      if(n_ev>100) htemp_met_mc->Draw("HIST SAME");
      leg5.AddEntry(htemp_met_mc, legname);
    }

    leg5.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    //tex_lumi->DrawLatex(0.6,0.91,"MC");
    c5->SaveAs(CorrectionObject::_outpath+"plots/control/METoverPt_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");
 


   TCanvas* c6 = new TCanvas();
    tdrCanvas(c6,"c6",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg6 = tdrLeg(0.22,0.6,0.88,0.79);
    leg6.SetNColumns(2);
    TH1D* htemp_met_data;
    for(int j=0; j<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
      TString name_met_data = "hist_data_METoverJetsPt_"+eta_name+"_"+pt_name;
      htemp_met_data = (TH1D*)f_mpf_data->Get(name_met_data);
      int n_ev = htemp_met_data->GetEntries();
      if(htemp_met_data->Integral() > 0)htemp_met_data->Scale(1/htemp_met_data->Integral());
      h->GetXaxis()->SetTitle("MET/#sum p_{T}");
      h->GetYaxis()->SetTitle("Norm. Entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,1.2);
      h->GetYaxis()->SetLimits(0,0.8);
      h->SetMaximum(0.3);
      if(j<9) htemp_met_data->SetLineColor(j+1);
      else    htemp_met_data->SetLineColor(j+31);
      htemp_met_data->SetLineWidth(3);
      if(n_ev>100) htemp_met_data->Draw("HIST SAME");
      leg6.AddEntry(htemp_met_data, legname);
    }
    leg6.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);

    //tex_lumi->DrawLatex(0.50,0.91,CorrectionObject::_lumitag+"(13TeV)");
    c6->SaveAs(CorrectionObject::_outpath+"plots/control/METoverPt_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + (abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") +".pdf");    

    ///END MET over sum pt

    TCanvas* c19 = new TCanvas();
    tdrCanvas(c19,"c19",hEF,4,10,kSquare,"MC");
    TLegend leg19 = tdrLeg(0.45,0.46,0.70,0.81);
    TH1D* htemp_probejet_phi_mc;
    for(int j=0; j<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){ 
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
      TString name_probejet_phi_mc = "hist_mc_probejet_phi_"+eta_name+"_"+pt_name;
      htemp_probejet_phi_mc = (TH1D*)f_mpf_mc->Get(name_probejet_phi_mc);
      //      htemp_probejet_phi_mc->Print();
      int n_ev =  htemp_probejet_phi_mc->GetEntries();
      if(htemp_probejet_phi_mc->Integral() > 0)htemp_probejet_phi_mc->Scale(1/htemp_probejet_phi_mc->Integral());
      hEF->GetXaxis()->SetTitle("probejet phi");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      //      hEF->GetXaxis()->SetLimits(0,1.5);
      hEF->GetXaxis()->SetLimits(-3.15,3.15);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      hEF->SetMaximum(0.1);
      //      hEF->SetMaximum(0.8);
      if(j<9) htemp_probejet_phi_mc->SetLineColor(j+1);
      else    htemp_probejet_phi_mc->SetLineColor(j+31);
      htemp_probejet_phi_mc->SetLineWidth(3);
      if(n_ev>100) htemp_probejet_phi_mc->Draw("HIST SAME");
      leg19.AddEntry(htemp_probejet_phi_mc, legname,"l");
    }

    leg19.Draw();
    tex->DrawLatex(0.47,0.85,"MC, " + text);
    tex->DrawLatex(0.15,0.95,txttag);
    c19->SaveAs(CorrectionObject::_outpath+"plots/control/probejet_phi_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")+ (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");

    TCanvas* c20 = new TCanvas();
    tdrCanvas(c20,"c20",hEF,4,10,kSquare,"DATA");
    TLegend leg20 = tdrLeg(0.45,0.46,0.70,0.81);
    TH1D* htemp_probejet_phi_data;
    for(int j=0; j<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); j++){ 
      TString pt_name = "pt_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      TString legname = "p_{T} #in [" + (eta_cut_bool?pt_range_HF:pt_range)[j] + "," + (eta_cut_bool?pt_range_HF:pt_range)[j+1] + "]";
      TString name_probejet_phi_data = "hist_data_probejet_phi_"+eta_name+"_"+pt_name;
      htemp_probejet_phi_data = (TH1D*)f_mpf_data->Get(name_probejet_phi_data);
      //      htemp_probejet_phi_data->Print();
      int n_ev =  htemp_probejet_phi_data->GetEntries();
      if(htemp_probejet_phi_data->Integral() > 0)htemp_probejet_phi_data->Scale(1/htemp_probejet_phi_data->Integral());
      hEF->GetXaxis()->SetTitle("probejet phi");
      hEF->GetYaxis()->SetTitle("Norm. Entries");
      hEF->GetYaxis()->SetTitleOffset(1.5);
      // h->SetMaximum(0.3);
      // hEF->GetXaxis()->SetLimits(0,1.5);
      hEF->GetXaxis()->SetLimits(-3.15,3.15);
      //      hEF->GetYaxis()->SetLimits(0,0.1);
      hEF->SetMaximum(0.1);
      //hEF->SetMaximum(0.8);
      if(j<9) htemp_probejet_phi_data->SetLineColor(j+1);
      else    htemp_probejet_phi_data->SetLineColor(j+31);
      htemp_probejet_phi_data->SetLineWidth(3);
      if(n_ev>100) htemp_probejet_phi_data->Draw("HIST SAME");
      leg20.AddEntry(htemp_probejet_phi_data, legname);
    }

    leg20.Draw();
    tex->DrawLatex(0.47,0.85,"Data, " + text);
    c20->SaveAs(CorrectionObject::_outpath+"plots/control/probejet_phi_DATA_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")+ (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");
   
    //END Different energy fractions

    delete tex;
    delete htemp_rel_data;
    delete htemp_met_data;
    delete c4;
    delete htemp_rel_mc;
    delete htemp_met_mc;
    delete c3;
    delete htemp_mpf_data;
    delete c2;
    delete htemp_mpf_mc;
    delete c1;
  } //END big n_eta loop



  //pT_ave for MC and data in bins of |eta|

  for(int i=0; i<n_eta-1; i++){
    TCanvas* c1 = new TCanvas();
    tdrCanvas(c1,"c1",h,4,10,kSquare,CorrectionObject::_lumitag);
  
    TLegend leg1 = tdrLeg(0.62,0.66,0.85,0.81);
    h->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
    h->GetXaxis()->SetLimits(0,2000);
    h->GetYaxis()->SetTitle("entries");
    h->GetYaxis()->SetTitleOffset(1.5);
    double maximum = std::max(hdata_pt_ave[i]->GetMaximum(), hmc_pt_ave[i]->GetMaximum());
    h->GetYaxis()->SetRangeUser(0,1.2*maximum);
    hdata_pt_ave[i]->SetMarkerColor(kBlack);
    hdata_pt_ave[i]->SetMarkerStyle(20);
    hdata_pt_ave[i]->Draw("SAME P");
    hmc_pt_ave[i]->SetLineColor(kBlue);
    hmc_pt_ave[i]->Draw("HIST SAME");
    leg1.AddEntry(hdata_pt_ave[i], "DATA");
    leg1.AddEntry(hmc_pt_ave[i], "MC");
    leg1.Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];
    tex->DrawLatex(0.52,0.85, text);

    c1->SaveAs(CorrectionObject::_outpath+"plots/control/Pt_ave_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")+ (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");    


    delete tex;
    delete c1;


    TCanvas* c2 = new TCanvas();
    tdrCanvas(c2,"c2",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg2 = tdrLeg(0.62,0.66,0.85,0.81);
    h->GetXaxis()->SetTitle("MET [GeV]");
    h->GetXaxis()->SetLimits(0,500);
    h->GetYaxis()->SetTitle("entries");
    h->GetYaxis()->SetTitleOffset(1.5);
    double maximum2 = std::max(hdata_MET[i]->GetMaximum(), hmc_MET[i]->GetMaximum());
    h->GetYaxis()->SetRangeUser(0,1.2*maximum2);
    hdata_MET[i]->SetMarkerColor(kBlack);
    hdata_MET[i]->SetMarkerStyle(20);
    hdata_MET[i]->Draw("SAME P");
    hmc_MET[i]->SetLineColor(kRed);
    hmc_MET[i]->Draw("L SAME");
    leg2.AddEntry(hdata_MET[i], "DATA");
    leg2.AddEntry(hmc_MET[i], "MC");
    leg2.Draw();

    TLatex *tex1 = new TLatex();
    tex1->SetNDC();
    tex1->SetTextSize(0.045); 
    TString text1 = eta_range[i] + " < |#eta| < " + eta_range[i+1];
    tex1->DrawLatex(0.52,0.85, text1);

    c2->SaveAs(CorrectionObject::_outpath+"plots/control/MET_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + (abs_asymmetry_cut ? "_wAsymCut":"")+ (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") +".pdf");    


    delete tex1;
    delete c2;

 TCanvas* c3 = new TCanvas();
    tdrCanvas(c3,"c3",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg3 = tdrLeg(0.62,0.66,0.85,0.81);
    h->SetMaximum(1);
    h->SetMinimum(0);
    h->GetXaxis()->SetTitle("alpha");
    h->GetXaxis()->SetLimits(0,1);
    h->GetYaxis()->SetTitle("entries");
    h->GetYaxis()->SetTitleOffset(1.5);
    double maximum3 = std::max(hdata_alpha[i]->GetMaximum(), hmc_alpha[i]->GetMaximum());
    h->GetYaxis()->SetRangeUser(0,1.2*maximum3);
    hdata_alpha[i]->SetMarkerColor(kBlack);
    hdata_alpha[i]->SetMarkerStyle(20);
    hdata_alpha[i]->DrawNormalized("SAME P");
    hmc_alpha[i]->SetLineColor(kBlue);
    hmc_alpha[i]->DrawNormalized("HIST SAME");
    leg3.AddEntry(hdata_alpha[i], "DATA");
    leg3.AddEntry(hmc_alpha[i], "MC");
    leg3.Draw();

    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.045); 
    TString text2 = eta_range[i] + " < |#eta| < " + eta_range[i+1];
    tex2->DrawLatex(0.52,0.85, text1);

    c3->SaveAs(CorrectionObject::_outpath+"plots/control/alpha_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"")+ (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "") + ".pdf");    


    delete tex2;
    delete c3;

   
TCanvas* c4 = new TCanvas();
    tdrCanvas(c4,"c4",h,4,10,kSquare,CorrectionObject::_lumitag);
    TLegend leg4 = tdrLeg(0.62,0.66,0.85,0.81);
    h->GetXaxis()->SetTitle("jet3 pt");
    h->GetXaxis()->SetLimits(0,300);
    h->GetYaxis()->SetTitle("entries");
    h->GetYaxis()->SetTitleOffset(1.5);
    double maximum4 = std::max(hdata_jet3_pt[i]->GetMaximum(), hmc_jet3_pt[i]->GetMaximum());
    h->GetYaxis()->SetRangeUser(0,1.2*maximum4);
    hdata_jet3_pt[i]->SetMarkerColor(kBlack);
    hdata_jet3_pt[i]->SetMarkerStyle(20);
    hdata_jet3_pt[i]->Draw("SAME P");
    hmc_jet3_pt[i]->SetLineColor(kBlue);
    hmc_jet3_pt[i]->Draw("HIST SAME");
    leg4.AddEntry(hdata_jet3_pt[i], "DATA");
    leg4.AddEntry(hmc_jet3_pt[i], "MC");
    leg4.Draw();

    TLatex *tex4 = new TLatex();
    tex4->SetNDC();
    tex4->SetTextSize(0.045); 
    TString text4 = eta_range[i] + " < |#eta| < " + eta_range[i+1];
    tex4->DrawLatex(0.52,0.85, text1);

    c4->SaveAs(CorrectionObject::_outpath+"plots/control/jet3_pt_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +(abs_asymmetry_cut ? "_wAsymCut":"") + (lumi_bin>=0 ? "_lumiBin" + to_string(lumi_bin)  : "")+ ".pdf");    

    delete tex4;
    delete c4;
  }

  delete c_0;
  delete h;
  
  for(int j=0; j<n_eta-1; j++){
    eta_cut_bool = fabs(eta_bins[j])>eta_cut; 
    for(int i=0; i< ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); i++){
      delete hdata_asymmetry[i][j];
      delete hmc_asymmetry[i][j];
      delete hdata_B[i][j];
      delete hmc_B[i][j];
    }
  }

  for(int i=0; i<n_eta-1; i++){
    delete hdata_pt_ave[i];
    delete hmc_pt_ave[i];
  }

  delete f_mpf_mc;
  delete f_mpf_data;
  delete f_rel_mc;
  delete f_rel_data;

}
