#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "TROOT.h"
#include "TMath.h"
#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/DiJetJERC/include/constants.h"
#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/PersonalCode/tdrstyle_all.C"


int shift_SM = 11;
int shift_FE = 3;
int shift_barrel = 1;

void LoadSF(std::vector<std::vector<std::vector<double>>> &SFs, TString filename) {
  std::vector<std::vector<double> > SF;
  std::vector<double> eta_bin_SM_center; 						//0
  std::vector<double> eta_bin_SM_error;    					//1
  std::vector<double> SF_uncorrelated_SM;       		//2
  std::vector<double> SF_uncorrelated_SM_error;			//3
  std::vector<double> SF_correlated_SM;        		 	//4
  std::vector<double> SF_correlated_SM_error;   		//5
  std::vector<double> SF_uncorrelated_SM_ptdep_min; //6
  std::vector<double> SF_uncorrelated_SM_ptdep_max;	//7
  std::vector<double> SF_correlated_SM_ptdep_min;   //8
  std::vector<double> SF_correlated_SM_ptdep_max;   //9

  // SF_uncorrelated_FE, SF_uncorrelated_FE_error, SF_correlated_FE, SF_correlated_FE_error, eta_bin_FE_center, eta_bin_FE_error;
  std::ifstream file(filename, ios::in);
  SF.push_back(eta_bin_SM_center);
  SF.push_back(eta_bin_SM_error);
  SF.push_back(SF_uncorrelated_SM);
  SF.push_back(SF_uncorrelated_SM_error);
  SF.push_back(SF_correlated_SM);
  SF.push_back(SF_correlated_SM_error);
  SF.push_back(SF_uncorrelated_SM_ptdep_min);
  SF.push_back(SF_uncorrelated_SM_ptdep_max);
  SF.push_back(SF_correlated_SM_ptdep_min);
  SF.push_back(SF_correlated_SM_ptdep_max);

  std::string line;
  getline(file, line);

  if (gSystem->AccessPathName(filename)) { std::cout << "check: " << filename << '\n'; return false;}

  while (!file.eof()) {
    getline(file, line);
    std::istringstream iss(line);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
    for (unsigned int i = 0; i < results.size(); i++) { SF.at(i).push_back(std::stold(results[i]));}
  }

  SFs.push_back(SF);
  SF.clear();
  file.close();
}

void Load_all_SF(std::vector<std::vector<std::vector<double>>> & SFs, TString path, TString central_SF_txt, std::vector<TString> systematics_name) {
  TString filename = path+central_SF_txt;
  //std::cout << filename << '\n';
  LoadSF(SFs, filename);
  for (unsigned int i = 0; i < systematics_name.size(); i++) {
    TString temp = central_SF_txt.Copy();
    TString sys = systematics_name.at(i);
    // pt_dep values are saved in all the folder. But the one in the nominal folder are calculated wrt to nominal values. It is saved only the Max and Min value of the difference of the SF
    if (sys.Contains("gaustails")) filename = path+temp.ReplaceAll("standard",sys);
    else if (sys.Contains("_")) filename = path+temp.ReplaceAll("standard",sys(0,sys.First("_"))+"/"+sys(sys.First("_")+1,sys.Length()));
    else if (sys.Contains("pTdep")) filename = path+central_SF_txt;
    else filename = path+temp.ReplaceAll("standard",sys);
    //std::cout << filename << '\n';
    std::cout << filename  << " " << sys << '\n';
    LoadSF(SFs, filename);
  }
}

void GetValuesAndUncertainties(std::vector <double> &SF, std::vector <double> &stat, std::vector <double> &eta_bin_center, std::vector <double> &eta_bin_error, std::vector <std::vector <double> > &systematics, std::vector<std::vector<std::vector<double>>> & SFs, std::vector<TString> systematics_name, double etabins, double method, double pt_dep_method) {
  for (unsigned int i = 0; i < 1+systematics_name.size(); i++){
    if (i==0) {
      // Getting Central Values & Statistical Uncertainties
      for (unsigned int j = 0; j < etabins; j++){
        eta_bin_center.push_back(SFs.at(0).at(0).at(j));
        eta_bin_error.push_back(SFs.at(0).at(1).at(j));
        SF.push_back(SFs.at(0).at(method).at(j));
        stat.push_back(SFs.at(0).at(method+1).at(j));
      }
    } else if (i < systematics_name.size()) {
      // Getting Systematics Uncertainties. Diff wrt the nominal value
      std::vector <double> temp;
      // for (unsigned int j = 0; j < etabins; j++) temp.push_back(TMath::Abs(SFs.at(0).at(method).at(j) - SFs.at(i).at(method).at(j)));
      for (unsigned int j = 0; j < etabins; j++) temp.push_back(SFs.at(i).at(method).at(j) - SFs.at(0).at(method).at(j));
      systematics.push_back(temp);
      temp.clear();
    }
    else{
      std::vector <double> temp;
      // Getting Systematics Uncertainties for pt_dep. It is saved only the Max and Min value
      for (unsigned int j = 0; j < etabins; j++) temp.push_back((SFs.at(i).at(method+pt_dep_method).at(j) + SFs.at(i).at(method+pt_dep_method+1).at(j))/2);
      systematics.push_back(temp);
      temp.clear();
    }
  }
}

void CalculateTotalError(std::vector <double> &total_error, std::vector <double> &systematics_all, std::vector <std::vector <double> > systematics, std::vector <double> &stat, std::vector<TString> systematics_name, double etabins, std::vector <double> SF, std::vector <double> eta_bin_center, std::vector <double> eta_bin_err) {
  TString text = " \t\t\t & stat ";
  for (unsigned int j = 0; j < systematics_name.size(); j++){
    text +=  " & " + systematics_name.at(j);
  }
  text += "\\\\";
  std::cout << text << '\n';
  for (unsigned int i = 0; i < etabins; i++) {
    double err = 0;
    text = Form("$[ %.2f-%.2f]$ & $%.2f $ & $ %.2f $", eta_bin_center.at(i)-eta_bin_err.at(i), eta_bin_center.at(i)+eta_bin_err.at(i), SF.at(i), stat.at(i)/SF.at(i)*100);
    for (unsigned int j = 0; j < systematics_name.size(); j++) {
      err += TMath::Power(systematics.at(j).at(i),2);
      // std::cout << "err FE " << i << " " << systematics_name.at(j) << " " << systematics.at(j).at(i)/SF.at(i)*100 <<  " " << TMath::Sqrt(err)/SF.at(i)*100 << '\n';
      text += Form(" & $ %.2f $", 100*systematics.at(j).at(i)/SF.at(i));
    }
    err = TMath::Sqrt(err);
    systematics_all.push_back(err);
    total_error.push_back(TMath::Sqrt(TMath::Power(err,2)+TMath::Power(stat.at(i),2)));
    text += Form("&$%.2f$", 100*TMath::Sqrt(TMath::Power(err,2)+TMath::Power(stat.at(i),2)-TMath::Power(systematics.at(0).at(i),2))/SF.at(i));
    text += Form(" & $ %.2f $\\\\", 100*TMath::Sqrt(TMath::Power(err,2)+TMath::Power(stat.at(i),2))/SF.at(i));
    std::cout << text << '\n';
  }
}

void Plot_Uncertainties(TString name_method, std::vector<double> eta_bins, std::vector<double> eta_bin_center, std::vector<double> SF, std::vector<double> eta_bin_error, std::vector<double> stat, std::vector<double> systematics_all, std::vector <std::vector <double> > systematics, std::vector<double> total_error, TString path, std::vector<TString> systematics_name,std::vector<double> SF_2) {
  std::vector<int> colors;
  colors.push_back(kBlue-4);
  colors.push_back(kGreen-2);
  colors.push_back(kOrange);
  colors.push_back(kViolet-3);
  colors.push_back(kCyan);
  colors.push_back(kMagenta);
  colors.push_back(kAzure+7);
  colors.push_back(kSpring);
  colors.push_back(kBlack);
  colors.push_back(kCyan+2);


  TCanvas* canv_SF = tdrCanvas("Statistics_"+name_method, eta_bins[0]-0.1, eta_bins[eta_bins.size()-1]+0.1, 0.5, 3.0, "#eta", "JER SF");
  canv_SF->SetTickx(0);
  canv_SF->SetTicky(0);
  TLegend *leg = tdrLeg(0.60,0.67,0.80,0.92, 0.025, 42, kBlack);
  tdrHeader(leg,"Uncertainties", 12);
  TGraphErrors* gr_stat = new TGraphErrors(eta_bins.size()-1, &(eta_bin_center[0]), &SF[0], &(eta_bin_error[0]), &stat[0]);
  TGraphErrors* gr_syst = new TGraphErrors(eta_bins.size()-1, &(eta_bin_center[0]), &SF[0], &(eta_bin_error[0]), &systematics_all[0]);
  TGraphErrors* gr_tot  = new TGraphErrors(eta_bins.size()-1, &(eta_bin_center[0]), &SF[0], &(eta_bin_error[0]), &total_error[0]);
  leg->AddEntry(gr_stat, "stat","f");
  leg->AddEntry(gr_syst, "syst","f");
  leg->AddEntry(gr_tot, "stat+syst","f");
  tdrDraw(gr_tot,  "P5", kFullDotLarge, kGreen-2, kSolid, kGreen-2, 3005, kGreen-2);
  tdrDraw(gr_syst, "P5", kFullDotLarge, kBlue-4,  kSolid, kBlue-4,  3005, kBlue-4);
  tdrDraw(gr_stat, "P5", kFullDotLarge, kRed+1,   kSolid, kRed+1,   3005, kRed+1);
  leg->Draw("same");

  canv_SF->Print(path+"SF_"+name_method+".pdf","pdf");
  canv_SF->Print(path+"SF_"+name_method+".png","png");

  TCanvas* canv_stat = tdrCanvas("Uncertainties_"+name_method,eta_bins[0]-0.1, eta_bins[eta_bins.size()-1]+0.1, 0.0001, 100.0, "#eta", "Uncertainties");
  canv_stat->SetTickx(0);
  canv_stat->SetTicky(0);
  canv_stat->SetLogy();
  leg = tdrLeg(0.55,0.67,0.80,0.92, 0.025, 42, kBlack);
  leg->SetNColumns(2);
  tdrHeader(leg,"Uncertainties", 12);
  TH1* h = new TH1F("", "", eta_bins.size()-1, &eta_bins[0] );
  for (unsigned int i = 0; i < eta_bins.size()-1; i++) { h->SetBinContent(i+1,stat.at(i));}
  leg->AddEntry(h, "stat","l");
  for (unsigned int j = 0; j <= systematics_name.size(); j++){
    TH1F* h = new TH1F("", "", eta_bins.size()-1, &eta_bins[0] );
    h->SetLineWidth(5);
    if (j != systematics_name.size()) {
      for (unsigned int i = 0; i < eta_bins.size()-1; i++) {h->SetBinContent(i+1,TMath::Abs(systematics.at(j).at(i)));}
      tdrDraw(h, "hist", kFullDotLarge, colors[j], kSolid, colors[j], 0, colors[j]);
      leg->AddEntry(h, systematics_name[j],"l");
    } else {
      int shift = shift_barrel;
      if (SF.size() > SF_2.size()) shift = - shift_barrel;
      for (unsigned int i = 0; i < eta_bins.size()-1; i++) {
        if (i < shift_FE && SF.size() < SF_2.size()) continue;
        if (i < shift_FE + shift && SF.size() > SF_2.size()) continue;
        if (i<SF_2.size()-1){h->SetBinContent(i+1, (TMath::Abs(SF.at(i)-SF_2.at(i-shift))/2));}
      }
      tdrDraw(h, "hist", kFullDotLarge, colors[j], kSolid, colors[j], 0, colors[j]);
      leg->AddEntry(h, "|SM-FE|/2","l");
    }
  }

  tdrDraw(h, "hist", kFullDotLarge, kRed+1, kSolid, kRed+1, 0, kRed+1);
  h->SetLineWidth(8); h->SetLineStyle(9);
  leg->Draw("same");

  canv_stat->Print(path+"Uncertainties_"+name_method+".pdf","pdf");
  canv_stat->Print(path+"Uncertainties_"+name_method+".png","png");

}


// SFs == (n_systematics, columns in files, eta_bins)
void plot_SF_systematics_(TString path = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_10_2_10/src/UHH2/JERSF/Analysis/JER/wide_eta_binning/file/Single/Fall17_17Nov2017_V10/AK4CHS/", TString QCD_DATA = "QCDPt/RunBCDEF/") {
  gPrintViaErrorHandler = kTRUE;
  gErrorIgnoreLevel = kFatal;

  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_13TeV = "[MC 102X] Run2018 41.53 fb^{-1}";
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  if (QCD_DATA.Contains("RunA"))     lumi_13TeV = "[MC 102X] Run2018A     14.00 fb^{-1}";
  if (QCD_DATA.Contains("RunB"))     lumi_13TeV = "[MC 102X] Run2018B      7.10 fb^{-1}";
  if (QCD_DATA.Contains("RunC"))     lumi_13TeV = "[MC 102X] Run2018C      6.94 fb^{-1}";
  if (QCD_DATA.Contains("RunD"))     lumi_13TeV = "[MC 102X] Run2018D     31.93 fb^{-1}";
  if (QCD_DATA.Contains("RunAB"))    lumi_13TeV = "[MC 102X] Run2018AB    14.10 fb^{-1}";
  if (QCD_DATA.Contains("RunABC"))   lumi_13TeV = "[MC 102X] Run2018ABC   28.04 fb^{-1}";
  if (QCD_DATA.Contains("RunABCD"))  lumi_13TeV = "[MC 102X] Run2018ABCD  59.97 fb^{-1}";

  //lumi_13TeV = "[MC Pythia8] RunII";


  std::vector<double> eta_bins_all(eta_bins_JER,                eta_bins_JER + sizeof(eta_bins_JER)/sizeof(double));
  std::vector<double> eta_bins_SM(eta_bins_JER,                 eta_bins_JER + sizeof(eta_bins_JER)/sizeof(double) - shift_SM);
  std::vector<double> eta_bins_FE(eta_bins_JER + shift_barrel,  eta_bins_JER + sizeof(eta_bins_JER)/sizeof(double));

  // std::vector<double> eta_bins_all(eta_bins2,                eta_bins2 + sizeof(eta_bins2)/sizeof(double));
  // std::vector<double> eta_bins_SM(eta_bins2,                 eta_bins2 + sizeof(eta_bins2)/sizeof(double) - shift_SM);
  // std::vector<double> eta_bins_FE(eta_bins2 + shift_barrel,  eta_bins2 + sizeof(eta_bins2)/sizeof(double));

  std::cout << "eta_bins_all " << eta_bins_all.size() << std::endl; for (size_t i = 0; i < eta_bins_all.size(); i++) std::cout << eta_bins_all[i] << " "; std::cout << '\n';
  std::cout << "eta_bins_SM " << eta_bins_SM.size() << std::endl; for (size_t i = 0; i < eta_bins_SM.size(); i++) std::cout << eta_bins_SM[i] << " "; std::cout << '\n';
  std::cout << "eta_bins_FE " << eta_bins_FE.size() << std::endl; for (size_t i = 0; i < eta_bins_FE.size(); i++) std::cout << eta_bins_FE[i] << " "; std::cout << '\n';

  int method = 4; //2-uncorr 4-corr
  int pt_dep_method = 4; //4-min value 5-max value

  std::vector<TString> systematics_name;

  // systematics_name.push_back("gaustails");
  systematics_name.push_back("gaustails_0.95");
  systematics_name.push_back("JEC_up");
  systematics_name.push_back("JEC_down");
  systematics_name.push_back("PU_up");
  systematics_name.push_back("PU_down");
  systematics_name.push_back("PLI_down");
  systematics_name.push_back("PLI_up");
  systematics_name.push_back("alpha");
  systematics_name.push_back("pTdep");
  std::vector<std::vector<std::vector<double>>> SFs_SM, SFs_FE;
  TString central_SM, central_FE, filename;

  // path = "AK4CHS/";
  central_SM = "standard/"+QCD_DATA+"output/scalefactors_ST.txt";
  central_FE = "standard/"+QCD_DATA+"output/scalefactors_FE.txt";

  /////////////////////////
  //   Standard Method   //
  /////////////////////////

  Load_all_SF(SFs_SM, path, central_SM, systematics_name);

  // for (size_t i = 0; i < SFs_SM.size(); i++) {
  //   for (size_t j = 0; j < SFs_SM[i].size(); j++) {
  //     for (size_t k = 0; k < SFs_SM[i][j].size(); k++) {
  //       std::cout << SFs_SM.size() << "\t" << SFs_SM[i].size() << "\t" << SFs_SM[i][j].size() << "\t" << i << " " << j << " " << k << "\t" << SFs_SM[i][j][k] << '\n';
  //     }
  //   }
  // }

  if (SFs_SM.size() != systematics_name.size()+1) {
    std::cout << "ERROR SM" << '\n';
    std::cout << SFs_SM.size() << '\n';
    return false;
  }

  std::vector <double> SF_SM, stat_SM, eta_bin_SM_center, eta_bin_SM_error;
  std::vector <std::vector <double> > systematics_SM;
  std::vector <double> total_error_SM, systematics_SM_all;
  GetValuesAndUncertainties(SF_SM, stat_SM, eta_bin_SM_center, eta_bin_SM_error, systematics_SM, SFs_SM, systematics_name, eta_bins_SM.size()-1, method, pt_dep_method);
  CalculateTotalError(total_error_SM, systematics_SM_all, systematics_SM, stat_SM, systematics_name, eta_bins_SM.size()-1, SF_SM, eta_bin_SM_center, eta_bin_SM_error);

  /////////////////////////
  //  Forward Extention  //
  /////////////////////////

  Load_all_SF(SFs_FE, path, central_FE, systematics_name);

  if (SFs_FE.size() != systematics_name.size()+1) {
    std::cout << "ERROR FE" << '\n';
    return false;
  }

  //
  // for (size_t i = 0; i < SFs_FE.size(); i++) {
  //   for (size_t j = 0; j < SFs_FE[i].size(); j++) {
  //     for (size_t k = 0; k < SFs_FE[i][j].size(); k++) {
  //       std::cout << SFs_FE.size() << "\t" << SFs_FE[i].size() << "\t" << SFs_FE[i][j].size() << "\t" << i << " " << j << " " << k << "\t" << SFs_FE[i][j][k] << '\n';
  //     }
  //   }
  // }

  std::vector <double> SF_FE, stat_FE, eta_bin_FE_center, eta_bin_FE_error;
  std::vector <std::vector <double> > systematics_FE;
  std::vector <double> total_error_FE, systematics_FE_all;
  GetValuesAndUncertainties(SF_FE, stat_FE, eta_bin_FE_center, eta_bin_FE_error, systematics_FE, SFs_FE, systematics_name, eta_bins_FE.size()-1, method, pt_dep_method);
  CalculateTotalError(total_error_FE, systematics_FE_all, systematics_FE, stat_FE, systematics_name, eta_bins_FE.size()-1, SF_FE, eta_bin_FE_center, eta_bin_FE_error);


  Plot_Uncertainties("SM", eta_bins_SM, eta_bin_SM_center, SF_SM, eta_bin_SM_error, stat_SM, systematics_SM_all, systematics_SM, total_error_SM, path+"standard/"+QCD_DATA, systematics_name, SF_FE);
  Plot_Uncertainties("FE", eta_bins_FE, eta_bin_FE_center, SF_FE, eta_bin_FE_error, stat_FE, systematics_FE_all, systematics_FE, total_error_FE, path+"standard/"+QCD_DATA, systematics_name, SF_SM);

  /////////////////////////
  //     Combination     //
  /////////////////////////

  std::vector <double> SF_final, SF_final_error, SF_final_error_stat, SF_final_error_syst, eta_bin_all_center, eta_bin_all_error;
  SF_final.clear();
  SF_final_error.clear();
  SF_final_error_stat.clear();
  SF_final_error_syst.clear();
  eta_bin_all_center.clear();
  eta_bin_all_error.clear();

  for (unsigned int i = 0; i < eta_bins_all.size()-1; i++) {
    if (i < shift_FE){
      /***/eta_bin_all_center.push_back(eta_bin_SM_center.at(i));
      /****/eta_bin_all_error.push_back(eta_bin_SM_error.at(i));
      /*************/SF_final.push_back(SF_SM.at(i));
      /*******/SF_final_error.push_back(total_error_SM.at(i));
      /**/SF_final_error_stat.push_back(stat_SM.at(i));
      /**/SF_final_error_syst.push_back(systematics_SM_all.at(i));
    } else if (i < eta_bins_all.size() - 1 - shift_SM){
      /***/eta_bin_all_center.push_back(eta_bin_SM_center.at(i));
      /****/eta_bin_all_error.push_back(eta_bin_SM_error.at(i));
      /*************/SF_final.push_back((SF_SM.at(i)+SF_FE.at(i-shift_FE+shift_barrel+1))/2);
      // /****/SF_final_error.push_back(TMath::Sqrt(TMath::Power((stat_SM.at(i)+stat_FE.at(i-shift_FE+shift_barrel+1))/2, 2)+TMath::Power((systematics_SM_all.at(i)+systematics_FE_all.at(i-shift_FE+shift_barrel+1))/2, 2)+TMath::Power((SF_SM.at(i)-SF_FE.at(i-shift_FE+shift_barrel+1))/2, 2)));
      /*******/SF_final_error.push_back(TMath::Sqrt(stat_SM.at(i)*stat_SM.at(i)+stat_FE.at(i-shift_FE+shift_barrel+1)*stat_FE.at(i-shift_FE+shift_barrel+1) +systematics_SM_all.at(i)*systematics_SM_all.at(i)+systematics_FE_all.at(i-shift_FE+shift_barrel+1)*systematics_FE_all.at(i-shift_FE+shift_barrel+1)+TMath::Power((SF_SM.at(i)-SF_FE.at(i-shift_FE+shift_barrel+1))/2, 2)));
      /**/SF_final_error_stat.push_back(TMath::Sqrt(stat_SM.at(i)*stat_SM.at(i)+stat_FE.at(i-shift_FE+shift_barrel+1)*stat_FE.at(i-shift_FE+shift_barrel+1) ));
      /**/SF_final_error_syst.push_back(TMath::Sqrt(TMath::Power(systematics_SM_all.at(i),2)+TMath::Power(systematics_FE_all.at(i-shift_FE+shift_barrel+1),2)+TMath::Power((SF_SM.at(i)-SF_FE.at(i-shift_FE+shift_barrel+1))/2, 2)));
    } else {
      /***/eta_bin_all_center.push_back(eta_bin_FE_center.at(i-1));
      /****/eta_bin_all_error.push_back(eta_bin_FE_error.at(i-1));
      /*************/SF_final.push_back(SF_FE.at(i-shift_barrel));
      /*******/SF_final_error.push_back(total_error_FE.at(i-shift_barrel));
      /**/SF_final_error_stat.push_back(stat_FE.at(i-shift_barrel));
      /**/SF_final_error_syst.push_back(systematics_FE_all.at(i-shift_barrel));
    }
  }


  TCanvas* canv_SF_final = tdrCanvas("SF_final", eta_bins_all.at(0)-0.1, eta_bins_all.at(eta_bins_all.size()-1)+0.5, 0., 3.0, "#eta", "JER SF");
  canv_SF_final->SetTickx(0);
  canv_SF_final->SetTicky(0);
  TLegend *leg_final = tdrLeg(0.64,0.67,0.79,0.92, 0.040, 42, kBlack);
  tdrHeader(leg_final,"", 12);

  double etaSummer16_25nsV1[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.853, 2.964, 3.139, 5.191};
  // double jerSummer16_25nsV1[13][2] = {{1.15953,0.0322243894},{1.19477,0.0337499818}, {1.1464, 0.0298593582},{1.1608608730818335, 0.0314114337}, {1.1277674503328965, 0.0227705413},{1.1000419290738948, 0.0448277854}, {1.1426190202932343, 0.0749862918},{1.1511576138738635, 0.0431038630}, {1.2962786307493799, 0.1221696907},{1.3418116992379743, 0.0896812121}, {1.77881,0.2007462079},{1.18695,0.1243701331}, {1.19218,0.1487939851}};
  double jerSummer16_25nsV1[13][2] = {{1.1595,0.0645},{1.1948,0.0652},{1.1464,0.0632},{1.1609,0.1025},{1.1278,0.0986},{1.1000,0.1079},{1.1426,0.1214},{1.1512,0.1140},{1.2963,0.2371},{1.3418,0.2091},{1.7788,0.2008},{1.1869,0.1243},{1.1922,0.1488}};
  // double jerSummer16_25nsV1[13][2] = {{1.17716, 0.04438},{1.21224, 0.04054},{1.14975, 0.082},{1.1395, 0.06360},{1.15978, 0.1112},{1.20805, 0.089},{1.39324, 0.302},{1.32341, 0.183},{1.58005, 0.598},{2.26888, 0.261},{2.65, 0.337},{1.27321, 0.103},{1.20094, 0.105}};

  // double etaSpring16_25nsV10[] = {0, 0.5, 0.8, 1.11, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0};
  // double etaSpring16_25nsV10[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.853, 2.964, 3.139, 5.191};
  // double jerSpring16_25nsV10[13][2] = {{1.109,0.008},{1.138,0.013},{1.114,0.013},{1.123,0.024},{1.084,0.011},{1.082,0.035},{1.140,0.047},{1.067,0.053},{1.177,0.041},{1.364,0.039},{1.857,0.071},{1.328,0.022},{1.16,0.029}};

  double etaFall17_V3[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.853, 2.964, 3.139, 5.191};
  double jerFall17_V3[13][2] = {{1.1432,0.0222},{1.1815,0.0484},{1.0989,0.0456},{1.1137,0.1397},{1.1307,0.1470},{1.1600,0.0976},{1.2393,0.1909},{1.2604,0.1501},{1.4085,0.2020},{1.9909,0.5684},{2.2923,0.3743},{1.2696,0.1089},{1.1542,0.1524}};

  double etaAutumn18_V1[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.853, 2.964, 3.139, 5.191};
  // RUNABCD
  //double jerAutumn18_V1[13][2] = {{1.15,0.043},{1.134,0.08},{1.102,0.052},{1.134,0.112},{1.104,0.211},{1.149,0.159},{1.148,0.209},{1.114,0.191},{1.347,0.274},{2.137,0.524},{1.65,0.941},{1.225,0.194},{1.082,0.198}};
  // double jerAutumn18_V1[13][2] = {{1.15,0.008},{1.134,0.018},{1.102,0.014},{1.134,0.033},{1.104,0.019},{1.149,0.02},{1.148,0.041},{1.114,0.025},{1.347,0.053},{2.137,0.102},{1.65,0.18},{1.225,0.054},{1.082,0.057}}; // stat
  // double jerAutumn18_V1[13][2] = {{1.15,0.042},{1.134,0.078},{1.102,0.05},{1.134,0.107},{1.104,0.21},{1.149,0.157},{1.148,0.205},{1.114,0.19},{1.347,0.269},{2.137,0.514},{1.65,0.923},{1.225,0.187},{1.082,0.189}}; // sys
  // RUND
  double jerAutumn18_V1[13][2] = {{1.1401,0.0323},{1.1370,0.0969},{1.1109,0.0676},{1.1585,0.1111},{1.0997,0.2897},{1.1428,0.2264},{1.1157,0.2216},{1.0903,0.2131},{1.4930,0.3106},{2.4518,0.5757},{1.8935,0.9319},{1.1826,0.2058},{1.0439,0.2040}};
  // double jerAutumn18_V1[13][2] = {{1.1401,0.0082},{1.1370,0.0183},{1.1109,0.0142},{1.1585,0.0344},{1.0997,0.0193},{1.1428,0.0210},{1.1157,0.0416},{1.0903,0.0257},{1.4930,0.0579},{2.4518,0.1171},{1.8935,0.2069},{1.1826,0.0518},{1.0439,0.0551}}; // stat
  // double jerAutumn18_V1[13][2] = {{1.1401,0.0313},{1.1370,0.0952},{1.1109,0.0660},{1.1585,0.1056},{1.0997,0.2890},{1.1428,0.2254},{1.1157,0.2177},{1.0903,0.2115},{1.4930,0.3051},{2.4518,0.5636},{1.8935,0.9087},{1.1826,0.1992},{1.0439,0.1964}}; // sys
  // RUNABC
  //double jerAutumn18_V1[13][2] = {{1.1609,0.0552},{1.1309,0.0610},{1.0918,0.0352},{1.1064,0.1131},{1.1097,0.1206},{1.1554,0.0820},{1.1843,0.1951},{1.1401,0.1669},{1.1818,0.2327},{1.7778,0.4659},{1.3718,0.9507},{1.2725,0.1810},{1.1255,0.1901}};
  // double jerAutumn18_V1[13][2] = {{1.1609,0.0077},{1.1309,0.0170},{1.0918,0.0129},{1.1064,0.0307},{1.1097,0.0177},{1.1554,0.0196},{1.1843,0.0412},{1.1401,0.0239},{1.1818,0.0467},{1.7778,0.0842},{1.3718,0.1489},{1.2725,0.0556},{1.1255,0.0589}}; // stat
  // double jerAutumn18_V1[13][2] = {{1.1609,0.0547},{1.1309,0.0586},{1.0918,0.0327},{1.1064,0.1089},{1.1097,0.1193},{1.1554,0.0797},{1.1843,0.1907},{1.1401,0.1652},{1.1818,0.2280},{1.7778,0.4582},{1.3718,0.9390},{1.2725,0.1723},{1.1255,0.1807}}; // sys


  double etaAutumn18_V3[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.853, 2.964, 3.139, 5.191};
  double jerAutumn18_V3[13][2] = {{1.1813,0.0639},{1.1136,0.0978},{1.1048,0.0328},{1.0741,0.0526},{1.0923,0.0953},{1.0779,0.0593},{1.0893,0.1941},{1.0755,0.1270},{1.4188,0.2210},{1.9206,0.4811},{2.0118,0.2916},{1.1904,0.1223},{1.0846,0.2697}};
  // double jerAutumn18_V3[13][2] = {{1.1813,0.0031},{1.1136,0.0068},{1.1048,0.0055},{1.0741,0.0133},{1.0923,0.0075},{1.0779,0.0081},{1.0893,0.0130},{1.0755,0.0101},{1.4188,0.0292},{1.9206,0.0214},{2.0118,0.0369},{1.1904,0.0122},{1.0846,0.0137}}; // stat
  // double jerAutumn18_V3[13][2] = {{1.1813,0.0638},{1.1136,0.0976},{1.1048,0.0324},{1.0741,0.0509},{1.0923,0.0950},{1.0779,0.0587},{1.0893,0.1936},{1.0755,0.1266},{1.4188,0.2191},{1.9206,0.4807},{2.0118,0.2893},{1.1904,0.1217},{1.0846,0.2694}}; // syst


  double etaAutumn18_V4[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.853, 2.964, 3.139, 5.191};
  // double jerAutumn18_V4[13][2] = {{1.1588,0.0392},{1.1504,0.0728},{1.1253,0.0343},{1.1217,0.0826},{1.1069,0.0890},{1.0916,0.0430},{1.0977,0.0951},{1.1177,0.0721},{1.4494,0.1502},{2.3588,0.6411},{2.2520,0.3542},{1.1759,0.0540},{1.0777,0.0542}}; // RunD
  double jerAutumn18_V4[13][2] = {{1.1677,0.0477},{1.1475,0.0520},{1.1029,0.0310},{1.0781,0.0950},{1.1006,0.0846},{1.1019,0.0386},{1.0459,0.1578},{1.1612,0.0646},{1.2299,0.1087},{1.6736,0.3792},{1.7292,0.2007},{1.2257,0.0452},{1.0733,0.0676}}; // RunABC
  //double jerAutumn18_V4[13][2] = {{1.1545,0.0308},{1.1481,0.0515},{1.0998,0.0386},{1.0929,0.0856},{1.1093,0.0718},{1.1005,0.0515},{1.0603,0.1301},{1.1287,0.0531},{1.3397,0.1147},{2.0325,0.5361},{2.0567,0.3060},{1.1868,0.0376},{1.0922,0.0489}}; // RunABCD

  double etaAutumn18_V4_wPUId[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.853, 2.964, 3.139, 5.191};
  double jerAutumn18_V4_wPUId[13][2] = {{1.1560,0.0161},{1.1588,0.0781},{1.1399,0.0870},{1.1223,0.0477},{1.1549,0.0776},{1.1023,0.0768},{1.1543,0.0952},{1.1626,0.0710},{1.3146,0.1141},{1.9192,0.3146},{2.0417,0.3043},{1.2588,0.0401},{1.1573,0.1365}}; // RunABCD
  double jerAutumn18_V4_noPUId[13][2] = {{1.1545,0.0308},{1.1481,0.0515},{1.0998,0.0386},{1.0929,0.0856},{1.1093,0.0718},{1.1005,0.0515},{1.0603,0.1301},{1.1287,0.0531},{1.3397,0.1147},{2.0325,0.5361},{2.0567,0.3060},{1.1868,0.0376},{1.0922,0.0489}}; // RunABCD


  double etaAutumn18_V5[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
  double jerAutumn18_V5[14][2] = {{1.0962,0.0708},{1.1206,0.0519},{1.0429,0.0357},{1.0626,0.0519},{1.1010,0.0573},{1.0391,0.0708},{1.0263,0.0714},{1.0883,0.0400},{1.1464,0.1022},{1.3271,0.0856},{1.7342,0.2126},{1.7644,0.2034},{1.2714,0.0433},{1.0248,0.1092}};//RunABC

  double etaAutumn18_V5_1[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
  double jerAutumn18_V5_1[14][2] = {{1.1631,0.0366},{1.1680,0.0564},{1.1300,0.0706},{1.1187,0.0441},{1.1393,0.0728},{1.1054,0.0764},{1.1242,0.0258},{1.1837,0.1065},{1.2067,0.0595},{1.5080,0.1476},{1.9058,0.2898},{1.8627,0.1835},{1.3025,0.0532},{1.1315,0.1257}};//RunABC

  double etaAutumn18_V6[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
  double jerAutumn18_V6[14][2] = {{1.1653,0.0273},{1.1929,0.0490},{1.1171,0.0474},{1.1071,0.0740},{1.1338,0.0567},{1.0909,0.0792},{1.1194,0.0491},{1.1889,0.1341},{1.1850,0.0401},{1.4596,0.1404},{1.8143,0.2356},{1.7414,0.2295},{1.2764,0.0606},{1.0986,0.1152}};//RunABC
  // double jerAutumn18_V6[14][2] = {{1.1639,0.0312},{1.1808,0.0168},{1.1340,0.0530},{1.1316,0.0738},{1.1389,0.0766},{1.1310,0.0600},{1.1388,0.0405},{1.1827,0.1309},{1.2874,0.0616},{1.8381,0.3179},{2.2087,0.2865},{2.1285,0.3096},{1.2063,0.0817},{1.0760,0.1492}};//RunABCD
  // double jerAutumn18_V6[14][2] = {{1.1742,0.0327},{1.1930,0.0372},{1.1451,0.0639},{1.1413,0.0601},{1.1396,0.0678},{1.1661,0.0905},{1.1581,0.0693},{1.1912,0.1294},{1.4015,0.1590},{2.1665,0.4095},{2.5335,0.3138},{2.3392,0.3996},{1.2284,0.1264},{1.0495,0.1725}};//RunD




  std::vector<double> etaSummer16_25nsV1_center, etaSummer16_25nsV1_err, SFSummer16_25nsV1, SFSummer16_25nsV1_Err;
  // std::vector<double> etaSpring16_25nsV10_center, etaSpring16_25nsV10_err, SFSpring16_25nsV10, SFSpring16_25nsV10_Err;
  std::vector<double> etaFall17_V3_center, etaFall17_V3_err, SFFall17_V3, SFFall17_V3_Err;
  std::vector<double> etaAutumn18_V1_center, etaAutumn18_V1_err, SFAutumn18_V1, SFAutumn18_V1_Err;
  std::vector<double> etaAutumn18_V3_center, etaAutumn18_V3_err, SFAutumn18_V3, SFAutumn18_V3_Err;
  std::vector<double> etaAutumn18_V4_center, etaAutumn18_V4_err, SFAutumn18_V4, SFAutumn18_V4_Err;
  std::vector<double> etaAutumn18_V4_wPUId_center, etaAutumn18_V4_wPUId_err, SFAutumn18_V4_wPUId, SFAutumn18_V4_wPUId_Err;
  std::vector<double> etaAutumn18_V4_noPUId_center, etaAutumn18_V4_noPUId_err, SFAutumn18_V4_noPUId, SFAutumn18_V4_noPUId_Err;
  std::vector<double> etaAutumn18_V5_center, etaAutumn18_V5_err, SFAutumn18_V5, SFAutumn18_V5_Err;
  std::vector<double> etaAutumn18_V5_1_center, etaAutumn18_V5_1_err, SFAutumn18_V5_1, SFAutumn18_V5_1_Err;
  std::vector<double> etaAutumn18_V6_center, etaAutumn18_V6_err, SFAutumn18_V6, SFAutumn18_V6_Err;


  for (unsigned int i = 0; i < 13; i++) {
    etaSummer16_25nsV1_center.push_back((etaSummer16_25nsV1[i+1]+etaSummer16_25nsV1[i])/2);
    etaSummer16_25nsV1_err.push_back((etaSummer16_25nsV1[i+1]-etaSummer16_25nsV1[i])/2);
    SFSummer16_25nsV1.push_back(jerSummer16_25nsV1[i][0]);
    SFSummer16_25nsV1_Err.push_back(jerSummer16_25nsV1[i][1]);

    // etaSpring16_25nsV10_center.push_back((etaSpring16_25nsV10[i+1]+etaSpring16_25nsV10[i])/2);
    // etaSpring16_25nsV10_err.push_back((etaSpring16_25nsV10[i+1]-etaSpring16_25nsV10[i])/2);
    // SFSpring16_25nsV10.push_back(jerSpring16_25nsV10[i][0]);
    // SFSpring16_25nsV10_Err.push_back(jerSpring16_25nsV10[i][1]);

    /****/etaFall17_V3_center.push_back((etaFall17_V3[i+1]+etaFall17_V3[i])/2);
    /*******/etaFall17_V3_err.push_back((etaFall17_V3[i+1]-etaFall17_V3[i])/2);
    /************/SFFall17_V3.push_back(jerFall17_V3[i][0]);
    /********/SFFall17_V3_Err.push_back(jerFall17_V3[i][1]);


    /**/etaAutumn18_V1_center.push_back((etaAutumn18_V1[i+1]+etaAutumn18_V1[i])/2);
    /*****/etaAutumn18_V1_err.push_back((etaAutumn18_V1[i+1]-etaAutumn18_V1[i])/2);
    /**********/SFAutumn18_V1.push_back(jerAutumn18_V1[i][0]);
    /******/SFAutumn18_V1_Err.push_back(jerAutumn18_V1[i][1]);

    /**/etaAutumn18_V3_center.push_back((etaAutumn18_V3[i+1]+etaAutumn18_V3[i])/2);
    /*****/etaAutumn18_V3_err.push_back((etaAutumn18_V3[i+1]-etaAutumn18_V3[i])/2);
    /**********/SFAutumn18_V3.push_back(jerAutumn18_V3[i][0]);
    /******/SFAutumn18_V3_Err.push_back(jerAutumn18_V3[i][1]);


    /**/etaAutumn18_V4_center.push_back((etaAutumn18_V4[i+1]+etaAutumn18_V4[i])/2);
    /*****/etaAutumn18_V4_err.push_back((etaAutumn18_V4[i+1]-etaAutumn18_V4[i])/2);
    /**********/SFAutumn18_V4.push_back(jerAutumn18_V4[i][0]);
    /******/SFAutumn18_V4_Err.push_back(jerAutumn18_V4[i][1]);


    /**/etaAutumn18_V4_wPUId_center.push_back((etaAutumn18_V4[i+1]+etaAutumn18_V4[i])/2);
    /*****/etaAutumn18_V4_wPUId_err.push_back((etaAutumn18_V4[i+1]-etaAutumn18_V4[i])/2);
    /**********/SFAutumn18_V4_wPUId.push_back(jerAutumn18_V4_wPUId[i][0]);
    /******/SFAutumn18_V4_wPUId_Err.push_back(jerAutumn18_V4_wPUId[i][1]);



    /**/etaAutumn18_V4_noPUId_center.push_back((etaAutumn18_V4[i+1]+etaAutumn18_V4[i])/2);
    /*****/etaAutumn18_V4_noPUId_err.push_back((etaAutumn18_V4[i+1]-etaAutumn18_V4[i])/2);
    /**********/SFAutumn18_V4_noPUId.push_back(jerAutumn18_V4_noPUId[i][0]);
    /******/SFAutumn18_V4_noPUId_Err.push_back(jerAutumn18_V4_noPUId[i][1]);

  }

  for (unsigned int i = 0; i < 14; i++) {
    /**/etaAutumn18_V5_center.push_back((etaAutumn18_V5[i+1]+etaAutumn18_V5[i])/2);
    /*****/etaAutumn18_V5_err.push_back((etaAutumn18_V5[i+1]-etaAutumn18_V5[i])/2);
    /**********/SFAutumn18_V5.push_back(jerAutumn18_V5[i][0]);
    /******/SFAutumn18_V5_Err.push_back(jerAutumn18_V5[i][1]);

    /**/etaAutumn18_V5_1_center.push_back((etaAutumn18_V5_1[i+1]+etaAutumn18_V5_1[i])/2);
    /*****/etaAutumn18_V5_1_err.push_back((etaAutumn18_V5_1[i+1]-etaAutumn18_V5_1[i])/2);
    /**********/SFAutumn18_V5_1.push_back(jerAutumn18_V5_1[i][0]);
    /******/SFAutumn18_V5_1_Err.push_back(jerAutumn18_V5_1[i][1]);

    /**/etaAutumn18_V6_center.push_back((etaAutumn18_V6[i+1]+etaAutumn18_V6[i])/2);
    /*****/etaAutumn18_V6_err.push_back((etaAutumn18_V6[i+1]-etaAutumn18_V6[i])/2);
    /**********/SFAutumn18_V6.push_back(jerAutumn18_V6[i][0]);
    /******/SFAutumn18_V6_Err.push_back(jerAutumn18_V6[i][1]);

  }


  // TGraphErrors* gr_SFSpring16_25nsV10 = new TGraphErrors(SFSpring16_25nsV10.size(), &(etaSpring16_25nsV10_center[0]), &SFSpring16_25nsV10[0], &(etaSpring16_25nsV10_err[0]), &SFSpring16_25nsV10_Err[0]);
  TGraphErrors* gr_SFSummer16_25nsV1  = new TGraphErrors(SFSummer16_25nsV1.size(), &(etaSummer16_25nsV1_center[0]), &SFSummer16_25nsV1[0], &(etaSummer16_25nsV1_err[0]), &SFSummer16_25nsV1_Err[0]);
  TGraphErrors* gr_SFFall17_V3        = new TGraphErrors(SFFall17_V3.size(), &(etaFall17_V3_center[0]), &SFFall17_V3[0], &(etaFall17_V3_err[0]), &SFFall17_V3_Err[0]);
  TGraphErrors* gr_SFAutumn18_V1      = new TGraphErrors(SFAutumn18_V1.size(), &(etaAutumn18_V1_center[0]), &SFAutumn18_V1[0], &(etaAutumn18_V1_err[0]), &SFAutumn18_V1_Err[0]);
  TGraphErrors* gr_SFAutumn18_V3      = new TGraphErrors(SFAutumn18_V3.size(), &(etaAutumn18_V3_center[0]), &SFAutumn18_V3[0], &(etaAutumn18_V3_err[0]), &SFAutumn18_V3_Err[0]);
  TGraphErrors* gr_SFAutumn18_V4      = new TGraphErrors(SFAutumn18_V4.size(), &(etaAutumn18_V4_center[0]), &SFAutumn18_V4[0], &(etaAutumn18_V4_err[0]), &SFAutumn18_V4_Err[0]);
  TGraphErrors* gr_SFAutumn18_V4_wPUId= new TGraphErrors(SFAutumn18_V4_wPUId.size(), &(etaAutumn18_V4_wPUId_center[0]), &SFAutumn18_V4_wPUId[0], &(etaAutumn18_V4_wPUId_err[0]), &SFAutumn18_V4_wPUId_Err[0]);
  TGraphErrors* gr_SFAutumn18_V4_noPUId= new TGraphErrors(SFAutumn18_V4_noPUId.size(), &(etaAutumn18_V4_noPUId_center[0]), &SFAutumn18_V4_noPUId[0], &(etaAutumn18_V4_noPUId_err[0]), &SFAutumn18_V4_noPUId_Err[0]);
  TGraphErrors* gr_SFAutumn18_V5      = new TGraphErrors(SFAutumn18_V5.size(), &(etaAutumn18_V5_center[0]), &SFAutumn18_V5[0], &(etaAutumn18_V5_err[0]), &SFAutumn18_V5_Err[0]);
  TGraphErrors* gr_SFAutumn18_V5_1    = new TGraphErrors(SFAutumn18_V5_1.size(), &(etaAutumn18_V5_1_center[0]), &SFAutumn18_V5_1[0], &(etaAutumn18_V5_1_err[0]), &SFAutumn18_V5_1_Err[0]);
  TGraphErrors* gr_SFAutumn18_V6      = new TGraphErrors(SFAutumn18_V6.size(), &(etaAutumn18_V6_center[0]), &SFAutumn18_V6[0], &(etaAutumn18_V6_err[0]), &SFAutumn18_V6_Err[0]);
  TGraphErrors* gr_final              = new TGraphErrors(SF_final.size(), &(eta_bin_all_center[0]), &SF_final[0], &(eta_bin_all_error[0]), &SF_final_error[0]); //tot
  // TGraphErrors* gr_final           = new TGraphErrors(SF_final.size(), &(eta_bin_all_center[0]), &SF_final[0], &(eta_bin_all_error[0]), &SF_final_error_stat[0]); //stat
  // TGraphErrors* gr_final           = new TGraphErrors(SF_final.size(), &(eta_bin_all_center[0]), &SF_final[0], &(eta_bin_all_error[0]), &SF_final_error_syst[0]); //sys
  // tdrDraw(gr_SFSpring16_25nsV10, "P5", kFullDotLarge, kOrange-1, kSolid, kOrange-1, 3004, kOrange-1);
  // tdrDraw(gr_SFSummer16_25nsV1, "P5", kFullDotLarge, kRed+1, kSolid, kRed+1, 3005, kRed+1);
  // tdrDraw(gr_SFFall17_V3, "P5", kFullDotLarge, kGreen-1, kSolid, kGreen-1, 3004, kGreen-1);
  // tdrDraw(gr_SFAutumn18_V1, "P5", kFullDotLarge, kGreen-1, kSolid, kGreen-1, 3005, kGreen-1);
  // tdrDraw(gr_SFAutumn18_V3, "P5", kFullDotLarge, kGreen-1, kSolid, kGreen-1, 3005, kGreen-1);
  // tdrDraw(gr_SFAutumn18_V4, "P5", kFullDotLarge, kRed+1, kSolid, kRed+1, 3005, kRed+1);
  tdrDraw(gr_SFAutumn18_V5_1, "P5", kFullDotLarge, kGreen-1, kSolid, kGreen-1, 3005, kGreen-1);
  // tdrDraw(gr_SFAutumn18_V5, "P5", kFullDotLarge, kRed+1, kSolid, kRed+1, 3005, kRed+1);
  tdrDraw(gr_SFAutumn18_V6, "P5", kFullDotLarge, kRed+1, kSolid, kRed+1, 3005, kRed+1);
  // tdrDraw(gr_SFAutumn18_V4_noPUId, "P5", kFullDotLarge, kRed+1, kSolid, kRed+1, 3005, kRed+1);
  tdrDraw(gr_final, "P5", kFullDotLarge, kBlue-4, kSolid, kBlue-4, 3005, kBlue-4);
  // tdrDraw(gr_SFAutumn18_V4_wPUId, "P5", kFullDotLarge, kGreen-1, kSolid, kGreen-1, 3005, kGreen-1);
  // // leg_final->AddEntry(gr_SFSpring16_25nsV10, "Spring16_25nsV10","f");
  // leg_final->AddEntry(gr_SFSummer16_25nsV1, "Summer16_25nsV1","f");
  // leg_final->AddEntry(gr_SFFall17_V3,       "Fall17_V3","f");
  // leg_final->AddEntry(gr_SFAutumn18_V1,  "Autumn18_V1RunABC","f");
  // leg_final->AddEntry(gr_SFAutumn18_V3,  "Autumn18_V3","f");
  // leg_final->AddEntry(gr_SFAutumn18_V4,  "Autumn18_V4RunABC","f");
  // leg_final->AddEntry(gr_SFAutumn18_V5,  "Autumn18_V5","f");
  leg_final->AddEntry(gr_SFAutumn18_V5_1,"Autumn18_JEC15","f");
  leg_final->AddEntry(gr_SFAutumn18_V6,"Autumn18_JEC16","f");
  // leg_final->AddEntry(gr_SFAutumn18_V4_noPUId,  "Autumn18_V4RunABC","f");
  leg_final->AddEntry(gr_final,             "Autumn18_JEC16h","f");
  // leg_final->AddEntry(gr_final,             "Autumn18_JEC16h"+QCD_DATA(QCD_DATA.Index("Run"), QCD_DATA.Length()-QCD_DATA.Index("Run")-1),"f");
  // leg_final->AddEntry(gr_final,             "Autumn18_NoCut","f");
  // leg_final->AddEntry(gr_SFAutumn18_V4_wPUId,  "Autumn18_V4_NoCut_PuId","f");
  // leg_final->AddEntry(gr_SFSummer16_25nsV1, "RunF_ECAL","f");
  // leg_final->AddEntry(gr_final,  "RunF","f");
  leg_final->Draw("same");

  canv_SF_final->Print(path+"standard/"+QCD_DATA+"SF_final.pdf","pdf");
  canv_SF_final->Print(path+"standard/"+QCD_DATA+"SF_final.png","png");

  ofstream SF_file_final;
  std::cout << path+"standard/"+QCD_DATA+"SF_final_tex.txt" << '\n';
  SF_file_final.open (path+"standard/"+QCD_DATA+"SF_final_tex.txt");
  for (unsigned int i = 0; i < SF_final.size(); i++) {
    SF_file_final << "{{" << eta_bin_all_center[i]+eta_bin_all_error[i] << ", " << SF_final[i] << ", " << SF_final[i]+SF_final_error[i] << ", " << SF_final[i]-SF_final_error[i] << "}},\n";
  }

  SF_file_final.close();
  std::cout << path+"standard/"+QCD_DATA+"SF_final.pdf" << '\n';

  std::cout << SF_final.size() << '\n';
  std::cout << path+"/standard/"+QCD_DATA << '\n';
  std::cout << "2017 vs 2018" << '\n';
  for (unsigned int i = 0; i < SFSummer16_25nsV1.size(); i++) {
    std::cout << Form("$[ %.2f-%.2f]$ & $%.2f \\pm %.2f $ \\% & $%.2f \\pm %.2f $ \\% \\\\", eta_bin_all_center[i]-eta_bin_all_error[i], eta_bin_all_center[i]+eta_bin_all_error[i], SFFall17_V3[i], SFFall17_V3_Err[i]/SFFall17_V3[i]*100, SF_final[i], SF_final_error[i]/SF_final[i]*100 ) << std::endl;
  }

  for (unsigned int i = 0; i < SFSummer16_25nsV1.size(); i++) {
    std::cout << Form("[ %.2f-%.2f] \t %.2f +- %.2f \% \t %.2f +- %.2f \% ", eta_bin_all_center[i]-eta_bin_all_error[i], eta_bin_all_center[i]+eta_bin_all_error[i], SFFall17_V3[i], SFFall17_V3_Err[i]/SFFall17_V3[i]*100, SF_final[i], SF_final_error[i]/SF_final[i]*100 ) << std::endl;
  }


  ofstream SF_file_twiki;
  SF_file_twiki.open (path+"standard/"+QCD_DATA+"SF_final_twiki.txt");

  std::cout << '\n' << "|  *abs(eta) region* |";
  SF_file_twiki << "|  *abs(eta) region* |";
  for (size_t i = 0; i < eta_bin_all_center.size(); i++) {
    std::cout << Form("|%.3f-%.3f|", eta_bin_all_center[i]-eta_bin_all_error[i], eta_bin_all_center[i]+eta_bin_all_error[i]);
    SF_file_twiki << Form("|%.3f-%.3f|", eta_bin_all_center[i]-eta_bin_all_error[i], eta_bin_all_center[i]+eta_bin_all_error[i]);
  }

  std::cout << '\n' << "SF_plot" << std::endl; for (size_t i = 0; i < eta_bin_all_center.size(); i++) std::cout << Form("{%.4f,%.4f},", SF_final[i],SF_final_error[i]);
  std::cout << '\n' << "SF_plot Stat" << std::endl; for (size_t i = 0; i < eta_bin_all_center.size(); i++) std::cout << Form("{%.4f,%.4f},", SF_final[i],SF_final_error_stat[i]);
  std::cout << '\n' << "SF_plot Sys" << std::endl; for (size_t i = 0; i < eta_bin_all_center.size(); i++) std::cout << Form("{%.4f,%.4f},", SF_final[i],SF_final_error_syst[i]);



  std::cout << '\n' << "|  *Data/MC SF*      |";
  SF_file_twiki << '\n' << "|  *Data/MC SF*      |";
  for (size_t i = 0; i < eta_bin_all_center.size(); i++) {
    std::cout << Form("|%.4f|", SF_final[i]);
    SF_file_twiki << Form("|%.4f|", SF_final[i]);
  }
  std::cout << '\n' << "|  *Stat.Unc*        |";
  SF_file_twiki << '\n' << "|  *Stat.Unc*        |";
  for (size_t i = 0; i < eta_bin_all_center.size(); i++) {
    std::cout << Form("|%.4f|", SF_final_error_stat[i]);
    SF_file_twiki << Form("|%.4f|", SF_final_error_stat[i]);
  }
  std::cout << '\n' << "|  *Syst.Unc*        |";
  SF_file_twiki << '\n' << "|  *Syst.Unc*        |";
  for (size_t i = 0; i < eta_bin_all_center.size(); i++) {
    std::cout << Form("|%.4f|", SF_final_error_syst[i]);
    SF_file_twiki << Form("|%.4f|", SF_final_error_syst[i]);
  }
  std::cout << '\n' << "|  *Total.Unc*       |";
  SF_file_twiki << '\n' << "|  *Total.Unc*       |";
  for (size_t i = 0; i < eta_bin_all_center.size(); i++) {
    std::cout << Form("|%.4f|", SF_final_error[i]);
    SF_file_twiki << Form("|%.4f|", SF_final_error[i]);
  }

  SF_file_twiki.close();

  ofstream SF_file_DB;
  SF_file_DB.open (path+"standard/"+QCD_DATA+"SF_final_DB.txt");

  std::cout << '\n' << "{1 JetEta 0 None ScaleFactor}" << '\n';
  SF_file_DB << '\n' << "{1 JetEta 0 None ScaleFactor}" << '\n';

  for (int i = SF_final.size()-1; i >= 0 ; i--) {
    std::cout << Form("%.3f %.3f 3 %.4f %.4f %.4f", -(eta_bin_all_center[i]+eta_bin_all_error[i]), -(eta_bin_all_center[i]-eta_bin_all_error[i]), SF_final[i], SF_final[i]-SF_final_error[i], SF_final[i]+SF_final_error[i] ) << std::endl;
    SF_file_DB << Form("%.3f %.3f 3 %.4f %.4f %.4f", -(eta_bin_all_center[i]+eta_bin_all_error[i]), -(eta_bin_all_center[i]-eta_bin_all_error[i]), SF_final[i], SF_final[i]-SF_final_error[i], SF_final[i]+SF_final_error[i] ) << std::endl;
  }
  for (unsigned int i = 0; i < SF_final.size(); i++) {
    std::cout << Form("%.3f %.3f 3 %.4f %.4f %.4f",   eta_bin_all_center[i]-eta_bin_all_error[i],    eta_bin_all_center[i]+eta_bin_all_error[i],  SF_final[i], SF_final[i]-SF_final_error[i], SF_final[i]+SF_final_error[i] ) << std::endl;
    SF_file_DB << Form("%.3f %.3f 3 %.4f %.4f %.4f",   eta_bin_all_center[i]-eta_bin_all_error[i],    eta_bin_all_center[i]+eta_bin_all_error[i],  SF_final[i], SF_final[i]-SF_final_error[i], SF_final[i]+SF_final_error[i] ) << std::endl;
  }

  SF_file_DB.close();

  std::vector<double> systematics_SM_rest(systematics_SM_all.size(),0);
  std::vector<double> systematics_FE_rest(systematics_FE_all.size(),0);
  std::vector<double> systematics_rest(eta_bins_all.size(),0);

  for (unsigned int i = 0; i < systematics_name.size(); i++) {
    TString sysName = systematics_name.at(i);
    std::vector<double> sys_SM(&(systematics_SM.at(i)[0]), &(systematics_SM.at(i)[0]) + systematics_SM.at(i).size());
    std::vector<double> sys_FE(&(systematics_FE.at(i)[0]), &(systematics_FE.at(i)[0]) + systematics_FE.at(i).size());
    std::vector<double> sys_comb(eta_bins_all.size(),0);
    for (size_t i = 0; i < sys_SM.size(); i++) sys_SM[i] /=  SF_SM[i];
    for (size_t i = 0; i < sys_FE.size(); i++) sys_FE[i] /=  SF_FE[i];
    for (unsigned int i = 0; i < eta_bins_all.size()-1; i++) {
      if (i < shift_FE) sys_comb[i] = sys_SM.at(i);
      else if (i < eta_bins_all.size() - 1 - shift_SM) sys_comb[i] = TMath::Sqrt(TMath::Power(sys_SM.at(i),2)+TMath::Power(sys_FE.at(i-shift_FE+shift_barrel+1),2));
      else sys_comb[i] = sys_FE.at(i-shift_barrel);
    }

    if (sysName!="JEC_up" && sysName!="JEC_down" && sysName!="gaustails_0.95") {
      std::cout << "Adding " << sysName << '\n';
      for (unsigned int j = 0; j < systematics_SM_rest.size(); j++) systematics_SM_rest.at(j) += TMath::Power(sys_SM.at(j)/SF_SM[j],2);
      for (unsigned int j = 0; j < systematics_FE_rest.size(); j++) systematics_FE_rest.at(j) += TMath::Power(sys_FE.at(j)/SF_FE[j],2);
    }

    std::vector<double> dummy_SM(eta_bin_SM_center.size(),0.001); std::vector<double> dummy_FE(eta_bin_FE_center.size(),0.001); std::vector<double> dummy(eta_bin_all_center.size(),0.001);
    TGraphErrors* gr_SM = new TGraphErrors(eta_bin_SM_center.size(), &(eta_bin_SM_center[0]), &(sys_SM[0]), &(eta_bin_SM_error[0]), &(dummy_SM[0]));
    TGraphErrors* gr_FE = new TGraphErrors(eta_bin_FE_center.size(), &(eta_bin_FE_center[0]), &(sys_FE[0]), &(eta_bin_FE_error[0]), &(dummy_FE[0]));
    TGraphErrors* gr_comb = new TGraphErrors(eta_bin_all_center.size(), &(eta_bin_all_center[0]), &(sys_comb[0]), &(eta_bin_all_error[0]), &(dummy[0]));
    double ymax = std::max(*std::max_element(sys_SM.begin(), sys_SM.end()), *std::max_element(sys_FE.begin(), sys_FE.end()));
    double ymin = std::min(*std::min_element(sys_SM.begin(), sys_SM.end()), *std::min_element(sys_FE.begin(), sys_FE.end()));
    ymax = std::max(ymax+0.02,1.2*ymax);
    ymin = std::min(ymin-0.02,1.2*ymin);
    ymax = 0.2; ymin=-0.2;
    TCanvas* canv_sys = tdrCanvas(sysName, eta_bins_all[0]-0.1, eta_bins_all[eta_bins_all.size()-1]+0.1, ymin, ymax, "#eta", "SF_{nominal} - SF_{"+sysName+"} [\%]");
    TLegend *leg_sys = tdrLeg(0.64,0.67,0.79,0.92, 0.040, 42, kBlack);
    tdrHeader(leg_sys,"", 12);
    tdrDraw(gr_comb, "", kFullDotLarge, kGreen-1, kSolid, kGreen-1, 3005, kGreen-1);
    tdrDraw(gr_SM, "", kFullDotLarge, kRed+1, kSolid, kRed+1, 3005, kRed+1);
    tdrDraw(gr_FE, "", kFullDotLarge, kBlue-4, kSolid, kBlue-4, 3005, kBlue-4);
    leg_sys->AddEntry(gr_SM,  "SM","l");
    leg_sys->AddEntry(gr_FE,  "FE","l");
    leg_sys->AddEntry(gr_comb,"comb","l");
    leg_sys->Draw("same");
    canv_sys->Print(path+"standard/"+QCD_DATA+"sys_"+sysName+".pdf","pdf");

  }

  for (unsigned int j = 0; j < systematics_SM_rest.size(); j++) systematics_SM_rest.at(j) = TMath::Sqrt(systematics_SM_rest.at(j));
  for (unsigned int j = 0; j < systematics_FE_rest.size(); j++) systematics_FE_rest.at(j) = TMath::Sqrt(systematics_FE_rest.at(j));

  for (unsigned int i = 0; i < eta_bins_all.size()-1; i++) {
    if (i < shift_FE) systematics_rest[i] = systematics_SM_rest.at(i);
    else if (i < eta_bins_all.size() - 1 - shift_SM) systematics_rest[i] = TMath::Sqrt(TMath::Power(systematics_SM_rest.at(i),2)+TMath::Power(systematics_FE_rest.at(i-shift_FE+shift_barrel+1),2)+TMath::Power((SF_SM.at(i)-SF_FE.at(i-shift_FE+shift_barrel+1))/2, 2));
    else systematics_rest[i] = systematics_FE_rest.at(i-shift_barrel);
  }


  std::vector<double> dummy_SM(eta_bin_SM_center.size(),0.001); std::vector<double> dummy_FE(eta_bin_FE_center.size(),0.001); std::vector<double> dummy(eta_bin_all_center.size(),0.001);
  TGraphErrors* gr_SM_rest = new TGraphErrors(eta_bin_SM_center.size(), &(eta_bin_SM_center[0]), &(systematics_SM_rest[0]), &(eta_bin_SM_error[0]), &(dummy_SM[0]));
  TGraphErrors* gr_FE_rest = new TGraphErrors(eta_bin_FE_center.size(), &(eta_bin_FE_center[0]), &(systematics_FE_rest[0]), &(eta_bin_FE_error[0]), &(dummy_FE[0]));
  TGraphErrors* gr_rest = new TGraphErrors(eta_bin_all_center.size(), &(eta_bin_all_center[0]), &(systematics_rest[0]), &(eta_bin_all_error[0]), &(dummy[0]));
  double ymax = std::max(*std::max_element(systematics_SM_rest.begin(), systematics_SM_rest.end()), *std::max_element(systematics_FE_rest.begin(), systematics_FE_rest.end()));
  double ymin = std::min(*std::min_element(systematics_SM_rest.begin(), systematics_SM_rest.end()), *std::min_element(systematics_FE_rest.begin(), systematics_FE_rest.end()));
  ymax = std::max(ymax+0.02,1.2*ymax);
  ymin = std::min(ymin-0.02,1.2*ymin);
  ymax = 0.2; ymin=-0.2;
  TCanvas* canv_sys_rest = tdrCanvas("rest", eta_bins_all[0]-0.1, eta_bins_all[eta_bins_all.size()-1]+0.1, ymin, ymax, "#eta", "SF_{nominal} - SF_{rest} [\%]");
  TLegend *leg_sys = tdrLeg(0.64,0.67,0.79,0.92, 0.040, 42, kBlack);
  tdrHeader(leg_sys,"", 12);
  tdrDraw(gr_rest, "", kFullDotLarge, kGreen-1, kSolid, kGreen-1, 3005, kGreen-1);
  tdrDraw(gr_SM_rest, "", kFullDotLarge, kRed+1, kSolid, kRed+1, 3005, kRed+1);
  tdrDraw(gr_FE_rest, "", kFullDotLarge, kBlue-4, kSolid, kBlue-4, 3005, kBlue-4);
  leg_sys->AddEntry(gr_SM_rest, "SM","l");
  leg_sys->AddEntry(gr_FE_rest, "FE","l");
  leg_sys->AddEntry(gr_rest,    "comb","l");
  leg_sys->Draw("same");
  canv_sys_rest->Print(path+"standard/"+QCD_DATA+"sys_rest.pdf","pdf");

  return true;

}



void plot_SF_systematics() {

  TString path ;
  TString path_ = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/";

  std::vector<TString> studies;
  studies.push_back("MergeL2Res");

  std::vector<TString> JECs;
  // JECs.push_back("Autumn18_V15");
  // JECs.push_back("Autumn18_V16");
  JECs.push_back("Autumn18_V16h");

  std::vector<TString> JETs;
  JETs.push_back("AK4CHS");

  std::vector<TString> QCDS;
  QCDS.push_back("QCDHT");

  std::vector<TString> DATAS;
  // DATAS.push_back("RunD");
  DATAS.push_back("RunABC");
  // DATAS.push_back("RunABCD");



  for(TString study : studies){
    for(TString JEC : JECs){
      for(TString JET : JETs){
        path = path_+study+"/"+JEC+"/"+JET+"/";
        if (!gSystem->AccessPathName(path)) {
          std::cout << path << '\n';
          for(TString QCD : QCDS){
            for(TString DATA : DATAS){
              TString QCD_DATA = QCD+"/"+DATA+"/";
              std::cout << "start: " << QCD_DATA << '\n';
              plot_SF_systematics_(path, QCD_DATA);
              std::cout << "end: " << QCD_DATA << '\n';
              sleep(5);
            }
          }
        }
      }
    }
  }
}
