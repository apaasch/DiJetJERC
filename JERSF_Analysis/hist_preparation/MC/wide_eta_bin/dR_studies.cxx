
#include <iostream>
#include <algorithm>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TSystem.h>
#include "MySelector.h"
#include "MyJet.h"
#include "../../../../include/contants.hpp"


void dR_plots(TString add = "wide_eta_bin", TString root_filename = "histograms_mc_incl_full_2D", TString histo_name = "dR_forward_probe10_pt2_alpha6"){
  TString filename = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/JER2017/Analysis/hist_preparation/MC/";
  filename+=add;
  filename += "/file/Single/PtBinned_full/";
  filename += root_filename;
  filename += ".root";
  TFile *file = new TFile( filename, "READ");
  std::cout << add << std::endl;
  std::cout << filename << std::endl;
  std::cout << histo_name << std::endl;
  TH2F *h = (TH2F*)file->Get(histo_name);
  // h->RebinX(4);
  // h->RebinY(4);
  TCanvas *c = new TCanvas("c","c",50, 50, 800, 600);
  h->Draw("colz");
  TString outname = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/JER2017/Analysis/hist_preparation/MC/wide_eta_bin/dRstudies/";
  c->Print(outname+histo_name+"_"+add+".pdf", "pdf");
  c->Print(outname+histo_name+"_"+add+".png", "png");
  delete c;
  delete h;
  delete file;
}



void check_asymmetry(TString histo_name, int &count1){
  TString add = "wide_eta_bin";
  TString root_filename = "histograms_mc_incl_full";
  TString filename = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_94/CMSSW_9_4_1/src/UHH2/JER2017/Analysis/hist_preparation/MC/";
  filename+=add;
  filename += "/file/Single/PtBinned_full/";
  filename += root_filename;
  filename += ".root";
  TFile *file = new TFile( filename, "READ");
  TH2F *h = (TH2F*)file->Get(histo_name);
  // h->RebinX(h->GetNbinsX()/2);
  h->RebinX(4);
  if (TMath::Abs(h->GetBinContent(h->FindBin(0))-h->GetBinContent(h->FindBin(0)-1))> TMath::Abs(h->GetBinError(h->FindBin(0))+h->GetBinError(h->FindBin(0)-1))) {
    count1++;
    std::cout << histo_name << std::endl;
    TCanvas *c = new TCanvas("c", "c", 50,50,800,600);
    h->Draw();
    c->Print("dRstudies/"+histo_name+".png", "png");
    delete c;
  }
  delete h;
  delete file;
}



TCanvas* fitSlices_inFile(TString filename = "file/Single/PtBinned_full/histograms_mc_incl_full_2D_dR.root",
TString histo_name = "asy_dR_barrel_forward_probe10_pt2_alpha6_dR_probe3") {
  TFile *file = new TFile( filename, "READ");
  TH2F *h2 = (TH2F*)file->Get(histo_name);
  h2->RebinX(4);
  h2->RebinY(4);
  // TObjArray *a = histo_name.Tokenize("_")
  // for (size_t i = 0; i < a->GetEntries(); i++) {
  //   /* code */
  // }
  // h2->SetTitle("test2")
  // Create a canvas and divide it
  TCanvas *c1 = new TCanvas("c1","c1",50, 50, 800, 600);
  // c1->SetFillColor(42);
  c1->Divide(2,1);
  TPad *leftPad = (TPad*)c1->cd(1);
  leftPad->Divide(1,2);
  // Draw 2-d original histogram
  leftPad->cd(1);
  gPad->SetTopMargin(0.12);
  // gPad->SetFillColor(33);
  h2->Draw("colz");
  h2->GetXaxis()->SetLabelSize(0.06);
  h2->GetYaxis()->SetLabelSize(0.06);
  // h2->SetMarkerColor(kYellow);
  // Fit slices projected along Y fron bins in X [7,32] with more than 20 bins  in Y filled
  // h2->FitSlicesX(0, 0, -1, 2);
  h2->FitSlicesX();
  // Show fitted "mean" for each slice
  leftPad->cd(2);
  // gPad->SetFillColor(33);
  TH2F *h2_0 = (TH2F*)file->Get(histo_name+"_0");
  h2_0->Draw();
  TPad *rightPad = (TPad*)c1->cd(2);
  rightPad->Divide(1,2);
  rightPad->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  // gPad->SetFillColor(33);
  TH2F *h2_1 = (TH2F*)file->Get(histo_name+"_1");
  h2_1->GetYaxis()->SetRangeUser(-0.3, 0.3);
  h2_1->Draw();
  // Show fitted "sigma" for each slice
  rightPad->cd(2);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  // gPad->SetFillColor(33);
  TH2F *h2_2 = (TH2F*)file->Get(histo_name+"_2");
  h2_2->GetYaxis()->SetRangeUser(-0.1, 0.4);
  h2_2->Draw();
  //attributes
  h2_0->SetLineColor(2);
  h2_1->SetLineColor(2);
  h2_2->SetLineColor(2);
  h2_0->SetMarkerColor(2);
  h2_1->SetMarkerColor(2);
  h2_2->SetMarkerColor(2);
  h2_0->SetMarkerStyle(21);
  h2_1->SetMarkerStyle(21);
  h2_2->SetMarkerStyle(21);
  h2_0->SetMarkerSize(0.6);
  h2_1->SetMarkerSize(0.6);
  h2_2->SetMarkerSize(0.6);
  return c1;
}

void fitSlices(){
  double dRbins [] = { 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0};
  std::vector<double> dR_bins(dRbins, dRbins + sizeof(dRbins)/sizeof(double));
  int min = 0;
  int max = dR_bins.size()-1;
  bool debug = false;


  TString filename = "file/Single/PtBinned_full/histograms_mc_incl_full_2D_dR.root";
  TString histo_name = "asy_dR_barrel_forward_probe10_pt2_alpha6_dR_probe";
  for (unsigned int i = min; i < max; i++) {
    histo_name = "asy_dR_barrel_forward_probe10_pt2_alpha6_dR_probe"; histo_name += i;
    TCanvas *c = fitSlices_inFile(filename, histo_name);
    TString outname = "dRstudies/asy_dR_barrel"; outname += i;
    c->Print(outname+".png", "png");
    c->Print(outname+".pdf", "pdf");
    delete c;
  }

  for (unsigned int i = min; i < max; i++) {
    histo_name = "asy_dR_probe_forward_probe10_pt2_alpha6_dR_barrel"; histo_name += i;
    TCanvas *c = fitSlices_inFile(filename, histo_name);
    TString outname = "dRstudies/asy_dR_probe"; outname += i;
    c->Print(outname+".png", "png");
    c->Print(outname+".pdf", "pdf");
    delete c;
  }

  if (debug){
    TString filename = "file/Single/PtBinned_full/histograms_mc_incl_full_2D_dR.root";
    TString histo_name = "asy_dR_barrel_forward_probe10_pt2_alpha6_dR_probe";
    int bin_min, bin_max;
    TH1D *h2_test = 0;
    filename = "file/Single/PtBinned_full/histograms_mc_incl_full_2D_dR.root";
    histo_name = "asy_dR_barrel_forward_probe10_pt2_alpha6_dR_probe"; histo_name += 2;
    TFile *file = new TFile( filename, "READ");
    TH2F *h2 = (TH2F*)file->Get(histo_name);
    h2->RebinX(4);
    h2->RebinY(4);
    new TCanvas;
    h2->Draw("colz");
    bin_min = h2->GetYaxis()->FindBin(2.4);
    bin_max = h2->GetYaxis()->FindBin(2.8);
    h2->ProjectionX("h2_test", bin_min+1, bin_max+1, "e");
    h2_test = (TH1D*)gDirectory->Get("h2_test");
    new TCanvas;
    h2_test->Draw();

    for (int i = 0; i < 16; i++) {
      TString name = i;
      h2->ProjectionX(name, i+1, i+2, "e");
      h2_test = (TH1D*)gDirectory->Get(name);
      new TCanvas;
      h2_test->Draw();
    }
  }
}



void dR_studies() {
  int count = 0;
  int count1 = 0;
  for (size_t i = 10; i < 13; i++) {
    for (size_t j = 2; j < 10; j++) {
      for (size_t k = 3; k < 7; k++) {
        TString histo_name = "forward_probe";
        histo_name += i;
        histo_name += "_pt";
        histo_name += j;
        histo_name += "_alpha";
        histo_name += k;
        count ++;
        check_asymmetry(histo_name, count1);
      }
    }
  }
  std::cout << count1 << " out of " << count << std::endl;


  count = 0;
  count1 = 0;
  for (size_t i = 1; i < 11; i++) {
    for (size_t j = 4; j < 10; j++) {
      for (size_t k = 3; k < 7; k++) {
        TString histo_name = "asymm_eta";
        histo_name += i;
        histo_name += "_pt";
        histo_name += j;
        histo_name += "_alpha";
        histo_name += k;
        count ++;
        // check_asymmetry(histo_name, count1);
      }
    }
  }
  std::cout << count1 << " out of " << count << std::endl;



  count = 0;
  count1 = 0;
  for (size_t i = 2; i < 10; i++) {
    for (size_t j = 4; j < 10; j++) {
      for (size_t k = 3; k < 7; k++) {
        TString histo_name = "forward_control_probe";
        histo_name += i;
        histo_name += "_pt";
        histo_name += j;
        histo_name += "_alpha";
        histo_name += k;
        count ++;
        // check_asymmetry(histo_name, count1);
      }
    }
  }
  std::cout << count1 << " out of " << count << std::endl;

  // dR_plots("wide_eta_bin", "histograms_mc_incl_full_2D", "dR_forward_probe10_pt2_alpha6");
  // dR_plots("wide_eta_bin", "histograms_mc_incl_full", "forward_probe10_pt2_alpha6");
  // for (int i = 1; i < 16; i++) {
  //   TString add = "wide_eta_bin_test"; add +=i;
  //   dR_plots(add, "histograms_mc_incl_full_control", "dR_forward_probe10_pt2_alpha6");
  //   dR_plots(add, "histograms_mc_incl_full", "forward_probe10_pt2_alpha6");
  // }

  // fitSlices();

}
