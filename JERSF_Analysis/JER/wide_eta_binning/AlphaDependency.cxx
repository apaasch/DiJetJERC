#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "TROOT.h"
#include "TString.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "../../../include/constants.h"
#include "../../../include/tdrstyle_all.h"
#include "Cosmetics.h"

TString year = "UL18";
TString study = "eta_common_fine_v5_highalpha__tmp";
TString corr = "Summer20UL18_V2";
TString jet = "AK4Puppi";

TString pt = "5";
TString eta = "3";

TString path = "/nfs/dust/cms/user/paaschal/UHH2_DiJet/CMSSW_10_6_28/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/";

vector<double> bins_eta = { 0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.566, 1.740, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 5.191};
vector<double> bins_pt = { 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120, 1240, 1400, 1500, 1600, 2000, 2500, 7000};

map<TString, TF1*> functions = {
  {"log", new TF1("log", "[0]+[1]*log(x)", 0, 1.05)},
  {"pol2", new TF1("pol2", "[0]+[1]*x+[2]*x*x", 0, 1.05)},
  // {"exp", new TF1("exp", "[0]+[1]*Exp([2]*x)", 0, 1.05)},
  {"lin03", new TF1("lin035", "[0]+[1]*x", 0, 0.35)},
  {"lin04", new TF1("lin45", "[0]+[1]*x", 0, 0.45)},
};

map<TString, int> f_colors = {
  {"log", kYellow+2},
  {"pol2", kBlue+2},
  {"exp", kAzure-3},
  {"lin03", kRed+2},
  {"lin04", kGreen+2},
};

map<TString, TString> f_legs = {
  {"log", "logarithmic"},
  {"pol2", "polynomial 2"},
  {"exp", "exponential"},
  {"lin03", "linear from #alpha #in (0,0.3)"},
  {"lin04", "linear from #alpha #in (0,0.4)"},
};

void FitHistogram(TH1F *h_, TF1* f_, int c_) {
  h_->Fit(f_, "QR");
  f_->SetLineColor(c_);
  f_->SetLineStyle(kSolid);
  f_->SetLineWidth(2);
}

TGraph *GetRatios(TH1F *h, TF1 *f){
  TGraph *gr = new TGraph();
  int i = 0;
  for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
    if (h->GetBinContent(bin) == 0) continue;
    double x = h->GetBinCenter(bin);
    double y = f->Eval(x) / h->GetBinContent(bin);
    gr->SetPoint(i++, x, y);
    // cout << setw(4) << i << setw(15) << x << setw(15) << y << setw(15) << f->Eval(x) << setw(15) << h->GetBinContent(bin) << endl;
  }
  gr->GetYaxis()->SetRangeUser(0.89, 1.11);
  gr->SetMarkerStyle(kFullCircle);
  gr->SetLineColor(f->GetLineColor());
  gr->SetLineStyle(f->GetLineStyle());
  gr->SetLineWidth(f->GetLineWidth());
  return gr;
}

void AlphaDependency() {

  TFile *file = new TFile(path+"file/"+study+"/"+year+"/"+corr+"/"+jet+"/standard/QCDHT/RunABCD/output/widths.root");

  TH1F *h_mc = (TH1F*) file->Get("widths_fe_eta"+eta+"_pt"+pt);
  h_mc->SetMarkerStyle(kFullCircle);

  map<TString, TF1*> f_mc;
  map<TString, TGraph*> r_mc;
  for(auto& f: functions){
    TString f_name = f.first;
    f_mc[f_name] = (TF1*) f.second->Clone();
    FitHistogram(h_mc, f_mc[f_name], f_colors[f_name]);
    r_mc[f_name] = GetRatios(h_mc, f_mc[f_name]);
  }

  TString nameXaxis = "#alpha_{max}";
  TString nameYaxis = "#sigma_{A}";
  extraText3.clear();
  extraText3.push_back(Form("%d GeV < p_{T}^{ave} < %d GeV", (int)bins_pt[pt.Atoi()], (int)bins_pt[pt.Atoi()+1]));
  extraText3.push_back(Form("%.1f < |#eta| < %.1f", bins_eta[eta.Atoi()], bins_eta[eta.Atoi()+1]));
  lumi_13TeV = "[MC 106X] "+year;
  setTDRStyle();
  TCanvas *canv = tdrDiCanvas2("alpha", 0, 1.05, 0.01, 0.31, 0.89, 1.11, nameXaxis, nameYaxis, "#frac{fit}{mc}", kRectangular, 0);
  canv->SetGrid();
  canv->cd(1);
  h_mc->Draw("same");
  TLegend *leg = tdrLeg(0.5,0.6,0.9,0.9);
  for(auto&f:functions){
    leg->AddEntry(f_mc[f.first],f_legs[f.first],"l");
    f_mc[f.first]->Draw("l same");
   }
  leg->Draw();
  canv->cd(2);
  tdrDraw(new TLine(0, 1, 1.05, 1));
  for(auto&r:r_mc) r.second->Draw("l same");
  canv->Print(path+"plots/alpha/alpha_eta"+eta+"_pt"+pt+".pdf");

  return;
}
