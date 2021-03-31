#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "TROOT.h"
#include "TString.h"
#include "TMath.h"
#include "../../../../DiJetJERC/include/constants.h"
#include "../../../../DiJetJERC/include/tdrstyle_all.h"

typedef std::vector<TString> VecS;
typedef std::vector<int> VecI;
typedef std::vector<double> VecD;
typedef std::vector<std::vector<double>> VecDD;

typedef std::map<TString, TGraphErrors*> MapG;
typedef std::map<TString, VecD> MapVD;
typedef std::map<TString, MapVD> MapVDD;

double zero_numerical = 1e-08;
double plotshift = 0.5;

TString lumiRunII = "137";

MapVD CreateTGraphSF(TString filename = "SF_final_tex.txt") {

  filename += "SF_final_DB_UncertaintiesSplit.txt";
  if (gSystem->AccessPathName(filename)) throw runtime_error("check: "+filename);

  MapVD SFs;
  for(const auto& mode: {"eta_center","eta_err","SF","SF_up","SF_stat","SF_err_stat","SF_err_syst"}) {SFs[mode] = {}; SFs[mode].clear();}

  string line;
  ifstream myfile(filename, ios::in);
  while (!myfile.eof()) {
    std::getline(myfile, line); if (line.length()==0) continue; if (line.find("JetEta")!=std::string::npos) continue;
    VecD values; TString tok; int from = 0;
    while (((TString)line).Tokenize(tok, from, " ")) {values.push_back(tok.Atof());}
    if (values[0]<0) continue;
    SFs["eta_center"].push_back((values[1]+values[0])/2);
    SFs["eta_err"].push_back((values[1]-values[0])/2);
    SFs["SF"].push_back(values[5]);
    SFs["SF_up"].push_back(values[7]);
    SFs["SF_stat"].push_back(values[9]);
    SFs["SF_err_stat"].push_back(values[9]-values[5]);
    SFs["SF_err_syst"].push_back(TMath::Sqrt(TMath::Power(values[7]-values[5],2) - TMath::Power(values[9]-values[5],2)));
  }
  myfile.close();
  return SFs;
}


void PlotBinPerYear(TString canvName, VecS years, MapVDD map_SF, double eta_min, double eta_max, double ymin = 0.8, double ymax = 1.5, TString xname = "", TString yname = "JER SF", bool square = false) {
  setTDRStyle();
  int W = (square ? 600 : 800); int H = (square ? 600 : 600);
  int W_ref = (square ? 600 : 800); int H_ref = (square ? 600 : 600);
  float T = (square ? 0.07*H_ref : 0.08*H_ref); float B = (square ? 0.13*H_ref : 0.12*H_ref); float L = (square ? 0.15*W_ref : 0.12*W_ref); float R = (square ? 0.05*W_ref : 0.04*W_ref);
  TCanvas *canv = new TCanvas("BinPerYear"+canvName+eta_min+eta_max,"BinPerYear"+canvName+eta_min+eta_max,50,50,W,H);
  canv->SetFillColor(0); canv->SetBorderMode(0); canv->SetFrameFillStyle(0); canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W ); canv->SetRightMargin( R/W ); canv->SetTopMargin( T/H ); canv->SetBottomMargin( B/H );

  TH1F* hframe = new TH1F("hframe","",years.size()+2,-1,years.size()+1);
  // for (int i = 0; i < years.size(); i++) {TString lab = years[i]; hframe->GetXaxis()->SetBinLabel(i+2, lab.ReplaceAll("_eta",""));};
  for (int i = 0; i < years.size(); i++) {TString lab = years[i]; hframe->GetXaxis()->SetBinLabel(i+2, lab);};

  hframe->SetBit(TH1::kNoStats); hframe->SetBit(kCanDelete);
  hframe->SetMinimum(ymin); hframe->SetMaximum(ymax); hframe->GetYaxis()->SetLimits(ymin,ymax);
  hframe->SetDirectory(0); hframe->Draw(" "); canv->Update();
  hframe->GetYaxis()->SetTitleOffset(square ? 1.25 : 1); hframe->GetXaxis()->SetTitleOffset(square ? 1.0 : 0.9);
  hframe->GetXaxis()->SetTitle(xname); hframe->GetYaxis()->SetTitle(yname);
  hframe->Draw("AXIS"); CMS_lumi(canv,4,11);
  canv->Update(); canv->RedrawAxis(); canv->GetFrame()->Draw();

  TLegend *leg = tdrLeg(0.65,0.67,0.79,0.92, 0.035, 42, kBlack);
  TGraphErrors* gr = new TGraphErrors(years.size());
  TGraphErrors* gr_stat = new TGraphErrors(years.size());

  for (int i = 0; i < years.size(); i++) {
    double val=0; double err=0;  double stat=0; double weights=0;
    for (int bin= 0; bin< map_SF[years[i]]["SF"].size(); bin++){
      double eta_center = map_SF[years[i]]["eta_center"].at(bin);
      double eta_err = map_SF[years[i]]["eta_err"].at(bin);
      if (((eta_center-eta_err)>=(eta_min-zero_numerical)) && ((eta_center+eta_err)<=(eta_max+zero_numerical))) {
        weights += 2*eta_err;
        val   += 2*eta_err*map_SF[years[i]]["SF"].at(bin);
        stat  += 2*eta_err*map_SF[years[i]]["SF_stat"].at(bin);
        err   += 2*eta_err*map_SF[years[i]]["SF_up"].at(bin);
      }
    }
    val  /= weights;
    stat /= weights;
    err  /= weights;
    gr_stat->SetPoint(i, i + 0.5,  val);
    gr->SetPoint(i, i + 0.5,  val);
    gr_stat->SetPointError(i, 0., stat-val);
    gr->SetPointError(i, 0.5, err-val);
  }
  tdrDraw(gr, "P", kFullDotLarge, kBlack, kSolid, kBlack, 0, kBlack);
  tdrDraw(gr_stat, "E", kFullDotLarge, kRed, kSolid, kRed, 0, kRed);
  leg->AddEntry(gr, TString::Format("#eta  #in  %.1f - %.1f", eta_min, eta_max),"p");
  leg->AddEntry(gr_stat, "Statistical","p");
  leg->AddEntry(gr, "Total","p");
  canv->SaveAs("studyTimeDependence/BinPerYear_"+canvName+TString::Format("_range_%.1f_%.1f", eta_min, eta_max)+".pdf");
}

void PlotBinPerYearAll(TString canvName, VecS years, MapVDD map_SF) {
  PlotBinPerYear(canvName, years, map_SF, 0.000, 0.522);
  PlotBinPerYear(canvName, years, map_SF, 0.522, 0.783);
  PlotBinPerYear(canvName, years, map_SF, 0.783, 1.131);
  PlotBinPerYear(canvName, years, map_SF, 1.131, 1.305);
  PlotBinPerYear(canvName, years, map_SF, 1.305, 1.740);
  PlotBinPerYear(canvName, years, map_SF, 1.740, 1.930);
  PlotBinPerYear(canvName, years, map_SF, 1.930, 2.043);
  PlotBinPerYear(canvName, years, map_SF, 2.043, 2.322);
  PlotBinPerYear(canvName, years, map_SF, 2.322, 2.500);

  // { 0.000,        0.522, 0.783,        1.131, 1.305,               1.740, 1.930, 2.043,        2.322, 2.500, 2.650, 2.853, 2.964, 3.139,               5.191};

  PlotBinPerYear(canvName, years, map_SF, 0.000, 1.305);
  PlotBinPerYear(canvName, years, map_SF, 1.305, 2.500);
  PlotBinPerYear(canvName, years, map_SF, 0.000, 2.500);
}

void PlotBand(TString canvName, MapVDD map_SF, VecS years, VecI colors, VecI shapes, bool isShape=true, bool isPoint=true) {
  int nbins = map_SF[years[0]]["eta_center"].size();
  VecD etabins = {0};
  for(int i=0; i<map_SF[years[0]]["eta_center"].size(); i++) etabins.push_back(map_SF[years[0]]["eta_center"].at(i)+ map_SF[years[0]]["eta_err"].at(i));
  double eta_min =  map_SF[years[0]]["eta_center"].at(0)- map_SF[years[0]]["eta_err"].at(0);
  double eta_max =  map_SF[years[0]]["eta_center"].at(nbins-1) + map_SF[years[0]]["eta_err"].at(nbins-1);
  TCanvas* canv_band = tdrCanvas("SF_band"+canvName, eta_min-plotshift, eta_max+plotshift, 0.85, 1.65, "|#eta|", "JER SF");
  TLegend *leg_band = tdrLeg(0.65,0.67,0.79,0.92, 0.035, 42, kBlack);
  tdrHeader(leg_band,"", 12);
  TLine* line_MC = new TLine(eta_min-plotshift, 1, eta_max+plotshift, 1);
  line_MC->SetLineWidth(1); line_MC->SetLineStyle(kDotted); line_MC->SetLineColor(kBlack); line_MC->Draw("same");

  VecD SF_MC_dummy(nbins, 1);
  MapG grs;

  // TH1D* h_unc_band = new TH1D("h_unc_band"+canvName, "h_unc_band"+canvName, etabins.size()-1, &etabins[0]);
  // for(int i=1; i<h_unc_band->GetNbinsX()+1; i++) h_unc_band->SetBinContent(i,0);
  for(int i=0; i<years.size(); i++){
    TString year = years[i];
    TString nameLeg = years[i];
    nameLeg = nameLeg.ReplaceAll("_eta","").ReplaceAll("UL","20").ReplaceAll("preVFP"," early").ReplaceAll("postVFP"," late");
    grs[year+"data"] = new TGraphErrors(nbins, &(map_SF[year]["eta_center"][0]), &map_SF[year]["SF"][0], &(map_SF[year]["eta_err"][0]), &map_SF[year]["SF_err_stat"][0]);
    if (isPoint) tdrDraw(grs[year+"data"], "P5", shapes[i], colors[i], kSolid, colors[i], 1001, colors[i], 0.15);
    if (isPoint) leg_band->AddEntry(grs[year+"data"], nameLeg,"fp");
    // grs[year+"ratio"] = new TGraphErrors(nbins, &(map_SF[year]["eta_center"][0]), &SF_MC_dummy[0], &(map_SF[year]["eta_err"][0]), &map_SF[year+"ratio"]["SF"][0]);
    // grs[year+"MC"] = new TGraphErrors(nbins, &(map_SF[year]["eta_center"][0]), &SF_MC_dummy[0], &(map_SF[year]["eta_err"][0]), &map_SF[year]["SF_err_syst"][0]);
    // if (isShape) tdrDraw(grs[year+"MC"], "E2", shapes[i], colors[i], kDashed, colors[i], 0, colors[i]);
    // if (isShape) tdrDraw(grs[year+"ratio"], "E3", shapes[i], colors[i], kSolid, colors[i], 0, colors[i]);

    TH1D* h_unc_line = new TH1D(year, year, etabins.size()-1, &etabins[0]);
    for(int i=1; i<h_unc_line->GetNbinsX()+1; i++) h_unc_line->SetBinContent(i,1+map_SF[year]["SF_err_syst"][i-1]);
    // for(int i=1; i<h_unc_line->GetNbinsX()+1; i++) h_unc_line->SetBinContent(i,1+map_SF[year+"ratio"]["SF"][i-1]);
    h_unc_line->SetBinContent(0,1);
    h_unc_line->SetBinContent(h_unc_line->GetNbinsX()+1,1);
    if (isShape) tdrDraw(h_unc_line, "hist ][", shapes[i], colors[i], kDashed, colors[i], 0, colors[i]);
    if (isShape) {
      if (! isPoint) leg_band->AddEntry(h_unc_line, nameLeg,"l");
    }

  }

  VecD v_unc_band(nbins,0);
  VecS years_RunII = {"UL16preVFP", "UL16postVFP", "UL17", "UL18"};
  for(int i=0; i<years_RunII.size(); i++){
    TString year = years_RunII[i];
    double lumi = 0; std::string check = year.Data();
    if (check.find("16preVFP")!=std::string::npos) lumi = 19.66;
    if (check.find("16post")!=std::string::npos) lumi = 16.23;
    if (check.find("17")!=std::string::npos) lumi = 41.53;
    if (check.find("18")!=std::string::npos) lumi = 59.74;
    for(int i=0; i<v_unc_band.size(); i++) v_unc_band.at(i) = v_unc_band.at(i)+lumi*(1+map_SF[year]["SF_err_syst"][i]);
  }

  if (isShape && isPoint) {
    TH1D* h_dummy = new TH1D("h_dummyForLeg", "h_dummyForLeg", 1, 0, 1);
    h_dummy->SetLineStyle(kDashed); h_dummy->SetLineColor(kBlack);
    leg_band->AddEntry(h_dummy, "Syst. per year","l");
  }

  std::string temp = canvName.Data();
  TString RunII_err = ((temp.find("eta")!=std::string::npos)?"RunII_eta": "RunII");
  TH1D* h_unc_RunII = new TH1D("h_unc_RunII"+canvName, "h_unc_RunII"+canvName, etabins.size()-1, &etabins[0]);
  for(int i=1; i<h_unc_RunII->GetNbinsX()+1; i++) h_unc_RunII->SetBinContent(i,1+map_SF[RunII_err]["SF_err_syst"].at(i-1));
  h_unc_RunII->SetLineWidth(2);
  if (isShape) tdrDraw(h_unc_RunII, "hist ][", kFullDotSmall, kBlack, kSolid, kBlack, 0, kBlack);
  if (isShape) leg_band->AddEntry(h_unc_RunII, "RunII Syst.","l");

  for(int i=0; i<v_unc_band.size(); i++) v_unc_band.at(i) = v_unc_band.at(i)/lumiRunII.Atof()-1;
  TGraphErrors* gr_unc_band = new TGraphErrors(nbins, &(map_SF[years[0]]["eta_center"][0]), &SF_MC_dummy[0], &(map_SF[years[0]]["eta_err"][0]), &v_unc_band[0]);
  tdrDraw(gr_unc_band, "E2", kFullDotLarge, kBlack, kSolid, kBlack, 1001, kBlack, 0.15);
  leg_band->AddEntry(gr_unc_band, "Systematic Uncertainties","f");


  TString fname = "studyTimeDependence/SF_band"+canvName;
  if (!isShape) fname += "_noerrors";
  if (!isPoint) fname += "_nopoints";

  canv_band->SaveAs(fname+".pdf");
}

void PlotBandAll(TString canvName, MapVDD map_SF, VecS years, VecI colors, VecI shapes) {
  PlotBand(canvName, map_SF, years, colors, shapes, true, true);
  PlotBand(canvName, map_SF, years, colors, shapes, true, false);
  PlotBand(canvName, map_SF, years, colors, shapes, false, true);
}


void studyTimeDependence() {
  TString path_ = std::getenv("CMSSW_BASE"); path_ += "/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/";
  TString outdir = "studyTimeDependence/";
  gPrintViaErrorHandler = kTRUE;
  gErrorIgnoreLevel = kFatal;
  extraText  = "Preliminary";
  lumi_sqrtS = "13 TeV";
  lumi_13TeV = "RunII Legacy, "+lumiRunII+" fb^{-1}";

  MapVDD map_SF;
  map_SF["UL17_DB"] = CreateTGraphSF("/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/Standard/UL17/Summer19UL17_V1_SimpleL1/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_SF["UL18_DB"] = CreateTGraphSF("/nfs/dust/cms/user/amalara//WorkingArea/File//DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/Standard/UL18/Summer19UL18_V2/AK4CHS/standard/QCDHT/RunABCD/");

  map_SF["UL16preVFP"] = CreateTGraphSF(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_SF["UL16postVFP"] = CreateTGraphSF(path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/");
  map_SF["UL17"] = CreateTGraphSF(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_SF["UL18"] = CreateTGraphSF(path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/");
  map_SF["RunII"] = CreateTGraphSF(path_+"eta_JER/Legacy/Summer19Legacy/AK4CHS/standard/QCDHT/RunII/");

  map_SF["UL16preVFP_eta"] = CreateTGraphSF(path_+"eta_simple/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_SF["UL16postVFP_eta"] = CreateTGraphSF(path_+"eta_simple/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/");
  map_SF["UL17_eta"] = CreateTGraphSF(path_+"eta_simple/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_SF["UL18_eta"] = CreateTGraphSF(path_+"eta_simple/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/");
  map_SF["RunII_eta"] = CreateTGraphSF(path_+"eta_simple/Legacy/Summer19Legacy/AK4CHS/standard/QCDHT/RunII/");


  map_SF["UL16preVFP_save"] = CreateTGraphSF(path_+"eta_JER_save/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_SF["UL16postVFP_save"] = CreateTGraphSF(path_+"eta_JER_save/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/");
  map_SF["UL17_save"] = CreateTGraphSF(path_+"eta_JER_save/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_SF["UL18_save"] = CreateTGraphSF(path_+"eta_JER_save/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunD/");
  map_SF["RunII_save"] = CreateTGraphSF(path_+"eta_JER_save/Legacy/Summer19Legacy/AK4CHS/standard/QCDHT/RunII/");



  map_SF["UL16preVFP_B"] = CreateTGraphSF(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunB/");
  map_SF["UL16preVFP_C"] = CreateTGraphSF(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunC/");
  map_SF["UL16preVFP_D"] = CreateTGraphSF(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunD/");
  map_SF["UL16preVFP_E"] = CreateTGraphSF(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunE/");
  map_SF["UL16preVFP_F"] = CreateTGraphSF(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunF/");

  map_SF["UL16preVFP_BCD"] = CreateTGraphSF(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCD/");
  map_SF["UL16preVFP_EF"] = CreateTGraphSF(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunEF/");


  map_SF["UL16postVFP_FG"] = CreateTGraphSF(path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFG/");
  map_SF["UL16postVFP_H"] = CreateTGraphSF(path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunH/");
  map_SF["UL17_B"] = CreateTGraphSF(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunB/");
  map_SF["UL17_C"] = CreateTGraphSF(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunC/");
  map_SF["UL17_D"] = CreateTGraphSF(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunD/");
  map_SF["UL17_E"] = CreateTGraphSF(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunE/");
  map_SF["UL17_F"] = CreateTGraphSF(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunF/");

  map_SF["UL18_A"] = CreateTGraphSF(path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunA/");
  map_SF["UL18_B"] = CreateTGraphSF(path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunB/");
  map_SF["UL18_C"] = CreateTGraphSF(path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunC/");
  map_SF["UL18_D"] = CreateTGraphSF(path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunD/");


  for (auto [name,sf]: map_SF) {
    TString ref = "UL16preVFP";
    std::string check = name.Data();
    if (check.find("_eta")!=std::string::npos) ref+="_eta";
    if (check.find("ratio")!=std::string::npos) continue;

    map_SF[name+"ratio"] = sf;
    for (int i=0; i < sf["SF"].size(); i++) map_SF[name+"ratio"]["SF"].at(i) = fabs(map_SF[name+"ratio"]["SF"].at(i)/map_SF[ref]["SF"].at(i) -1);
  }

  // PlotBinPerYearAll("UL16",  {"UL16preVFP", "UL16postVFP", "UL16preVFP_eta", "UL16postVFP_eta"}, map_SF);
  // PlotBinPerYearAll("UL17",  {"UL17_DB", "UL17_eta", "UL17"}, map_SF);
  // PlotBinPerYearAll("UL18",  {"UL18_DB", "UL18_eta", "UL18"}, map_SF);
  // PlotBinPerYearAll("UL_DB", {"UL16preVFP", "UL16postVFP", "UL17_DB", "UL18_DB"}, map_SF);
  // PlotBinPerYearAll("UL_eta",{"UL16preVFP_eta", "UL16postVFP_eta", "UL17_eta", "UL18_eta"}, map_SF);
  PlotBinPerYearAll("UL",    {"UL16preVFP", "UL16postVFP", "UL17", "UL18", "RunII"}, map_SF);

  // PlotBinPerYearAll("UL_check",    {"UL16preVFP_save", "UL16preVFP", "UL16postVFP_save", "UL16postVFP", "UL17_save", "UL17", "UL18_save", "UL18", "RunII_save", "RunII"}, map_SF);
  // PlotBinPerYearAll("UL_check",    {"UL16preVFP_save","UL16postVFP_save","UL17_save","UL18_save", "UL16preVFP", "UL16postVFP", "UL17", "UL18", "UL16preVFP_eta", "UL16postVFP_eta", "UL17_eta", "UL18_eta"}, map_SF);

  PlotBinPerYearAll("UL_check",    {"UL16preVFP_save", "UL16preVFP", "UL16preVFP_eta", "UL16postVFP_save", "UL16postVFP", "UL16postVFP_eta", "UL17_save", "UL17", "UL17_eta", "UL18_save", "UL18", "UL18_eta"}, map_SF);

  PlotBandAll("",map_SF, {"UL16preVFP", "UL16postVFP", "UL17", "UL18"}, {kOrange+1, kRed+1,kAzure+2,kGreen-1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare});
  PlotBandAll("_eta",map_SF, {"UL16preVFP_eta", "UL16postVFP_eta", "UL17_eta", "UL18_eta"}, {kOrange+1,kRed+1,kAzure+2,kGreen-1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare});



  // PlotBandAll("_UL16preVFP",map_SF, {"UL16preVFP_BCD","UL16preVFP_EF","UL16preVFP"}, {kOrange+1, kRed+1,kAzure+2,kGreen-1, kBlue+1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare, kFullCircle});
  // PlotBandAll("_UL16postVFP",map_SF, {"UL16postVFP_FG","UL16postVFP_H","UL16postVFP"}, {kOrange+1, kRed+1,kAzure+2,kGreen-1, kBlue+1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare, kFullCircle});
  // PlotBandAll("_UL17",map_SF, {"UL17_B","UL17_C","UL17_D","UL17_E","UL17_F"}, {kOrange+1, kRed+1,kAzure+2,kGreen-1, kBlue+1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare, kFullCircle});
  // PlotBandAll("_UL18",map_SF, {"UL18_A","UL18_B","UL18_C","UL18_D","UL18"}, {kOrange+1, kRed+1,kAzure+2,kGreen-1, kBlue+1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare, kFullCircle});

}
