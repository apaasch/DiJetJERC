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

typedef std::map<TString, double> MapD;
typedef std::map<TString, MapD> MapDD;
typedef std::map<TString, TGraphErrors*> MapG;
typedef std::map<TString, VecD> MapVD;
typedef std::map<TString, MapVD> MapVDD;
typedef std::map<TString, TString> MapTS;
typedef std::map<TString, MapTS> MapDTS;

TString lumiRunII = "137";

// ---------------------------------------------------------------------------------
MapVD GetFitParameters(TString path_, TString filename_) {

  TString filename = path_+"pdfy/kValues/"+filename_+".txt";
  if (gSystem->AccessPathName(filename)) throw runtime_error("check: "+filename);

  MapVD parameters;
  for(const auto& mode: {"eta_center","eta_err","N","Nerr","S","Serr","C","Cerr","kNS","kNSerr","kC","kCerr"}) {parameters[mode] = {}; parameters[mode].clear();}

  string line;
  ifstream myfile(filename, ios::in);
  while (!myfile.eof()) {
    std::getline(myfile, line);
    if (line.length()==0) continue;
    if (line.find("kNS")!=std::string::npos) continue; // skip first line
    VecD values; TString tok; int from = 0;
    while (((TString)line).Tokenize(tok, from, " ")) {values.push_back(tok.Atof());} // cuts string
    parameters["eta_center"].push_back((values[1]+values[0])/2);
    parameters["eta_err"].push_back((values[1]-values[0])/2);
    parameters["N"].push_back(values[2]);
    parameters["Nerr"].push_back(values[3]);
    parameters["S"].push_back(values[4]);
    parameters["Serr"].push_back(values[5]);
    parameters["C"].push_back(values[6]);
    parameters["Cerr"].push_back(values[7]);
    parameters["kNS"].push_back(values[8]);
    parameters["kNSerr"].push_back(values[9]);
    parameters["kC"].push_back(values[10]);
    parameters["kCerr"].push_back(values[11]);
  }
  myfile.close();
  return parameters;
}

// ---------------------------------------------------------------------------------
MapVDD GetPtStudies(MapVDD parameters) {

  MapVD  single_studies;
  MapVDD all_studies;

  for(map<TString,MapVD>::iterator year = parameters.begin(); year != parameters.end(); ++year)
  {
    TString key_year = year->first;
    VecD vNS    = year->second["kNS"];
    VecD vNSerr = year->second["kNSerr"];
    VecD vC     = year->second["kC"];
    VecD vCerr  = year->second["kCerr"];

    for(unsigned int value=0; value<vNS.size(); value++){
      single_studies["Ratio"].push_back(vNS[value]/vC[value]);
      single_studies["RatioErr"].push_back(sqrt(pow(vNSerr[value]/vNS[value],2)+pow(vCerr[value]/vC[value],2)));

      single_studies["Difference"].push_back(vNS[value]-vC[value]);
      single_studies["DifferenceErr"].push_back(sqrt(pow(vNSerr[value],2)+pow(vCerr[value],2)));
    }

    all_studies[key_year] = single_studies;
  }

  return all_studies;
}

// ---------------------------------------------------------------------------------
void PlotValue(TString yname, TString year, TString method, TString corr, TString save, VecD eta, VecD eta_err, VecD values, VecD values_err, double eta_max = 5.2, double ymin = 0.5, double ymax = 1.5) {
  setTDRStyle();
  // int W = (square ? 600 : 800); int H = (square ? 600 : 600);
  // int W = 600; int H = 600;
  // TCanvas *canv = new TCanvas(canvName,canvName,50,50,W,H);
  // // float T = (square ? 0.07*H_ref : 0.08*H_ref); float B = (square ? 0.13*H_ref : 0.12*H_ref); float L = (square ? 0.15*W_ref : 0.12*W_ref); float R = (square ? 0.05*W_ref : 0.04*W_ref);
  // float T = 0.08*H_ref; float B = 0.12*H_ref; float L = 0.12*W_ref; float R = 0.04*W_ref;
  // canv->SetFillColor(0); canv->SetBorderMode(0); canv->SetFrameFillStyle(0); canv->SetFrameBorderMode(0);
  // canv->SetLeftMargin( L/W ); canv->SetRightMargin( R/W ); canv->SetTopMargin( T/H ); canv->SetBottomMargin( B/H );

  TGraphErrors *gr = new TGraphErrors(eta.size(), &(eta[0]), &(values[0]), &(eta_err[0]), &(values_err[0]));

  TCanvas* canv = tdrCanvas(yname+year, 0, eta_max, ymin, ymax, "#eta", yname);
  canv->SetTickx(0);
  canv->SetTicky(0);
  tdrDraw(gr, "P", kFullDotLarge);
  canv -> Print(save+"/"+yname+"_"+corr+"_"+method+".pdf","pdf");
}

// ---------------------------------------------------------------------------------

// ------------------------------
//          MAIN PROGRAM
// ------------------------------

void studyFitParameters(){
  TString path_ = std::getenv("CMSSW_BASE"); path_ += "/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/";
  TString outdir = "studyTimeDependence/";
  gPrintViaErrorHandler = kTRUE;
  gErrorIgnoreLevel = kFatal;
  extraText  = "Preliminary";
  lumi_sqrtS = "13 TeV";
  lumi_13TeV = "RunII Legacy, "+lumiRunII+" fb^{-1}";

  MapVDD map_FIT;
  map_FIT["UL16preVFP"]  = GetFitParameters(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/", "correlated_FE");
  map_FIT["UL16postVFP"] = GetFitParameters(path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/", "correlated_FE");
  map_FIT["UL17"]        = GetFitParameters(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/", "correlated_FE");
  map_FIT["UL18"]        = GetFitParameters(path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/", "correlated_FE");

  MapVDD map_PtStudies = GetPtStudies(map_FIT);

  MapTS map_SAVE;
  map_SAVE["UL16preVFP"]  = path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/pdfy/kValues";
  map_SAVE["UL16postVFP"] = path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/pdfy/kValues";
  map_SAVE["UL17"]        = path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/pdfy/kValues";
  map_SAVE["UL18"]        = path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/pdfy/kValues";

  PlotValue("Ratio", "UL18", "FE", "correlated", map_SAVE["UL18"], map_FIT["UL18"]["eta_center"], map_FIT["UL18"]["eta_err"], map_PtStudies["UL18"]["Ratio"], map_PtStudies["UL18"]["RatioErr"]);

}
