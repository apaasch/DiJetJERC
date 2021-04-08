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
bool debug = false;

// ---------------------------------------------------------------------------------
MapVD GetFitParameters(TString path_, TString filename_) {

  TString filename = path_+"pdfy/kValues/"+filename_+".txt";
  if (gSystem->AccessPathName(filename)) throw runtime_error("check: "+filename);

  MapVD parameters;
  for(const auto& mode: {"eta","eta_err","N","Nerr","S","Serr","C","Cerr","kNS","kNSerr","kC","kCerr"}) {parameters[mode] = {}; parameters[mode].clear();}

  string line;
  ifstream myfile(filename, ios::in);
  while (!myfile.eof()) {
    std::getline(myfile, line);
    if(debug) cout << line << endl;
    if (line.length()==0) continue;
    if (line.find("kNS")!=std::string::npos) continue; // skip first line
    VecD values; TString tok; int from = 0;
    while (((TString)line).Tokenize(tok, from, " ")) {values.push_back(tok.Atof());} // cuts string
    parameters["eta"].push_back((values[1]+values[0])/2);
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

  MapVDD studies;
  for(map<TString,MapVD>::iterator year = parameters.begin(); year != parameters.end(); ++year)
  {
    TString key_year = year->first;
    VecD kNS    = year->second["kNS"];
    VecD kNSerr = year->second["kNSerr"];
    VecD kC     = year->second["kC"];
    VecD kCerr  = year->second["kCerr"];
    for(unsigned int value=0; value<kNS.size(); value++){
      // Store again for PlotValue()
      studies[key_year]["eta"].push_back(parameters[key_year]["eta"][value]);
      studies[key_year]["eta_err"].push_back(parameters[key_year]["eta_err"][value]);

      studies[key_year]["Ratio"].push_back(kNS[value]/kC[value]);
      studies[key_year]["Ratioerr"].push_back(sqrt(pow(kNSerr[value]/kNS[value],2)+pow(kCerr[value]/kC[value],2)));

      studies[key_year]["Diff"].push_back(kNS[value]-kC[value]);
      studies[key_year]["Differr"].push_back(sqrt(pow(kNSerr[value],2)+pow(kCerr[value],2)));
    }
  }
  return studies;
}


// ---------------------------------------------------------------------------------
// Purpose: Debugging
void printVector(VecD vector) {
  cout << "\nSize: " << vector.size() << endl;
  for(double value: vector) cout << value << endl;
  cout << endl;
}

// ---------------------------------------------------------------------------------
void PlotValue(TString study, TString year, TString method, MapTS save, MapVDD map, double ymin = 0.5, double ymax = 1.5, double eta_max = 5.2) {

  TGraphErrors *gr = new TGraphErrors(map[year]["eta"].size(), &(map[year]["eta"][0]), &(map[year][study][0]), &(map[year]["eta_err"][0]), &(map[year][study+"err"][0]));
  TCanvas* canv = tdrCanvas(study+year, 0, eta_max, ymin, ymax, "#eta", study);
  canv->SetTickx(0);
  canv->SetTicky(0);
  tdrDraw(gr, "P", kFullDotLarge);
  canv->Print(save[year]+"/"+study+"_"+method+".pdf","pdf");
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

  if(debug) cout << "Getting Parameters for correlated FE ..." << endl;
  MapVDD map_FIT_FE_corr;
  map_FIT_FE_corr["UL16preVFP"]  = GetFitParameters(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/", "correlated_FE");
  map_FIT_FE_corr["UL16postVFP"] = GetFitParameters(path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/", "correlated_FE");
  map_FIT_FE_corr["UL17"]        = GetFitParameters(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/", "correlated_FE");
  map_FIT_FE_corr["UL18"]        = GetFitParameters(path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/", "correlated_FE");

  if(debug) cout << "Getting Parameters for correlated SM ..." << endl;
  MapVDD map_FIT_SM_corr;
  map_FIT_SM_corr["UL16preVFP"]  = GetFitParameters(path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/", "correlated_SM");
  map_FIT_SM_corr["UL16postVFP"] = GetFitParameters(path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/", "correlated_SM");
  map_FIT_SM_corr["UL17"]        = GetFitParameters(path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/", "correlated_SM");
  map_FIT_SM_corr["UL18"]        = GetFitParameters(path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/", "correlated_SM");

  MapVDD map_PtStudies_FE_corr = GetPtStudies(map_FIT_FE_corr);
  MapVDD map_PtStudies_SM_corr = GetPtStudies(map_FIT_SM_corr);

  MapTS map_SAVE;
  map_SAVE["UL16preVFP"]  = path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/pdfy/kValues";
  map_SAVE["UL16postVFP"] = path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/pdfy/kValues";
  map_SAVE["UL17"]        = path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/pdfy/kValues";
  map_SAVE["UL18"]        = path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/pdfy/kValues";

  cout << "Starting with UL16preVFP ... " << endl;
  PlotValue("Ratio", "UL16preVFP", "correlated_FE", map_SAVE, map_PtStudies_FE_corr,   0.25,  1.75);
  PlotValue("N",     "UL16preVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,       -10.00, 10.50);
  PlotValue("S",     "UL16preVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,         0.50,  1.50);
  PlotValue("C",     "UL16preVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,        -0.05,  0.16);
  PlotValue("kNS",   "UL16preVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,         1.00,  1.50);
  PlotValue("kC",    "UL16preVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,         0.70,  2.00);

  PlotValue("Ratio", "UL16preVFP", "correlated_SM", map_SAVE, map_PtStudies_SM_corr,   0.25,  1.75, 2.322);
  PlotValue("N",     "UL16preVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,       -10.00, 10.50, 2.322);
  PlotValue("S",     "UL16preVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,         0.50,  1.50, 2.322);
  PlotValue("C",     "UL16preVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,        -0.05,  0.16, 2.322);
  PlotValue("kNS",   "UL16preVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,         1.00,  1.50, 2.322);
  PlotValue("kC",    "UL16preVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,         0.70,  2.00, 2.322);

  cout << "Starting with UL16postVFP ... " << endl;
  PlotValue("Ratio", "UL16postVFP", "correlated_FE", map_SAVE, map_PtStudies_FE_corr,   0.25,  1.75);
  PlotValue("N",     "UL16postVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,       -10.00, 10.50);
  PlotValue("S",     "UL16postVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,         0.50,  1.50);
  PlotValue("C",     "UL16postVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,        -0.05,  0.16);
  PlotValue("kNS",   "UL16postVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,         1.00,  1.50);
  PlotValue("kC",    "UL16postVFP", "correlated_FE", map_SAVE, map_FIT_FE_corr,         0.70,  2.00);

  PlotValue("Ratio", "UL16postVFP", "correlated_SM", map_SAVE, map_PtStudies_SM_corr,   0.25,  1.75, 2.322);
  PlotValue("N",     "UL16postVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,       -10.00, 10.50, 2.322);
  PlotValue("S",     "UL16postVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,         0.50,  1.50, 2.322);
  PlotValue("C",     "UL16postVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,        -0.05,  0.16, 2.322);
  PlotValue("kNS",   "UL16postVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,         1.00,  1.50, 2.322);
  PlotValue("kC",    "UL16postVFP", "correlated_SM", map_SAVE, map_FIT_SM_corr,         0.70,  2.00, 2.322);

  cout << "Starting with UL17 ... " << endl;
  PlotValue("Ratio", "UL17", "correlated_FE", map_SAVE, map_PtStudies_FE_corr,   0.25,  1.75);
  PlotValue("N",     "UL17", "correlated_FE", map_SAVE, map_FIT_FE_corr,       -10.00, 10.50);
  PlotValue("S",     "UL17", "correlated_FE", map_SAVE, map_FIT_FE_corr,         0.50,  1.50);
  PlotValue("C",     "UL17", "correlated_FE", map_SAVE, map_FIT_FE_corr,        -0.05,  0.16);
  PlotValue("kNS",   "UL17", "correlated_FE", map_SAVE, map_FIT_FE_corr,         1.00,  1.50);
  PlotValue("kC",    "UL17", "correlated_FE", map_SAVE, map_FIT_FE_corr,         0.70,  2.00);

  PlotValue("Ratio", "UL17", "correlated_SM", map_SAVE, map_PtStudies_SM_corr,   0.25,  1.75, 2.322);
  PlotValue("N",     "UL17", "correlated_SM", map_SAVE, map_FIT_SM_corr,       -10.00, 10.50, 2.322);
  PlotValue("S",     "UL17", "correlated_SM", map_SAVE, map_FIT_SM_corr,         0.50,  1.50, 2.322);
  PlotValue("C",     "UL17", "correlated_SM", map_SAVE, map_FIT_SM_corr,        -0.05,  0.16, 2.322);
  PlotValue("kNS",   "UL17", "correlated_SM", map_SAVE, map_FIT_SM_corr,         1.00,  1.50, 2.322);
  PlotValue("kC",    "UL17", "correlated_SM", map_SAVE, map_FIT_SM_corr,         0.70,  2.00, 2.322);

  cout << "Starting with UL18 ... " << endl;
  PlotValue("Ratio", "UL18", "correlated_FE", map_SAVE, map_PtStudies_FE_corr,   0.25,  1.75);
  PlotValue("N",     "UL18", "correlated_FE", map_SAVE, map_FIT_FE_corr,       -10.00, 10.50);
  PlotValue("S",     "UL18", "correlated_FE", map_SAVE, map_FIT_FE_corr,         0.50,  1.50);
  PlotValue("C",     "UL18", "correlated_FE", map_SAVE, map_FIT_FE_corr,        -0.05,  0.16);
  PlotValue("kNS",   "UL18", "correlated_FE", map_SAVE, map_FIT_FE_corr,         1.00,  1.50);
  PlotValue("kC",    "UL18", "correlated_FE", map_SAVE, map_FIT_FE_corr,         0.70,  2.00);

  PlotValue("Ratio", "UL18", "correlated_SM", map_SAVE, map_PtStudies_SM_corr,   0.25,  1.75, 2.322);
  PlotValue("N",     "UL18", "correlated_SM", map_SAVE, map_FIT_SM_corr,       -10.00, 10.50, 2.322);
  PlotValue("S",     "UL18", "correlated_SM", map_SAVE, map_FIT_SM_corr,         0.50,  1.50, 2.322);
  PlotValue("C",     "UL18", "correlated_SM", map_SAVE, map_FIT_SM_corr,        -0.05,  0.16, 2.322);
  PlotValue("kNS",   "UL18", "correlated_SM", map_SAVE, map_FIT_SM_corr,         1.00,  1.50, 2.322);
  PlotValue("kC",    "UL18", "correlated_SM", map_SAVE, map_FIT_SM_corr,         0.70,  2.00, 2.322);

}
