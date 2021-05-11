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

typedef std::map<TString, VecD> MapVD;
typedef std::map<TString, MapVD> MapVDD;

typedef std::map<TString, TString> MapTS;
typedef std::map<TString, MapTS> MapDTS;

typedef std::map<TString, TH1F*> MapH;
typedef std::map<TString, MapH> MapMH;
typedef std::map<TString, TGraphErrors*> MapG;

TString lumiRunII = "137";
VecS samples = {"standard", "alpha", "gaustails_0.95", "JEC", "JER", "PLI", "PU"};
VecS syst    = {"","/nominal","/up","/down"};

bool debug      = false;
bool same_range = false;

////////////////////////////////////////////////////////////////////////////
//    General Functions                                                   //
////////////////////////////////////////////////////////////////////////////
// =======================================================================================================
int Round(double wert)
{
  return (int) (wert + ((wert < 0)? - 0.5 : 0.5));
}

// =======================================================================================================
TString dtos(double number, int precision)
{
  stringstream stream;
  stream << std::fixed << std::setprecision(precision) << number;
  return stream.str();
}

// =======================================================================================================
void print_bin_content(MapVD content, TString name, bool isUnc=false){
  // function for debugging purpose only
  TString con = "content";
  if(isUnc) con = "uncertainty";
  cout << setprecision(5) << endl;
  cout << "==============================================================" << endl;
  cout << " ------ "+name+"; Entries: " << content["eta"].size() << " " << content["pT"].size() << " " << content[con].size() << endl;
  cout << "==============================================================" << endl;
  cout << right << setw(8) << "eta" << setw(8) << "pt" << setw(13) << con << "\n";
  cout << "--------------------------------------------------------------" << endl;
  for(unsigned int i=0; i<content["eta"].size(); i++) // size of binEOY and binsUL equal
  {

    cout << setw(8) << content["eta"][i] << setw(8) << content["pT"][i] << setw(13) << content[con][i] << "\n";
  }
  cout << "==============================================================" << endl;
  cout << left << endl;
}

// =======================================================================================================
// Purpose: Debugging
void printVector(VecD vector) {
  cout << "\nSize: " << vector.size() << endl;
  for(double value: vector) cout << value << ", ";
  cout << endl;
}

////////////////////////////////////////////////////////////////////////////
//    Fit Studies                                                         //
////////////////////////////////////////////////////////////////////////////
// =======================================================================================================
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


// =======================================================================================================
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

// =======================================================================================================
void PlotValue(TString study, TString year, TString method, MapTS save, MapVDD map, double ymin = 0.5, double ymax = 1.5, double eta_max = 5.2) {
  TGraphErrors *gr = new TGraphErrors(map[year]["eta"].size(), &(map[year]["eta"][0]), &(map[year][study][0]), &(map[year]["eta_err"][0]), &(map[year][study+"err"][0]));
  TCanvas* canv = tdrCanvas(study+year, 0, eta_max, ymin, ymax, "#eta", study);
  canv->SetTickx(0);
  canv->SetTicky(0);
  tdrDraw(gr, "P", kFullDotLarge);
  canv->Print(save[year]+"/pdfy/kValues/"+study+"_"+method+".pdf","pdf");
}

////////////////////////////////////////////////////////////////////////////
//    Uncertainty studies                                                 //
////////////////////////////////////////////////////////////////////////////
// =======================================================================================================
bool skipSample(TString sample, TString sys){
  bool skip = false;
  if((sample.EqualTo("standard")     && (!sys.EqualTo("")))                           ||
  (sample.EqualTo("alpha")           && (!sys.EqualTo("")))                           ||
  (sample.EqualTo("gaustails_0.95")  && (!sys.EqualTo("")))                           ||
  (sample.EqualTo("JER")             && (!sys.EqualTo("/nominal")))                   ||
  (sample.EqualTo("JEC")             && (!sys.EqualTo("/up")&&!sys.EqualTo("/down"))) ||
  (sample.EqualTo("PU")              && (!sys.EqualTo("/up")&&!sys.EqualTo("/down"))) ||
  (sample.EqualTo("PLI")             && (!sys.EqualTo("/up")&&!sys.EqualTo("/down")))) skip=true;
  return skip;
}


// =======================================================================================================
MapVD ReadFile(TString filename_){

  MapVD content;
  for(const auto& mode: {"eta","eta_err","pT", "content"}) {content[mode] = {}; content[mode].clear();}
  if (gSystem->AccessPathName(filename_)) throw runtime_error("check: "+filename_);
  string line;
  ifstream myfile(filename_, ios::in);
  while (!myfile.eof()) {
    std::getline(myfile, line);
    if (line.length()==0) continue;
    if (line.find("pT")!=std::string::npos) continue; // skip first line
    VecD values; TString tok; int from = 0;
    while (((TString)line).Tokenize(tok, from, " ")){values.push_back(tok.Atof());}
    content["eta"].push_back((values[1]+values[0])/2);
    content["eta_err"].push_back((values[1]-values[0])/2);
    content["pT"].push_back(values[2]);
    content["content"].push_back(values[3]);
  }
  myfile.close();

  return content;
}

// =======================================================================================================
MapVDD GetBinContent(TString filename_, MapTS save, TString year_){
  MapVDD content;
  for(TString sample: samples){
    for(TString sys: syst){
      if(skipSample(sample, sys)) continue;
      TString path = save[year_];
      TString filename = path+"/pdfy/JERs/"+filename_+".txt";
      filename.ReplaceAll("standard", sample+sys);
      if(debug) cout << filename << endl;
      content[sample+sys] = ReadFile(filename);
    }
  }
  return content;
}

// =======================================================================================================
MapVDD CalculateUncertainty(MapVDD content){
  MapVDD uncertainties;
  for(TString sample: samples){
    for(TString sys: syst){
      if(skipSample(sample, sys)) continue;
      if(debug) cout << sample+sys << ": " << endl;
      for(unsigned int i=0; i<content["standard"]["eta"].size(); i++){
        uncertainties[sample+sys]["eta"].push_back(content["standard"]["eta"][i]);
        uncertainties[sample+sys]["pT"].push_back(content["standard"]["pT"][i]);
        uncertainties[sample+sys]["uncertainty"].push_back(content["standard"]["content"][i]-content[sample+sys]["content"][i]);
        if(debug) if(abs(content["standard"]["content"][i]-content[sample+sys]["content"][i])==0) cout << sample+sys << " is 0 at " <<  content["standard"]["eta"][i] << " and " << content["standard"]["pT"][i] << endl;
      }
    }
  }
  return uncertainties;
}

// =======================================================================================================
MapD SetSettings(double xmin, double xmax, double binning, double l_gauss = -0.05, double r_gauss = 0.05) {
  MapD hist;
  hist["xmin"]    = xmin;
  hist["xmax"]    = xmax;
  hist["binning"] = binning;
  hist["lgauss"]  = l_gauss;
  hist["rgauss"]  = r_gauss;
  return hist;
}

// =======================================================================================================
TH1F* FillHist(MapD settings, TString name, VecD content, VecD eta, VecD etaRange={0,5}) {
  TString h_name;
  double xmin = settings["xmin"]; double xmax = settings["xmax"]; double binning = settings["binning"];
  if(etaRange[0]==0&&etaRange[1]==5) h_name = name;
  else                               h_name = name+"_eta["+to_string(Round(etaRange[0]))+","+to_string(Round(etaRange[1]))+"]";
  TH1F* hist = new TH1F(h_name, "", binning, xmin, xmax);
  hist->SetTitle(h_name);
  for(unsigned int i=0; i<content.size(); i++){
    if(eta[i]<etaRange[0]||eta[i]>etaRange[1]) continue;
    hist->Fill(content[i]);
  }
  return hist;
}

// =======================================================================================================
MapVD CombineUncert(MapVDD map, VecS uncerts){
  MapVD combined;
  combined["eta"] = {}; combined["pT"] = {}; combined["uncertainty"] = {};
  for(unsigned int i=0; i<uncerts.size(); i++){
    combined["eta"].insert(combined["eta"].end(), map[uncerts[i]]["eta"].begin(), map[uncerts[i]]["eta"].end());
    combined["pT"].insert(combined["pT"].end(), map[uncerts[i]]["pT"].begin(), map[uncerts[i]]["pT"].end());
    combined["uncertainty"].insert(combined["uncertainty"].end(), map[uncerts[i]]["uncertainty"].begin(), map[uncerts[i]]["uncertainty"].end());
  }
  if(debug) cout << "Size new vector: " << combined["uncertainty"].size() << " " << combined["uncertainty"].size()/uncerts.size() << endl;
  return combined;
}

// =======================================================================================================
MapH SetHists(MapVD uncert, MapD settings, TString name) {
  MapH h_range;
  h_range["full"]       = FillHist(settings, name, uncert["uncertainty"], uncert["eta"], {0,5});
  h_range["eta0-2"]     = FillHist(settings, name, uncert["uncertainty"], uncert["eta"], {0,2});
  h_range["eta2-5"]     = FillHist(settings, name, uncert["uncertainty"], uncert["eta"], {2,5});
  h_range["eta0.0-1.3"] = FillHist(settings, name, uncert["uncertainty"], uncert["eta"], {0.0,1.3});
  h_range["eta1.3-2.5"] = FillHist(settings, name, uncert["uncertainty"], uncert["eta"], {1.3,2.5});
  h_range["eta2.5-3.1"] = FillHist(settings, name, uncert["uncertainty"], uncert["eta"], {2.5,3.1});
  h_range["eta3.1-5.2"] = FillHist(settings, name, uncert["uncertainty"], uncert["eta"], {3.1,5.2});
  return h_range;
}

// =======================================================================================================
void FitGaussian(TH1F* hist, double left, double right){
  TF1 *gauss = new TF1("gauss","gaus",left,right);
  gauss->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS() );
  hist->Fit("gauss", "RMQ+");
  cout << left << " " << right << endl;
}

// =======================================================================================================
void DrawAdditionalHist(TH1F* hist, TString option="", int Lcolor=kBlack, int Lstyle=1, int Lwidth=2, int Mstyle=kFullCircle, int Mcolor=kBlack){
  hist->SetLineStyle(Lstyle);
  hist->SetLineColor(Lcolor);
  hist->SetLineWidth(Lwidth);
  hist->SetMarkerStyle(Mstyle);
  hist->SetMarkerColor(Mcolor);
  hist->Draw(option+" SAME");
}

// =======================================================================================================
void PlotUncert(MapH hists, TString path, MapD settings, TString uncert) {

  TH1F* h_main = hists["full"];
  FitGaussian(h_main, settings["lgauss"], settings["rgauss"]);
  TF1* gaus = h_main->GetFunction("gauss");
  double gmean      = gaus->GetParameter(1);
  double gmean_err  = gaus->GetParError(1);
  double gwidth     = gaus->GetParameter(2);
  double gwidth_err = gaus->GetParError(2);
  double hmean          = h_main->GetMean();
  double hwidth         = h_main->GetRMS();
  double hmean_err      = h_main->GetMeanError();
  double hwidth_err     = h_main->GetRMSError();

  TString h_name = h_main->GetTitle();
  TCanvas* canv = tdrCanvas(h_name, settings["xmin"], settings["xmax"], 0, h_main->GetMaximum()*1.1, "Uncertainty", "Events");
  canv->SetTickx(0);
  canv->SetTicky(0);
  tdrDraw(h_main, "HIST");
  gaus->Draw("SAME");
  DrawAdditionalHist(hists["eta0.0-1.3"], "hist", kRed+2,     1);
  DrawAdditionalHist(hists["eta1.3-2.5"], "hist", kBlue+2,    1);
  DrawAdditionalHist(hists["eta2.5-3.1"], "hist", kGreen+2,   1);
  DrawAdditionalHist(hists["eta3.1-5.2"], "hist", kMagenta+2, 1);
  TLegend *legend, *gleg, *hleg;
  legend = tdrLeg(0.65,0.65,0.90,0.9, 0.04);
  legend->AddEntry((TObject*) 0,    uncert,           "");
  legend->AddEntry(hists["full"],       "full range",    "l");
  legend->AddEntry(hists["eta0.0-1.3"], "#eta in [0.0,1.3]", "l");
  legend->AddEntry(hists["eta1.3-2.5"], "#eta in [1.3,2.5]", "l");
  legend->AddEntry(hists["eta2.5-3.1"], "#eta in [2.5,3.1]", "l");
  legend->AddEntry(hists["eta3.1-5.2"], "#eta in [3.1,5.2]", "l");
  gleg = tdrLeg(0.65,0.55,0.90,0.65, 0.025);
  gleg->AddEntry((TObject*) 0, "Gaussian: ", "");
  char line[100];
  sprintf(line, "Mean = %.5f #pm %.5f", gmean, gmean_err);    gleg->AddEntry((TObject*)0, line, "");
  sprintf(line, "Width = %.5f #pm %.5f", gwidth, gwidth_err); gleg->AddEntry((TObject*)0, line, "");
  hleg = tdrLeg(0.65,0.45,0.90,0.55, 0.025);
  hleg->AddEntry((TObject*) 0, "Hist: ", "");
  sprintf(line, "Mean = %.5f #pm %.5f", hmean, hmean_err);    hleg->AddEntry((TObject*)0, line, "");
  sprintf(line, "Width = %.5f #pm %.5f", hwidth, hwidth_err); hleg->AddEntry((TObject*)0, line, "");

  canv->Print(path+"/pdfy/kValues/uncert_"+h_name+".pdf","pdf");
}


// ===========================================================================
//                               MAIN PROGRAM
// ===========================================================================

void studyFitParameters(){

  ////////////////////////////////////////////////////////////////////////////
  //    I load all parameters I will need                                   //
  ////////////////////////////////////////////////////////////////////////////

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
  map_SAVE["UL16preVFP"]  = path_+"eta_JER/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF";
  map_SAVE["UL16postVFP"] = path_+"eta_JER/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH";
  map_SAVE["UL17"]        = path_+"eta_JER/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF";
  map_SAVE["UL18"]        = path_+"eta_JER/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD";


  ////////////////////////////////////////////////////////////////////////////
  //    I plot all values of the fit parameters                             //
  ////////////////////////////////////////////////////////////////////////////

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

  ////////////////////////////////////////////////////////////////////////////
  //    I load all values for the Error studies                             //
  ////////////////////////////////////////////////////////////////////////////

  cout << "\nStarting with uncertainty study ... " << endl;
  // MapVDD GetBinContent(TString filename_, MapTS save, TString year_)
  MapVDD content = GetBinContent("values_FE_correlated", map_SAVE, "UL18");
  MapVDD uncertainty = CalculateUncertainty(content);
  uncertainty["JEC"] = CombineUncert(uncertainty, {"JEC/up", "JEC/down"});
  uncertainty["PU"]  = CombineUncert(uncertainty, {"PU/up", "PU/down"});
  uncertainty["PLI"] = CombineUncert(uncertainty, {"PLI/up", "PLI/down"});
  uncertainty["all"] = CombineUncert(uncertainty, {"alpha", "gaustails_0.95", "JEC", "PLI", "PU"});

  // Set settings for histograms
  // MapD SetSettings(double xmin, double xmax, double binning, double l_gauss, double r_gauss)
  MapDD settings;
  settings["all"]            = SetSettings(-0.03,   0.03, 30,  -0.01,  0.01);
  settings["alpha"]          = SetSettings(-0.03,   0.03, 30,  -0.01,  0.01);
  settings["gaustails_0.95"] = SetSettings(-0.05,   0.05, 20, -0.005, 0.025);
  settings["JEC"]            = SetSettings(-0.02,   0.02, 20,  -0.01,  0.01);
  settings["PU"]             = SetSettings(-0.02,   0.02, 40,  -0.05,  0.05);
  settings["PLI"]            = SetSettings(-0.003, 0.003, 30,  -0.02,  0.02);

  settings["JEC/up"]         = SetSettings(-0.02,   0.02, 20,  -0.005,  0.005);
  settings["JEC/down"]       = SetSettings(-0.02,   0.02, 20,  -0.005,  0.005);
  settings["JER/nominal"]    = SetSettings(-0.05,   0.05, 25,   -0.02,   0.01);
  settings["PU/up"]          = SetSettings(-0.02,   0.02, 40,  -0.005,  0.005);
  settings["PU/down"]        = SetSettings(-0.02,   0.02, 40,  -0.005,  0.005);
  settings["PLI/up"]         = SetSettings(-0.003, 0.003, 30,     0.0, 0.0015);
  settings["PLI/down"]       = SetSettings(-0.003, 0.003, 30, -0.0015,    0.0);

  if(same_range){
    settings["all"]            = SetSettings(-0.05,   0.05, 30,  -0.01,  0.01);
    settings["alpha"]          = SetSettings(-0.05,   0.05, 30,  -0.01,  0.01);
    settings["gaustails_0.95"] = SetSettings(-0.05,   0.05, 30, -0.005, 0.025);
    settings["JEC"]            = SetSettings(-0.05,   0.05, 30,  -0.01,  0.01);
    settings["PU"]             = SetSettings(-0.05,   0.05, 30,  -0.05,  0.05);
    settings["PLI"]            = SetSettings(-0.05,   0.05, 30,  -0.02,  0.02);

    settings["JEC/up"]         = SetSettings(-0.05,   0.05, 30, -0.005,  0.005);
    settings["JEC/down"]       = SetSettings(-0.05,   0.05, 30, -0.005,  0.005);
    settings["JER/nominal"]    = SetSettings(-0.05,   0.05, 30,  -0.02,   0.01);
    settings["PU/up"]          = SetSettings(-0.05,   0.05, 30, -0.005,  0.005);
    settings["PU/down"]        = SetSettings(-0.05,   0.05, 30, -0.005,  0.005);
    settings["PLI/up"]         = SetSettings(-0.05,   0.05, 30,    0.0, 0.0015);
    settings["PLI/down"]       = SetSettings(-0.05,   0.05, 30, -0.005,  0.005);
  }

  cout << "Set Hists ... " << endl;
  // MapH SetHists(MapVD uncert, MapD settings, TString name)
  MapMH hists;
  hists["all"]            = SetHists(uncertainty["all"],            settings["all"],            "all");
  hists["alpha"]          = SetHists(uncertainty["alpha"],          settings["alpha"],          "alpha");
  hists["JEC"]            = SetHists(uncertainty["JEC"],            settings["JEC"],            "JEC");
  hists["PU"]             = SetHists(uncertainty["PU"],             settings["PU"],             "PU");
  hists["PLI"]            = SetHists(uncertainty["PLI"],            settings["PLI"],            "PLI");
  hists["JEC/up"]         = SetHists(uncertainty["JEC/up"],         settings["JEC/up"],         "JECup");
  hists["JEC/down"]       = SetHists(uncertainty["JEC/down"],       settings["JEC/down"],       "JECdown");
  hists["gaustails_0.95"] = SetHists(uncertainty["gaustails_0.95"], settings["gaustails_0.95"], "gaustails95");
  hists["JER/nominal"]    = SetHists(uncertainty["JER/nominal"],    settings["JER/nominal"],    "JERnominal");
  hists["PU/up"]          = SetHists(uncertainty["PU/up"],          settings["PU/up"],          "PUup");
  hists["PU/down"]        = SetHists(uncertainty["PU/down"],        settings["PU/down"],        "PUdown");
  hists["PLI/up"]         = SetHists(uncertainty["PLI/up"],         settings["PLI/up"],         "PLIup");
  hists["PLI/down"]       = SetHists(uncertainty["PLI/down"],       settings["PLI/down"],       "PLIdown");


  cout << "Plot Uncetainties ... " << endl;
  // void PlotUncert(MapH hists, TString path, MapD settings, TString uncert)
  PlotUncert(hists["all"],            map_SAVE["UL18"], settings["all"],            "all");
  PlotUncert(hists["alpha"],          map_SAVE["UL18"], settings["alpha"],          "alpha");
  PlotUncert(hists["JEC"],            map_SAVE["UL18"], settings["JEC"],            "JEC");
  PlotUncert(hists["JEC/up"],         map_SAVE["UL18"], settings["JEC/up"],         "JECup");
  PlotUncert(hists["JEC/down"],       map_SAVE["UL18"], settings["JEC/down"],       "JECdown");
  PlotUncert(hists["gaustails_0.95"], map_SAVE["UL18"], settings["gaustails_0.95"], "gaustails95");
  PlotUncert(hists["JER/nominal"],    map_SAVE["UL18"], settings["JER/nominal"],    "JERnominal");
  PlotUncert(hists["PU"],             map_SAVE["UL18"], settings["PU"],             "PU");
  PlotUncert(hists["PU/up"],          map_SAVE["UL18"], settings["PU/up"],          "PUup");
  PlotUncert(hists["PU/down"],        map_SAVE["UL18"], settings["PU/down"],        "PUdown");
  PlotUncert(hists["PLI"],            map_SAVE["UL18"], settings["PLI"],            "PLI");
  PlotUncert(hists["PLI/up"],         map_SAVE["UL18"], settings["PLI/up"],         "PLIup");
  PlotUncert(hists["PLI/down"],       map_SAVE["UL18"], settings["PLI/down"],       "PLIdown");

}
