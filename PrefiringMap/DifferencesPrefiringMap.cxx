#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TChain.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TString.h>

typedef std::vector<double> VecD;
typedef std::vector<vector<double>> VecDD;

typedef std::map<TString, double> MapD;
typedef std::map<TString, VecD>   MapVD;
typedef std::map<TString, MapVD>  MapMVD;

typedef std::map<int, VecD>       MapiVD;
typedef std::map<TString, MapiVD> MapMiVD;

typedef std::map<TString, TH2F*>  MapH2;
typedef std::map<TString, MapH2>  MapMH2;

bool debug = false;

VecD jet_pt_center       = { 12.5,  17.5,  22.5,  27.5,  32.5,  37.5,    45,    60,    85,   125,   175,   250,   400}; // equal for jet & jetem and EOY & UL.
VecD photonEOY_pt_center = {               22.5,  27.5,  32.5,  37.5,    45,    60,    85,   150,   350};               // Note the different binning
VecD photonUL_pt_center  = {        17.5,  22.5,  27.5,  32.5,  37.5,    45,    60,    85,   150,   350};               // Note the different binning

// same for all three versions (jet, jetem, photon)
VecD EOY_eta_center = { 0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.625, 2.875,  3.05,   3.3,  3.75,   4.5}; // also negative
VecD UL_eta_center  = {    -4,-2.875,-2.625,-2.375,-2.125, 0, 2.125, 2.375, 2.625, 2.875,     4};

VecD wanted_eta = {-2.875, -2.625, -2.375, -2.125, 2.125, 2.375, 2.625, 2.875};

// MapVD same_eta_bins; // Calculated via Code
// same_eta_bins["EOY"] = {  5,  6,  7,  8, 25, 26, 27, 28};
// same_eta_bins["UL"]  = {  2,  3,  4,  5,  7,  8,  9, 10};

// =======================================================================================================
void print_bin_content(TH2F* hist, TString name){
  // function for debugging purpose only
  // bin with x or y = 0 does not contribute to histogram
  cout << setprecision(5) << endl;
  cout << "================================================================================================" << endl;
  cout << " ------ "+name  << endl;
  cout << "================================================================================================" << endl;
  cout << right << setw(2) << "x" << setw(4) << "y" << setw(5) << "bin" << setw(9) << "eta" << setw(6) << "e.w." << setw(9) << "pt" << setw(6) << "p.w." << setw(13) << "content" << setw(13) << "contentBin"<< "\n";
  cout << "------------------------------------------------------------------------------------------------" << endl;

  for(unsigned int x=1; x<=hist->GetNbinsX(); x++) // size of binEOY and binsUL equal
  {
    for(unsigned int y=1; y<=hist->GetNbinsY(); y++)
    {
      double eta        = hist->GetXaxis()->GetBinCenter(x);
      double etaWidth   = hist->GetXaxis()->GetBinWidth(x);
      double pt         = hist->GetYaxis()->GetBinCenter(y);
      double ptWidth    = hist->GetYaxis()->GetBinWidth(y);
      double bin        = hist->GetBin(x,y);
      double content    = hist->GetBinContent(x,y);
      double contentBin = hist->GetBinContent(bin);
      cout << setw(2) << x << setw(4) << y << setw(5) << bin << setw(9) << eta << setw(6) << etaWidth << setw(9) << pt << setw(6) << ptWidth << setw(13) << content << setw(13) << contentBin << "\n";
    }
  }
  cout << "================================================================================================" << endl;
  cout << left << endl;
}

// =======================================================================================================
void print_vector(VecD vec, TString name="", int space=12){
  bool noName = name.EqualTo("");
  if(!noName && name.Contains("EOY")) cout << name+" = {";
  if(!noName && name.Contains("UL"))  cout << name+"  = {";
  for(unsigned int i=0; i<vec.size(); i++)
  {
    if(vec[i]<0) continue; // since negative and positiv bins are equal
    cout << setw(space) << vec[i];
    if(!(i==vec.size()-1)) cout << ", ";
  }
  if(!noName) cout << "}" << endl;
  else        cout << endl;
}

// =======================================================================================================
void draw_map(TH2F* map, TString name, TString ztitle="weight", double ymin=0, double ymax=500, double zmin=0, double zmax=1){
  gStyle->SetOptStat(0);
  map->GetXaxis()->SetRangeUser(-4, 4);
  map->GetYaxis()->SetRangeUser(ymin, ymax);
  map->GetZaxis()->SetRangeUser(zmin, zmax);

  map->SetTitle(name);
  map->GetXaxis()->SetTitle("#eta");
  map->GetYaxis()->SetTitle("p_{jet}");
  map->GetZaxis()->CenterTitle();
  map->GetZaxis()->SetTitle(ztitle);
  map->GetZaxis()->SetTitleOffset(0.8);

  TCanvas *A = new TCanvas(name, name, 800, 600);
  A->SetRightMargin(0.15);
  map->Draw("COLZ");
  A->SaveAs("/nfs/dust/cms/user/paaschal/CMSSW_106X_v1/CMSSW_10_6_13/src/UHH2/DiJetJERC/PrefiringMap/maps/"+name+".pdf");
}

// =======================================================================================================
VecD save_bins(VecD vec, double ValueUp, double ValueDown){
  VecD bins = {};
  for(unsigned int bin=0; bin<vec.size(); bin++)
  {
    if(abs(vec[bin])>ValueUp && abs(vec[bin])<ValueDown) bins.push_back(bin+1);
  }
  return bins;
}

// =======================================================================================================
MapVD get_same_bins(VecD binsEOY, VecD binsUL, double ValueUp, double ValueDown){
  MapVD same_bins;
  same_bins["EOY"] = save_bins(binsEOY, ValueUp, ValueDown);
  same_bins["UL"]  = save_bins(binsUL,  ValueUp, ValueDown);
  if(debug) print_vector(same_bins["EOY"], "EOY", 3);
  if(debug) print_vector(same_bins["UL"],  "UL",  3);
  return same_bins;
}

// =======================================================================================================
bool value_in_vec(VecD vec, double value){
  cout << "\nIn value_in_vec()" << endl;
  bool isIn = false;
  for(unsigned int i=0; i<vec.size(); i++)
  {
    if(vec[i]==value) isIn = true;
  }
  return isIn;
}

// =======================================================================================================
bool within_uncertainties(double contentEOY, double contentUL, double uncert){
  bool EOYisGreater = false;
  if(contentEOY>contentUL) EOYisGreater = true;

  bool isInUncert = false;
  double above, below;
  if(EOYisGreater)
  {
    above      = contentEOY - contentEOY*uncert;
    below      = contentUL  + contentUL*uncert;
    isInUncert = (above - below)<=0;
  }
  else
  {
    below      = contentEOY + contentEOY*uncert;
    above      = contentUL  - contentUL*uncert;
    isInUncert = (above - below)<=0;
  }
  return isInUncert;
}

// =======================================================================================================
void set_bin_content(TH2F* hist, MapiVD content, MapVD eta_bins, VecD pt_UL){
  for(unsigned int x=0; x<eta_bins["UL"].size(); x++) // size of binEOY and binsUL equal
  {
    for(unsigned int y=1; y<=pt_UL.size(); y++)
    {
      hist->SetBinContent(eta_bins["UL"][x], y, content[eta_bins["UL"][x]][y-1]);
    }
  }
}

// =======================================================================================================
MapMiVD get_bin_content(MapH2 hists, MapVD eta_bins, VecD pt_EOY, VecD pt_UL, TString name, bool isPhoton=false){

  int all       = 0;
  int inUncert  = 0;
  int inUncertlow  = 0;
  int outUncert = 0;

  ofstream output;
  output.open("/nfs/dust/cms/user/paaschal/CMSSW_106X_v1/CMSSW_10_6_13/src/UHH2/DiJetJERC/PrefiringMap/"+name+".txt");

  output << right << fixed << setprecision(5);
  output << setw(10) << "eta" << setw(10) << "pt" << setw(10) << "EOY" << setw(10) << "UL" << setw(10) << "diff" << setw(10) << "ratio" << "\n";


  MapMiVD bin_content;
  for(unsigned int x=0; x<eta_bins["UL"].size(); x++) // size of binEOY and binsUL equal
  {
    for(unsigned int y=1; y<=pt_UL.size(); y++)
    {
      int yEOY = y; int yUL = y;
      if(isPhoton&&y>1) yEOY = y-1;

      double contentEOY  = hists["EOY"]->GetBinContent(eta_bins["EOY"][x],yEOY);
      if(isPhoton&&y==1) contentEOY = 0.; // because EOY has one pt bin less than UL for photon.
      double contentUL   = hists["UL"]->GetBinContent(eta_bins["UL"][x],yUL);

      double diff   = abs(contentEOY-contentUL);
      double ratio_ = contentEOY/contentUL;
      double ratio;
      if(ratio_>100.) ratio = 10.;
      else            ratio = ratio_;

      bin_content["EOY"][eta_bins["EOY"][x]].push_back(contentEOY);
      bin_content["UL"][eta_bins["UL"][x]].push_back(contentUL);
      bin_content["diff"][eta_bins["UL"][x]].push_back(diff);   // Clone of UL
      bin_content["ratio"][eta_bins["UL"][x]].push_back(ratio); // Clone of UL
      output << setw(10) << wanted_eta[x] << setw(10) << pt_UL[y-1] << setw(10) << contentEOY << setw(10) << contentUL << setw(10) << diff << setw(10) << ratio << "\n";

      // count bins
      bool withinUncert = within_uncertainties(contentEOY, contentUL, 0.25);
      bool below100     = hists["UL"]->GetYaxis()->GetBinCenter(yUL)<100;
      bool both      = !withinUncert&&below100;
      bool checkboth = withinUncert&&below100;
      all++;
      if(withinUncert) inUncert++;
      if(checkboth)    outUncert++;
      if(both)         inUncertlow++;
    }
  }
  cout << setprecision(3) << "number bins: " << setw(4) << all << " | within Uncert: " << setw(3) << inUncert << " (" << inUncert/(double)all*100 << "%)" << endl;;
  return bin_content;
}

// =======================================================================================================
MapMVD get_bin_info(TH2F* hist, TString name=""){
  MapVD X_bin_info, Y_bin_info;
  MapMVD hist_bin_info;
  VecD BinUp, BinDown, BinCenter;

  double nx = hist->GetNbinsX();
  double ny = hist->GetNbinsY();

  TAxis *axis;
  for(int x=1; x<=nx; x++)
  {
    axis = hist->GetXaxis();
    X_bin_info["low"].push_back(axis->GetBinLowEdge(x));
    X_bin_info["up"].push_back(axis->GetBinUpEdge(x));
    X_bin_info["center"].push_back(axis->GetBinCenter(x));
    X_bin_info["width"].push_back(axis->GetBinWidth(x));
  }

  for(int y=1; y<=ny; y++)
  {
    axis = hist->GetYaxis();
    Y_bin_info["low"].push_back(axis->GetBinLowEdge(y));
    Y_bin_info["up"].push_back(axis->GetBinUpEdge(y));
    Y_bin_info["center"].push_back(axis->GetBinCenter(y));
    Y_bin_info["width"].push_back(axis->GetBinWidth(y));
  }

  hist_bin_info["eta"] = X_bin_info;
  hist_bin_info["pt"] = Y_bin_info;

  if(debug) print_vector(hist_bin_info["eta"]["center"], name+"_eta_center", 5);
  if(debug) print_vector(hist_bin_info["pt"]["center"],  name+"_pt_center ", 5);
  if(debug) print_vector(hist_bin_info["pt"]["low"],     name+"_pt_low    ", 5);

  return hist_bin_info;

}


// =======================================================================================================
//                                  MAIN
// =======================================================================================================
void DifferencesPrefiringMap(){
  TFile *maps = new TFile("L1PrefiringMaps_WithUL17.root");

  TString jetemptvseta_EOY  = "L1prefiring_jetemptvseta_2017BtoF";
  TString jetptvseta_EOY    = "L1prefiring_jetptvseta_2017BtoF";
  TString photonptvseta_EOY = "L1prefiring_photonptvseta_2017BtoF";

  TString jetemptvseta_UL   = "L1prefiring_jetemptvseta_UL2017BtoF";
  TString jetptvseta_UL     = "L1prefiring_jetptvseta_UL2017BtoF";
  TString photonptvseta_UL  = "L1prefiring_photonptvseta_UL2017BtoF";

  TH2F* h_jetemptvseta_EOY  = (TH2F*)maps->Get(jetemptvseta_EOY);
  TH2F* h_jetptvseta_EOY    = (TH2F*)maps->Get(jetptvseta_EOY);
  TH2F* h_photonptvseta_EOY = (TH2F*)maps->Get(photonptvseta_EOY);

  TH2F* h_jetemptvseta_UL   = (TH2F*)maps->Get(jetemptvseta_UL);
  TH2F* h_jetptvseta_UL     = (TH2F*)maps->Get(jetptvseta_UL);
  TH2F* h_photonptvseta_UL  = (TH2F*)maps->Get(photonptvseta_UL);

  // Clone UL; lesser eta bins
  TH2F* h_jetptvseta_diff    = (TH2F*)h_jetptvseta_UL->Clone();
  TH2F* h_jetemptvseta_diff  = (TH2F*)h_jetemptvseta_UL->Clone();
  TH2F* h_photonptvseta_diff = (TH2F*)h_photonptvseta_UL->Clone();

  TH2F* h_jetptvseta_ratio    = (TH2F*)h_jetptvseta_UL->Clone();
  TH2F* h_jetemptvseta_ratio  = (TH2F*)h_jetemptvseta_UL->Clone();
  TH2F* h_photonptvseta_ratio = (TH2F*)h_photonptvseta_UL->Clone();

  MapMH2 hists;
  hists["jet"]["EOY"]     = h_jetptvseta_EOY;
  hists["jetem"]["EOY"]   = h_jetemptvseta_EOY;
  hists["photon"]["EOY"]  = h_photonptvseta_EOY;

  hists["jet"]["UL"]      = h_jetptvseta_UL;
  hists["jetem"]["UL"]    = h_jetemptvseta_UL;
  hists["photon"]["UL"]   = h_photonptvseta_UL;

  hists["jet"]["diff"]    = h_jetptvseta_diff;
  hists["jetem"]["diff"]  = h_jetemptvseta_diff;
  hists["photon"]["diff"] = h_photonptvseta_diff;

  hists["jet"]["ratio"]    = h_jetptvseta_ratio;
  hists["jetem"]["ratio"]  = h_jetemptvseta_ratio;
  hists["photon"]["ratio"] = h_photonptvseta_ratio;

  // Get Bin Inofrmation: Center, Lower Edge, Upper Edge, etc.. Main purpose: compare both 2D Hists
  if(debug) cout << "\nIn get_bin_info() with jet" << endl;
  MapMVD bin_jet_EOY    = get_bin_info(h_jetptvseta_EOY, "jetEOY");
  MapMVD bin_jet_UL     = get_bin_info(h_jetptvseta_UL,  "jetUL");

  if(debug) cout << "\nIn get_bin_info() with jet" << endl;
  MapMVD bin_jetem_EOY  = get_bin_info(h_jetemptvseta_EOY, "jetemEOY");
  MapMVD bin_jetem_UL   = get_bin_info(h_jetemptvseta_UL,  "jetemUL");

  if(debug) cout << "\nIn get_bin_info() with photon" << endl;
  MapMVD bin_photon_EOY = get_bin_info(h_photonptvseta_EOY, "photonEOY");
  MapMVD bin_photon_UL  = get_bin_info(h_photonptvseta_UL,  "photonUL");

  if(debug) cout << "\nIn save_bins()" << endl;
  // pT bins are equal for jet and jetem. Need a different treatment for photons.
  // Eta bins are the same for all the cases in each year.
  MapVD same_eta_bins = get_same_bins(bin_jet_EOY["eta"]["center"], bin_jet_UL["eta"]["center"], 2, 3);

  // MapMiVD get_bin_content(MapH2 hists, MapVD eta_bins, VecD pt_EOY, VecD pt_UL, TString name, bool isPhoton=false)
  MapMiVD content_jet    = get_bin_content(hists["jet"],    same_eta_bins, bin_jet_EOY["pt"]["center"],    bin_jet_UL["pt"]["center"],    "jetptvseta",    false);
  MapMiVD content_jetem  = get_bin_content(hists["jetem"],  same_eta_bins, bin_jetem_EOY["pt"]["center"],  bin_jetem_UL["pt"]["center"],  "jetemptvseta",  false);
  MapMiVD content_photon = get_bin_content(hists["photon"], same_eta_bins, bin_photon_EOY["pt"]["center"], bin_photon_UL["pt"]["center"], "photonptvseta", true);

  // set_bin_content(TH2F* hist, MapiVD content, MapVD eta_bins, VecD pt_UL)
  // set_bin_content(hists["jet"]["diff"],     content_jet["diff"],     same_eta_bins, bin_jet_UL["pt"]["center"]);
  // set_bin_content(hists["jetem"]["diff"],   content_jetem["diff"],   same_eta_bins, bin_jetem_UL["pt"]["center"]);
  // set_bin_content(hists["photon"]["diff"],  content_photon["diff"],  same_eta_bins, bin_photon_UL["pt"]["center"]);

  set_bin_content(hists["jet"]["ratio"],    content_jet["ratio"],    same_eta_bins, bin_jet_UL["pt"]["center"]);
  set_bin_content(hists["jetem"]["ratio"],  content_jetem["ratio"],  same_eta_bins, bin_jetem_UL["pt"]["center"]);
  set_bin_content(hists["photon"]["ratio"], content_photon["ratio"], same_eta_bins, bin_photon_UL["pt"]["center"]);

  if(debug) cout << "\nIn draw_map() with EOY" << endl;
  draw_map(hists["jet"]["EOY"],    "jetptvseta_EOY",    20);
  draw_map(hists["jetem"]["EOY"],  "jetemptvseta_EOY",  20);
  draw_map(hists["photon"]["EOY"], "photonptvseta_EOY", 20);

  if(debug) cout << "\nIn draw_map() with UL" << endl;
  draw_map(hists["jet"]["UL"],   "jetptvseta_UL",    20);
  draw_map(hists["jetem"]["UL"], "jetemptvseta_UL",  20);
  draw_map(hists["photon"]["UL"],"photonptvseta_UL", 20);

  if(debug) cout << "\nIn draw_map() with Diff" << endl;
  // draw_map(hists["jet"]["diff"],    "jetptvseta_diff",    "|EOY-UL|", 20, 500, 0, 0.07);
  // draw_map(hists["jetem"]["diff"],  "jetemptvseta_diff",  "|EOY-UL|", 20, 500, 0, 0.2);
  // draw_map(hists["photon"]["diff"], "photonptvseta_diff", "|EOY-UL|", 20, 500, 0, 0.35);

  if(debug) cout << "\nIn draw_map() with ratio" << endl;
  draw_map(hists["jet"]["ratio"],    "jetptvseta_ratio",    "EOY/UL", 20, 500, 0.5, 1.5);
  draw_map(hists["jetem"]["ratio"],  "jetemptvseta_ratio",  "EOY/UL", 20, 500, 0.5, 1.5);
  draw_map(hists["photon"]["ratio"], "photonptvseta_ratio", "EOY/UL", 20, 500, 0.5, 1.5);



}
