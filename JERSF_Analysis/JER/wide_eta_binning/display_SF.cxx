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

TGraphErrors* CreateTGraphSF(TString filename = "SF_final_tex.txt") {

  filename += "SF_final_tex.txt";

  if (gSystem->AccessPathName(filename)) throw runtime_error("check: "+filename);

  VecD eta, eta_center, eta_err, SF, SF_Err;
  eta.push_back(0.); string line;
  ifstream myfile(filename, ios::in);
  while (!myfile.eof()) {
    std::getline(myfile, line); if (line.length()==0) continue;
    line = line.substr(line.find("{{")+2,line.find("}}")-2);
    VecD values; TString tok; int from = 0;
    while (((TString)line).Tokenize(tok, from, ", ")) {values.push_back(tok.Atof());}
    eta.push_back(values[0]);
    SF.push_back(values[1]);
    SF_Err.push_back(values[2]-values[1]);
  }
  myfile.close();

  for (size_t i = 0; i < eta.size()-1; i++) { eta_center.push_back((eta[i+1]+eta[i])/2); eta_err.push_back((eta[i+1]-eta[i])/2); }
  TGraphErrors* gr = new TGraphErrors(SF.size(), &(eta_center[0]), &SF[0], &(eta_err[0]), &SF_Err[0]);
  return gr;
}

TGraphErrors* CreateTGraphSF(VecD eta, std::vector<std::vector<double>> jer) {
  VecD eta_center, eta_err, SF, SF_Err;
  for (unsigned int i = 0; i < eta.size()-1; i++) {
    eta_center.push_back((eta[i+1]+eta[i])/2);
    eta_err.push_back((eta[i+1]-eta[i])/2);
    SF.push_back(jer[i][0]);
    SF_Err.push_back(jer[i][1]);
  }
  TGraphErrors* gr = new TGraphErrors(SF.size(), &(eta_center[0]), &SF[0], &(eta_err[0]), &SF_Err[0]);
  return gr;
}

TGraphErrors* CreateTGraphSFfromFile(TString ver, TString fname="Source" ) {
  VecD eta_center, eta_err, SF, SF_Err;
  TString filename = std::getenv("CMSSW_BASE"); filename += "/src/UHH2/JRDatabase/textFiles/"+ver+"_MC/"+fname+".txt";
  std::ifstream infile(filename);
  std::string line;
  std::getline(infile, line); //skip header
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    if (line=="" || line=="\n") continue;
    TString tok;
    int from = 0;
    double eta_max=0, eta_min=0, sf=0, sf_down=0, sf_err=0;
    while (((TString)line).Tokenize(tok, from, "[ |]")) {
      if (atof(tok)<0) break;
      if (from == 6)  eta_min = atof(tok);
      if (from == 12) eta_max = atof(tok);
      if (from == 21) sf = atof(tok);
      if (from >= 27) sf_down = atof(tok);
    }
    if (sf<=0) continue;
    eta_center.push_back((eta_max+eta_min)/2);
    eta_err.push_back((eta_max-eta_min)/2);
    SF.push_back(sf);
    SF_Err.push_back(sf-sf_down);
  }
  infile.close();

  TGraphErrors* gr = new TGraphErrors(SF.size(), &(eta_center[0]), &SF[0], &(eta_err[0]), &SF_Err[0]);
  return gr;
}


void PlotCanvas(TString cName, MapG map_gr, VecS names, VecI colors, VecI shapes){
  lumi_13TeV = "RunII Legacy";

  //TCanvas* canv = tdrCanvas(cName, -0.5, 5.191+0.5, 0., 3.0, "|#eta|", "JER SF");
  TCanvas* canv = tdrCanvas(cName, -0.5, 5.191+0.5, 0.85, 1.8, "|#eta|", "JER SF");
  TLegend* leg = tdrLeg(0.65,0.67,0.79,0.90, 0.035, 42, kBlack);
  for(int i=0; i<names.size(); i++){
    TString name = names[i];
    int color = colors[i];
    int shape = shapes[i];
    tdrDraw(map_gr[name], "P5", shape, color, kSolid, color, 1001, color, 0.15);
    leg->AddEntry(map_gr[name], name,"fp");
  }

  leg->Draw("same");
  canv->SaveAs("plots/"+cName+".pdf");
}

void display_SF() {
  TString path_ = std::getenv("CMSSW_BASE"); path_ += "/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/";
  TString outdir = "plots/";
  gPrintViaErrorHandler = kTRUE;
  gErrorIgnoreLevel = kFatal;

  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  // TString lumi;
  //
  // if (year.Contains("18"))        lumi = "59.74";
  // if (year.Contains("17"))        lumi = "41.53";
  // if (year.Contains("16"))        lumi = "";
  // if (year.Contains("16preVFP"))  lumi = "";
  // if (year.Contains("16postVFP")) lumi = "";

  // lumi_13TeV = lumi+" fb^{-1}";
  lumi_13TeV = "RunII Legacy";


  MapG map_gr;
  map_gr["EOY16"] = CreateTGraphSFfromFile("Summer16_25nsV1", "Summer16_25nsV1_MC_SF_AK4PFchs");
  map_gr["EOY17"] = CreateTGraphSFfromFile("Fall17_V3", "Fall17_V3_MC_SF_AK4PFchs");
  // map_gr["EOY18"] = CreateTGraphSFfromFile("Autumn18_V7b");
  map_gr["EOY18"] = CreateTGraphSF("/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/MergeL2Res/Autumn18_V17_save/AK4CHS/standard/QCDHT/RunABCD/");

  map_gr["UL16"] = CreateTGraphSFfromFile("Summer19UL16_JRV1", "Summer19UL16_JRV1_MC_SF_AK4PFchs");
  map_gr["UL17"] = CreateTGraphSFfromFile("Summer19UL17_JRV3", "Summer19UL17_JRV3_MC_SF_AK4PFchs");
  map_gr["UL18"] = CreateTGraphSFfromFile("Summer19UL18_JRV2", "Summer19UL18_JRV2_MC_SF_AK4PFchs");

  map_gr["UL16preVFP"]  = CreateTGraphSF(path_+"eta_JER_tot/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_gr["UL16postVFP"] = CreateTGraphSF(path_+"eta_JER_tot/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/");

  map_gr["UL16_new"] = CreateTGraphSF(path_+"eta_JER_tot/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/");
  map_gr["UL17_new"] = CreateTGraphSF(path_+"eta_JER_tot/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_gr["UL18_new"] = CreateTGraphSF(path_+"eta_JER_tot/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/");
  map_gr["RunII_new"] = CreateTGraphSF(path_+"eta_JER_tot/Legacy/Summer19Legacy/AK4CHS/standard/QCDHT/RunII/");

  map_gr["UL16preVFP_eta"]  = CreateTGraphSF(path_+"eta_simple/UL16preVFP/Summer19UL16APV_V3/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_gr["UL16postVFP_eta"] = CreateTGraphSF(path_+"eta_simple/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/");
  map_gr["UL16_eta"] = CreateTGraphSF(path_+"eta_simple/UL16postVFP/Summer19UL16_V2/AK4CHS/standard/QCDHT/RunFGH/");
  map_gr["UL17_eta"] = CreateTGraphSF(path_+"eta_simple/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/");
  map_gr["UL18_eta"] = CreateTGraphSF(path_+"eta_simple/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/");
  map_gr["RunII_eta"] = CreateTGraphSF(path_+"eta_simple/Legacy/Summer19Legacy/AK4CHS/standard/QCDHT/RunII/");

  map_gr["UL17_B"] = CreateTGraphSF(path_+"eta_JER_tot/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunB/");
  map_gr["UL17_C"] = CreateTGraphSF(path_+"eta_JER_tot/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunC/");
  map_gr["UL17_D"] = CreateTGraphSF(path_+"eta_JER_tot/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunD/");
  map_gr["UL17_E"] = CreateTGraphSF(path_+"eta_JER_tot/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunE/");
  map_gr["UL17_F"] = CreateTGraphSF(path_+"eta_JER_tot/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunF/");
  map_gr["UL17_BCDEF"] = CreateTGraphSF(path_+"eta_JER_tot/UL17/Summer19UL17_V5/AK4CHS/standard/QCDHT/RunBCDEF/");

  map_gr["UL18_A"] = CreateTGraphSF(path_+"eta_JER_tot/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunA/");
  map_gr["UL18_B"] = CreateTGraphSF(path_+"eta_JER_tot/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunB/");
  map_gr["UL18_C"] = CreateTGraphSF(path_+"eta_JER_tot/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunC/");
  map_gr["UL18_D"] = CreateTGraphSF(path_+"eta_JER_tot/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunD/");
  map_gr["UL18_ABCD"] = CreateTGraphSF(path_+"eta_JER_tot/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/");

  map_gr["UL18_A_eta"] = CreateTGraphSF(path_+"eta_simple/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunA/");
  map_gr["UL18_B_eta"] = CreateTGraphSF(path_+"eta_simple/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunB/");
  map_gr["UL18_C_eta"] = CreateTGraphSF(path_+"eta_simple/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunC/");
  map_gr["UL18_D_eta"] = CreateTGraphSF(path_+"eta_simple/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunD/");
  map_gr["UL18_ABCD_eta"] = CreateTGraphSF(path_+"eta_simple/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/");

  map_gr["UL18_ABCD_JER"] = CreateTGraphSF("/nfs/dust/cms/user/paaschal/CMSSW_106X_v1/CMSSW_10_6_13/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/eta_JER_original_temp/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/");
  map_gr["UL18_ABCD_common"] = CreateTGraphSF(path_+"eta_common_fine_aNew/UL18/Summer19UL18_V5/AK4CHS/standard/QCDHT/RunABCD/");

  for (auto [name,gr]: map_gr) gr->SetMarkerSize(0.7);

  int col2016 = kRed+1;
  int col2017 = kRed+1;
  int colUL17 = kGreen-1;
  int col2018 = kBlue+1;
  int colUL18 = kOrange+1;
  int colUL16 = kAzure+2;

  TString cName = "SF_final";
  TCanvas* canv_SF_final = tdrCanvas(cName, -0.5, 5.191+0.5, 0., 3.0, "|#eta|", "JER SF");
  TLegend *leg_final = tdrLeg(0.65,0.67,0.79,0.92, 0.035, 42, kBlack);
  tdrHeader(leg_final,"", 12);

  // tdrDraw(map_gr["EOY16"], "P5", kFullTriangleUp, col2016, kSolid, col2016, 1001, col2016, 0.15);
  // tdrDraw(map_gr["UL17"], "P5", kFullTriangleUp, colUL17, kSolid, colUL17, 1001, colUL17, 0.15);
  tdrDraw(map_gr["UL18_ABCD_common"], "P5", kFullTriangleUp, colUL18, kSolid, colUL18, 1001, colUL18, 0.15);
  // tdrDraw(map_gr["Fall17_V3"], "P5", kFullTriangleDown, col2017, kSolid, col2017, 1001, col2017, 0.15);
  // tdrDraw(map_gr["Autumn18_V7b"], "P5", kFullTriangleUp, col2018, kSolid, col2018, 1001, col2018, 0.15);
  // tdrDraw(map_gr["SFAutumn18_V7_RunABCD"], "P5", kFullSquare, kGreen-1, kSolid, kGreen-1, 1001, kGreen-1, 0.15);
  // tdrDraw(map_gr["SFAutumn18_V8"+DATA], "P5", kFullDotLarge, kGreen-1, kSolid, kGreen-1, 1001, kGreen-1, 0.15);
  // tdrDraw(map_gr["SFAutumn18_V8_AK4Puppi"+DATA], "P5", kFullDotLarge, kRed+1, kSolid, kRed+1, 3004, kRed+1, 0.15);
  // leg_final->AddEntry(map_gr["EOY16"], "EOY16","f");
  // leg_final->AddEntry(map_gr["EOY16"], "EOY16","f");
  // leg_final->AddEntry(map_gr["Fall17_V3"], "EOY17","f");
  // leg_final->AddEntry(map_gr["Autumn18_V7b"], "EOY18","f");
  // leg_final->AddEntry(map_gr["UL17"], "UL17","f");
  leg_final->AddEntry(map_gr["UL18_ABCD_common"], "UL18_ABCD","f");
  // leg_final->AddEntry(map_gr["SFAutumn18_V7_RunABCD"],"Autumn18_V7","f");
  // leg_final->AddEntry(map_gr["SFAutumn18_V7"+DATA],"JEC_V17_"+DATA,"f");
  //leg_final->AddEntry(map_gr["SFAutumn18_V8"+DATA],"JEC_V19_AK4CHS_"+DATA,"f");
  //leg_final->AddEntry(map_gr["SFAutumn18_V8_AK4Puppi"+DATA],"JEC_V19_AK4Puppi_"+DATA,"f");
  // leg_final->AddEntry(gr_final, "UL16postVFP","f");
  leg_final->Draw("same");

  canv_SF_final->SaveAs(outdir+cName+"_common_fine_aNew.pdf");


  lumi_13TeV = "RunII Legacy";

  PlotCanvas("SF_EOY", map_gr, {"EOY16","EOY17","EOY18"}, {kGreen-1,kRed+1,kBlue+1}, {kFullTriangleUp,kFullTriangleDown,kFullSquare});
  PlotCanvas("SF_UL", map_gr, {"UL16","UL17","UL18"}, {kGreen-1,kRed+1,kBlue+1}, {kFullTriangleUp,kFullTriangleDown,kFullSquare});
  PlotCanvas("SF_UL_new", map_gr, {"UL16_new","UL17_new","UL18_new"}, {kGreen-1,kRed+1,kBlue+1}, {kFullTriangleUp,kFullTriangleDown,kFullSquare});

  PlotCanvas("SF_16", map_gr, {"EOY16","UL16", "UL16_new"}, {kRed+1,kBlue+1,kGreen-1}, {kFullTriangleDown,kFullTriangleUp,kFullSquare});
  PlotCanvas("SF_17", map_gr, {"EOY17","UL17", "UL17_new"}, {kRed+1,kBlue+1,kGreen-1}, {kFullTriangleDown,kFullTriangleUp,kFullSquare});
  PlotCanvas("SF_18", map_gr, {"EOY18","UL18", "UL18_new"}, {kRed+1,kBlue+1,kGreen-1}, {kFullTriangleDown,kFullTriangleUp,kFullSquare});


  PlotCanvas("SF_16_eta", map_gr, {"UL16_new", "UL16_eta"}, {kRed+1, kOrange+1}, {kFullSquare,kFullCircle});
  PlotCanvas("SF_17_eta", map_gr, {"UL17_new", "UL17_eta"}, {kRed+1, kOrange+1}, {kFullSquare,kFullCircle});
  PlotCanvas("SF_18_eta", map_gr, {"UL18_new", "UL18_eta"}, {kRed+1, kOrange+1}, {kFullSquare,kFullCircle});


  PlotCanvas("SF_UL16", map_gr, {"EOY16","UL16preVFP", "UL16postVFP"}, {kRed+1,kBlue+1,kGreen-1}, {kFullTriangleDown,kFullTriangleUp,kFullSquare});
  PlotCanvas("SF_UL16_eta", map_gr, {"UL16preVFP", "UL16preVFP_eta"}, {kRed+1, kOrange+1}, {kFullSquare,kFullCircle});

  PlotCanvas("SF_RunII_eta", map_gr, {"RunII_new", "RunII_eta"}, {kRed+1, kOrange+1}, {kFullSquare,kFullCircle});

  PlotCanvas("SF_UL_RunII", map_gr, {"UL16_new","UL17_new","UL18_new","RunII_new"}, {kRed+1,kBlue+1,kGreen-1,kOrange+1}, {kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});
  PlotCanvas("SF_UL_RunII_eta", map_gr, {"UL16_eta","UL17_eta","UL18_eta","RunII_eta"}, {kRed+1,kBlue+1,kGreen-1,kOrange+1}, {kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});



  PlotCanvas("SF_UL_RunII_full",     map_gr, {"UL16preVFP", "UL16_new","UL17_new","UL18_new","RunII_new"}, {kAzure+2, kRed+1,kBlue+1,kGreen-1,kOrange+1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});
  PlotCanvas("SF_UL_RunII_eta_full", map_gr, {"UL16preVFP_eta", "UL16_eta","UL17_eta","UL18_eta","RunII_eta"}, {kAzure+2, kRed+1,kBlue+1,kGreen-1,kOrange+1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});

  PlotCanvas("minsuk", map_gr, {"UL16preVFP", "UL16postVFP","UL18"}, {kRed+1,kBlue+1,kGreen-1}, {kFullTriangleDown,kFullTriangleUp,kFullSquare});
  PlotCanvas("minsuk2", map_gr, {"UL16preVFP", "UL16postVFP","UL18_new"}, {kRed+1,kBlue+1,kGreen-1}, {kFullTriangleDown,kFullTriangleUp,kFullSquare});
  PlotCanvas("mikko", map_gr, {"UL16preVFP", "UL16postVFP","UL17_new", "UL18_new"}, {kRed+1,kBlue+1,kOrange+1,kGreen-1}, {kFullTriangleDown,kFullTriangleUp,kFullCircle, kFullSquare});


  PlotCanvas("UL18_run", map_gr, {"UL18_A", "UL18_B","UL18_C", "UL18_D"}, {kRed+1,kBlue+1,kGreen-1,kOrange+1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});
  PlotCanvas("UL18_run_all", map_gr, {"UL18_A", "UL18_B","UL18_C", "UL18_D", "UL18_ABCD"}, {kRed+1,kBlue+1,kGreen-1,kOrange+1,kAzure+2}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});

  PlotCanvas("UL18_eta_run", map_gr, {"UL18_A_eta", "UL18_B_eta","UL18_C_eta", "UL18_D_eta"}, {kRed+1,kBlue+1,kGreen-1,kOrange+1}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});
  PlotCanvas("UL18_eta_run_all", map_gr, {"UL18_A_eta", "UL18_B_eta","UL18_C_eta", "UL18_D_eta", "UL18_ABCD_eta"}, {kRed+1,kBlue+1,kGreen-1,kOrange+1,kAzure+2}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});


  PlotCanvas("UL17_run", map_gr, {"UL17_B", "UL17_C","UL17_D", "UL17_E", "UL17_F"}, {kRed+1,kBlue+1,kGreen-1,kOrange+1,kAzure+2}, {kFullStar, kFullTriangleDown,kFullTriangleUp,kFullSquare,kFullCircle});
  PlotCanvas("Compare_JER_common", map_gr, {"UL18_ABCD_common", "UL18_ABCD_JER"}, {kRed+1, kGreen-1}, {kFullCircle, kFullSquare});

}
