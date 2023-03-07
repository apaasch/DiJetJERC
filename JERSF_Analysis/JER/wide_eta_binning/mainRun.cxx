#include <iostream>
#include <iomanip>
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
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TString.h>

#include "TMinuit.h"
#include "TFitter.h"

#include "constants.h"
#include "functions.C"
#include "tdrstyle_all.h"

// double min_fit = 100.;
// double max_fit = 1200.;
double min_fit = 150.;
double max_fit = 2700.;
double min_fit_lin = 150.;
double max_fit_lin = 2700.;

bool useBoth = true; // NSxPC for eta <2.5 and NSC for eta >2.5
bool useP = false; // decides if default is NSC or NSxPC
bool useOriginal = false; // Fit before pt studies

// Code by Andrea Malara
// Based on code by Marek Niedziela, Matthias Schröder, Kristin Goebel

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void setParameters(TF1* func, std::vector<double> initial, std::vector<double> N, std::vector<double> S,std::vector<double> C,std::vector<double> P={0,0}){
  if(initial.size()==3){
    func->SetParameters(initial[0], initial[1], initial[2]);
    func->SetParLimits(0, N[0], N[1]);
    func->SetParLimits(1, S[0], S[1]);
    func->SetParLimits(2, C[0], C[1]);
  }
  else if(initial.size()==4){
    func->SetParameters(initial[0], initial[1], initial[2], initial[3]);
    func->SetParLimits(0, N[0], N[1]);
    func->SetParLimits(1, S[0], S[1]);
    func->SetParLimits(2, C[0], C[1]);
    func->SetParLimits(3, P[0], P[1]);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void DrawComparison(double eta, bool inRatio = false){
  TString name = to_string(eta);
  cout << name << endl;
  TF1* func = new TF1( "Comparison"+name, fNSxPC, 70, max_fit);
  if(eta==0.6525) func->SetParameters(2.738, 0.8313, 0.04231, -0.9279);
  if(eta==0.957)  func->SetParameters(3.152, 0.8212, 0.04697, -0.8962);
  if(eta==1.218)  func->SetParameters(3.531, 0.9372, 0.05745, -0.9065);
  if(eta==1.5225) func->SetParameters(3.589, 1.327, 0.06284, -1.05);
  if(eta==1.835)  func->SetParameters(3.705, 1.324, 0.05218, -1.051);
  if(eta==1.9865) func->SetParameters(3.989, 1.427, 0.05677, -1.084);
  if(eta==2.1825) func->SetParameters(5.831, 0.9054, 0.05601, -0.887);
  if(eta==2.411)  func->SetParameters(6.041, 1.288, 0.05329, -1.07);
  if(eta==2.575)  func->SetParameters(9.419, 0.2345, 4.562e-07, -0.4472);
  if(eta==2.7515) func->SetParameters(9.419, 0.2345, 4.562e-07, -0.4472);
  if(eta==2.9085) func->SetParameters(10.52, 0.1805, 5.12e-06, -0.1993);
  if(eta==3.0515) func->SetParameters(9.814, 0.00034, 0.1636, -1.239);
  if(eta==4.165)  func->SetParameters(5.616, 0.2159, 0.05147, -0.3009);
  func->SetLineColor(kGreen+2);
  func->Draw("SAME");

  char line[100];
  TLegend *legend;
  legend = tdrLeg(0.30,0.55,0.50,0.7, 0.025, 42, kGreen+2);
  legend->AddEntry((TObject*)0, "Mesurement", "");
  sprintf(line, "N = %.5f", func->GetParameter(0));  legend->AddEntry((TObject*)0, line, "");
  sprintf(line, "S = %.5f", func->GetParameter(1));  legend->AddEntry((TObject*)0, line, "");
  sprintf(line, "C = %.5f", func->GetParameter(2));  legend->AddEntry((TObject*)0, line, "");
  if(useP){ sprintf(line, "P = %.5f", func->GetParameter(3));  legend->AddEntry((TObject*)0, line, ""); }
  legend->Draw("same");
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void MCFIT(TH1* hist, TH1F* mcERR, double &N, double &S, double &C, double &P, double &Nerr, double &Serr, double &Cerr, double& Perr, double &mcChi, int &mcNDF, double eta) {
  TF1* mcFIT;
  if (useBoth){
    useP = eta<2.5?true:false;
    useOriginal = eta<2.5?false:true;
    if( (eta<2.5 && !useP) || (eta>2.5 && !useOriginal) ){
      cout << eta << " " << useP << " " << useOriginal << endl;
      throw runtime_error("Set of NSxPC and NSC funtions are wrong in MCFIT");
    }
  }
  if(useP)
  {
    mcFIT = new TF1( "mcFIT", fNSxPC, min_fit, max_fit);
    // if(eta==2.7515) mcFIT = new TF1( "mcFIT", fNSxPC, minfit, 430); // exclude last point
    // if(eta==2.9085) mcFIT = new TF1( "mcFIT", fNSxPC, 110, max_fit); // exclude first point

    setParameters(mcFIT, {5, 1, 0.05,-0.8}, {0.,20.}, {0.,2.}, {0.015,1.}, {-3.,0.});
    // setParameters(mcFIT, {5, 1, 0.05,-0.8}, {-20.,20.}, {0.,2.}, {0.015,1.}, {-3.,0.});
    if(eta>2.5) setParameters(mcFIT, {5, 1, 0.05,-0.8}, {0.,20.}, {0.,2.}, {0.015,1.}, {-3.,0.});
    // if(eta>2.5) setParameters(mcFIT, {5, 1, 0.05,-0.8}, {-20.,20.}, {0.,2.}, {0.015,1.}, {-3.,0.});
  }
  else if(useOriginal){
    mcFIT = new TF1( "mcFIT", fNSC, min_fit, max_fit);
    setParameters(mcFIT, {0.00015, 0.8, 0.04}, {0.,20.}, {0.,2.}, {0.015,1.});
  }
  else{
    mcFIT = new TF1( "mcFIT", fNSC, 70, max_fit);
    if(eta==2.7515) mcFIT = new TF1( "mcFIT", fNSC, min_fit, 430); // exclude last point
    if(eta==2.9085) mcFIT = new TF1( "mcFIT", fNSC, 110, max_fit); // exclude first point

    setParameters(mcFIT, {0.00015, 0.8, 0.04}, {0.,20.}, {0.,2.}, {0.015,1.});
    if(eta>2.5) setParameters(mcFIT, {0.00015, 0.8, 0.04}, {0.,20.}, {0.,2.}, {0.,1.});

  }

  hist-> Fit("mcFIT", "RMQ+"); // at last for CL/mcERR
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mcERR,0.68);

  N = mcFIT -> GetParameter( 0 );
  S = mcFIT -> GetParameter( 1 );
  C = mcFIT -> GetParameter( 2 );
  if(useP) P = mcFIT -> GetParameter( 3 );
  Nerr = mcFIT -> GetParError( 0 );
  Serr = mcFIT -> GetParError( 1 );
  Cerr = mcFIT -> GetParError( 2 );
  if(useP) Perr = mcFIT -> GetParError( 3 );
  mcChi = mcFIT->GetChisquare();
  mcNDF = mcFIT->GetNDF();
  mcFIT->SetLineColor(kBlue+2);
  mcERR->SetStats(kFALSE);
  mcERR->GetXaxis()->SetRange(min_fit,max_fit);

}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void DTFIT(TH1* hist, TH1F* dtERR, double &N, double &S, double &C, double &P, double &kNS, double &kC, double &Nerr, double &Serr, double &Cerr, double& Perr, double &kNSerr, double &kCerr, double &dtChi, int &dtNDF, double eta, int m, bool isFE) {
  TF1 * dtFIT;
  if (useBoth){
    useP = eta<2.5?true:false;
    useOriginal = eta<2.5?false:true;
    if( (eta<2.5 && !useP) || (eta>2.5 && !useOriginal) ){
      cout << eta << " " << useP << " " << useOriginal << endl;
      throw runtime_error("Set of NSxPC and NSC funtions are wrong in DTFIT");
    }
  }
  if(useP)
  {
    dtFIT = new TF1( "dtFIT", "TMath::Sqrt( ([3]*[3]*[0]*abs([0])/(x*x))+[3]*[3]*[1]*[1]*pow(x,[5])+[4]*[4]*[2]*[2] )", min_fit, max_fit);
    // if(eta==2.7515) dtFIT = new TF1( "mdtFIT", "TMath::Sqrt( ([3]*[3]*[0]*abs([0])/(x*x))+[3]*[3]*[1]*[1]*pow(x,[5])+[4]*[4]*[2]*[2] )", min_fit, 430); // exclude last point
    // if(eta==2.9085) dtFIT = new TF1( "mdtFIT", "TMath::Sqrt( ([3]*[3]*[0]*abs([0])/(x*x))+[3]*[3]*[1]*[1]*pow(x,[5])+[4]*[4]*[2]*[2] )", 110, max_fit); // exclude first point
  }
  else if(useOriginal)
  {
    dtFIT = new TF1( "dtFIT", "TMath::Sqrt( ([3]*[3]*[0]*abs([0])/(x*x))+[3]*[3]*[1]*[1]/x+[4]*[4]*[2]*[2] )", min_fit, max_fit);
  }
  else
  {
    dtFIT = new TF1( "dtFIT", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[3]*[3]*[1]*[1]/x+[4]*[4]*[2]*[2] )", min_fit, max_fit);
    if(eta==2.7515) dtFIT = new TF1( "mdtFIT", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[3]*[3]*[1]*[1]/x+[4]*[4]*[2]*[2] )", min_fit, 430); // exclude last point
    if(eta==2.9085) dtFIT = new TF1( "mdtFIT", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[3]*[3]*[1]*[1]/x+[4]*[4]*[2]*[2] )", 110, max_fit); // exclude first point
  }
  // TF1 * dtFIT = new TF1( "dtFIT", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[3]*[3]*[1]*[1]/x+[4]*[4]*[2]*[2] )", min_fit, max_fit);
  dtFIT -> FixParameter(0, N);
  dtFIT -> FixParameter(1, S);
  dtFIT -> FixParameter(2, C);
  if(useP) dtFIT -> FixParameter(5, P);
  dtFIT -> SetParameter(3, 1.2);
  dtFIT -> SetParameter(4, 1.2);
  dtFIT -> SetParLimits(3, 0.5, 3);
  dtFIT -> SetParLimits(4, 0.5, 3);
  // if (!isFE && m == 3) { dtFIT -> SetParLimits(3, 0.0, 15); dtFIT -> SetParLimits(4, 0.0, 15); }
  // if (!isFE && m == 3) { dtFIT -> SetParameter(0, N); dtFIT -> SetParameter(1, S); dtFIT -> SetParameter(2, C); dtFIT -> FixParameter(3,1); dtFIT -> FixParameter(4, 1); }
  if (hist->GetEntries() != 0) hist -> Fit("dtFIT", "RMQ+");
  else hist->GetListOfFunctions()->Add(dtFIT);

  kNS = dtFIT -> GetParameter( 3 );
  kC = dtFIT -> GetParameter( 4 );
  kNSerr = dtFIT -> GetParError( 3 );
  kCerr = dtFIT -> GetParError( 4 );
  dtChi = dtFIT->GetChisquare();
  dtNDF = dtFIT->GetNDF();
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(dtERR,0.68);
  dtERR->SetStats(kFALSE);
  dtERR->GetXaxis()->SetRange(min_fit,max_fit);
}

void DTFIT2(TH1* hist, TH1F* dtERR, double &N, double &S, double &C, double &P, double &kN, double &kS, double &kC, double &Nerr, double &Serr, double &Cerr, double& Perr, double &kNerr, double &kSerr, double &kCerr, double &dtChi, int &dtNDF, double eta, int m, bool isFE) {
  TF1 * dtFIT;
  if (useBoth){
    useP = eta<2.5?true:false;
    useOriginal = eta<2.5?false:true;
    if( (eta<2.5 && !useP) || (eta>2.5 && !useOriginal) ){
      cout << eta << " " << useP << " " << useOriginal << endl;
      throw runtime_error("Set of NSxPC and NSC funtions are wrong in DTFIT");
    }
  }
  if(useP)
  {
    dtFIT = new TF1( "dtFIT", "TMath::Sqrt( ([3]*[3]*[0]*abs([0])/(x*x))+[4]*[4]*[1]*[1]*pow(x,[6])+[5]*[5]*[2]*[2] )", min_fit, max_fit);
    // if(eta==2.7515) dtFIT = new TF1( "mdtFIT", "TMath::Sqrt( ([3]*[3]*[0]*abs([0])/(x*x))+[3]*[3]*[1]*[1]*pow(x,[5])+[4]*[4]*[2]*[2] )", min_fit, 430); // exclude last point
    // if(eta==2.9085) dtFIT = new TF1( "mdtFIT", "TMath::Sqrt( ([3]*[3]*[0]*abs([0])/(x*x))+[3]*[3]*[1]*[1]*pow(x,[5])+[4]*[4]*[2]*[2] )", 110, max_fit); // exclude first point
  }
  else if(useOriginal)
  {
    dtFIT = new TF1( "dtFIT", "TMath::Sqrt( ([3]*[3]*[0]*abs([0])/(x*x))+[4]*[4]*[1]*[1]/x+[5]*[5]*[2]*[2] )", min_fit, max_fit);
  }
  else
  {
    dtFIT = new TF1( "dtFIT", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[4]*[4]*[1]*[1]/x+[5]*[5]*[2]*[2] )", min_fit, max_fit);
    if(eta==2.7515) dtFIT = new TF1( "mdtFIT", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[4]*[4]*[1]*[1]/x+[5]*[5]*[2]*[2] )", min_fit, 430); // exclude last point
    if(eta==2.9085) dtFIT = new TF1( "mdtFIT", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[4]*[4]*[1]*[1]/x+[5]*[5]*[2]*[2] )", 110, max_fit); // exclude first point
  }
  // TF1 * dtFIT = new TF1( "dtFIT", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[3]*[3]*[1]*[1]/x+[4]*[4]*[2]*[2] )", min_fit, max_fit);
  dtFIT -> FixParameter(0, N);
  dtFIT -> FixParameter(1, S);
  dtFIT -> FixParameter(2, C);
  if(useP) dtFIT -> FixParameter(6, P);
  dtFIT -> SetParameter(3, 1.2);
  dtFIT -> SetParameter(4, 1.2);
  dtFIT -> SetParameter(5, 1.2);
  dtFIT -> SetParLimits(3, 0.5, 3);
  dtFIT -> SetParLimits(4, 0.5, 3);
  dtFIT -> SetParLimits(5, 0.5, 3);
  // if (!isFE && m == 3) { dtFIT -> SetParLimits(3, 0.0, 15); dtFIT -> SetParLimits(4, 0.0, 15); }
  // if (!isFE && m == 3) { dtFIT -> SetParameter(0, N); dtFIT -> SetParameter(1, S); dtFIT -> SetParameter(2, C); dtFIT -> FixParameter(3,1); dtFIT -> FixParameter(4, 1); }
  if (hist->GetEntries() != 0) hist -> Fit("dtFIT", "RMQ+");
  else hist->GetListOfFunctions()->Add(dtFIT);

  kN = dtFIT -> GetParameter( 3 );
  kS = dtFIT -> GetParameter( 4 );
  kC = dtFIT -> GetParameter( 5 );
  kNerr = dtFIT -> GetParError( 3 );
  kSerr = dtFIT -> GetParError( 4 );
  kCerr = dtFIT -> GetParError( 5 );
  dtChi = dtFIT->GetChisquare();
  dtNDF = dtFIT->GetNDF();
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(dtERR,0.68);
  dtERR->SetStats(kFALSE);
  dtERR->GetXaxis()->SetRange(min_fit,max_fit);
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PLOT_ALLWIDTHS(std::vector< std::vector< TH1F* > > h_data, std::vector< std::vector< TH1F* > > h_MC, std::vector< std::vector< TH1F* > > h_gen, TString outdir) {
  int color_data = kBlue;
  int color_MC = kRed;
  int color_gen = kGreen+2;
  for( unsigned int m = 0; m < h_data.size(); m++ ){
    for( unsigned int p = 0; p < h_data.at(m).size(); p++ ){
      TString canvName  = h_data.at(m).at(p)->GetTitle();
      TString nameXaxis = h_data.at(m).at(p)->GetXaxis()->GetTitle();
      TString nameYaxis = h_data.at(m).at(p)->GetYaxis()->GetTitle();
      std::vector<TH1*> vec;
      vec.push_back(h_data.at(m).at(p));
      vec.push_back(h_MC.at(m).at(p));
      vec.push_back(h_gen.at(m).at(p));
      double x_min, x_max, y_min, y_max;
      findExtreme2(vec, &x_min, &x_max, &y_min, &y_max);
      TCanvas* canv = tdrCanvas(canvName, 0, 0.35, 0, y_max*1.6, nameXaxis, nameYaxis);
      canv->SetTickx(0);
      canv->SetTicky(0);
      tdrDraw(h_data.at(m).at(p), "", kFullCircle, color_data );
      tdrDraw(h_MC.at(m).at(p),   "", kFullCircle, color_MC );
      tdrDraw(h_gen.at(m).at(p),  "", kFullCircle, color_gen );
      canv -> Update();
      canv -> Print(outdir+canvName+".pdf","pdf");
      delete canv;
      // canv -> Print(outdir+canvName+".png","png");
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PLOT_MCT(std::vector< TH1F* > h_MCTruth, std::vector< TH1F* > h_uncor, std::vector< TH1F* > h_cor, std::vector< TH1F* > h_015 , TString outdir, std::vector<double> eta_bins, bool isFE) {
  int MCTruth      = kMagenta;
  int color_uncor  = kRed;
  int color_cor    = kBlue;
  int color_015    = kGreen+2;
  for( unsigned int m = 0; m < h_cor.size(); m++ ){
    TString canvName  = h_cor.at(m)->GetTitle();
    TString nameXaxis = h_cor.at(m)->GetXaxis()->GetTitle();
    TString nameYaxis = h_cor.at(m)->GetYaxis()->GetTitle();
    nameXaxis = "p_{T} [GeV]";
    std::vector<TH1*> vec;

    if (isFE) vec.push_back(h_MCTruth.at(m+1));
    else vec.push_back(h_MCTruth.at(m));
    vec.push_back(h_MCTruth.at(m));
    vec.push_back(h_uncor.at(m));
    vec.push_back(h_cor.at(m));
    vec.push_back(h_015.at(m));
    double x_min, x_max, y_min, y_max;
    findExtreme(vec, &x_min, &x_max, &y_min, &y_max);
    TCanvas* canv = tdrCanvas(canvName, x_min, x_max, y_min, y_max, nameXaxis, nameYaxis);
    canv->SetTickx(0);
    canv->SetTicky(0);
    tdrDraw(h_MCTruth.at(m),  "", kFullCircle, MCTruth );
    tdrDraw(h_uncor.at(m),    "", kFullCircle, color_uncor );
    tdrDraw(h_cor.at(m),      "", kFullCircle, color_cor );
    tdrDraw(h_015.at(m),      "", kFullCircle, color_015 );

    char legTitle[100];
    TLegend *leg = tdrLeg(0.6,0.7,0.9,0.9);
    sprintf(legTitle,     "#eta #in [%.3f,%.3f]", eta_bins[m], eta_bins[m+1]);
    leg->AddEntry(h_MCTruth.at(m),  "MC Truth",         "lep");
    leg->AddEntry(h_uncor.at(m),    "MC Uncorrelated",  "lep");
    leg->AddEntry(h_cor.at(m),      "MC Correlated",    "lep");
    leg->AddEntry(h_015.at(m),      "MC 015",           "lep");

    if(isFE){
      if (m<9) canvName = canvName(0, canvName.Length()-1);
      else canvName = canvName(0, canvName.Length()-2);
      canvName += (m+2);
    }
    canv->Print(outdir+canvName+".pdf","pdf");
    delete canv;
    // canv->Print(outdir+canvName+".png","png");
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PLOT_ASY(std::vector< std::vector< std::vector< TH1F* > > > h_data, std::vector< std::vector< std::vector< TH1F* > > > h_MC, std::vector< std::vector< std::vector< TH1F* > > > h_gen, std::vector< std::vector< std::vector< double > > > h_data_width, std::vector< std::vector< std::vector< double > > > h_MC_width, std::vector< std::vector< std::vector< double > > > h_gen_width, std::vector< std::vector< std::vector< double > > > h_data_width_err, std::vector< std::vector< std::vector< double > > > h_MC_width_err, std::vector< std::vector< std::vector< double > > > h_gen_width_err, TString outdir, std::vector<double> eta_bins, std::vector<double> Pt_bins_Central, std::vector<double> Pt_bins_HF, std::vector<double> alpha, std::vector< std::vector< std::vector< double > > > lower_x_data, std::vector< std::vector< std::vector< double > > > upper_x_data, std::vector< std::vector< std::vector< double > > > lower_x, std::vector< std::vector< std::vector< double > > > upper_x, std::vector< std::vector< std::vector< double > > > gen_lower_x, std::vector< std::vector< std::vector< double > > > gen_upper_x) {
  int color_data = kBlue;
  int color_MC = kRed;
  int color_gen = kGreen+2;
  for( unsigned int m = 0; m < h_data.size(); m++ ){
    // if(m!=3) continue;
    std::vector<double> Pt_bins;
    if (eta_bins[m]< eta_cut) for (size_t i = 0; i <= h_data.at(m).size(); i++) { Pt_bins.push_back(Pt_bins_Central[i]); }
    else for (size_t i = 0; i <= h_data.at(m).size(); i++) { Pt_bins.push_back(Pt_bins_HF[i]); }
    for( unsigned int p = 0; p < h_data.at(m).size(); p++ ){
      // if(p!=9) continue;
      for( unsigned int r = 0; r < h_data.at(m).at(p).size(); r++ ){
        // if(r!=2) continue;
        h_data.at(m).at(p).at(r)-> Scale(1./h_data.at(m).at(p).at(r) -> Integral());
        h_MC.at(m).at(p).at(r)  -> Scale(1./h_MC.at(m).at(p).at(r) -> Integral());
        h_gen.at(m).at(p).at(r) -> Scale(1./h_gen.at(m).at(p).at(r) -> Integral());

        TString canvName  = h_data.at(m).at(p).at(r)->GetTitle();
        TString nameXaxis = "Asymmetry"; // h_data.at(m).at(p).at(r)->GetXaxis()->GetTitle();
        TString nameYaxis = "a.u."; // h_data.at(m).at(p).at(r)->GetYaxis()->GetTitle();
        std::vector<TH1*> vec;
        vec.push_back(h_data.at(m).at(p).at(r));
        vec.push_back(h_MC.at(m).at(p).at(r));
        vec.push_back(h_gen.at(m).at(p).at(r));
        double x_min, x_max, y_min, y_max;
        findExtreme(vec, &x_min, &x_max, &y_min, &y_max);

        extraText3.clear();
        extraText3.push_back(Form("%d GeV < p_{T}^{ave} < %d GeV", (int)Pt_bins[p], (int)Pt_bins[p+1]));
        extraText3.push_back(Form("%.1f < |#eta| < %.1f", eta_bins[m], eta_bins[m+1]));
        extraText3.push_back(Form("#alpha_{max} < %.2f", alpha.at(r)));

        TCanvas* canv = tdrCanvas(canvName, -0.5, 0.5, 0.00001, 100, nameXaxis, nameYaxis);
        canv->SetLogy();

        tdrDraw(h_data.at(m).at(p).at(r), "", kFullCircle, color_data );
        tdrDraw(h_MC.at(m).at(p).at(r),   "", kFullCircle, color_MC );
        tdrDraw(h_gen.at(m).at(p).at(r),  "", kFullCircle, color_gen );

        if (h_data.at(m).at(p).at(r)!=0 && h_data.at(m).at(p).at(r)->GetEntries()>50) {
          TF1* gaus_data = new TF1("gauss_data", "gaus", -0.6, 0.6); TF1* gaus_data2 = new TF1("gauss_data2", "gaus", lower_x_data.at(m).at(p).at(r), upper_x_data.at(m).at(p).at(r));
          gaus_data->SetLineColor(color_data); gaus_data->SetLineStyle(kDashed); gaus_data2->SetLineColor(color_data); gaus_data2->SetLineStyle(kSolid);
          h_data.at(m).at(p).at(r)->Fit("gauss_data", "MQ", "", lower_x_data.at(m).at(p).at(r), upper_x_data.at(m).at(p).at(r));
          h_data.at(m).at(p).at(r)->Fit("gauss_data2", "MQ", "", lower_x_data.at(m).at(p).at(r), upper_x_data.at(m).at(p).at(r));
          gaus_data->Draw("same"); gaus_data2->Draw("same");
        }

        if (h_MC.at(m).at(p).at(r)!=0 && h_MC.at(m).at(p).at(r)->GetEntries()>50) {
          TF1* gaus_mc = new TF1("gauss_MC", "gaus", -0.6, 0.6); TF1* gaus_mc2 = new TF1("gauss_MC2", "gaus", lower_x.at(m).at(p).at(r), upper_x.at(m).at(p).at(r));
          gaus_mc->SetLineColor(color_MC); gaus_mc->SetLineStyle(kDashed); gaus_mc2->SetLineColor(color_MC); gaus_mc2->SetLineStyle(kSolid);
          h_MC.at(m).at(p).at(r)->Fit("gauss_MC", "MQ", "", lower_x.at(m).at(p).at(r), upper_x.at(m).at(p).at(r));
          h_MC.at(m).at(p).at(r)->Fit("gauss_MC2", "MQ", "", lower_x.at(m).at(p).at(r), upper_x.at(m).at(p).at(r));
          gaus_mc->Draw("same"); gaus_mc2->Draw("same");
        }

        if (h_gen.at(m).at(p).at(r)!=0 && h_gen.at(m).at(p).at(r)->GetEntries()>50) {
          TF1* gaus_gen = new TF1("gauss_gen", "gaus", -0.6, 0.6); TF1* gaus_gen2 = new TF1("gauss_gen2", "gaus", gen_lower_x.at(m).at(p).at(r), gen_upper_x.at(m).at(p).at(r));
          gaus_gen->SetLineColor(color_gen); gaus_gen->SetLineStyle(kDashed); gaus_gen2->SetLineColor(color_gen); gaus_gen2->SetLineStyle(kSolid);
          h_gen.at(m).at(p).at(r)->Fit("gauss_gen", "MQ", "", gen_lower_x.at(m).at(p).at(r), gen_upper_x.at(m).at(p).at(r));
          h_gen.at(m).at(p).at(r)->Fit("gauss_gen2", "MQ", "", gen_lower_x.at(m).at(p).at(r), gen_upper_x.at(m).at(p).at(r));
          gaus_gen->Draw("same"); gaus_gen2->Draw("same");
        }

        TLine* line_lower_data = new TLine(lower_x_data.at(m).at(p).at(r), 0.00001, lower_x_data.at(m).at(p).at(r), 1);
        TLine* line_upper_data = new TLine(upper_x_data.at(m).at(p).at(r), 0.00001, upper_x_data.at(m).at(p).at(r), 1);
        TLine* line_lower_MC = new TLine(lower_x.at(m).at(p).at(r), 0.00001, lower_x.at(m).at(p).at(r), 1);
        TLine* line_upper_MC = new TLine(upper_x.at(m).at(p).at(r), 0.00001, upper_x.at(m).at(p).at(r), 1);
        TLine* line_lower_gen = new TLine(gen_lower_x.at(m).at(p).at(r), 0.00001, gen_lower_x.at(m).at(p).at(r), 1);
        TLine* line_upper_gen = new TLine(gen_upper_x.at(m).at(p).at(r), 0.00001, gen_upper_x.at(m).at(p).at(r), 1);

        line_lower_data->SetLineWidth(1); line_lower_data->SetLineColor(color_data); line_lower_data->Draw("same");
        line_upper_data->SetLineWidth(1); line_upper_data->SetLineColor(color_data); line_upper_data->Draw("same");
        line_lower_MC->SetLineWidth(1); line_lower_MC->SetLineColor(color_MC); line_lower_MC->Draw("same");
        line_upper_MC->SetLineWidth(1); line_upper_MC->SetLineColor(color_MC); line_upper_MC->Draw("same");
        line_lower_gen->SetLineWidth(1); line_lower_gen->SetLineColor(color_gen); line_lower_gen->Draw("same");
        line_upper_gen->SetLineWidth(1); line_upper_gen->SetLineColor(color_gen); line_upper_gen->Draw("same");

        TLegend *leg = tdrLeg(0.45,0.7,0.95,0.9);
        leg->AddEntry(h_data.at(m).at(p).at(r), Form("#sigma_{A}^{data} = %.4f #pm %.4f", h_data_width.at(m).at(p).at(r), h_data_width_err.at(m).at(p).at(r)), "lep");
        leg->AddEntry(h_MC.at(m).at(p).at(r),   Form("#sigma_{A}^{MC}  = %.4f #pm %.4f",  h_MC_width.at(m).at(p).at(r),   h_MC_width_err.at(m).at(p).at(r)),   "lep");
        leg->AddEntry(h_gen.at(m).at(p).at(r),  Form("#sigma_{A}^{gen}  = %.4f #pm %.4f", h_gen_width.at(m).at(p).at(r),  h_gen_width_err.at(m).at(p).at(r)),  "lep");
        leg->Draw("same");

        canv->Update();
        canv -> Print(outdir+canvName+".pdf","pdf");
        delete canv;
        // canv -> Print(outdir+canvName+".png","png");
      }
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PLOT_WIDTH(std::vector< std::vector< TH1F* > > h_data, std::vector< std::vector< TH1F* > > h_MC, std::vector< std::vector< TH1F* > > h_gen, TString outdir, std::vector<double> eta_bins, std::vector<double> Pt_bins_Central, std::vector<double> Pt_bins_HF) {
  int color_data = kBlue;
  int color_MC = kRed;
  int color_gen = kGreen+2;
  for( unsigned int m = 0; m < h_data.size(); m++ ){
    std::vector<double> Pt_bins;
    if (eta_bins[m]< eta_cut) for (size_t i = 0; i <= h_data.at(m).size(); i++) { Pt_bins.push_back(Pt_bins_Central[i]); }
    else for (size_t i = 0; i <= h_data.at(m).size(); i++) { Pt_bins.push_back(Pt_bins_HF[i]); }
    for( unsigned int p = 0; p < h_data.at(m).size(); p++ ){
      TString canvName  = h_data.at(m).at(p)->GetTitle();
      TString nameXaxis = h_data.at(m).at(p)->GetXaxis()->GetTitle();
      TString nameYaxis = h_data.at(m).at(p)->GetYaxis()->GetTitle();
      std::vector<TH1*> vec;
      vec.push_back(h_data.at(m).at(p));
      vec.push_back(h_MC.at(m).at(p));
      vec.push_back(h_gen.at(m).at(p));
      double x_min, x_max, y_min, y_max;
      findExtreme2(vec, &x_min, &x_max, &y_min, &y_max);
      TCanvas* canv = tdrCanvas(canvName, 0, 0.35, 0, y_max*1.6, nameXaxis, nameYaxis);

      tdrDraw(h_data.at(m).at(p), "", kFullCircle, color_data );
      tdrDraw(h_MC.at(m).at(p),   "", kFullCircle, color_MC );
      tdrDraw(h_gen.at(m).at(p),  "", kFullCircle, color_gen );

      if( h_MC.at(m).at(p)  -> GetEntries() != 0 ) h_MC.at(m).at(p)   -> GetFunction("linfit")->SetLineColor(color_MC);
      if( h_gen.at(m).at(p) -> GetEntries() != 0 ) h_gen.at(m).at(p)  -> GetFunction("linfit")->SetLineColor(color_gen);
      if( h_data.at(m).at(p)-> GetEntries() != 0 ) h_data.at(m).at(p) -> GetFunction("linfit")->SetLineColor(color_data);

      TF1* f;
      char line[100];
      TLegend *legend;

      if (h_MC.at(m).at(p) -> GetFunction("linfit")) {
        f = h_MC.at(m).at(p) -> GetFunction("linfit");
        legend = tdrLeg(0.50,0.70,0.70,0.85, 0.025, 42, color_MC);
        tdrHeader(legend,"MC", 22);
        sprintf(line, "#chi^{2}/ndf = %.2f/%d", f->GetChisquare(),  f->GetNDF());       legend->AddEntry((TObject*)0, line, "");
        sprintf(line, "p0 = %.5f #pm %.5f",     f->GetParameter(0), f->GetParError(0)); legend->AddEntry((TObject*)0, line, "");
        sprintf(line, "p1 = %.5f #pm %.5f",     f->GetParameter(1), f->GetParError(1)); legend->AddEntry((TObject*)0, line, "");
        legend->Draw("same");
      }

      if (h_gen.at(m).at(p) -> GetFunction("linfit")) {
        f = h_gen.at(m).at(p) -> GetFunction("linfit");
        legend = tdrLeg(0.70,0.70,0.90,0.85, 0.025, 42, color_gen);
        tdrHeader(legend,"gen", 22);
        sprintf(line, "#chi^{2}/ndf = %.2f/%d", f->GetChisquare(), f->GetNDF());      legend->AddEntry((TObject*)0, line, "");
        sprintf(line, "p0 = %.5f #pm %.5f", f->GetParameter(0), f->GetParError(0));   legend->AddEntry((TObject*)0, line, "");
        sprintf(line, "p1 = %.5f #pm %.5f", f->GetParameter(1), f->GetParError(1));    legend->AddEntry((TObject*)0, line, "");
        legend->Draw("same");
      }

      if (h_data.at(m).at(p) -> GetFunction("linfit")) {
        f = h_data.at(m).at(p) -> GetFunction("linfit");
        legend = tdrLeg(0.30,0.70,0.50,0.85, 0.025, 42, color_data);
        tdrHeader(legend,"Data", 22);
        sprintf(line, "#chi^{2}/ndf = %.2f/%d", f->GetChisquare(), f->GetNDF());      legend->AddEntry((TObject*)0, line, "");
        sprintf(line, "p0 = %.5f #pm %.5f", f->GetParameter(0), f->GetParError(0));   legend->AddEntry((TObject*)0, line, "");
        sprintf(line, "p1 = %.5f #pm %.5f", f->GetParameter(1), f->GetParError(1));    legend->AddEntry((TObject*)0, line, "");
        legend->Draw("same");
      }

      TLegend *leg = tdrLeg(0.35,0.85,0.9,0.89, 0.025, 42, kBlack);
      tdrHeader(leg, Form("#eta #in [%.3f,%.3f], p_{T} #in [%.0f,%.0f] GeV", eta_bins[m], eta_bins[m+1], Pt_bins[p], Pt_bins[p+1]));
      leg->Draw("same");

      canv -> Update();
      canv -> Print(outdir+canvName+".pdf","pdf");
      delete canv;
      // canv -> Print(outdir+canvName+".png","png");
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PLOT_WIDTH_gr(std::vector< std::vector< TGraphErrors* > > h_data, std::vector< std::vector< TGraphErrors* > > h_MC, std::vector< std::vector< TGraphErrors* > > h_gen, TString outdir, std::vector<double> eta_bins, std::vector<double> Pt_bins_Central, std::vector<double> Pt_bins_HF, bool isFE) {
  int color_data = kBlue;
  int color_MC = kRed;
  int color_gen = kGreen+2;

  for( unsigned int m = 0; m < h_data.size(); m++ ){
    std::vector<double> Pt_bins;
    if (eta_bins[m]< eta_cut) for (size_t i = 0; i <= h_data.at(m).size(); i++) { Pt_bins.push_back(Pt_bins_Central[i]); }
    else for (size_t i = 0; i <= h_data.at(m).size(); i++) { Pt_bins.push_back(Pt_bins_HF[i]); }

    for( unsigned int p = 0; p < h_data.at(m).size(); p++ ){

      if(extenda){
        if( m<4 && !extendAlpha(p) ) continue;
      }

      double range = std::max(0.05,removePointsforAlphaExtrapolation(isFE, eta_bins.at(m), p+1));

      TString canvName  = h_data.at(m).at(p)->GetTitle();
      TString nameXaxis = "#alpha_{max}";
      TString nameYaxis = "#sigma_{A}";
      std::vector<TGraphErrors*> vec;
      double x_min, x_max, y_min, y_max;
      vec.push_back(h_data.at(m).at(p));
      vec.push_back(h_MC.at(m).at(p));
      vec.push_back(h_gen.at(m).at(p));
      findExtreme2(vec, &x_min, &x_max, &y_min, &y_max);
      extraText3.clear();
      extraText3.push_back(Form("%d GeV < p_{T}^{ave} < %d GeV", (int)Pt_bins[p], (int)Pt_bins[p+1]));
      extraText3.push_back(Form("%.1f < |#eta| < %.1f", eta_bins[m], eta_bins[m+1]));

      double amax = 0.35;
      if(outdir.Contains("highalpha")) amax = 1.05;
      TCanvas* canv = tdrCanvas(canvName, 0, amax, 0, y_max*1.6, nameXaxis, nameYaxis);
      canv->SetGrid();

      tdrDraw(h_data.at(m).at(p), "p ", kFullCircle, color_data );
      tdrDraw(h_MC.at(m).at(p),   "p ", kFullCircle, color_MC );
      tdrDraw(h_gen.at(m).at(p),  "p ", kFullCircle, color_gen );

      TF1* f;
      TF1* Temp = new TF1();
      char line[100];
      TLegend *legend;

      if (h_MC.at(m).at(p) -> GetFunction("lin_extrapol_mc")) {
        f = h_MC.at(m).at(p) -> GetFunction("lin_extrapol_mc"); f->SetLineStyle(2); f->SetLineColor(color_MC);
        Temp = (TF1*) h_MC.at(m).at(p)->GetFunction("lin_extrapol_mc")->Clone(); gStyle -> SetOptFit(0000);
        Temp->SetRange(range,0.33); Temp->SetLineStyle(1); Temp->Draw("same");

        if(alpha_new) f = h_MC.at(m).at(p) -> GetFunction("lin_extrapol_mc_first"); f->SetLineStyle(3); f->SetLineColor(color_MC);
        if(alpha_new) f = h_MC.at(m).at(p) -> GetFunction("lin_extrapol_mc_last"); f->SetLineStyle(4); f->SetLineColor(color_MC);
      }

      if (h_gen.at(m).at(p) -> GetFunction("lin_extrapol_mc")) {

        f = h_gen.at(m).at(p) -> GetFunction("lin_extrapol_mc"); f->SetLineStyle(2); f->SetLineColor(color_gen);
        Temp = (TF1*) h_gen.at(m).at(p)->GetFunction("lin_extrapol_mc")->Clone(); gStyle -> SetOptFit(0000);
        Temp->SetRange(range,0.33); Temp->SetLineStyle(1); Temp->Draw("same");

        if(alpha_new) f = h_gen.at(m).at(p) -> GetFunction("lin_extrapol_mc_first"); f->SetLineStyle(3); f->SetLineColor(color_gen);
        if(alpha_new) f = h_gen.at(m).at(p) -> GetFunction("lin_extrapol_mc_last"); f->SetLineStyle(4); f->SetLineColor(color_gen);

        legend = tdrLeg(0.70,0.70,0.90,0.85, 0.025, 42, color_gen);
      }

      if (h_data.at(m).at(p) -> GetFunction("lin_extrapol_mc")) {

        f = h_data.at(m).at(p) -> GetFunction("lin_extrapol_mc"); f->SetLineStyle(2); f->SetLineColor(color_data);
        Temp = (TF1*) h_data.at(m).at(p)->GetFunction("lin_extrapol_mc")->Clone(); gStyle -> SetOptFit(0000);
        Temp->SetRange(range,0.33); Temp->SetLineStyle(1); Temp->Draw("same");

        if(alpha_new) f = h_data.at(m).at(p) -> GetFunction("lin_extrapol_mc_first"); f->SetLineStyle(3); f->SetLineColor(color_data);
        if(alpha_new) f = h_data.at(m).at(p) -> GetFunction("lin_extrapol_mc_last"); f->SetLineStyle(4); f->SetLineColor(color_data);
      }

      TLegend *leg = tdrLeg(0.50,0.70,0.85,0.90, 0.045, 42, kBlack);

      TLine* first = new TLine(0, 0, 0, 0); first->SetLineColor(kGray+1); first->SetLineStyle(3);
      TLine* last = new TLine(0, 0, 0, 0); last->SetLineColor(kGray+1); last->SetLineStyle(4);

      if(alpha_new) leg->SetNColumns(2);
      leg->AddEntry(h_data.at(m).at(p), "data", "lp"); if(alpha_new) leg->AddEntry(first, "remove first point", "l");
      leg->AddEntry(h_MC.at(m).at(p), "MC", "lp"); if(alpha_new) leg->AddEntry(last, "remove last point", "l");
      leg->AddEntry(h_gen.at(m).at(p), "gen", "lp");
      leg->Draw("same");

      TString add_name = "Extrapol_";
      if (isFE) add_name = "Extrapol_FE_";

      canv->Print(outdir+add_name+Form("Eta%i_pt%i.pdf", m+1, p+1));
      delete canv;
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void SFtoTXT(std::ofstream& texfile, std::vector< TH1F* > h_JER, std::vector< std::vector< std::vector< double > > > width_pt, std::vector<double> &SF, std::vector<double> &SF_err, std::vector<double> &SF_ptdep_min, std::vector<double> &SF_ptdep_max, std::vector<double> &eta_bin_center, std::vector<double> &eta_bin_err, std::vector<double> eta_bins, int shift, bool isFE, bool isCorr ) {
  for( unsigned int m = 0; m < h_JER.size(); m++ ){
    TF1 * constfit = h_JER.at(m) -> GetFunction("constfit");
    TF1 * NSC_ratio = h_JER.at(m) -> GetFunction("NSC_ratio");
    h_JER.at(m)->GetFunction("NSC_ratio")->SetBit(TF1::kNotDraw);
    h_JER.at(m) -> Write();
    int diff = 0;
    if (isFE) diff = shift - h_JER.size();
    if (!isCorr) {
      eta_bin_center.push_back((eta_bins[diff+m+1]+eta_bins[diff+m]) /2);
      eta_bin_err.push_back((eta_bins[diff+m+1]-eta_bins[diff+m]) /2);
    }
    SF_ptdep_min.push_back(findMinMax(h_JER.at(m), width_pt.at(m), NSC_ratio, constfit, 1));
    SF_ptdep_max.push_back(findMinMax(h_JER.at(m), width_pt.at(m), NSC_ratio, constfit, 0));
    SF.push_back(constfit -> GetParameter( 0 ));
    SF_err.push_back(constfit -> GetParError( 0 ));
    texfile << constfit -> GetParameter( 0 ) << " \\pm " << constfit -> GetParError( 0 ) << " & ";
    if(m == h_JER.size()-1 ) texfile << "\\\\";
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void SIGtoTXT(std::vector< TH1F* > h_MC, std::vector<double> eta_bins, TString outdir, TString name){
  ofstream output_values;
  output_values.open(outdir+"pdfy/JERs/values_"+name+".txt");
  output_values << right << fixed << setprecision(5);
  output_values << setw(7) << "eta[" << " " << setw(7) << "eta]" << " ";
  output_values << setw(12) << "pT"    << " " << setw(8) << "sigma\n";
  for( unsigned int m = 0; m < h_MC.size(); m++ ){
    for(unsigned int bin = 1; bin <= h_MC.at(m)->GetNbinsX(); bin++){
      double pt    = h_MC.at(m)->GetXaxis()->GetBinCenter(bin);
      double width = h_MC.at(m)->GetBinContent(bin);
      if(width!=0){
        output_values << setw(7) << eta_bins[m] << " " << setw(7) << eta_bins[m+1] << " ";
        output_values << setw(12) << pt         << " " << setw(7) << width         << "\n";
      }
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PLOT_SF(std::vector< TH1F* > h_uncor, std::vector< TH1F* > h_cor, std::vector< TH1F* > h_015, std::vector< TH1F* > h_fitERR, TString outdir, std::vector<double> eta_bins, bool isFE, std::vector< TH1F* > h_JER_corr_NSC) {
  if(debug) cout << "In PLOT_SF " << ((isFE)?"FE":"SM") << endl;
  int color_uncor = kBlue;
  int color_cor   = kRed;
  int color_015   = kGreen+2;
  for( unsigned int m = 0; m < h_uncor.size(); m++ ){
    TString canvName  = h_uncor.at(m)->GetTitle();
    TString nameXaxis = "p_{T}^{ave} [GeV]";
    TString nameYaxis = h_uncor.at(m)->GetYaxis()->GetTitle();
    std::vector<TH1*> vec;
    vec.push_back(h_uncor.at(m));
    vec.push_back(h_cor.at(m));
    vec.push_back(h_015.at(m));
    double x_min, x_max, y_min, y_max;
    findExtreme(vec, &x_min, &x_max, &y_min, &y_max);
    extraText3.clear();
    extraText3.push_back(Form("%.1f < |#eta| < %.1f", eta_bins[m], eta_bins[m+1]));

    TCanvas* canv = tdrCanvas(canvName, 50, 2700, 0.99, 1.4, nameXaxis, nameYaxis);

    tdrDraw(h_uncor.at(m), "", kFullCircle, color_uncor );
    tdrDraw(h_cor.at(m), "", kFullCircle, color_cor );
    tdrDraw(h_015.at(m), "", kFullCircle, color_015 );

    TF1* f;
    char line[100];
    TLegend *legend;

    // if(debug) cout << "In before fitERR " << m << " " << h_fitERR.size() << " " << h_cor.size() << endl;
    h_fitERR.at(m)->Draw("E4 same");

    if(h_uncor.at(m)->GetFunction("constfit")){
      f = h_uncor.at(m) -> GetFunction("constfit"); f->SetLineColor(color_uncor);f->SetLineWidth(2);f->Draw("same");
      f->SetParameter(0,f->GetParameter(0)-0.002);
    } else { std::cout << "Fit uncor function at bin " << m << " was not found" << std::endl; }

    if(h_cor.at(m)->GetFunction("constfit")){
      f = h_cor.at(m) -> GetFunction("constfit"); f->SetLineColor(color_cor);f->SetLineWidth(2);f->Draw("same");
    } else { std::cout << "Fit cor function at bin " << m << " was not found" << std::endl; }

    f = new TF1( "constfit", "pol0", min_fit_lin, max_fit_lin );
    f->SetParameter(0,1);
    if (h_015.at(m)->GetEntries() > 1) h_015.at(m) -> Fit("constfit","RMQ+");
    else h_015.at(m)->GetListOfFunctions()->Add(f);
    if (h_015.at(m) -> GetFunction("constfit")==0) h_015.at(m)->GetListOfFunctions()->Add(f);
    f->SetLineColor(color_015); f->SetLineWidth(2);f->Draw("same");

    f = h_JER_corr_NSC.at(m) -> GetFunction("NSC_ratio");
    f->SetLineColor(kOrange+1); f->SetLineWidth(2); f->Draw("same");


    TLegend *leg = tdrLeg(0.55,0.55,0.9,0.85, 0.04, 42, kBlack);
    leg->AddEntry(h_uncor.at(m),  "#sigma_{JER}^{data}/#sigma_{JER}^{MC} uncorrelated","lep");
    leg->AddEntry(h_cor.at(m),    "#sigma_{JER}^{data}/#sigma_{JER}^{MC} correlated","lep");
    leg->AddEntry(h_015.at(m),    "#sigma_{0.15}^{data}/#sigma_{0.15}^{MC}","lep");
    leg->AddEntry(f, "NSC ratio", "l");
    leg->Draw("same");

    if (isFE) {
      if (m<9) canvName = canvName(0, canvName.Length()-1);
      else canvName = canvName(0, canvName.Length()-2);
      canvName += (m+2);
    }
    canv -> Print(outdir+canvName+".pdf","pdf");
    delete canv;
    // canv -> Print(outdir+canvName+".png","png");
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PrintFitInfo(TF1* fit, TString function, std::vector<double> range,  std::vector<double> initial, std::vector<double> N, std::vector<double> S, std::vector<double> C, std::vector<double> P={0,0}, TString addition=""){
  cout << "Draw " << addition << function << " in [" << setw(4) << dtos(range[0], 0) << ", " << setw(4) << dtos(range[1], 0) << "] with ";
  cout << "N in [" << setw(5) << dtos(N[0], 2) << ", " << setw(4) << dtos(N[1], 2) << "] (N0: " << setw(4) << dtos(initial[0], 2) << ") " << "and ";
  cout << "S in [" << setw(5) << dtos(S[0], 2) << ", " << setw(4) << dtos(S[1], 2) << "] (S0: " << setw(4) << dtos(initial[1], 2) << ") " << "and ";
  cout << "C in [" << setw(5) << dtos(C[0], 2) << ", " << setw(4) << dtos(C[1], 2) << "] (C0: " << setw(4) << dtos(initial[2], 2) << ") ";
  if(function.Contains("NSxPC")) cout << "and P in [" << setw(5) << dtos(P[0], 2) << ", " << setw(4) << dtos(P[1], 2) << "] (P0: " << setw(4) << dtos(initial[3], 2) << ") ";
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void DrawLegFit(TF1* fit, TLegend* leg, double vN, double vNerr, double vS, double vSerr, double vC, double vCerr, double vCh, int vNDF, int color, TString line1, TString line2, int place){
  bool isP = line1.Contains("NSxPC");
  double shift = (isP)?0.025:0.;
  char line[100];
  if(place<0 || 5<place) cerr << "You selected the wrong legend position for the ratio plots" << endl;
  if(place==0) leg = tdrLeg(0.30,0.65-shift,0.45,0.8,  0.025, 42, color);
  if(place==1) leg = tdrLeg(0.45,0.65-shift,0.65,0.8,  0.025, 42, color);
  if(place==2) leg = tdrLeg(0.65,0.65-shift,0.80,0.8,  0.025, 42, color);
  if(place==3) leg = tdrLeg(0.45,0.47-shift,0.65,0.62, 0.025, 42, color);
  if(place==4) leg = tdrLeg(0.65,0.47-shift,0.80,0.62, 0.025, 42, color);
  if(place==5) leg = tdrLeg(0.30,0.47-shift,0.45,0.62, 0.025, 42, color);
  leg->AddEntry((TObject*)0, line1, "");
  leg->AddEntry((TObject*)0, line2, "");
  vN  = fit->GetParameter(0); vNerr = fit->GetParError(0);
  vS  = fit->GetParameter(1); vSerr = fit->GetParError(1);
  vC  = fit->GetParameter(2); vCerr = fit->GetParError(2);
  vCh = fit->GetChisquare();  vNDF  = fit->GetNDF();
  sprintf(line, "#chi^{2}/ndf = %.2f/%d", vCh,vNDF); leg->AddEntry((TObject*)0, line, "");
  sprintf(line, "N = %.5f #pm %.5f", vN, vNerr);     leg->AddEntry((TObject*)0, line, "");
  sprintf(line, "S = %.5f #pm %.5f", vS, vSerr);     leg->AddEntry((TObject*)0, line, "");
  sprintf(line, "C = %.5f #pm %.5f", vC, vCerr);     leg->AddEntry((TObject*)0, line, "");
  if(isP){sprintf(line, "P = %.5f #pm %.5f", fit->GetParameter(3), fit->GetParError(3)); leg->AddEntry((TObject*)0, line, "");}
  leg->Draw("same");

}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void AddFit(bool draw, int color, int place, TCanvas* canv, TString name, TString function, TString info1,  TString info2, TH1F* hist, double r1, double r2, std::vector<double> initial, std::vector<double> N, std::vector<double> S,std::vector<double> C,std::vector<double> P={0,0}){
  TF1* fit= new TF1( name, function, r1, r2);
  setParameters(fit, initial, N, S, C, P);
  hist->Fit(name, "0RMQ+");

  if(draw)
  {
    canv->cd(1);
    TF1* func = hist->GetFunction(name);
    func->SetLineColor(color);
    func->Draw("same");

    double vN, vNerr, vS, vSerr, vC, vCerr, vCh;
    int vNDF;
    char line[100];
    TLegend *rleg;
    DrawLegFit(func, rleg, vN, vNerr, vS, vSerr, vC, vCerr, vCh, vNDF, color, info1, info2, place); // Place==0 for default value

    canv->cd(2);

    TF1* ratio = GetFitsRatio(func, hist->GetFunction("mcFIT"), r1, color, function);
    ratio->Draw("SAME");
  }
}

void AddFit(bool draw, int color, int place, TCanvas* canv, TString name, TString function, TString info1,  TString info2, TGraphErrors* hist, TF1* func2, double r1, double r2, std::vector<double> initial, std::vector<double> N, std::vector<double> S,std::vector<double> C,std::vector<double> P={0,0}){
  TF1* fit= new TF1( name, function, r1, r2);
  setParameters(fit, initial, N, S, C, P);
  hist->Fit(name, "0RMQ+");

  if(draw)
  {
    canv->cd(1);
    TF1* func = hist->GetFunction(name);
    func->SetLineColor(color);
    func->Draw("same");
    double vN, vNerr, vS, vSerr, vC, vCerr, vCh;
    int vNDF;
    char line[100];
    TLegend *rleg;
    DrawLegFit(func, rleg, vN, vNerr, vS, vSerr, vC, vCerr, vCh, vNDF, color, info1, info2, place); // Place==0 for default value
    canv->cd(2);
    TF1* ratio = GetFitsRatio(func, func2, r1, color, function);
    ratio->Draw("SAME");
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PLOT_NSC(std::vector< TH1F* > h_data, std::vector< TH1F* > h_MC, std::vector< TH1F* > h_SF, std::vector< TH1F* > &h_fitERR, TString outdir, std::vector<double> eta_bins, bool isFE, bool isCorr) {

  int color_data  = kBlue;
  int color_MC    = kRed;
  int color_NSC   = kGreen+2;

  TString method, corr;
  if(isFE) method = "FE";
  else     method = "SM";
  if(isCorr) corr = "correlated";
  else       corr = "uncorrelated";
  ofstream output_fit;
  output_fit.open(outdir+"pdfy/kValues/"+corr+"_"+method+".txt");

  output_fit << right;
  output_fit << setw(7) << "eta[" << " " << setw(7) << "eta]" << " ";
  output_fit << setw(7) << "N" << " " << setw(7) << "Nerr" << " " << setw(7) << "S" << " " << setw(7) << "Serr" << " ";
  output_fit << setw(7) << "C" << " " << setw(7) << "Cerr" << " ";
  if(useP || useBoth) output_fit << setw(7) << "P" << " " << setw(7) << "Perr" << " ";
  output_fit << setw(7) << "kNS"  << " " << setw(7) << "kNSerr" << " " << setw(7) << "kC" << " " << setw(8) << "kCerr\n";

  ofstream FitPar;
  FitPar.open(outdir+"pdfy/kValues/MultipleFitPar.txt");
  FitPar << fixed << setprecision(5) << right;

  for( unsigned int m = 0; m < h_data.size(); m++ ){
    double eta = (eta_bins[m]+eta_bins[m+1])/2;
    // if(m==1) RemoveAsymBinsFromJER(h_MC.at(m), 5);
    h_data.at(m)->SetStats(kFALSE);
    h_MC.at(m)->SetStats(kFALSE);
    // double N, S, C, P, kN, kNS, Nerr, Serr, Cerr, Perr, kNSerr, kCerr, mcChi, dtChi;
    double N, S, C, P, kN, kS, kC, Nerr, Serr, Cerr, Perr, kNerr, kSerr, kCerr, mcChi, dtChi;
    int mcNDF, dtNDF;

    TH1F* mcERR = (TH1F*)h_data.at(m)->Clone();
    TH1F* dtERR = (TH1F*)h_data.at(m)->Clone();
    MCFIT(h_MC.at(m), mcERR, N, S, C, P, Nerr, Serr, Cerr, Perr, mcChi, mcNDF, eta);
    // DTFIT(h_data.at(m), dtERR, N, S, C, P, kNS, kC, Nerr, Serr, Cerr, Perr, kNSerr, kCerr, dtChi, dtNDF, eta, m, isFE);
    DTFIT2(h_data.at(m), dtERR, N, S, C, P, kN, kS, kC, Nerr, Serr, Cerr, Perr, kNerr, kSerr, kCerr, dtChi, dtNDF, eta, m, isFE);

    output_fit << fixed << setprecision(5);
    output_fit << setw(7) << eta_bins[m] << " " << setw(7) << eta_bins[m+1] << " ";
    output_fit << setw(7) << N << " " << setw(7) << Nerr << " " << setw(7) << S    << " " << setw(7) << Serr   << " ";
    output_fit << setw(7) << C << " " << setw(7) << Cerr << " ";
    if(useP) output_fit << setw(7) << P << " " << setw(7) << Perr << " ";
    // output_fit << setw(7) << kNS  << " " << setw(7) << kNSerr << " " << setw(7) << kC          << " " << setw(7) << kCerr         << "\n";
    output_fit << setw(7) << kN  << " " << setw(7) << kNerr << setw(7) << kS  << " " << setw(7) << kSerr << " " << setw(7) << kC          << " " << setw(7) << kCerr         << "\n";

    if (h_data.at(m)-> GetFunction("dtFIT")!=0) h_data.at(m)-> GetFunction("dtFIT")->SetLineColor(color_data+2);
    if (h_MC.at(m)-> GetFunction("mcFIT")!=0) h_MC.at(m)-> GetFunction("mcFIT")->SetLineColor(color_MC+2);
    mcERR->SetFillColorAlpha(color_MC+2,0.35);
    dtERR->SetFillColorAlpha(color_data+2,0.35);

    TString canvName  = h_data.at(m)->GetTitle();
    TString nameXaxis = "p_{T}^{ave} [GeV]"; // h_data.at(m)->GetXaxis()->GetTitle();
    TString nameYaxis = h_data.at(m)->GetYaxis()->GetTitle();

    std::vector<TH1*> vec;
    vec.push_back(h_data.at(m));
    vec.push_back(h_MC.at(m));

    double x_min, x_max, y_min, y_max;
    findExtreme(vec, &x_min, &x_max, &y_min, &y_max);
    extraText3.clear();
    extraText3.push_back(Form("%.1f < |#eta| < %.1f", eta_bins[m], eta_bins[m+1]));

    // TCanvas* canv = tdrCanvas(canvName, x_min, x_max, y_min, y_max, nameXaxis, nameYaxis);
    TCanvas* canv = tdrCanvas(canvName, 49, x_max, y_min, y_max, nameXaxis, nameYaxis);
    canv->SetLogx(1);

    mcERR->Draw("E4 SAME");
    dtERR->Draw("E4 SAME");
    tdrDraw(h_data.at(m), "", kFullCircle, color_data );
    tdrDraw(h_MC.at(m), "", kFullCircle, color_MC );
    TLegend *leg = tdrLeg(0.6,0.7,0.9,0.9);
    leg->SetNColumns(2);
    leg->AddEntry(h_MC.at(m),"MC","lp");
    leg->AddEntry(h_data.at(m),"data","lp");
    char line[100];
    TLegend *legend;
    legend = tdrLeg(0.50,0.55,0.70,0.7, 0.025, 42, color_MC);
    sprintf(line, "N = %.5f #pm %.5f", N, Nerr);  legend->AddEntry((TObject*)0, line, "");
    sprintf(line, "S = %.5f #pm %.5f", S, Serr);  legend->AddEntry((TObject*)0, line, "");
    sprintf(line, "C = %.5f #pm %.5f", C, Cerr);  legend->AddEntry((TObject*)0, line, "");
    if(useP){ sprintf(line, "P = %.5f #pm %.5f", P, Perr);  legend->AddEntry((TObject*)0, line, ""); }
    legend->Draw("same");
    legend = tdrLeg(0.70,0.55,0.85,0.7, 0.025, 42, color_data);
    // sprintf(line, "k_{NS} = %.5f #pm %.5f", kNS,kNSerr);legend->AddEntry((TObject*)0, line, "");
    sprintf(line, "k_{N} = %.5f #pm %.5f", kN, kNerr); legend->AddEntry((TObject*)0, line, "");
    sprintf(line, "k_{S} = %.5f #pm %.5f", kS, kSerr); legend->AddEntry((TObject*)0, line, "");
    sprintf(line, "k_{C} = %.5f #pm %.5f", kC, kCerr); legend->AddEntry((TObject*)0, line, "");
    legend->AddEntry((TObject*)0, "", "");
    legend->Draw("same");
    if(isFE){
      if (m<9) canvName = canvName(0, canvName.Length()-1);
      else canvName = canvName(0, canvName.Length()-2);
      canvName += (m+2);
    }
    canv -> Print(outdir+"pdfy/JERs/"+canvName+".pdf","pdf");

    if(isFE&&isCorr){
      // if(debug) cout << "Start comparing plots" << endl;
      TF1* mcFIT = h_MC.at(m)->GetFunction("mcFIT"); mcFIT->SetLineColor(color_MC+2);

      bool drawEXT = false; int colorEXT = kGray+2;
      bool drawIND = false; int colorIND = kCyan+1;

      TGraphErrors* MC_graph = TH1toTGraphError(h_MC.at(m), 0);
      MC_graph->GetListOfFunctions()->Add(mcFIT);
      MC_graph->SetMarkerColor(colorIND);
      MC_graph->SetLineColor(colorIND);

      TGraphErrors* error_ext = TH1toTGraphError(h_MC.at(m), ((eta_bins[m+1]>2.5)?0.001:0.001));
      MC_graph->GetListOfFunctions()->Add(mcFIT);
      error_ext->SetLineColor(colorEXT);

      // Exclude points in eta/pt-bins: 5/2 6/2 7/2 8/2,6 10/11 11/9 12/2,3 | Bin starts at 0 for TGraph
      // Bins are inserted by Hand. Way to do it automatically?
      double x,y; MC_graph->GetPoint(0, x, y);
      // if(m==8||m==9||m==10) MC_graph->RemovePoint(0);// For bin 10-13
      // if(m==3)  MC_graph->RemovePoint(1);
      // if(m==4)  MC_graph->RemovePoint(1);
      // if(m==5)  MC_graph->RemovePoint(1);
      // if(m==6){ MC_graph->RemovePoint(1); MC_graph->RemovePoint(4);}
      // if(m==8)  MC_graph->RemovePoint(10);
      // if(m==9)  MC_graph->RemovePoint(9);
      // if(m==10) MC_graph->RemovePoint(0);

      // if(m==11) MC-Graph->RemovePoint(9); // eta 2.6 < 2.9 and pt 440 ?
      if(m==10) MC_graph->RemovePoint(2);  // eta 2.9 < 3.0 and pt 440 ?
      // if(m==12) MC_graph->RemovePoint(9); // eta 3.1 < 5.2 and pt 440 ?

      TCanvas*               canv2 = tdrDiCanvas2(canvName, x_min, x_max, 0., 0.2, 0.85, 1.15, "p_{T} [GeV]", nameYaxis,"Ratio", false, 4, 11);
      if(eta_bins[m]==2.964) canv2 = tdrDiCanvas2(canvName, x_min, x_max, 0., 0.3, 0.85, 1.15, "p_{T} [GeV]", nameYaxis,"Ratio", false, 4, 11);
      canv2->SetTickx(0);
      canv2->SetTicky(0);
      canv2->cd(1); // JER
      tdrDraw(h_MC.at(m), "", kFullCircle, color_MC );
      mcERR->Draw("E4 same");
      // mcFIT->Draw("SAME");

      if(drawEXT) error_ext->Draw("ERROR SAME");
      if(drawEXT) h_MC.at(m)->Draw("ERROR SAME");
      if(drawIND) MC_graph->Draw("P SAME");

      // if(debug) cout << "Start Legend for Ratio" << endl;
      double vN, vNerr, vS, vSerr, vC, vCerr, vCh;
      int vNDF;
      char line[100];
      TLegend *rleg;
      rleg = tdrLeg(0.45,0.8,0.65,0.9);

      rleg->Draw("same");
      TString default_text = "Default; "+(TString) (useP?"NSxPC":"NSC");
      DrawLegFit(mcFIT, rleg, vN, vNerr, vS, vSerr, vC, vCerr, vCh, vNDF, color_MC, default_text, "pT in [70, 1200] GeV",  1); // Place==0 for default value

      canv2->cd(2); // Ratio Plots
      TH1F* ratio_MC  = GetMCFitRatio(h_MC.at(m), mcFIT, color_MC, 1);
      TLine* fit_default = new TLine(0, 1, 1100, 1);
      TH1F* ratio_err = (TH1F*) mcERR->Clone();

      for(unsigned int i=1; i<=mcERR->GetNbinsX(); i++){
        double bin_content = ratio_err->GetBinContent(i);
        double bin_error = ratio_err->GetBinError(i);
        double error = bin_error/bin_content;
        ratio_err->SetBinContent(i, 1);
        ratio_err->SetBinError(i, error);
      }
      ratio_err->SetMarkerSize(0);

      // add_function
      fit_default->SetLineColor(color_MC+2);
      fit_default->Draw();
      ratio_err->Draw("E4 SAME");
      ratio_MC->Draw("P SAME");

      // Add new fit functions and compare to default. For pT-dependent SFs studies
      // void AddFit(bool draw, color, place, canv, TString name, TString function, info1,  info2, hist, r1, r2, initial, N, S, C, P={0,0})
      AddFit(false, kOrange+2,  1, canv2, "mcFITold",    fNSC,      "NSC; old",           "pT in [70, 1200] GeV",  h_MC.at(m), 70,  1200, {0.00015, 0.8, 0.04},    {0.,10.},   {0.,2.}, {0.,1.});
      AddFit(false, kGreen+2,   3, canv2, "mcFIT100",    fNSC,      "NSC; ",              "pT in [100, 1200] GeV", h_MC.at(m), 100, 1200, {0.00015, 0.8, 0.04},    {0.,10.},   {0.,2.}, {0.,1.});
      AddFit(false, kGreen+2,   3, canv2, "mcFIT200",    fNSC,      "NSC; ",              "pT in [200, 1200] GeV", h_MC.at(m), 200, 1200, {0.00015, 0.8, 0.04},    {0.,10.},   {0.,2.}, {0.,1.});
      AddFit(false, kGray+2,    4, canv2, "mcFITv1",     fNSC,      "NSC; v1",            "pT in [70, 1200] GeV",  h_MC.at(m), 70,  1200, {5, 0.8, 0.04},          {0.,20.},   {0.,2.}, {0.02,1.});
      AddFit(false, kBlue+2,    9, canv2, "mcFITsig",    fNSCsign,  "sign(N)NSC; ",       "pT in [70, 1200] GeV",  h_MC.at(m), 70,  1200, {0.00015, 0.8, 0.04},    {-10.,10.}, {0.,2.}, {0.,1.});
      AddFit(false, kMagenta+2, 9, canv2, "mcFITpowCut", fNSxPC,    "NSxPC; ",            "pT in [110, 1200] GeV", h_MC.at(m), 110, 1200, {0.00015, 0.8, 0.06,-1}, {-20.,20.}, {0.,2.}, {0.02,1.}, {-3.,0.});
      AddFit(false, kOrange+2,  2, canv2, "mcFITpow",    fNSxPC,    "NSxPC; ",            "pT in [70, 1200] GeV",  h_MC.at(m), 70,  1200, {0.00015, 0.8, 0.06,-1}, {-20.,20.}, {0.,6.}, {0.015,1.}, {-3.,0.});
      AddFit(false, kOrange+2,  3, canv2, "mcFITpowV1",  fNSxPC,    "NSxPC; new initial", "pT in [70, 1200] GeV",  h_MC.at(m), 70,  1200, {5, 1, 0.05,-0.8},       {-20.,20.}, {0.,2.}, {0.015,1.}, {-3.,0.});
      AddFit(false, kOrange+2,  4, canv2, "mcFITpowN",   fNSxPC,    "NSxPC; new N limit", "pT in [70, 1200] GeV",  h_MC.at(m), 70,  1200, {1, 1, 0.05,-1.0},       {-5.,5.},   {0.,2.}, {0.015,1.}, {-3.,0.});
      AddFit(false, kCyan+2,    1, canv2, "mcFITpowdiv", fNSxPCdiv, "NSxPC; C*C/100",     "pT in [70, 1200] GeV",  h_MC.at(m), 70,  1200, {5, 1, 5,-0.8},          {-20.,20.}, {0.,2.}, {0.015,10.}, {-3.,0.});

      // void AddFit(bool draw, color, place, canv, TString name, TString function, info1,  info2, graph, ratio, r1, r2, initial, N, S, C, P={0,0}){
      AddFit(drawIND, kCyan+1, 3, canv2, "mcFITind", fNSC, "NSC; ", "Individual", MC_graph, mcFIT, 70, 1200, {0.00015, 0.8, 0.04}, {0.,10.}, {0.,2.}, {0.02,1.});
      AddFit(drawIND, kBlue+2, 2, canv2, "mcFITindpow", fNSxPC, "NSxPC; ", "Individual", MC_graph, mcFIT, 70, 1200, {0.5, 0.8, 0.04, -1}, {-20.,20.}, {0.,2.}, {0.02,1.}, {-3.,0.});
      AddFit(drawEXT, colorEXT, 4, canv2, "mcFIText", fNSC, "NSC; ", (TString) "syst. error = "+((eta_bins[m+1]>2.5)?"0.001":"0.001"), error_ext, mcFIT, 70, 1200, {0.00015, 0.8, 0.04}, {0.,10.}, {0.,2.}, {0.,1.});


      canv2 -> cd(0);
      canv2 -> Print(outdir+"pdfy/JERs/ratio_"+canvName+".pdf","pdf");
      delete canv2;

      if(isFE&&isCorr){ // add_function
        // if(debug) cout << "Start extracting FitParameters" << endl;
        FitPar << "====================\n";
        FitPar << "eta in ["+dtos(eta_bins[m],3)+","+dtos(eta_bins[m+1],3)+"]\n";
        FitPar << "====================\n";
        TString info = useP?"NSxPC":"NSC";
        ExtractFitParameters(FitPar, h_MC.at(m)->GetFunction("mcFIT"), info+"; 70-1200 GeV ");
        ExtractFitParameters(FitPar, h_MC.at(m)->GetFunction("mcFITv1"), "NSC; v1; 70-1200 GeV ");
        ExtractFitParameters(FitPar, h_MC.at(m)->GetFunction("mcFIT200"), "NSC; 200-1200 GeV ");
        ExtractFitParameters(FitPar, h_MC.at(m)->GetFunction("mcFIT100"), "NSC; 100-1200 GeV ");
        ExtractFitParameters(FitPar, h_MC.at(m)->GetFunction("mcFITsig"), "sign(N)NSC; 70-1200 GeV ");
        ExtractFitParameters(FitPar, h_MC.at(m)->GetFunction("mcFITold"), "NSC; old; 70-1200 GeV ");
        ExtractFitParameters(FitPar, h_MC.at(m)->GetFunction("mcFITpow"), "NSxPC; 70-1200 GeV ");
        ExtractFitParameters(FitPar, h_MC.at(m)->GetFunction("mcFITpowCut"), "NSxPC; 100-1200 GeV ");
        ExtractFitParameters(FitPar, MC_graph->GetFunction("mcFITind"), "NSC; Individual ");
        ExtractFitParameters(FitPar, MC_graph->GetFunction("mcFITindpow"), "NSxPC; Individual ");
        ExtractFitParameters(FitPar, error_ext->GetFunction("mcFIText"), "NSC; 70-1200 GeV; large errors");
      }
    }

    canvName  = h_SF.at(m)->GetTitle();
    nameXaxis = h_SF.at(m)->GetXaxis()->GetTitle();
    nameYaxis = h_SF.at(m)->GetYaxis()->GetTitle();
    vec.clear();
    vec.push_back(h_SF.at(m));
    findExtreme(vec, &x_min, &x_max, &y_min, &y_max);
    extraText3.clear();
    extraText3.push_back(Form("%.1f < |#eta| < %.1f", eta_bins[m], eta_bins[m+1]));
    canv = tdrCanvas(canvName, 150, 800, 0.9, 1.4, nameXaxis, nameYaxis);
    tdrDraw(h_SF.at(m), "", kFullCircle, color_data );

    TF1 * constfit = new TF1( "constfit", "pol0", min_fit_lin, max_fit_lin );
    constfit->SetParameter(0,1);
    constfit->SetLineColor(color_MC);
    constfit->SetLineWidth(2);
    TH1F* fitERR = (TH1F*)h_SF.at(m)->Clone();
    if (h_SF.at(m)->GetEntries() > 1){
      h_SF.at(m) -> Fit("constfit","RMQ+");
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(fitERR,0.68);
      fitERR->SetStats(kFALSE);
      fitERR->GetXaxis()->SetRange(min_fit_lin,max_fit_lin);
      fitERR->SetFillColorAlpha(kRed,0.35); // defined in PLOT_SF color_cor
      h_fitERR.push_back(fitERR);
    }
    else h_SF.at(m)->GetListOfFunctions()->Add(constfit);
    if (h_SF.at(m) -> GetFunction("constfit")==0) h_SF.at(m)->GetListOfFunctions()->Add(constfit);

    // TF1 *NSC_ratio = new TF1("NSC_ratio", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[3]*[3]*[1]*[1]/x+[4]*[4]*[2]*[2] )/TMath::Sqrt( ([0]*[0]/(x*x))+[1]*[1]/x+[2]*[2] )",min_fit,max_fit);
    // NSC_ratio -> SetParameters(N, S, C, kNS, kC);
    TF1 *NSC_ratio = new TF1("NSC_ratio", "TMath::Sqrt( ([3]*[3]*[0]*[0]/(x*x))+[4]*[4]*[1]*[1]/x+[5]*[5]*[2]*[2] )/TMath::Sqrt( ([0]*[0]/(x*x))+[1]*[1]/x+[2]*[2] )",min_fit,max_fit);
    NSC_ratio -> SetParameters(N, S, C, kN, kS, kC);
    NSC_ratio->SetLineColor(color_NSC);
    NSC_ratio->SetLineWidth(2);
    NSC_ratio->Draw("same");
    h_SF.at(m)->GetListOfFunctions()->Add(NSC_ratio);
    leg = tdrLeg(0.55,0.7,0.9,0.9);
    leg->AddEntry(h_SF.at(m),"#sigma_{JER}^{data}/#sigma_{JER}^{mc} correlated","lep");
    leg->AddEntry(constfit,"Constant Fit","l");
    leg->AddEntry(NSC_ratio,"Ratio of NSC-Fits","l");
    if(isFE){
      if (m<9) canvName = canvName(0, canvName.Length()-1);
      else canvName = canvName(0, canvName.Length()-2);
      canvName += (m+2);
    }
    canv -> Print(outdir+"pdfy/NSC_SFs/NSC"+canvName+".pdf","pdf");
    // canv -> Print(outdir+"pdfy/NSC_SFs/NSC"+canvName+".png","png");
    canv->Write();
    delete canv;
  }
}

// ------------------------------
//          MAIN PROGRAM
// ------------------------------
//
// bool data_ = false;
// bool real_data = true;
// const char* filename = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/JERSF/Analysis/hist_preparation/MC/wide_eta_bin/file/Single/Fall17_17Nov2017_V10/AK4CHS/histograms_mc_incl_full.root"
// const char* filename_data = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/JERSF/Analysis/hist_preparation/data/wide_eta_bin/file/Single/Fall17_17Nov2017_V10/AK4CHS/RunBCDEF/histograms_data_incl_full.root"
//
// TString Trigger = "Single"
//
// double gaustails = 0.985
// float shiftForPLI = 0.0
// int ref_shift = 3
// int shift = ref_shift

int mainRun(std::string year, bool data_, const char* filename, const char* filename_data, TString lumi, TString label_mc, TString label_data, TString Trigger, TString outdir, double gaustails = 0.985, float shiftForPLI = 0.0, int ref_shift = 3){

  time_t start, end; // TIME
  time_t startp, endp; // TIME
  double time_taken;

  time(&startp);

  gPrintViaErrorHandler = kTRUE;
  gErrorIgnoreLevel = kFatal;
  double hist_max_value = 0.3;
  double hist_max_value_SF = 3.0;

  g_year = year;
  g_study = Trigger;

  if(useBoth) {useP=false; useOriginal=false;};
  if(useP) {useBoth=false; useOriginal=false;};
  if(useOriginal) {useBoth=false; useP=false;};

  ////////////////////////////////////////////////////////////////////////////
  //    I load all histograms I will need                                   //
  ////////////////////////////////////////////////////////////////////////////

  bool real_data = true; // If false: compare MC to MC:
  // read out genlevel info from both MC and data.
  // If true: compare MC to data:
  // use MC gen info for data as well.
  //bool data = false;     // MC is not data and gets this boolean.

  if ( !real_data && (label_data == "Data" || label_data == "data")){
    std::cout << "WARNING: are you using data?" << std::endl;
    std::cout << "Program has set real_data=" << real_data << std::endl;
    std::cout << "PLI will be corrected with MC gen info only: " << real_data <<std::endl;
  }

  // ------------------------------
  //           bin values
  // ------------------------------

  if(Trigger.Contains("eta_common")) ref_shift = 5;
  if(Trigger.Contains("eta_calo")) ref_shift = 15;
  if(debug) cout << "Reference shift is " << ref_shift << endl;

  bool isAK8 = outdir.Contains("AK8");
  if (debug) std::cout << outdir << "\t" << isAK8 << "\t" << year << "\n";

  std::vector<double> eta_bins;

  if (Trigger.Contains("eta_narrow"))      {binHF=11; eta_bins = std::vector<double>(eta_bins_narrow, eta_bins_narrow + n_eta_bins_narrow);}
  else if (Trigger.Contains("eta_common")) {binHF=15; eta_bins = std::vector<double>(eta_bins_common, eta_bins_common + n_eta_bins_common);}
  else if (Trigger.Contains("eta_calo"))   {binHF=29; eta_bins = std::vector<double>(eta_bins_calo, eta_bins_calo + n_eta_bins_calo);}
  else if (Trigger.Contains("eta_simple")) {binHF=11; eta_bins = std::vector<double>(eta_bins_simple, eta_bins_simple + n_eta_bins_simple);}
  else if (Trigger.Contains("eta_L2R"))    {binHF=11; eta_bins = std::vector<double>(eta_bins_L2R, eta_bins_L2R + n_eta_bins_L2R);}
  else                                     {binHF=11; eta_bins = std::vector<double>(eta_bins_JER, eta_bins_JER + n_eta_bins_JER);}

  int n_eta_bins = eta_bins.size();

  int EtaBins_SM            = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return i<eta_cut; });; // st method bins
  int EtaBins_SM_control    = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return i>eta_cut; });; // st method bins control
  int EtaBins_FE_reference  = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return i<s_eta_barr;});; // fe method bins reference
  int EtaBins_FE_control    = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return (i>=s_eta_barr)&&(i<eta_cut);});; // fe method bins control
  int EtaBins_FE            = std::count_if(&eta_bins[0], &eta_bins[0]+n_eta_bins, [](double i) { return i>eta_cut; });; // fe method bins

  if (debug) {
    std::cout << "EtaBins_SM " << EtaBins_SM << '\n';
    std::cout << "EtaBins_SM_control " << EtaBins_SM_control << '\n';
    std::cout << "EtaBins_FE_reference " << EtaBins_FE_reference << '\n';
    std::cout << "EtaBins_FE_control " << EtaBins_FE_control << '\n';
    std::cout << "EtaBins_FE " << EtaBins_FE << '\n';
  }

  int shift_ = EtaBins_SM + EtaBins_FE;// TODO recheck
  if (debug) std::cout << "shift_: " << shift_ << '\n';

  int etaShift_SM           = 0;
  int etaShift_SM_control   = EtaBins_SM;
  int etaShift_FE_reference = 0;
  int etaShift_FE_control   = EtaBins_FE_reference;
  int etaShift_FE           = EtaBins_FE_reference + EtaBins_FE_control;

  if(debug) cout << "General: " << etaShift_FE << "\t | Control: " << etaShift_FE_control << "\t | Ref.: " << etaShift_FE_reference << endl;

  int PtBins_Central = 9, PtBins_HF = 6, AlphaBins = 6;

  std::vector<double> Pt_bins_Central, Pt_bins_HF;

  std::string triggerName = isAK8? "SingleJet" : "DiJet";
  std::string name_pt_bin = triggerName+"_central_";
  if (isAK8) name_pt_bin += "AK8_";
  name_pt_bin += year+"_ptbins";
  if(Trigger.Contains("default")) name_pt_bin += "_default";
  if(Trigger.Contains("fine_v1")) name_pt_bin += "_fine_v1";
  if(Trigger.Contains("fine_v2")) name_pt_bin += "_fine_v2";
  if(Trigger.Contains("fine_v3")) name_pt_bin += "_fine_v3";
  if(Trigger.Contains("fine_v4")) name_pt_bin += "_fine_v4";
  if(Trigger.Contains("fine_v5")) name_pt_bin += "_fine_v5";
  if(Trigger.Contains("fine_v6")) name_pt_bin += "_fine_v6";
  cout << year << " " << name_pt_bin << endl;
  PtBins_Central = pt_trigger_thr.at(name_pt_bin).size();
  for (auto &pt: pt_trigger_thr.at(name_pt_bin)) Pt_bins_Central.push_back(pt);

  name_pt_bin = triggerName+"_forward_";
  if (isAK8) name_pt_bin += "AK8_";
  name_pt_bin += year+"_ptbins";
  if(Trigger.Contains("default")) name_pt_bin += "_default";
  if(Trigger.Contains("fine_v1")) name_pt_bin += "_fine_v1";
  if(Trigger.Contains("fine_v2")) name_pt_bin += "_fine_v2";
  if(Trigger.Contains("fine_v3")) name_pt_bin += "_fine_v3";
  if(Trigger.Contains("fine_v4")) name_pt_bin += "_fine_v4";
  if(Trigger.Contains("fine_v5")) name_pt_bin += "_fine_v5";
  if(Trigger.Contains("fine_v6")) name_pt_bin += "_fine_v6";
  cout << year << " " << name_pt_bin << endl;
  PtBins_HF = pt_trigger_thr.at(name_pt_bin).size();
  for (auto &pt: pt_trigger_thr.at(name_pt_bin)) Pt_bins_HF.push_back(pt);
  Pt_bins_Central.push_back(7000);
  Pt_bins_HF.push_back(7000);
  if(year.find("UL16") != std::string::npos && isAK8 ) Pt_bins_HF = Pt_bins_Central;

  if(Trigger.Contains("fine")) sizeHist = 3100;

  usedPtTrigger_central = Pt_bins_Central; // global variable used in function.C (correctForRef)
  usedPtTrigger_forward = Pt_bins_HF; // global variable used in function.C (correctForRef)

  if (debug) {
    std::cout << "Trigger " << Trigger << " " << Pt_bins_Central.size() << " " << Pt_bins_HF.size() << '\n';
    std::cout << "Pt_bins_Central\t"; for (size_t i = 0; i < Pt_bins_Central.size(); i++) std::cout << Pt_bins_Central[i] << '\t'; std::cout << '\n';
    std::cout << "Pt_bins_HF\t";  for (size_t i = 0; i < Pt_bins_HF.size(); i++) std::cout << Pt_bins_HF[i] << '\t'; std::cout << '\n';
  }

  //std::vector<double> eta_bins_edge_SM(eta_bins, eta_bins + sizeof(eta_bins)/sizeof(double));
  // std::vector<double> eta_bins_edge_FE(eta_bins+1, eta_bins + sizeof(eta_bins)/sizeof(double));

  std::vector<double> eta_bins_edge_SM(&eta_bins[0], &eta_bins[0]+n_eta_bins);
  std::vector<double> eta_bins_edge_FE(&eta_bins[1], &eta_bins[0]+n_eta_bins);

  // std::vector<double> eta_bins_edge_SM(eta_bins2, eta_bins2 + sizeof(eta_bins2)/sizeof(double));
  // std::vector<double> eta_bins_edge_FE(eta_bins2, eta_bins2 + sizeof(eta_bins2)/sizeof(double));

  if (debug) {
    std::cout << "Eta bins " << Trigger << " " << eta_bins_edge_SM.size() << " " << eta_bins_edge_FE.size() << '\n';
    for (size_t i = 0; i < eta_bins_edge_SM.size(); i++) std::cout << eta_bins_edge_SM[i] << '\t'; std::cout << '\n';
    for (size_t i = 0; i < eta_bins_edge_FE.size(); i++) std::cout << eta_bins_edge_FE[i] << '\t'; std::cout << '\n';
  }

  // EtaBins_FE   = 3;
  // EtaBins_SM       = n_eta_bins - EtaBins_FE - 1;
  // etaShift_FE_control  = eta_bins_edge_FE.size() -1 ;

  alpha = {0.05,0.10,0.15,0.20,0.25,0.30};
  if(Trigger.Contains("finealpha")) alpha = {0.05,0.075,0.10,0.125,0.15,0.175,0.20,0.225,0.25,0.275,0.30};
  if(Trigger.Contains("highalpha")) alpha = {0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.};
  AlphaBins = alpha.size();
  cout << "Consider " << AlphaBins << " alpha bins;";
  for(double a:alpha){printf(" %1.3f",a);} cout << endl;

  TFile *f, *f_data;
  f = new TFile( filename, "READ");
  f_data = new TFile( filename_data, "READ");

  if (debug) {
    std::cout << "filename MC " << filename << '\n';
    std::cout << "filename DATA " << filename_data << '\n';
  }

  // ------------------------------
  //      loading histograms
  // ------------------------------

  time(&start); // TIME

  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_SM, gen_asymmetries_SM, asymmetries_data_SM, MC_Truth_asymmetries_SM;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_pt_SM, gen_asymmetries_pt_SM, asymmetries_pt_data_SM;

  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_FE, gen_asymmetries_FE, asymmetries_data_FE, MC_Truth_asymmetries_FE;
  std::vector< std::vector< std::vector< TH1F* > > > asymmetries_pt_FE, gen_asymmetries_pt_FE, asymmetries_pt_data_FE;

  std::vector< std::vector< std::vector< TH1F* > > > gen_asymmetries_data_SM, gen_asymmetries_data_FE, gen_asymmetries_pt_data_SM, gen_asymmetries_pt_data_FE;

  std::vector< std::vector< std::vector< TH1F* > > > dummy_hists;

  // real_data = true, data = false
  histLoadAsym( *f_data,  real_data,  "asymm_SM",           asymmetries_data_SM,  gen_asymmetries_data_SM,  EtaBins_SM,           PtBins_Central, AlphaBins, etaShift_SM);
  histLoadAsym( *f,       data_,      "asymm_SM",           asymmetries_SM,       gen_asymmetries_SM,       EtaBins_SM,           PtBins_Central, AlphaBins, etaShift_SM);
  histLoadAsym( *f_data,  real_data,  "asymm_SM_control",   asymmetries_data_SM,  gen_asymmetries_data_SM,  EtaBins_SM_control,   PtBins_HF,      AlphaBins, etaShift_SM_control);
  histLoadAsym( *f,       data_,      "asymm_SM_control",   asymmetries_SM,       gen_asymmetries_SM,       EtaBins_SM_control,   PtBins_HF,      AlphaBins, etaShift_SM_control);
  histLoadAsym( *f_data,  real_data,  "asymm_FE_reference", asymmetries_data_FE,  gen_asymmetries_data_FE,  EtaBins_FE_reference, PtBins_Central, AlphaBins, etaShift_FE_reference);
  histLoadAsym( *f,       data_,      "asymm_FE_reference", asymmetries_FE,       gen_asymmetries_FE,       EtaBins_FE_reference, PtBins_Central, AlphaBins, etaShift_FE_reference);
  histLoadAsym( *f_data,  real_data,  "asymm_FE_control",   asymmetries_data_FE,  gen_asymmetries_data_FE,  EtaBins_FE_control,   PtBins_Central, AlphaBins, etaShift_FE_control);
  histLoadAsym( *f,       data_,      "asymm_FE_control",   asymmetries_FE,       gen_asymmetries_FE,       EtaBins_FE_control,   PtBins_Central, AlphaBins, etaShift_FE_control);
  histLoadAsym( *f_data,  real_data,  "asymm_FE",           asymmetries_data_FE,  gen_asymmetries_data_FE,  EtaBins_FE,           PtBins_HF,      AlphaBins, etaShift_FE);
  histLoadAsym( *f,       data_,      "asymm_FE",           asymmetries_FE,       gen_asymmetries_FE,       EtaBins_FE,           PtBins_HF,      AlphaBins, etaShift_FE);

  if (debug) {
    std::cout << "asymmetries_data_SM " << asymmetries_data_SM.size() << "=" << EtaBins_SM << "+" << EtaBins_SM_control << '\n';
    std::cout << "asymmetries_data_FE " << asymmetries_data_FE.size() << "=" << EtaBins_FE_reference << "+" << EtaBins_FE_control << "+" << EtaBins_FE << '\n';
  }

  histLoadAsym( *f_data,  real_data,  "asymmpt_SM",           asymmetries_pt_data_SM,  gen_asymmetries_pt_data_SM,  EtaBins_SM,           PtBins_Central, AlphaBins, etaShift_SM);
  histLoadAsym( *f,       data_,      "asymmpt_SM",           asymmetries_pt_SM,       gen_asymmetries_pt_SM,       EtaBins_SM,           PtBins_Central, AlphaBins, etaShift_SM);
  histLoadAsym( *f_data,  real_data,  "asymmpt_SM_control",   asymmetries_pt_data_SM,  gen_asymmetries_pt_data_SM,  EtaBins_SM_control,   PtBins_HF,      AlphaBins, etaShift_SM_control);
  histLoadAsym( *f,       data_,      "asymmpt_SM_control",   asymmetries_pt_SM,       gen_asymmetries_pt_SM,       EtaBins_SM_control,   PtBins_HF,      AlphaBins, etaShift_SM_control);
  histLoadAsym( *f_data,  real_data,  "asymmpt_FE_reference", asymmetries_pt_data_FE,  gen_asymmetries_pt_data_FE,  EtaBins_FE_reference, PtBins_Central, AlphaBins, etaShift_FE_reference);
  histLoadAsym( *f,       data_,      "asymmpt_FE_reference", asymmetries_pt_FE,       gen_asymmetries_pt_FE,       EtaBins_FE_reference, PtBins_Central, AlphaBins, etaShift_FE_reference);
  histLoadAsym( *f_data,  real_data,  "asymmpt_FE_control",   asymmetries_pt_data_FE,  gen_asymmetries_pt_data_FE,  EtaBins_FE_control,   PtBins_Central, AlphaBins, etaShift_FE_control);
  histLoadAsym( *f,       data_,      "asymmpt_FE_control",   asymmetries_pt_FE,       gen_asymmetries_pt_FE,       EtaBins_FE_control,   PtBins_Central, AlphaBins, etaShift_FE_control);
  histLoadAsym( *f_data,  real_data,  "asymmpt_FE",           asymmetries_pt_data_FE,  gen_asymmetries_pt_data_FE,  EtaBins_FE,           PtBins_HF,      AlphaBins, etaShift_FE);
  histLoadAsym( *f,       data_,      "asymmpt_FE",           asymmetries_pt_FE,       gen_asymmetries_pt_FE,       EtaBins_FE,           PtBins_HF,      AlphaBins, etaShift_FE);

  histLoadAsym( *f,       data_,       "mctruth_SM",           MC_Truth_asymmetries_SM,  dummy_hists,       EtaBins_SM,           PtBins_Central, AlphaBins, etaShift_SM);
  histLoadAsym( *f,       data_,       "mctruth_SM_control",   MC_Truth_asymmetries_SM,  dummy_hists,       EtaBins_SM_control,   PtBins_HF,      AlphaBins, etaShift_SM_control);
  histLoadAsym( *f,       data_,       "mctruth_FE_reference", MC_Truth_asymmetries_FE,  dummy_hists,       EtaBins_FE_reference, PtBins_Central, AlphaBins, etaShift_FE_reference);
  histLoadAsym( *f,       data_,       "mctruth_FE_control",   MC_Truth_asymmetries_FE,  dummy_hists,       EtaBins_FE_control,   PtBins_Central, AlphaBins, etaShift_FE_control);
  histLoadAsym( *f,       data_,       "mctruth_FE",           MC_Truth_asymmetries_FE,  dummy_hists,       EtaBins_FE,           PtBins_HF,      AlphaBins, etaShift_FE);

  if (debug) {
    std::cout << "asymmpt_SM " << asymmetries_pt_data_SM.size() << '\n';
    std::cout << "asymmpt_FE " << asymmetries_pt_data_FE.size() << '\n';
  }

  // std::vector < TH2F* > Map_mean_data;
  //
  // Fill_Map3D(asymmetries_data_SM, Map_mean_data, eta_bins_edge_SM, Pt_bins_Central);
  //
  // TFile maps("pdfy/maps/maps.root","RECREATE");
  // for (unsigned int r = 0; r < alpha.size(); r++) {
  //   TString canvName = "Map_mean_"; canvName += r;
  //   TString nameXaxis = "#eta";
  //   TString nameYaxis = "p_{T} (GeV)";
  //   TCanvas* canv = tdrCanvas(canvName, 0, 6, 0, 1100, nameXaxis, nameYaxis);
  //   canv->SetTickx(0);
  //   canv->SetTicky(0);
  //   gStyle->SetPalette(kLightTemperature);
  //   Map_mean_data.at(r)->Draw("colz");
  //   canv->Update();
  //   Map_mean_data.at(r)->Write();
  //   canv -> Print(outdir+"pdfy/maps/"+canvName+".pdf","pdf");
  //   canv -> Print(outdir+"pdfy/maps/"+canvName+".png","png");
  // }
  // maps.Close();

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by histLoadAsym is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I calculate pt_mean for each alpha and pt bin and eta bin.          //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME

  std::vector< std::vector< std::vector< double > > > width_pt_SM, gen_width_pt_SM, width_pt_data_SM, gen_width_pt_data_SM;
  std::vector< std::vector< std::vector< double > > > width_pt_FE, gen_width_pt_FE, width_pt_data_FE, gen_width_pt_data_FE;

  histMeanPt( asymmetries_pt_SM,          width_pt_SM);
  histMeanPt( gen_asymmetries_pt_SM,      gen_width_pt_SM);
  histMeanPt( asymmetries_pt_data_SM,     width_pt_data_SM);
  histMeanPt( gen_asymmetries_pt_data_SM, gen_width_pt_data_SM);

  histMeanPt( asymmetries_pt_FE,          width_pt_FE );
  histMeanPt( gen_asymmetries_pt_FE,      gen_width_pt_FE );
  histMeanPt( asymmetries_pt_data_FE,     width_pt_data_FE );
  histMeanPt( gen_asymmetries_pt_data_FE, gen_width_pt_data_FE );

  if (debug) {
    std::cout << "histMeanPt_SM " << width_pt_SM.size() << '\n';
    std::cout << "histMeanPt_FE " << width_pt_FE.size() << '\n';
  }

  // gen_all_data is for running a cross check with smeared MC.

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by histMeanPt is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I calculate width of asymmetry distributions only for               //
  //    alpha bins above 10 GeV thresholds (too soft contriubtions)         //
  //    e.g. for bin p_T_ave (55-75) alpha 0.1 corresponds to 57 GeV jets   //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< std::vector< double > > > asymmetries_width_SM, gen_asymmetries_width_SM, asymmetries_width_data_SM, gen_asymmetries_width_data_SM, lower_x_SM, upper_x_SM, gen_lower_x_SM, gen_upper_x_SM, lower_x_data_SM, upper_x_data_SM, gen_lower_x_data_SM, gen_upper_x_data_SM;
  std::vector< std::vector< std::vector< double > > > asymmetries_width_SM_error, gen_asymmetries_width_SM_error, asymmetries_width_data_SM_error, gen_asymmetries_width_data_SM_error;
  std::vector< std::vector< std::vector< double > > > asymmetries_width_FE, gen_asymmetries_width_FE, asymmetries_width_data_FE, gen_asymmetries_width_data_FE, lower_x_FE, upper_x_FE, gen_lower_x_FE, gen_upper_x_FE, lower_x_data_FE, upper_x_data_FE, gen_lower_x_data_FE, gen_upper_x_data_FE;
  std::vector< std::vector< std::vector< double > > > asymmetries_width_FE_error, gen_asymmetries_width_FE_error, asymmetries_width_data_FE_error, gen_asymmetries_width_data_FE_error;

  histWidthAsym( asymmetries_SM , asymmetries_width_SM, asymmetries_width_SM_error, false, gaustails, 0, lower_x_SM, upper_x_SM, eta_bins, alpha);
  histWidthAsym( gen_asymmetries_SM , gen_asymmetries_width_SM, gen_asymmetries_width_SM_error, false, gaustails, 0, gen_lower_x_SM, gen_upper_x_SM, eta_bins, alpha);
  histWidthAsym( asymmetries_data_SM , asymmetries_width_data_SM, asymmetries_width_data_SM_error, false, gaustails, 0, lower_x_data_SM, upper_x_data_SM, eta_bins, alpha);
  histWidthAsym( gen_asymmetries_data_SM , gen_asymmetries_width_data_SM, gen_asymmetries_width_data_SM_error, false, gaustails, 0, gen_lower_x_data_SM, gen_upper_x_data_SM, eta_bins, alpha);

  histWidthAsym( asymmetries_FE , asymmetries_width_FE, asymmetries_width_FE_error, false, gaustails, 1, lower_x_FE, upper_x_FE, eta_bins, alpha);
  histWidthAsym( gen_asymmetries_FE , gen_asymmetries_width_FE, gen_asymmetries_width_FE_error, false, gaustails, 1, gen_lower_x_FE, gen_upper_x_FE, eta_bins, alpha);
  histWidthAsym( asymmetries_data_FE , asymmetries_width_data_FE, asymmetries_width_data_FE_error, false, gaustails, 1, lower_x_data_FE, upper_x_data_FE, eta_bins, alpha);
  histWidthAsym( gen_asymmetries_data_FE , gen_asymmetries_width_data_FE, gen_asymmetries_width_data_FE_error, false, gaustails, 1, gen_lower_x_data_FE, gen_upper_x_data_FE, eta_bins, alpha);

  if (debug) {
    std::cout << "asymmetries_width_SM " << asymmetries_width_SM.size() << '\n';
    std::cout << "asymmetries_width_Fe " << asymmetries_width_FE.size() << '\n';
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by histWidthAsym is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I calculate widths, this time also including                        //
  //    alpha bins below 10GeV threshold                                    //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< std::vector< double > > > soft_asymmetries_width_SM, soft_gen_asymmetries_width_SM, soft_asymmetries_width_data_SM, soft_gen_asymmetries_width_data_SM;
  std::vector< std::vector< std::vector< double > > > soft_asymmetries_width_SM_error, soft_gen_asymmetries_width_SM_error, soft_asymmetries_width_data_SM_error, soft_gen_asymmetries_width_data_SM_error;

  std::vector< std::vector< std::vector< double > > > soft_asymmetries_width_FE, gen_soft_asymmetries_width_FE, soft_asymmetries_width_data_FE, gen_soft_asymmetries_width_data_FE;
  std::vector< std::vector< std::vector< double > > > soft_asymmetries_width_FE_error, gen_soft_asymmetries_width_FE_error, soft_asymmetries_width_data_FE_error, gen_soft_asymmetries_width_data_FE_error;
  std::vector< std::vector< std::vector< double > > > dummy_vec_x, dummy_vec_y;

  histWidthAsym( asymmetries_SM , soft_asymmetries_width_SM, soft_asymmetries_width_SM_error, true, gaustails, 0, dummy_vec_x, dummy_vec_y, eta_bins, alpha);
  histWidthAsym( gen_asymmetries_SM , soft_gen_asymmetries_width_SM, soft_gen_asymmetries_width_SM_error, true, gaustails, 0, dummy_vec_x, dummy_vec_y, eta_bins, alpha);
  histWidthAsym( asymmetries_data_SM , soft_asymmetries_width_data_SM, soft_asymmetries_width_data_SM_error, true, gaustails, 0, dummy_vec_x, dummy_vec_y, eta_bins, alpha);
  histWidthAsym( gen_asymmetries_data_SM , soft_gen_asymmetries_width_data_SM, soft_gen_asymmetries_width_data_SM_error, true, gaustails, 0, dummy_vec_x, dummy_vec_y, eta_bins, alpha);

  histWidthAsym( asymmetries_FE , soft_asymmetries_width_FE, soft_asymmetries_width_FE_error, true, gaustails, 1, dummy_vec_x, dummy_vec_y, eta_bins, alpha);
  histWidthAsym( gen_asymmetries_FE , gen_soft_asymmetries_width_FE, gen_soft_asymmetries_width_FE_error, true, gaustails, 1, dummy_vec_x, dummy_vec_y, eta_bins, alpha);
  histWidthAsym( asymmetries_data_FE , soft_asymmetries_width_data_FE, soft_asymmetries_width_data_FE_error, true, gaustails, 1, dummy_vec_x, dummy_vec_y, eta_bins, alpha);
  histWidthAsym( gen_asymmetries_data_FE , gen_soft_asymmetries_width_data_FE, gen_soft_asymmetries_width_data_FE_error, true, gaustails, 1, dummy_vec_x, dummy_vec_y, eta_bins, alpha);

  if (debug) {
    std::cout << "asymmetries_SM " << asymmetries_SM.size() << '\n';
    std::cout << "asymmetries_FE " << asymmetries_FE.size() << '\n';
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by histWidthAsym add is : " << fixed << setprecision(2) << time_taken << " sec " << endl;


  ////////////////////////////////////////////////////////////////////////////
  //    Calculate mcTruth resolution for cross check with dijet calculation //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< std::vector< double > > > mcTruth_res_SM, mcTruth_res_FE;
  std::vector< std::vector< std::vector< double > > > mcTruth_res_SM_error, mcTruth_res_FE_error;

  histWidthMCTruth( MC_Truth_asymmetries_SM, mcTruth_res_SM, mcTruth_res_SM_error);
  histWidthMCTruth( MC_Truth_asymmetries_FE, mcTruth_res_FE, mcTruth_res_FE_error);

  if (debug) {
    std::cout << "MC_Truth_asymmetries_SM " << MC_Truth_asymmetries_SM.size() << '\n';
    std::cout << "MC_Truth_asymmetries_FE " << MC_Truth_asymmetries_FE.size() << '\n';
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by histWidthMCTruth is : " << fixed << setprecision(2) << time_taken << " sec " << endl;


  ////////////////////////////////////////////////////////////////////////////
  //     I fill width(alpha_max) histograms                                 //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< TH1F* > > widths_hist_SM, gen_widths_hist_SM, widths_hist_data_SM, gen_widths_hist_data_SM;
  std::vector< std::vector< TH1F* > > widths_hist_FE, gen_widths_hist_FE, widths_hist_data_FE, gen_widths_hist_data_FE;

  fill_widths_hists( "widths", widths_hist_SM , asymmetries_width_SM, asymmetries_width_SM_error );
  fill_widths_hists( "widths_gen", gen_widths_hist_SM , gen_asymmetries_width_SM, gen_asymmetries_width_SM_error );
  fill_widths_hists( "widths_data", widths_hist_data_SM , asymmetries_width_data_SM, asymmetries_width_data_SM_error );
  fill_widths_hists( "widths_gen_data", gen_widths_hist_data_SM , gen_asymmetries_width_data_SM, gen_asymmetries_width_data_SM_error );

  fill_widths_hists( "widths_fe", widths_hist_FE , asymmetries_width_FE, asymmetries_width_FE_error );
  fill_widths_hists( "widths_gen_fe", gen_widths_hist_FE , gen_asymmetries_width_FE, gen_asymmetries_width_FE_error );
  fill_widths_hists( "widths_data_fe", widths_hist_data_FE , asymmetries_width_data_FE, asymmetries_width_data_FE_error );
  fill_widths_hists( "widths_gen_data_fe", gen_widths_hist_data_FE , gen_asymmetries_width_data_FE, gen_asymmetries_width_data_FE_error );

  if (debug) {
    std::cout << "widths_hist_SM " << widths_hist_SM.size() << '\n';
    std::cout << "widths_hist_FE " << widths_hist_FE.size() << '\n';
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by fill_widths_hists is : " << fixed << setprecision(2) << time_taken << " sec " << endl;
  
  ////////////////////////////////////////////////////////////////////////////
  //    I do same for alpha unconstrained widths                            //
  //    one needs these plots to prove which points should be rejected!     //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< TH1F* > > soft_widths_hist_SM, soft_gen_widths_hist_SM, soft_widths_hist_data_SM, soft_gen_widths_hist_data_SM;
  std::vector< std::vector< TH1F* > > soft_widths_hist_FE, soft_gen_widths_hist_FE, soft_widths_hist_data_FE, soft_gen_widths_hist_data_FE;

  fill_widths_hists( "all_widths", soft_widths_hist_SM, soft_asymmetries_width_SM, soft_asymmetries_width_SM_error );
  fill_widths_hists( "all_widths_gen", soft_gen_widths_hist_SM, soft_gen_asymmetries_width_SM, soft_gen_asymmetries_width_SM_error );
  fill_widths_hists( "all_widths_data", soft_widths_hist_data_SM, soft_asymmetries_width_data_SM, soft_asymmetries_width_data_SM_error );
  fill_widths_hists( "all_widths_gen_data", soft_gen_widths_hist_data_SM, soft_gen_asymmetries_width_data_SM, soft_gen_asymmetries_width_data_SM_error );

  fill_widths_hists( "all_widths_fe", soft_widths_hist_FE, soft_asymmetries_width_FE, soft_asymmetries_width_FE_error );
  fill_widths_hists( "all_widths_gen_fe", soft_gen_widths_hist_FE, gen_soft_asymmetries_width_FE, gen_soft_asymmetries_width_FE_error );
  fill_widths_hists( "all_widths_data_fe", soft_widths_hist_data_FE, soft_asymmetries_width_data_FE, soft_asymmetries_width_data_FE_error );
  fill_widths_hists( "all_widths_gen_data_fe", soft_gen_widths_hist_data_FE, gen_soft_asymmetries_width_data_FE, gen_soft_asymmetries_width_data_FE_error ); 
  
  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by fill_widths_hists 2 is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I fit line or const to width(alpha_max)                             //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< double > > extrapolated_widths_SM, extrapolated_gen_widths_SM, extrapolated_widths_data_SM, extrapolated_gen_widths_data_SM;
  std::vector< std::vector< double > > extrapolated_widths_SM_error, extrapolated_gen_widths_SM_error, extrapolated_widths_data_SM_error, extrapolated_gen_widths_data_SM_error;

  std::vector< std::vector< double > > extrapolated_widths_FE, extrapolated_gen_widths_FE, extrapolated_widths_data_FE, extrapolated_gen_widths_data_FE;
  std::vector< std::vector< double > > extrapolated_widths_FE_error, extrapolated_gen_widths_FE_error, extrapolated_widths_data_FE_error, extrapolated_gen_widths_data_FE_error;

  histLinFit( widths_hist_SM , extrapolated_widths_SM, extrapolated_widths_SM_error, false );
  histLinFit( gen_widths_hist_SM , extrapolated_gen_widths_SM, extrapolated_gen_widths_SM_error, false );
  histLinFit( widths_hist_data_SM , extrapolated_widths_data_SM, extrapolated_widths_data_SM_error, false );
  histLinFit( gen_widths_hist_data_SM , extrapolated_gen_widths_data_SM, extrapolated_gen_widths_data_SM_error, false );

  histLinFit( widths_hist_FE , extrapolated_widths_FE, extrapolated_widths_FE_error, true );
  histLinFit( gen_widths_hist_FE , extrapolated_gen_widths_FE, extrapolated_gen_widths_FE_error, true );
  histLinFit( widths_hist_data_FE , extrapolated_widths_data_FE, extrapolated_widths_data_FE_error, true );
  histLinFit( gen_widths_hist_data_FE , extrapolated_gen_widths_data_FE, extrapolated_gen_widths_data_FE_error, true );

  if (debug) {
    std::cout << "extrapolated_widths_SM " << extrapolated_widths_SM.size() << '\n';
    std::cout << "extrapolated_widths_FE " << extrapolated_widths_FE.size() << '\n';
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by histLinFit is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    Correlated fit                                                      //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< double > > extrapolated_widths_correlated_SM, extrapolated_gen_widths_correlated_SM, extrapolated_widths_correlated_data_SM, extrapolated_gen_widths_correlated_data_SM;
  std::vector< std::vector< double > > extrapolated_widths_correlated_SM_error, extrapolated_gen_widths_correlated_SM_error, extrapolated_widths_correlated_data_SM_error, extrapolated_gen_widths_correlated_data_SM_error;

  std::vector< std::vector< double > > extrapolated_widths_correlated_FE, extrapolated_gen_widths_correlated_FE, extrapolated_widths_correlated_data_FE, extrapolated_gen_widths_correlated_data_FE;
  std::vector< std::vector< double > > extrapolated_widths_correlated_FE_error, extrapolated_gen_widths_correlated_FE_error, extrapolated_widths_correlated_data_FE_error, extrapolated_gen_widths_correlated_data_FE_error;

  std::vector< std::vector< TGraphErrors *> > MC_correlated_graphs_SM, data_correlated_graphs_SM, gen_correlated_graphs_SM, gen_data_correlated_graphs_SM;
  std::vector< std::vector< TGraphErrors *> > MC_correlated_graphs_FE, data_correlated_graphs_FE, gen_correlated_graphs_FE, gen_data_correlated_graphs_FE;
  TH1F* h_chi2_tot = new TH1F("h_chi2_tot","h_chi2_tot",100, -10,10);

  histLinCorFit(asymmetries_width_SM, asymmetries_width_SM_error, MC_correlated_graphs_SM, extrapolated_widths_correlated_SM, extrapolated_widths_correlated_SM_error, false, true, h_chi2_tot);
  histLinCorFit(asymmetries_width_data_SM, asymmetries_width_data_SM_error, data_correlated_graphs_SM, extrapolated_widths_correlated_data_SM, extrapolated_widths_correlated_data_SM_error, false, false, h_chi2_tot);
  histLinCorFit(gen_asymmetries_width_SM, gen_asymmetries_width_SM_error, gen_correlated_graphs_SM, extrapolated_gen_widths_correlated_SM, extrapolated_gen_widths_correlated_SM_error, false, false, h_chi2_tot);
  histLinCorFit(gen_asymmetries_width_data_SM, gen_asymmetries_width_data_SM_error, gen_data_correlated_graphs_SM, extrapolated_gen_widths_correlated_data_SM, extrapolated_gen_widths_correlated_data_SM_error, false, false, h_chi2_tot);

  histLinCorFit(asymmetries_width_FE, asymmetries_width_FE_error, MC_correlated_graphs_FE, extrapolated_widths_correlated_FE, extrapolated_widths_correlated_FE_error, true, true, h_chi2_tot);
  histLinCorFit(asymmetries_width_data_FE, asymmetries_width_data_FE_error, data_correlated_graphs_FE, extrapolated_widths_correlated_data_FE, extrapolated_widths_correlated_data_FE_error, true, false, h_chi2_tot);
  histLinCorFit(gen_asymmetries_width_FE, gen_asymmetries_width_FE_error, gen_correlated_graphs_FE, extrapolated_gen_widths_correlated_FE, extrapolated_gen_widths_correlated_FE_error, true, false, h_chi2_tot);
  histLinCorFit(gen_asymmetries_width_data_FE, gen_asymmetries_width_data_FE_error, gen_data_correlated_graphs_FE, extrapolated_gen_widths_correlated_data_FE, extrapolated_gen_widths_correlated_data_FE_error, true, false, h_chi2_tot);

  if (debug) {
    std::cout << "extrapolated_widths_correlated_SM " << extrapolated_widths_correlated_SM.size() << '\n';
    std::cout << "extrapolated_widths_correlated_FE " << extrapolated_widths_correlated_FE.size() << '\n';
  }

  if (debug){
    std::cout << "\n === extrapolated widths: (ex. data | ex. data err | ex. MC | ex. MC err)\n" << endl;
    for(unsigned int i=1; i<extrapolated_widths_correlated_data_FE[1].size(); i++ ){
      double content = extrapolated_widths_correlated_data_FE[1][i];
      if(content==0) continue;
      cout << setw(4) << i << setw(12) << extrapolated_widths_correlated_data_FE[1][i] << setw(12) << extrapolated_widths_correlated_data_FE_error[1][i] << setw(12) << extrapolated_widths_correlated_FE[1][i] << setw(12) << extrapolated_widths_correlated_FE_error[1][i] << endl;
    }
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by histLinCorFit is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I make histograms ratio of widths(alpha=0.15)                       //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< TH1F* > widths_015_ratios_SM, widths_015_ratios_FE;

  widths_015_ratios( "widths_015_SM_ratios", widths_015_ratios_SM, asymmetries_width_data_SM, asymmetries_width_data_SM_error, asymmetries_width_SM, asymmetries_width_SM_error, width_pt_SM );
  widths_015_ratios( "widths_015_FE_ratios", widths_015_ratios_FE, asymmetries_width_data_FE, asymmetries_width_data_FE_error, asymmetries_width_FE, asymmetries_width_FE_error, width_pt_FE );

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by widths_015_ratios is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I correct for PLI using b parameter from line fit to sigma_A(alpha) //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< double > > JER_uncorrelated_corrected_MC_SM, JER_uncorrelated_corrected_data_SM;
  std::vector< std::vector< double > > JER_uncorrelated_corrected_MC_SM_error, JER_uncorrelated_corrected_data_SM_error;

  std::vector< std::vector< double > > JER_uncorrelated_corrected_MC_FE_ref, JER_uncorrelated_corrected_data_FE_ref;
  std::vector< std::vector< double > > JER_uncorrelated_corrected_MC_FE_ref_error, JER_uncorrelated_corrected_data_FE_ref_error;

  correctJERwithPLI( JER_uncorrelated_corrected_MC_SM, JER_uncorrelated_corrected_MC_SM_error, extrapolated_widths_SM, extrapolated_widths_SM_error, extrapolated_gen_widths_SM, extrapolated_gen_widths_SM_error, shiftForPLI);
  if( real_data ){
    // correct data with PLI from MC
    correctJERwithPLI( JER_uncorrelated_corrected_data_SM, JER_uncorrelated_corrected_data_SM_error, extrapolated_widths_data_SM, extrapolated_widths_data_SM_error, extrapolated_gen_widths_SM, extrapolated_gen_widths_SM_error, shiftForPLI);
  } else {
    // correct 'data' with own PLI
    correctJERwithPLI( JER_uncorrelated_corrected_data_SM, JER_uncorrelated_corrected_data_SM_error, extrapolated_widths_data_SM, extrapolated_widths_data_SM_error, extrapolated_gen_widths_data_SM, extrapolated_gen_widths_data_SM_error, shiftForPLI);
  }

  correctJERwithPLI( JER_uncorrelated_corrected_MC_FE_ref, JER_uncorrelated_corrected_MC_FE_ref_error, extrapolated_widths_FE, extrapolated_widths_FE_error, extrapolated_gen_widths_FE, extrapolated_gen_widths_FE_error, shiftForPLI);
  if( real_data ){
    // correct data with PLI from MC
    correctJERwithPLI( JER_uncorrelated_corrected_data_FE_ref, JER_uncorrelated_corrected_data_FE_ref_error, extrapolated_widths_data_FE, extrapolated_widths_data_FE_error, extrapolated_gen_widths_FE, extrapolated_gen_widths_FE_error, shiftForPLI);
  } else {
    // correct 'data' with own PLI
    correctJERwithPLI( JER_uncorrelated_corrected_data_FE_ref, JER_uncorrelated_corrected_data_FE_ref_error, extrapolated_widths_data_FE, extrapolated_widths_data_FE_error, extrapolated_gen_widths_data_FE, extrapolated_gen_widths_data_FE_error, shiftForPLI);
  }

  if (debug) {
    std::cout << "JER_uncorrelated_corrected_MC_SM " << JER_uncorrelated_corrected_MC_SM.size() << '\n';
    std::cout << "JER_uncorrelated_corrected_MC_FE_ref " << JER_uncorrelated_corrected_MC_FE_ref.size() << '\n';
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by correctJERwithPLI is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    PLI corrected using b parameters                                    //
  ////////////////////////////////////////////////////////////////////////////
  //    Same correction but for correlated fit results                      //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< std::vector< double > > JER_correlated_corrected_MC_SM, JER_correlated_corrected_data_SM;
  std::vector< std::vector< double > > JER_correlated_corrected_MC_SM_error, JER_correlated_corrected_data_SM_error;

  std::vector< std::vector< double > > JER_correlated_corrected_MC_FE_ref, JER_correlated_corrected_data_FE_ref;
  std::vector< std::vector< double > > JER_correlated_corrected_MC_FE_ref_error, JER_correlated_corrected_data_FE_ref_error;

  correctJERwithPLI( JER_correlated_corrected_MC_SM, JER_correlated_corrected_MC_SM_error, extrapolated_widths_correlated_SM, extrapolated_widths_correlated_SM_error, extrapolated_gen_widths_correlated_SM, extrapolated_gen_widths_correlated_SM_error, shiftForPLI);
  if( real_data ){
    // correct data with PLI from MC
    correctJERwithPLI( JER_correlated_corrected_data_SM, JER_correlated_corrected_data_SM_error, extrapolated_widths_correlated_data_SM, extrapolated_widths_correlated_data_SM_error, extrapolated_gen_widths_correlated_SM, extrapolated_gen_widths_correlated_SM_error, shiftForPLI);
  } else {
    // correct 'data' with own PLI
    correctJERwithPLI( JER_correlated_corrected_data_SM, JER_correlated_corrected_data_SM_error, extrapolated_widths_correlated_data_SM, extrapolated_widths_correlated_data_SM_error, extrapolated_gen_widths_correlated_data_SM, extrapolated_gen_widths_correlated_data_SM_error, shiftForPLI);
  }

  correctJERwithPLI( JER_correlated_corrected_MC_FE_ref, JER_correlated_corrected_MC_FE_ref_error, extrapolated_widths_correlated_FE, extrapolated_widths_correlated_FE_error, extrapolated_gen_widths_correlated_FE, extrapolated_gen_widths_correlated_FE_error, shiftForPLI);
  if( real_data ){
    // correct data with PLI from MC
    correctJERwithPLI( JER_correlated_corrected_data_FE_ref, JER_correlated_corrected_data_FE_ref_error, extrapolated_widths_correlated_data_FE, extrapolated_widths_correlated_data_FE_error, extrapolated_gen_widths_correlated_FE, extrapolated_gen_widths_correlated_FE_error, shiftForPLI);
  } else {
    // correct 'data' with own PLI
    correctJERwithPLI( JER_correlated_corrected_data_FE_ref, JER_correlated_corrected_data_FE_ref_error, extrapolated_widths_correlated_data_FE, extrapolated_widths_correlated_data_FE_error, extrapolated_gen_widths_correlated_data_FE, extrapolated_gen_widths_correlated_data_FE_error, shiftForPLI);
  }

  if (debug) {
    std::cout << "JER_correlated_corrected_MC_SM " << JER_correlated_corrected_MC_SM.size() << '\n';
    std::cout << "JER_correlated_corrected_MC_FE_ref " << JER_correlated_corrected_MC_FE_ref.size() << '\n';
  }

  if (debug){
    std::cout << "\n === Corrected with PLI: (ex. data | ex. data err | ex. MC | ex. MC err)\n" << endl;
    for(unsigned int i=1; i<extrapolated_widths_correlated_data_FE[1].size(); i++ ){
      double content = extrapolated_widths_correlated_data_FE[1][i];
      if(content==0) continue;
      cout << setw(4) << i << setw(12) << JER_correlated_corrected_data_FE_ref[1][i] << setw(12) << JER_correlated_corrected_data_FE_ref_error[1][i] << setw(12) << JER_correlated_corrected_MC_FE_ref[1][i] << setw(12) << JER_correlated_corrected_MC_FE_ref_error[1][i] << endl;
    }
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by correctJERwithPLI 2 is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I do the same for widths at alpha = 0.15                            //
  ////////////////////////////////////////////////////////////////////////////
  if(debug) cout << "LINE " << __LINE__ << endl;

  time(&start); // TIME
  
  std::vector< std::vector< double > > JER015_MC_SM, JER015_data_SM;
  std::vector< std::vector< double > > JER015_MC_SM_error, JER015_data_SM_error;

  std::vector< std::vector< double > > JER015_MC_FE_ref, JER015_data_FE_ref;
  std::vector< std::vector< double > > JER015_MC_FE_ref_error, JER015_data_FE_ref_error;

  correctJERwithPLI015(JER015_MC_SM, JER015_MC_SM_error, asymmetries_width_SM, asymmetries_width_SM_error, gen_asymmetries_width_SM, gen_asymmetries_width_SM_error, shiftForPLI);
  if( real_data ){
    // correct data with PLI from MC
    correctJERwithPLI015(JER015_data_SM, JER015_data_SM_error, asymmetries_width_data_SM, asymmetries_width_data_SM_error, gen_asymmetries_width_SM, gen_asymmetries_width_SM_error, shiftForPLI);
  } else {
    // correct 'data' with own PLI
    correctJERwithPLI015(JER015_data_SM, JER015_data_SM_error, asymmetries_width_data_SM, asymmetries_width_data_SM_error, gen_asymmetries_width_data_SM, gen_asymmetries_width_data_SM_error, shiftForPLI);
  }
  correctJERwithPLI015(JER015_MC_FE_ref, JER015_MC_FE_ref_error, asymmetries_width_FE, asymmetries_width_FE_error, gen_asymmetries_width_FE, gen_asymmetries_width_FE_error, shiftForPLI);
  if( real_data ){
    // correct data with PLI from MC
    correctJERwithPLI015(JER015_data_FE_ref, JER015_data_FE_ref_error, asymmetries_width_data_FE, asymmetries_width_data_FE_error, gen_asymmetries_width_FE, gen_asymmetries_width_FE_error, shiftForPLI);
  } else {
    // correct 'data' with own PLI
    correctJERwithPLI015(JER015_data_FE_ref, JER015_data_FE_ref_error, asymmetries_width_data_FE, asymmetries_width_data_FE_error, gen_asymmetries_width_data_FE, gen_asymmetries_width_data_FE_error, shiftForPLI);
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by correctJERwithPLI015 is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I corrected alpha = 0.15 widhs for PLI correct way                  //
  ////////////////////////////////////////////////////////////////////////////
  //    I correct forward widths for Ref region.                            //
  ////////////////////////////////////////////////////////////////////////////
  if(debug) cout << "LINE " << __LINE__ << endl;

  time(&start); // TIME
  
  std::vector< std::vector< double > > JER_uncorrelated_corrected_MC_FE, JER_uncorrelated_corrected_data_FE;
  std::vector< std::vector< double > > JER_uncorrelated_corrected_MC_FE_error, JER_uncorrelated_corrected_data_FE_error;

  correctForRef( "mccorrected", JER_uncorrelated_corrected_MC_FE,   JER_uncorrelated_corrected_MC_FE_error,   JER_uncorrelated_corrected_MC_FE_ref,   JER_uncorrelated_corrected_MC_FE_ref_error,   JER_uncorrelated_corrected_MC_SM,   JER_uncorrelated_corrected_MC_SM_error,   width_pt_FE, ref_shift, outdir);
  correctForRef( "datacorrect", JER_uncorrelated_corrected_data_FE, JER_uncorrelated_corrected_data_FE_error, JER_uncorrelated_corrected_data_FE_ref, JER_uncorrelated_corrected_data_FE_ref_error, JER_uncorrelated_corrected_data_SM, JER_uncorrelated_corrected_data_SM_error, width_pt_FE, ref_shift, outdir);

  if (debug) std::cout << "JER_uncorrelated_corrected_MC_FE " << JER_uncorrelated_corrected_MC_FE.size() << '\n';
  if (debug) std::cout << "JER_uncorrelated_corrected_data_FE " << JER_uncorrelated_corrected_data_FE.size() << '\n';

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by correctForRef is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    forward widths corrected for Ref widths!                            //
  ////////////////////////////////////////////////////////////////////////////
  //    same correction for correlated fit results                          //
  ////////////////////////////////////////////////////////////////////////////
  if(debug) cout << "LINE " << __LINE__ << endl;

  time(&start); // TIME
  
  std::vector< std::vector< double > > JER_correlated_corrected_MC_FE, JER_correlated_corrected_data_FE;
  std::vector< std::vector< double > > JER_correlated_corrected_MC_FE_error, JER_correlated_corrected_data_FE_error;

  correctForRef( "mc_cor_corrected", JER_correlated_corrected_MC_FE,   JER_correlated_corrected_MC_FE_error,   JER_correlated_corrected_MC_FE_ref,   JER_correlated_corrected_MC_FE_ref_error,   JER_correlated_corrected_MC_SM,   JER_correlated_corrected_MC_SM_error,   width_pt_FE, ref_shift, outdir);
  correctForRef( "data_cor_correct", JER_correlated_corrected_data_FE, JER_correlated_corrected_data_FE_error, JER_correlated_corrected_data_FE_ref, JER_correlated_corrected_data_FE_ref_error, JER_correlated_corrected_data_SM, JER_correlated_corrected_data_SM_error, width_pt_FE, ref_shift, outdir);

  // if (debug) std::cout << "JER_correlated_corrected_MC_FE " << JER_correlated_corrected_MC_FE.size() << '\n';
  if(debug) std::cout << "JER_correlated_corrected_MC_FE " << JER_correlated_corrected_MC_FE.size() << '\n';
  if(debug) std::cout << "JER_correlated_corrected_data_FE " << JER_correlated_corrected_data_FE.size() << '\n';

  if (debug){
    std::cout << "\n === Corrected for Ref: (ex. data | ex. data err | ex. MC | ex. MC err)\n" << endl;
    for(unsigned int i=1; i<extrapolated_widths_correlated_data_FE[1].size(); i++ ){
      double content = extrapolated_widths_correlated_data_FE[1][i];
      if(content==0) continue;
      cout << setw(4) << i << setw(12) << JER_correlated_corrected_data_FE[1][i] << setw(12) << JER_correlated_corrected_data_FE_error[1][i] << setw(12) << JER_correlated_corrected_MC_FE[1][i] << setw(12) << JER_correlated_corrected_MC_FE_error[1][i] << endl;
    }
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by correctForRef 2 is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    ref region corrected for correlated fit                             //
  ////////////////////////////////////////////////////////////////////////////
  //    and again, Ref region for widths at alpha = 0.15                    //
  ////////////////////////////////////////////////////////////////////////////
  if(debug) cout << "LINE " << __LINE__ << endl;

  time(&start); // TIME
  
  std::vector< std::vector< double > > JER015_MC_FE, JER015_data_FE;
  std::vector< std::vector< double > > JER015_MC_FE_error, JER015_data_FE_error;

  correctForRef( "mccorrected015", JER015_MC_FE,   JER015_MC_FE_error,   JER015_MC_FE_ref,   JER015_MC_FE_ref_error,   JER015_MC_SM,   JER015_MC_SM_error,    width_pt_FE, ref_shift, outdir);
  correctForRef( "datacorrect015", JER015_data_FE, JER015_data_FE_error, JER015_data_FE_ref, JER015_data_FE_ref_error, JER015_data_SM, JER015_data_SM_error,  width_pt_FE, ref_shift, outdir);

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by correctForRef 3 is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    Ref reg corrected for widths at alpha = 0.15                        //
  ////////////////////////////////////////////////////////////////////////////
  //    I make make vectors with ratios of widths                           //
  ////////////////////////////////////////////////////////////////////////////
  if(debug) cout << "LINE " << __LINE__ << endl;

  time(&start); // TIME
  
  std::vector< std::vector< double > > scales_uncorrelated_SM, scales_uncorrelated_FE, scales_uncorrelated_FE_control;
  std::vector< std::vector< double > > scales_uncorrelated_SM_error, scales_uncorrelated_FE_error, scales_uncorrelated_FE_control_error;

  makeScales( scales_uncorrelated_SM,         scales_uncorrelated_SM_error,         JER_uncorrelated_corrected_data_SM,     JER_uncorrelated_corrected_data_SM_error,     JER_uncorrelated_corrected_MC_SM,     JER_uncorrelated_corrected_MC_SM_error );
  makeScales( scales_uncorrelated_FE,         scales_uncorrelated_FE_error,         JER_uncorrelated_corrected_data_FE,     JER_uncorrelated_corrected_data_FE_error,     JER_uncorrelated_corrected_MC_FE,     JER_uncorrelated_corrected_MC_FE_error );
  makeScales( scales_uncorrelated_FE_control, scales_uncorrelated_FE_control_error, JER_uncorrelated_corrected_data_FE_ref, JER_uncorrelated_corrected_data_FE_ref_error, JER_uncorrelated_corrected_MC_FE_ref, JER_uncorrelated_corrected_MC_FE_ref_error );

  // same thing for correlated fit results
  std::vector< std::vector< double > > scales_correlated_SM, scales_correlated_FE, scales_correlated_FE_control;
  std::vector< std::vector< double > > scales_correlated_SM_error, scales_correlated_FE_error, scales_correlated_FE_control_error;

  makeScales( scales_correlated_SM,         scales_correlated_SM_error,         JER_correlated_corrected_data_SM,          JER_correlated_corrected_data_SM_error,          JER_correlated_corrected_MC_SM,          JER_correlated_corrected_MC_SM_error );
  makeScales( scales_correlated_FE,         scales_correlated_FE_error,         JER_correlated_corrected_data_FE,       JER_correlated_corrected_data_FE_error,       JER_correlated_corrected_MC_FE,       JER_correlated_corrected_MC_FE_error );
  makeScales( scales_correlated_FE_control, scales_correlated_FE_control_error, JER_correlated_corrected_data_FE_ref, JER_correlated_corrected_data_FE_ref_error, JER_correlated_corrected_MC_FE_ref, JER_correlated_corrected_MC_FE_ref_error );  // uncorrected for reference region width, as a cross check. (i think)

  if (debug){
    std::cout << "\n === SFs: (SF FE | SF FE err | SF FE control | SF FE control err )\n" << endl;
    for(unsigned int i=1; i<JER_correlated_corrected_data_FE[1].size(); i++ ){
      double content = JER_correlated_corrected_data_FE[1][i];
      if(content==0) continue;
      cout << setw(4) << i << setw(12) << scales_correlated_FE[1][i] << setw(12) << scales_correlated_FE_error[1][i] << setw(12) << scales_correlated_FE_control[1][i] << setw(12) << scales_correlated_FE_control_error[1][i] << setw(12) << scales_correlated_FE_error[1][i] << endl;
    }
  }

  std::vector< std::vector< double > > scales015_SM, scales015_FE;
  std::vector< std::vector< double > > scales015_SM_error, scales015_FE_error;

  makeScales( scales015_SM, scales015_SM_error, JER015_data_SM, JER015_data_SM_error, JER015_MC_SM, JER015_MC_SM_error );
  makeScales( scales015_FE, scales015_FE_error, JER015_data_FE, JER015_data_FE_error, JER015_MC_FE, JER015_MC_FE_error );

  if (debug) {
    std::cout << "scales_uncorrelated_SM " << scales_uncorrelated_SM.size() << '\n';
    std::cout << "scales_uncorrelated_FE " << scales_uncorrelated_FE.size() << '\n';
    std::cout << "scales_uncorrelated_FE_control " << scales_uncorrelated_FE_control.size() << '\n';
    std::cout << "scales_correlated_FE " << scales_correlated_FE.size() << " " << scales_correlated_FE[0].size() << '\n';
    std::cout << "scales_correlated_FE_control " << scales_correlated_FE_control.size() << " " << scales_correlated_FE_control[0].size() << '\n';
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by makeScales is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I make plots with MCTruth: Res from dijet                           //
  ////////////////////////////////////////////////////////////////////////////

  time(&start); // TIME
  
  std::vector< TH1F* > JER_MC_Truth_SM, JER_MC_Truth_FE;

  fill_mctruth_hist( "MC_Truth"    , JER_MC_Truth_SM, mcTruth_res_SM, mcTruth_res_SM_error, width_pt_SM, hist_max_value);
  fill_mctruth_hist( "MC_Truth_Fwd", JER_MC_Truth_FE, mcTruth_res_FE, mcTruth_res_FE_error, width_pt_FE, hist_max_value);

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by fill_mctruth_hist is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  ////////////////////////////////////////////////////////////////////////////
  //    I make histograms with  JERs and scale factors                      //
  ////////////////////////////////////////////////////////////////////////////
     
  time(&start); // TIME
        
  std::vector< TH1F* > JER_uncorrelated_MC_hist_SM,         JER_uncorrelated_data_hist_SM,          JER_uncorrelated_scale_hist_SM,         JER015_uncorrelated_MC_hist_SM;
  std::vector< TH1F* > JER_uncorrelated_MC_hist_FE,         JER_uncorrelated_data_hist_FE,          JER_uncorrelated_scale_hist_FE,         JER015_uncorrelated_MC_hist_FE;
  std::vector< TH1F* > JER_uncorrelated_MC_hist_FE_control, JER_uncorrelated_data_hist_FE_control,  JER_FE_uncorrelated_scale_control_hist, JER015_uncorrelated_MC_hist_SM_control;

  fill_hist( "MC_JER_uncorrelated_SM",    JER_uncorrelated_MC_hist_SM,    JER_uncorrelated_corrected_MC_SM,   JER_uncorrelated_corrected_MC_SM_error,   width_pt_SM, hist_max_value);
  fill_hist( "data_JER_uncorrelated_SM",  JER_uncorrelated_data_hist_SM , JER_uncorrelated_corrected_data_SM, JER_uncorrelated_corrected_data_SM_error, width_pt_SM, hist_max_value);
  fill_hist( "SF_uncorrelated_SM",        JER_uncorrelated_scale_hist_SM, scales_uncorrelated_SM,             scales_uncorrelated_SM_error,             width_pt_SM, hist_max_value_SF);

  fill_hist( "MC_JER_uncorrelated_FE",    JER_uncorrelated_MC_hist_FE,    JER_uncorrelated_corrected_MC_FE,   JER_uncorrelated_corrected_MC_FE_error,   width_pt_FE, hist_max_value,1);
  fill_hist( "data_JER_uncorrelated_FE",  JER_uncorrelated_data_hist_FE,  JER_uncorrelated_corrected_data_FE, JER_uncorrelated_corrected_data_FE_error, width_pt_FE, hist_max_value,1);
  fill_hist( "SF_uncorrelated_FE",        JER_uncorrelated_scale_hist_FE, scales_uncorrelated_FE,             scales_uncorrelated_FE_error,             width_pt_FE, hist_max_value_SF,1);

  fill_hist( "MC_JER_uncorrelated_FE_control",    JER_uncorrelated_MC_hist_FE_control,    JER_uncorrelated_corrected_MC_FE_ref,   JER_uncorrelated_corrected_MC_FE_ref_error,   width_pt_FE, hist_max_value);
  fill_hist( "data_JER_uncorrelated_FE_control",  JER_uncorrelated_data_hist_FE_control,  JER_uncorrelated_corrected_data_FE_ref, JER_uncorrelated_corrected_data_FE_ref_error, width_pt_FE, hist_max_value);
  fill_hist( "SF_uncorrelated_FE_control",        JER_FE_uncorrelated_scale_control_hist, scales_uncorrelated_FE_control,         scales_uncorrelated_FE_control_error,           width_pt_FE, hist_max_value_SF);

  fill_hist( "MC_JER015_uncorrelated_SM", JER015_uncorrelated_MC_hist_SM, JER015_MC_SM,                       JER015_MC_SM_error,                       width_pt_SM,      hist_max_value);
  fill_hist( "MC_JER015_uncorrelated_FE", JER015_uncorrelated_MC_hist_FE, JER015_MC_FE,                       JER015_MC_FE_error,                       width_pt_FE,      hist_max_value,1);

  std::vector< TH1F* > JER015_scale_hist_SM, JER015_scale_hist_FE;

  fill_hist( "SF_SM015", JER015_scale_hist_SM,  scales015_SM, scales015_SM_error, width_pt_SM, hist_max_value_SF);
  fill_hist( "SF_FE015", JER015_scale_hist_FE,  scales015_FE, scales015_FE_error, width_pt_FE, hist_max_value_SF,1);

  std::vector< TH1F* > JER_correlated_MC_hist_SM,         JER_correlated_data_hist_SM,         JER_correlated_scale_hist_SM;
  std::vector< TH1F* > JER_correlated_MC_hist_FE,         JER_correlated_data_hist_FE,         JER_correlated_scale_hist_FE;
  std::vector< TH1F* > JER_correlated_MC_hist_FE_control, JER_correlated_data_hist_FE_control, JER_correlated_scale_hist_FE_control;

  fill_hist( "MC_JER_correlated_SM",   JER_correlated_MC_hist_SM,    JER_correlated_corrected_MC_SM,         JER_correlated_corrected_MC_SM_error,         width_pt_SM, hist_max_value);
  fill_hist( "data_JER_correlated_SM", JER_correlated_data_hist_SM,  JER_correlated_corrected_data_SM,       JER_correlated_corrected_data_SM_error,       width_pt_SM, hist_max_value);
  fill_hist( "SF_correlated_SM",       JER_correlated_scale_hist_SM, scales_correlated_SM, scales_correlated_SM_error, width_pt_SM, hist_max_value_SF);

  fill_hist( "MC_JER_correlated_FE",   JER_correlated_MC_hist_FE,    JER_correlated_corrected_MC_FE,   JER_correlated_corrected_MC_FE_error,   width_pt_FE, hist_max_value,1);
  fill_hist( "data_JER_correlated_FE", JER_correlated_data_hist_FE,  JER_correlated_corrected_data_FE, JER_correlated_corrected_data_FE_error, width_pt_FE, hist_max_value,1);
  fill_hist( "SF_correlated_FE",       JER_correlated_scale_hist_FE, scales_correlated_FE,   scales_correlated_FE_error,   width_pt_FE, hist_max_value_SF,1);

  fill_hist( "MC_JER_correlated_FE_control",   JER_correlated_MC_hist_FE_control,    JER_correlated_corrected_MC_FE_ref,   JER_correlated_corrected_MC_FE_ref_error,   width_pt_FE, hist_max_value);
  fill_hist( "data_JER_correlated_FE_control", JER_correlated_data_hist_FE_control,  JER_correlated_corrected_data_FE_ref, JER_correlated_corrected_data_FE_ref_error, width_pt_FE, hist_max_value);
  fill_hist( "SF_correlated_FE_control",       JER_correlated_scale_hist_FE_control, scales_correlated_FE_control, scales_correlated_FE_control_error, width_pt_FE, hist_max_value_SF);

  if (debug) {
    std::cout << "JER_correlated_MC_hist_SM " << JER_correlated_MC_hist_SM.size() << '\n';
    std::cout << "JER_correlated_data_hist_SM " << JER_correlated_data_hist_SM.size() << '\n';
    std::cout << "JER_correlated_scale_hist_SM " << JER_correlated_scale_hist_SM.size() << '\n';

    std::cout << "JER_correlated_MC_hist_FE " << JER_correlated_MC_hist_FE.size() << '\n';
    std::cout << "JER_correlated_data_hist_FE " << JER_correlated_data_hist_FE.size() << '\n';
    std::cout << "JER_correlated_scale_hist_FE " << JER_correlated_scale_hist_FE.size() << '\n';

    std::cout << "JER_correlated_MC_hist_FE_control " << JER_correlated_MC_hist_FE_control.size() << '\n';
    std::cout << "JER_correlated_data_hist_FE_control " << JER_correlated_data_hist_FE_control.size() << '\n';
    std::cout << "JER_correlated_scale_hist_FE_control " << JER_correlated_scale_hist_FE_control.size() << '\n';
  }

  if (debug){
    std::cout << "\n === Data JERs: (MC JER | Data JER | Data/MC | SF)\n" << endl;
    for(unsigned int i=1; i<JER_correlated_MC_hist_FE[1]->GetNbinsX(); i++ ){
      double content = JER_correlated_MC_hist_FE.at(1)->GetBinContent(i);
      if(content==0) continue;
      cout << setw(4) << i << setw(12) << JER_correlated_MC_hist_FE[1]->GetBinContent(i) << setw(12) << JER_correlated_data_hist_FE[1]->GetBinContent(i) << setw(12) << JER_correlated_data_hist_FE[1]->GetBinContent(i)/JER_correlated_MC_hist_FE[1]->GetBinContent(i) << setw(12) << JER_correlated_scale_hist_FE[1]->GetBinContent(i) << endl;
    }
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by fill_hist is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  //////////////////////////////////////////////////////////////////////////////////////////
  //    store SFs as 2D (pt,eta) histograms                                               //
  //////////////////////////////////////////////////////////////////////////////////////////
  if(false){ // Not studied each time. Keep out for due to time consumption
    time(&start);
    TH2Poly* JER_SF_uncorrelated_SM_2D = fill_2Dhist( "2D_SF_SM", scales_uncorrelated_SM, scales_uncorrelated_SM_error, Pt_bins_Central, Pt_bins_HF, eta_bins_edge_SM,eta_cut);
    TH2Poly* JER_SF_uncorrelated_FE_2D = fill_2Dhist( "2D_SF_FE", scales_uncorrelated_FE, scales_uncorrelated_FE_error, Pt_bins_Central, Pt_bins_HF, eta_bins_edge_FE,eta_cut);
    TH2Poly* JER_SF_correlated_SM_2D   = fill_2Dhist( "2D_SF_SM", scales_correlated_SM,   scales_correlated_SM_error,   Pt_bins_Central, Pt_bins_HF, eta_bins_edge_SM,eta_cut);
    TH2Poly* JER_SF_correlated_FE_2D   = fill_2Dhist( "2D_SF_FE", scales_correlated_FE,   scales_correlated_FE_error,   Pt_bins_Central, Pt_bins_HF, eta_bins_edge_FE,eta_cut);

    TFile JERSF2Droot(outdir+"output/DijetJERSF2D.root","RECREATE");
    JER_SF_uncorrelated_FE_2D->Write();
    JER_SF_uncorrelated_SM_2D->Write();
    JER_SF_correlated_FE_2D->Write();
    JER_SF_correlated_SM_2D->Write();
    JERSF2Droot.Close();

    TCanvas* canv_2D_SF = new TCanvas();
    gStyle->SetPaintTextFormat("5.2f");
    canv_2D_SF->SetTickx(0);
    canv_2D_SF->SetTicky(0);
    JER_SF_uncorrelated_FE_2D->Draw("colz TEXT0");
    canv_2D_SF->SaveAs(outdir+"output/DijetJERSF2D_FE_uncorrelated.pdf");
    JER_SF_uncorrelated_SM_2D->Draw("colz TEXT0");
    canv_2D_SF->SaveAs(outdir+"output/DijetJERSF2D_SM_uncorrelated.pdf");
    JER_SF_correlated_FE_2D->Draw("colz TEXT0");
    canv_2D_SF->SaveAs(outdir+"output/DijetJERSF2D_FE_correlated.pdf");
    JER_SF_correlated_SM_2D->Draw("colz TEXT0");
    canv_2D_SF->SaveAs(outdir+"output/DijetJERSF2D_SM_correlated.pdf");
    delete JER_SF_uncorrelated_SM_2D;
    delete JER_SF_uncorrelated_FE_2D;
    delete JER_SF_correlated_SM_2D;
    delete JER_SF_correlated_FE_2D;

    time(&end); // TIME
    time_taken = double(end - start);
    if(debug) cout << "Time taken by 2D hists is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  }
  //////////////////////////////////////////////////////////////////////////////////////////
  //    resolution cross check with mcTruth                                               //
  //////////////////////////////////////////////////////////////////////////////////////////
  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_13TeV = lumi; //"2.1 fb^{-1}";
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  time(&start); // TIME

  std::cout << "START PLOTTING" << '\n';
  if(false){
    TFile fMCTruth("output/MCTruth.root","RECREATE");
    PLOT_MCT(JER_MC_Truth_SM,JER_uncorrelated_MC_hist_SM,JER_correlated_MC_hist_SM,JER015_uncorrelated_MC_hist_SM,outdir+"pdfy/MCTruth/",eta_bins_edge_SM, false);
    PLOT_MCT(JER_MC_Truth_FE,JER_uncorrelated_MC_hist_FE,JER_correlated_MC_hist_FE,JER015_uncorrelated_MC_hist_FE,outdir+"pdfy/MCTruth/",eta_bins_edge_FE, true);
    fMCTruth.Close();
  }

  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by PLOT_MCT is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  bool plot_all = false;
  if (plot_all) {
    ////////////////////////////////////////////////////////////////////////////
    //  Plots Asymmetries                                                     //
    ////////////////////////////////////////////////////////////////////////////
    


    if(true){ // Ignore to save time consumption
    
	    time(&start); // TIME
    
    	PLOT_ASY(asymmetries_data_SM,  asymmetries_SM, gen_asymmetries_SM,  asymmetries_width_data_SM, asymmetries_width_SM, gen_asymmetries_width_SM, asymmetries_width_data_SM_error, asymmetries_width_SM_error, gen_asymmetries_width_SM_error, outdir+"output/asymmetries/", eta_bins_edge_SM, Pt_bins_Central, Pt_bins_HF, alpha, lower_x_data_SM, upper_x_data_SM, lower_x_SM, upper_x_SM, gen_lower_x_SM, gen_upper_x_SM);
    	PLOT_ASY(asymmetries_data_FE,  asymmetries_FE, gen_asymmetries_FE,  asymmetries_width_data_FE, asymmetries_width_FE, gen_asymmetries_width_FE, asymmetries_width_data_FE_error, asymmetries_width_FE_error, gen_asymmetries_width_FE_error, outdir+"output/asymmetries/", eta_bins_edge_SM, Pt_bins_Central, Pt_bins_HF, alpha, lower_x_data_FE, upper_x_data_FE, lower_x_FE, upper_x_FE, gen_lower_x_FE, gen_upper_x_FE);
      
    	time(&end); // TIME
    	time_taken = double(end - start);
    	if(debug) cout << "Time taken by PLOT_ASY is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

      	time(&start); // TIME
      	TFile asyroot(outdir+"output/asymmetries.root","RECREATE");
      	for( unsigned int m = 0; m < asymmetries_SM.size(); m++){
          for( unsigned int p = 0; p < asymmetries_SM.at(m).size(); p++){
            for( unsigned int r = 0; r < asymmetries_SM.at(m).at(p).size(); r++){
                asymmetries_SM.at(m).at(p).at(r) -> Write();
                gen_asymmetries_SM.at(m).at(p).at(r) -> Write();
                asymmetries_data_SM.at(m).at(p).at(r) -> Write();
              }
            }
        } 
        for( unsigned int m = 0; m < asymmetries_FE.size(); m++){
          for( unsigned int p = 0; p < asymmetries_FE.at(m).size(); p++){
            for( unsigned int r = 0; r < asymmetries_FE.at(m).at(p).size(); r++){
              asymmetries_FE.at(m).at(p).at(r) -> Write();
              gen_asymmetries_FE.at(m).at(p).at(r) -> Write();
              asymmetries_data_FE.at(m).at(p).at(r) -> Write();
            }
          }
        } 
        asyroot.Close();

        time(&end); // TIME
        time_taken = double(end - start);
        if(debug) cout << "Time taken by STORE_ASY is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

    }

    ////////////////////////////////////////////////////////////////////////////
    //  Plots with widths(alpha)                                              //
    ////////////////////////////////////////////////////////////////////////////

    time(&start); // TIME

    // TFile WIDTHroot(outdir+"output/Width.root","RECREATE");

    // PLOT_WIDTH(widths_hist_data_SM, widths_hist_SM, gen_widths_hist_SM, outdir+"pdfy/widths/", eta_bins_edge_SM,Pt_bins_Central,Pt_bins_HF);
    // PLOT_WIDTH(widths_hist_data_FE, widths_hist_FE, gen_widths_hist_FE, outdir+"pdfy/widths/", eta_bins_edge_SM,Pt_bins_Central,Pt_bins_HF);

    if(true){
      PLOT_WIDTH_gr(data_correlated_graphs_SM, MC_correlated_graphs_SM, gen_correlated_graphs_SM, outdir+"ClosureTest/", eta_bins_edge_SM, Pt_bins_Central, Pt_bins_HF, false);
      PLOT_WIDTH_gr(data_correlated_graphs_FE, MC_correlated_graphs_FE, gen_correlated_graphs_FE, outdir+"ClosureTest/", eta_bins_edge_SM, Pt_bins_Central, Pt_bins_HF, true);
    }
    // WIDTHroot.Close();

    time(&end); // TIME
    time_taken = double(end - start);
    if(debug) cout << "Time taken by PLOT_WIDTH is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Plots with all points of widths(alpha)                                                                                               //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // PLOT_ALLWIDTHS(soft_widths_hist_data_SM, soft_widths_hist_SM, soft_gen_widths_hist_SM,outdir+"pdfy/widths/allPoints_");
    // PLOT_ALLWIDTHS(soft_widths_hist_data_FE, soft_widths_hist_FE, soft_gen_widths_hist_FE,outdir+"pdfy/widths/allPoints_");
       
    time(&start); // TIME

    TFile widthroot(outdir+"output/widths.root","RECREATE");
    if(false){ // Ignore to save time consumption
      for( unsigned int m = 0; m < widths_hist_SM.size(); m++ ){
        for( unsigned int p = 0; p < widths_hist_SM.at(m).size(); p++ ){
          widths_hist_SM.at(m).at(p) -> Write();
          gen_widths_hist_SM.at(m).at(p) -> Write();
          widths_hist_data_SM.at(m).at(p) -> Write();
        }
      }
    }
    for( unsigned int m = 0; m < widths_hist_FE.size(); m++ ){
      for( unsigned int p = 0; p < widths_hist_FE.at(m).size(); p++ ){
        widths_hist_FE.at(m).at(p) -> Write();
        gen_widths_hist_FE.at(m).at(p) -> Write();
        widths_hist_data_FE.at(m).at(p) -> Write();
      }
    }
    h_chi2_tot->Write();
    widthroot.Close();
       
    time(&end); // TIME
    time_taken = double(end - start);
    if(debug) cout << "Time taken by STORE_WIDTHS is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  }

  std::cout << "plot_all : " << plot_all << '\n';

  /////////////////////////////////////////////////////////////////////////////////////////
  // plot with JERs with NSC fit
  /////////////////////////////////////////////////////////////////////////////////////////
  std::vector< TH1F* > h_fitERR_uncorrelated_SM, h_fitERR_correlated_SM, h_fitERR_uncorrelated_FE, h_fitERR_correlated_FE;
  // TFile NSCroot(outdir+"output/NSC.root","RECREATE");
  PLOT_NSC(JER_uncorrelated_data_hist_SM,JER_uncorrelated_MC_hist_SM,JER_uncorrelated_scale_hist_SM,h_fitERR_uncorrelated_SM,outdir,eta_bins_edge_SM,false, false);
  PLOT_NSC(JER_correlated_data_hist_SM,JER_correlated_MC_hist_SM,JER_correlated_scale_hist_SM,h_fitERR_correlated_SM,outdir,eta_bins_edge_SM,false, true);
  PLOT_NSC(JER_uncorrelated_data_hist_FE,JER_uncorrelated_MC_hist_FE,JER_uncorrelated_scale_hist_FE,h_fitERR_uncorrelated_FE,outdir,eta_bins_edge_FE, true, false);
  PLOT_NSC(JER_correlated_data_hist_FE,JER_correlated_MC_hist_FE,JER_correlated_scale_hist_FE,h_fitERR_correlated_FE,outdir,eta_bins_edge_FE, true, true);
  // NSCroot.Close();

  //////////////////////////////////////////////////////////////////////////////////////////
  // Extract Widths
  //////////////////////////////////////////////////////////////////////////////////////////

  SIGtoTXT(JER_correlated_MC_hist_FE, eta_bins_edge_FE, outdir, "FE_correlated");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  time(&start); // TIME

  TFile Dijetroot(outdir+"output/dijet_balance_"+year+".root","RECREATE");

  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  for( unsigned int m = 0; m < JER_correlated_MC_hist_SM.size(); m++ ) cout << m << " " << JER_correlated_MC_hist_SM.at(m)->GetTitle() << endl;
  for( unsigned int m = 0; m < JER_correlated_MC_hist_SM.size(); m++ ){   g = TH1toTGraphAsymmErrors(m, eta_bins, JER_correlated_MC_hist_SM.at(m),   "MC",   "SM", "nominal"); g -> Write();}
  for( unsigned int m = 0; m < JER_correlated_data_hist_SM.size(); m++ ){ g = TH1toTGraphAsymmErrors(m, eta_bins, JER_correlated_data_hist_SM.at(m), "Data", "SM", "nominal"); g -> Write();}
  for( unsigned int m = 0; m < JER_correlated_MC_hist_FE.size(); m++ ){   g = TH1toTGraphAsymmErrors(m, eta_bins, JER_correlated_MC_hist_FE.at(m),   "MC",   "FE", "nominal"); g -> Write();}
  for( unsigned int m = 0; m < JER_correlated_data_hist_FE.size(); m++ ){ g = TH1toTGraphAsymmErrors(m, eta_bins, JER_correlated_data_hist_FE.at(m), "Data", "FE", "nominal"); g -> Write();}
  Dijetroot.Close();
       
  time(&end); // TIME
  time_taken = double(end - start);
  if(debug) cout << "Time taken by STORE_DIJET is : " << fixed << setprecision(2) << time_taken << " sec " << endl;

  if(false){ // Replaced by dijet_balance file
    TFile JERroot(outdir+"output/JERs.root","RECREATE");
    for( unsigned int m = 0; m < JER_uncorrelated_MC_hist_SM.size(); m++ ){           JER_uncorrelated_MC_hist_SM.at(m)           -> Write();}
    for( unsigned int m = 0; m < JER_uncorrelated_data_hist_SM.size(); m++ ){         JER_uncorrelated_data_hist_SM.at(m)         -> Write();}
    for( unsigned int m = 0; m < JER_uncorrelated_MC_hist_FE.size(); m++ ){           JER_uncorrelated_MC_hist_FE.at(m)           -> Write();}
    for( unsigned int m = 0; m < JER_uncorrelated_data_hist_FE.size(); m++ ){         JER_uncorrelated_data_hist_FE.at(m)         -> Write();}
    for( unsigned int m = 0; m < JER_uncorrelated_MC_hist_FE_control.size(); m++ ){   JER_uncorrelated_MC_hist_FE_control.at(m)   -> Write();}
    for( unsigned int m = 0; m < JER_uncorrelated_data_hist_FE_control.size(); m++ ){ JER_uncorrelated_data_hist_FE_control.at(m) -> Write();}
    for( unsigned int m = 0; m < JER_correlated_MC_hist_SM.size(); m++ ){             JER_correlated_MC_hist_SM.at(m)             -> Write();}
    for( unsigned int m = 0; m < JER_correlated_data_hist_SM.size(); m++ ){           JER_correlated_data_hist_SM.at(m)           -> Write();}
    for( unsigned int m = 0; m < JER_correlated_MC_hist_FE.size(); m++ ){             JER_correlated_MC_hist_FE.at(m)             -> Write();}
    for( unsigned int m = 0; m < JER_correlated_data_hist_FE.size(); m++ ){           JER_correlated_data_hist_FE.at(m)           -> Write();}
    for( unsigned int m = 0; m < JER_correlated_MC_hist_FE_control.size(); m++ ){     JER_correlated_MC_hist_FE_control.at(m)     -> Write();}
    for( unsigned int m = 0; m < JER_correlated_data_hist_FE_control.size(); m++ ){   JER_correlated_data_hist_FE_control.at(m)   -> Write();}
    JERroot.Close();
  }

  TFile gbis(outdir+"output/SFs.root","RECREATE");
  
  std::vector<double> SF_uncorrelated_SM, SF_uncorrelated_SM_error, SF_correlated_SM, SF_correlated_SM_error, SF_uncorrelated_FE, SF_uncorrelated_FE_error, SF_correlated_FE, SF_correlated_FE_error, eta_bin_SM_center, eta_bin_SM_error, eta_bin_FE_center, eta_bin_FE_error;
  std::vector<double> SF_uncorrelated_SM_ptdep_min, SF_correlated_SM_ptdep_min, SF_uncorrelated_FE_ptdep_min, SF_correlated_FE_ptdep_min;
  std::vector<double> SF_uncorrelated_SM_ptdep_max, SF_correlated_SM_ptdep_max, SF_uncorrelated_FE_ptdep_max, SF_correlated_FE_ptdep_max;
  ofstream texfile;
  texfile.open ("output/scalefactors_tex.txt");
  texfile << "standard method\n";
  SFtoTXT(texfile, JER_uncorrelated_scale_hist_SM, width_pt_SM, SF_uncorrelated_SM, SF_uncorrelated_SM_error, SF_uncorrelated_SM_ptdep_min, SF_uncorrelated_SM_ptdep_max, eta_bin_SM_center, eta_bin_SM_error, eta_bins_edge_SM, shift_, false, false);

  texfile << "\n forward extension \n";
  SFtoTXT(texfile, JER_uncorrelated_scale_hist_FE, width_pt_FE, SF_uncorrelated_FE, SF_uncorrelated_FE_error, SF_uncorrelated_FE_ptdep_min, SF_uncorrelated_FE_ptdep_max, eta_bin_FE_center, eta_bin_FE_error, eta_bins_edge_SM, shift_, true, false);

  texfile << "\n standard method, correlated fit\n";
  SFtoTXT(texfile, JER_correlated_scale_hist_SM, width_pt_SM, SF_correlated_SM, SF_correlated_SM_error, SF_correlated_SM_ptdep_min, SF_correlated_SM_ptdep_max, eta_bin_SM_center, eta_bin_SM_error, eta_bins_edge_SM, shift_, false, true);

  texfile << "\n forward extension, correlated fit \n";
  SFtoTXT(texfile, JER_correlated_scale_hist_FE, width_pt_FE, SF_correlated_FE, SF_correlated_FE_error, SF_correlated_FE_ptdep_min, SF_correlated_FE_ptdep_max, eta_bin_FE_center, eta_bin_FE_error, eta_bins_edge_SM, shift_, true, true);

  texfile << "\n";
  texfile.close();

  gbis.Close();
  ofstream txt_ST, txt_FE;
  txt_ST.open ("output/scalefactors_ST.txt");
  txt_FE.open ("output/scalefactors_FE.txt");
  // std::cout << "#eta_bin_SM_center" << " " << "eta_bin_SM_error" << " " << "SF_uncorrelated_SM" << " " << "SF_uncorrelated_SM_error" << " " << "SF_correlated_SM" << " " << "SF_correlated_SM_error" << "\n";
  txt_ST << "#eta_bin_SM_center" << " " << "eta_bin_SM_error" << " " << "SF_uncorrelated_SM" << " " << "SF_uncorrelated_SM_error" << " " << "SF_correlated_SM" << " " << "SF_correlated_SM_error" << "\n";
  txt_FE << "#eta_bin_FE_center" << " " << "eta_bin_FE_error" << " " << "SF_uncorrelated_FE" << " " << "SF_uncorrelated_FE_error" << " " << "SF_correlated_FE" << " " << "SF_correlated_FE_error" << "\n";

  // std::cout << eta_bin_SM_center.size() << " " << eta_bin_SM_error.size() << " " << SF_uncorrelated_SM.size() << " " << SF_uncorrelated_SM_error.size() << " " << SF_correlated_SM.size() << " " << SF_correlated_SM_error.size() << " " << SF_uncorrelated_SM_ptdep_min.size() << " " << SF_uncorrelated_SM_ptdep_max.size() << " " << SF_correlated_SM_ptdep_min.size() << " " << SF_correlated_SM_ptdep_max.size() << "\n";
  if (debug) std::cout << "\n === Print ScaleFactors SM:\n" << endl;
  for (unsigned int i = 0; i < JER_uncorrelated_scale_hist_SM.size(); i++) {
    txt_ST << eta_bin_SM_center[i] << " " << eta_bin_SM_error[i] << " " << SF_uncorrelated_SM[i] << " " << SF_uncorrelated_SM_error[i] << " " << SF_correlated_SM[i] << " " << SF_correlated_SM_error[i] << " " << SF_uncorrelated_SM_ptdep_min[i] << " " << SF_uncorrelated_SM_ptdep_max[i] << " " << SF_correlated_SM_ptdep_min[i] << " " << SF_correlated_SM_ptdep_max[i] << "\n";
    if (debug) std::cout << eta_bin_SM_center[i] << " " << eta_bin_SM_error[i] << " " << SF_uncorrelated_SM[i] << " " << SF_uncorrelated_SM_error[i] << " " << SF_correlated_SM[i] << " " << SF_correlated_SM_error[i] << " " << SF_uncorrelated_SM_ptdep_min[i] << " " << SF_uncorrelated_SM_ptdep_max[i] << " " << SF_correlated_SM_ptdep_min[i] << " " << SF_correlated_SM_ptdep_max[i] << "\n";
  }
  if (debug) std::cout << "\n === Print ScaleFactors FE:\n" << endl;
  for (unsigned int i = 0; i < JER_uncorrelated_scale_hist_FE.size(); i++) {
    txt_FE << eta_bin_FE_center[i] << " " << eta_bin_FE_error[i] << " " << SF_uncorrelated_FE[i] << " " << SF_uncorrelated_FE_error[i] << " " << SF_correlated_FE[i] << " " << SF_correlated_FE_error[i] << " " << SF_uncorrelated_FE_ptdep_min[i] << " " << SF_uncorrelated_FE_ptdep_max[i] << " " << SF_correlated_FE_ptdep_min[i] << " " << SF_correlated_FE_ptdep_max[i] << "\n";
    if (debug) std::cout << eta_bin_FE_center[i] << " " << eta_bin_FE_error[i] << " " << SF_uncorrelated_FE[i] << " " << SF_uncorrelated_FE_error[i] << " " << SF_correlated_FE[i] << " " << SF_correlated_FE_error[i] << " " << SF_uncorrelated_FE_ptdep_min[i] << " " << SF_uncorrelated_FE_ptdep_max[i] << " " << SF_correlated_FE_ptdep_min[i] << " " << SF_correlated_FE_ptdep_max[i] << "\n";
  }
  txt_ST.close();
  txt_FE.close();

  TCanvas* canv_SF = tdrCanvas("JER SF",eta_bins_edge_SM[0]-0.1, eta_bins_edge_SM[EtaBins_SM + EtaBins_FE]+0.5, 0.8, 3.0, "#eta", "JER SF");
  canv_SF->SetTickx(0);
  canv_SF->SetTicky(0);
  TGraphErrors *gr_st = new TGraphErrors(JER_correlated_scale_hist_SM.size(), &(eta_bin_SM_center[0]), &(SF_correlated_SM[0]), &(eta_bin_SM_error[0]), &(SF_correlated_SM_error[0]));
  TGraphErrors *gr_fe = new TGraphErrors(JER_correlated_scale_hist_FE.size(), &(eta_bin_FE_center[0]), &(SF_correlated_FE[0]), &(eta_bin_FE_error[0]), &(SF_correlated_FE_error[0]));
  tdrDraw(gr_st, "P", kFullDotLarge, kRed);
  tdrDraw(gr_fe, "P", kFullDotLarge, kBlue);

  TLegend *leg = tdrLeg(0.70,0.75,0.9,0.9);
  leg->AddEntry(gr_st, "Standard","lep");
  leg->AddEntry(gr_fe, "Forward ext", "lep");
  leg->Draw("same");

  canv_SF->SaveAs(outdir+"output/SF.pdf");
  // canv_SF->SaveAs(outdir+"output/SF.png");
  // canv_SF->SaveAs(outdir+"output/SF.root");

  //////////////////////////////////////////////////////////////////////////////////////////
  //    SF plots overlayed with ...                                                       //
  //////////////////////////////////////////////////////////////////////////////////////////
  // TFile SFsoverlayed(outdir+"output/SFsoverlayed.root","RECREATE");

  PLOT_SF(JER_uncorrelated_scale_hist_SM, JER_correlated_scale_hist_SM, JER015_scale_hist_SM, h_fitERR_correlated_SM, outdir+"pdfy/SFs/", eta_bins_edge_SM, false, JER_correlated_scale_hist_SM);
  PLOT_SF(JER_uncorrelated_scale_hist_FE, JER_correlated_scale_hist_FE, JER015_scale_hist_FE, h_fitERR_correlated_FE, outdir+"pdfy/SFs/", eta_bins_edge_FE, true, JER_correlated_scale_hist_FE);

  // Program takes long to close after return. Ran multiple cleaning options and calculated the time taken for the full program to estimate single steps. Keep.
  time(&endp); // TIME
  time_taken = double(endp - startp);
  printf("Close Program; time taken %5.2f\n",time_taken);
  return true;

}
