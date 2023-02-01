#include <iostream>
#include <fstream>
#include <TChain.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>
#include <TStyle.h>
#include <string>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TProfile.h>

#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMinuit.h>
#include <TMatrixD.h>
#include "constants.h"
#include "common_binning.hpp"

#include <bits/stdc++.h>

bool debug = true;
bool alpha_new = true;
bool extenda = true;
bool useRMS = true;

int binHF = 15;

int printEta = 15;
int printPt = 6;

int sizeHist = 3000;

TString g_year = "";
TString g_study = "";

static std::vector<double> usedPtTrigger_central; // Set in mainRun.cxx
static std::vector<double> usedPtTrigger_forward; // Set in mainRun.cxx

// Code of Andrea Malara
// Based on code by Marek Niedziela, Matthias Schröder, Kristin Goebel

// Smeared MC instead of data:
//bool smeared = true;
TString labeldata = "Data";
//TString labeldata = "Smeared MC";
TString sep = "|";

// Functions
TString fNSC      = "TMath::Sqrt( ([0]*[0]/(x*x))+[1]*[1]/x+[2]*[2] )";
TString fNSCsign  = "TMath::Sqrt( (TMath::Sign(1, [0])*[0]*[0]/(x*x))+[1]*[1]/x+[2]*[2] )";
TString fNSxPC    = "TMath::Sqrt([0]*abs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])";
TString fNSxPCdiv = "TMath::Sqrt([0]*abs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2]/100)";

struct fit_data{
  // x, y values:
  std::vector<double> x_val, y_val;
  // variance-covariance matrix for y values:
  TMatrixD y_cov;
  // inverted cov matrix; calculated by chi2_linear "on demand".
  TMatrixD y_cov_inv;

  void reset() {
    x_val.clear();
    y_val.clear();
    y_cov.ResizeTo(0,0);
    y_cov_inv.ResizeTo(0,0);
  }

  void CheckPoints() {
    std::vector<int> RemovedPoints;
    TMatrixD y_cov_new;
    int j = 0;

    for(unsigned int i = 0; i < y_val.size(); i++) {
      if ( y_val.at(i) == 0) {
        x_val.erase(x_val.begin()+i);
        y_val.erase(y_val.begin()+i);
        RemovedPoints.push_back(j);
        i = i-1;
      }
      j++;
    }

    y_cov_new.ResizeTo(x_val.size(),x_val.size());
    for(unsigned int i=0; i < x_val.size(); i++) {
      for(unsigned int k= 0; k < x_val.size(); k++) {
        y_cov_new(i,k) = y_cov(i+RemovedPoints.size(),k+RemovedPoints.size());
      }
    }
    y_cov.ResizeTo(0,0);
    y_cov.ResizeTo(x_val.size(),x_val.size());
    y_cov = y_cov_new;
  }
};

fit_data data_;

TString dtos(double number, int precision);

void SetBinsToZero(TGraphErrors* graph, double bin);
void ExtractFitParameters(ofstream& FitPar, TF1* fit, TString info);
void RemoveAsymBinsFromJER(TH1F* hist, double bin);
void histLoadAsym( TFile &f, bool data, TString text, std::vector< std::vector< std::vector< TH1F* > > > &Asymmetry, std::vector< std::vector< std::vector< TH1F* > > > &GenAsymmetry, int etaBins, int ptBins, int AlphaBins, int etaShift);
void histMeanPt( std::vector< std::vector< std::vector< TH1F* > > > &Asymmetry , std::vector< std::vector< std::vector< double > > > &Widths );
void histWidthAsym_old( std::vector<std::vector<std::vector<TH1F*> > > &Asymmetry , std::vector<std::vector<std::vector<double> > > &Widths, std::vector<std::vector<std::vector<double> > > &WidthsError , bool fill_all );
void histWidthAsym( std::vector< std::vector< std::vector< TH1F* > > > &Asymmetry , std::vector< std::vector< std::vector< double > > > &Widths, std::vector< std::vector< std::vector< double > > > &WidthsError , bool fill_all, double alpha, bool isFE, std::vector< std::vector< std::vector< double > > > &lower_x, std::vector< std::vector< std::vector< double > > > &upper_x, std::vector <double> eta_bins, std::vector <double> alpha_bins);
void histWidthMCTruth( std::vector<std::vector<std::vector<TH1F*> > > &Asymmetry , std::vector<std::vector<std::vector<double> > > &Widths, std::vector<std::vector<std::vector<double> > > &WidthsError );
void fill_widths_hists( TString name1, std::vector< std::vector< TH1F* > > &widths , std::vector< std::vector< std::vector< double > > > Widths, std::vector< std::vector< std::vector< double > > > WidthsError);
void histLinFit( std::vector< std::vector< TH1F* > > widths_hist_all , std::vector< std::vector< double > > &Widths, std::vector< std::vector< double > > &WidthsError, bool isFE );
void histLinCorFit( std::vector< std::vector< std::vector< double > > > Widths, std::vector< std::vector< std::vector< double > > > WidthsError, std::vector< std::vector< TGraphErrors* > > &output_graph, std::vector< std::vector< double > > &output, std::vector< std::vector< double > > &output_error, bool isFE, bool isMC, TH1F* h_chi2_tot);
void widths_015_ratios( TString name1, std::vector<TH1F*> &widths, std::vector<std::vector<std::vector<double> > > Widths, std::vector<std::vector<std::vector<double> > > WidthsError, std::vector<std::vector<std::vector<double> > > WidthsTwo, std::vector<std::vector<std::vector<double> > > WidthsTwoError, std::vector<std::vector<std::vector<double> > > forward_width_pt );
void correctJERwithPLI(std::vector< std::vector< double > > &Output, std::vector< std::vector< double > > &OutputError, std::vector< std::vector< double > > Widths, std::vector< std::vector< double > > WidthsError, std::vector< std::vector< double > > PLI, std::vector< std::vector< double > > PLIError, float shift = 0.0);
void correctJERwithPLI015(std::vector<std::vector<double> > &Output, std::vector<std::vector<double> > &OutputError, std::vector<std::vector<std::vector<double> > > Widths, std::vector<std::vector<std::vector<double> > > WidthsError, std::vector<std::vector<std::vector<double> > > PLI, std::vector<std::vector< std::vector< double > > > PLIError, float shift = 0.0);
// void correctForRef( TString name1, std::vector<std::vector<double> > &Output, std::vector<std::vector<double> > &OutputError, std::vector<std::vector<double> > Input, std::vector<std::vector<double> > InputError, std::vector<std::vector<std::vector<double> > > width_pt, int shift, TString outdir);
void correctForRef( TString name1, std::vector<std::vector<double> > &Output, std::vector<std::vector<double> > &OutputError, std::vector<std::vector<double> > Input, std::vector<std::vector<double> > InputError, std::vector<std::vector<double> > Input2, std::vector<std::vector<double> > InputError2, std::vector<std::vector<std::vector<double> > > width_pt, int shift, TString outdir);
void makeScales( std::vector< std::vector< double > > &Output, std::vector< std::vector< double > > &OutputError, std::vector< std::vector< double > > Input1, std::vector< std::vector< double > > Input1Error, std::vector< std::vector< double > > Input2, std::vector< std::vector< double > > Input2Error );
void fill_mctruth_hist( TString name1, std::vector< TH1F* > &output, std::vector< std::vector< std::vector< double > > > Widths, std::vector< std::vector< std::vector< double > > > WidthsError, std::vector< std::vector< std::vector< double > > > pt_binning, double range);
void fill_hist( TString name1, std::vector< TH1F* > &output, std::vector< std::vector< double > > Widths, std::vector< std::vector< double > > WidthsError, std::vector< std::vector< std::vector< double > > > pt_binning, double range, int shift2 = 0);
void Fill_Map3D(std::vector< std::vector < std::vector < TH1F* > > > &Asymmetry, std::vector < TH2F* > &Map, std::vector < double > &eta_bins, std::vector < double > &pt_bins );
void make_lin_fit(double & slope, double & d_slope, double & offset, double & d_offset, double min_slope, double max_slope, double min_offset, double max_offset, double & chi2);
void chi2_linear(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* p, Int_t status);
void chi2_calculation(Double_t& fval, Double_t* p);
void fitLin( TH1F &hist, double &width, double &error );
double sumSquare(double a, double b);
double findMinMax(TH1F* JER, std::vector< std::vector< double > > pt_width, TF1* NSC_ratio, TF1* constfit, bool isMin, bool print=false);
double removePointsforAlphaExtrapolation(bool isFE, double eta, int p);
std::vector<double> Confidence(TH1F* hist, double confLevel);
bool removePointsforFit(bool isFE, int m, int p);
inline bool extendAlpha(int p);

TGraphErrors* TH1toTGraphError(TH1* hist, double extend_err);
TGraphAsymmErrors* TH1toTGraphAsymmErrors(int m, std::vector <double> eta_bins, TH1* hist, TString sample, TString method, TString var="nominal");
TH1F* GetMCFitRatio(TH1F* hist, TF1* fit, int color, int MarkerSize);
TF1* GetFitsRatio(TF1* fit1, TF1* fit2, double ptmin, int color, TString func);
string BinToString(float low, float high, int precision=3);

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TString dtos(double number, int precision)
{
  stringstream stream;
  stream << std::fixed << std::setprecision(precision) << number;
  return stream.str();
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TGraphErrors* TH1toTGraphError(TH1* hist, double extend_err){
  vector<double> xvalues = {};
  vector<double> yvalues = {};
  vector<double> yerrors = {};
  vector<double> dummy = {};
  for(unsigned int i=1; i<=hist->GetNbinsX(); i ++){
    if(hist->GetBinContent(i)==0) continue;
    // cout << hist->GetTitle() << "   " << hist->GetBinContent(i) << endl;
    xvalues.push_back(hist->GetXaxis()->GetBinCenter(i));
    yvalues.push_back(hist->GetBinContent(i));
    yerrors.push_back(hist->GetBinError(i)+extend_err);
    dummy.push_back(0);
  }
  // cout << " " << xvalues.size() << " " << yvalues.size() << " " << yerrors.size() << " " << dummy.size() << endl;
  TGraphErrors* MC_graph = new TGraphErrors(xvalues.size(), &xvalues[0], &yvalues[0], &dummy[0], &yerrors[0]);
  return MC_graph;
}

string BinToString(float low, float high, int precision=3) {
  stringstream ss;
  ss <<  left << setfill('0') << setw(precision) << low;
  string low_ = ss.str();
  if (low==0) low_ = "0."+low_;
  replace(low_.begin(), low_.end(), '.', 'p');
  ss.str("");
  ss.clear();
  ss <<  left << setfill('0') << setw(precision) << high;
  string high_ = ss.str();
  if (high==0) high_ = "0."+high_;
  replace(high_.begin(), high_.end(), '.', 'p');
  return low_+"_"+high_;
}

TGraphAsymmErrors* TH1toTGraphAsymmErrors(int m, std::vector <double> eta_bins, TH1* hist, TString sample, TString method, TString var){
  vector<double> xvalues = {};
  vector<double> xerrors_lo = {};
  vector<double> xerrors_hi = {};
  vector<double> yvalues = {};
  vector<double> yerrors_lo = {};
  vector<double> yerrors_hi = {};
  vector<double> dummy = {};
  std::vector<double> usedPtBinning = m<14?usedPtTrigger_central:usedPtTrigger_forward;
  for(unsigned int i=2; i<=hist->GetNbinsX(); i ++){
    if(hist->GetBinContent(i)==0) continue;
    // cout << hist->GetTitle() << "   " << hist->GetBinContent(i) << endl;
    double pT    = hist->GetXaxis()->GetBinCenter(i);
    double pTerr = hist->GetBinError(i);
    xvalues.push_back(pT);
    yvalues.push_back(hist->GetBinContent(i));
    yerrors_lo.push_back(pTerr);
    yerrors_hi.push_back(pTerr);
  }
  // if(debug) cout << " x: " << setw(8) << xvalues.size() << " " << setw(8) << xerrors_lo.size() << " " << setw(8) << xerrors_hi.size() << " | y: " << setw(12) << yvalues.size() << " " << setw(12) << yerrors_lo.size() << " " << setw(12) << yerrors_hi.size() << endl;
  for(unsigned int p=0; p<xvalues.size();p++){
    xerrors_lo.push_back(xvalues[p]-usedPtBinning[p]);
    xerrors_hi.push_back(usedPtBinning[p+1]-xvalues[p]);
    // if(debug) cout << " x:" << setw(8) << xvalues[p] << setw(8) << xerrors_lo[p] << setw(8) << xerrors_hi[p] << " | y:" << setw(12) << yvalues[p] << setw(12) << yerrors_lo[p] << setw(12) << yerrors_hi[p] << " | pt:" << setw(8) << usedPtBinning[p] << setw(8) << usedPtBinning[p+1] << endl;
  }
  TGraphAsymmErrors* MC_graph = new TGraphAsymmErrors(xvalues.size(), &xvalues[0], &yvalues[0], &xerrors_lo[0], &xerrors_hi[0], &yerrors_lo[0], &yerrors_hi[0]);
  // TString names = "TGraph_"+(TString) hist->GetName();
  TString bins = BinToString(eta_bins[m], eta_bins[m+1], 5);
  TString names = "dijet_balance_jer_"+sample+"_";
  names += bins+"_";
  names += method+"_";
  names += var;
  // cout << names << endl;
  MC_graph->SetNameTitle(names, names);
  return MC_graph;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TH1F* GetMCFitRatio(TH1F* hist, TF1* fit, int color, int MarkerSize = 1){
  TH1F* h_ratio = (TH1F*) hist->Clone();
  // delete h_ratio->GetListOfFunctions()->FindObject("mcFIT");
  // delete h_ratio->GetFunction("mcFIT");
  for(unsigned int bin = 1; bin <= hist->GetNbinsX(); bin++){
    double pt    = hist->GetXaxis()->GetBinCenter(bin);
    double width = hist->GetBinContent(bin);
    double value = fit->Eval(pt);
    double ratio = width/value;
    h_ratio->SetBinContent(bin, ratio);
  }
  h_ratio->SetMarkerSize(MarkerSize);
  h_ratio->SetMarkerStyle(kFullCircle);

  return h_ratio;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void RemoveAsymBinsFromJER(TH1F* hist, double bin){
  int count_bins = 0;
  for(unsigned int i=0; i<hist->GetNbinsX(); i++)
  {
    if(hist->GetBinContent(i)==0) continue;
    // cout << i << " " << hist->GetXaxis()->GetBinCenter(i) << endl;
    count_bins++;
  }
  // cout << count_bins << endl;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void ExtractFitParameters(ofstream& FitPar, TF1* fit, TString info){
  FitPar << info << "\n";
  if(fit==0) return;
  FitPar << "N = " << setw(8) << fit->GetParameter(0) << " \u00B1 " << setw(8) << fit->GetParError(0) << "\n";
  FitPar << "S = " << setw(8) << fit->GetParameter(1) << " \u00B1 " << setw(8) << fit->GetParError(1) << "\n";
  FitPar << "C = " << setw(8) << fit->GetParameter(2) << " \u00B1 " << setw(8) << fit->GetParError(2) << "\n";
  if(info.Contains("NSxPC")) FitPar << "P = " << setw(8) << fit->GetParameter(3) << " \u00B1 " << setw(8) << fit->GetParError(3) << "\n";
  FitPar << "\u03A7 = " << setw(8) << fit->GetChisquare() << "\n";
  FitPar << "NDF = " << setw(3) << fit->GetNDF() << "\n";
  FitPar << "\n";
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void histLoadAsym( TFile &f, bool data, TString text, std::vector< std::vector< std::vector< TH1F* > > > &Asymmetry, std::vector< std::vector< std::vector< TH1F* > > > &GenAsymmetry, int etaBins, int ptBins, int AlphaBins, int etaShift) {
  for( int m = etaShift; m < etaBins+etaShift; m++ ) {
    std::vector< std::vector< TH1F* > > temp2;
    std::vector< std::vector< TH1F* > > temp2gen;
    for( int p = 0; p < ptBins; p++ ) {
      std::vector< TH1F* > temp1;
      std::vector< TH1F* > temp1gen;
      int alpha = 6;
      if(extendAlpha(p)) alpha = 13;
      // for( int r = 0; r < AlphaBins; r++ ) {
      for( int r = 0; r < alpha; r++ ) {
        TString name = text;                        name    += "_eta"; name     += m+1; name    += "_pt"; name    += p+1; name    += "_alpha"; name     += r+1;
        TString namegen = "gen_"; namegen += text;  namegen += "_eta"; namegen  += m+1; namegen += "_pt"; namegen += p+1; namegen += "_alpha"; namegen  += r+1;
        TH1F* h = (TH1F*)f.Get(name);
        if (h) h->Rebin(2);
        temp1.push_back(h);
        if (( TH1F* )f.Get( name ) == 0) {
          std::cout << "[ERROR] " << name << std::endl;
        }
        if ( data == false ) {
          TH1F* h_gen = (TH1F*)f.Get(namegen);
          if (h_gen) h_gen->Rebin(2);
          temp1gen.push_back(h_gen);
        }
      }
      temp2.push_back(temp1);
      if ( data == false ) temp2gen.push_back(temp1gen);
    }
    Asymmetry.push_back(temp2);
    if ( data == false ) GenAsymmetry.push_back(temp2gen);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================


void histMeanPt( std::vector< std::vector< std::vector< TH1F* > > > &Asymmetry , std::vector< std::vector< std::vector< double > > > &Widths ) {
  for( unsigned int m = 0; m < Asymmetry.size(); m++ ) {
    std::vector< std::vector< double > > temp1;
    for( unsigned int p = 0; p < Asymmetry.at(m).size(); p++ ) {
      std::vector< double > temp2;
      for( unsigned int r = 0; r < Asymmetry.at(m).at(p).size(); r++ ) {
        temp2.push_back( (*Asymmetry.at(m).at(p).at(r)).GetMean() );
      }
      temp1.push_back(temp2);
    }
    Widths.push_back(temp1);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

std::vector<double> Confidence(TH1F* hist, double confLevel)
{
  // sigmas = [0.99, 0.98, 0.95, 0.87, 0.68]
  // conversions = [2.576,2.326,1.960,1.514,0.9945]
  // for perc, conv in zip(sigmas,conversions):
  // width = Confidence(hist, confLevel = perc)/(2* conv)

  // 99% confidence is [-2.5sigma,+2.5sigma] so divide by 2*2.576
  // 98% confidence is [-2.3sigma,+2.3sigma] so divide by 2*2.326
  // 95% confidence is [-2sigma,+2sigma] so divide by 2*1.960
  // 87% confidence is [-1.5*sigma,+1.5sigma] so divide by 2*1.514
  // 68% confidence is [-sigma,+sigma] so divide by 2*0.9945

  double ix = hist->GetXaxis()->FindBin(hist->GetMean());
  double ixlow = ix;
  double ixhigh = ix;
  double nb = hist->GetNbinsX();
  double ntot = hist->Integral();
  double nsum = hist->GetBinContent(ix);
  double width = hist->GetBinWidth(ix);
  double error = 0;
  if (ntot==0) return {0,0,0,0};
  while (nsum < confLevel * ntot){
    double nlow = ixlow>0?hist->GetBinContent(ixlow-1):0;
    double nhigh = ixhigh<nb?hist->GetBinContent(ixhigh+1):0;
    if (nsum+max(nlow,nhigh) < confLevel * ntot){
      if (nlow>=nhigh && ixlow>0){
        nsum += nlow;
        ixlow -=1;
        width += hist->GetBinWidth(ixlow);
      }
      else if (ixhigh<nb){
        nsum += nhigh;
        ixhigh+=1;
        width += hist->GetBinWidth(ixhigh);
      }
      else{
        throw std::runtime_error("BOOM");
      }
    }
    else{
      if (nlow>nhigh){
        width += hist->GetBinWidth(ixlow-1) * (confLevel * ntot - nsum) / nlow;
      }
      else{
        width += hist->GetBinWidth(ixhigh+1) * (confLevel * ntot - nsum) / nhigh;
      }
      nsum = ntot;
    }
  }
  std::vector<double> results = {width, hist->GetRMSError(), ixlow, ixhigh};
  return results;
  // return width;
}

void histWidthAsym( std::vector< std::vector< std::vector< TH1F* > > > &Asymmetry , std::vector< std::vector< std::vector< double > > > &Widths, std::vector< std::vector< std::vector< double > > > &WidthsError , bool fill_all, double alpha, bool isFE, std::vector< std::vector< std::vector< double > > > &lower_x, std::vector< std::vector< std::vector< double > > > &upper_x, std::vector <double> eta_bins, std::vector <double> alpha_bins) {
  std::vector<double> asyms;
  double asym;
  double asymerr;
  std::map<double, double> conversions = { {0.99, 2.576}, {0.985, 2.446}/*TODO*/, {0.98, 2.326}, {0.95, 1.960}, {0.87, 1.514}, {0.68, 0.9945} };
  // if(isFE && !fill_all) cout << endl << " === NEW ASYMMETRY ===" << endl;

  for( unsigned int m = 0; m < Asymmetry.size(); m++ ) {
    std::vector< std::vector< double > > temp2;
    std::vector< std::vector< double > > temp_error2;
    std::vector< std::vector< double > > temp2_lower_x, temp2_upper_x;
    for( unsigned int p = 0; p < Asymmetry.at(m).size() ; p++ ) {
      std::vector< double > temp1;
      std::vector< double > temp_error1;
      std::vector< double > temp1_lower_x, temp1_upper_x;
      // std::cout << setw(5) << m << setw(5) << p << setw(5);
      for( unsigned int r = 0; r < Asymmetry.at(m).at(p).size(); r++ ) {
        TH1F* temp_hist = (TH1F*) Asymmetry.at(m).at(p).at(r)->Clone();
        double low, up;
        if (temp_hist->Integral() > 0) {
          for (int i = 0; i <= temp_hist->GetNbinsX(); i++) {
            if (i < temp_hist->FindBin(-0.5) || i > temp_hist->FindBin(0.5)) {
              temp_hist->SetBinContent(i,0);
            }
          }

          if(useRMS){
            temp_hist->ComputeIntegral();
            Double_t xq[2], yq[2];
            xq[0] = std::min(alpha, 1.-alpha);
            xq[1] = std::max(alpha, 1.-alpha);
            temp_hist->GetQuantiles(2, yq, xq);
            temp_hist->GetXaxis()->SetRange(temp_hist->FindBin(yq[0]), temp_hist->FindBin(yq[1]));
            asym = temp_hist->GetRMS();
            asymerr = temp_hist->GetRMSError();
            low = temp_hist->GetBinCenter(temp_hist->FindBin(yq[0])); up = temp_hist->GetBinCenter(temp_hist->FindBin(yq[1]));
          }
          else{
            asyms = Confidence(temp_hist, alpha);
            asym = asyms[0]/(2*conversions[alpha]);
            asymerr = asyms[1];
            low = temp_hist->GetBinCenter(temp_hist->FindBin(asyms[2]));
            up = temp_hist->GetBinCenter(temp_hist->FindBin(asyms[3]));
          }

        } else { asym = 0.; asymerr = 0.; low = 10.0; up = 10.0;};

        if (!fill_all) {
          if (alpha_bins.at(r)<removePointsforAlphaExtrapolation(isFE, eta_bins.at(m), p+1)) {asym = 0.; asymerr = 0.;}
        }

        // cout << setw(10) << asym << " (" << setw(2) << r << ")";

        temp1.push_back( asym );
        temp_error1.push_back( asymerr );
        temp1_lower_x.push_back(low);
        temp1_upper_x.push_back(up);
      }
      // cout << endl;
      temp2.push_back(temp1);
      temp_error2.push_back(temp_error1);
      temp2_lower_x.push_back(temp1_lower_x);
      temp2_upper_x.push_back(temp1_upper_x);
    }
    Widths.push_back(temp2);
    WidthsError.push_back(temp_error2);
    lower_x.push_back(temp2_lower_x);
    upper_x.push_back(temp2_upper_x);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void histWidthMCTruth( std::vector<std::vector<std::vector<TH1F*> > > &Asymmetry , std::vector<std::vector<std::vector<double> > > &Widths, std::vector<std::vector<std::vector<double> > > &WidthsError ) {
  double asym;
  double asymerr;

  for( unsigned int m = 0; m < Asymmetry.size(); m++ ) {
    std::vector< std::vector< double > > temp2;
    std::vector< std::vector< double > > temp_error2;
    for( unsigned int p = 0; p < Asymmetry.at(m).size() ; p++ ) {
      std::vector< double > temp1;
      std::vector< double > temp_error1;
      for( unsigned int r = 0; r < Asymmetry.at(m).at(p).size(); r++ ) {
        asym    = ((*Asymmetry.at(m).at(p).at(r)).GetRMS())/TMath::Sqrt(2);
        asymerr = ((*Asymmetry.at(m).at(p).at(r)).GetRMSError())/TMath::Sqrt(2);
        temp1.push_back( asym );
        temp_error1.push_back( asymerr );
      }
      temp2.push_back(temp1);
      temp_error2.push_back(temp_error1);
    }
    Widths.push_back(temp2);
    WidthsError.push_back(temp_error2);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void fill_widths_hists( TString name1, std::vector< std::vector< TH1F* > > &widths , std::vector< std::vector< std::vector< double > > > Widths, std::vector< std::vector< std::vector< double > > > WidthsError) {
  double aMax[] = { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 , 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1. };
  double temp;
  for( unsigned int m = 0; m < Widths.size(); m++ ) {
    std::vector< TH1F* > temp1;
    for( unsigned int p = 0; p < Widths.at(m).size(); p++ ) {
      TString name_width = name1;
      name_width += "_eta"; name_width += m+1; name_width += "_pt"; name_width += p+1;
      // if (name1.Contains("_fe")) { name_width += m+2; }
      // else { name_width += m+1; }
      TH1F *h1 = new TH1F( name_width, name_width, 105, 0, 1.05);
      h1 ->GetYaxis()->SetTitle("#sigma_{A}");	h1 ->GetXaxis()->SetTitle("#alpha_{max}"); h1 ->GetYaxis()->SetTitleOffset(1.);
      h1 -> Sumw2();
      // cout << setw(3) << m << setw(3) << p;
      for( unsigned int r = 0; r < Widths.at(m).at(p).size(); r++ ) {
        // cout << setw(3) << r << " (" << setw(3) << h1 -> FindBin( aMax[ r ] ) << ", " << setw(10) << Widths.at(m).at(p).at(r) << ")";
        temp = Widths.at(m).at(p).at(r);
        if ( !(TMath::IsNaN(temp)) && temp != 0) h1 -> SetBinContent( h1 -> FindBin( aMax[ r ] ), Widths.at(m).at(p).at(r) );
        if ( !(TMath::IsNaN(temp)) && temp != 0) h1 -> SetBinError( h1 -> FindBin( aMax[ r ] ), WidthsError.at(m).at(p).at(r) );
      }
      // cout << endl;
      h1 ->GetYaxis()-> SetRangeUser( 0., 0.3 );
      temp1.push_back(h1);
    }
    widths.push_back(temp1);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void histLinFit( std::vector< std::vector< TH1F* > > widths_hist_all , std::vector< std::vector< double > > &Widths, std::vector< std::vector< double > > &WidthsError, bool isFE ) {
  for( unsigned int m = 0; m < widths_hist_all.size(); m++ ) {
    std::vector<double> temp2;
    std::vector<double> temp_error2;
    for( unsigned int p = 0; p < widths_hist_all.at(m).size(); p++ ) {
      double value, error;
      if ( widths_hist_all.at(m).at(p)->GetEntries() != 0 ) {
        fitLin( *( widths_hist_all.at(m).at(p) ), value, error );

        // if (removePointsforFit(isFE, m, p)) {value = 0; error = 0;}

        temp2.push_back(value);
        temp_error2.push_back(error);
      } else {
        temp2.push_back(0.);
        temp_error2.push_back(0.);
      }
    }
    Widths.push_back(temp2);
    WidthsError.push_back(temp_error2);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void createCov( std::vector< double > Widths, std::vector< double > WidthsError, std::vector<float> alpha, TMatrixD &y_cov_mc){
  for(unsigned int ialpha=0; ialpha < alpha.size(); ++ialpha) {
    for (unsigned int jalpha =0; jalpha < alpha.size(); jalpha++) {
      if ( ialpha <= jalpha ) {
        double n1_mc = pow(Widths.at(ialpha),2)/(2*pow(WidthsError.at(ialpha),2));
        double n2_mc = pow(Widths.at(jalpha),2)/(2*pow(WidthsError.at(jalpha),2));
        y_cov_mc(ialpha, jalpha) = pow(WidthsError.at(ialpha),2) * pow((n1_mc/n2_mc),2)*
        (Widths.at(ialpha)/Widths.at(jalpha));
      } else {
        double n1_mc = pow(Widths.at(jalpha),2)/(2*pow(WidthsError.at(jalpha),2));
        double n2_mc = pow(Widths.at(ialpha),2)/(2*pow(WidthsError.at(ialpha),2));
        y_cov_mc(ialpha, jalpha) = pow(WidthsError.at(jalpha),2) * pow((n1_mc/n2_mc),2)*
        (Widths.at(jalpha)/Widths.at(ialpha));
      }
    }
  }
}

int CheckNumberAlphaPoints( std::vector< double > Widths){
  int length = 0;
  for(auto a: Widths) length = length+((a==0)?0:1);
  return length;
}

bool IsGreaterZero(double i){
  return 0.<i;
}

void SetFit(vector<double> Widths, vector<double> alpha, TMatrixD y_cov_, TF1 *lin_extrapol, int m, int p, bool isFE, bool isMC){
  data_.reset();
  data_.x_val = alpha;
  data_.y_val = Widths;
  data_.y_cov.ResizeTo(alpha.size(), alpha.size());
  data_.y_cov = y_cov_;
  data_.CheckPoints();
  // choose start values for the fit
  // double slope = (Widths.at(Widths.size()-1) - Widths.at(Widths.size()-3))/(x.at(x.size()-1) - x.at(x.size()-3));
  // double offset = Widths.at(Widths.size()-1) - (slope*x.at(x.size()-1));
  double chi2 = -10.0;
  int ndf = Widths.size()-2;
  double slope  = 0.15;
  double offset = 0.05;
  double d_slope = slope;
  double d_offset = offset;
  double min_slope = 0.05;
  double max_slope = 0.5;
  double min_offset = 0.001;
  double max_offset = 0.15;
  if ( !isFE && m == 7 && p == 2) { min_slope = 0.15; max_offset = 0.08;}
  if ( !isFE && m == 7 && p == 6) { min_slope = 0.01; max_offset = 0.08;}
  if ( !isFE && m == 8 && p == 6) { min_slope = 0.09; max_offset = 0.045;}
  if ( !isFE && m >= 10&& p == 2) { max_offset = 0.4;}
  if (isFE&&isMC&&m==10&& p == 6) { min_slope = 0.13; max_offset = 0.06;}
  if (!isFE&&isMC&&m==7&& p == 6) { min_slope = 0.01; max_offset = 0.08; min_offset = 0.04;}
  if (!isFE&&isMC&&m==7&& p == 2) { min_slope = 0.15; max_offset = 0.08; min_offset = 0.06;}

  //         std::cout << "eta: " << m <<  ", p_T: " << p << std::endl;
  //         std::cout << "fit start values: " << "slope: " << slope << " offset: " << offset << std::endl;
  make_lin_fit(slope, d_slope, offset, d_offset, min_slope, max_slope, min_offset, max_offset, chi2);
  //         std::cout << "fit values: " << "slope: " << slope << " offset: " << offset << std::endl;
  // if(isFE&&isMC) cout << "During: " << offset << endl;
  lin_extrapol->SetParameter(0, offset);
  lin_extrapol->SetParError(0, d_offset);
  lin_extrapol->SetParameter(1, slope);
  lin_extrapol->SetParError(1, d_slope);
  lin_extrapol->SetChisquare(chi2);//TODO: set the correct chi/2
  lin_extrapol->SetNDF(ndf);
}

inline bool extendAlpha(int p){
  // return ( extenda && (p<5) );
  return ( extenda );
  // return ( extenda && (p<10 || (20<p && p<25) || 30<p) );
}

void histLinCorFit( std::vector< std::vector< std::vector< double > > > Widths, std::vector< std::vector< std::vector< double > > > WidthsError, std::vector< std::vector< TGraphErrors* > > &output_graph, std::vector< std::vector< double > > &output, std::vector< std::vector< double > > &output_error, bool isFE, bool isMC, TH1F* h_chi2_tot) {
  std::vector<float> alpha;
  // std::vector<float> alpha;
  // alpha.push_back(0.05); alpha.push_back(0.1); alpha.push_back(0.15); alpha.push_back(0.20); alpha.push_back(0.25); alpha.push_back(0.3);
  for( unsigned int m = 0; m < Widths.size(); m++ ) {
    // eta loop
    std::vector<TGraphErrors*> temp2_graph;
    std::vector<double> temp2;
    std::vector<double> temp_error2;

    for( unsigned int p = 0; p < Widths.at(m).size(); p++ ) {
      // p_T loop
      std::vector<double> x,x_e;

      alpha = { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3  };
      if( extendAlpha(p) ) alpha = { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1. };

      // cout << setw(5) << m << setw(5) << p;
      for(int ialpha=0; ialpha < alpha.size(); ++ialpha) {
      // for(int ialpha=0; ialpha < AlphaBins; ++ialpha) {
      // cout << setw(5) << ialpha;
        x.push_back(alpha.at(ialpha));
        x_e.push_back(0.);
      }
      // cout << endl;

      vector<double> widths     = Widths.at(m).at(p);
      vector<double> widths_err = WidthsError.at(m).at(p);

      TMatrixD y_cov_mc;
      // y_cov_mc.ResizeTo(AlphaBins, AlphaBins);
      y_cov_mc.ResizeTo(alpha.size(), alpha.size());

      TMatrixD y_cov_mc_first;
      vector<double> widths_first = widths; vector<double> widths_first_err = widths_err;
      if(2<CheckNumberAlphaPoints(widths_first)){
        vector<double>::iterator it = std::find_if(widths_first.begin(), widths_first.end(), IsGreaterZero);
        int index = it-widths_first.begin();
        widths_first[index] = 0.; widths_first_err[index] = 0.;
      }
      // y_cov_mc_first.ResizeTo(AlphaBins, AlphaBins);
      y_cov_mc_first.ResizeTo(alpha.size(), alpha.size());

      TMatrixD y_cov_mc_last;
      vector<float> alpha_last = alpha; vector<double> x_last = x;
      vector<double> widths_last = widths; vector<double> widths_last_err = widths_err;
      if(2<CheckNumberAlphaPoints(widths_last)){
        widths_last.pop_back(); widths_last_err.pop_back(); alpha_last.pop_back(); x_last.pop_back();
      }
      // y_cov_mc_last.ResizeTo(AlphaBins-1, AlphaBins-1);
      y_cov_mc_last.ResizeTo(alpha_last.size(), alpha_last.size());

      createCov(widths, widths_err, alpha, y_cov_mc);
      createCov(widths_first, widths_first_err, alpha, y_cov_mc_first);
      createCov(widths_last, widths_last_err, alpha_last, y_cov_mc_last);


      //create TGraphErrors from previously defined vectors
      TGraphErrors* extrapol_MC = new TGraphErrors(alpha.size(),&x[0],&widths.at(0),&x_e[0],&widths_err.at(0));
      TString name = "Graph_SM_eta";
      if (isFE) name = "Graph_FE_eta";
      name += m+1; name+="_pt"; name += p+1;
      extrapol_MC->SetName(name);

      // fit linear extrapolation function
      TF1 *lin_extrapol_mc = new TF1("lin_extrapol_mc","[0]+[1]*x",0, 0.35); // alpha.back()+0.05
      TF1 *lin_extrapol_mc_first = new TF1("lin_extrapol_mc_first","[0]+[1]*x",0, 0.35);
      TF1 *lin_extrapol_mc_last = new TF1("lin_extrapol_mc_last","[0]+[1]*x",0, 0.35);

      TF1 *log_extrapol_mc = new TF1("log_extrapol_mc","[0]+[1]*log(x)",0,alpha.back()+0.05);
      TF1 *pol2_extrapol_mc = new TF1("pol2_extrapol_mc","pol2(0)",0,alpha.back()+0.05);
      TF1 *pol3_extrapol_mc = new TF1("pol3_extrapol_mc","pol3(0)",0,alpha.back()+0.05);

      // void SetFit(vector<double> Widths, vector<double> alpha, TMatrixD y_cov_, TF1 *lin_extrapol, int m, int p, bool isFE, bool isMC){
      SetFit(widths,       x,      y_cov_mc,       lin_extrapol_mc,       m, p, isFE, isMC);
      SetFit(widths_first, x,      y_cov_mc_first, lin_extrapol_mc_first, m, p, isFE, isMC);
      SetFit(widths_last,  x_last, y_cov_mc_last,  lin_extrapol_mc_last,  m, p, isFE, isMC);

      SetFit(widths,       x,      y_cov_mc,       log_extrapol_mc,       m, p, isFE, isMC);
      SetFit(widths,       x,      y_cov_mc,       pol2_extrapol_mc,      m, p, isFE, isMC);
      SetFit(widths,       x,      y_cov_mc,       pol3_extrapol_mc,      m, p, isFE, isMC);

      extrapol_MC->GetListOfFunctions()->Add(lin_extrapol_mc);
      extrapol_MC->GetListOfFunctions()->Add(log_extrapol_mc);
      extrapol_MC->GetListOfFunctions()->Add(pol2_extrapol_mc);
      extrapol_MC->GetListOfFunctions()->Add(pol3_extrapol_mc);

      double offset_nom = lin_extrapol_mc->GetParameter(0);
      double offset_first = lin_extrapol_mc_first->GetParameter(0);
      double offset_last = lin_extrapol_mc_last->GetParameter(0);
      double offset = ( offset_nom + offset_first + offset_last)/3;

      double diff_nom = abs(offset-offset_nom);
      double diff_first = abs(offset-offset_first);
      double diff_last = abs(offset-offset_last);

      double err_stat = lin_extrapol_mc->GetParError(0);
      double err_offset = 0.;
      if(diff_nom<diff_first && diff_last<diff_first) err_offset = diff_first;
      else if(diff_first<diff_nom && diff_last<diff_nom) err_offset = diff_nom;
      else if(diff_nom<diff_last && diff_first<diff_last) err_offset = diff_last;
      else err_offset = err_stat; // When errors are really small down to 0; They are not considered anyway

      double err_total = sumSquare(err_offset, err_stat);

      h_chi2_tot->Fill(lin_extrapol_mc->GetChisquare()/lin_extrapol_mc->GetNDF()); // TODO for new method

      if(alpha_new)
      {
        extrapol_MC->GetListOfFunctions()->Add(lin_extrapol_mc_first);
        extrapol_MC->GetListOfFunctions()->Add(lin_extrapol_mc_last);
      }

      temp2.push_back(alpha_new?offset:offset_nom);
      temp_error2.push_back(alpha_new?err_total:err_stat);
      temp2_graph.push_back(extrapol_MC);
    }
    output.push_back(temp2);
    output_error.push_back(temp_error2);
    output_graph.push_back(temp2_graph);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void widths_015_ratios( TString name1, std::vector<TH1F*> &widths, std::vector<std::vector<std::vector<double> > > Widths, std::vector<std::vector<std::vector<double> > > WidthsError, std::vector<std::vector<std::vector<double> > > WidthsTwo, std::vector<std::vector<std::vector<double> > > WidthsTwoError, std::vector<std::vector<std::vector<double> > > forward_width_pt ) {
  double temp, tempError;
  for( unsigned int m = 0; m < Widths.size(); m++ ) {
    TString name_width = name1;
    name_width += "_eta"; name_width += m+1;
    TH1F *hist = new TH1F( name_width, name_width, sizeHist, 0, sizeHist );
    hist ->GetYaxis()->SetTitle();
    hist ->GetXaxis()->SetTitle("p_{T} [GeV]");
    hist -> GetYaxis()->SetRangeUser( 0, 2. );
    for( unsigned int p = 0; p < Widths.at(m).size(); p++ ) {
      if (WidthsTwo.at(m).at(p).at(2) != 0.) {
        temp = Widths.at(m).at(p).at(2)/WidthsTwo.at(m).at(p).at(2);
        // tempError = WidthsError.at(m).at(p).at(2)/WidthsTwo.at(m).at(p).at(2) + ( Widths.at(m).at(p).at(2) * WidthsTwoError.at(m).at(p).at(2) ) / ( WidthsTwo.at(m).at(p).at(2) * WidthsTwo.at(m).at(p).at(2) ) ;
        tempError = sumSquare(WidthsError.at(m).at(p).at(2)/WidthsTwo.at(m).at(p).at(2) , ( Widths.at(m).at(p).at(2) * WidthsTwoError.at(m).at(p).at(2) ) / ( WidthsTwo.at(m).at(p).at(2) * WidthsTwo.at(m).at(p).at(2) ) );
        hist -> SetBinContent( hist -> FindBin( forward_width_pt.at(m).at(p).at(2) ), temp );
        hist -> SetBinError( hist -> FindBin( forward_width_pt.at(m).at(p).at(2) ), tempError );
      }
    }
    widths.push_back(hist);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void correctJERwithPLI(std::vector< std::vector< double > > &Output, std::vector< std::vector< double > > &OutputError, std::vector< std::vector< double > > Widths, std::vector< std::vector< double > > WidthsError, std::vector< std::vector< double > > PLI, std::vector< std::vector< double > > PLIError, float shift) {


  for( unsigned int i = 0; i < Widths.size(); i++ ) {
    std::vector< double > temp2;
    std::vector< double > temp_error2;
    for( unsigned int j = 0; j < Widths.at(i).size(); j++ ) {
      double temp;
      double temp_error;
      if (Widths.at(i).at(j) != 0. ) {
        // With PLI correction (change also alpha=015):
        temp = TMath::Sqrt(2)*TMath::Sqrt( Widths.at(i).at(j) * Widths.at(i).at(j) - (1.+shift)*PLI.at(i).at(j) * PLI.at(i).at(j) );
        // temp_error = ( Widths.at(i).at(j) * WidthsError.at(i).at(j) + (1.+shift)* PLI.at(i).at(j) * PLIError.at(i).at(j) )/temp;
        temp_error = TMath::Sqrt(2)*sumSquare( Widths.at(i).at(j) * WidthsError.at(i).at(j), (1.+shift)*PLI.at(i).at(j) * PLIError.at(i).at(j) )/temp;
      } else {
        temp = 0.;
        temp_error = 0.;
      }
      if ( TMath::IsNaN(temp) ) { temp = 0. ; temp_error = 0.; }

      temp2.push_back(temp);
      temp_error2.push_back(temp_error);
    }
    Output.push_back(temp2);
    OutputError.push_back(temp_error2);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void correctJERwithPLI015(std::vector<std::vector<double> > &Output, std::vector<std::vector<double> > &OutputError, std::vector<std::vector<std::vector<double> > > Widths, std::vector<std::vector<std::vector<double> > > WidthsError, std::vector<std::vector<std::vector<double> > > PLI, std::vector<std::vector< std::vector< double > > > PLIError, float shift) {
  for( unsigned int i = 0; i < Widths.size(); i++ ) {
    std::vector< double > temp2;
    std::vector< double > temp_error2;
    for( unsigned int j = 0; j < Widths.at(i).size(); j++ ) {
      double temp;
      double temp_error;
      // With PLI correction (change also alpha=0):
      temp = TMath::Sqrt(2)*TMath::Sqrt( Widths.at(i).at(j).at(2) * Widths.at(i).at(j).at(2) - (1.+shift)* PLI.at(i).at(j).at(2) * PLI.at(i).at(j).at(2) );
      // temp_error = ( Widths.at(i).at(j).at(2) * WidthsError.at(i).at(j).at(2) + (1.+shift) * PLI.at(i).at(j).at(2) * PLIError.at(i).at(j).at(2) )/temp;
      temp_error = TMath::Sqrt(2)*sumSquare( Widths.at(i).at(j).at(2) * WidthsError.at(i).at(j).at(2) , (1.+shift) * PLI.at(i).at(j).at(2) * PLIError.at(i).at(j).at(2) )/temp;
      if ( TMath::IsNaN(temp) ) { temp = 0 ; temp_error = 0; }
      temp2.push_back(temp);
      temp_error2.push_back(temp_error);
    }
    Output.push_back(temp2);
    OutputError.push_back(temp_error2);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

// void correctForRef( TString name1, std::vector<std::vector<double> > &Output, std::vector<std::vector<double> > &OutputError, std::vector<std::vector<double> > Input, std::vector<std::vector<double> > InputError, std::vector<std::vector<std::vector<double> > > width_pt, int shift, TString outdir) {
void correctForRef( TString name1, std::vector<std::vector<double> > &Output, std::vector<std::vector<double> > &OutputError, std::vector<std::vector<double> > Input, std::vector<std::vector<double> > InputError, std::vector<std::vector<double> > Input2, std::vector<std::vector<double> > InputError2, std::vector<std::vector<std::vector<double> > > width_pt, int shift, TString outdir) {
  double Ref, Probe, RefError, ProbeError, pT;
  if(debug) cout << "In correctForRef with " << name1 << endl;

  TH1F *hist = new TH1F( name1+"_hist", name1+"_hist", sizeHist, 0, sizeHist );
  hist ->GetYaxis()->SetTitle("#sigma_{A}");
  hist ->GetXaxis()->SetTitle("p_{T}");

  for( unsigned int p = 0; p < Input.at(0).size(); p++ ) {
    double temp = 0;
    double temp_error = 0;
    for (unsigned int i = 0; i < shift; i++) {
      // temp += Input.at(i).at(p)/shift;
      // temp_error += InputError.at(i).at(p)/shift;
      temp += Input2.at(i).at(p)/shift;
      temp_error += InputError2.at(i).at(p)/shift;
    }
    if ( temp != 0 && !(TMath::IsNaN(temp))) {
      pT = (double)(*std::max_element(width_pt.at(0).at(p).begin(),width_pt.at(0).at(p).end()));
      hist->SetBinContent(hist->FindBin(pT), temp);
      hist->SetBinError(  hist->FindBin(pT), temp_error );
    }
  }

  if(debug) cout << "In correctForRef | LINE " << __LINE__ << endl;
  TCanvas* canv = new TCanvas("nscplot","nscplot",50,50,800,600);
  hist -> Draw();
  canv -> Print(outdir+"pdfy/JERs/reference"+name1+".pdf","pdf");
  canv -> Print(outdir+"pdfy/JERs/reference"+name1+".png","png");

  for( unsigned int m = 1; m < shift; m++ ) {
    std::vector< double > temp2;
    std::vector< double > temp_error2;

    for( unsigned int p = 0; p < Input.at(m).size(); p++ ) {
      double temp;
      double temp_error;
      if ( Input.at(m).at(p) != 0. ) {
        pT = (double)(*std::max_element(width_pt.at(0).at(p).begin(),width_pt.at(0).at(p).end()));
        // Other pT value ? average? minimum? ...

        temp = Input.at(m).at(p);
        temp_error = InputError.at(m).at(p);

        if ( !(TMath::IsNaN(temp)) ) temp2.push_back(temp);
        else temp2.push_back(0.);
        if ( !(TMath::IsNaN(temp)) ) temp_error2.push_back(temp_error);
        else temp_error2.push_back(0.);
      }
      else { temp2.push_back(0.); temp_error2.push_back(0.); }
    }
    Output.push_back(temp2);
    OutputError.push_back(temp_error2);
  }

  if(debug) cout << "In correctForRef | LINE " << __LINE__ << endl;
  for( unsigned int m = shift; m < Input.size() ; m++ ) {
    vector<double> pTbinsValue = (m<binHF)?usedPtTrigger_central:usedPtTrigger_forward;
    vector<double> pTbinsRef = usedPtTrigger_central;
    std::vector< double > temp2;
    std::vector< double > temp_error2;

    for( unsigned int p = 0; p < Input.at(m).size(); p++ ) {
      double temp;
      double temp_error;

      if ( Input.at(m).at(p) != 0. ) {
        // pT = width_pt.at(m).at(p).at(5) ;
        int index_p = -1; // In order to throw break
        double ptvalue = (pTbinsValue[p]+pTbinsValue[p+1])/2; // 1st: Get mid of pT bin of HF bins
        for(unsigned int z=0; z<pTbinsRef.size()-1; z++){
          if(pTbinsRef[z]<ptvalue&&ptvalue<pTbinsRef[z+1]){index_p = z; break;}
        }
        pT = (double)(*std::max_element(width_pt.at(0).at(index_p).begin(),width_pt.at(0).at(index_p).end()));

        Ref = hist->GetBinContent(hist->FindBin(pT));
        RefError = hist->GetBinError(hist->FindBin(pT));
        Probe = Input.at(m).at(p);
        ProbeError = InputError.at(m).at(p);

        temp = TMath::Sqrt( 2*Probe*Probe - Ref*Ref);
        temp_error = sumSquare(2*Probe*ProbeError, Ref*RefError )/temp;

        // Remove ptbins by Hand; How to automate? TODO
        // bool removeByHand = (
        //   (m==1&&p==0)|| // 0.3-0.5
        //   (m==8&&p==27)|| // 1.9-2.0
        //   (m==9&&(p==0||p==1))|| // 2.0-2.2
        //   (m==13&&p==17)|| // 2.7-2.9
        //   (m==16&&p==12) // 3.1-3.5
        // );
        bool removeByHand = false;
        if( g_year.Contains("UL16post") && g_study.Contains("common_fine_v1") ){
          if(m==6){
            int dummy = 0;
            // if(p==30) removeByHand = true;
          } else if (m==7){
            if(p==30||p==29) removeByHand = true;
          } else if (m==8){
            if(p==30||p==29||p==28) removeByHand = true;
          } else if (m==9){
            if(p==30||p==29||p==28) removeByHand = true;
          } else if (m==10){
            if(p==30||p==29||p==28||p==27) removeByHand = true;
          } else if (m==11){
            if(p==30||p==29||p==28||p==27||p==26||p==25) removeByHand = true;
          } else if (m==12){
            if(p==30||p==29||p==28||p==27||p==26||p==25||p==24) removeByHand = true;
          } else if (m==13){
            if(p==30||p==29||p==28||p==27||p==26||p==25||p==24||p==23||p==22||p==21) removeByHand = true;
          } else if (m==16){
            if(p==13) removeByHand = true;
          }
        }
        if( g_year.Contains("UL18") && g_study.Contains("common_fine_v1") ){
          if(m==6){
            if(p==30) removeByHand = true;
          } else if (m==7){
            if(p==30) removeByHand = true;
          } else if (m==8){
            if(p==30||p==29||p==28) removeByHand = true;
          } else if (m==9){
            if(p==30||p==29||p==28) removeByHand = true;
          } else if (m==10){
            if(p==30||p==29||p==28||p==27) removeByHand = true;
          } else if (m==11){
            if(p==30||p==29||p==28||p==27||p==26||p==25) removeByHand = true;
          } else if (m==12){
            if(p==30||p==29||p==28||p==27||p==26||p==25||p==24) removeByHand = true;
          } else if (m==13){
            if(p==30||p==29||p==28||p==27||p==26||p==25||p==24||p==23||p==22||p==21) removeByHand = true;
          } else if (m==16){
            if(p==13) removeByHand = true;
          }
        }

        if(removeByHand){
          // cout << "Remove eta bin " << m+1 << " | " << p+1 << " from JER by Hand" << endl;
          if(debug) cout << "Remove eta bin " << m << " | " << p+1 << " from JER by Hand" << endl;
          temp2.push_back(0.);
          temp_error2.push_back(0.);
        }
        else{
          if ( !(TMath::IsNaN(temp)) ) temp2.push_back(temp);
          else temp2.push_back(0.);
          if ( !(TMath::IsNaN(temp)) ) temp_error2.push_back(temp_error);
          else temp_error2.push_back(0.);
        }
      }
      else { temp2.push_back(0.); temp_error2.push_back(0.); }
    }
    Output.push_back(temp2);
    OutputError.push_back(temp_error2);
  }
  delete canv;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void makeScales( std::vector< std::vector< double > > &Output, std::vector< std::vector< double > > &OutputError, std::vector< std::vector< double > > Input1, std::vector< std::vector< double > > Input1Error, std::vector< std::vector< double > > Input2, std::vector< std::vector< double > > Input2Error ) {
  for( unsigned int i = 0; i < Input1.size(); i++ ) {
    std::vector< double > temp2;
    std::vector< double > temp_error2;
    for( unsigned int j = 0; j < Input1.at(i).size(); j++ ) {
      double temp;
      double temp_error;
      // Cut on scale factors not between 0.7 and 2
      //if ( (Input1.at(i).at(j)!=0.) && (Input2.at(i).at(j)!=0.) && Input1.at(i).at(j) / Input2.at(i).at(j)>0.7 && Input1.at(i).at(j) / Input2.at(i).at(j)<2. )
      if ( (Input1.at(i).at(j)!=0.) && (Input2.at(i).at(j)!=0.) ) {
        temp = Input1.at(i).at(j) / Input2.at(i).at(j);
        // temp_error = Input1Error.at(i).at(j) / Input2.at(i).at(j) + ( Input1.at(i).at(j) * Input2Error.at(i).at(j)) / ( Input2.at(i).at(j) * Input2.at(i).at(j) ) ;
        temp_error = sumSquare( Input1Error.at(i).at(j)/Input2.at(i).at(j) , ( Input1.at(i).at(j) * Input2Error.at(i).at(j)) / ( Input2.at(i).at(j) * Input2.at(i).at(j) ) ) ;
      } else { temp = 0.; temp_error = 0.; }
      if ( TMath::IsNaN(temp) ) { temp = 0.; temp_error = 0.; }
      temp2.push_back(temp);
      temp_error2.push_back(temp_error);
    }
    Output.push_back(temp2);
    OutputError.push_back(temp_error2);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void fill_mctruth_hist( TString name1, std::vector< TH1F* > &output, std::vector< std::vector< std::vector< double > > > Widths, std::vector< std::vector< std::vector< double > > > WidthsError, std::vector< std::vector< std::vector< double > > > pt_binning, double range) {
  for( unsigned int m = 0; m <  Widths.size(); m++ ) {
    TString name_width = name1; name_width += m+1;
    TH1F *h1 = new TH1F( name_width, name_width, sizeHist, 0, sizeHist );
    h1 ->GetYaxis()->SetTitle("#sigma_{MCTruth}");	h1 ->GetXaxis()->SetTitle("p_{T}");	h1 -> Sumw2();

    for( unsigned int p = 0; p <  Widths.at(m).size(); p++ ) {
      double pT = (double)(*std::max_element(pt_binning.at(m).at(p).begin(),pt_binning.at(m).at(p).end()));
      if ( ( !(TMath::IsNaN(Widths.at(m).at(p).at(5))) )      && Widths.at(m).at(p).at(5)!= 0. ) h1 -> SetBinContent( h1 -> FindBin(pT), Widths.at(m).at(p).at(5)      );
      if ( ( !(TMath::IsNaN(WidthsError.at(m).at(p).at(5))) ) && Widths.at(m).at(p).at(5)!= 0. ) h1 -> SetBinError(   h1 -> FindBin(pT), WidthsError.at(m).at(p).at(5) );
    }
    h1 ->GetYaxis()-> SetRangeUser( 0., range );
    output.push_back(h1);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void fill_hist( TString name1, std::vector< TH1F* > &output, std::vector< std::vector< double > > Widths, std::vector< std::vector< double > > WidthsError, std::vector< std::vector< std::vector< double > > > pt_binning, double range, int shift2) {
  // int shift = 0;
  // if (name1.Contains("SF_")) {
  //   shift = 1;
  //   if (Widths.size() > 10) shift = 3;
  // }
  for( unsigned int m = 0; m <  Widths.size(); m++ ) {
    TString name = name1;
    name += m+1;
    TH1F *h1 = new TH1F( name, name, sizeHist, 0, sizeHist );
    // std::cout << m << " " << Widths.at(m).size() << std::endl;;
    if (name1.Contains("SF_")) h1 ->GetYaxis()->SetTitle("Scale factor");
    else                       h1 ->GetYaxis()->SetTitle("#sigma_{JER}");
    h1 ->GetXaxis()->SetTitle("p_{T}");	h1 -> Sumw2();

    for( unsigned int p = 0; p <  Widths.at(m).size(); p++ ) {
      double pT = (double)(*std::max_element(pt_binning.at(m+shift2).at(p).begin(),pt_binning.at(m+shift2).at(p).end()));
      if ( ( !(TMath::IsNaN(Widths.at(m).at(p))) ) && Widths.at(m).at(p)!= 0. )      h1 -> SetBinContent(h1 -> FindBin(pT), Widths.at(m).at(p) );
      if ( ( !(TMath::IsNaN(WidthsError.at(m).at(p))) ) && Widths.at(m).at(p)!= 0. ) h1 -> SetBinError(  h1 -> FindBin(pT), WidthsError.at(m).at(p) );

      if (name1.Contains("SF_") && h1->GetBinContent(h1->FindBin(pT))<1) h1->SetBinContent(h1->FindBin(pT),1);
    }
    h1 ->GetYaxis()-> SetRangeUser( 0., range );
    output.push_back(h1);
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void Fill_Map3D(std::vector< std::vector < std::vector < TH1F* > > > &Asymmetry, std::vector < TH2F* > &Map, std::vector < double > &eta_bins, std::vector < double > &pt_bins ) {
  for (unsigned int r = 0; r < 6; r++) {
    TString name = "Map_mean_"; name += r;
    TH2F* temp = new TH2F(name,name, pt_bins.size(), pt_bins.at(0), pt_bins.at(pt_bins.size()-1), eta_bins.size(), eta_bins.at(0), eta_bins.at(eta_bins.size()-1));
    // TH2F* temp = new TH2F(name,name, eta_bins.size(), &eta_bins[0], pt_bins.size(), &pt_bins[0]);
    Map.push_back(temp);
  }

  for (unsigned int m = 0; m < Asymmetry.size(); m++) {
    for (unsigned int p = 0; p < Asymmetry.at(m).size(); p++) {
      for (unsigned int r = 0; r < Asymmetry.at(m).at(p).size(); r++) {
        Map.at(r) -> SetBinContent( Map.at(r) -> FindBin(pt_bins[p]), Map.at(r) -> FindBin(eta_bins[m]), Asymmetry.at(m).at(p).at(r)->GetMean());
      }
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void make_lin_fit(double & slope, double & d_slope, double & offset, double & d_offset, double min_slope, double max_slope, double min_offset, double max_offset, double & chi2) {
  TMinuit min;
  min.SetPrintLevel(-1);
  //min.SetPrintLevel(0);
  if (slope  < 0.05  || slope  > 0.5) slope  = 0.15;
  if (offset < 0.001 || offset > 0.1) offset = 0.05;
  int err = min.DefineParameter(0, "slope", slope, d_slope, min_slope, max_slope);
  assert(err==0);
  err = min.DefineParameter(1, "offset", offset, d_offset, min_offset, max_offset);
  assert(err==0);
  min.SetFCN(chi2_linear);
  min.mnmigr();
  min.GetParameter(0, slope, d_slope);
  min.GetParameter(1, offset, d_offset);
  Double_t par[2]; par[0] = slope; par[1] = offset;
  chi2_calculation(chi2, par);
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void chi2_linear(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* p, Int_t status) {
  if (data_.y_cov_inv.GetNcols()==0) {
    double dummy;
    int ncols = data_.y_cov.GetNcols();
    data_.y_cov_inv.ResizeTo(ncols, ncols);
    data_.y_cov_inv = data_.y_cov.Invert(&dummy);
  }
  const size_t ndata = data_.x_val.size(); // number of data points in x,y graph to fit to
  std::vector<double> delta_y(ndata);
  for(size_t i=0; i<ndata; ++i) {
    delta_y[i] = data_.x_val[i]*p[0] + p[1] - data_.y_val[i];
  }
  // now calculate the chi2, i.e.
  //  dy^T * C^{-1} * dy
  // where C is the variance--covariance matrix and dy = (y_data - y_pred)
  // This could probably be implemented in ROOT, but it's so simple, we just do it here:
  fval = 0.0;
  for(size_t i=0; i<ndata; ++i) {
    for(size_t j=0; j<ndata; ++j) {
      fval += delta_y[i] * delta_y[j] * data_.y_cov_inv(i,j);
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void chi2_calculation(Double_t& fval, Double_t* p) {
  if (data_.y_cov_inv.GetNcols()==0) {
    double dummy;
    int ncols = data_.y_cov.GetNcols();
    data_.y_cov_inv.ResizeTo(ncols, ncols);
    data_.y_cov_inv = data_.y_cov.Invert(&dummy);
  }
  const size_t ndata = data_.x_val.size(); // number of data points in x,y graph to fit to
  std::vector<double> delta_y(ndata);
  for(size_t i=0; i<ndata; ++i) {
    delta_y[i] = data_.x_val[i]*p[0] + p[1] - data_.y_val[i];
  }
  // now calculate the chi2, i.e.
  //  dy^T * C^{-1} * dy
  // where C is the variance--covariance matrix and dy = (y_data - y_pred)
  // This could probably be implemented in ROOT, but it's so simple, we just do it here:
  fval = 0.0;
  for(size_t i=0; i<ndata; ++i) {
    for(size_t j=0; j<ndata; ++j) {
      fval += delta_y[i] * delta_y[j] * data_.y_cov_inv(i,j);
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

double sumSquare(double a, double b) {
  return TMath::Sqrt(a*a + b*b);
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void findExtreme(std::vector<TH1*> vec, double *x_min, double *x_max, double *y_min, double *y_max) {
  std::vector<double> x;
  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetXaxis()->GetXmin());}
  *x_min = *std::min_element(x.begin(), x.end());
  x.clear();

  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetXaxis()->GetXmax());}
  *x_max = *std::max_element(x.begin(), x.end());
  x.clear();

  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetMinimum());}
  *y_min = *std::min_element(x.begin(), x.end());
  x.clear();

  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetMaximum());}
  *y_max = *std::max_element(x.begin(), x.end());
  x.clear();

  if ((*x_min) == (*x_max)) {
    *x_min = (*x_min)*0.9;
    *x_max = (*x_max)*1.2;
  }

  if ((*y_min) == (*y_max)) {
    *y_min = *y_min*0.9;
    *y_max = *y_max*1.2;
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void findExtreme2(std::vector<TH1*> vec, double *x_min, double *x_max, double *y_min, double *y_max) {
  std::vector<double> x;
  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetXaxis()->GetXmin());}
  *x_min = *std::min_element(x.begin(), x.end());
  x.clear();

  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetXaxis()->GetXmax());}
  *x_max = *std::max_element(x.begin(), x.end());
  x.clear();

  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetMinimum());}
  *y_min = *std::min_element(x.begin(), x.end());
  x.clear();

  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetBinContent(vec.at(i)->GetMaximumBin()));}
  *y_max = *std::max_element(x.begin(), x.end());
  x.clear();

  if ((*x_min) == (*x_max)) {
    *x_min = (*x_min)*0.9;
    *x_max = (*x_max)*1.2;
  }

  if ((*y_min) == (*y_max)) {
    *y_min = *y_min*0.9;
    *y_max = *y_max*1.2;
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void findExtreme2(std::vector<TGraphErrors*> vec, double *x_min, double *x_max, double *y_min, double *y_max) {
  std::vector<double> x;
  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(TMath::MinElement(vec.at(i)->GetN(),vec.at(i)->GetX()));}
  *x_min = *std::min_element(x.begin(), x.end());
  x.clear();

  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(TMath::MaxElement(vec.at(i)->GetN(),vec.at(i)->GetX()));}
  *x_max = *std::max_element(x.begin(), x.end());
  x.clear();

  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(TMath::MinElement(vec.at(i)->GetN(),vec.at(i)->GetY()));}
  *y_min = *std::min_element(x.begin(), x.end());
  x.clear();

  // for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetBinContent(vec.at(i)->GetMaximumBin()));}
  // for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetHistogram()->GetMaximum());}
  // for (int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetHistogram()->GetBinContent(vec.at(i)->GetHistogram()->GetMaximumBin()));}
  for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(TMath::MaxElement(vec.at(i)->GetN(),vec.at(i)->GetY()));}
  *y_max = *std::max_element(x.begin(), x.end());
  x.clear();

  if ((*x_min) == (*x_max)) {
    *x_min = (*x_min)*0.9;
    *x_max = (*x_max)*1.2;
  }

  if ((*y_min) == (*y_max)) {
    *y_min = *y_min*0.9;
    *y_max = *y_max*1.2;
  }

  if (*x_min >= *x_max ||  *y_min >= *y_max)
  {
    cout << "Aiuto" << endl;
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

template<typename TT>
void findExtreme_gr(std::vector<TT*> vec, double *x_min, double *x_max, double *y_min, double *y_max) {
  std::vector<double> x;
  // if (std::is_same<TT, TH1F>::value) for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetXaxis()->GetXmin());}
  if (std::is_same<TT, TGraphErrors>::value) for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(TMath::MinElement(vec.at(i)->GetN(),vec.at(i)->GetX()));}
  *x_min = *std::min_element(x.begin(), x.end());
  x.clear();

  // if (std::is_same<TT, TH1F>::value) for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetXaxis()->GetXmax());}
  if (std::is_same<TT, TGraphErrors>::value) for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(TMath::MaxElement(vec.at(i)->GetN(),vec.at(i)->GetX()));}
  *x_max = *std::max_element(x.begin(), x.end());
  x.clear();


  // if (std::is_same<TT, TH1F>::value) for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetMinimum());}
  if (std::is_same<TT, TGraphErrors>::value) for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(TMath::MinElement(vec.at(i)->GetN(),vec.at(i)->GetY()));}
  *y_min = *std::min_element(x.begin(), x.end());
  x.clear();

  // if (std::is_same<TT, TH1F>::value) for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(vec.at(i)->GetBinContent(vec.at(i)->GetMaximumBin()));}
  if (std::is_same<TT, TGraphErrors>::value) for (unsigned int i = 0; i < vec.size(); i++) {x.push_back(TMath::MaxElement(vec.at(i)->GetN(),vec.at(i)->GetY()));}
  *y_max = *std::max_element(x.begin(), x.end());
  x.clear();

  if ((*x_min) == (*x_max)) {
    *x_min = (*x_min)*0.9;
    *x_max = (*x_max)*1.2;
  }

  if ((*y_min) == (*y_max)) {
    *y_min = *y_min*0.9;
    *y_max = *y_max*1.2;
  }

}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

double findMinMax(TH1F* JER, std::vector< std::vector< double > > pt_width, TF1* NSC_ratio, TF1* constfit, bool isMin, bool print) {
  double min = 10000;
  double max = 0;
  for (unsigned int p = 2; p < pt_width.size(); p++) { // TODO It's set to 2 just because in the following steps pt>2"bin are used
  double pT = (double)(*std::max_element(pt_width.at(p).begin(),pt_width.at(p).end()));
  if (JER->GetBinContent(JER->FindBin(pT))== 0.) continue;
  min = std::min(min, TMath::Abs(NSC_ratio->Eval(pT) - constfit->Eval(pT)));
  max = std::max(min, TMath::Abs(NSC_ratio->Eval(pT) - constfit->Eval(pT)));
  if (print){
    cout << "Max pT " << pT << " and JER " << JER->GetBinContent(JER->FindBin(pT));
    if(isMin) cout << " min " << min << endl;
    else      cout << " max " << max << endl;
  }
}
if (debug) cout << endl;
if (isMin) return min;
else return max;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void fitLin( TH1F &hist, double &width, double &error ) {
  TF1 * linfit = new TF1( "linfit", "[0]+x*[1]", 0, 0.3 );
  linfit -> SetParameter( 0, 0.1 );
  linfit -> SetParameter( 1, 0.1 );
  linfit -> SetParLimits( 0, 0., 1. );
  hist.Fit( "linfit", "QM+" );
  width = linfit -> GetParameter(0);
  error = 1. * ( linfit -> GetParError(0) );
  delete linfit;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

double removePointsforAlphaExtrapolation(bool isFE, double eta, int p) {
  // Cut to skip from fitting alpha point with not enough statistics due to pt_3 thr.
  // PT is meant to be the bin as in the pdf! aka p+1.
  // Values checked for UL
  double check = 0.;
  if (p>=1)  check = (eta>=eta_cut)? 0.15 : 0.20;
  if (p>=4)  check = (eta>=eta_cut)? 0.1 : 0.15;
  if (p>=7)  check = 0.1;
  if (p>=10) check = (eta>=eta_cut)? 0.1 : 0.05;

  // if(eta==2.043 && p==26) check = 0.1;

  return check;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

bool removePointsforFit(bool isFE, int m, int p) {
  bool check = false;
  //TODO FOR UL17 ONLY!
  // if ( isFE && p>=6 ) check = true;
  if ( isFE && m==3  && p>=6 ) check = true;
  if ( isFE && m==5  && p==6 ) check = true;
  if ( isFE && m==6  && p==8 ) check = true;
  if ( isFE && m==7  && p>=7 ) check = true;
  if ( isFE && m==8  && p==6 ) check = true;

  // if ( p<=1 ) check = true;
  // if ( !isFE && m==5  && p==8 ) check = true;
  // if ( !isFE && m==6  && p==8 ) check = true;
  // if ( !isFE && m==7  && p>=8 ) check = true;
  // if ( !isFE && m==8  && p>=5 ) check = true;
  // if ( !isFE && m==9  && p>=5 ) check = true;
  // if (  isFE && m==5  && p==2 ) check = true;
  // if (  isFE && m==6  && p==8 ) check = true;
  // if (  isFE && m==7  && p==7 ) check = true;
  // if (  isFE && m==8  && p==2 ) check = true;
  // if (  isFE && m==8  && p==7 ) check = true;
  // if (  isFE && m==10 && p==6 ) check = true;
  // if (  isFE && m==10 && p>=8 ) check = true;
  // if (  isFE && m==11 && p>=7 ) check = true;
  // if (  isFE && m==12 && p>=7 ) check = true;
  return check;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================


bool removePointsforFit2(bool isFE, int m, int p) {
  bool check = false;
  if ( p<=1 ) check = true;
  if ( !isFE && m==5  && p==8 ) check = true;
  if ( !isFE && m==6  && p==8 ) check = true;
  if ( !isFE && m==7  && p>=8 ) check = true;
  if ( !isFE && m==8  && p>=5 ) check = true;
  if ( !isFE && m==9  && p>=5 ) check = true;
  if (  isFE && m==5  && p==2 ) check = true;
  if (  isFE && m==6  && p==8 ) check = true;
  if (  isFE && m==7  && p==7 ) check = true;
  if (  isFE && m==8  && p==2 ) check = true;
  if (  isFE && m==8  && p==7 ) check = true;
  if (  isFE && m==11 && p==6 ) check = true;
  if (  isFE && m==11 && p>=8 ) check = true;
  if (  isFE && m==12 && p>=7 ) check = true;
  if (  isFE && m==13 && p>=5 ) check = true;
  return check;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TH2Poly* fill_2Dhist( TString name1, std::vector< std::vector< double > > SFs, std::vector< std::vector< double > > SFsError, std::vector<double> Pt_bins_Central, std::vector<double> Pt_bins_HF,std::vector<double> eta_bins, double eta_cut) {

  TH2Poly* h2 = new TH2Poly();
  for( unsigned int m = 0; m <  SFs.size(); m++ ) {
    //    for( unsigned int p = 0; p <  SFs.at(m).size(); p++ ) {
    for( unsigned int p = 0; p <  SFs.at(m).size()-1; p++ ) {
      bool eta_cut_bool = (fabs(eta_bins[m])+1e-3)>eta_cut;
      double pt_low = (eta_cut_bool?Pt_bins_HF:Pt_bins_Central)[p];
      double pt_up = (eta_cut_bool?Pt_bins_HF:Pt_bins_Central)[p+1];
      //      if(p==0) cout<<eta_cut_bool<<" eta = "<<eta_bins[m]<<" "<<eta_bins[m+1]<<" pt_low = "<<pt_low<<" pt_up = "<<pt_up<<" SFs.at(m).size() = "<<SFs.at(m).size()<<" eta_cut = "<<eta_cut<<endl;
      int realbinnumber = h2->AddBin(pt_low, eta_bins[m], pt_up, eta_bins[m+1]);
    }
  }

  h2->SetName(name1);
  h2->SetTitle(name1);
  h2->GetXaxis()->SetTitle("p_{T}^{ave}");
  h2->GetYaxis()->SetTitle("|#eta|");
  h2->GetZaxis()->SetRangeUser(0,3.0);
  h2->SetStats(kFALSE);
  if (name1.Contains("SF_")) h2 ->GetZaxis()->SetTitle("Scale factor");
  else                       h2 ->GetZaxis()->SetTitle("#sigma_{JER}");


  // for( unsigned int m = 0; m <  SFs.size(); m++ ) {
  //   for( unsigned int p = 0; p <  SFs.at(m).size(); p++ ) {
  for( unsigned int m = 0; m < SFs.size(); m++ ) {
    //      cout<<"SFs.at(m).size() = "<<SFs.at(m).size()<<endl;
    //      for( unsigned int p = 0; p < SFs.at(m).size(); p++ ) {
    for( unsigned int p = 0; p < SFs.at(m).size()-1; p++ ) {
      bool eta_cut_bool = (fabs(eta_bins[m])+1e-3)>eta_cut;
      double eta = fabs(eta_bins[m])+1e-3;
      double pT = (eta_cut_bool?Pt_bins_HF:Pt_bins_Central)[p]+1e-3;
      //	cout<<"m (eta) = "<<m<<" p (pT) = "<<p<<" pT = "<<pT<<" eta = "<<eta<<" Bin #"<<h2 -> FindBin(pT,eta)<<" SFs.at(m).at(p) = "<<SFs.at(m).at(p)<<" Error = "<<SFsError.at(m).at(p)<<endl;
      if ( ( !(TMath::IsNaN(SFs.at(m).at(p))) ) && SFs.at(m).at(p)!= 0. )      h2 -> SetBinContent(h2 -> FindBin(pT,eta), SFs.at(m).at(p) );
      if ( ( !(TMath::IsNaN(SFsError.at(m).at(p))) ) && SFs.at(m).at(p)!= 0. ) h2 -> SetBinError(h2 -> FindBin(pT,eta), SFsError.at(m).at(p) );
    }
  }
  //    cout<<"END pof 2D fill"<<endl;
  return h2;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TF1* GetFitsRatio(TF1* fit1, TF1* fit2, double ptmin, int color, TString func){
  bool isPOW = func.EqualTo(fNSxPC);
  bool isSIG = func.EqualTo(fNSCsign);
  TString f1 = "TMath::Sqrt( [0]*[0]/(x*x)+[1]*[1]/x+[2]*[2] )";
  TString f2 = "TMath::Sqrt( [3]*[3]/(x*x)+[4]*[4]/x+[5]*[5] )";
  if(isSIG) f1 = "TMath::Sqrt( TMath::Sign(1, [0])*[0]*[0]/(x*x)+[1]*[1]/x+[2]*[2] )";
  if(isPOW) f1 = "TMath::Sqrt([0]*abs([0])/(x*x)+[1]*[1]*pow(x,[6])+[2]*[2])";

  TF1* ratio = new TF1("ratio_", f1+"/"+f2, ptmin, 1200);
  ratio->SetParameter(0, fit1->GetParameter(0));
  ratio->SetParameter(1, fit1->GetParameter(1));
  ratio->SetParameter(2, fit1->GetParameter(2));
  ratio->SetParameter(3, fit2->GetParameter(0));
  ratio->SetParameter(4, fit2->GetParameter(1));
  ratio->SetParameter(5, fit2->GetParameter(2));
  if(isPOW) ratio->SetParameter(6, fit1->GetParameter(3));
  ratio->SetLineColor(color);

  return ratio;
}
