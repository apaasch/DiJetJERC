#pragma once

#include <iostream>
#include <utility>
#include <cmath>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

std::pair<double,double> Rmc_to_Rdata(std::pair<double,double> mc, std::pair<double,double> data);

std::pair<double,double> GetValueAndError(TH1D *hin);

double ErrorPropagation_AoverB(std::pair<double,double> Ap, std::pair<double,double> Bp);
double ErrorPropagation_AB(std::pair<double,double> Ap, std::pair<double,double> Bp);

TGraphErrors* BuildRatio(TGraphErrors* input, double ave, double err_ave);
TGraphErrors* CleanEmptyPoints(TGraphErrors* input);

TH1D* GetHist(TFile *rootfile, TString selection, TString varName, int nbins, double low, double up);
Double_t SmoothFit(Double_t *v, Double_t *par);
