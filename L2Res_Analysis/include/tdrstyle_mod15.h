#pragma once

#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"
// Collect all the #includes here
#include "TROOT.h"
#include "TGraph.h"
#include "TLegend.h"

#include "TStyle.h"

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

#include "TFrame.h"

////////////////////////////////////
// Useful small macros (by Mikko) //
////////////////////////////////////

#include "TCanvas.h"
#include "TH1D.h"

#include <iostream>
#include <assert.h>

const bool kSquare = true;
const bool kRectangular = false;

using namespace std;

void tdrDraw(TH1* h, string opt,
	     int marker=kFullCircle, int mcolor = kBlack,
	     int lstyle=kSolid, int lcolor=-1,
	     int fstyle=1001, int fcolor=kYellow+1);

void tdrDraw(TGraph* g, string opt,
	     int marker=kFullCircle, int mcolor = kBlack,
	     int lstyle=kSolid, int lcolor=-1,
	     int fstyle=1001, int fcolor=kYellow+1);

TLegend tdrLeg(double x1, double y1, double x2, double y2);

//////////////////////////////////////////
// New CMS Style from 2014              //
// https://ghm.web.cern.ch/ghm/plots/   //
// Merged all macros into one
//////////////////////////////////////////

////////////////
// tdrstyle.C //
////////////////


// tdrGrid: Turns the grid lines on (true) or off (false)

void tdrGrid(bool gridOn);

// fixOverlay: Redraws the axis

void fixOverlay();

void setTDRStyle();

////////////////
// CMS_lumi.h //
////////////////

//
// Global variables
//

const TString cmsText = "CMS";;
const float cmsTextFont = 61;;  // default is helvetic-bold

const bool writeExtraText = true;//false;
const TString extraText   = "Simulation"; //Simulation, Preliminary
const TString extraText2  = ""; // For Simulation Preliminary on two lines
const float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
const float lumiTextSize     = 0.6;
const float lumiTextOffset   = 0.2;
const float cmsTextSize      = 0.75;
const float cmsTextOffset    = 0.1;  // only used in outOfFrame version

const float relPosX    = 0.045;
const float relPosY    = 0.035;
const float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
const float extraOverCmsTextSize  = 0.76;

const TString lumi_8TeV  = "19.7 fb^{-1}";
const TString lumi_7TeV  = "5.1 fb^{-1}";

const bool drawLogo      = false;


void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10, TString lumi_13TeV = "27.2 fb^{-1}" );



///////////////
// myMacro.C //
///////////////

// Give the macro an empty histogram for h->Draw("AXIS");
// Create h after calling setTDRStyle to get all the settings right
void tdrCanvas(TCanvas *&canv_in, const char* canvName, TH1D *&h,
		   int iPeriod = 2, int iPos = 11,
		   bool square = kRectangular, TString lumi_13TeV = "27.2 fb^{-1}");



// Give the macro empty histograms for h->Draw("AXIS");
// Create h after calling setTDRStyle to get all the settings right
// Created by: Mikko Voutilainen (HIP)
TCanvas* tdrDiCanvas(const char* canvName, TH1D *hup, TH1D *hdw,
		     int iPeriod = 2, int iPos = 11, TString lumi_13TeV = "27.2 fb^{-1}");
