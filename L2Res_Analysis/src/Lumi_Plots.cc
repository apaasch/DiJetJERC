#include "../include/parameters.h"
#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"

#include <iostream>
#include <fstream>
#include <string> 

#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TFile.h>
#include <TProfile.h>
#include <TString.h>

using namespace std;

void CorrectionObject::Lumi_Plots(){
  cout << "--------------- StartingLumi_Plots () ---------------" << endl << endl;
  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/lumi/");
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(000);

  int n_trigger = n_triggerSi;
  vector<int> triggerVal(triggerValSi, triggerValSi + sizeof(triggerValSi) / sizeof(triggerValSi[0]));
  
  // fill the histos for pt average in bins of eta
  TH1D* lumi_data[n_trigger+5];

  int countPt = 0;
  TString namehist = "lumi_";

  TString dirNames[n_trigger+5];
  for(int i = 0; i<n_trigger; i++){
    dirNames[i] = TString("Lumi_Trig") + to_string(triggerVal[i]); 
  }
  dirNames[n_trigger] = "Lumi_noCuts";
  dirNames[n_trigger+1] = "Lumi_Selection";
  dirNames[n_trigger+2] = "Lumi_diJet";
  dirNames[n_trigger+3] = "Lumi_JetMatching";
  dirNames[n_trigger+4] = "Lumi_Final";

  TString var1 = "luminosity";
  
    
  for(int i=0; i<n_trigger+5; i++){
    TString dirName = dirNames[i];
  
    TH1D* hist;
    
    cout << "looking for "<<dirName << endl;
    
    lumi_data[i] =(TH1D*) CorrectionObject::_DATAFile->Get(dirName+"/"+var1)->Clone();
    // lumi_data[i]->Rebin(6);
    // lumi_data[i]->SaveAs(CorrectionObject::_outpath+"plots/lumi/"+dirName+".root");
    }

  cout << "fit the thresholds" << endl << endl;
    
  cout << "Draw plots to "<<CorrectionObject::_outpath+"plots/lumi" << endl;

  for(int i=0; i<n_trigger+5; i++){
      TCanvas* c1 = new TCanvas("c1");
      lumi_data[i]->Draw();
      c1->Print(CorrectionObject::_outpath+"plots/lumi/"+dirNames[i]+".eps","eps");

  }

  
}
