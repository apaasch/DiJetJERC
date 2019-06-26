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

void CorrectionObject::OnOffResp_Plots(){
  cout << "--------------- StartingOnOffResp_Plots () ---------------" << endl << endl;
  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/onoffresp/");
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(000);

  TH1D* oor1 = (TH1D*) CorrectionObject::_DATAFile->Get("Final/pt_1_onoff_Resp")->Clone();
  TH1D* oor2 = (TH1D*) CorrectionObject::_DATAFile->Get("Final/pt_2_onoff_Resp")->Clone();
  TCanvas* c11 = new TCanvas("c11");
  oor1->Draw();
  c11->Print(CorrectionObject::_outpath+"plots/onoffresp/"+"final"+"_pt_resp_0.eps","eps");

  TCanvas* c12 = new TCanvas("c12");
  oor2->Draw();
  c12->Print(CorrectionObject::_outpath+"plots/onoffresp/"+"final"+"_pt_resp_1.eps","eps");
     
  int n_trigger = n_triggerSi;
  vector<int> triggerVal(triggerValSi, triggerValSi + sizeof(triggerValSi) / sizeof(triggerValSi[0]));
  
  // fill the histos for pt average in bins of eta
  TH1D* trg_data0[n_trigger];
  TH1D* trg_data1[n_trigger];
  
  int countPt = 0;
  TString namehist = "trg_";

  TString dirNames[n_trigger];
  for(int i = 0; i<n_trigger; i++){
    dirNames[i] = TString("HLT_PFJet") + to_string(triggerVal[i]); 
  }

  TString var1 = "pt_1_onoff_Resp";
  TString var2 = "pt_2_onoff_Resp";  
    
  for(int i=0; i<n_trigger; i++){
    TString dirName = dirNames[i];
  
    TH1D* hist;
    
    cout << "looking for "<<dirName << endl;
    
    trg_data0[i] =(TH1D*) CorrectionObject::_DATAFile->Get(dirName+"/"+var1)->Clone();
    trg_data1[i] =(TH1D*) CorrectionObject::_DATAFile->Get(dirName+"/"+var2)->Clone();
    // trg_data[i]->Rebin(6);
    // trg_data[i]->SaveAs(CorrectionObject::_outpath+"plots/trg/"+dirName+".root");
    }

  cout << "Draw plots to "<<CorrectionObject::_outpath+"plots/onoffresp" << endl;

  for(int i=0; i<n_trigger; i++){
      TCanvas* c1 = new TCanvas("c1");
      trg_data0[i]->Draw();
      c1->Print(CorrectionObject::_outpath+"plots/onoffresp/"+dirNames[i]+"_pt_resp_0.eps","eps");
      
      TCanvas* c2 = new TCanvas("c2");
      trg_data1[i]->Draw();
      c2->Print(CorrectionObject::_outpath+"plots/onoffresp/"+dirNames[i]+"_pt_resp_1.eps","eps");
  }

  
}
