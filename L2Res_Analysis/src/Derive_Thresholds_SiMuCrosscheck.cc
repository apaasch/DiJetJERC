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

void CorrectionObject::Derive_Thresholds_SiMuCrosscheck(TString siMu_dirName_ = "HLT_IsoMu20"){
  cout << "--------------- Starting Derive_Thresholds_SiMuCrosscheck() ---------------" << endl << endl;
  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/thresholds/");
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(000);

  int n_trigger = n_triggerSi;
  vector<int> triggerVal(triggerValSi, triggerValSi + sizeof(triggerValSi) / sizeof(triggerValSi[0]));
  
  // fill the histos for pt average in bins of eta
  TH1D* ptave_data[n_trigger];
  TH1D* ptave_data_eff[n_trigger-1];  
  int countPt = 0;
  TString namehist = "ptave_";

  TString siMu_dirName = siMu_dirName_ ;

  TString var1 = "pt_ave";
  
  TH1D* ptave_data_SiMu = (TH1D*) CorrectionObject::_DATAFile->Get(siMu_dirName+"/"+var1)->Clone();
  ptave_data_SiMu->SetName(siMu_dirName);
  ptave_data_SiMu->Print();
  
  for(int i=0; i<n_trigger; i++){
    TString dirName = "HLT_PFJet";
    dirName += to_string(triggerVal[i]); 
  
    TH1D* hist;
    
    cout << "looking for "<<dirName << endl;
    
    countPt++;

    ptave_data_eff[i] =(TH1D*) CorrectionObject::_DATAFile->Get(dirName+"/"+var1)->Clone();
    ptave_data_eff[i]->SetName(dirName);
    ptave_data_eff[i]->Print();
    ptave_data_eff[i]->Divide(ptave_data_SiMu);
    ptave_data_eff[i]->Rebin(6);
    ptave_data_eff[i]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"SiMuCrosscheck_HLT_PFJet"+to_string(triggerVal[i])+".root");
    }

  cout << "fit the thresholds" << endl << endl;

  TF1 *func[n_trigger];
  double thresholds[n_trigger];
  double thresholds09[n_trigger];
  double thresholds_err[n_trigger];
  double thresholds09_err[n_trigger];
  double thresholds_errUp[n_trigger];
  double thresholds09_errUp[n_trigger];
  double thresholds_errDown[n_trigger];
  double thresholds09_errDown[n_trigger];
  double fitrange[n_trigger] = {50, 50, 50, 50, 40, 10, 0 , 0, 0, 0};
  for(int i=0; i<n_trigger; i++){
    TString fitname = "fit";
    fitname +=  to_string(triggerVal[i]); 
    func[i] = new TF1(fitname,SmoothFit,(triggerVal[i-1]?i>0:50)-fitrange[i],
		      triggerVal[i]+100+fitrange[i],3);
    func[i]->SetParameters(triggerVal[i], 40., 1.);
    func[i]->SetParNames("p0", "p1", "N");
    ptave_data_eff[i]->Fit(func[i],"R");

    thresholds[i] = func[i]->GetX(0.95*func[i]->GetParameter(2), triggerVal[i]-10, 800);
    thresholds09[i] = func[i]->GetX(0.9*func[i]->GetParameter(2), triggerVal[i]-10, 800);
 

        cout<< std::fixed << std::setprecision(2) << "Trigger value: "<< triggerVal[i] << "; 0.95 threshold: "<< thresholds[i]  << endl;;
    cout<< std::fixed << std::setprecision(2) << "Trigger value: "<< triggerVal[i] << "; 0.90 threshold: "<< thresholds09[i]  << endl;
    cout << "Old Threshold: " << pt_bins[i] <<endl <<endl;
  }
  
    
  cout << "Draw plots to "<<CorrectionObject::_outpath+"plots/thresholds" << endl;

  ofstream myfile;
  myfile.open (CorrectionObject::_outpath+"plots/thresholds/"+"/thresholds.txt",ios::trunc);
  myfile << "trigger  old  \t  0.90 \t \t  0.95  \n";
  for(int i=0; i<n_trigger; i++){
    myfile<< std::fixed << std::setprecision(2) << triggerVal[i]<< " \t " << pt_bins[i]<< " \t "<< thresholds09[i]<< " \t "<< thresholds[i] << " \n"; 
  }
  myfile.close();
  
  for(int i=0; i<n_trigger; i++){
      TCanvas* c1 = new TCanvas("c1");
      ptave_data_eff[i]->Draw();
      TLine *line = new TLine(thresholds[i],0,thresholds[i],
			      1.05*func[i]->GetParameter(2));
      line->SetLineColor(kBlack);
      line->Draw();
       TLine *line09 = new TLine(thresholds09[i],0,thresholds09[i],
				 1.05*func[i]->GetParameter(2));
      line09->SetLineColor(kBlack);
      line09->SetLineStyle(2);
      line09->Draw();     
      c1->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[i])+".eps","eps");

  }

  
}
