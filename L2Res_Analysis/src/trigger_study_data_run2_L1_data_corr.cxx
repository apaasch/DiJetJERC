#include "TLatex.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TFile.h"
#include <iostream>
#include <TPad.h>
#include "TH1F.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TFitResultPtr.h"
#include <cmath>
#include <iomanip>
#include "TMinuit.h"
#include <algorithm>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "../include/parameters.h"

using namespace std;
double n_denom;

class TriggerThresholdDetermination {

    private:
        TFile * f_dijet;
        bool firstHist;
        TCanvas *c1, *c11, *c2;
        TGraphErrors *graph;

    public:
    TriggerThresholdDetermination(TString name = "");
    ~TriggerThresholdDetermination();

    void SetFile(TString fileName);

    void BuildRatio(TH1F* nominator, TH1F* denominator, TString nomin, TString denomin);

    void GetOfflineThreshold(double triggerThreshold);

};

TriggerThresholdDetermination::TriggerThresholdDetermination(TString name)
{
    firstHist = true;
    gStyle -> SetFrameBorderMode(0);

    gStyle -> SetNdivisions(10);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(1);
    gStyle -> SetOptTitle(1);
    gStyle -> SetStatFont(42);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetPadColor(0);
    gStyle -> SetTitleFont(62,"xy");
    gStyle -> SetLabelFont(62,"xy");
    gStyle -> SetTitleFontSize(0.06);
    gStyle -> SetTitleSize(0.08,"xy");
    gStyle -> SetLabelSize(0.06,"xy");
    gStyle -> SetHistFillStyle(3001);
    gStyle -> SetHistFillColor(0);
    gStyle -> SetHistLineStyle(1);
    gStyle -> SetHistLineWidth(2);
    gStyle -> SetHistLineColor(2);
    gStyle -> SetOptStat(1110);
    gStyle -> SetOptStat(kFALSE);
    gStyle -> SetOptFit(0111);
    gStyle -> SetStatH(0.1);
    gStyle -> SetErrorX(0.001);
    gStyle -> SetLabelOffset(0.02,"xy");
    gStyle -> SetTitleOffset(1.0,"y");
    gStyle -> SetTitleOffset(0.99,"x");

    c1 = new TCanvas(name+"c1","c1",0,0,800,600);
    gStyle->SetOptStat(0);
    c1->Divide(1,1,0,0);
    c1->SetFrameFillColor(0);

    c11 = new TCanvas(name+"c11","c11",0,0,800,600);
    gStyle->SetOptStat(0);
    c11->Divide(1,1,0,0);
    c11->SetFrameFillColor(0);


    c2 = new TCanvas(name+"c2","c2",0,0,800,600);
    c2->Divide(1,1);
    c2->SetFrameFillColor(0);

    graph = new TGraphErrors();
    graph->SetMarkerStyle(20);
}

TriggerThresholdDetermination::~TriggerThresholdDetermination()
{
}

void TriggerThresholdDetermination::SetFile(TString fileName){
    f_dijet = new TFile (fileName, "read");
}

Double_t SmoothFit_(Double_t *x, Double_t *par) {
   // cout << "n_denom = " << n_denom << endl;
    if (x[0] < n_denom) {
        TF1::RejectPoint();
        //  return 0;
    }
    Double_t p0 = par[0];
    Double_t p1 = par[1];
    Double_t N  = par[2];

  //Double_t fitval = N * (0.5 * (TMath::Erf(p0 * (x[0] - p1)) +1));
    Double_t fitval = 0.5 * N * (1. + TMath::Erf((x[0] - p0)/(pow(2,0.7)*p1)));

    return fitval;
}

void TriggerThresholdDetermination::BuildRatio(TH1F* nominator, TH1F* denominator, TString nomin, TString denomin)
{
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);

    TH1F * nom   = nominator;
    TH1F * denom = denominator;

   // TH1F * eff_norm = new TH1F("eff_norm","eff_norm", 500, 0, 500);

    TH1F * eff = (TH1F*) nom -> Clone("eff");
     eff -> Divide(denom);

    double offlineThreshold = eff->GetBinLowEdge(eff->FindFirstBinAbove(0.99,1));

    TString sThreshold = nomin;
    int startIdx = sThreshold.Index("DiPFJetAve", 10, 0, TString::kExact);
    sThreshold  = sThreshold(startIdx+10,4);
    double triggerThreshold = sThreshold.Atof();
    std::cout << triggerThreshold << " -> " << offlineThreshold << std::endl;

    int oldSize = graph->GetN();
    cout<<"Old  Size: "<<oldSize<<endl;

    TString sThreshold_fit = denomin;
    int startIdx_fit = sThreshold_fit.Index("DiPFJetAve", 10, 0, TString::kExact);
    sThreshold_fit  = sThreshold_fit(startIdx_fit+10,4);
    double triggerThreshold_fit = sThreshold_fit.Atof();
    n_denom = triggerThreshold_fit;

    TString sThreshold_tf1 = nomin;
    //TF1 * fEff_fi = sThreshold_tf1.Atof();

    cout<<"sThreshold: "<<sThreshold_tf1<<endl;
    cout<<"N Denom: "<<n_denom<<endl;
    TF1 * fEff_fit;
    fEff_fit = new TF1(sThreshold_tf1, SmoothFit_, n_denom, 250., 3);

    TF1 * fEff_fit_max = new TF1(sThreshold_tf1, SmoothFit_, n_denom, 250., 3);
    TF1 * fEff_fit_min = new TF1(sThreshold_tf1, SmoothFit_, n_denom, 250., 3);

   //Double_t Parameters[3] = {214.1630, 11.26654, 4.912097};
///    Double_t Parameters[3] = {84.1630, 7.26654, 27.912097};
   // Double_t Parameters[3] = {426.247, 11.26654, 0.912097};

    ////Double_t Parameters[3] = {426.247, 16.6424, 2.54782};
  
 if(firstHist) {
   //  Double_t Parameters[3] = {63, 2, 18};
  Double_t Parameters[3] = {80, 10, 1.3};
    fEff_fit -> SetParameters(Parameters);
    fEff_fit -> SetLineColor(4);
    fEff_fit -> SetParNames ("p0","p1","N");

 }
 else{
    cout<<"Trigger Threshold: "<<triggerThreshold<<endl;
   Double_t Parameters[3] = {triggerThreshold, 10, 3.7};
    fEff_fit -> SetParameters(Parameters);
    fEff_fit -> SetLineColor(4);
    fEff_fit -> SetParNames ("p0","p1","N");
 }



    c1->cd();
    eff -> SetTitle(0);
    eff->GetYaxis()->SetTitle("Efficiency");
    eff->GetXaxis()->SetTitle("p_{T}^{ave} (GeV)");
    eff->GetXaxis()->SetTitleOffset(1.25);
    eff->GetYaxis()->SetTitleOffset(1.25);
    eff->SetMarkerStyle(20);
    if(oldSize < 5)    eff->SetMarkerColor(kRed+oldSize);
    else    eff->SetMarkerColor(kBlue-oldSize+1);

    TLegend *leg = new TLegend(0.15,0.48+(oldSize*0.045),0.3,0.53+(oldSize*0.045));
    leg -> SetBorderSize(0);
    leg -> AddEntry(eff,nomin,"p");
    leg -> SetTextSize(0.035);
    leg -> SetFillColor(0);


    if(firstHist) {
        eff -> Draw("E");
        eff -> Fit(fEff_fit,"RLI","",0.5*triggerThreshold,3*triggerThreshold);

    }
    else {
        eff -> Draw("same E");
        eff -> Fit(fEff_fit,"RLI","same",0.5*triggerThreshold,1.5*triggerThreshold);
    }
    leg -> Draw("same");

    double Norm = fEff_fit -> GetParameter(2);
    cout << "Norm = "<< Norm<<endl;

    TH1F * eff_n = (TH1F*) eff -> Clone("eff_n");
    eff_n -> Scale(1/Norm);

 if(firstHist) {
  Double_t Parameters2[3] = {63, 2, 1};
    fEff_fit -> SetParameters(Parameters2);
    fEff_fit -> SetLineColor(4);
    fEff_fit -> SetParNames ("p0","p1","N");

 }
 else{
    cout<<"Trigger Threshold: "<<triggerThreshold<<endl;
   Double_t Parameters2[3] = {triggerThreshold, 10, 1};
    fEff_fit -> SetParameters(Parameters2);
    fEff_fit -> SetLineColor(4);
    fEff_fit -> SetParNames ("p0","p1","N");
 }


    c11->cd();
    if(firstHist) {
        eff_n->Draw("E");
    eff_n -> Fit(fEff_fit,"RLI","",0.5*triggerThreshold,3*triggerThreshold);        firstHist = false;
    } else {
        eff_n->Draw("same E");
   eff_n -> Fit(fEff_fit,"RLI","same",0.5*triggerThreshold,1.5*triggerThreshold);
    }
 

    fEff_fit_min -> SetParameter(0,(fEff_fit -> GetParameter(0) - fEff_fit -> GetParError(0)));
    fEff_fit_min -> SetParameter(1,(std::max(0.,fEff_fit -> GetParameter(1) - fEff_fit -> GetParError(1))));
    fEff_fit_min -> SetParameter(2,(fEff_fit -> GetParameter(2)));
    //   fEff_fit_min -> Draw("same");
    fEff_fit_max -> SetParameter(0,(fEff_fit -> GetParameter(0) + fEff_fit -> GetParError(0)));
    fEff_fit_max -> SetParameter(1,(fEff_fit -> GetParameter(1) + fEff_fit -> GetParError(1)));
    fEff_fit_max -> SetParameter(2,(fEff_fit -> GetParameter(2)));
    //    fEff_fit_max -> Draw("same");

      leg -> Draw("same");


    double offlineThreshold_n_h = eff_n->GetBinLowEdge(eff_n->FindFirstBinAbove(0.985,1));
    double offlineThreshold_n = fEff_fit ->GetX(0.99, 0, 600);
    double offlineThreshold_n_max = fEff_fit_max ->GetX(0.99, 0, 600);
    double offlineThreshold_n_min = fEff_fit_min ->GetX(0.99, 0, 600);

    cout << "min "<< offlineThreshold_n_min <<" nom "<< offlineThreshold_n << " max "<< offlineThreshold_n_max<<endl;

    int nBins = eff_n -> GetNbinsX(); // # of bins from histo
    int bin_n = eff_n-> GetBin(offlineThreshold_n);

    double offlineThreshold_n_err = offlineThreshold_n_max - offlineThreshold_n;
    cout << "bin num "<< bin_n <<""<<endl;



    std::cout <<"NORMALIZED: "<< triggerThreshold << " -> " << offlineThreshold_n_h << std::endl;
    cout <<"Thresholds from the fit for: "<< triggerThreshold << " -> " <<offlineThreshold_n<<" +/- "<< offlineThreshold_n_err<<endl;

    graph->SetPoint(oldSize, triggerThreshold, offlineThreshold_n);
    graph->SetPointError(oldSize, 0, offlineThreshold_n_err);

    c2->cd();
    graph->GetXaxis()->SetLabelOffset(0.01);
    graph->GetYaxis()->SetLabelOffset(0.01);
    graph->GetXaxis()->SetLabelSize(0.035);
    graph->GetYaxis()->SetLabelSize(0.035);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitle("99% efficiency threshold (GeV)");
    graph->GetXaxis()->SetTitle("Nominal trigger threshold (GeV)");
    graph->Draw("AP");

}

void TriggerThresholdDetermination::GetOfflineThreshold(double triggerThreshold)
{
    TF1* pol1 = new TF1("pol1","[0]+[1]*x",0,600);
    graph->Fit(pol1);
    std::cout << "extrapolated: " << 40 << " -> " << pol1->Eval(40) << std::endl;
    std::cout << "extrapolated: " << 60 << " -> " << pol1->Eval(60) << std::endl;
    std::cout << "extrapolated: " << 80 << " -> " << pol1->Eval(80) << std::endl;

}

void trigger_study_data_run2_L1_data_corr()
{
  cout << "--------------- Starting Derive_Thresholds_alternativeWay() ---------------" << endl << endl;

  TH1F *hdata_pt_ave[n_triggerDi];
  TH1F *hdata_pt_ave_next[n_triggerDi];

  TH1F *hdata_pt_ave_yield[n_triggerDi];
  TH1F *hdata_pt_ave_HF_yield[n_triggerDi_HF];

  TH1F *hdata_pt_ave_evt[n_triggerDi];
  TH1F *hdata_pt_ave_HF_evt[n_triggerDi_HF];

  TH1F *hdata_pt_ave_eta[n_triggerDi];
  TH1F *hdata_pt_ave_HF_eta[n_triggerDi_HF];

  TH1F *hdata_pt_ave_HF[n_triggerDi_HF];
  TH1F *hdata_pt_ave_HF_next[n_triggerDi_HF];

  int triggerThres[n_triggerDi+1];
  int triggerThres_HF[n_triggerDi_HF+1];

  for(int i = 0; i<n_triggerDi+1;i++){
    triggerThres[i] = pt_bins[i];
  }
  for(int i = 0; i<n_triggerDi_HF+1;i++){
    triggerThres_HF[i] = pt_bins_HF[i];
  }
  
  TString trigger_name[n_triggerDi];
  for(int i = 0; i<n_triggerDi;i++){
    trigger_name[i]="HLT_DiPFJetAve";
    trigger_name[i]+=to_string(triggerValDi[i]);
  }
  
  TString trigger_name_HF[n_triggerDi_HF] ;  
  for(int i = 0; i<n_triggerDi_HF;i++){
    trigger_name_HF[i]="HLT_DiPFJetAve";
    trigger_name_HF[i]+=to_string(triggerValDi_HF[i]);
    trigger_name_HF[i]+="_HFJEC";  
  }  

  for(int j=0; j<n_triggerDi; j++){
    TString name = "pt_ave_trg"+to_string(triggerValDi[j]);
    TString name2 = "pt_ave_next_trg"+to_string(triggerValDi[j]);
    TString name3 = "pt_ave_trig_yield"+to_string(triggerValDi[j]);
    TString name4 = "pt_ave_trig_evt"+to_string(triggerValDi[j]);
    TString name5 = "eta_trig_evt"+to_string(triggerValDi[j]);

    hdata_pt_ave[j]= new TH1F(name,"",nResponseBins,0,700);
    hdata_pt_ave[j]->Sumw2();
    hdata_pt_ave_yield[j]=new TH1F(name3,"",200,0,1000);
    hdata_pt_ave_yield[j]->Sumw2();
    hdata_pt_ave_evt[j]=new TH1F(name4,"",200,0,1000);
    hdata_pt_ave_evt[j]->Sumw2();
    hdata_pt_ave_next[j]= new TH1F(name2,"",nResponseBins,0,700);
    hdata_pt_ave_next[j]->Sumw2();
    hdata_pt_ave_eta[j] = new TH1F(name5,"",100,-5.191,5.191);
    hdata_pt_ave_eta[j] ->Sumw2();
  }
   
 for(int j=0; j<n_triggerDi_HF; j++){
    TString name = "pt_ave_HF_trg"+to_string(triggerValDi_HF[j]);
    TString name2 = "pt_ave_HF_next_trg"+to_string(triggerValDi_HF[j]);
    TString name3 = "pt_ave_HF_trg_yield"+to_string(triggerValDi_HF[j]);
    TString name4 = "pt_ave_trig_HF_evt"+to_string(triggerValDi_HF[j]);
    TString name5 = "eta_trig_HF_evt"+to_string(triggerValDi_HF[j]);

    hdata_pt_ave_HF[j]= new TH1F(name,"",60,0,450);
    hdata_pt_ave_HF[j]->Sumw2();
    hdata_pt_ave_HF_next[j]= new TH1F(name2,"",60,0,450);
    hdata_pt_ave_HF_next[j]->Sumw2();
    hdata_pt_ave_HF_yield[j]= new TH1F(name3,"",200,0,1000);
    hdata_pt_ave_HF_yield[j]->Sumw2();
    hdata_pt_ave_HF_evt[j]= new TH1F(name4,"",200,0,1000);
    hdata_pt_ave_HF_evt[j]->Sumw2();
    hdata_pt_ave_HF_eta[j]= new TH1F(name5,"",100,-5.191,5.191);
    hdata_pt_ave_HF_eta[j]->Sumw2();
  }
  
 
     TFile* _DATAFile = new TFile("/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_DeriveThresholds_woHOTVR_woXCONE/RunD_17Nov17_2017.root","READ");
 
  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", _DATAFile);
  TTreeReaderValue<int> trg40(myReader_DATA, "triggerDi40");
  TTreeReaderValue<int> trg60(myReader_DATA, "triggerDi60");
  TTreeReaderValue<int> trg80(myReader_DATA, "triggerDi80");
  TTreeReaderValue<int> trg140(myReader_DATA, "triggerDi140");
  TTreeReaderValue<int> trg200(myReader_DATA, "triggerDi200");
  TTreeReaderValue<int> trg260(myReader_DATA, "triggerDi260");
  TTreeReaderValue<int> trg320(myReader_DATA, "triggerDi320");
  TTreeReaderValue<int> trg400(myReader_DATA, "triggerDi400");
  TTreeReaderValue<int> trg500(myReader_DATA, "triggerDi500");
 
  
  TTreeReaderValue<int> trg60_HF(myReader_DATA, "trigger60_HF");
  TTreeReaderValue<int> trg80_HF(myReader_DATA, "trigger80_HF");
  TTreeReaderValue<int> trg100_HF(myReader_DATA, "trigger100_HF");
  TTreeReaderValue<int> trg160_HF(myReader_DATA, "trigger160_HF");
  TTreeReaderValue<int> trg220_HF(myReader_DATA, "trigger220_HF");
  TTreeReaderValue<int> trg300_HF(myReader_DATA, "trigger300_HF");
  

  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Int_t>   n_muon_data(myReader_DATA, "Nmuon");


  int myCount = 0;
  int myCount_notX = 0;
  bool allExclusive = true;
  while (myReader_DATA.Next()) {
    // if(*n_muon_data > 0) continue;
    
    int trg_arr[n_triggerDi] = {*trg40,*trg60,*trg80,*trg140,*trg200,*trg260,*trg320,*trg400,*trg500};
    int trg_arr_HF[n_triggerDi_HF] = {*trg60_HF,*trg80_HF,*trg100_HF,*trg160_HF,*trg220_HF,*trg300_HF};
    
    bool exclusive = true;
    bool exclusive_HF = true;
    exclusive = (*trg40)^(*trg60)^(*trg80)^(*trg140)^(*trg200)^(*trg260)^(*trg320)^(*trg400)^(*trg500);
    exclusive_HF = (*trg60_HF)^(*trg80_HF)^(*trg100_HF)^(*trg160_HF)^(*trg220_HF)^(*trg300_HF);
    
    if(fabs(*probejet_eta_data)<eta_cut){
      for(int j=0; j<n_triggerDi; j++){
	if((trg_arr[j])){
	  hdata_pt_ave[j]->Fill(*pt_ave_data);
	  if(*pt_ave_data>triggerThres[j]){
	    hdata_pt_ave_yield[j]->Fill(*pt_ave_data, *weight_data);
	    if(*pt_ave_data<triggerThres[j+1]){
	      hdata_pt_ave_evt[j]->Fill(*pt_ave_data, *weight_data);
	      hdata_pt_ave_eta[j]->Fill(*probejet_eta_data, *weight_data);
	    }
	  }
	}
	
	if((trg_arr[j+1])){
	  hdata_pt_ave_next[j]->Fill(*pt_ave_data);
	} 
      }
    }
    else{
      for(int j=0; j<n_triggerDi_HF; j++){
	if((trg_arr_HF[j])){
	  hdata_pt_ave_HF[j]->Fill(*pt_ave_data);
	  if(*pt_ave_data>triggerThres_HF[j]){
	    hdata_pt_ave_HF_yield[j]->Fill(*pt_ave_data, *weight_data);
	    if(*pt_ave_data<triggerThres_HF[j+1]){
	      hdata_pt_ave_HF_evt[j]->Fill(*pt_ave_data, *weight_data);
	      hdata_pt_ave_HF_eta[j]->Fill(*probejet_eta_data, *weight_data);
	    }
	  }
	  
	}
	if((trg_arr_HF[j+1])){
	  hdata_pt_ave_HF_next[j]->Fill(*pt_ave_data);
	}
	
      }
    }
    myCount++;
    if(!exclusive){
      myCount_notX++;
    }
  }
  
  
    TriggerThresholdDetermination trig;
    
    trig.BuildRatio(hdata_pt_ave_next[0],    hdata_pt_ave[0], "HLT_DiPFJetAve60", "HLT_DiPFJetAve40");
    
    trig.BuildRatio(hdata_pt_ave_next[1],    hdata_pt_ave[1], "HLT_DiPFJetAve80", "HLT_DiPFJetAve60");
    trig.BuildRatio(hdata_pt_ave_next[2],    hdata_pt_ave[2], "HLT_DiPFJetAve140", "HLT_DiPFJetAve80");
    trig.BuildRatio(hdata_pt_ave_next[3],    hdata_pt_ave[3], "HLT_DiPFJetAve200", "HLT_DiPFJetAve140");
    
    trig.BuildRatio(hdata_pt_ave_next[4],    hdata_pt_ave[4], "HLT_DiPFJetAve260", "HLT_DiPFJetAve200");
    
    trig.BuildRatio(hdata_pt_ave_next[5],    hdata_pt_ave[5], "HLT_DiPFJetAve320", "HLT_DiPFJetAve260");
    trig.BuildRatio(hdata_pt_ave_next[6],    hdata_pt_ave[6], "HLT_DiPFJetAve400", "HLT_DiPFJetAve320");
    trig.BuildRatio(hdata_pt_ave_next[7],    hdata_pt_ave[7], "HLT_DiPFJetAve500", "HLT_DiPFJetAve400");
  
    
    TriggerThresholdDetermination trigHF("HF"); 
    trigHF.BuildRatio(hdata_pt_ave_HF_next[0],    hdata_pt_ave_HF[0], "HLT_DiPFJetAve80_HFJEC", "HLT_DiPFJetAve60_HFJEC");    
    trigHF.BuildRatio(hdata_pt_ave_HF_next[1],    hdata_pt_ave_HF[1], "HLT_DiPFJetAve100_HFJEC", "HLT_DiPFJetAve80_HFJEC");
    trigHF.BuildRatio(hdata_pt_ave_HF_next[2],    hdata_pt_ave_HF[2], "HLT_DiPFJetAve160_HFJEC", "HLT_DiPFJetAve100_HFJEC");
    trigHF.BuildRatio(hdata_pt_ave_HF_next[3],    hdata_pt_ave_HF[3], "HLT_DiPFJetAve220_HFJEC", "HLT_DiPFJetAve160_HFJEC");
    trigHF.BuildRatio(hdata_pt_ave_HF_next[4],    hdata_pt_ave_HF[4], "HLT_DiPFJetAve300_HFJEC", "HLT_DiPFJetAve220_HFJEC");
    
    trig.GetOfflineThreshold(100);
    trigHF.GetOfflineThreshold(100);
    
    TLegend *leg1 = new TLegend(0.6,0.48,0.8,0.85);
    leg1 -> SetBorderSize(0);
    leg1 -> SetTextSize(0.035);
    leg1 -> SetFillColor(0);   

    TCanvas* c3 = new TCanvas("trig_yield","",0,0,800,600);
    gStyle->SetOptStat(0);
    c3->SetFrameFillColor(0);
    double maximum = 0;
    for(int i=0;i<n_triggerDi;i++){
      maximum = max(maximum,hdata_pt_ave_yield[i]->GetMaximum());
    }
    for(int i=0;i<n_triggerDi;i++){
      if(i==0){
	hdata_pt_ave_yield[i]->GetYaxis()->SetTitle("Events");
	hdata_pt_ave_yield[i]->GetYaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_yield[i]->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
	hdata_pt_ave_yield[i]->GetXaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_yield[i]->SetMaximum(maximum*1.2);
      }
      hdata_pt_ave_yield[i]->SetMarkerStyle(20);
      hdata_pt_ave_yield[i]->SetMarkerColor(kRed-i*5);
      hdata_pt_ave_yield[i]->Draw("P SAME");

      TLine * line = new TLine(triggerThres[i], 0,triggerThres[i], maximum*1.2);
      line->SetLineStyle(2);
      line->Draw("SAME");

      leg1 -> AddEntry(hdata_pt_ave_yield[i],trigger_name[i],"p");
    }
    leg1->Draw("SAME");
    c3->Print("/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_DeriveThresholds_woHOTVR_woXCONE/Trigger_yield.pdf");


    TLegend *leg2 = new TLegend(0.5,0.48,0.8,0.8);
    leg2 -> SetBorderSize(0);
    leg2 -> SetTextSize(0.035);
    leg2 -> SetFillColor(0); 

    TCanvas* c4 = new TCanvas("trig_yield_HF","",0,0,800,600);
    gStyle->SetOptStat(0);
    c4->SetFrameFillColor(0);  
    maximum = 0;
    for(int i=0;i<n_triggerDi_HF;i++){
      maximum = max(maximum,hdata_pt_ave_HF_yield[i]->GetMaximum());
    }
    for(int i=0;i<n_triggerDi_HF;i++){
      if(i==0){
	hdata_pt_ave_HF_yield[i]->GetYaxis()->SetTitle("Events");
	hdata_pt_ave_HF_yield[i]->GetYaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_HF_yield[i]->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
	hdata_pt_ave_HF_yield[i]->GetXaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_HF_yield[i]->SetMaximum(maximum*1.2);
      }
      hdata_pt_ave_HF_yield[i]->SetMarkerStyle(20);
      hdata_pt_ave_HF_yield[i]->SetMarkerColor(kBlue+i*3);
      hdata_pt_ave_HF_yield[i]->Draw("P SAME");

      TLine * line = new TLine(triggerThres_HF[i], 0,triggerThres_HF[i], maximum*1.2);
      line->SetLineStyle(2);
      line->Draw("SAME");

      leg2 -> AddEntry(hdata_pt_ave_HF_yield[i],trigger_name_HF[i],"p");
    }
    leg2->Draw("SAME");
    c4->Print("/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_DeriveThresholds_woHOTVR_woXCONE/Trigger_HF_yield.pdf");
    
    TLegend *leg3 = new TLegend(0.6,0.15,0.8,0.45);
    leg3 -> SetBorderSize(0);
    leg3 -> SetTextSize(0.035);
    leg3 -> SetFillColor(0); 

    TCanvas* c5 = new TCanvas("trig_evt","",0,0,800,600);
    c5->SetLogy();
    gStyle->SetOptStat(0);
    c5->SetFrameFillColor(0);
    
     maximum = 0;
    for(int i=0;i<n_triggerDi;i++){
      maximum = max(maximum,hdata_pt_ave_evt[i]->GetMaximum());
    }
    for(int i=0;i<n_triggerDi;i++){
      if(i==0){
	hdata_pt_ave_evt[i]->GetYaxis()->SetTitle("Events");
	hdata_pt_ave_evt[i]->GetYaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_evt[i]->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
	hdata_pt_ave_evt[i]->GetXaxis()->SetTitleOffset(1.1);
      	hdata_pt_ave_evt[i]->SetMaximum(maximum*1.2);
	hdata_pt_ave_evt[i]->SetMinimum(80);
      }
      hdata_pt_ave_evt[i]->SetMarkerStyle(20);
      hdata_pt_ave_evt[i]->SetMarkerColor(kRed-i*5);
      hdata_pt_ave_evt[i]->Draw("P SAME");

      TLine * line = new TLine(triggerThres[i], 0,triggerThres[i], maximum*1.2);
      line->SetLineStyle(2);
      line->Draw("SAME");
      leg3 -> AddEntry(hdata_pt_ave_evt[i],trigger_name[i],"p");
    }
    leg3->Draw("SAME");
    c5->Print("/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_DeriveThresholds_woHOTVR_woXCONE/Trigger_Events.pdf");


    TLegend *leg4 = new TLegend(0.5,0.5,0.8,0.8);
    leg4 -> SetBorderSize(0);
    leg4 -> SetTextSize(0.035);
    leg4 -> SetFillColor(0); 


    TCanvas* c6 = new TCanvas("trig_evt_HF","",0,0,800,600);
    c6->SetLogy();
    gStyle->SetOptStat(0);
    c6->SetFrameFillColor(0); 
     maximum = 0;
    for(int i=0;i<n_triggerDi_HF;i++){
      maximum = max(maximum,hdata_pt_ave_HF_evt[i]->GetMaximum());
    } 
    for(int i=0;i<n_triggerDi_HF;i++){
      if(i==0){
	hdata_pt_ave_HF_evt[i]->GetYaxis()->SetTitle("Events");
	hdata_pt_ave_HF_evt[i]->GetYaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_HF_evt[i]->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
	hdata_pt_ave_HF_evt[i]->GetXaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_HF_evt[i]->SetMaximum(maximum*1.2);
	hdata_pt_ave_HF_evt[i]->SetMinimum(0.8);
      }
      hdata_pt_ave_HF_evt[i]->SetMarkerStyle(20);
      hdata_pt_ave_HF_evt[i]->SetMarkerColor(kBlue+i*3);
      hdata_pt_ave_HF_evt[i]->Draw("P SAME");

      TLine * line = new TLine(triggerThres_HF[i], 0,triggerThres_HF[i], maximum*1.2);
      line->SetLineStyle(2);
      line->Draw("SAME");
      leg4 -> AddEntry(hdata_pt_ave_HF_evt[i],trigger_name_HF[i],"p");
    }
    leg4->Draw("SAME");
    c6->Print("/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_DeriveThresholds_woHOTVR_woXCONE/Trigger/Trigger_HF_Events.pdf");



    TCanvas* c7 = new TCanvas("trig_eta","",0,0,800,600);
    c7->SetLogy();
    gStyle->SetOptStat(0);
    c7->SetFrameFillColor(0);
    
    maximum = 0;
    for(int i=0;i<n_triggerDi;i++){
      maximum = max(maximum,hdata_pt_ave_eta[i]->GetMaximum());
    }
    for(int i=0;i<n_triggerDi;i++){
      if(i==0){
	hdata_pt_ave_eta[i]->GetYaxis()->SetTitle("Events");
	hdata_pt_ave_eta[i]->GetYaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_eta[i]->GetXaxis()->SetTitle("probejet #eta");
	hdata_pt_ave_eta[i]->GetXaxis()->SetTitleOffset(1.1);
      	hdata_pt_ave_eta[i]->SetMaximum(maximum*1.2);
	hdata_pt_ave_eta[i]->SetMinimum(0.8);
      }
      hdata_pt_ave_eta[i]->SetMarkerStyle(20);
      hdata_pt_ave_eta[i]->SetMarkerColor(kRed-i*5);
      hdata_pt_ave_eta[i]->Draw("P SAME");
 
      //leg3 -> AddEntry(hdata_pt_ave_eta[i],trigger_name[i],"p");
    }
    //  leg3->Draw("SAME");
    TLine * line1 = new TLine(-2.853, 0, -2.8, maximum*1.2);
      line1->SetLineStyle(2);
      line1->Draw("SAME");
      line1 = new TLine(2.853, 0, 2.8, maximum*1.2);
      line1->SetLineStyle(2);
      line1->Draw("SAME");
    c7->Print("/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_DeriveThresholds_woHOTVR_woXCONE/Trigger_eta.pdf");


    TCanvas* c8 = new TCanvas("trig_eta_HF","",0,0,800,600);
    c8->SetLogy();
    gStyle->SetOptStat(0);
    c8->SetFrameFillColor(0);
    
   maximum = 0;
    for(int i=0;i<n_triggerDi_HF;i++){
      maximum = max(maximum,hdata_pt_ave_HF_eta[i]->GetMaximum());
    }
    for(int i=0;i<n_triggerDi_HF;i++){
      if(i==0){
	hdata_pt_ave_HF_eta[i]->GetYaxis()->SetTitle("Events");
	hdata_pt_ave_HF_eta[i]->GetYaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_HF_eta[i]->GetXaxis()->SetTitle("probejet #eta");
	hdata_pt_ave_HF_eta[i]->GetXaxis()->SetTitleOffset(1.1);
	hdata_pt_ave_HF_eta[i]->SetMaximum(maximum*1.2);
	hdata_pt_ave_HF_eta[i]->SetMinimum(0.08);
      }
      hdata_pt_ave_HF_eta[i]->SetMarkerStyle(20);
      hdata_pt_ave_HF_eta[i]->SetMarkerColor(kBlue+i*3);
      hdata_pt_ave_HF_eta[i]->Draw("P SAME");

      //     leg4 -> AddEntry(hdata_pt_ave_HF_eta[i],trigger_name_HF[i],"p");
    }
    // leg4->Draw("SAME");
     line1 = new TLine(-2.853, 0,-2.8, maximum*1.2);
      line1->SetLineStyle(2);
      line1->Draw("SAME");
      line1 = new TLine(2.853, 0, 2.8, maximum*1.2);
      line1->SetLineStyle(2);
      line1->Draw("SAME");
      c8->Print("/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_DeriveThresholds_woHOTVR_woXCONE/Trigger_HF_eta.pdf");

}
