// Read histograms filled in JECAnalysisPUjetsHists (by BaconJet/AnalysisModule module) 
// to study jets not matched to GEN jets (most probable PU jets)
// created on 10.04.2019 by A. Karavdina

#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include <TH1F.h>
#include <TH2D.h>
#include <TH3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <THStack.h>
#include <TLatex.h>
using namespace std;

void doPlot_2D_PIid_PF(TH2D* Jet_PU_id_chEM_BB1_MC,TH2D* Jet_PU_id_chEM_BB1_bad,TH2D* Jet_PU_id_chEM_BB1_DATA, TString PFlabel,TString ptbinLabel,TString txttag,TString SavePlots, int eta_bin){
  TLegend* leg =  new TLegend(0.17,0.67,0.52,0.89);
  leg->SetTextSize(0.045);
  TString eta_bin_txt; 
  TString eta_bin_label; 
  if(eta_bin==0){
    eta_bin_txt="PU ID (|#eta|<1.3)";
    eta_bin_label="BB1";
  }
  if(eta_bin==1){ 
    eta_bin_txt="PU ID (1.3<|#eta|<2.1)";
    eta_bin_label="BB2";
  }
  if(eta_bin==2){
    eta_bin_txt="PU ID (2.1<|#eta|<2.5)";
    eta_bin_label="EC1";
  }
  if(eta_bin==3){
    eta_bin_txt="PU ID (2.5<|#eta|<3.0)";
    eta_bin_label="EC2";
  }
  if(eta_bin==4){
    eta_bin_txt="PU ID (|#eta|>3.0)";
    eta_bin_label="HF";
  }


  TCanvas* cJet_PU_id_chEM_BB1_MC = new TCanvas();
  Jet_PU_id_chEM_BB1_MC->GetXaxis()->SetLabelSize(0.045);
  Jet_PU_id_chEM_BB1_MC->GetYaxis()->SetLabelSize(0.045);
  Jet_PU_id_chEM_BB1_MC->GetYaxis()->SetTitleSize(0.045);
  Jet_PU_id_chEM_BB1_MC->SetTitle("");

  Jet_PU_id_chEM_BB1_MC->SetMarkerColor(kRed);
  Jet_PU_id_chEM_BB1_MC->SetLineColor(kRed);
  Jet_PU_id_chEM_BB1_MC->SetMarkerSize(1.1);
  //  Jet_PU_id_chEM_BB1_MC->GetZaxis()->SetRangeUser(1e-3,1e3);
  Jet_PU_id_chEM_BB1_MC->GetXaxis()->SetTitle(eta_bin_txt);
  Jet_PU_id_chEM_BB1_MC->GetYaxis()->SetTitle(PFlabel);
  Jet_PU_id_chEM_BB1_MC->Draw("P");

  Jet_PU_id_chEM_BB1_bad->SetMarkerColor(kBlack);
  Jet_PU_id_chEM_BB1_bad->SetLineColor(kBlack);
  Jet_PU_id_chEM_BB1_bad->SetMarkerSize(1.1);
  Jet_PU_id_chEM_BB1_bad->Draw("SAME");
  // Jet_PU_id_chEM_BB1_DATA->SetMarkerColor(kGreen);
  // Jet_PU_id_chEM_BB1_DATA->SetLineColor(kGreen);
  // //  Jet_PU_id_chEM_BB1_DATA->SetMarkerStyle(20);
  // Jet_PU_id_chEM_BB1_DATA->SetMarkerSize(1.1);
  // Jet_PU_id_chEM_BB1_DATA->Draw("SAME");

  // //  cJet_PU_id_chEM_BB1_MC->SetLogy();
  // // cJet_PU_id_chEM_BB1_MC->SetLogz();
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045); 
  tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  leg->AddEntry(Jet_PU_id_chEM_BB1_MC,"MC, matched","pl");
  leg->AddEntry(Jet_PU_id_chEM_BB1_bad,"MC, not matched","pl");
  //  leg->AddEntry(Jet_PU_id_chEM_BB1_DATA,"DATA","pl");
  leg->Draw();
  cJet_PU_id_chEM_BB1_MC->Print(SavePlots + "_Jet_PU_id_"+PFlabel+"_"+eta_bin_label+".pdf","pdf");
}

void doPlot_PUid(TH1F* PU_id_BB1_MC, TH1F* PU_id_BB1_notmatchedMC, TH1F* PU_id_BB1_DATA, int eta_bin, bool is3rd,TString ptbinLabel,TString txttag,TString SavePlots, bool isParton){
  TLegend* leg =  new TLegend(0.17,0.67,0.52,0.89);
  leg->SetTextSize(0.045);
  TString legend_label = "1st and 2nd jets";
  if(is3rd) legend_label = "3rd jet";
  leg->SetHeader(legend_label);
  TString eta_bin_txt; 
  TString eta_bin_label; 
  if(eta_bin==0){
    eta_bin_txt="PU ID (|#eta|<1.3)";
    eta_bin_label="BB1";
  }
  if(eta_bin==1){ 
    eta_bin_txt="PU ID (1.3<|#eta|<2.1)";
    eta_bin_label="BB2";
  }
  if(eta_bin==2){
    eta_bin_txt="PU ID (2.1<|#eta|<2.5)";
    eta_bin_label="EC1";
  }
  if(eta_bin==3){
    eta_bin_txt="PU ID (2.5<|#eta|<3.0)";
    eta_bin_label="EC2";
  }
  if(eta_bin==4){
    eta_bin_txt="PU ID (|#eta|>3.0)";
    eta_bin_label="HF";
  }
  
   TCanvas* cPU_id_BB1_MC = new TCanvas();
   PU_id_BB1_MC->GetXaxis()->SetTitle(eta_bin_txt);
   PU_id_BB1_MC->GetXaxis()->SetLabelSize(0.045);
   PU_id_BB1_MC->GetYaxis()->SetLabelSize(0.045);
   PU_id_BB1_MC->GetYaxis()->SetTitleSize(0.045);
   PU_id_BB1_MC->SetMarkerColor(kRed);
   PU_id_BB1_MC->SetMarkerStyle(20);
   PU_id_BB1_DATA->SetMarkerColor(kGreen);
   PU_id_BB1_DATA->SetMarkerStyle(22);

   PU_id_BB1_MC->SetTitle("");
   PU_id_BB1_notmatchedMC->SetMarkerColor(kBlack);
   PU_id_BB1_notmatchedMC->SetMarkerStyle(21);
   PU_id_BB1_MC->Print();
   PU_id_BB1_notmatchedMC->Print();
   Double_t scale = 1/PU_id_BB1_MC->Integral();
   if(scale>0)
     PU_id_BB1_MC->Scale(scale);
   PU_id_BB1_MC->GetYaxis()->SetRangeUser(0.0,0.2);
   PU_id_BB1_MC->Draw();
   scale = 1/PU_id_BB1_DATA->Integral();
   //   PU_id_BB1_DATA->Print();
   if(scale>0){
     PU_id_BB1_DATA->Scale(scale);
     PU_id_BB1_DATA->Draw("same");
   }

   scale = 1/PU_id_BB1_notmatchedMC->Integral();
   if(scale>0){
     PU_id_BB1_notmatchedMC->Scale(scale);
     PU_id_BB1_notmatchedMC->Draw("same");
   }
   //   TString txttag=CorrectionObject::_generator_tag; 
   TLatex *tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.045);  
   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   if(isParton){
     leg->AddEntry(PU_id_BB1_MC,"MC, ptcl matched","p");
     leg->AddEntry(PU_id_BB1_notmatchedMC,"MC, ptcl not matched","p");
   }
   else{
     leg->AddEntry(PU_id_BB1_MC,"MC, matched","p");
     leg->AddEntry(PU_id_BB1_notmatchedMC,"MC, not matched","p");
   }
   leg->AddEntry(PU_id_BB1_DATA,"DATA","p");

   leg->Draw();
   if(isParton){
     if(is3rd){
       cPU_id_BB1_MC->Print(SavePlots + "_PU_id_parton_3rdJet_"+eta_bin_label+".pdf","pdf");
     }
     else
       cPU_id_BB1_MC->Print(SavePlots + "_PU_id_parton_"+eta_bin_label+".pdf","pdf");
     
   }
   else{
     if(is3rd){
       cPU_id_BB1_MC->Print(SavePlots + "_PU_id_3rdJet_"+eta_bin_label+".pdf","pdf");
     }
     else
       cPU_id_BB1_MC->Print(SavePlots + "_PU_id_"+eta_bin_label+".pdf","pdf");
   }
   delete cPU_id_BB1_MC;
}

void CorrectionObject::ControlPlotsPUjets(bool forEveryPtBin){
  cout << "--------------- Starting ControlPlotsPUjets() ---------------" << endl;
  bool forEveryPtBin_ = forEveryPtBin;
  if(forEveryPtBin_){
    cout<<"ControlPlotsPUjets are gonne be made for each pt bin seperatly"<<endl;
  }
  TString txttag=CorrectionObject::_generator_tag; 
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"x");  
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"x");  
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleYOffset(0.9);
  gStyle->SetTitleXOffset(0.8);
  gStyle->SetPalette(55);
  gROOT->ForceStyle();



  // TString dirNames[30] = {"goodRECO_all", 
// 			  "goodRECO_10_15","goodRECO_15_20","goodRECO_20_25","goodRECO_25_30","goodRECO_30_40",
// 			  "goodRECO_40_50","goodRECO_50_70","goodRECO_70_120","goodRECO_120_250",
// 			  "goodRECO_250_400","goodRECO_400_650","goodRECO_650_1000",
// 			  "goodRECO_1000_1500","goodRECO_1500_2000",
// 			  "badRECO_all",
// 			  "badRECO_10_15","badRECO_15_20","badRECO_20_25","badRECO_25_30","badRECO_30_40",
// 			  "badRECO_40_50","badRECO_50_70","badRECO_70_120","badRECO_120_250",
// 			  "badRECO_250_400","badRECO_400_650","badRECO_650_1000",
// 			  "badRECO_1000_1500","badRECO_1500_2000",
//   };

//   TString ptbinsLabels[30] = {"10 #leqp^{RECO}_{T}#leq #infty",
// 			     "10 #leqp^{RECO}_{T}#leq 15",
// 			     "15 #leqp^{RECO}_{T}#leq 20",
// 			     "20 #leqp^{RECO}_{T}#leq 25",
// 			     "25 #leqp^{RECO}_{T}#leq 30",
// 			     "30 #leqp^{RECO}_{T}#leq 40",
// 			     "40 #leqp^{RECO}_{T}#leq 50", "50 #leqp^{RECO}_{T}#leq 70", "70 #leqp^{RECO}_{T}#leq 120", 
// 			     "120 #leqp^{RECO}_{T}#leq 250", "250 #leqp^{RECO}_{T}#leq 400", "400 #leqp^{RECO}_{T}#leq 650",
// 			     "650 #leqp^{RECO}_{T}#leq 1000", "1000 #leqp^{RECO}_{T}#leq 1500",
// 			     "1500#leqp^{RECO}_{T}#leq 2000",
// 			     "10 #leqp^{RECO}_{T}#leq #infty",
// 			     "10 #leqp^{RECO}_{T}#leq 15",
// 			     "15 #leqp^{RECO}_{T}#leq 20",
// 			     "20 #leqp^{RECO}_{T}#leq 25",
// 			     "25 #leqp^{RECO}_{T}#leq 30",
// 			     "30 #leqp^{RECO}_{T}#leq 40",
// 			     "40 #leqp^{RECO}_{T}#leq 50", "50 #leqp^{RECO}_{T}#leq 70", "70 #leqp^{RECO}_{T}#leq 120", 
// 			     "120 #leqp^{RECO}_{T}#leq 250", "250 #leqp^{RECO}_{T}#leq 400", "400 #leqp^{RECO}_{T}#leq 650",
// 			     "650 #leqp^{RECO}_{T}#leq 1000", "1000 #leqp^{RECO}_{T}#leq 1500",
// 			     "1500#leqp^{RECO}_{T}#leq 2000",
// };
 TString dirNames[15] = {"RECO_all", 
			  "RECO_10_15","RECO_15_20","RECO_20_25","RECO_25_30","RECO_30_40",
			  "RECO_40_50","RECO_50_70","RECO_70_120","RECO_120_250",
			  "RECO_250_400","RECO_400_650","RECO_650_1000",
			  "RECO_1000_1500","RECO_1500_2000",
  };

  TString ptbinsLabels[15] = {"10 #leqp^{RECO}_{T}#leq #infty",
			     "10 #leqp^{RECO}_{T}#leq 15",
			     "15 #leqp^{RECO}_{T}#leq 20",
			     "20 #leqp^{RECO}_{T}#leq 25",
			     "25 #leqp^{RECO}_{T}#leq 30",
			     "30 #leqp^{RECO}_{T}#leq 40",
			     "40 #leqp^{RECO}_{T}#leq 50", "50 #leqp^{RECO}_{T}#leq 70", "70 #leqp^{RECO}_{T}#leq 120", 
			     "120 #leqp^{RECO}_{T}#leq 250", "250 #leqp^{RECO}_{T}#leq 400", "400 #leqp^{RECO}_{T}#leq 650",
			     "650 #leqp^{RECO}_{T}#leq 1000", "1000 #leqp^{RECO}_{T}#leq 1500",
			     "1500#leqp^{RECO}_{T}#leq 2000",
};

  int loop_inter = 1;
  if(forEveryPtBin) loop_inter=15;
  //  if(forEveryPtBin) loop_inter=1; //TEST
  
  for(int i = 0; i<loop_inter ; i++){
  
   TString dirName = (TString)dirNames[i];  
   TString ptbinLabel = (TString)ptbinsLabels[i];
   cout<<endl;  
   cout<<dirName<<endl;
      
   TString MCtitle = "MC";
   TString SavePlots = CorrectionObject::_outpath + "plots/control/ControlPlotsPUjets/CP_RecoGenMatched_" + dirName + "_" + CorrectionObject::_generator_tag;
   CorrectionObject::make_path(std::string((_outpath + "plots/").Data()));
   CorrectionObject::make_path(std::string((_outpath + "plots/control/").Data()));
   CorrectionObject::make_path(std::string((_outpath + "plots/control/ControlPlotsPUjets/").Data()));

   cout<<"++++++++++++++++ Collect all histograms ++++++++++++++++\n";
   TH2D* Jet_eta_pt_MC = (TH2D*)CorrectionObject::_MCFile->Get("good"+dirName+"/Jet_eta_pt");
   TH2D* Jet_eta_PUid_MC = (TH2D*)CorrectionObject::_MCFile->Get("good"+dirName+"/Jet_eta_PUid");


   TH2D* Jet_eta_pt_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data"+dirName+"/Jet_eta_pt");
   TH2D* Jet_eta_PUid_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data"+dirName+"/Jet_eta_PUid");

   TH2D* Jet_eta_pt_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad"+dirName+"/Jet_eta_pt");
   TH2D* Jet_eta_PUid_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad"+dirName+"/Jet_eta_PUid");

   
   TH1F* nPu_MC = (TH1F*)CorrectionObject::_MCFile->Get("good"+dirName+"/nPu");
   TH1F* N_PV_MC = (TH1F*)CorrectionObject::_MCFile->Get("good"+dirName+"/N_PV");
   TH1F* PU_id_MC = (TH1F*)CorrectionObject::_MCFile->Get("good"+dirName+"/PU_id");
   TH1F* PU_id_BB1_MC = (TH1F*)CorrectionObject::_MCFile->Get("good"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_MC = (TH1F*)CorrectionObject::_MCFile->Get("good"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_MC = (TH1F*)CorrectionObject::_MCFile->Get("good"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_MC = (TH1F*)CorrectionObject::_MCFile->Get("good"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_MC = (TH1F*)CorrectionObject::_MCFile->Get("good"+dirName+"/PU_id_HF");

   TH1F* nPu_notmatchedMC = (TH1F*)CorrectionObject::_MCFile->Get("bad"+dirName+"/nPu");
   TH1F* N_PV_notmatchedMC = (TH1F*)CorrectionObject::_MCFile->Get("bad"+dirName+"/N_PV");
   TH1F* PU_id_notmatchedMC = (TH1F*)CorrectionObject::_MCFile->Get("bad"+dirName+"/PU_id");
   TH1F* PU_id_BB1_notmatchedMC = (TH1F*)CorrectionObject::_MCFile->Get("bad"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_notmatchedMC = (TH1F*)CorrectionObject::_MCFile->Get("bad"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_notmatchedMC = (TH1F*)CorrectionObject::_MCFile->Get("bad"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_notmatchedMC = (TH1F*)CorrectionObject::_MCFile->Get("bad"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_notmatchedMC = (TH1F*)CorrectionObject::_MCFile->Get("bad"+dirName+"/PU_id_HF");

   TH1F* nPu_MC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/nPu");
   TH1F* N_PV_MC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/N_PV");
   TH1F* PU_id_MC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/PU_id");
   TH1F* PU_id_BB1_MC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_MC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_MC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_MC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_MC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/PU_id_HF");

   TH1F* nPu_notmatchedMC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/nPu");
   TH1F* N_PV_notmatchedMC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/N_PV");
   TH1F* PU_id_notmatchedMC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/PU_id");
   TH1F* PU_id_BB1_notmatchedMC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_notmatchedMC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_notmatchedMC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_notmatchedMC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_notmatchedMC_3rd = (TH1F*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/PU_id_HF");

   TH1F* PU_id_MCparton= (TH1F*)CorrectionObject::_MCFile->Get("goodparton"+dirName+"/PU_id");
   TH1F* PU_id_BB1_MCparton= (TH1F*)CorrectionObject::_MCFile->Get("goodparton"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_MCparton= (TH1F*)CorrectionObject::_MCFile->Get("goodparton"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_MCparton= (TH1F*)CorrectionObject::_MCFile->Get("goodparton"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_MCparton= (TH1F*)CorrectionObject::_MCFile->Get("goodparton"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_MCparton= (TH1F*)CorrectionObject::_MCFile->Get("goodparton"+dirName+"/PU_id_HF");

   TH1F* nPu_notmatchedMCparton = (TH1F*)CorrectionObject::_MCFile->Get("badparton"+dirName+"/nPu");
   TH1F* N_PV_notmatchedMCparton = (TH1F*)CorrectionObject::_MCFile->Get("badparton"+dirName+"/N_PV");
   TH1F* PU_id_notmatchedMCparton = (TH1F*)CorrectionObject::_MCFile->Get("badparton"+dirName+"/PU_id");
   TH1F* PU_id_BB1_notmatchedMCparton = (TH1F*)CorrectionObject::_MCFile->Get("badparton"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_notmatchedMCparton = (TH1F*)CorrectionObject::_MCFile->Get("badparton"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_notmatchedMCparton = (TH1F*)CorrectionObject::_MCFile->Get("badparton"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_notmatchedMCparton = (TH1F*)CorrectionObject::_MCFile->Get("badparton"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_notmatchedMCparton = (TH1F*)CorrectionObject::_MCFile->Get("badparton"+dirName+"/PU_id_HF");


   TH1F* PU_id_MCparton_3rd= (TH1F*)CorrectionObject::_MCFile->Get("goodparton3rd"+dirName+"/PU_id");
   TH1F* PU_id_BB1_MCparton_3rd= (TH1F*)CorrectionObject::_MCFile->Get("goodparton3rd"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_MCparton_3rd= (TH1F*)CorrectionObject::_MCFile->Get("goodparton3rd"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_MCparton_3rd= (TH1F*)CorrectionObject::_MCFile->Get("goodparton3rd"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_MCparton_3rd= (TH1F*)CorrectionObject::_MCFile->Get("goodparton3rd"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_MCparton_3rd= (TH1F*)CorrectionObject::_MCFile->Get("goodparton3rd"+dirName+"/PU_id_HF");

   TH1F* nPu_notmatchedMCparton_3rd = (TH1F*)CorrectionObject::_MCFile->Get("badparton3rd"+dirName+"/nPu");
   TH1F* N_PV_notmatchedMCparton_3rd = (TH1F*)CorrectionObject::_MCFile->Get("badparton3rd"+dirName+"/N_PV");
   TH1F* PU_id_notmatchedMCparton_3rd = (TH1F*)CorrectionObject::_MCFile->Get("badparton3rd"+dirName+"/PU_id");
   TH1F* PU_id_BB1_notmatchedMCparton_3rd = (TH1F*)CorrectionObject::_MCFile->Get("badparton3rd"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_notmatchedMCparton_3rd = (TH1F*)CorrectionObject::_MCFile->Get("badparton3rd"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_notmatchedMCparton_3rd = (TH1F*)CorrectionObject::_MCFile->Get("badparton3rd"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_notmatchedMCparton_3rd = (TH1F*)CorrectionObject::_MCFile->Get("badparton3rd"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_notmatchedMCparton_3rd = (TH1F*)CorrectionObject::_MCFile->Get("badparton3rd"+dirName+"/PU_id_HF");


   TH1F* nPu_DATA = (TH1F*)CorrectionObject::_DATAFile->Get("data"+dirName+"/nPu");
   TH1F* N_PV_DATA = (TH1F*)CorrectionObject::_DATAFile->Get("data"+dirName+"/N_PV");
   TH1F* PU_id_DATA = (TH1F*)CorrectionObject::_DATAFile->Get("data"+dirName+"/PU_id");
   TH1F* PU_id_BB1_DATA = (TH1F*)CorrectionObject::_DATAFile->Get("data"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_DATA = (TH1F*)CorrectionObject::_DATAFile->Get("data"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_DATA = (TH1F*)CorrectionObject::_DATAFile->Get("data"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_DATA = (TH1F*)CorrectionObject::_DATAFile->Get("data"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_DATA = (TH1F*)CorrectionObject::_DATAFile->Get("data"+dirName+"/PU_id_HF");

   TH1F* nPu_DATA_3rd = (TH1F*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/nPu");
   TH1F* N_PV_DATA_3rd = (TH1F*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/N_PV");
   TH1F* PU_id_DATA_3rd = (TH1F*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/PU_id");
   TH1F* PU_id_BB1_DATA_3rd = (TH1F*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/PU_id_BB1");
   TH1F* PU_id_BB2_DATA_3rd = (TH1F*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/PU_id_BB2");
   TH1F* PU_id_EC1_DATA_3rd = (TH1F*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/PU_id_EC1");
   TH1F* PU_id_EC2_DATA_3rd = (TH1F*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/PU_id_EC2");
   TH1F* PU_id_HF_DATA_3rd = (TH1F*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/PU_id_HF");

   TH2D* Jet_eta_chEM_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_eta_chEM");
   TH2D* Jet_eta_neuEM_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_eta_neuEM");
   TH2D* Jet_eta_chHAD_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_eta_chHAD");
   TH2D* Jet_eta_neuHAD_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_eta_neuHAD");

   TH2D* Jet_eta_chEM_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_eta_chEM");
   TH2D* Jet_eta_neuEM_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_eta_neuEM");
   TH2D* Jet_eta_chHAD_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_eta_chHAD");
   TH2D* Jet_eta_neuHAD_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_eta_neuHAD");

   TH2D* Jet_eta_chEM_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_eta_chEM");
   TH2D* Jet_eta_neuEM_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_eta_neuEM");
   TH2D* Jet_eta_chHAD_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_eta_chHAD");
   TH2D* Jet_eta_neuHAD_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_eta_neuHAD");


   //PI id vs PFs
   TH2D* Jet_PU_id_chEM_BB1_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_chEM_BB1");
   TH2D* Jet_PU_id_neuEM_BB1_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_neuEM_BB1");
   TH2D* Jet_PU_id_chHAD_BB1_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_chHAD_BB1");
   TH2D* Jet_PU_id_neuHAD_BB1_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_neuHAD_BB1");

   TH2D* Jet_PU_id_chEM_BB1_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_chEM_BB1");
   TH2D* Jet_PU_id_neuEM_BB1_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_neuEM_BB1");
   TH2D* Jet_PU_id_chHAD_BB1_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_chHAD_BB1");
   TH2D* Jet_PU_id_neuHAD_BB1_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_neuHAD_BB1");

   TH2D* Jet_PU_id_chEM_BB1_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_chEM_BB1");
   TH2D* Jet_PU_id_neuEM_BB1_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_neuEM_BB1");
   TH2D* Jet_PU_id_chHAD_BB1_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_chHAD_BB1");
   TH2D* Jet_PU_id_neuHAD_BB1_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_neuHAD_BB1");

   TH2D* Jet_PU_id_chEM_BB2_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_chEM_BB2");
   TH2D* Jet_PU_id_neuEM_BB2_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_neuEM_BB2");
   TH2D* Jet_PU_id_chHAD_BB2_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_chHAD_BB2");
   TH2D* Jet_PU_id_neuHAD_BB2_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_neuHAD_BB2");

   TH2D* Jet_PU_id_chEM_BB2_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_chEM_BB2");
   TH2D* Jet_PU_id_neuEM_BB2_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_neuEM_BB2");
   TH2D* Jet_PU_id_chHAD_BB2_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_chHAD_BB2");
   TH2D* Jet_PU_id_neuHAD_BB2_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_neuHAD_BB2");

   TH2D* Jet_PU_id_chEM_BB2_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_chEM_BB2");
   TH2D* Jet_PU_id_neuEM_BB2_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_neuEM_BB2");
   TH2D* Jet_PU_id_chHAD_BB2_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_chHAD_BB2");
   TH2D* Jet_PU_id_neuHAD_BB2_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_neuHAD_BB2");

   TH2D* Jet_PU_id_chEM_EC1_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_chEM_EC1");
   TH2D* Jet_PU_id_neuEM_EC1_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_neuEM_EC1");
   TH2D* Jet_PU_id_chHAD_EC1_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_chHAD_EC1");
   TH2D* Jet_PU_id_neuHAD_EC1_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_neuHAD_EC1");

   TH2D* Jet_PU_id_chEM_EC1_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_chEM_EC1");
   TH2D* Jet_PU_id_neuEM_EC1_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_neuEM_EC1");
   TH2D* Jet_PU_id_chHAD_EC1_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_chHAD_EC1");
   TH2D* Jet_PU_id_neuHAD_EC1_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_neuHAD_EC1");

   TH2D* Jet_PU_id_chEM_EC1_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_chEM_EC1");
   TH2D* Jet_PU_id_neuEM_EC1_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_neuEM_EC1");
   TH2D* Jet_PU_id_chHAD_EC1_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_chHAD_EC1");
   TH2D* Jet_PU_id_neuHAD_EC1_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_neuHAD_EC1");

   TH2D* Jet_PU_id_chEM_EC2_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_chEM_EC2");
   TH2D* Jet_PU_id_neuEM_EC2_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_neuEM_EC2");
   TH2D* Jet_PU_id_chHAD_EC2_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_chHAD_EC2");
   TH2D* Jet_PU_id_neuHAD_EC2_MC = (TH2D*)CorrectionObject::_MCFile->Get("good3rd"+dirName+"/Jet_PU_id_neuHAD_EC2");

   TH2D* Jet_PU_id_chEM_EC2_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_chEM_EC2");
   TH2D* Jet_PU_id_neuEM_EC2_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_neuEM_EC2");
   TH2D* Jet_PU_id_chHAD_EC2_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_chHAD_EC2");
   TH2D* Jet_PU_id_neuHAD_EC2_notmatchedMC = (TH2D*)CorrectionObject::_MCFile->Get("bad3rd"+dirName+"/Jet_PU_id_neuHAD_EC2");

   TH2D* Jet_PU_id_chEM_EC2_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_chEM_EC2");
   TH2D* Jet_PU_id_neuEM_EC2_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_neuEM_EC2");
   TH2D* Jet_PU_id_chHAD_EC2_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_chHAD_EC2");
   TH2D* Jet_PU_id_neuHAD_EC2_DATA = (TH2D*)CorrectionObject::_DATAFile->Get("data3rd"+dirName+"/Jet_PU_id_neuHAD_EC2");


   Jet_eta_pt_MC->SetName("Jet_eta_pt_MC");
   Jet_eta_PUid_MC->SetName("Jet_eta_PUid_MC");


   nPu_MC->SetName("nPu_MC");
   N_PV_MC->SetName("N_PV_MC");
   PU_id_MC->SetName("PU_id_MC");
   PU_id_BB1_MC->SetName("PU_id_BB1_MC");
   PU_id_BB2_MC->SetName("PU_id_BB2_MC");
   PU_id_EC1_MC->SetName("PU_id_EC1_MC");
   PU_id_EC2_MC->SetName("PU_id_EC2_MC");
   PU_id_HF_MC->SetName("PU_id_HF_MC");

   nPu_DATA->SetName("nPu_DATA");
   N_PV_DATA->SetName("N_PV_DATA");
   PU_id_DATA->SetName("PU_id_DATA");
   PU_id_BB1_DATA->SetName("PU_id_BB1_DATA");
   PU_id_BB2_DATA->SetName("PU_id_BB2_DATA");
   PU_id_EC1_DATA->SetName("PU_id_EC1_DATA");
   PU_id_EC2_DATA->SetName("PU_id_EC2_DATA");
   PU_id_HF_DATA->SetName("PU_id_HF_DATA");

   nPu_DATA_3rd->SetName("nPu_DATA_3rd");
   N_PV_DATA_3rd->SetName("N_PV_DATA_3rd");
   PU_id_DATA_3rd->SetName("PU_id_DATA_3rd");
   PU_id_BB1_DATA_3rd->SetName("PU_id_BB1_DATA_3rd");
   PU_id_BB2_DATA_3rd->SetName("PU_id_BB2_DATA_3rd");
   PU_id_EC1_DATA_3rd->SetName("PU_id_EC1_DATA_3rd");
   PU_id_EC2_DATA_3rd->SetName("PU_id_EC2_DATA_3rd");
   PU_id_HF_DATA_3rd->SetName("PU_id_HF_DATA_3rd");

   nPu_notmatchedMC->SetName("nPu_notmatchedMC");
   N_PV_notmatchedMC->SetName("N_PV_notmatchedMC");
   PU_id_notmatchedMC->SetName("PU_id_notmatchedMC");
   PU_id_BB1_notmatchedMC->SetName("PU_id_BB1_notmatchedMC");
   PU_id_BB2_notmatchedMC->SetName("PU_id_BB2_notmatchedMC");
   PU_id_EC1_notmatchedMC->SetName("PU_id_EC1_notmatchedMC");
   PU_id_EC2_notmatchedMC->SetName("PU_id_EC2_notmatchedMC");
   PU_id_HF_notmatchedMC->SetName("PU_id_HF_notmatchedMC");

   nPu_notmatchedMC_3rd->SetName("nPu_notmatchedMC_3rd");
   N_PV_notmatchedMC_3rd->SetName("N_PV_notmatchedMC_3rd");
   PU_id_notmatchedMC_3rd->SetName("PU_id_notmatchedMC_3rd");
   PU_id_BB1_notmatchedMC_3rd->SetName("PU_id_BB1_notmatchedMC_3rd");
   PU_id_BB2_notmatchedMC_3rd->SetName("PU_id_BB2_notmatchedMC_3rd");
   PU_id_EC1_notmatchedMC_3rd->SetName("PU_id_EC1_notmatchedMC_3rd");
   PU_id_EC2_notmatchedMC_3rd->SetName("PU_id_EC2_notmatchedMC_3rd");
   PU_id_HF_notmatchedMC_3rd->SetName("PU_id_HF_notmatchedMC_3rd");

   nPu_MC_3rd->SetName("nPu_MC_3rd");
   N_PV_MC_3rd->SetName("N_PV_MC_3rd");
   PU_id_MC_3rd->SetName("PU_id_MC_3rd");
   PU_id_BB1_MC_3rd->SetName("PU_id_BB1_MC_3rd");
   PU_id_BB2_MC_3rd->SetName("PU_id_BB2_MC_3rd");
   PU_id_EC1_MC_3rd->SetName("PU_id_EC1_MC_3rd");
   PU_id_EC2_MC_3rd->SetName("PU_id_EC2_MC_3rd");
   PU_id_HF_MC_3rd->SetName("PU_id_HF_MC_3rd");
   

   Jet_PU_id_chEM_BB1_MC->SetName("Jet_PU_id_chEM_BB1_good"); 
   Jet_PU_id_neuEM_BB1_MC->SetName("Jet_PU_id_neuEM_BB1_good");
   Jet_PU_id_chHAD_BB1_MC->SetName("Jet_PU_id_chHAD_BB1_good");
   Jet_PU_id_neuHAD_BB1_MC->SetName("Jet_PU_id_neuHAD_BB1_good");
   Jet_PU_id_chEM_BB1_notmatchedMC->SetName("Jet_PU_id_chEM_BB1_bad"); 
   Jet_PU_id_neuEM_BB1_notmatchedMC->SetName("Jet_PU_id_neuEM_BB1_bad");
   Jet_PU_id_chHAD_BB1_notmatchedMC->SetName("Jet_PU_id_chHAD_BB1_bad");
   Jet_PU_id_neuHAD_BB1_notmatchedMC->SetName("Jet_PU_id_neuHAD_BB1_bad");
   Jet_PU_id_chEM_BB1_DATA->SetName("Jet_PU_id_chEM_BB1_data");
   Jet_PU_id_neuEM_BB1_DATA->SetName("Jet_PU_id_neuEM_BB1_data");
   Jet_PU_id_chHAD_BB1_DATA->SetName("Jet_PU_id_chHAD_BB1_data");
   Jet_PU_id_neuHAD_BB1_DATA->SetName("Jet_PU_id_neuHAD_BB1_data");

   Jet_PU_id_chEM_BB2_MC->SetName("Jet_PU_id_chEM_BB2_good"); 
   Jet_PU_id_neuEM_BB2_MC->SetName("Jet_PU_id_neuEM_BB2_good");
   Jet_PU_id_chHAD_BB2_MC->SetName("Jet_PU_id_chHAD_BB2_good");
   Jet_PU_id_neuHAD_BB2_MC->SetName("Jet_PU_id_neuHAD_BB2_good");
   Jet_PU_id_chEM_BB2_notmatchedMC->SetName("Jet_PU_id_chEM_BB2_bad"); 
   Jet_PU_id_neuEM_BB2_notmatchedMC->SetName("Jet_PU_id_neuEM_BB2_bad");
   Jet_PU_id_chHAD_BB2_notmatchedMC->SetName("Jet_PU_id_chHAD_BB2_bad");
   Jet_PU_id_neuHAD_BB2_notmatchedMC->SetName("Jet_PU_id_neuHAD_BB2_bad");
   Jet_PU_id_chEM_BB2_DATA->SetName("Jet_PU_id_chEM_BB2_data");
   Jet_PU_id_neuEM_BB2_DATA->SetName("Jet_PU_id_neuEM_BB2_data");
   Jet_PU_id_chHAD_BB2_DATA->SetName("Jet_PU_id_chHAD_BB2_data");
   Jet_PU_id_neuHAD_BB2_DATA->SetName("Jet_PU_id_neuHAD_BB2_data");


   Jet_PU_id_chEM_EC1_MC->SetName("Jet_PU_id_chEM_EC1_good"); 
   Jet_PU_id_neuEM_EC1_MC->SetName("Jet_PU_id_neuEM_EC1_good");
   Jet_PU_id_chHAD_EC1_MC->SetName("Jet_PU_id_chHAD_EC1_good");
   Jet_PU_id_neuHAD_EC1_MC->SetName("Jet_PU_id_neuHAD_EC1_good");
   Jet_PU_id_chEM_EC1_notmatchedMC->SetName("Jet_PU_id_chEM_EC1_bad"); 
   Jet_PU_id_neuEM_EC1_notmatchedMC->SetName("Jet_PU_id_neuEM_EC1_bad");
   Jet_PU_id_chHAD_EC1_notmatchedMC->SetName("Jet_PU_id_chHAD_EC1_bad");
   Jet_PU_id_neuHAD_EC1_notmatchedMC->SetName("Jet_PU_id_neuHAD_EC1_bad");
   Jet_PU_id_chEM_EC1_DATA->SetName("Jet_PU_id_chEM_EC1_data");
   Jet_PU_id_neuEM_EC1_DATA->SetName("Jet_PU_id_neuEM_EC1_data");
   Jet_PU_id_chHAD_EC1_DATA->SetName("Jet_PU_id_chHAD_EC1_data");
   Jet_PU_id_neuHAD_EC1_DATA->SetName("Jet_PU_id_neuHAD_EC1_data");

   Jet_PU_id_chEM_EC2_MC->SetName("Jet_PU_id_chEM_EC2_good"); 
   Jet_PU_id_neuEM_EC2_MC->SetName("Jet_PU_id_neuEM_EC2_good");
   Jet_PU_id_chHAD_EC2_MC->SetName("Jet_PU_id_chHAD_EC2_good");
   Jet_PU_id_neuHAD_EC2_MC->SetName("Jet_PU_id_neuHAD_EC2_good");
   Jet_PU_id_chEM_EC2_notmatchedMC->SetName("Jet_PU_id_chEM_EC2_bad"); 
   Jet_PU_id_neuEM_EC2_notmatchedMC->SetName("Jet_PU_id_neuEM_EC2_bad");
   Jet_PU_id_chHAD_EC2_notmatchedMC->SetName("Jet_PU_id_chHAD_EC2_bad");
   Jet_PU_id_neuHAD_EC2_notmatchedMC->SetName("Jet_PU_id_neuHAD_EC2_bad");
   Jet_PU_id_chEM_EC2_DATA->SetName("Jet_PU_id_chEM_EC2_data");
   Jet_PU_id_neuEM_EC2_DATA->SetName("Jet_PU_id_neuEM_EC2_data");
   Jet_PU_id_chHAD_EC2_DATA->SetName("Jet_PU_id_chHAD_EC2_data");
   Jet_PU_id_neuHAD_EC2_DATA->SetName("Jet_PU_id_neuHAD_EC2_data");

   TLatex *tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.045);  

   cout<<"+++++++++++++++++++++ Plot extremely nice histograms ++++++++++++++++++++++++++++\n";
   TCanvas* cJet_eta_pt_MC = new TCanvas();
   //  Jet_eta_pt_MC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_pt_MC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_pt_MC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_pt_MC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_pt_MC->SetTitle("");
   Jet_eta_pt_MC->Draw("colz");
   cJet_eta_pt_MC->SetLogy();
   cJet_eta_pt_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_pt_MC->Print(SavePlots + "_Jet_eta_pt.pdf","pdf");

   TCanvas* cJet_eta_PUid_MC = new TCanvas();
   //  Jet_eta_PUid_MC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_PUid_MC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_PUid_MC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_PUid_MC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_PUid_MC->SetTitle("");
   Jet_eta_PUid_MC->Draw("colz");
   cJet_eta_PUid_MC->SetLogy();
   cJet_eta_PUid_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_PUid_MC->Print(SavePlots + "_Jet_eta_PUid.pdf","pdf");


   //GOOD
   TCanvas* cJet_eta_chEM_MC = new TCanvas();
   //  Jet_eta_chEM_MC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_chEM_MC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_chEM_MC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_chEM_MC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_chEM_MC->SetTitle("");
   Jet_eta_chEM_MC->Draw("colz");
   cJet_eta_chEM_MC->SetLogy();
   cJet_eta_chEM_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_chEM_MC->Print(SavePlots + "_good_Jet_eta_chEM.pdf","pdf");

   TCanvas* cJet_eta_neuEM_MC = new TCanvas();
   //  Jet_eta_neuEM_MC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_neuEM_MC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_neuEM_MC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_neuEM_MC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_neuEM_MC->SetTitle("");
   Jet_eta_neuEM_MC->Draw("colz");
   cJet_eta_neuEM_MC->SetLogy();
   cJet_eta_neuEM_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_neuEM_MC->Print(SavePlots + "_good_Jet_eta_neuEM.pdf","pdf");

   TCanvas* cJet_eta_chHAD_MC = new TCanvas();
   //  Jet_eta_chHAD_MC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_chHAD_MC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_chHAD_MC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_chHAD_MC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_chHAD_MC->SetTitle("");
   Jet_eta_chHAD_MC->Draw("colz");
   cJet_eta_chHAD_MC->SetLogy();
   cJet_eta_chHAD_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_chHAD_MC->Print(SavePlots + "_good_Jet_eta_chHAD.pdf","pdf");

   TCanvas* cJet_eta_neuHAD_MC = new TCanvas();
   //  Jet_eta_neuHAD_MC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_neuHAD_MC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_neuHAD_MC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_neuHAD_MC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_neuHAD_MC->SetTitle("");
   Jet_eta_neuHAD_MC->Draw("colz");
   cJet_eta_neuHAD_MC->SetLogy();
   cJet_eta_neuHAD_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_neuHAD_MC->Print(SavePlots + "_good_Jet_eta_neuHAD.pdf","pdf");

   //BAD
   TCanvas* cJet_eta_chEM_notmatchedMC = new TCanvas();
   //  Jet_eta_chEM_notmatchedMC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_chEM_notmatchedMC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_chEM_notmatchedMC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_chEM_notmatchedMC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_chEM_notmatchedMC->SetTitle("");
   Jet_eta_chEM_notmatchedMC->Draw("colz");
   cJet_eta_chEM_notmatchedMC->SetLogy();
   cJet_eta_chEM_notmatchedMC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_chEM_notmatchedMC->Print(SavePlots + "_bad_Jet_eta_chEM.pdf","pdf");

   TCanvas* cJet_eta_neuEM_notmatchedMC = new TCanvas();
   //  Jet_eta_neuEM_notmatchedMC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_neuEM_notmatchedMC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_neuEM_notmatchedMC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_neuEM_notmatchedMC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_neuEM_notmatchedMC->SetTitle("");
   Jet_eta_neuEM_notmatchedMC->Draw("colz");
   cJet_eta_neuEM_notmatchedMC->SetLogy();
   cJet_eta_neuEM_notmatchedMC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_neuEM_notmatchedMC->Print(SavePlots + "_bad_Jet_eta_neuEM.pdf","pdf");

   TCanvas* cJet_eta_chHAD_notmatchedMC = new TCanvas();
   //  Jet_eta_chHAD_notmatchedMC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_chHAD_notmatchedMC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_chHAD_notmatchedMC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_chHAD_notmatchedMC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_chHAD_notmatchedMC->SetTitle("");
   Jet_eta_chHAD_notmatchedMC->Draw("colz");
   cJet_eta_chHAD_notmatchedMC->SetLogy();
   cJet_eta_chHAD_notmatchedMC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_chHAD_notmatchedMC->Print(SavePlots + "_bad_Jet_eta_chHAD.pdf","pdf");

   TCanvas* cJet_eta_neuHAD_notmatchedMC = new TCanvas();
   //  Jet_eta_neuHAD_notmatchedMC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_neuHAD_notmatchedMC->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_neuHAD_notmatchedMC->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_neuHAD_notmatchedMC->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_neuHAD_notmatchedMC->SetTitle("");
   Jet_eta_neuHAD_notmatchedMC->Draw("colz");
   cJet_eta_neuHAD_notmatchedMC->SetLogy();
   cJet_eta_neuHAD_notmatchedMC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag+", 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_neuHAD_notmatchedMC->Print(SavePlots + "_bad_Jet_eta_neuHAD.pdf","pdf");

   // DATA
   TCanvas* cJet_eta_chEM_DATA = new TCanvas();
   //  Jet_eta_chEM_DATA->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_chEM_DATA->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_chEM_DATA->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_chEM_DATA->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_chEM_DATA->SetTitle("");
   Jet_eta_chEM_DATA->Draw("colz");
   cJet_eta_chEM_DATA->SetLogy();
   cJet_eta_chEM_DATA->SetLogz();
   tex->DrawLatex(0.15,0.95,"DATA, 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_chEM_DATA->Print(SavePlots + "_data_Jet_eta_chEM.pdf","pdf");

   TCanvas* cJet_eta_neuEM_DATA = new TCanvas();
   //  Jet_eta_neuEM_DATA->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_neuEM_DATA->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_neuEM_DATA->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_neuEM_DATA->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_neuEM_DATA->SetTitle("");
   Jet_eta_neuEM_DATA->Draw("colz");
   cJet_eta_neuEM_DATA->SetLogy();
   cJet_eta_neuEM_DATA->SetLogz();
   tex->DrawLatex(0.15,0.95,"DATA, 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_neuEM_DATA->Print(SavePlots + "_data_Jet_eta_neuEM.pdf","pdf");

   TCanvas* cJet_eta_chHAD_DATA = new TCanvas();
   //  Jet_eta_chHAD_DATA->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_chHAD_DATA->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_chHAD_DATA->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_chHAD_DATA->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_chHAD_DATA->SetTitle("");
   Jet_eta_chHAD_DATA->Draw("colz");
   cJet_eta_chHAD_DATA->SetLogy();
   cJet_eta_chHAD_DATA->SetLogz();
   tex->DrawLatex(0.15,0.95,"DATA, 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_chHAD_DATA->Print(SavePlots + "_data_Jet_eta_chHAD.pdf","pdf");

   TCanvas* cJet_eta_neuHAD_DATA = new TCanvas();
   //  Jet_eta_neuHAD_DATA->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
   Jet_eta_neuHAD_DATA->GetXaxis()->SetLabelSize(0.045);
   Jet_eta_neuHAD_DATA->GetYaxis()->SetLabelSize(0.045);
   Jet_eta_neuHAD_DATA->GetYaxis()->SetTitleSize(0.045);
   Jet_eta_neuHAD_DATA->SetTitle("");
   Jet_eta_neuHAD_DATA->Draw("colz");
   cJet_eta_neuHAD_DATA->SetLogy();
   cJet_eta_neuHAD_DATA->SetLogz();
   tex->DrawLatex(0.15,0.95,"DATA, 3rd jet");
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cJet_eta_neuHAD_DATA->Print(SavePlots + "_data_Jet_eta_neuHAD.pdf","pdf");


   //PU id vs PFs
  
   doPlot_2D_PIid_PF(Jet_PU_id_chEM_BB1_MC, Jet_PU_id_chEM_BB1_notmatchedMC, Jet_PU_id_chEM_BB1_DATA, "chEM", ptbinLabel, txttag, SavePlots,0);
   doPlot_2D_PIid_PF(Jet_PU_id_chEM_BB2_MC, Jet_PU_id_chEM_BB2_notmatchedMC, Jet_PU_id_chEM_BB2_DATA, "chEM", ptbinLabel, txttag, SavePlots,1);
   doPlot_2D_PIid_PF(Jet_PU_id_chEM_EC1_MC, Jet_PU_id_chEM_EC1_notmatchedMC, Jet_PU_id_chEM_EC1_DATA, "chEM", ptbinLabel, txttag, SavePlots,2);
   doPlot_2D_PIid_PF(Jet_PU_id_chEM_EC2_MC, Jet_PU_id_chEM_EC2_notmatchedMC, Jet_PU_id_chEM_EC2_DATA, "chEM", ptbinLabel, txttag, SavePlots,3);

   doPlot_2D_PIid_PF(Jet_PU_id_neuEM_BB1_MC, Jet_PU_id_neuEM_BB1_notmatchedMC, Jet_PU_id_neuEM_BB1_DATA, "neuEM", ptbinLabel, txttag, SavePlots,0);
   doPlot_2D_PIid_PF(Jet_PU_id_neuEM_BB2_MC, Jet_PU_id_neuEM_BB2_notmatchedMC, Jet_PU_id_neuEM_BB2_DATA, "neuEM", ptbinLabel, txttag, SavePlots,1);
   doPlot_2D_PIid_PF(Jet_PU_id_neuEM_EC1_MC, Jet_PU_id_neuEM_EC1_notmatchedMC, Jet_PU_id_neuEM_EC1_DATA, "neuEM", ptbinLabel, txttag, SavePlots,2);
   doPlot_2D_PIid_PF(Jet_PU_id_neuEM_EC2_MC, Jet_PU_id_neuEM_EC2_notmatchedMC, Jet_PU_id_neuEM_EC2_DATA, "neuEM", ptbinLabel, txttag, SavePlots,3);

   doPlot_2D_PIid_PF(Jet_PU_id_chHAD_BB1_MC, Jet_PU_id_chHAD_BB1_notmatchedMC, Jet_PU_id_chHAD_BB1_DATA, "chHAD", ptbinLabel, txttag, SavePlots,0);
   doPlot_2D_PIid_PF(Jet_PU_id_chHAD_BB2_MC, Jet_PU_id_chHAD_BB2_notmatchedMC, Jet_PU_id_chHAD_BB2_DATA, "chHAD", ptbinLabel, txttag, SavePlots,1);
   doPlot_2D_PIid_PF(Jet_PU_id_chHAD_EC1_MC, Jet_PU_id_chHAD_EC1_notmatchedMC, Jet_PU_id_chHAD_EC1_DATA, "chHAD", ptbinLabel, txttag, SavePlots,2);
   doPlot_2D_PIid_PF(Jet_PU_id_chHAD_EC2_MC, Jet_PU_id_chHAD_EC2_notmatchedMC, Jet_PU_id_chHAD_EC2_DATA, "chHAD", ptbinLabel, txttag, SavePlots,3);

   doPlot_2D_PIid_PF(Jet_PU_id_neuHAD_BB1_MC, Jet_PU_id_neuHAD_BB1_notmatchedMC, Jet_PU_id_neuHAD_BB1_DATA, "neuHAD", ptbinLabel, txttag, SavePlots,0);
   doPlot_2D_PIid_PF(Jet_PU_id_neuHAD_BB2_MC, Jet_PU_id_neuHAD_BB2_notmatchedMC, Jet_PU_id_neuHAD_BB2_DATA, "neuHAD", ptbinLabel, txttag, SavePlots,1);
   doPlot_2D_PIid_PF(Jet_PU_id_neuHAD_EC1_MC, Jet_PU_id_neuHAD_EC1_notmatchedMC, Jet_PU_id_neuHAD_EC1_DATA, "neuHAD", ptbinLabel, txttag, SavePlots,2);
   doPlot_2D_PIid_PF(Jet_PU_id_neuHAD_EC2_MC, Jet_PU_id_neuHAD_EC2_notmatchedMC, Jet_PU_id_neuHAD_EC2_DATA, "neuHAD", ptbinLabel, txttag, SavePlots,3);
   


   //////////////////////


   TCanvas* cnPu_MC = new TCanvas();
   nPu_MC->GetXaxis()->SetTitle("N_{PU}");
   nPu_MC->GetXaxis()->SetLabelSize(0.045);
   nPu_MC->GetYaxis()->SetLabelSize(0.045);
   nPu_MC->GetYaxis()->SetTitleSize(0.045);
   nPu_MC->SetMarkerColor(kRed);
   nPu_MC->SetMarkerStyle(20);
   nPu_MC->SetTitle("");
   nPu_MC->Draw();
   nPu_notmatchedMC->SetMarkerColor(kBlack);
   nPu_notmatchedMC->SetMarkerStyle(21);
   nPu_notmatchedMC->Draw("same");
   //   cnPu_MC->SetLogy();
   //   cnPu_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cnPu_MC->Print(SavePlots + "_nPu.pdf","pdf");

   TCanvas* cN_PV_MC = new TCanvas();
   N_PV_MC->GetXaxis()->SetTitle("N_{PV}");
   N_PV_MC->GetXaxis()->SetLabelSize(0.045);
   N_PV_MC->GetYaxis()->SetLabelSize(0.045);
   N_PV_MC->GetYaxis()->SetTitleSize(0.045);
   N_PV_MC->SetMarkerColor(kRed);
   N_PV_MC->SetMarkerStyle(20);
   N_PV_MC->SetTitle("");
   N_PV_MC->Draw();
   N_PV_notmatchedMC->SetMarkerColor(kBlack);
   N_PV_notmatchedMC->SetMarkerStyle(21);
   N_PV_notmatchedMC->Draw("same");
   //   cN_PV_MC->SetLogy();
   //   cN_PV_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   cN_PV_MC->Print(SavePlots + "_N_PV.pdf","pdf");

   //   TLegend* leg =  new TLegend(0.70,0.67,0.99,0.87);
   TLegend* leg =  new TLegend(0.17,0.67,0.52,0.89);
   leg->SetTextSize(0.045);
   TCanvas* cPU_id_MC = new TCanvas();
   PU_id_MC->GetXaxis()->SetTitle("PU ID");
   PU_id_MC->GetXaxis()->SetLabelSize(0.045);
   PU_id_MC->GetYaxis()->SetLabelSize(0.045);
   PU_id_MC->GetYaxis()->SetTitleSize(0.045);
   PU_id_MC->SetMarkerColor(kRed);
   PU_id_MC->SetMarkerStyle(20);
   PU_id_DATA->SetMarkerColor(kGreen);
   PU_id_DATA->SetMarkerStyle(22);
   PU_id_DATA_3rd->SetMarkerColor(kGreen+2);
   PU_id_DATA_3rd->SetMarkerStyle(21);

   PU_id_MC->SetTitle("");
   PU_id_MC->GetYaxis()->SetRangeUser(0.0,0.2);
   PU_id_MC->Draw();
   PU_id_notmatchedMC->SetMarkerColor(kBlack);
   PU_id_notmatchedMC->SetMarkerStyle(21);
   PU_id_notmatchedMC->Draw("same");
   leg->AddEntry(PU_id_MC,"MC, matched","p");
   leg->AddEntry(PU_id_notmatchedMC,"MC, not matched","p");
   leg->AddEntry(PU_id_DATA,"DATA, [1st,2nd] jets","p");
   //   leg->AddEntry(PU_id_DATA_3rd,"DATA, 3rd jet","p");
   //   cPU_id_MC->SetLogy();
   //   cPU_id_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   leg->Draw();
   cPU_id_MC->Print(SavePlots + "_PU_id.pdf","pdf");




   doPlot_PUid(PU_id_BB1_MC,PU_id_BB1_notmatchedMC, PU_id_BB1_DATA,0,false,ptbinLabel,txttag,SavePlots,false);
   doPlot_PUid(PU_id_BB2_MC,PU_id_BB2_notmatchedMC, PU_id_BB2_DATA,1,false,ptbinLabel,txttag,SavePlots,false);
   doPlot_PUid(PU_id_EC1_MC,PU_id_EC1_notmatchedMC, PU_id_EC1_DATA,2,false,ptbinLabel,txttag,SavePlots,false);
   doPlot_PUid(PU_id_EC2_MC,PU_id_EC2_notmatchedMC, PU_id_EC2_DATA,3,false,ptbinLabel,txttag,SavePlots,false);
   doPlot_PUid(PU_id_HF_MC,PU_id_HF_notmatchedMC, PU_id_HF_DATA,4,false,ptbinLabel,txttag,SavePlots,false);

   doPlot_PUid(PU_id_BB1_MC_3rd,PU_id_BB1_notmatchedMC_3rd, PU_id_BB1_DATA_3rd,0,true,ptbinLabel,txttag,SavePlots,false);
   doPlot_PUid(PU_id_BB2_MC_3rd,PU_id_BB2_notmatchedMC_3rd, PU_id_BB2_DATA_3rd,1,true,ptbinLabel,txttag,SavePlots,false);
   doPlot_PUid(PU_id_EC1_MC_3rd,PU_id_EC1_notmatchedMC_3rd, PU_id_EC1_DATA_3rd,2,true,ptbinLabel,txttag,SavePlots,false);
   doPlot_PUid(PU_id_EC2_MC_3rd,PU_id_EC2_notmatchedMC_3rd, PU_id_EC2_DATA_3rd,3,true,ptbinLabel,txttag,SavePlots,false);
   doPlot_PUid(PU_id_HF_MC_3rd,PU_id_HF_notmatchedMC_3rd, PU_id_HF_DATA_3rd,4,true,ptbinLabel,txttag,SavePlots,false);

   doPlot_PUid(PU_id_BB1_MCparton,PU_id_BB1_notmatchedMCparton, PU_id_BB1_DATA,0,false,ptbinLabel,txttag,SavePlots,true);
   doPlot_PUid(PU_id_BB2_MCparton,PU_id_BB2_notmatchedMCparton, PU_id_BB2_DATA,1,false,ptbinLabel,txttag,SavePlots,true);
   doPlot_PUid(PU_id_EC1_MCparton,PU_id_EC1_notmatchedMCparton, PU_id_EC1_DATA,2,false,ptbinLabel,txttag,SavePlots,true);
   doPlot_PUid(PU_id_EC2_MCparton,PU_id_EC2_notmatchedMCparton, PU_id_EC2_DATA,3,false,ptbinLabel,txttag,SavePlots,true);
   doPlot_PUid(PU_id_HF_MCparton,PU_id_HF_notmatchedMCparton, PU_id_HF_DATA,4,false,ptbinLabel,txttag,SavePlots,true);

   doPlot_PUid(PU_id_BB1_MCparton_3rd,PU_id_BB1_notmatchedMCparton_3rd, PU_id_BB1_DATA_3rd,0,true,ptbinLabel,txttag,SavePlots,true);
   doPlot_PUid(PU_id_BB2_MCparton_3rd,PU_id_BB2_notmatchedMCparton_3rd, PU_id_BB2_DATA_3rd,1,true,ptbinLabel,txttag,SavePlots,true);
   doPlot_PUid(PU_id_EC1_MCparton_3rd,PU_id_EC1_notmatchedMCparton_3rd, PU_id_EC1_DATA_3rd,2,true,ptbinLabel,txttag,SavePlots,true);
   doPlot_PUid(PU_id_EC2_MCparton_3rd,PU_id_EC2_notmatchedMCparton_3rd, PU_id_EC2_DATA_3rd,3,true,ptbinLabel,txttag,SavePlots,true);
   doPlot_PUid(PU_id_HF_MCparton_3rd,PU_id_HF_notmatchedMCparton_3rd, PU_id_HF_DATA_3rd,4,true,ptbinLabel,txttag,SavePlots,true);

   /*
   TCanvas* cPU_id_BB1_MC = new TCanvas();
   PU_id_BB1_MC->GetXaxis()->SetTitle("PU ID (|#eta|<1.3)");
   PU_id_BB1_MC->GetXaxis()->SetLabelSize(0.045);
   PU_id_BB1_MC->GetYaxis()->SetLabelSize(0.045);
   PU_id_BB1_MC->GetYaxis()->SetTitleSize(0.045);
   PU_id_BB1_MC->SetMarkerColor(kRed);
   PU_id_BB1_MC->SetMarkerStyle(20);
   PU_id_BB1_DATA->SetMarkerColor(kGreen);
   PU_id_BB1_DATA->SetMarkerStyle(22);
   PU_id_BB1_DATA_3rd->SetMarkerColor(kGreen+2);
   PU_id_BB1_DATA_3rd->SetMarkerStyle(25);

   PU_id_BB1_MC->SetTitle("");
   PU_id_BB1_notmatchedMC->SetMarkerColor(kBlack);
   PU_id_BB1_notmatchedMC->SetMarkerStyle(21);

   Double_t scale = 1/PU_id_BB1_MC->Integral();
   if(scale>0)
     PU_id_BB1_MC->Scale(scale);
   PU_id_BB1_MC->GetYaxis()->SetRangeUser(0.0,0.2);
   PU_id_BB1_MC->Draw();
   scale = 1/PU_id_BB1_DATA->Integral();
   PU_id_BB1_DATA->Print();
   if(scale>0){
     PU_id_BB1_DATA->Scale(scale);
     PU_id_BB1_DATA->Draw("same");
   }

   scale = 1/PU_id_BB1_notmatchedMC->Integral();
   if(scale>0){
     PU_id_BB1_notmatchedMC->Scale(scale);
     PU_id_BB1_notmatchedMC->Draw("same");
   }

   scale = 1/PU_id_BB1_DATA_3rd->Integral();
   if(scale>0){
     PU_id_BB1_DATA_3rd->Scale(scale);
     PU_id_BB1_DATA_3rd->Draw("same");
   }

   PU_id_BB1_DATA_3rd->Print();

   //   cPU_id_BB1_MC->SetLogy();
   //   cPU_id_BB1_MC->SetLogz();
   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   leg->Draw();
   cPU_id_BB1_MC->Print(SavePlots + "_PU_id_BB1.pdf","pdf");

   TCanvas* cPU_id_BB2_MC = new TCanvas();
   PU_id_BB2_MC->GetXaxis()->SetTitle("PU ID (1.3<|#eta|<2.1)");
   PU_id_BB2_MC->GetXaxis()->SetLabelSize(0.045);
   PU_id_BB2_MC->GetYaxis()->SetLabelSize(0.045);
   PU_id_BB2_MC->GetYaxis()->SetTitleSize(0.045);
   PU_id_BB2_MC->SetMarkerColor(kRed);
   PU_id_BB2_MC->SetMarkerStyle(20);
   PU_id_BB2_DATA->SetMarkerColor(kGreen);
   PU_id_BB2_DATA->SetMarkerStyle(22);
   PU_id_BB2_DATA_3rd->SetMarkerColor(kGreen+2);
   PU_id_BB2_DATA_3rd->SetMarkerStyle(25);

   PU_id_BB2_MC->SetTitle("");
   PU_id_BB2_notmatchedMC->SetMarkerColor(kBlack);
   PU_id_BB2_notmatchedMC->SetMarkerStyle(21);
   //   cPU_id_BB2_MC->SetLogy();
   //   cPU_id_BB2_MC->SetLogz();
   scale = 1/PU_id_BB2_MC->Integral();
   if(scale>0)
     PU_id_BB2_MC->Scale(scale);
   PU_id_BB2_MC->GetYaxis()->SetRangeUser(0.0,0.2);
   PU_id_BB2_MC->Draw();
   scale = 1/PU_id_BB2_notmatchedMC->Integral();
   if(scale>0){
     PU_id_BB2_notmatchedMC->Scale(scale);
     PU_id_BB2_notmatchedMC->Draw("same");
   }
   scale = 1/PU_id_BB2_DATA->Integral();
   if(scale>0){
     PU_id_BB2_DATA->Scale(scale);
     PU_id_BB2_DATA->Draw("same");
   }
   scale = 1/PU_id_BB2_DATA_3rd->Integral();
   if(scale>0){
     PU_id_BB2_DATA_3rd->Scale(scale);
     PU_id_BB2_DATA_3rd->Draw("same");
   }

   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   leg->Draw();
   cPU_id_BB2_MC->Print(SavePlots + "_PU_id_BB2.pdf","pdf");

   TCanvas* cPU_id_EC1_MC = new TCanvas();
   PU_id_EC1_MC->GetXaxis()->SetTitle("PU ID (2.1<|#eta|<2.5)");
   PU_id_EC1_MC->GetXaxis()->SetLabelSize(0.045);
   PU_id_EC1_MC->GetYaxis()->SetLabelSize(0.045);
   PU_id_EC1_MC->GetYaxis()->SetTitleSize(0.045);
   PU_id_EC1_MC->SetMarkerColor(kRed);
   PU_id_EC1_MC->SetMarkerStyle(20);
   PU_id_EC1_DATA->SetMarkerColor(kGreen);
   PU_id_EC1_DATA->SetMarkerStyle(22);
   PU_id_EC1_DATA_3rd->SetMarkerColor(kGreen+2);
   PU_id_EC1_DATA_3rd->SetMarkerStyle(25);

   PU_id_EC1_MC->SetTitle("");
   PU_id_EC1_notmatchedMC->SetMarkerColor(kBlack);
   PU_id_EC1_notmatchedMC->SetMarkerStyle(21);

   //   cPU_id_EC1_MC->SetLogy();
   //   cPU_id_EC1_MC->SetLogz();
   scale = 1/PU_id_EC1_MC->Integral();
   if(scale>0)
     PU_id_EC1_MC->Scale(scale);
   PU_id_EC1_MC->GetYaxis()->SetRangeUser(0.0,0.2);
   PU_id_EC1_MC->Draw();
   scale = 1/PU_id_EC1_notmatchedMC->Integral();
   if(scale>0){
     PU_id_EC1_notmatchedMC->Scale(scale);
     PU_id_EC1_notmatchedMC->Draw("same");
   }
   scale = 1/PU_id_EC1_DATA->Integral();
   if(scale>0){
     PU_id_EC1_DATA->Scale(scale);
     PU_id_EC1_DATA->Draw("same");
   }
   scale = 1/PU_id_EC1_DATA_3rd->Integral();
   if(scale>0){
     PU_id_EC1_DATA_3rd->Scale(scale);
     PU_id_EC1_DATA_3rd->Draw("same");
   }

   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   leg->Draw();
   cPU_id_EC1_MC->Print(SavePlots + "_PU_id_EC1.pdf","pdf");


   TCanvas* cPU_id_EC2_MC = new TCanvas();
   PU_id_EC2_MC->GetXaxis()->SetTitle("PU ID (2.5<|#eta|<3.0)");
   PU_id_EC2_MC->GetXaxis()->SetLabelSize(0.045);
   PU_id_EC2_MC->GetYaxis()->SetLabelSize(0.045);
   PU_id_EC2_MC->GetYaxis()->SetTitleSize(0.045);
   PU_id_EC2_MC->SetMarkerColor(kRed);
   PU_id_EC2_MC->SetMarkerStyle(20);
   PU_id_EC2_MC->SetTitle("");
   PU_id_EC2_DATA->SetMarkerColor(kGreen);
   PU_id_EC2_DATA->SetMarkerStyle(22);
   PU_id_EC2_DATA_3rd->SetMarkerColor(kGreen+2);
   PU_id_EC2_DATA_3rd->SetMarkerStyle(25);

   PU_id_EC2_notmatchedMC->SetMarkerColor(kBlack);
   PU_id_EC2_notmatchedMC->SetMarkerStyle(21);
   //   cPU_id_EC2_MC->SetLogy();
   //   cPU_id_EC2_MC->SetLogz();
   scale = 1/PU_id_EC2_MC->Integral();
   if(scale>0)
     PU_id_EC2_MC->Scale(scale);
   PU_id_EC2_MC->GetYaxis()->SetRangeUser(0.0,0.2);
   PU_id_EC2_MC->Draw();
   scale = 1/PU_id_EC2_DATA->Integral();
   if(scale>0){
     PU_id_EC2_DATA->Scale(scale);
     PU_id_EC2_DATA->Draw("same");
   }
   scale = 1/PU_id_EC2_notmatchedMC->Integral();
   if(scale>0){
     PU_id_EC2_notmatchedMC->Scale(scale);
     PU_id_EC2_notmatchedMC->Draw("same");
   }
   scale = 1/PU_id_EC2_DATA_3rd->Integral();
   if(scale>0){
     PU_id_EC2_DATA_3rd->Scale(scale);
     PU_id_EC2_DATA_3rd->Draw("same");
   }

   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   leg->Draw();
   cPU_id_EC2_MC->Print(SavePlots + "_PU_id_EC2.pdf","pdf");

   TCanvas* cPU_id_HF_MC = new TCanvas();
   PU_id_HF_MC->GetXaxis()->SetTitle("PU ID (|#eta|>3.0)");
   PU_id_HF_MC->GetXaxis()->SetLabelSize(0.045);
   PU_id_HF_MC->GetYaxis()->SetLabelSize(0.045);
   PU_id_HF_MC->GetYaxis()->SetTitleSize(0.045);
   PU_id_HF_MC->SetMarkerColor(kRed);
   PU_id_HF_MC->SetMarkerStyle(20);
   PU_id_HF_MC->SetTitle("");
   PU_id_HF_DATA->SetMarkerColor(kGreen);
   PU_id_HF_DATA->SetMarkerStyle(22);
   PU_id_HF_DATA_3rd->SetMarkerColor(kGreen+2);
   PU_id_HF_DATA_3rd->SetMarkerStyle(25);

   PU_id_HF_notmatchedMC->SetMarkerColor(kBlack);
   PU_id_HF_notmatchedMC->SetMarkerStyle(21);
   //   cPU_id_HF_MC->SetLogy();
   //   cPU_id_HF_MC->SetLogz();
   scale = 1/PU_id_HF_MC->Integral();
   if(scale>0)
     PU_id_HF_MC->Scale(scale);
   PU_id_HF_MC->GetYaxis()->SetRangeUser(0.0,0.2);
   PU_id_HF_MC->Draw();
   scale = 1/PU_id_HF_DATA->Integral();
   if(scale>0){
     PU_id_HF_DATA->Scale(scale);
     PU_id_HF_DATA->Draw("same");
   }

   scale = 1/PU_id_HF_notmatchedMC->Integral();
   if(scale>0){
     PU_id_HF_notmatchedMC->Scale(scale);
     PU_id_HF_notmatchedMC->Draw("same");
   }
   scale = 1/PU_id_HF_DATA_3rd->Integral();
   if(scale>0){
     PU_id_HF_DATA_3rd->Scale(scale);
     PU_id_HF_DATA_3rd->Draw("same");
   }

   tex->DrawLatex(0.15,0.95,txttag);
   tex->DrawLatex(0.68,0.95,ptbinLabel);
   leg->Draw();
   cPU_id_HF_MC->Print(SavePlots + "_PU_id_HF.pdf","pdf");

   */
   delete cJet_eta_pt_MC;
   delete cJet_eta_PUid_MC;
   delete cJet_eta_chEM_MC;
   delete cJet_eta_neuEM_MC;
   delete cJet_eta_chHAD_MC;
   delete cJet_eta_neuHAD_MC;

   delete cnPu_MC;
   delete cN_PV_MC;
   delete cPU_id_MC;
  };
}
