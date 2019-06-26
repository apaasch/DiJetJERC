// Read histograms filled in JECAnalysisRecoGenMatchedHists (by BaconJet/AnalysisModule module) for hadron composition of GEN probe jet
// created on 19.12.2018 by A. Karavdina

#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TH3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <THStack.h>
#include <TLatex.h>
using namespace std;

void CorrectionObject::ControlPlotsRecoGenCompare(bool forEveryPtBin){
  cout << "--------------- Starting ControlPlotsRecoGenCompare() ---------------" << endl;
  bool forEveryPtBin_ = forEveryPtBin;
  if(forEveryPtBin_){
    cout<<"ControlPlotsRecoGenCompare are gonne be made for each pt bin seperatly"<<endl;
  }
  TString txttag=CorrectionObject::_generator_tag; 
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.06,"x");  
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"x");  
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleYOffset(0.9);
  gStyle->SetTitleXOffset(0.85);
  gStyle->SetPalette(55);
  gROOT->ForceStyle();



  TString dirNames[15] = {"MatchedRecoGen_all", "MatchedRecoGen_10_15", "MatchedRecoGen_15_20","MatchedRecoGen_20_25", "MatchedRecoGen_25_30",
			  "MatchedRecoGen_30_40","MatchedRecoGen_40_50",
			  "MatchedRecoGen_50_70","MatchedRecoGen_70_120","MatchedRecoGen_120_250",
			  "MatchedRecoGen_250_400","MatchedRecoGen_400_650","MatchedRecoGen_650_1000",
			  "MatchedRecoGen_1000_1500","MatchedRecoGen_1500_2000"};

  TString ptbinsLabels[15] = {"10 #leqp^{GEN}_{T}#leq #infty",
			      "10 #leqp^{GEN}_{T}#leq 15",
			      "15 #leqp^{GEN}_{T}#leq 20",
			      "20 #leqp^{GEN}_{T}#leq 25",
			      "25 #leqp^{GEN}_{T}#leq 30",
			      "30 #leqp^{GEN}_{T}#leq 40",
			      "40 #leqp^{GEN}_{T}#leq 50", "50 #leqp^{GEN}_{T}#leq 70", "70 #leqp^{GEN}_{T}#leq 120", 
			      "120 #leqp^{GEN}_{T}#leq 250", "250 #leqp^{GEN}_{T}#leq 400", "400 #leqp^{GEN}_{T}#leq 650",
			      "650 #leqp^{GEN}_{T}#leq 1000", "1000 #leqp^{GEN}_{T}#leq 1500",
			      "1500#leqp^{GEN}_{T}#leq 2000"};
  int loop_inter = 1;
  if(forEveryPtBin) loop_inter=12;
  //  if(forEveryPtBin) loop_inter=1;//TEST
  
  for(int i = 0; i<loop_inter ; i++){
  
   TString dirName = (TString)dirNames[i];  
   TString ptbinLabel = (TString)ptbinsLabels[i];
   cout<<endl;  
   cout<<dirName<<endl;
      
  TString MCtitle = "MC";
  TString SavePlots = CorrectionObject::_outpath + "plots/control/ControlPlotsRecoGenCompare/CP_RecoGenMatched_" + dirName + "_" + CorrectionObject::_generator_tag;
  CorrectionObject::make_path(std::string((_outpath + "plots/").Data()));
  CorrectionObject::make_path(std::string((_outpath + "plots/control/").Data()));
  CorrectionObject::make_path(std::string((_outpath + "plots/control/ControlPlotsRecoGenCompare/").Data()));

  cout<<"++++++++++++++++ Collect all histograms ++++++++++++++++\n";
  TH2F* Response_eta_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/Response_eta");
  TH2F* Hadrons_energy_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/Hadrons_energy");
  TH2F* Hadrons_energy_rel_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/Hadrons_energy_rel");
  TH2F* Hadrons_genjeteta_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/Hadrons_genjeteta");
  TH2F* Hadrons_energy_rel_event_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/Hadrons_energy_rel_event");
  TH2F* Hadrons_count_event_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/Hadrons_count_event");
  TH1F* Hadrons_total_energy_rel_event_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/Hadrons_total_energy_rel_event");
  TH1F* Hadrons_total_count_event_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/Hadrons_total_count_event");
  TH3D* Hadrons_energy_rel_event_eta_MC = (TH3D*)CorrectionObject::_MCFile->Get(dirName+"/Hadrons_energy_rel_event_eta");
  TH3D* PF_frac_event_eta_MC = (TH3D*)CorrectionObject::_MCFile->Get(dirName+"/PF_frac_event_eta");
  TH3D* PF_to_HAD_event_eta_MC = (TH3D*)CorrectionObject::_MCFile->Get(dirName+"/PF_to_HAD_event_eta");
  TH3D* N_PF_to_HAD_event_eta_MC = (TH3D*)CorrectionObject::_MCFile->Get(dirName+"/N_PF_to_HAD_event_eta");
  cout<<" +++ read DATA hists +++ "<<endl;
  TH2F* Response_eta_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/Response_eta");
  TH2F* Hadrons_energy_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/Hadrons_energy");
  TH2F* Hadrons_energy_rel_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/Hadrons_energy_rel");
  TH2F* Hadrons_genjeteta_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/Hadrons_genjeteta");
  TH2F* Hadrons_energy_rel_event_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/Hadrons_energy_rel_event");
  TH2F* Hadrons_count_event_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/Hadrons_count_event");
  TH1F* Hadrons_total_energy_rel_event_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/Hadrons_total_energy_rel_event");
  TH1F* Hadrons_total_count_event_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/Hadrons_total_count_event");
  TH3D* Hadrons_energy_rel_event_eta_DATA = (TH3D*)CorrectionObject::_DATAFile->Get(dirName+"/Hadrons_energy_rel_event_eta");
  TH3D* PF_frac_event_eta_DATA = (TH3D*)CorrectionObject::_DATAFile->Get(dirName+"/PF_frac_event_eta");
  TH3D* PF_to_HAD_event_eta_DATA = (TH3D*)CorrectionObject::_DATAFile->Get(dirName+"/PF_to_HAD_event_eta");
  TH3D* N_PF_to_HAD_event_eta_DATA = (TH3D*)CorrectionObject::_DATAFile->Get(dirName+"/N_PF_to_HAD_event_eta");


  PF_frac_event_eta_MC->SetName("PF_frac_event_eta_MC");
  PF_frac_event_eta_DATA->SetName("PF_frac_event_eta_DATA");
  Response_eta_MC->SetName("Response_eta_MC");
  Response_eta_DATA->SetName("Response_eta_DATA");
  Hadrons_energy_rel_MC->SetName("Hadrons_energy_rel_MC");
  Hadrons_genjeteta_MC->SetName("Hadrons_genjeteta_MC");
  Hadrons_energy_rel_event_MC->SetName("Hadrons_energy_rel_event_MC");
  Hadrons_count_event_MC->SetName("Hadrons_count_event_MC");
  Hadrons_total_energy_rel_event_MC->SetName("Hadrons_total_energy_rel_event_MC");
  Hadrons_total_count_event_MC->SetName("Hadrons_total_count_event_MC");
  Hadrons_energy_rel_event_eta_MC->SetName("Hadrons_energy_rel_event_eta_MC");
  PF_to_HAD_event_eta_MC->SetName("PF_to_HAD_event_eta_MC");
  N_PF_to_HAD_event_eta_MC->SetName("N_PF_to_HAD_event_eta_MC");
  Hadrons_energy_rel_DATA->SetName("Hadrons_energy_rel_DATA");
  Hadrons_genjeteta_DATA->SetName("Hadrons_genjeteta_DATA");
  Hadrons_energy_rel_event_DATA->SetName("Hadrons_energy_rel_event_DATA");
  Hadrons_count_event_DATA->SetName("Hadrons_count_event_DATA");
  Hadrons_total_energy_rel_event_DATA->SetName("Hadrons_total_energy_rel_event_DATA");
  Hadrons_total_count_event_DATA->SetName("Hadrons_total_count_event_DATA");
  Hadrons_energy_rel_event_eta_DATA->SetName("Hadrons_energy_rel_event_eta_DATA");
  PF_to_HAD_event_eta_DATA->SetName("PF_to_HAD_event_eta_DATA");
  N_PF_to_HAD_event_eta_DATA->SetName("N_PF_to_HAD_event_eta_DATA");



 
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);  

  cout<<"+++++++++++++++++++++ Plot extremely nice histograms ++++++++++++++++++++++++++++\n";
  TCanvas* cResponse_eta_MC = new TCanvas();
  Response_eta_MC->GetYaxis()->SetTitle("p^{RECO}_{T}/p^{GEN}_{T}, GeV");
  //  Response_eta_MC->GetXaxis()->SetTitle("Hadrons");
  Response_eta_MC->GetXaxis()->SetLabelSize(0.07);
  Response_eta_MC->GetYaxis()->SetLabelSize(0.045);
  Response_eta_MC->GetYaxis()->SetTitleSize(0.045);
  Response_eta_MC->SetTitle("");
  Response_eta_MC->Draw("colz");
  cResponse_eta_MC->SetLogy();
  cResponse_eta_MC->SetLogz();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cResponse_eta_MC->Print(SavePlots + "_Response_eta.pdf","pdf");

  Response_eta_MC->FitSlicesY();
  Response_eta_DATA->FitSlicesY();
  TCanvas* cResponse_eta_MC_mean = new TCanvas();
  TH1D* Response_eta_MC_mean = (TH1D*)gDirectory->Get("Response_eta_MC_1");
  Response_eta_MC_mean->SetMarkerColor(2);
  Response_eta_MC_mean->SetMarkerStyle(20);
  Response_eta_MC_mean->GetXaxis()->SetLabelSize(0.045);
  Response_eta_MC_mean->GetYaxis()->SetLabelSize(0.045);
  Response_eta_MC_mean->GetYaxis()->SetTitleSize(0.045);
  Response_eta_MC_mean->GetYaxis()->SetRangeUser(0.7,1.3);
  Response_eta_MC_mean->GetYaxis()->SetTitle("<p^{RECO}_{T}/p^{GEN}_{T}>");
  Response_eta_MC_mean->GetXaxis()->SetTitle("jet #eta");
  Response_eta_MC_mean->SetTitle("");
  Response_eta_MC_mean->Draw();
  TH1D* Response_eta_DATA_mean = (TH1D*)gDirectory->Get("Response_eta_DATA_1");
  Response_eta_DATA_mean->SetMarkerColor(1);
  Response_eta_DATA_mean->SetMarkerStyle(21);
  Response_eta_DATA_mean->Draw("same");
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  TLegend* legRes =  new TLegend(0.19,0.70,0.35,0.87);
  legRes->AddEntry(Response_eta_MC_mean,"old","p");
  legRes->AddEntry(Response_eta_DATA_mean,"new","p");
  legRes->Draw();
  cResponse_eta_MC_mean->Print(SavePlots + "_Response_eta_mean.pdf","pdf");

  TCanvas* cResponse_eta_MC_std = new TCanvas();
  TH1D* Response_eta_MC_std = (TH1D*)gDirectory->Get("Response_eta_MC_2");
  Response_eta_MC_std->Divide(Response_eta_MC_mean);
  Response_eta_MC_std->SetMarkerColor(2);
  Response_eta_MC_std->SetMarkerStyle(20);
  Response_eta_MC_std->GetXaxis()->SetLabelSize(0.045);
  Response_eta_MC_std->GetYaxis()->SetLabelSize(0.045);
  Response_eta_MC_std->GetYaxis()->SetTitleSize(0.045);
  Response_eta_MC_std->GetYaxis()->SetRangeUser(0.0,0.35);
  Response_eta_MC_std->GetYaxis()->SetTitle("#sigma(p^{RECO}_{T}/p^{GEN}_{T})/<p^{RECO}_{T}/p^{GEN}_{T}>");
  Response_eta_MC_std->GetXaxis()->SetTitle("jet #eta");
  Response_eta_MC_std->SetTitle("");
  Response_eta_MC_std->Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  TH1D* Response_eta_DATA_std = (TH1D*)gDirectory->Get("Response_eta_DATA_2");
  Response_eta_DATA_std->Divide(Response_eta_DATA_mean);
  Response_eta_DATA_std->SetMarkerColor(1);
  Response_eta_DATA_std->SetMarkerStyle(21);
  Response_eta_DATA_std->Draw("same");
  legRes->Draw();
  cResponse_eta_MC_std->Print(SavePlots + "_Response_eta_std.pdf","pdf");


 bool use_mc=true;
  if (Hadrons_energy_MC->GetEntries()==0 && Hadrons_genjeteta_MC->GetEntries()==0){
    cout<<"Hadrons are not stored in MC, skip"<<endl;
    use_mc=false;
  }
  if(!use_mc) continue;

  TCanvas* cHadrons_energy_MC = new TCanvas();
  Hadrons_energy_MC->GetYaxis()->SetTitle("hadron E, GeV");
  //  Hadrons_energy_MC->GetXaxis()->SetTitle("Hadrons");
  Hadrons_energy_MC->GetXaxis()->SetLabelSize(0.07);
  Hadrons_energy_MC->GetYaxis()->SetLabelSize(0.045);
  Hadrons_energy_MC->GetYaxis()->SetTitleSize(0.045);
  Hadrons_energy_MC->SetTitle("");
  Hadrons_energy_MC->Draw("colz");
  cHadrons_energy_MC->SetLogy();
  cHadrons_energy_MC->SetLogz();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cHadrons_energy_MC->Print(SavePlots + "_Hadrons_energy.pdf","pdf");

  TCanvas* cHadrons_energy_rel_MC = new TCanvas();
  Hadrons_energy_rel_MC->GetYaxis()->SetTitle("single hadron E/GEN jet E");
  //  Hadrons_energy_rel_MC->GetXaxis()->SetTitle("Hadrons");
  Hadrons_energy_rel_MC->GetXaxis()->SetLabelSize(0.07);
  Hadrons_energy_rel_MC->GetYaxis()->SetLabelSize(0.045);
  Hadrons_energy_rel_MC->GetYaxis()->SetTitleSize(0.045);
  Hadrons_energy_rel_MC->SetTitle("");
  // Double_t scale = 1/Hadrons_energy_rel_MC->Integral();
  // Hadrons_energy_rel_MC->Scale(scale);
  // Hadrons_energy_rel_MC->GetZaxis()->SetRangeUser(0,1.0);
  Hadrons_energy_rel_MC->Draw("colz");
  cHadrons_energy_rel_MC->SetLogy();
  // cHadrons_energy_rel_MC->SetLogz();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cHadrons_energy_rel_MC->Print(SavePlots + "_Hadrons_energy_rel.pdf","pdf");

  TCanvas* cHadrons_genjeteta_MC = new TCanvas();
  //  Hadrons_genjeteta_MC->GetXaxis()->SetTitle("Hadrons");
  Hadrons_genjeteta_MC->GetXaxis()->SetLabelSize(0.07);
  Hadrons_genjeteta_MC->GetYaxis()->SetLabelSize(0.045);
  Hadrons_genjeteta_MC->GetYaxis()->SetTitleSize(0.045);
  Hadrons_genjeteta_MC->SetTitle("");
  Hadrons_genjeteta_MC->Draw("colz");
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cHadrons_genjeteta_MC->Print(SavePlots + "_Hadrons_genjeteta.pdf","pdf");

  TCanvas* cHadrons_energy_rel_event_MC = new TCanvas();
  Hadrons_energy_rel_event_MC->GetYaxis()->SetTitle("#sum (hadron E)/GEN jet E");
  //  Hadrons_energy_rel_event_MC->GetXaxis()->SetTitle("Hadrons");
  Hadrons_energy_rel_event_MC->GetXaxis()->SetLabelSize(0.07);
  Hadrons_energy_rel_event_MC->GetYaxis()->SetLabelSize(0.045);
  Hadrons_energy_rel_event_MC->GetYaxis()->SetTitleSize(0.045);
  Hadrons_energy_rel_event_MC->SetTitle("");
  Hadrons_energy_rel_event_MC->Draw("colz");
  cHadrons_energy_rel_event_MC->SetLogy();
  cHadrons_energy_rel_event_MC->SetLogz();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cHadrons_energy_rel_event_MC->Print(SavePlots + "_Hadrons_energy_rel_event.pdf","pdf");

  TCanvas* cHadrons_count_event_MC = new TCanvas();
  Hadrons_count_event_MC->GetYaxis()->SetTitle("Number of hadrons per jet");
  //  Hadrons_count_event_MC->SetTitleYOffset(0.9);
  Hadrons_count_event_MC->GetXaxis()->SetLabelSize(0.07);
  Hadrons_count_event_MC->GetYaxis()->SetLabelSize(0.045);
  Hadrons_count_event_MC->GetYaxis()->SetTitleSize(0.045);
  Hadrons_count_event_MC->SetTitle("");
  Hadrons_count_event_MC->Draw("colz");
  // cHadrons_count_event_MC->SetLogy();
  // cHadrons_count_event_MC->SetLogz();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);

  cHadrons_count_event_MC->Print(SavePlots + "_Hadrons_count_event.pdf","pdf");


  TCanvas* cHadrons_total_energy_rel_event_MC = new TCanvas();
  Hadrons_total_energy_rel_event_MC->GetXaxis()->SetLabelSize(0.045);
  Hadrons_total_energy_rel_event_MC->GetYaxis()->SetLabelSize(0.045);
  Hadrons_total_energy_rel_event_MC->GetXaxis()->SetRangeUser(0.7,1.3);
  Hadrons_total_energy_rel_event_MC->GetXaxis()->SetTitleSize(0.045);
  Hadrons_total_energy_rel_event_MC->SetMarkerStyle(20);
  Hadrons_total_energy_rel_event_MC->SetMarkerColor(kRed);
  Hadrons_total_energy_rel_event_MC->SetTitle("");
  Hadrons_total_energy_rel_event_MC->Draw("P");
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cHadrons_total_energy_rel_event_MC->Print(SavePlots + "_Hadrons_total_energy_rel_event.pdf","pdf");

  TCanvas* cHadrons_total_count_event_MC = new TCanvas();
  Hadrons_total_count_event_MC->GetXaxis()->SetLabelSize(0.045);
  Hadrons_total_count_event_MC->GetYaxis()->SetLabelSize(0.045);
  Hadrons_total_count_event_MC->GetXaxis()->SetTitleSize(0.045);
  Hadrons_total_count_event_MC->SetMarkerStyle(20);
  Hadrons_total_count_event_MC->SetMarkerColor(kRed);
  Hadrons_total_count_event_MC->SetTitle("");
  Hadrons_total_count_event_MC->Draw("P");
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);

  TString textMean ="Mean=";textMean+=Form("%.1f",Hadrons_total_count_event_MC->GetMean());
  tex->DrawLatex(0.70,0.85,textMean);
  cHadrons_total_count_event_MC->Print(SavePlots + "_Hadrons_total_count_event.pdf","pdf");



  double hadron_energy_fractions[n_eta_common-1][8]; //pt fraction ratio in eta bins for each hadron
  double hadron_energy_fractions_median[n_eta_common-1][8]; //pt fraction ratio in eta bins for each hadron
  double pf_had_energy_fractions[n_eta_common-1][2]; //pf fraction/had fraction ratio in eta bins for ch/neut. hadrons
  double n_pf_had_energy_fractions[n_eta_common-1][2]; //pf fraction/had fraction ratio in eta bins for ch/neut. hadrons
  double pf_energy_fractions[n_eta_common-1][5]; //pf in eta bins for each hadron

  double hadron_energy_fractions_DATA[n_eta_common-1][8]; //pt fraction ratio in eta bins for each hadron
  double hadron_energy_fractions_median_DATA[n_eta_common-1][8]; //pt fraction ratio in eta bins for each hadron
  double pf_had_energy_fractions_DATA[n_eta_common-1][2]; //pf fraction/had fraction ratio in eta bins for ch/neut. hadrons
  double n_pf_had_energy_fractions_DATA[n_eta_common-1][2]; //pf fraction/had fraction ratio in eta bins for ch/neut. hadrons
  double pf_energy_fractions_DATA[n_eta_common-1][5]; //pf in eta bins for each hadron

  TH1D *ptfrac[8]; 
  THStack *hs_energyfrac = new THStack("hadron_energy_fractions_hs","; GEN jet #eta;<hadron E/GEN jet E> (mean)");
  TLegend* leg2 =  new TLegend(0.89,0.42,0.99,0.99);
  TString name_ext[8];
  for(int k=0;k<8;k++){
    name_ext[k] = Hadrons_energy_rel_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
    TString name = "ptfrac_"+name_ext[k]; 
    name +="_GLbinN"; name +=i;
    //    cout<<name.Data()<<endl;
    ptfrac[k] = new TH1D(name,"",n_eta-1, eta_bins);
    int markerCol=k+1;
    if(markerCol>4) markerCol=k+2;
    ptfrac[k]->SetLineColor(markerCol);
    ptfrac[k]->SetFillColor(markerCol);
    leg2 ->AddEntry(ptfrac[k],name_ext[k],"f");
    //    ptfrac[k]->SetLineWidth(2);
  }

 
  TH1D *ptfrac_median[8]; 
  THStack *hs_energyfrac_median = new THStack("hadron_energy_fractions_hs_median","; GEN jet #eta;<hadron E/GEN jet E> (median)");
  TLegend* leg2_median =  new TLegend(0.89,0.42,0.99,0.99);
  for(int k=0;k<8;k++){
    TString name = "ptfrac_median_"+name_ext[k]; 
    name +="_GLbinN"; name +=i;
    //    cout<<name.Data()<<endl;
    ptfrac_median[k] = new TH1D(name,"",n_eta-1, eta_bins);
    int markerCol=k+1;
    if(markerCol>4) markerCol=k+2;
    ptfrac_median[k]->SetLineColor(markerCol);
    ptfrac_median[k]->SetFillColor(markerCol);
    leg2_median ->AddEntry(ptfrac_median[k],name_ext[k],"f");
  }


  //PF to HAD fractions
  TString name_ext2[2];
  TH1D *pfhadfrac[2]; 
  THStack *hs_pfhadfrac = new THStack("pf_to_hadron_energy_fractions_hs_MC","; RECO jet #eta;<PF fraction/HAD fraction>");
  TH1D *pfhadfrac_DATA[2]; 
  THStack *hs_pfhadfrac_DATA = new THStack("pf_to_hadron_energy_fractions_hs_DATA","; RECO jet #eta;<PF fraction/HAD fraction>");

  TLegend* leg3 =  new TLegend(0.40,0.68,0.88,0.88);
  for(int k=0;k<2;k++){
    name_ext2[k] = PF_to_HAD_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
    TString name = "pffrac_to_hadfrac_"+name_ext2[k];
    name +="_GLbinN_MC"; name +=i;  
    pfhadfrac[k] = new TH1D(name,"",n_eta-1, eta_bins);
    name = "pffrac_to_hadfrac_"+name_ext2[k];
    name +="_GLbinN_DATA"; name +=i;  
    pfhadfrac_DATA[k] = new TH1D(name,"",n_eta-1, eta_bins);
    int markerCol=k+1;
    if(markerCol>4) markerCol=k+2;
    pfhadfrac[k]->SetLineColor(markerCol);
    pfhadfrac[k]->SetFillColorAlpha(markerCol,0.8-0.4*k);
    pfhadfrac_DATA[k]->SetMarkerColor(10*markerCol);
    pfhadfrac_DATA[k]->SetMarkerStyle(27);
    for(int j=0; j<n_eta_common-1; j++){
      pfhadfrac[k]->SetBinContent(j,0.0);
      pfhadfrac_DATA[k]->SetBinContent(j,0.0);
    }
    leg3 ->AddEntry(pfhadfrac[k],name_ext2[k],"f");
  }
  TLegend* legMCDATA = new TLegend(0.28,0.89,0.37,0.98);
  legMCDATA->AddEntry(pfhadfrac[0],"Old","f");
  legMCDATA->AddEntry(pfhadfrac_DATA[0],"New","p");

  //N PF to HAD constituence
  TString n_name_ext2[2];
  TH1D *n_pfhadfrac[2]; 
  THStack *hs_n_pfhadfrac = new THStack("n_pf_to_hadron_energy_fractions_hs_MC","; RECO jet #eta;<N PFs/N hadrons>");
  TH1D *n_pfhadfrac_DATA[2]; 
  THStack *hs_n_pfhadfrac_DATA = new THStack("n_pf_to_hadron_energy_fractions_hs_DATA","; RECO jet #eta;<N PFs/N hadrons>");

  TLegend* leg5 =  new TLegend(0.40,0.68,0.88,0.88);
  for(int k=0;k<2;k++){
    n_name_ext2[k] = N_PF_to_HAD_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
    TString name = "n_pffrac_to_hadfrac_"+name_ext2[k];
    name +="_GLbinN_MC"; name +=i;  
    n_pfhadfrac[k] = new TH1D(name,"",n_eta-1, eta_bins);
    name = "n_pffrac_to_hadfrac_"+name_ext2[k];
    name +="_GLbinN_DATA"; name +=i;  
    n_pfhadfrac_DATA[k] = new TH1D(name,"",n_eta-1, eta_bins);

    int markerCol=k+1;
    if(markerCol>4) markerCol=k+2;
    n_pfhadfrac[k]->SetLineColor(markerCol);
    n_pfhadfrac[k]->SetFillColorAlpha(markerCol,0.8-0.4*k);
    n_pfhadfrac_DATA[k]->SetMarkerColor(10*markerCol);
    n_pfhadfrac_DATA[k]->SetMarkerStyle(27);
    for(int j=0; j<n_eta_common-1; j++){
      n_pfhadfrac[k]->SetBinContent(j,0.0);
      n_pfhadfrac_DATA[k]->SetBinContent(j,0.0);
    }
    leg5 ->AddEntry(n_pfhadfrac[k],name_ext2[k],"f");
  }

  //PF  fractions
  TString name_ext3[5];
  TH1D *pffrac[5]; 
  THStack *hs_pffrac = new THStack("pf_energy_fractions_hs_MC","; RECO jet #eta;<PF fraction>");
  TH1D *pffrac_DATA[5]; 
  THStack *hs_pffrac_DATA = new THStack("pf_energy_fractions_hs_DATA","; RECO jet #eta;<PF fraction>");

  TLegend* leg4 =  new TLegend(0.85,0.30,0.99,0.99);
  for(int k=0;k<5;k++){
    name_ext3[k] = PF_frac_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
    TString name = "pffrac_"+name_ext3[k]; 
    name +="_GLbinN_MC"; name +=i;
    pffrac[k] = new TH1D(name,"",n_eta-1, eta_bins);
    name = "pffrac_"+name_ext3[k]; 
    name +="_GLbinN_DATA"; name +=i;
    pffrac_DATA[k] = new TH1D(name,"",n_eta-1, eta_bins);

    int markerCol=k+1;
    if(markerCol>4) markerCol=k+2;
    pffrac[k]->SetLineColor(markerCol);
    pffrac[k]->SetFillColorAlpha(markerCol,0.8-0.1*k);
    //pffrac_DATA[k]->SetMarkerColor(10*markerCol);
    pffrac_DATA[k]->SetMarkerColor(1);
    pffrac_DATA[k]->SetMarkerStyle(27);
    leg4 ->AddEntry(pffrac[k],name_ext3[k],"f");
  }

  double tot_frac[n_eta_common];   double tot_frac_mean[n_eta_common];  double tot_frac_median[n_eta_common];  double tot_frac_pf[n_eta_common];
  double tot_frac_DATA[n_eta_common];   double tot_frac_mean_DATA[n_eta_common];  double tot_frac_median_DATA[n_eta_common];  double tot_frac_pf_DATA[n_eta_common];

  //  double tot_frac_pf_had[n_eta_common];
  for(int j=0; j<n_eta_common-1; j++){

    Hadrons_energy_rel_event_eta_MC->GetXaxis()->SetRange(j,j+1);
    Hadrons_energy_rel_event_eta_DATA->GetXaxis()->SetRange(j,j+1);

    PF_frac_event_eta_MC->GetXaxis()->SetRange(j,j+1);
    PF_to_HAD_event_eta_MC->GetXaxis()->SetRange(j,j+1);
    N_PF_to_HAD_event_eta_MC->GetXaxis()->SetRange(j,j+1);
    PF_frac_event_eta_DATA->GetXaxis()->SetRange(j,j+1);
    PF_to_HAD_event_eta_DATA->GetXaxis()->SetRange(j,j+1);
    N_PF_to_HAD_event_eta_DATA->GetXaxis()->SetRange(j,j+1);

    TH2D *Hadrons_energy_rel_event_eta_MC_pyz = (TH2D*)Hadrons_energy_rel_event_eta_MC->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *Hadrons_energy_rel_event_eta_DATA_pyz = (TH2D*)Hadrons_energy_rel_event_eta_DATA->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *PF_frac_event_eta_MC_pyz = (TH2D*)PF_frac_event_eta_MC->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *PF_frac_event_eta_DATA_pyz = (TH2D*)PF_frac_event_eta_DATA->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *PF_to_HAD_event_eta_MC_pyz = (TH2D*)PF_to_HAD_event_eta_MC->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *N_PF_to_HAD_event_eta_MC_pyz = (TH2D*)N_PF_to_HAD_event_eta_MC->Project3D("yz");//fraction vs hadron in 1 eta bin

    TH2D *PF_to_HAD_event_eta_DATA_pyz = (TH2D*)PF_to_HAD_event_eta_DATA->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *N_PF_to_HAD_event_eta_DATA_pyz = (TH2D*)N_PF_to_HAD_event_eta_DATA->Project3D("yz");//fraction vs hadron in 1 eta bin

    tot_frac[j] = 0; tot_frac_mean[j] = 0;
    for(int k=0;k<8;k++){
      //      TString tmpName ="GLbinN"; tmpName +=i;tmpName +="_Hadron";tmpName +=name_ext[k];
      TString tmpName ="GLbinN_MC_"; tmpName +=i;tmpName +="_Hadron";tmpName +=k;
      TH1D *Hadrons_energy_rel_event_eta_MC_pyz_1D = (TH1D*)Hadrons_energy_rel_event_eta_MC_pyz->ProjectionY(tmpName,k+1,k+1);
      tmpName ="GLbinN_DATA_"; tmpName +=i;tmpName +="_Hadron";tmpName +=k;
      TH1D *Hadrons_energy_rel_event_eta_DATA_pyz_1D = (TH1D*)Hadrons_energy_rel_event_eta_DATA_pyz->ProjectionY(tmpName,k+1,k+1);


      hadron_energy_fractions[j][k] = 0.0;
      hadron_energy_fractions_median[j][k] = 0.0;
      //      int netn = Hadrons_energy_rel_event_eta_MC_pyz->GetEntries();
      //      if(netn<200) continue;
     
      // // //Mean
      hadron_energy_fractions[j][k] = Hadrons_energy_rel_event_eta_MC_pyz_1D->GetMean();//<fraction> for hadron   
      // //Median
      // Double_t x, q;
      // q = 0.5; // 0.5 for "median"
      // Hadrons_energy_rel_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      // Hadrons_energy_rel_event_eta_MC_pyz_1D->GetQuantiles(1, &x, &q);
      // //      std::cout << "median = " << x << std::endl;
      // hadron_energy_fractions[j][k] = x;
      tot_frac[j] +=hadron_energy_fractions[j][k];
      tot_frac_mean[j] += Hadrons_energy_rel_event_eta_MC_pyz_1D->GetMean();

      // //DEBUG ------------------------------
      TCanvas* cHadrons_energy_rel_event_eta_MC_pyz_1D = new TCanvas();
      Hadrons_energy_rel_event_eta_MC_pyz_1D ->SetTitle("");
      TString nameTitle = Hadrons_energy_rel_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
      //      Hadrons_energy_rel_event_eta_MC_pyz_1D ->SetTitle(nameTitle);
      Hadrons_energy_rel_event_eta_MC_pyz_1D ->SetMarkerStyle(20);
      int markerCol=k+1;
      if(markerCol>4) markerCol=k+2;
      Hadrons_energy_rel_event_eta_MC_pyz_1D ->SetMarkerColor(markerCol);
      Hadrons_energy_rel_event_eta_MC_pyz_1D ->Draw("P");

      Hadrons_energy_rel_event_eta_DATA_pyz_1D ->SetMarkerStyle(27);
      Hadrons_energy_rel_event_eta_DATA_pyz_1D ->SetMarkerColor(markerCol);
      Hadrons_energy_rel_event_eta_DATA_pyz_1D ->Draw("same");

      tex->DrawLatex(0.15,0.95,txttag);
      TString eta_lab = ""; eta_lab+=eta_range[j]; eta_lab +=" < |#eta| < "; eta_lab += eta_range[j+1];
      tex->DrawLatex(0.65,0.96,ptbinLabel);
      tex->DrawLatex(0.65,0.92,eta_lab);
      tex->DrawLatex(0.15,0.85,nameTitle);
      TString textMean ="Mean(old)=";textMean+=Form("%.3f",Hadrons_energy_rel_event_eta_MC_pyz_1D->GetMean());
      // // //Median
      Double_t median, q;
      q = 0.5; // 0.5 for "median"
      Hadrons_energy_rel_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      Hadrons_energy_rel_event_eta_MC_pyz_1D->GetQuantiles(1, &median, &q);
      TString textMedian ="Median(old)=";textMedian+=Form("%.3f",median);
      tex->DrawLatex(0.63,0.85,textMean);
      tex->DrawLatex(0.63,0.80,textMedian);
      textMean ="Mean(new)=";textMean+=Form("%.3f",Hadrons_energy_rel_event_eta_DATA_pyz_1D->GetMean());
      Hadrons_energy_rel_event_eta_DATA_pyz_1D->ComputeIntegral(); // just a precaution
      Hadrons_energy_rel_event_eta_DATA_pyz_1D->GetQuantiles(1, &median, &q);
      textMedian ="Median(new)=";textMedian+=Form("%.3f",median);
      tex->DrawLatex(0.63,0.71,textMean);
      tex->DrawLatex(0.63,0.66,textMedian);
      TString textNevents = "Events(old): ";textNevents+=Form("%.0f",Hadrons_energy_rel_event_eta_MC_pyz_1D->GetEntries());
      tex->DrawLatex(0.30,0.85,textNevents);
      textNevents="Events(new): ";textNevents+=Form("%.0f",Hadrons_energy_rel_event_eta_DATA_pyz_1D->GetEntries());
      tex->DrawLatex(0.30,0.79,textNevents);

      TLegend* legHADfrac =  new TLegend(0.19,0.72,0.27,0.84);
      legHADfrac->AddEntry(Hadrons_energy_rel_event_eta_MC_pyz_1D,"old","p");
      legHADfrac->AddEntry(Hadrons_energy_rel_event_eta_DATA_pyz_1D,"new","p");
      legHADfrac->Draw();

      TString CanvasHistName= SavePlots;
      CanvasHistName +="_HadronN"; CanvasHistName+=k; 
      CanvasHistName+="_eta_";CanvasHistName+=eta_range2[j]; CanvasHistName+= "_"; CanvasHistName+= eta_range2[j+1]; 
      CanvasHistName+= "_Hadrons_energy_rel_event_eta_MC_pyz_1D.pdf";
      cHadrons_energy_rel_event_eta_MC_pyz_1D->SetLogy();
      cHadrons_energy_rel_event_eta_MC_pyz_1D->Print(CanvasHistName,"pdf");
      // //[END] DEBUG ------------------------------
      hadron_energy_fractions_median[j][k] = median;
      tot_frac_median[j] += median;
    }
    for(int k=0;k<8;k++){
      if(tot_frac[j]>0)
	hadron_energy_fractions[j][k] = hadron_energy_fractions[j][k]/tot_frac[j];
      else 
	hadron_energy_fractions[j][k] = 0;
      if(tot_frac_median[j]>0)
	hadron_energy_fractions_median[j][k] = hadron_energy_fractions_median[j][k]/tot_frac_median[j];
      else
	hadron_energy_fractions_median[j][k] = 0;
      ptfrac[k]->SetBinContent(j,hadron_energy_fractions[j][k]);
      ptfrac_median[k]->SetBinContent(j,hadron_energy_fractions_median[j][k]);
    }
    delete Hadrons_energy_rel_event_eta_MC_pyz;

    for(int k=0;k<2;k++){
      //      TString tmpName ="GLbinN"; tmpName +=i;tmpName +="_Hadron";tmpName +=name_ext2[k];
      TString tmpName ="GLbinN_MC"; tmpName +=i;tmpName +="_Hadron";tmpName +=k;
      TH1D *PF_to_HAD_event_eta_MC_pyz_1D = (TH1D*)PF_to_HAD_event_eta_MC_pyz->ProjectionY(tmpName,k+1,k+1);
      tmpName ="GLbinN_DATA"; tmpName +=i;tmpName +="_Hadron";tmpName +=k;
      TH1D *PF_to_HAD_event_eta_DATA_pyz_1D = (TH1D*)PF_to_HAD_event_eta_DATA_pyz->ProjectionY(tmpName,k+1,k+1);

      int nentpf = PF_to_HAD_event_eta_MC_pyz->GetEntries();
      pf_had_energy_fractions[j][k] = 0.0;
      pf_had_energy_fractions_DATA[j][k] = 0.0;
      if(nentpf<200) continue;
      // //DEBUG ------------------------------
      TCanvas* cPF_to_HAD_event_eta_MC_pyz_1D = new TCanvas();
      PF_to_HAD_event_eta_MC_pyz_1D ->SetTitle("");
      TString nameTitle = PF_to_HAD_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
      //      PF_to_HAD_event_eta_MC_pyz_1D ->SetTitle(nameTitle);
      PF_to_HAD_event_eta_MC_pyz_1D ->SetMarkerStyle(20);
      PF_to_HAD_event_eta_MC_pyz_1D ->SetMarkerColor(k+1);
      PF_to_HAD_event_eta_MC_pyz_1D ->Draw("P");
      PF_to_HAD_event_eta_DATA_pyz_1D ->SetMarkerStyle(27);
      PF_to_HAD_event_eta_DATA_pyz_1D ->SetMarkerColor(k+1);
      PF_to_HAD_event_eta_DATA_pyz_1D ->Draw("same");

      tex->DrawLatex(0.15,0.95,txttag);
      TString eta_lab = ""; eta_lab+=eta_range[j]; eta_lab +=" < |#eta| < "; eta_lab += eta_range[j+1];
      tex->DrawLatex(0.65,0.85,ptbinLabel);
      tex->DrawLatex(0.65,0.80,eta_lab);
      tex->DrawLatex(0.15,0.85,nameTitle);
      TString textMean ="Mean(old)=";textMean+=Form("%.3f",PF_to_HAD_event_eta_MC_pyz_1D->GetMean());
      // // //Median
      Double_t median, q;
      q = 0.5; // 0.5 for "median"
      PF_to_HAD_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      PF_to_HAD_event_eta_MC_pyz_1D->GetQuantiles(1, &median, &q);
      TString textMedian ="Median(old)=";textMedian+=Form("%.3f",median);
      tex->DrawLatex(0.63,0.68,textMean);
      tex->DrawLatex(0.63,0.63,textMedian);
      textMean ="Mean(new)=";textMean+=Form("%.3f",PF_to_HAD_event_eta_DATA_pyz_1D->GetMean());
      PF_to_HAD_event_eta_DATA_pyz_1D->ComputeIntegral(); // just a precaution
      PF_to_HAD_event_eta_DATA_pyz_1D->GetQuantiles(1, &median, &q);
      textMedian ="Median(new)=";textMedian+=Form("%.3f",median);
      tex->DrawLatex(0.63,0.53,textMean);
      tex->DrawLatex(0.63,0.48,textMedian);


      TString CanvasHistName= SavePlots;
      CanvasHistName +="_HadronN"; CanvasHistName+=k; 
      CanvasHistName+="_eta_";CanvasHistName+=eta_range2[j]; CanvasHistName+= "_"; CanvasHistName+= eta_range2[j+1]; 
      CanvasHistName+= "_PF_to_HAD_event_eta_MC_pyz_1D.pdf";
      cPF_to_HAD_event_eta_MC_pyz_1D->Print(CanvasHistName,"pdf");
      // //[END] DEBUG ------------------------------


      // // // //Mean
      pf_had_energy_fractions[j][k] = PF_to_HAD_event_eta_MC_pyz_1D->GetMean();
      pf_had_energy_fractions_DATA[j][k] = PF_to_HAD_event_eta_DATA_pyz_1D->GetMean();
      //      pf_had_energy_fractions[j][k] = median;
      //      tot_frac_pf_had[j] +=pf_had_energy_fractions[j][k];
    }
    for(int k=0;k<2;k++){
      //      cout<<"pf_had_energy_fractions[j][k]="<<pf_had_energy_fractions[j][k]<<" eta="<<j<<" GLbin="<<i<<endl;
      pfhadfrac[k]->SetBinContent(j,pf_had_energy_fractions[j][k]);
      pfhadfrac_DATA[k]->SetBinContent(j,pf_had_energy_fractions_DATA[j][k]);
    }
    delete PF_to_HAD_event_eta_MC_pyz;

    //N PFs to N hadrons
    for(int k=0;k<2;k++){
      //      TString tmpName ="GLbinN"; tmpName +=i;tmpName +="_Hadron";tmpName +=n_name_ext2[k];
      TString tmpName ="GLbinN_MC_"; tmpName +=i;tmpName +="_Hadron";tmpName +=k;
      TH1D *N_PF_to_HAD_event_eta_MC_pyz_1D = (TH1D*)N_PF_to_HAD_event_eta_MC_pyz->ProjectionY(tmpName,k+1,k+1);
      tmpName ="GLbinN_DATA_"; tmpName +=i;tmpName +="_Hadron";tmpName +=k;
      TH1D *N_PF_to_HAD_event_eta_DATA_pyz_1D = (TH1D*)N_PF_to_HAD_event_eta_DATA_pyz->ProjectionY(tmpName,k+1,k+1);
      int nentpf = N_PF_to_HAD_event_eta_MC_pyz->GetEntries();
      n_pf_had_energy_fractions[j][k] = 0.0;
      n_pf_had_energy_fractions_DATA[j][k] = 0.0;
      if(nentpf<200) continue;
      // //DEBUG ------------------------------
      TCanvas* cN_PF_to_HAD_event_eta_MC_pyz_1D = new TCanvas();
      N_PF_to_HAD_event_eta_MC_pyz_1D ->SetTitle("");
      TString nameTitle = N_PF_to_HAD_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
      //      N_PF_to_HAD_event_eta_MC_pyz_1D ->SetTitle(nameTitle);
      N_PF_to_HAD_event_eta_MC_pyz_1D ->SetMarkerStyle(20);
      N_PF_to_HAD_event_eta_MC_pyz_1D ->SetMarkerColor(k+1);
      N_PF_to_HAD_event_eta_MC_pyz_1D ->Draw("P");
      N_PF_to_HAD_event_eta_DATA_pyz_1D ->SetMarkerStyle(27);
      N_PF_to_HAD_event_eta_DATA_pyz_1D ->SetMarkerColor(k+1);
      N_PF_to_HAD_event_eta_DATA_pyz_1D ->Draw("same");

      tex->DrawLatex(0.15,0.95,txttag);
      TString eta_lab = ""; eta_lab+=eta_range[j]; eta_lab +=" < |#eta| < "; eta_lab += eta_range[j+1];
      tex->DrawLatex(0.65,0.85,ptbinLabel);
      tex->DrawLatex(0.65,0.80,eta_lab);
      tex->DrawLatex(0.15,0.85,nameTitle);
      TString textMean ="Mean(old)=";textMean+=Form("%.3f",N_PF_to_HAD_event_eta_MC_pyz_1D->GetMean());
      // // //Median
      Double_t median, q;
      q = 0.5; // 0.5 for "median"
      N_PF_to_HAD_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      N_PF_to_HAD_event_eta_MC_pyz_1D->GetQuantiles(1, &median, &q);
      TString textMedian ="Median(old)=";textMedian+=Form("%.3f",median);
      tex->DrawLatex(0.63,0.68,textMean);
      tex->DrawLatex(0.63,0.63,textMedian);

      textMean ="Mean(new)=";textMean+=Form("%.3f",N_PF_to_HAD_event_eta_DATA_pyz_1D->GetMean());
      N_PF_to_HAD_event_eta_DATA_pyz_1D->ComputeIntegral(); // just a precaution
      N_PF_to_HAD_event_eta_DATA_pyz_1D->GetQuantiles(1, &median, &q);
      textMedian ="Median(new)=";textMedian+=Form("%.3f",median);
      tex->DrawLatex(0.63,0.53,textMean);
      tex->DrawLatex(0.63,0.48,textMedian);


      TString CanvasHistName= SavePlots;
      CanvasHistName +="_HadronN"; CanvasHistName+=k; 
      CanvasHistName+="_eta_";CanvasHistName+=eta_range2[j]; CanvasHistName+= "_"; CanvasHistName+= eta_range2[j+1]; 
      CanvasHistName+= "_N_PF_to_HAD_event_eta_MC_pyz_1D.pdf";
      cN_PF_to_HAD_event_eta_MC_pyz_1D->Print(CanvasHistName,"pdf");
      // //[END] DEBUG ------------------------------


      // // // //Mean
      n_pf_had_energy_fractions[j][k] = N_PF_to_HAD_event_eta_MC_pyz_1D->GetMean();
      n_pf_had_energy_fractions_DATA[j][k] = N_PF_to_HAD_event_eta_DATA_pyz_1D->GetMean();
      cout<<"AAAAAAAAA n_pf_had_energy_fractions[j][k] = "<<n_pf_had_energy_fractions[j][k]<<" nentpf ="<<nentpf<<" n_pf_had_energy_fractions_DATA[j][k] = "<<n_pf_had_energy_fractions_DATA[j][k]<<endl;
      //      pf_had_energy_fractions[j][k] = median;
      //  tot_frac_pf_had[j] +=n_pf_had_energy_fractions[j][k];
    }
    for(int k=0;k<2;k++){
      //      cout<<"pf_had_energy_fractions[j][k]="<<pf_had_energy_fractions[j][k]<<" eta="<<j<<" GLbin="<<i<<endl;
      n_pfhadfrac[k]->SetBinContent(j,n_pf_had_energy_fractions[j][k]);
      n_pfhadfrac_DATA[k]->SetBinContent(j,n_pf_had_energy_fractions_DATA[j][k]);
    }
    delete N_PF_to_HAD_event_eta_MC_pyz;

    //PF fractions
    tot_frac_pf[j] = 0;
    tot_frac_pf_DATA[j] = 0;
    for(int k=0;k<5;k++){
      TString tmpName ="GLbinN_MC_"; tmpName +=i;tmpName +="_Hadron";tmpName +=k;
      cout<<tmpName.Data()<<endl;
      TH1D *PF_frac_event_eta_MC_pyz_1D = (TH1D*)PF_frac_event_eta_MC_pyz->ProjectionY(tmpName,k+1,k+1); //The projection is made from the channels along the X axis ranging from firstxbin to lastxbin included
      tmpName ="GLbinN_DATA_"; tmpName +=i;tmpName +="_Hadron";tmpName +=k;
      cout<<tmpName.Data()<<endl;
      TH1D *PF_frac_event_eta_DATA_pyz_1D = (TH1D*)PF_frac_event_eta_DATA_pyz->ProjectionY(tmpName,k+1,k+1); //The projection is made from the channels along the X axis ranging from firstxbin to lastxbin included


      int nentpf = PF_frac_event_eta_MC_pyz_1D->GetEntries();
      pf_energy_fractions[j][k]= 0.0;
      pf_energy_fractions_DATA[j][k]= 0.0;
      //      cout<<"total nentpf = "<<nentpf<<", this one = "<<PF_frac_event_eta_MC_pyz_1D->GetEntries()<<endl;
      if(nentpf<200) continue;
      // // //DEBUG ------------------------------
      TCanvas* cPF_frac_event_eta_MC_pyz_1D = new TCanvas();
      PF_frac_event_eta_MC_pyz_1D ->SetTitle("");
      TString nameTitle = PF_frac_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
      //      PF_frac_event_eta_MC_pyz_1D->GetXaxis()->SetRangeUser(-0.01,0.5);
      PF_frac_event_eta_MC_pyz_1D ->SetMarkerStyle(20);
      int markerCol=k+1;
      if(markerCol>4) markerCol=k+2;
      PF_frac_event_eta_MC_pyz_1D ->SetMarkerColor(markerCol);
      PF_frac_event_eta_MC_pyz_1D ->Draw("P");

      PF_frac_event_eta_DATA_pyz_1D ->SetMarkerStyle(27);
      PF_frac_event_eta_DATA_pyz_1D ->SetMarkerColor(markerCol);
      PF_frac_event_eta_DATA_pyz_1D ->Draw("same");
      TLegend* legPFfrac =  new TLegend(0.19,0.72,0.27,0.84);
      legPFfrac->AddEntry(PF_frac_event_eta_MC_pyz_1D,"old","p");
      legPFfrac->AddEntry(PF_frac_event_eta_DATA_pyz_1D,"new","p");
      legPFfrac->Draw();

      tex->DrawLatex(0.15,0.95,txttag);
      TString eta_lab = ""; eta_lab+=eta_range[j]; eta_lab +=" < |#eta| < "; eta_lab += eta_range[j+1];
      tex->DrawLatex(0.65,0.96,ptbinLabel);
      tex->DrawLatex(0.65,0.91,eta_lab);
      tex->DrawLatex(0.15,0.85,nameTitle);
      TString textMean ="Mean(old)=";textMean+=Form("%.3f",PF_frac_event_eta_MC_pyz_1D->GetMean());
      // // //Median
      Double_t median, q;
      q = 0.5; // 0.5 for "median"
      PF_frac_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      PF_frac_event_eta_MC_pyz_1D->GetQuantiles(1, &median, &q);
      TString textMedian ="Median(old)=";textMedian+=Form("%.3f",median);
      tex->DrawLatex(0.63,0.84,textMean);
      tex->DrawLatex(0.63,0.78,textMedian);

      TString textMeanNew ="Mean(new)=";textMeanNew+=Form("%.3f",PF_frac_event_eta_DATA_pyz_1D->GetMean());
      // // //Median
      PF_frac_event_eta_DATA_pyz_1D->ComputeIntegral(); // just a precaution
      PF_frac_event_eta_DATA_pyz_1D->GetQuantiles(1, &median, &q);
      TString textMedianNew ="Median(new)=";textMedianNew+=Form("%.3f",median);
      tex->DrawLatex(0.63,0.68,textMeanNew);
      tex->DrawLatex(0.63,0.62,textMedianNew);

      TString CanvasHistName= SavePlots;
      CanvasHistName +="_FracN"; CanvasHistName+=k; 
      CanvasHistName+="_eta_";CanvasHistName+=eta_range2[j]; CanvasHistName+= "_"; CanvasHistName+= eta_range2[j+1]; 
      CanvasHistName+= "_PF_frac_event_eta_MC_pyz_1D.pdf";
      cPF_frac_event_eta_MC_pyz_1D->SetLogy();
      cPF_frac_event_eta_MC_pyz_1D->Print(CanvasHistName,"pdf");
      // CanvasHistName+= "_PF_frac_event_eta_MC_pyz_1D.root";
      // cPF_frac_event_eta_MC_pyz_1D->Print(CanvasHistName,"root");
      // //[END] DEBUG ------------------------------


      // // // //Mean
      pf_energy_fractions[j][k] = PF_frac_event_eta_MC_pyz_1D->GetMean();
      pf_energy_fractions_DATA[j][k] = PF_frac_event_eta_DATA_pyz_1D->GetMean();
      //      pf_energy_fractions[j][k] = median;
      tot_frac_pf[j] +=pf_energy_fractions[j][k];
      tot_frac_pf_DATA[j] +=pf_energy_fractions_DATA[j][k];
      //      delete PF_frac_event_eta_MC_pyz_1D;//TEST
      //      cout<<"BEFORE norm: mean PF fractions["<<k<<"]="<<pf_energy_fractions[j][k]<<" norm = "<<tot_frac_pf[j]<<" label="<<name_ext3[k]<<" nentpf="<<nentpf<<endl;
    }
    for(int k=0;k<5;k++){
      if(tot_frac_pf[j]>0)
      	pf_energy_fractions[j][k] = pf_energy_fractions[j][k]/tot_frac_pf[j];
      else 
	pf_energy_fractions[j][k] = 0;
      pffrac[k]->SetBinContent(j,pf_energy_fractions[j][k]);
      if(tot_frac_pf_DATA[j]>0)
      	pf_energy_fractions_DATA[j][k] = pf_energy_fractions_DATA[j][k]/tot_frac_pf_DATA[j];
      else 
	pf_energy_fractions_DATA[j][k] = 0;
      pffrac_DATA[k]->SetBinContent(j,pf_energy_fractions_DATA[j][k]);


    }
    delete PF_frac_event_eta_MC_pyz;
    delete PF_frac_event_eta_DATA_pyz;
  }

  for(int k=0;k<8;k++){
    hs_energyfrac_median->Add(ptfrac_median[k]);
    hs_energyfrac->Add(ptfrac[k]);
    if(k<2){ 
      pfhadfrac[k]->GetYaxis()->SetRangeUser(0,3.1);
      hs_pfhadfrac->Add(pfhadfrac[k]);
      pfhadfrac_DATA[k]->GetYaxis()->SetRangeUser(0,3.1);
      hs_pfhadfrac_DATA->Add(pfhadfrac_DATA[k]);

      n_pfhadfrac[k]->GetYaxis()->SetRangeUser(0,4.1);
      hs_n_pfhadfrac->Add(n_pfhadfrac[k]);

      n_pfhadfrac_DATA[k]->GetYaxis()->SetRangeUser(0,4.1);
      hs_n_pfhadfrac_DATA->Add(n_pfhadfrac_DATA[k]);

    }
    if(k<5){
      hs_pffrac->Add(pffrac[k]);
      hs_pffrac_DATA->Add(pffrac_DATA[k]);
    }
  }
  

  TCanvas* cHadrons_energy_rel_event_eta_MC = new TCanvas();
  hs_energyfrac->Draw("HIST");
  leg2 -> Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cHadrons_energy_rel_event_eta_MC->Print(SavePlots + "_Hadrons_energy_rel_vs_eta_MC.pdf","pdf");

  TCanvas* cHadrons_energy_rel_event_eta_MC_median = new TCanvas();
  hs_energyfrac_median->Draw("HIST");
  leg2_median -> Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cHadrons_energy_rel_event_eta_MC_median->Print(SavePlots + "_Hadrons_energy_rel_vs_eta_MC_median.pdf","pdf");


  TCanvas* cPF_event_eta_MC = new TCanvas();
  hs_pffrac->Draw("HIST");
  hs_pffrac_DATA->Draw("HIST P SAME");
  leg4 -> Draw();
  legMCDATA -> Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);

  cPF_event_eta_MC->Print(SavePlots + "_PF_frac_energy_vs_eta_MC.pdf","pdf");

  TCanvas* cPF_to_HAD_event_eta_MC = new TCanvas();
  // // Upper plot will be in pad1
  // TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  // pad1->SetBottomMargin(0); // Upper and lower plot are joined
  // pad1->SetGridx();         // Vertical grid
  // pad1->Draw();             // Draw the upper pad: pad1
  // pad1->cd();               // pad1 becomes the current pad
  hs_pfhadfrac->Draw("NOSTACK");
  hs_pfhadfrac_DATA->Draw("NOSTACK P SAME");
  leg3 -> Draw();
  legMCDATA -> Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  // // Do not draw the Y axis label on the upper plot and redraw a small
  // // axis instead, in order to avoid the first label (0) to be clipped.
  // hs_pfhadfrac->GetYaxis()->SetLabelSize(0.);
  // TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  // axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  // axis->SetLabelSize(15);
  // axis->Draw();
  // cPF_to_HAD_event_eta_MC->cd();          // Go back to the main canvas before defining pad2
  // TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  // pad2->SetTopMargin(0);
  // pad2->SetBottomMargin(0.2);
  // pad2->SetGridx(); // vertical grid
  // pad2->Draw();
  // pad2->cd();       // pad2 becomes the current pad
  cPF_to_HAD_event_eta_MC->Print(SavePlots + "_PF_to_HAD_energy_vs_eta_MC.pdf","pdf");

  TCanvas* cN_PF_to_HAD_event_eta_MC = new TCanvas();
  hs_n_pfhadfrac->Draw("NOSTACK");
  hs_n_pfhadfrac_DATA->Draw("NOSTACK P SAME");
  leg5 -> Draw();
  legMCDATA -> Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cN_PF_to_HAD_event_eta_MC->Print(SavePlots + "_N_PF_to_HAD_energy_vs_eta_MC.pdf","pdf");

 


  for(int k=0;k<8;k++){
    delete ptfrac[k];
    delete ptfrac_median[k];
    if(k<2){
      delete pfhadfrac[k];
      delete pfhadfrac_DATA[k];
    }
    if(k<5){
      delete pffrac[k];
      delete pffrac_DATA[k];
    }
  }

  delete cHadrons_energy_MC;
  delete cHadrons_energy_rel_MC;
  delete cHadrons_energy_rel_event_MC;
  delete cHadrons_count_event_MC;
  delete cHadrons_genjeteta_MC;
  delete cHadrons_total_energy_rel_event_MC;
  delete cHadrons_total_count_event_MC;
  delete cHadrons_energy_rel_event_eta_MC;
  delete cResponse_eta_MC;
  delete cResponse_eta_MC_mean;
  delete cResponse_eta_MC_std;
  };
}
