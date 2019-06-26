// Read histograms filled by BaconJet/AnalysisModule module for hadron composition of GEN probe jet
//created on 19.12.2018 by A. Karavdina

#include "../include/CorrectionObject.h"
//#include "../include/parameters.h"
#include "../../include/constants.h"
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

void CorrectionObject::ControlPlotsHadrons(bool forEverySingleTrigger){
  cout << "--------------- Starting ControlPlotsHadrons() ---------------" << endl;
  bool forEverySingleTrigger_ = forEverySingleTrigger;
  if(forEverySingleTrigger_){
    cout<<"ControlPlotsHadrons are gonne be made for all trigger seperatly"<<endl;
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



 TString dirNames[21] = {"Hadrons", "Hadrons_3rdJet", "Hadrons_BB", "Hadrons_EC1", "Hadrons_EC2", "Hadrons_HF",
"Hadrons_HLT_DiPFJetAve40", "Hadrons_HLT_DiPFJetAve60", "Hadrons_HLT_DiPFJetAve80","Hadrons_HLT_DiPFJetAve140","Hadrons_HLT_DiPFJetAve200","Hadrons_HLT_DiPFJetAve260","Hadrons_HLT_DiPFJetAve320","Hadrons_HLT_DiPFJetAve400","Hadrons_HLT_DiPFJetAve500", "Hadrons_HLT_DiPFJetAve60_HFJEC", "Hadrons_HLT_DiPFJetAve80_HFJEC", "Hadrons_HLT_DiPFJetAve100_HFJEC", "Hadrons_HLT_DiPFJetAve160_HFJEC", "Hadrons_HLT_DiPFJetAve220_HFJEC", "Hadrons_HLT_DiPFJetAve300_HFJEC"};


 TString ptbinsLabels[21] = {"73#leqp^{RECO}_{T}#leq#infty", "3rd jet, 15#leqp^{RECO}_{T}#leq#infty", 
			     "0.261#leq#eta#leq0.522","2.172#leq#eta#leq2.322",
			     "2.50#leq#eta#leq2.65","3.4892#leq#eta#leq3.839",
			     "73#leqp^{RECO}_{T}#leq85", "85#leqp^{RECO}_{T}#leq97","97#leqp^{RECO}_{T}#leq179","179#leqp^{RECO}_{T}#leq307",
			     "307#leqp^{RECO}_{T}#leq370","370#leqp^{RECO}_{T}#leq434","434#leqp^{RECO}_{T}#leq520","520#leqp^{RECO}_{T}#leq649", "649#leqp^{RECO}_{T}#leq#infty", 
			     "73#leqp^{RECO}_{T}#leq93", "93#leqp^{RECO}_{T}#leq113", "113#leqp^{RECO}_{T}#leq176", "176#leqp^{RECO}_{T}#leq239", "239#leqp^{RECO}_{T}#leq318", "318#leqp^{RECO}_{T}#leq#infty"};

  int loop_inter = 2;
  if(forEverySingleTrigger) loop_inter=21;
  
  for(int i = 0; i<loop_inter ; i++){
  //for(int i = 0; i<1 ; i++){//TEST
  
   TString dirName = (TString)dirNames[i];  
   TString ptbinLabel = (TString)ptbinsLabels[i];
   cout<<endl;  
   cout<<dirName<<endl;
      
  TString MCtitle = "MC";
  TString SavePlots = CorrectionObject::_outpath + "plots/control/ControlPlotsHadrons/CPHadrons_" + dirName + "_" + CorrectionObject::_generator_tag;
  CorrectionObject::make_path(std::string((_outpath + "plots/").Data()));
  CorrectionObject::make_path(std::string((_outpath + "plots/control/").Data()));
  CorrectionObject::make_path(std::string((_outpath + "plots/control/ControlPlotsHadrons/").Data()));

  cout<<"++++++++++++++++ Collect all histograms ++++++++++++++++\n";
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



  bool use_mc=true;
  if (Hadrons_energy_MC->GetEntries()==0 && Hadrons_genjeteta_MC->GetEntries()==0){
    cout<<"no MC found, wont be used"<<endl;
    use_mc=false;
  }
  if(!use_mc) continue;
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);  

  cout<<"+++++++++++++++++++++ Plot extremely nice histograms ++++++++++++++++++++++++++++\n";
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
    ptfrac[k]->SetLineColor(k+1);
    ptfrac[k]->SetFillColor(k+1);
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
    ptfrac_median[k]->SetLineColor(k+1);
    ptfrac_median[k]->SetFillColor(k+1);
    leg2_median ->AddEntry(ptfrac_median[k],name_ext[k],"f");
  }


  //PF to HAD fractions
  TH1D *pfhadfrac[2]; 
  TString name_ext2[2];
  THStack *hs_pfhadfrac = new THStack("pf_to_hadron_energy_fractions_hs","; RECO jet #eta;<PF fraction/HAD fraction>");
  TLegend* leg3 =  new TLegend(0.40,0.68,0.88,0.88);
  for(int k=0;k<2;k++){
    name_ext2[k] = PF_to_HAD_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
    TString name = "pffrac_to_hadfrac_"+name_ext2[k];
    name +="_GLbinN"; name +=i;  
    pfhadfrac[k] = new TH1D(name,"",n_eta-1, eta_bins);
    pfhadfrac[k]->SetLineColor(k+1);
    pfhadfrac[k]->SetFillColorAlpha(k+1,1-0.4*k);
    for(int j=0; j<n_eta_common-1; j++)
      pfhadfrac[k]->SetBinContent(j,0.0);
    leg3 ->AddEntry(pfhadfrac[k],name_ext2[k],"f");
  }
  //N PF to HAD constituence
  TH1D *n_pfhadfrac[2]; 
  TString n_name_ext2[2];
  THStack *hs_n_pfhadfrac = new THStack("n_pf_to_hadron_energy_fractions_hs","; RECO jet #eta;<N PFs/N hadrons>");
  TLegend* leg5 =  new TLegend(0.40,0.68,0.88,0.88);
  for(int k=0;k<2;k++){
    n_name_ext2[k] = N_PF_to_HAD_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
    TString name = "n_pffrac_to_hadfrac_"+name_ext2[k];
    name +="_GLbinN"; name +=i;  
    n_pfhadfrac[k] = new TH1D(name,"",n_eta-1, eta_bins);
    n_pfhadfrac[k]->SetLineColor(k+1);
    n_pfhadfrac[k]->SetFillColorAlpha(k+1,1-0.4*k);
    for(int j=0; j<n_eta_common-1; j++)
      n_pfhadfrac[k]->SetBinContent(j,0.0);
    leg5 ->AddEntry(n_pfhadfrac[k],name_ext2[k],"f");
  }

  //PF  fractions
  TH1D *pffrac[5]; 
  TString name_ext3[5];
  THStack *hs_pffrac = new THStack("pf_energy_fractions_hs","; RECO jet #eta;<PF fraction>");
  TLegend* leg4 =  new TLegend(0.85,0.30,0.99,0.99);
  for(int k=0;k<5;k++){
    name_ext3[k] = PF_frac_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
    TString name = "pffrac_"+name_ext3[k]; 
    name +="_GLbinN"; name +=i;
    pffrac[k] = new TH1D(name,"",n_eta-1, eta_bins);
    pffrac[k]->SetLineColor(k+1);
    pffrac[k]->SetFillColorAlpha(k+1,1-0.1*k);
    leg4 ->AddEntry(pffrac[k],name_ext3[k],"f");
  }

  double tot_frac[n_eta_common];   double tot_frac_mean[n_eta_common];  double tot_frac_median[n_eta_common];
  double tot_frac_pf[n_eta_common];
  //  double tot_frac_pf_had[n_eta_common];
  for(int j=0; j<n_eta_common-1; j++){

    Hadrons_energy_rel_event_eta_MC->GetXaxis()->SetRange(j,j+1);

    PF_frac_event_eta_MC->GetXaxis()->SetRange(j,j+1);
    PF_to_HAD_event_eta_MC->GetXaxis()->SetRange(j,j+1);

    TH2D *Hadrons_energy_rel_event_eta_MC_pyz = (TH2D*)Hadrons_energy_rel_event_eta_MC->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *PF_frac_event_eta_MC_pyz = (TH2D*)PF_frac_event_eta_MC->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *PF_to_HAD_event_eta_MC_pyz = (TH2D*)PF_to_HAD_event_eta_MC->Project3D("yz");//fraction vs hadron in 1 eta bin
    TH2D *N_PF_to_HAD_event_eta_MC_pyz = (TH2D*)N_PF_to_HAD_event_eta_MC->Project3D("yz");//fraction vs hadron in 1 eta bin
    tot_frac[j] = 0; tot_frac_mean[j] = 0;
    for(int k=0;k<8;k++){
      TString tmpName ="GLbinN"; tmpName +=i;tmpName +="_Hadron";tmpName +=name_ext[k];
      //      TH1D *Hadrons_energy_rel_event_eta_MC_pyz_1D = (TH1D*)Hadrons_energy_rel_event_eta_MC_pyz->ProjectionY(tmpName,k,k+1);
      TH1D *Hadrons_energy_rel_event_eta_MC_pyz_1D = (TH1D*)Hadrons_energy_rel_event_eta_MC_pyz->ProjectionY(tmpName,k+1,k+1);//TEST
      int netn = Hadrons_energy_rel_event_eta_MC_pyz->GetEntries();
      hadron_energy_fractions[j][k] = 0.0;
      hadron_energy_fractions_median[j][k] = 0.0;
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
      Hadrons_energy_rel_event_eta_MC_pyz_1D ->SetMarkerColor(k+1);
      Hadrons_energy_rel_event_eta_MC_pyz_1D ->Draw("P");

      tex->DrawLatex(0.15,0.95,txttag);
      TString eta_lab = ""; eta_lab+=eta_range[j]; eta_lab +=" < |#eta| < "; eta_lab += eta_range[j+1];
      tex->DrawLatex(0.65,0.85,ptbinLabel);
      tex->DrawLatex(0.65,0.80,eta_lab);
      tex->DrawLatex(0.15,0.85,nameTitle);
      TString textMean ="Mean=";textMean+=Form("%.3f",Hadrons_energy_rel_event_eta_MC_pyz_1D->GetMean());
      // // //Median
      Double_t median, q;
      q = 0.5; // 0.5 for "median"
      Hadrons_energy_rel_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      Hadrons_energy_rel_event_eta_MC_pyz_1D->GetQuantiles(1, &median, &q);
      TString textMedian ="Median=";textMedian+=Form("%.3f",median);
      tex->DrawLatex(0.67,0.68,textMean);
      tex->DrawLatex(0.67,0.58,textMedian);
      TString textNevents = "Events: ";textNevents+=Form("%.0f",Hadrons_energy_rel_event_eta_MC_pyz_1D->GetEntries());
      tex->DrawLatex(0.30,0.85,textNevents);
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
    //    delete Hadrons_energy_rel_event_eta_MC_pyz;

    for(int k=0;k<2;k++){
      TString tmpName ="GLbinN"; tmpName +=i;tmpName +="_Hadron";tmpName +=name_ext2[k];
      // TH1D *PF_to_HAD_event_eta_MC_pyz_1D = (TH1D*)PF_to_HAD_event_eta_MC_pyz->ProjectionY(tmpName,k,k+1);
      TH1D *PF_to_HAD_event_eta_MC_pyz_1D = (TH1D*)PF_to_HAD_event_eta_MC_pyz->ProjectionY(tmpName,k+1,k+1); //TEST
      int nentpf = PF_to_HAD_event_eta_MC_pyz->GetEntries();
      pf_had_energy_fractions[j][k] = 0.0;
      if(nentpf<200) continue;
      // // //DEBUG ------------------------------
      // TCanvas* cPF_to_HAD_event_eta_MC_pyz_1D = new TCanvas();
      // PF_to_HAD_event_eta_MC_pyz_1D ->SetTitle("");
      // TString nameTitle = PF_to_HAD_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
      // //      PF_to_HAD_event_eta_MC_pyz_1D ->SetTitle(nameTitle);
      // PF_to_HAD_event_eta_MC_pyz_1D ->SetMarkerStyle(20);
      // PF_to_HAD_event_eta_MC_pyz_1D ->SetMarkerColor(k+1);
      // PF_to_HAD_event_eta_MC_pyz_1D ->Draw("P");
      // tex->DrawLatex(0.15,0.95,txttag);
      // TString eta_lab = ""; eta_lab+=eta_range[j]; eta_lab +=" < |#eta| < "; eta_lab += eta_range[j+1];
      // tex->DrawLatex(0.65,0.85,ptbinLabel);
      // tex->DrawLatex(0.65,0.80,eta_lab);
      // tex->DrawLatex(0.15,0.85,nameTitle);
      // TString textMean ="Mean=";textMean+=Form("%.3f",PF_to_HAD_event_eta_MC_pyz_1D->GetMean());
      // // // //Median
      // Double_t median, q;
      // q = 0.5; // 0.5 for "median"
      // PF_to_HAD_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      // PF_to_HAD_event_eta_MC_pyz_1D->GetQuantiles(1, &median, &q);
      // TString textMedian ="Median=";textMedian+=Form("%.3f",median);
      // tex->DrawLatex(0.67,0.68,textMean);
      // tex->DrawLatex(0.67,0.58,textMedian);
      // TString CanvasHistName= SavePlots;
      // CanvasHistName +="_HadronN"; CanvasHistName+=k; 
      // CanvasHistName+="_eta_";CanvasHistName+=eta_range2[j]; CanvasHistName+= "_"; CanvasHistName+= eta_range2[j+1]; 
      // CanvasHistName+= "_PF_to_HAD_event_eta_MC_pyz_1D.pdf";
      // cPF_to_HAD_event_eta_MC_pyz_1D->Print(CanvasHistName,"pdf");
      // // //[END] DEBUG ------------------------------


      // // // //Mean
      pf_had_energy_fractions[j][k] = PF_to_HAD_event_eta_MC_pyz_1D->GetMean();
      //      pf_had_energy_fractions[j][k] = median;
      //      tot_frac_pf_had[j] +=pf_had_energy_fractions[j][k];
    }
    for(int k=0;k<2;k++){
      //      cout<<"pf_had_energy_fractions[j][k]="<<pf_had_energy_fractions[j][k]<<" eta="<<j<<" GLbin="<<i<<endl;
      pfhadfrac[k]->SetBinContent(j,pf_had_energy_fractions[j][k]);
    }
    delete PF_to_HAD_event_eta_MC_pyz;

    //N PFs to N hadrons
    for(int k=0;k<2;k++){
      TString tmpName ="GLbinN"; tmpName +=i;tmpName +="_Hadron";tmpName +=n_name_ext2[k];
      TH1D *N_PF_to_HAD_event_eta_MC_pyz_1D = (TH1D*)N_PF_to_HAD_event_eta_MC_pyz->ProjectionY(tmpName,k+1,k+1);
      int nentpf = N_PF_to_HAD_event_eta_MC_pyz->GetEntries();
      n_pf_had_energy_fractions[j][k] = 0.0;
      if(nentpf<200) continue;
      // // //DEBUG ------------------------------
      // TCanvas* cN_PF_to_HAD_event_eta_MC_pyz_1D = new TCanvas();
      // N_PF_to_HAD_event_eta_MC_pyz_1D ->SetTitle("");
      // TString nameTitle = N_PF_to_HAD_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
      // //      N_PF_to_HAD_event_eta_MC_pyz_1D ->SetTitle(nameTitle);
      // N_PF_to_HAD_event_eta_MC_pyz_1D ->SetMarkerStyle(20);
      // N_PF_to_HAD_event_eta_MC_pyz_1D ->SetMarkerColor(k+1);
      // N_PF_to_HAD_event_eta_MC_pyz_1D ->Draw("P");
      // tex->DrawLatex(0.15,0.95,txttag);
      // TString eta_lab = ""; eta_lab+=eta_range[j]; eta_lab +=" < |#eta| < "; eta_lab += eta_range[j+1];
      // tex->DrawLatex(0.65,0.85,ptbinLabel);
      // tex->DrawLatex(0.65,0.80,eta_lab);
      // tex->DrawLatex(0.15,0.85,nameTitle);
      // TString textMean ="Mean=";textMean+=Form("%.3f",N_PF_to_HAD_event_eta_MC_pyz_1D->GetMean());
      // // // //Median
      // Double_t median, q;
      // q = 0.5; // 0.5 for "median"
      // N_PF_to_HAD_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      // N_PF_to_HAD_event_eta_MC_pyz_1D->GetQuantiles(1, &median, &q);
      // TString textMedian ="Median=";textMedian+=Form("%.3f",median);
      // tex->DrawLatex(0.67,0.68,textMean);
      // tex->DrawLatex(0.67,0.58,textMedian);
      // TString CanvasHistName= SavePlots;
      // CanvasHistName +="_HadronN"; CanvasHistName+=k; 
      // CanvasHistName+="_eta_";CanvasHistName+=eta_range2[j]; CanvasHistName+= "_"; CanvasHistName+= eta_range2[j+1]; 
      // CanvasHistName+= "_N_PF_to_HAD_event_eta_MC_pyz_1D.pdf";
      // cN_PF_to_HAD_event_eta_MC_pyz_1D->Print(CanvasHistName,"pdf");
      // // //[END] DEBUG ------------------------------


      // // // //Mean
      n_pf_had_energy_fractions[j][k] = N_PF_to_HAD_event_eta_MC_pyz_1D->GetMean();
      //      pf_had_energy_fractions[j][k] = median;
      //  tot_frac_pf_had[j] +=n_pf_had_energy_fractions[j][k];
    }
    for(int k=0;k<2;k++){
      //      cout<<"pf_had_energy_fractions[j][k]="<<pf_had_energy_fractions[j][k]<<" eta="<<j<<" GLbin="<<i<<endl;
      n_pfhadfrac[k]->SetBinContent(j,n_pf_had_energy_fractions[j][k]);
    }
    delete N_PF_to_HAD_event_eta_MC_pyz;

    //PF fractions
    tot_frac_pf[j] = 0;
    for(int k=0;k<5;k++){
      TString tmpName ="GLbinN"; tmpName +=i;tmpName +="_Hadron";tmpName +=name_ext2[k];
      TH1D *PF_frac_event_eta_MC_pyz_1D = (TH1D*)PF_frac_event_eta_MC_pyz->ProjectionY(tmpName,k+1,k+1); //The projection is made from the channels along the X axis ranging from firstxbin to lastxbin included
      //      TH2D *PF_frac_event_eta_MC_pyz_1D = PF_frac_event_eta_MC_pyz;//TEST
      int nentpf = PF_frac_event_eta_MC_pyz_1D->GetEntries();
      pf_energy_fractions[j][k]= 0.0;
      //      cout<<"total nentpf = "<<nentpf<<", this one = "<<PF_frac_event_eta_MC_pyz_1D->GetEntries()<<endl;
      if(nentpf<200) continue;
      // // //DEBUG ------------------------------
      // TCanvas* cPF_frac_event_eta_MC_pyz_1D = new TCanvas();
      // PF_frac_event_eta_MC_pyz_1D ->SetTitle("");
      // TString nameTitle = PF_frac_event_eta_MC->GetZaxis()->GetBinLabel(k+1);
      // //      PF_frac_event_eta_MC_pyz_1D->GetXaxis()->SetRangeUser(-0.01,0.5);
      // PF_frac_event_eta_MC_pyz_1D ->SetMarkerStyle(20);
      // PF_frac_event_eta_MC_pyz_1D ->SetMarkerColor(k+1);
      // PF_frac_event_eta_MC_pyz_1D ->Draw("P");
      // tex->DrawLatex(0.15,0.95,txttag);
      // TString eta_lab = ""; eta_lab+=eta_range[j]; eta_lab +=" < |#eta| < "; eta_lab += eta_range[j+1];
      // tex->DrawLatex(0.65,0.85,ptbinLabel);
      // tex->DrawLatex(0.65,0.80,eta_lab);
      // tex->DrawLatex(0.15,0.85,nameTitle);
      // TString textMean ="Mean=";textMean+=Form("%.3f",PF_frac_event_eta_MC_pyz_1D->GetMean());
      // // // //Median
      // Double_t median, q;
      // q = 0.5; // 0.5 for "median"
      // PF_frac_event_eta_MC_pyz_1D->ComputeIntegral(); // just a precaution
      // PF_frac_event_eta_MC_pyz_1D->GetQuantiles(1, &median, &q);
      // TString textMedian ="Median=";textMedian+=Form("%.3f",median);
      // tex->DrawLatex(0.67,0.68,textMean);
      // tex->DrawLatex(0.67,0.58,textMedian);
      // TString CanvasHistName= SavePlots;
      // CanvasHistName +="_HadronN"; CanvasHistName+=k; 
      // CanvasHistName+="_eta_";CanvasHistName+=eta_range2[j]; CanvasHistName+= "_"; CanvasHistName+= eta_range2[j+1]; 
      // CanvasHistName+= "_PF_frac_event_eta_MC_pyz_1D.pdf";
      // cPF_frac_event_eta_MC_pyz_1D->Print(CanvasHistName,"pdf");
      // // CanvasHistName+= "_PF_frac_event_eta_MC_pyz_1D.root";
      // // cPF_frac_event_eta_MC_pyz_1D->Print(CanvasHistName,"root");
      // // //[END] DEBUG ------------------------------


      // // // //Mean
      pf_energy_fractions[j][k] = PF_frac_event_eta_MC_pyz_1D->GetMean();
      //      pf_energy_fractions[j][k] = median;
      tot_frac_pf[j] +=pf_energy_fractions[j][k];
      //      delete PF_frac_event_eta_MC_pyz_1D;//TEST
      //      cout<<"BEFORE norm: mean PF fractions["<<k<<"]="<<pf_energy_fractions[j][k]<<" norm = "<<tot_frac_pf[j]<<" label="<<name_ext3[k]<<" nentpf="<<nentpf<<endl;
    }
    for(int k=0;k<5;k++){
      if(tot_frac_pf[j]>0)
      	pf_energy_fractions[j][k] = pf_energy_fractions[j][k]/tot_frac_pf[j];
      else 
	pf_energy_fractions[j][k] = 0;
      //      cout<<"AFTER norm: mean PF fractions ["<<k<<"]="<<pf_energy_fractions[j][k]<<" norm = "<<tot_frac_pf[j]<<" label="<<name_ext3[k]<<endl;
      pffrac[k]->SetBinContent(j,pf_energy_fractions[j][k]);
    }
    delete PF_frac_event_eta_MC_pyz;
  }

  for(int k=0;k<8;k++){
    hs_energyfrac_median->Add(ptfrac_median[k]);
    hs_energyfrac->Add(ptfrac[k]);
    if(k<2){ 
      pfhadfrac[k]->GetYaxis()->SetRangeUser(0,3.1);
      hs_pfhadfrac->Add(pfhadfrac[k]);
      n_pfhadfrac[k]->GetYaxis()->SetRangeUser(0,3.1);
      hs_n_pfhadfrac->Add(n_pfhadfrac[k]);
    }
    if(k<5) hs_pffrac->Add(pffrac[k]);
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
  leg4 -> Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cPF_event_eta_MC->Print(SavePlots + "_PF_frac_energy_vs_eta_MC.pdf","pdf");

  TCanvas* cPF_to_HAD_event_eta_MC = new TCanvas();
  //  hs_pfhadfrac->Draw("HIST");
  hs_pfhadfrac->Draw("NOSTACK");
  leg3 -> Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cPF_to_HAD_event_eta_MC->Print(SavePlots + "_PF_to_HAD_energy_vs_eta_MC.pdf","pdf");

  TCanvas* cN_PF_to_HAD_event_eta_MC = new TCanvas();
  hs_n_pfhadfrac->Draw("NOSTACK");
  leg5 -> Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  tex->DrawLatex(0.68,0.95,ptbinLabel);
  cN_PF_to_HAD_event_eta_MC->Print(SavePlots + "_N_PF_to_HAD_energy_vs_eta_MC.pdf","pdf");



  for(int k=0;k<8;k++){
    delete ptfrac[k];
    delete ptfrac_median[k];
    if(k<2) delete pfhadfrac[k];
    if(k<5) delete pffrac[k];
  }

  delete cHadrons_energy_MC;
  delete cHadrons_energy_rel_MC;
  delete cHadrons_energy_rel_event_MC;
  delete cHadrons_count_event_MC;
  delete cHadrons_genjeteta_MC;
  delete cHadrons_total_energy_rel_event_MC;
  delete cHadrons_total_count_event_MC;
  delete cHadrons_energy_rel_event_eta_MC;
  };
}
