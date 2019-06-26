#include "../include/CorrectionObject.h"
// #include "../../include/helper_fkt.h"
#include "../include/parameters.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>

using namespace std;

void CorrectionObject::ControlPlots(bool forEverySingleTrigger){
  cout << "--------------- Starting ControlPlots() ---------------" << endl;
  bool forEverySingleTrigger_ = forEverySingleTrigger;
  if(forEverySingleTrigger_){
    cout<<"ControlPlots are gonne be made for all trigger seperatly"<<endl;
  }
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);
  
  TString dirNames[16] = {"Selection", "HLT_DiPFJetAve40", "HLT_DiPFJetAve60", "HLT_DiPFJetAve80","HLT_DiPFJetAve140","HLT_DiPFJetAve200","HLT_DiPFJetAve260","HLT_DiPFJetAve320","HLT_DiPFJetAve400","HLT_DiPFJetAve500", "HLT_DiPFJetAve60_HFJEC", "HLT_DiPFJetAve80_HFJEC", "HLT_DiPFJetAve100_HFJEC", "HLT_DiPFJetAve160_HFJEC", "HLT_DiPFJetAve220_HFJEC", "HLT_DiPFJetAve300_HFJEC"};

  int loop_inter = 1;
  if(forEverySingleTrigger) loop_inter=16;
  
  for(int i = 0; i<loop_inter ; i++){

   TString dirName = (TString)dirNames[i];
   
    cout<<endl;  
    cout<<"DirName: "<<dirName<<endl;
      
  TString DATAtitle = "DATA";
  TString MCtitle = "MC";
  TString SavePlots = CorrectionObject::_outpath + "plots/control/ControlPlots_" + dirName + "_" + CorrectionObject::_generator_tag;
 CorrectionObject::make_path(std::string((_outpath + "plots/control/").Data()));

 cout<<"++++++++++++++++ Collect all histograms ++++++++++++++++\n";
  TH1F* pt_jet1_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/pt_1");
  TH1F* pt_jet2_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/pt_2");
  TH1F* pt_jet3_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/pt_3");
  TH1F* pt_jet1_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/pt_1");
  TH1F* pt_jet2_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/pt_2");
  TH1F* pt_jet3_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/pt_3");
  TH1F* eta_jet1_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/eta_1");
  TH1F* eta_jet2_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/eta_2");
  TH1F* eta_jet1_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/eta_1");
  TH1F* eta_jet2_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/eta_2");
  TH1F *Njets_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/N_jets");
  TH1F *Njets_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/N_jets");
  TH1F *pt_ave_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/pt_ave");
  TH1F *pt_ave_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/pt_ave");
  TH1F *asym_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/asym");
  TH1F *asym_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/asym");
  TH1F *r_dijet_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/r_rel");
  TH1F *r_dijet_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/r_rel");
  TH1F *r_mpf_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/mpf");
  TH1F *r_mpf_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/mpf");
  TH1F *N_PV_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/N_PV");
  TH1F *N_PV_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/N_PV");
  TH1F *Eta_pos_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/eta_probe_pos");
  TH1F *Eta_neg_DATA = (TH1F*)CorrectionObject::_DATAFile->Get(dirName+"/eta_probe_neg");
  TH1F *Eta_assym_top_DATA = (TH1F*)Eta_pos_DATA->Clone();
  Eta_assym_top_DATA->Add(Eta_neg_DATA,-1);
  TH1F *Eta_assym_bot_DATA = (TH1F*)Eta_pos_DATA->Clone();
  Eta_assym_bot_DATA->Add(Eta_neg_DATA,+1);
  Eta_assym_top_DATA->Divide(Eta_assym_bot_DATA);
  TH1F *Eta_pos_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/eta_probe_pos");
  TH1F *Eta_neg_MC = (TH1F*)CorrectionObject::_MCFile->Get(dirName+"/eta_probe_neg");
  TH1F *Eta_assym_top_MC = (TH1F*)Eta_pos_MC->Clone();
  Eta_assym_top_MC->Add(Eta_neg_MC,-1);
  TH1F *Eta_assym_bot_MC = (TH1F*)Eta_pos_MC->Clone();
  Eta_assym_bot_MC->Add(Eta_neg_MC,+1);
  Eta_assym_top_MC->Divide(Eta_assym_bot_MC);

  bool use_mc=true;
  if (pt_jet1_MC->GetEntries()==0 && pt_jet2_MC->GetEntries()==0){
    cout<<"no MC found, wont be used"<<endl;
    use_mc=false;
  }
  else{
    cout<<"Number of events in MC: "<<pt_jet1_MC->GetEntries()<<endl;
  }

  cout<<"+++++++++++++++++++++ Plot extremely nice histograms ++++++++++++++++++++++++++++\n";
  TCanvas* a = new TCanvas();
  a->Divide(3,2);

  if(use_mc) cout<<"mc is used"<<endl;
  
  a->cd(1);
  pt_jet1_DATA->SetMarkerStyle(20);
  pt_jet1_DATA->SetMarkerSize(0.5);
  pt_jet1_DATA->SetMarkerColor(1);
  pt_jet1_DATA->SetLineColor(1);
  pt_jet1_DATA->GetXaxis()->SetTitle("p_{T}, GeV");
  pt_jet1_DATA->Scale(1/pt_jet1_DATA->Integral());
if(use_mc)     pt_jet1_DATA->GetYaxis()->SetRangeUser(0,0.1);
 // cout<<pt_jet1_DATA<<endl;
  pt_jet1_DATA->Draw();
  pt_jet1_MC->SetMarkerStyle(22);
  pt_jet1_MC->SetMarkerSize(0.5);
  pt_jet1_MC->SetMarkerColor(2);
  pt_jet1_MC->SetLineColor(2);
  if(use_mc){
    cout<<"pt_jet1_MC->Integral() = "<<pt_jet1_MC->Integral()<<endl;
    cout<<"pt_jet2_MC->Integral() = "<<pt_jet2_MC->Integral()<<endl;
    cout<<"pt_jet3_MC->Integral() = "<<pt_jet3_MC->Integral()<<endl;
    cout<<"pt_jet1_DATA->Integral() = "<<pt_jet1_DATA->Integral()<<endl;
    cout<<"pt_jet2_DATA->Integral() = "<<pt_jet2_DATA->Integral()<<endl;
    cout<<"pt_jet3_DATA->Integral() = "<<pt_jet3_DATA->Integral()<<endl;

  }
  if(use_mc)  pt_jet1_MC->Scale(1/pt_jet1_MC->Integral());
  if(use_mc)    pt_jet1_MC->Draw("same");

  TLegend *leg1;
  leg1 = new TLegend(0.63,0.68,0.9,0.85,"","brNDC");//x+0.1
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetFillColor(10);
  leg1->SetLineColor(1);
  leg1->SetTextFont(42);
  leg1->AddEntry(pt_jet1_DATA,DATAtitle,"lp");
  if(use_mc)    leg1->AddEntry(pt_jet1_MC,MCtitle,"lp");
  leg1->Draw();

  a->cd(2);
  pt_jet2_DATA->SetMarkerStyle(20);
  pt_jet2_DATA->SetMarkerSize(0.5);
  pt_jet2_DATA->SetMarkerColor(1);
  pt_jet2_DATA->SetLineColor(1);
  pt_jet2_DATA->GetXaxis()->SetTitle("p_{T}, GeV");
  pt_jet2_DATA->Scale(1/pt_jet2_DATA->Integral());
if(use_mc)     pt_jet2_DATA->GetYaxis()->SetRangeUser(0,0.1);
 // cout<<pt_jet2_DATA<<endl;
  pt_jet2_DATA->Draw();
  pt_jet2_MC->SetMarkerStyle(22);
  pt_jet2_MC->SetMarkerSize(0.5);
  pt_jet2_MC->SetMarkerColor(2);
  pt_jet2_MC->SetLineColor(2);
  if(use_mc)    pt_jet2_MC->Scale(1/pt_jet2_MC->Integral());
  if(use_mc)    pt_jet2_MC->Draw("same");

  a->cd(3);
  pt_jet3_DATA->SetMarkerStyle(20);
  pt_jet3_DATA->SetMarkerSize(0.5);
  pt_jet3_DATA->SetMarkerColor(1);
  pt_jet3_DATA->SetLineColor(1);
  pt_jet3_DATA->GetXaxis()->SetTitle("p_{T}, GeV");
  pt_jet3_DATA->Scale(1/pt_jet3_DATA->Integral());
  if(use_mc)   pt_jet3_DATA->GetYaxis()->SetRangeUser(0,0.1);
  // cout<<pt_jet3_DATA<<endl;
  pt_jet3_DATA->Draw();
  pt_jet3_MC->SetMarkerStyle(22);
  pt_jet3_MC->SetMarkerSize(0.5);
  pt_jet3_MC->SetMarkerColor(2);
  pt_jet3_MC->SetLineColor(2);
  if(use_mc)    pt_jet3_MC->Scale(1/pt_jet3_MC->Integral());
  if(use_mc)    pt_jet3_MC->Draw("same");

  a->cd(4);
  eta_jet1_DATA->SetMarkerStyle(20);
  eta_jet1_DATA->SetMarkerSize(0.5);
  eta_jet1_DATA->SetMarkerColor(1);
  eta_jet1_DATA->SetLineColor(1);
  eta_jet1_DATA->GetXaxis()->SetTitle("#eta");
  eta_jet1_DATA->Scale(1./eta_jet1_DATA->Integral());
  if(use_mc)   eta_jet1_DATA->GetYaxis()->SetRangeUser(0,0.04);
  // cout<<eta_jet1_DATA<<endl;
  eta_jet1_DATA->Draw();
  eta_jet1_MC->SetMarkerStyle(22);
  eta_jet1_MC->SetMarkerSize(0.5);
  eta_jet1_MC->SetMarkerColor(2);
  eta_jet1_MC->SetLineColor(2);
  if(use_mc)     eta_jet1_MC->Scale(1./eta_jet1_MC->Integral());
  if(use_mc)     eta_jet1_MC->Draw("same");

  a->cd(5);
  eta_jet2_DATA->SetMarkerStyle(20);
  eta_jet2_DATA->SetMarkerSize(0.5);
  eta_jet2_DATA->SetMarkerColor(1);
  eta_jet2_DATA->SetLineColor(1);
  eta_jet2_DATA->GetXaxis()->SetTitle("#eta");
  eta_jet2_DATA->Scale(1./eta_jet2_DATA->Integral());
  if(use_mc)   eta_jet2_DATA->GetYaxis()->SetRangeUser(0,0.04);
  // cout<<eta_jet2_DATA<<endl;
  eta_jet2_DATA->Draw();
  eta_jet2_MC->SetMarkerStyle(22);
  eta_jet2_MC->SetMarkerSize(0.5);
  eta_jet2_MC->SetMarkerColor(2);
  eta_jet2_MC->SetLineColor(2);
  if(use_mc)    eta_jet2_MC->Scale(1./eta_jet2_MC->Integral());
  if(use_mc)    eta_jet2_MC->Draw("same");
 
  a->cd(6);
  Njets_DATA->SetMarkerStyle(20);
  Njets_DATA->SetMarkerSize(0.5);
  Njets_DATA->SetMarkerColor(1);
  Njets_DATA->SetLineColor(1);
  Njets_DATA->GetXaxis()->SetTitle("Number of jets");
  Njets_DATA->Scale(1./Njets_DATA->Integral());
  if(use_mc)   Njets_DATA->GetYaxis()->SetRangeUser(0,0.1);
  // cout<<Njets_DATA<<endl;
  Njets_DATA->Draw();
  Njets_MC->SetMarkerStyle(22);
  Njets_MC->SetMarkerSize(0.5);
  Njets_MC->SetMarkerColor(2);
  Njets_MC->SetLineColor(2);  
  if(use_mc)    Njets_MC->Scale(1./Njets_MC->Integral());
  if(use_mc)    Njets_MC->Draw("same");

  cout << "Printing plots to " << SavePlots + "_jets.pdf" << endl;
    // cout<<"1\n";
    // cout<<SavePlots<<endl;
    // cout<<a<<endl;
    // a->Draw();
         // cout<<"1.2\n"; 
  a->Print(SavePlots + "_jets.pdf","pdf");
  cout<<" Random message 1 \n";

  TCanvas* b = new TCanvas();
  b->Divide(3,2);

  b->cd(1);
  N_PV_DATA->SetMarkerStyle(20);
  N_PV_DATA->SetMarkerSize(0.5);
  N_PV_DATA->SetMarkerColor(1);
  N_PV_DATA->SetLineColor(1);
  N_PV_DATA->GetXaxis()->SetTitle("N good PV");
  N_PV_DATA->Scale(1./N_PV_DATA->Integral());
  if(use_mc)   N_PV_DATA->GetYaxis()->SetRangeUser(0,0.1);
  N_PV_DATA->Draw();
  N_PV_MC->SetMarkerStyle(22);
  N_PV_MC->SetMarkerSize(0.5);
  N_PV_MC->SetMarkerColor(2);
  N_PV_MC->SetLineColor(2);
   if(use_mc)   N_PV_MC->Scale(1./N_PV_MC->Integral());
  if(use_mc)    N_PV_MC->Draw("same");

  TLegend *leg2;
  leg2 = new TLegend(0.63,0.68,0.87,0.86,"","brNDC");//x+0.1
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  leg2->SetFillColor(10);
  leg2->SetLineColor(1);
  leg2->SetTextFont(42);
  leg2->AddEntry(Njets_DATA,DATAtitle,"lp");
  if(use_mc)    leg2->AddEntry(Njets_MC,MCtitle,"lp");
  leg2->Draw();

  b->cd(2);
  leg2->Draw();

  b->cd(3);
  pt_ave_DATA->SetMarkerStyle(20);
  pt_ave_DATA->SetMarkerSize(0.5);
  pt_ave_DATA->SetMarkerColor(1);
  pt_ave_DATA->SetLineColor(1);
  pt_ave_DATA->Scale(1./pt_ave_DATA->Integral());
 if(use_mc)    pt_ave_DATA->GetYaxis()->SetRangeUser(0,0.03);
  pt_ave_DATA->Draw();
  pt_ave_MC->SetMarkerStyle(22);
  pt_ave_MC->SetMarkerSize(0.5);
  pt_ave_MC->SetMarkerColor(2);
  pt_ave_MC->SetLineColor(2);
  if(use_mc)    pt_ave_MC->Scale(1./pt_ave_MC->Integral());
  if(use_mc)    pt_ave_MC->Draw("same");

    cout<<"2\n";
  b->cd(4);
  gPad->SetLogy();
  asym_DATA->SetMarkerStyle(20);
  asym_DATA->SetMarkerSize(0.5);
  asym_DATA->SetMarkerColor(1);
  asym_DATA->SetLineColor(1);
  asym_DATA->Scale(1./asym_DATA->Integral());
  asym_DATA->Draw();
  asym_MC->SetMarkerStyle(22);
  asym_MC->SetMarkerSize(0.5);
  asym_MC->SetMarkerColor(2);
  asym_MC->SetLineColor(2);
  asym_MC->Scale(1./asym_MC->Integral());
  if(use_mc)    asym_MC->Draw("same");
  leg1->Draw();

  b->cd(5);
  r_dijet_DATA->SetMarkerStyle(20);
  r_dijet_DATA->SetMarkerSize(0.5);
  r_dijet_DATA->SetMarkerColor(1);
  r_dijet_DATA->SetLineColor(1);
  r_dijet_DATA->Scale(1./r_dijet_DATA->Integral());
 if(use_mc)    r_dijet_DATA->GetYaxis()->SetRangeUser(0,0.03);
  r_dijet_DATA->Draw();
  r_dijet_MC->SetMarkerStyle(22);
  r_dijet_MC->SetMarkerSize(0.5);
  r_dijet_MC->SetMarkerColor(2);
  r_dijet_MC->SetLineColor(2);
  if(use_mc)    r_dijet_MC->Scale(1./r_dijet_MC->Integral());
  if(use_mc)    r_dijet_MC->Draw("same");
  leg1->Draw();

  b->cd(6);
  r_mpf_DATA->SetMarkerStyle(20);
  r_mpf_DATA->SetMarkerSize(0.5);
  r_mpf_DATA->SetMarkerColor(1);
  r_mpf_DATA->SetLineColor(1);
  r_mpf_DATA->Scale(1./r_mpf_DATA->Integral());
 if(use_mc)    r_mpf_DATA->GetYaxis()->SetRangeUser(0,0.04);
  r_mpf_DATA->Draw();
  r_mpf_MC->SetMarkerStyle(22);
  r_mpf_MC->SetMarkerSize(0.5);
  r_mpf_MC->SetMarkerColor(2);
  r_mpf_MC->SetLineColor(2);
  if(use_mc)    r_mpf_MC->Scale(1./r_mpf_MC->Integral());
  if(use_mc)    r_mpf_MC->Draw("same");

  b->Print(SavePlots + "_dijet.pdf");


    cout<<"3\n";

  TCanvas* f = new TCanvas();
  TH2F *mpf_vs_etaProbe_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/mpf_vs_etaProbe");
  TH2F *mpf_vs_etaProbe_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/mpf_vs_etaProbe");
  TH2F *r_rel_vs_etaProbe_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/r_rel_vs_etaProbe");
  TH2F *r_rel_vs_etaProbe_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/r_rel_vs_etaProbe");

  f->Divide(2,2);								    
  f->cd(1);									    
  r_rel_vs_etaProbe_DATA->SetTitle(DATAtitle);				    
  r_rel_vs_etaProbe_DATA->Draw("colz");					    
  f->cd(2);									    
  r_rel_vs_etaProbe_MC->SetTitle(MCtitle);					    
  r_rel_vs_etaProbe_MC->Draw("colz");						    
  f->cd(3);									    
  mpf_vs_etaProbe_DATA->SetTitle(DATAtitle);					    
  mpf_vs_etaProbe_DATA->Draw("colz");						    
  f->cd(4);									    
  mpf_vs_etaProbe_MC->SetTitle(MCtitle);					    
  if(use_mc)    mpf_vs_etaProbe_MC->Draw("colz");						    
  
  f->Print(SavePlots + "_ResponsesVsEta.pdf");
  
  TCanvas* fpeprobe = new TCanvas();    
  TH2F *phi_vs_etaProbe_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/phi_vs_etaProbe");
  phi_vs_etaProbe_DATA->SetTitle(DATAtitle);
  phi_vs_etaProbe_DATA->DrawNormalized("colz");
  fpeprobe->Print(SavePlots + "_PhiVsEtaProbe.pdf");
 
  TCanvas* fpe = new TCanvas();    
  TH2F *phi_vs_eta_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/phi_vs_etaAll");
  phi_vs_eta_DATA->SetTitle(DATAtitle);
  phi_vs_eta_DATA->DrawNormalized("colz");
  fpe->Print(SavePlots + "_PhiVsEtaAll.pdf");

    TCanvas* fpec = new TCanvas();    
    TH2F *phi_vs_etaC_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/phi_vs_etaAllCut30");
    if(phi_vs_etaC_DATA!=NULL){
      phi_vs_etaC_DATA->SetTitle(DATAtitle);
      phi_vs_etaC_DATA->DrawNormalized("colz");
      fpec->Print(SavePlots + "_PhiVsEtaAllCut30.pdf");
  }
  
  TCanvas* g = new TCanvas();
  TH2F *pt_ave_vs_etaProbe_DATA = (TH2F*)CorrectionObject::_DATAFile->Get(dirName+"/pt_ave_vs_etaProbe");
  TH2F *pt_ave_vs_etaProbe_MC = (TH2F*)CorrectionObject::_MCFile->Get(dirName+"/pt_ave_vs_etaProbe");

  use_mc = !(pt_ave_vs_etaProbe_MC->GetEntries()==0);
  
  TLine *lineEta[2*n_eta];
  for(int i=0;i<n_eta;i++){
    lineEta[i] = new TLine(-1*eta_bins[i],0,-1*eta_bins[i],1000);
    lineEta[i]->SetLineColor(2);
    lineEta[i]->SetLineStyle(2);
  }
  for(int i=n_eta;i<2*n_eta;i++){
    lineEta[i] = new TLine(eta_bins[i-n_eta],0,eta_bins[i-n_eta],1000);
    lineEta[i]->SetLineColor(2);
    lineEta[i]->SetLineStyle(2);
    lineEta[i]->SetLineWidth(2);
  }

  TLine *linePt[n_pt];
  for(int i=0;i<n_pt;i++){
    linePt[i] = new TLine(-5.2,pt_bins[i],5.2,pt_bins[i]);
    linePt[i]->SetLineColor(2);
    linePt[i]->SetLineStyle(2);
    linePt[i]->SetLineWidth(2);
  }

  g->Divide(2,1);
  g->cd(1);
  pt_ave_vs_etaProbe_DATA->SetTitle(DATAtitle);
  pt_ave_vs_etaProbe_DATA->Draw("colz");
  for(int i=0;i<2*n_eta;i++)
    lineEta[i]->Draw("same");
  for(int i=0;i<n_pt;i++)
    linePt[i]->Draw("same");
  g->cd(2);
  pt_ave_vs_etaProbe_MC->SetTitle(MCtitle);
 if(use_mc)  pt_ave_vs_etaProbe_MC->Draw("colz");
  for(int i=0;i<2*n_eta;i++)
    lineEta[i]->Draw("same");
  for(int i=0;i<n_pt;i++)
    linePt[i]->Draw("same");

  g->Print(SavePlots + "_Pt_aveVsEtaProbe.pdf");




  //delete everything

  for(int i=0; i<n_pt; i++) delete linePt[i];
  for(int i=0; i<2*n_eta; i++) delete lineEta[i];
  delete g;
  delete f;
  delete leg2;
  delete b;
  delete leg1;
  delete a;

  };
}
