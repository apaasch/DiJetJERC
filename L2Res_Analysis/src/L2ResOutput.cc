#include<TStyle.h>
#include<TLine.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include "../include/tdrstyle_mod15.h"
#include <TFile.h>
#include <iostream>
#include <assert.h>
#include "../include/CorrectionObject.h"
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include "../include/parameters.h"

using namespace std;

void CorrectionObject::L2ResOutput(){
  cout << "--------------- Starting L2ResOutput() ---------------" << endl << endl;

  TFile* f_Res_mpf = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");
  TFile* f_Res_dijet = new TFile(CorrectionObject::_outpath+"Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");  

  // TFile* f_Res_mpf_old   = new TFile("/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V7/AK4CHS/MC_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor_w_JER_noEtaPhiCleaning/RunBCDEFGH/Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");
  // TFile* f_Res_dijet_old = new TFile("/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V7/AK4CHS/MC_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor_w_JER_noEtaPhiCleaning/RunBCDEFGH/Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");  
  TString txttag=CorrectionObject::_generator_tag; 
  TString JetDescrib;                                                                                                                            
  if (CorrectionObject::_collection=="AK4CHS") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";
  if (CorrectionObject::_collection=="AK4Puppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";
  if (CorrectionObject::_collection=="AK8CHS") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";
  if (CorrectionObject::_collection=="AK8Puppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 

  //plot results for "nominal" variation
 

  // get the (R_{MC}/R_{DATA}) hists for MPF and pt balance
  TH1D* pt_depend_const_mpf = (TH1D*)f_Res_mpf->Get("ptave_const_mpf");
  TH1D* pt_depend_logpt_mpf = (TH1D*)f_Res_mpf->Get("ptave_logpt_mpf");
  TH1D* pt_depend_const_dijet = (TH1D*)f_Res_dijet->Get("ptave_const_dijet");
  TH1D* pt_depend_logpt_dijet = (TH1D*)f_Res_dijet->Get("ptave_logpt_dijet");

  // get the kFSR hists for MPF and pt balance
  TH1D* kfsr_mpf = (TH1D*)f_Res_mpf->Get("kfsr_mpf");
  TH1D* kfsr_mpf_fit = (TH1D*)f_Res_mpf->Get("hist_kfsr_fit_mpf");
  TH1D* kfsr_dijet = (TH1D*)f_Res_dijet->Get("kfsr_dijet");
  TH1D* kfsr_dijet_fit = (TH1D*)f_Res_dijet->Get("hist_kfsr_fit_dijet");

  //get L2Res hists for MPF and pt balance
  TH1D* res_const_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_const_mpf");
  TH1D* res_const_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_const_dijet");
  TH1D* res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
  TH1D* res_logpt_dijet_kfsrfit = (TH1D*)f_Res_dijet->Get("res_logpt_dijet");

  //************************* get L2Res hists for MPF and pt balance WITHOUT FIT VALUES****************************
  TH1D* res_const_mpf_kfsrfit_val = (TH1D*)f_Res_mpf->Get("res_const_mpf_val");
  TH1D* res_const_dijet_kfsrfit_val = (TH1D*)f_Res_dijet->Get("res_const_dijet_val");
  TH1D* res_logpt_mpf_kfsrfit_val = (TH1D*)f_Res_mpf->Get("res_logpt_mpf_val");
  TH1D* res_logpt_dijet_kfsrfit_val = (TH1D*)f_Res_dijet->Get("res_logpt_dijet_val");
  //***************************************************************************************************************  

  // //get L2Res hists for MPF and pt balance
  // TH1D* res_const_mpf_kfsrfit_old = (TH1D*)f_Res_mpf_old->Get("res_const_mpf");
  // TH1D* res_const_dijet_kfsrfit_old = (TH1D*)f_Res_dijet_old->Get("res_const_dijet");
  // TH1D* res_logpt_mpf_kfsrfit_old = (TH1D*)f_Res_mpf_old->Get("res_logpt_mpf");
  // TH1D* res_logpt_dijet_kfsrfit_old = (TH1D*)f_Res_dijet_old->Get("res_logpt_dijet");

  TH1D* res_mpf_pT_diff = (TH1D*) res_logpt_mpf_kfsrfit->Clone("res_mpf_pT_diff");
  res_mpf_pT_diff->Add(res_logpt_dijet_kfsrfit,-1);

 
  // TH1D* res_const_mpf_diff = (TH1D*) res_const_mpf_kfsrfit->Clone("res_const_mpf_diff");
  // res_const_mpf_diff->Add(res_const_mpf_kfsrfit_old,-1);

  // TH1D* res_const_dijet_diff = (TH1D*) res_const_dijet_kfsrfit->Clone("res_const_dijet_diff");
  // res_const_dijet_diff->Add(res_const_dijet_kfsrfit_old,-1);

  // TH1D* res_logpt_mpf_diff = (TH1D*) res_logpt_mpf_kfsrfit->Clone("res_logpt_mpf_diff");
  // res_logpt_mpf_diff->Add(res_logpt_mpf_kfsrfit_old,-1);

  // TH1D* res_logpt_dijet_diff = (TH1D*) res_logpt_dijet_kfsrfit->Clone("res_logpt_dijet_diff");
  // res_logpt_dijet_diff->Add(res_logpt_dijet_kfsrfit_old,-1);
  


  //default Canvas
  TCanvas* c1 = new TCanvas();

  // create horizontal line for plotting ("ideal value")
  TLine *line = new TLine(0,1,5.191,1);
  TLine *line2 = new TLine(0,0,5.191,0);
  res_const_mpf_kfsrfit->SetLineWidth(2);
  res_const_mpf_kfsrfit->SetLineColor(kRed+1);
  res_const_dijet_kfsrfit->SetLineWidth(2);
  res_const_dijet_kfsrfit->SetLineColor(kBlue+1);
  res_const_mpf_kfsrfit->SetLineStyle(2);
  res_const_dijet_kfsrfit->SetLineStyle(2);
  res_logpt_mpf_kfsrfit->SetLineWidth(2);
  res_logpt_mpf_kfsrfit->SetLineColor(kRed+1);
  res_logpt_dijet_kfsrfit->SetLineWidth(2);
  res_logpt_dijet_kfsrfit->SetLineColor(kBlue+1);

  //**************************  kFSR NO FIT Value **********************************
  res_const_mpf_kfsrfit_val->SetLineWidth(2);
  res_const_mpf_kfsrfit_val->SetLineColor(kRed+1);
  res_const_dijet_kfsrfit_val->SetLineWidth(2);
  res_const_dijet_kfsrfit_val->SetLineColor(kBlue+1);
  res_const_mpf_kfsrfit_val->SetLineStyle(2);
  res_const_dijet_kfsrfit_val->SetLineStyle(2);
  res_logpt_mpf_kfsrfit_val->SetLineWidth(2);
  res_logpt_mpf_kfsrfit_val->SetLineColor(kRed+1);
  res_logpt_dijet_kfsrfit_val->SetLineWidth(2);
  res_logpt_dijet_kfsrfit_val->SetLineColor(kBlue+1);
  //*********************************************************************************

  pt_depend_const_mpf->SetLineWidth(2);
  pt_depend_const_mpf->SetLineColor(kRed+1);
  pt_depend_const_dijet->SetLineWidth(2);
  pt_depend_const_dijet->SetLineColor(kBlue+1);
  pt_depend_const_mpf->SetLineStyle(2);
  pt_depend_const_dijet->SetLineStyle(2);
  pt_depend_logpt_mpf->SetLineWidth(2);
  pt_depend_logpt_mpf->SetLineColor(kRed+1);
  pt_depend_logpt_dijet->SetLineWidth(2);
  pt_depend_logpt_dijet->SetLineColor(kBlue+1);

  kfsr_mpf->SetLineWidth(2);
  kfsr_mpf->SetLineColor(kRed+1);
  kfsr_dijet->SetLineWidth(2);
  kfsr_dijet->SetLineColor(kBlue+1);
  
  bool kSquare = true;
  // Draw results
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);
  //h->GetXaxis()->SetTitleSize(0.05);
  //h->GetYaxis()->SetTitleSize(0.05);
  //lumi_13TeV = CorrectionObject::_lumitag;
  tdrCanvas(c1,"c1",h,4,10,kSquare,CorrectionObject::_lumitag);
  
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);  
  tex->DrawLatex(0.45,0.87,JetDescrib);

  TLatex *tex_bars = new TLatex();
  tex_bars->SetNDC();
  tex_bars->SetTextSize(0.039);

  TLegend leg1 = tdrLeg(0.20,0.19,0.43,0.40);
  leg1.AddEntry(pt_depend_const_mpf, "MPF Flat","L");
  leg1.AddEntry(pt_depend_const_dijet, "Pt Flat","L");
  leg1.AddEntry(pt_depend_logpt_mpf, "MPF Loglin","L");
  leg1.AddEntry(pt_depend_logpt_dijet, "Pt Loglin","L"); 
  
  TCanvas* c2 = new TCanvas();
  tdrCanvas(c2,"c2",h,4,10,kSquare,CorrectionObject::_lumitag);
  
  TString alVal;
  alVal.Form("%0.2f\n",alpha_cut);
  TString altitle = "{#alpha<"+alVal+"}";
  TString axistitle = "(R^{MC}/R^{data})_";
  axistitle +=altitle;
  h->GetYaxis()->SetTitle(axistitle);
  h->GetYaxis()->SetRangeUser(0.81,1.15); 

  pt_depend_const_mpf->SetMarkerStyle(1);
  pt_depend_const_dijet->SetMarkerStyle(1);
  pt_depend_logpt_mpf->SetMarkerStyle(1);
  pt_depend_logpt_dijet->SetMarkerStyle(1);
  pt_depend_const_mpf->Draw("E1 SAME");
  pt_depend_const_dijet->Draw("E1 SAME");
  pt_depend_logpt_mpf->Draw("E1 SAME");
  pt_depend_logpt_dijet->Draw("E1 SAME");

  line->Draw("SAME");
  leg1.Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  tex->DrawLatex(0.15,0.95,txttag);
  c2->SaveAs(CorrectionObject::_outpath+"plots/Ratio_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");

  
  TCanvas *c3 = new TCanvas();
  tdrCanvas(c3,"c3",h,4,10,kSquare,CorrectionObject::_lumitag);
  h->GetYaxis()->SetTitle("k_{FSR}");
  //  h->GetYaxis()->SetRangeUser(0.81,1.15);
  h->GetYaxis()->SetRangeUser(0.94,1.06);
  kfsr_dijet_fit->SetMarkerStyle(1);
  kfsr_mpf_fit->SetMarkerStyle(1);
  kfsr_mpf->SetMarkerStyle(1);
  kfsr_dijet->SetMarkerStyle(1);
  kfsr_dijet_fit->Draw("E3 SAME");
  kfsr_mpf_fit->Draw("E3 SAME");
  kfsr_mpf->Draw("E1 SAME");
  kfsr_dijet->Draw("E1 SAME");
  line->Draw("SAME");

  TLegend leg2 = tdrLeg(0.20,0.19,0.43,0.30);
  leg2 . AddEntry(kfsr_mpf, "MPF","L");
  leg2 . AddEntry(kfsr_dijet, "Pt","L");
  leg2.Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  tex->DrawLatex(0.15,0.95,txttag);
  c3->SaveAs(CorrectionObject::_outpath+"plots/kFSR_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");


  TCanvas *c4 = new TCanvas();
  tdrCanvas(c4,"L2res_kFSRfit",h,4,10,kSquare,CorrectionObject::_lumitag);
  h->GetYaxis()->SetTitle("Relative correction");
  h->GetYaxis()->SetRangeUser(0.81,1.15);
  res_const_mpf_kfsrfit->SetMarkerStyle(1);
  res_const_dijet_kfsrfit->SetMarkerStyle(1);
  res_logpt_mpf_kfsrfit->SetMarkerStyle(1);
  res_logpt_dijet_kfsrfit->SetMarkerStyle(1);
  res_const_mpf_kfsrfit->Draw("E1 SAME");
  res_const_dijet_kfsrfit->Draw("E1 SAME");
  res_logpt_mpf_kfsrfit->Draw("E1 SAME");
  res_logpt_dijet_kfsrfit->Draw("E1 SAME");
  leg1.Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  line->Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  c4->SaveAs(CorrectionObject::_outpath+"plots/L2Res_kFSRfit_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");

 TCanvas *c10 = new TCanvas();
  tdrCanvas(c10,"L2res_kFSRfit",h,4,10,kSquare,CorrectionObject::_lumitag);
  h->GetYaxis()->SetTitle("Relative correction");
  h->GetYaxis()->SetRangeUser(0.81,1.15);
  res_const_mpf_kfsrfit_val->SetMarkerStyle(1);
  res_const_dijet_kfsrfit_val->SetMarkerStyle(1);
  res_logpt_mpf_kfsrfit_val->SetMarkerStyle(1);
  res_logpt_dijet_kfsrfit_val->SetMarkerStyle(1);
  res_const_mpf_kfsrfit_val->Draw("E1 SAME");
  res_const_dijet_kfsrfit_val->Draw("E1 SAME");
  res_logpt_mpf_kfsrfit_val->Draw("E1 SAME");
  res_logpt_dijet_kfsrfit_val->Draw("E1 SAME");
  leg1.Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);
  line->Draw();
  tex->DrawLatex(0.15,0.95,txttag);
  c10->SaveAs(CorrectionObject::_outpath+"plots/L2Res_kFSRval_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");


  //pt-dependence of L2Res (Money plot)

  TString var="";
  TH1D* res_logpt_mpf_kfsrfit_var[4];
  TH1D* res_logpt_dijet_kfsrfit_var[4];
  
  // TH1D* res_logpt_mpf_kfsrfit_old_var[4];
  // TH1D* res_logpt_dijet_kfsrfit_old_var[4];

  TH1D* res_mpf_pT_diff_var[4];

  // TH1D* res_logpt_mpf_diff_var[4];
  // TH1D* res_logpt_dijet_diff_var[4];

  TFile* f_Res_mpf_var;
  TFile* f_Res_dijet_var;

  // TFile* f_Res_mpf_diff_var;
  // TFile* f_Res_dijet_diff_var;

  TFile* f_Res_mpf_pT_var;

  // TFile* f_Res_mpf_old_var;
  // TFile* f_Res_dijet_old_var;

  for(int i=0;i<4;i++){
    if(i==0) var="central"; 
    if(i==1) var="down";
    if(i==2) var="up";
    if(i==3) var="doubleup";

    f_Res_mpf_var = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ");
    f_Res_dijet_var = new TFile(CorrectionObject::_outpath+"Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ");

    f_Res_mpf_pT_var = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ");

    // f_Res_mpf_diff_var = new TFile(CorrectionObject::_outpath+"Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ");
    // f_Res_dijet_diff_var = new TFile(CorrectionObject::_outpath+"Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_"+var+".root","READ"); 


    // f_Res_mpf_old_var   = new TFile("/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V7/AK4CHS/MC_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor_w_JER_noEtaPhiCleaning/RunBCDEFGH/Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");
    // f_Res_dijet_old_var = new TFile("/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V7/AK4CHS/MC_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor_w_JER_noEtaPhiCleaning/RunBCDEFGH/Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");  


    res_logpt_mpf_kfsrfit_var[i] = (TH1D*)f_Res_mpf_var->Get("res_logpt_mpf");
    res_logpt_dijet_kfsrfit_var[i] = (TH1D*)f_Res_dijet_var->Get("res_logpt_dijet");

    // res_logpt_mpf_kfsrfit_old_var[i] = (TH1D*)f_Res_mpf_old_var->Get("res_logpt_mpf");
    // res_logpt_dijet_kfsrfit_old_var[i] = (TH1D*)f_Res_dijet_old_var->Get("res_logpt_dijet");
    
    res_mpf_pT_diff_var[i] =  (TH1D*)f_Res_mpf_pT_var->Get("res_logpt_mpf");
    res_mpf_pT_diff_var[i] ->Add(res_logpt_dijet_kfsrfit_var[i],-1);
    
    // res_logpt_mpf_diff_var[i] = (TH1D*)f_Res_mpf_diff_var->Get("res_logpt_mpf");
    // res_logpt_mpf_diff_var[i] -> Add(res_logpt_mpf_kfsrfit_old_var[i],-1);
    
    // res_logpt_dijet_diff_var[i] = (TH1D*)f_Res_dijet_diff_var->Get("res_logpt_dijet");
    // res_logpt_dijet_diff_var[i] -> Add(res_logpt_dijet_kfsrfit_old_var[i],-1);
    
  }

  TCanvas *c5 = new TCanvas();
  tdrCanvas(c5,"L2res_logpt_MPF_kFSRfit_ptDepend",h,4,10,kSquare,CorrectionObject::_lumitag);
  res_logpt_mpf_kfsrfit->SetLineColor(kBlack);
  res_logpt_mpf_kfsrfit->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_logpt_mpf_kfsrfit_var[i]->SetLineColor(kRed-3*i);
    res_logpt_mpf_kfsrfit_var[i]->SetMarkerColor(kRed-3*i);
    res_logpt_mpf_kfsrfit_var[i]->SetMarkerStyle(20+i);
    res_logpt_mpf_kfsrfit_var[i]->Draw("E1 SAME"); 
  }

  TLegend leg3 = tdrLeg(0.20,0.19,0.43,0.42); 
  leg3 . AddEntry(res_logpt_mpf_kfsrfit_var[1] , "60 GeV","LP");  
  leg3 . AddEntry(res_logpt_mpf_kfsrfit_var[0] , "120 GeV","LP");
  leg3 . AddEntry(res_logpt_mpf_kfsrfit_var[2] , "240 GeV","LP"); 
  leg3 . AddEntry(res_logpt_mpf_kfsrfit_var[3] , "480 GeV","LP");
  leg3 . AddEntry(res_logpt_mpf_kfsrfit, "Nominal","LP");
  leg3.Draw();              
  tex->DrawLatex(0.45,0.87,JetDescrib);      
  tex->DrawLatex(0.15,0.95,txttag);
  c5->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_MPF_kFSRfit_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");

  
  TCanvas *c7 = new TCanvas();
  tdrCanvas(c7,"L2res_logpt_DiJet_kFSRfit_ptDepend",h,4,10,kSquare,CorrectionObject::_lumitag);
  res_logpt_dijet_kfsrfit->SetLineColor(kBlack);
  res_logpt_dijet_kfsrfit->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_logpt_dijet_kfsrfit_var[i]->SetLineColor(kBlue+3*i);
    res_logpt_dijet_kfsrfit_var[i]->SetMarkerColor(kBlue+3*i);
    res_logpt_dijet_kfsrfit_var[i]->SetMarkerStyle(20+i);
    res_logpt_dijet_kfsrfit_var[i]->Draw("E1 SAME");  

  }
  TLegend leg5 = tdrLeg(0.20,0.19,0.43,0.42); 
  leg5 . AddEntry(res_logpt_dijet_kfsrfit_var[1] , "60 GeV","LP");  
  leg5 . AddEntry(res_logpt_dijet_kfsrfit_var[0] , "120 GeV","LP");
  leg5 . AddEntry(res_logpt_dijet_kfsrfit_var[2] , "240 GeV","LP"); 
  leg5 . AddEntry(res_logpt_dijet_kfsrfit_var[3] , "480 GeV","LP");
  leg5 . AddEntry(res_logpt_dijet_kfsrfit, "Nominal","LP");
  leg5.Draw();
  tex->DrawLatex(0.45,0.87,JetDescrib);                    
  tex->DrawLatex(0.15,0.95,txttag);
  c7->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_DiJet_kFSRfit_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");
  


  
 // TCanvas *c6 = new TCanvas();
 //  tdrCanvas(c6,"L2res_logpt_MPF_kFSRfit_ptDepend_diff",h,4,10,kSquare,CorrectionObject::_lumitag);
 //  h->GetYaxis()->SetTitle("new - prev. iteration");
 //  h->GetYaxis()->SetRangeUser(-0.3,0.3);
 //   res_logpt_mpf_diff->SetLineColor(kBlack);
 //  res_logpt_mpf_diff->Draw("E1 SAME");
 //  for(int i=0;i<4;i++){
 //    res_logpt_mpf_diff_var[i]->SetLineColor(kRed-3*i);
 //    res_logpt_mpf_diff_var[i]->SetMarkerColor(kRed-3*i);
 //    res_logpt_mpf_diff_var[i]->SetMarkerStyle(20+i);
 //    res_logpt_mpf_diff_var[i]->Draw("E1 SAME"); 
 //  }

 //  TLegend leg4 = tdrLeg(0.20,0.19,0.43,0.42); 
 //  leg4 . AddEntry(res_logpt_mpf_kfsrfit_var[1] , "60 GeV","LP");  
 //  leg4 . AddEntry(res_logpt_mpf_kfsrfit_var[0] , "120 GeV","LP");
 //  leg4 . AddEntry(res_logpt_mpf_kfsrfit_var[2] , "240 GeV","LP"); 
 //  leg4 . AddEntry(res_logpt_mpf_kfsrfit_var[3] , "480 GeV","LP");
 //  leg4 . AddEntry(res_logpt_mpf_kfsrfit, "Nominal","LP");
 //  leg4.Draw();     
 //  line2->Draw("SAME");
 //  tex->DrawLatex(0.45,0.87,JetDescrib);      
 //  c6->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_MPF_kFSRfit_diff"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");
  
  
 //  TCanvas *c8 = new TCanvas();
 //  tdrCanvas(c8,"L2res_logpt_DiJet_kFSRfit_ptDepend",h,4,10,kSquare,CorrectionObject::_lumitag);
 //  h->GetYaxis()->SetTitle("new - prev. iteration");
 //  h->GetYaxis()->SetRangeUser(-0.3,0.3);
 //  res_logpt_dijet_diff->SetLineColor(kBlack);
 //  res_logpt_dijet_diff->Draw("E1 SAME");
 //  for(int i=0;i<4;i++){
 //    res_logpt_dijet_diff_var[i]->SetLineColor(kBlue+3*i);
 //    res_logpt_dijet_diff_var[i]->SetMarkerColor(kBlue+3*i);
 //    res_logpt_dijet_diff_var[i]->SetMarkerStyle(20+i);
 //    res_logpt_dijet_diff_var[i]->Draw("E1 SAME");  

 //  }
 //  TLegend leg6 = tdrLeg(0.20,0.19,0.43,0.42); 
 //  leg6 . AddEntry(res_logpt_dijet_kfsrfit_var[1] , "60 GeV","LP");  
 //  leg6 . AddEntry(res_logpt_dijet_kfsrfit_var[0] , "120 GeV","LP");
 //  leg6 . AddEntry(res_logpt_dijet_kfsrfit_var[2] , "240 GeV","LP"); 
 //  leg6 . AddEntry(res_logpt_dijet_kfsrfit_var[3] , "480 GeV","LP");
 //  leg6 . AddEntry(res_logpt_dijet_kfsrfit, "Nominal","LP");
 //  leg6.Draw();
 //  line2->Draw("SAME");
 //  tex->DrawLatex(0.45,0.87,JetDescrib);                    
 //  c8->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_DiJet_kFSRfit_diff"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");
  
  TCanvas *c9 = new TCanvas();
  tdrCanvas(c9,"L2res_logpt_MPF_DiJet_kFSRfit_ptDepend",h,4,10,kSquare,CorrectionObject::_lumitag);
  h->GetYaxis()->SetTitle("MPF - pT-bal");
  h->GetYaxis()->SetRangeUser(-0.12,0.12);
  res_mpf_pT_diff->SetLineColor(kBlack);
  res_mpf_pT_diff->Draw("E1 SAME");
  for(int i=0;i<4;i++){
    res_mpf_pT_diff_var[i]->SetLineColor(kOrange+3*i);
    res_mpf_pT_diff_var[i]->SetMarkerColor(kOrange+3*i);
    res_mpf_pT_diff_var[i]->SetMarkerStyle(20+i);
    res_mpf_pT_diff_var[i]->Draw("E1 SAME");  

  }
  TLegend leg7 = tdrLeg(0.20,0.19,0.43,0.42); 
  leg7 . AddEntry(res_mpf_pT_diff_var[1] , "60 GeV","LP");  
  leg7 . AddEntry(res_mpf_pT_diff_var[0] , "120 GeV","LP");
  leg7 . AddEntry(res_mpf_pT_diff_var[2] , "240 GeV","LP"); 
  leg7 . AddEntry(res_mpf_pT_diff_var[3] , "480 GeV","LP");
  leg7 . AddEntry(res_mpf_pT_diff, "Nominal","LP");
  leg7.Draw();
  line2->Draw("SAME");
  tex->DrawLatex(0.45,0.87,JetDescrib);                    
  tex->DrawLatex(0.15,0.95,txttag);
  c9->SaveAs(CorrectionObject::_outpath+"plots/L2Res_logpt_MPF_DiJet_kFSRfit"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".pdf");
  



  cout << " ------------------- Finished L2Res(). Going to delete everything. ------------------------ " << endl;

  //delete everything
  
  delete c9;
  //  delete c8;
  delete c7;
  //  delete c6;
  delete c5;
  delete c4;
  delete c3;
  delete c2;
  delete c1;

 
  for(int i=0; i<4; i++){
    // delete res_logpt_mpf_kfsrfit_old_var[i];
    // delete res_logpt_dijet_kfsrfit_old_var[i];
    delete res_logpt_mpf_kfsrfit_var[i];
    delete res_logpt_dijet_kfsrfit_var[i];
    //delete res_logpt_mpf_diff_var[i];
    //delete res_logpt_dijet_diff_var[i];
   }


  delete f_Res_dijet_var;
  delete f_Res_mpf_var;
  // delete f_Res_mpf_old_var;
  // delete f_Res_dijet_old_var;

  //delete leg2;
  //delete leg1;
  delete tex_bars;
  delete tex;
  
  delete h;
  delete line;
  delete line2;
 

  delete res_logpt_dijet_kfsrfit;
  delete res_logpt_mpf_kfsrfit;
  delete res_const_dijet_kfsrfit;
  delete res_const_mpf_kfsrfit;
  delete kfsr_dijet_fit;
  delete kfsr_dijet;
  delete kfsr_mpf_fit;
  delete kfsr_mpf;
  delete pt_depend_logpt_dijet;
  delete pt_depend_const_dijet;
  delete pt_depend_logpt_mpf;
  delete pt_depend_const_mpf;
  delete f_Res_dijet;
  delete f_Res_mpf;

}
