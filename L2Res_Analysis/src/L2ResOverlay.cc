#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>

#include "../include/tdrstyle_mod15.h"






using namespace std;

//Function to overlay current results with another set of corrections
void CorrectionObject::L2ResOverlay(bool is_MPF){
  cout << "--------------- Starting L2ResOverlay() ---------------" << endl << endl;

  TCanvas* c1 = new TCanvas();
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h->SetMaximum(1.2); //1.2
  h->SetMinimum(0.8); //0.8
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  tdrCanvas(c1,"c1",h,4,10,true,CorrectionObject::_lumitag);

  TLegend leg1 = tdrLeg(0.17,0.14,0.35,0.45);
  TLine *line = new TLine(0.,1,5.191,1);              

  const int n_runs = 4;
  TString runnr_v[n_runs]={"BCD","EFearly","FlateG","H"};
  vector<TLegend> leg_var, leg_var_diff;
  for(int i=0; i<n_runs; i++){
    TLegend tmp = tdrLeg(0.17,0.14,0.35,0.45);
    leg_var.push_back(tmp);
    TLegend tmp_diff =  tdrLeg(0.17,0.14,0.35,0.45);
    leg_var_diff.push_back(tmp_diff);
  }
  TString names_var[4] = {"60 GeV", "120 GeV", "240 GeV", "480 GeV"};
  //TString runnr_v[n_runs]={"EFearly"};
  TCanvas* c_var[n_runs], *c_var_diff[n_runs];
  TCanvas* c2 = new TCanvas;
  tdrCanvas(c2,"c2",h,4,10,true,"V4 vs. V5");

  TFile* f_Res_nom;
  TFile* f_Res_var[4];
  TFile* f_Res_nom_old;
  TFile* f_Res_var_old[4];
  TH1D* res_kfsrfit_nom;
  TH1D* res_kfsrfit_var[4];
  TH1D* res_kfsrfit_nom_old;
  TH1D* res_kfsrfit_var_old[4];
  TH1D* res_diff_nom[n_runs];
  TH1D* res_diff_var[n_runs][4];
  for(int j=0;j<n_runs;j++){
    TString runnr=runnr_v[j];
    TString path, path_old;
    if(is_MPF){
      path = CorrectionObject::_outpath + "../Run" + runnr + "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag;
      path_old = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Summer16_23Sep2016_V1/AK4CHS/AveResponse/Run" + runnr + "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag;
    }
    else{
      path = CorrectionObject::_outpath + "../Run" + runnr + "/" + "Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag;
      path_old = "/nfs/dust/cms/user/reimersa/JEC/2016ReReco/Residuals/Summer16_23Sep2016_V1/AK4CHS/AveResponse/Run" + runnr + "/" + "Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag;
    }

    f_Res_nom = new TFile(path+".root","READ");   
    f_Res_nom_old = new TFile(path_old+".root","READ"); 
    f_Res_var[0] = new TFile(path+"_down.root","READ");   
    f_Res_var_old[0] = new TFile(path_old+"_down.root","READ");  
    f_Res_var[1] = new TFile(path+"_central.root","READ");   
    f_Res_var_old[1] = new TFile(path_old+"_central.root","READ");  
    f_Res_var[2] = new TFile(path+"_up.root","READ");   
    f_Res_var_old[2] = new TFile(path_old+"_up.root","READ");  
    f_Res_var[3] = new TFile(path+"_doubleup.root","READ");   
    f_Res_var_old[3] = new TFile(path_old+"_doubleup.root","READ");  
    
 
    if(is_MPF){
      res_kfsrfit_nom = (TH1D*)f_Res_nom->Get("res_logpt_mpf");
      res_kfsrfit_nom_old = (TH1D*)f_Res_nom_old->Get("res_logpt_mpf");
      for(int i=0; i<4; i++){
	res_kfsrfit_var[i] = (TH1D*)f_Res_var[i]->Get("res_logpt_mpf");
	res_kfsrfit_var_old[i] = (TH1D*)f_Res_var_old[i]->Get("res_logpt_mpf");
      }
    }
    else{
      res_kfsrfit_nom = (TH1D*)f_Res_nom->Get("res_logpt_dijet");
      res_kfsrfit_nom_old = (TH1D*)f_Res_nom_old->Get("res_logpt_dijet");
      for(int i=0; i<4; i++){
	res_kfsrfit_var[i] = (TH1D*)f_Res_var[i]->Get("res_logpt_dijet");
	res_kfsrfit_var_old[i] = (TH1D*)f_Res_var_old[i]->Get("res_logpt_dijet");
      }
    }

    res_diff_nom[j] = (TH1D*)res_kfsrfit_nom_old->Clone();
    res_diff_nom[j]->Add(res_kfsrfit_nom,-1);
    for(int i=0; i<4; i++){
      res_diff_var[j][i] = (TH1D*)res_kfsrfit_var_old[i]->Clone();
      res_diff_var[j][i]->Add(res_kfsrfit_var[i],-1);
    }

    

    //style for nominal
    res_kfsrfit_nom->SetMarkerStyle(1);
    res_kfsrfit_nom_old->SetMarkerStyle(1); 
    res_diff_nom[j]->SetMarkerStyle(1); 
    res_kfsrfit_nom_old->SetLineStyle(2); // dashed line 
    res_kfsrfit_nom->SetLineWidth(2);
    res_kfsrfit_nom_old->SetLineWidth(2);
    res_diff_nom[j]->SetLineWidth(2);
    res_kfsrfit_nom->SetLineColor(1+j); 
    res_kfsrfit_nom_old->SetLineColor(1+j); 
    res_diff_nom[j]->SetLineColor(1+j); 
    if(j == 2){
      res_kfsrfit_nom->SetLineColor(kGreen-2); 
      res_kfsrfit_nom_old->SetLineColor(kGreen-2);
      res_diff_nom[j]->SetLineColor(kGreen-2);
    }
    leg1.AddEntry(res_kfsrfit_nom,runnr+" V5","L");
    leg1.AddEntry(res_kfsrfit_nom_old,runnr+" V4","L");
    res_kfsrfit_nom->Draw("E SAME");  
    res_kfsrfit_nom_old->Draw("E SAME");  
  }
  line->Draw("SAME");
  leg1.Draw();



  TString save_as;
  if(is_MPF) save_as = CorrectionObject::_outpath + "../plots/L2Res_MPF_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V5_vs_CorrectFormulae_MCReweighting.pdf";
  else       save_as = CorrectionObject::_outpath + "../plots/L2Res_DiJet_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V5_vs_CorrectFormulae_MCReweighting.pdf";
  cout << "save-as = " << save_as << endl;
  c2->SaveAs(save_as); 

  for(int j=0;j<n_runs;j++){
    for(int i=0; i<4; i++){

      c_var[j] = new TCanvas();
      tdrCanvas(c_var[j],"c2"+runnr_v[j],h,4,10,true,"Run"+runnr_v[j]+" V4 vs. V5");

      //style for variations
      res_kfsrfit_var[i]->SetMarkerStyle(1);
      res_kfsrfit_var_old[i]->SetMarkerStyle(1); 
      res_diff_var[j][i]->SetMarkerStyle(20+i); 
      res_kfsrfit_var_old[i]->SetLineStyle(2); // dashed line 
      res_kfsrfit_var[i]->SetLineWidth(2);
      res_kfsrfit_var_old[i]->SetLineWidth(2);
      res_diff_var[j][i]->SetLineWidth(2);
      res_kfsrfit_var[i]->SetLineColor(1+j); 
      res_kfsrfit_var_old[i]->SetLineColor(1+j); 
      if(is_MPF){
	res_diff_var[j][i]->SetLineColor(kRed-3*i); 
	res_diff_var[j][i]->SetMarkerColor(kRed-3*i); 
      }
      else{
	res_diff_var[j][i]->SetLineColor(kBlue+3*i); 
	res_diff_var[j][i]->SetMarkerColor(kBlue+3*i); 
      }
      if(j == 2){
	res_kfsrfit_var[i]->SetLineColor(kGreen-2); 
	res_kfsrfit_var_old[i]->SetLineColor(kGreen-2);
      }
      leg_var[j].AddEntry(res_kfsrfit_nom,"nominal V5","L");
      leg_var[j].AddEntry(res_kfsrfit_nom_old,"nominal V4","L");
      leg_var[j].AddEntry(res_kfsrfit_var[i],names_var[i]+" V5","L");
      leg_var[j].AddEntry(res_kfsrfit_var_old[i],names_var[i]+" V4","L");
      res_kfsrfit_var[i]->Draw("E SAME");  
      res_kfsrfit_var_old[i]->Draw("E SAME");  
    }

    line->Draw("SAME");
    leg_var[j].Draw();

    TString save_as_var;
    if(is_MPF) save_as_var = CorrectionObject::_outpath + "../plots/L2Res_MPF_LOGLIN_PtDepend_Run" + runnr_v[j]  + "_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V5_vs_CorrectFormulae_MCReweighting.pdf";
    else       save_as_var = CorrectionObject::_outpath + "../plots/L2Res_DiJet_LOGLIN_PtDepend_Run" + runnr_v[j]  + "_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V5_vs_CorrectFormulae_MCReweighting.pdf";
    cout << "save-as_var = " << save_as_var << endl;
    c_var[j]->SaveAs(save_as_var); 
  }

  TCanvas* c3 = new TCanvas;
  h->SetMaximum(0.3); 
  h->SetMinimum(-0.3); 
  h->GetYaxis()->SetTitle("V4 - V5");
  tdrCanvas(c3,"c3",h,4,10,true,"V4 - V5");
  TLegend leg2 = tdrLeg(0.17,0.14,0.35,0.45);
  for(int i=0; i<n_runs; i++){
    TString runnr=runnr_v[i];
    leg2.AddEntry(res_diff_nom[i],runnr+": (V4-V5)","L");
    res_diff_nom[i]->Draw("E SAME");
  }
  line->Draw("SAME");
  leg2.Draw();

  TString save_as_diff;
  if(is_MPF) save_as_diff = CorrectionObject::_outpath + "../plots/L2Res_MPF_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V5_vs_CorrectFormulae_MCReweighting_diff.pdf";
  else       save_as_diff = CorrectionObject::_outpath + "../plots/L2Res_DiJet_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V5_vs_CorrectFormulae_MCReweighting_diff.pdf";
  cout << "save-as_diff = " << save_as_diff << endl;
  c3->SaveAs(save_as_diff); 



  for(int j=0;j<n_runs;j++){

    c_var_diff[j] = new TCanvas();
    tdrCanvas(c_var_diff[j],"c3"+runnr_v[j],h,4,10,true,"Run"+runnr_v[j]+" V4 - V5");
    res_diff_nom[j]->SetLineColor(kBlack);
    res_diff_nom[j]->Draw("E SAME");
    leg_var_diff[j].AddEntry(res_diff_nom[j],"nominal (V4-V5)","L");
    for(int i=0; i<4; i++){
      leg_var_diff[j].AddEntry(res_diff_var[j][i],names_var[i]+" (V4-V5)","LP");
      res_diff_var[j][i]->Draw("E SAME");
    }


    line->Draw("SAME");
    leg_var_diff[j].Draw();

    TString save_as_var_diff;
    if(is_MPF) save_as_var_diff = CorrectionObject::_outpath + "../plots/L2Res_MPF_LOGLIN_PtDepend_Run" + runnr_v[j]  + "_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V5_vs_CorrectFormulae_MCReweighting_diff.pdf";
    else       save_as_var_diff = CorrectionObject::_outpath + "../plots/L2Res_DiJet_LOGLIN_PtDepend_Run" + runnr_v[j]  + "_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_VersionComparison_Summer16_23Sep2016_V5_vs_CorrectFormulae_MCReweighting_diff.pdf";
    cout << "save-as_var_diff = " << save_as_var_diff << endl;
    c_var_diff[j]->SaveAs(save_as_var_diff); 
  }




  //delete everything
  for(int i=0; i<n_runs; i++) delete res_diff_nom[i];
  for(int j=0; j<n_runs; j++){
    for(int i=0; i<4; i++){
      delete res_diff_var[j][i];
    }
  }
  delete res_kfsrfit_nom;
  delete res_kfsrfit_nom_old;
  for(int i=0; i<4; i++){
    delete res_kfsrfit_var[i];
    delete res_kfsrfit_var_old[i];
  }
  delete f_Res_nom;
  delete f_Res_nom_old;
  for(int i=0; i<4; i++){
    delete f_Res_var[i];
    delete f_Res_var_old[i];
  }
  delete line;
  delete h;
  delete c1;


}
