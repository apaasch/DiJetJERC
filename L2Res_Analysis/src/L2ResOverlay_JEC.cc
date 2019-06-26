#include "../include/CorrectionObject.h"
//#include "../include/parameters.h"
#include "../../include/constants.h"
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
void CorrectionObject::L2ResOverlay_JEC(){
  cout << "--------------- Starting L2ResOverlay_JEC() ---------------" << endl << endl;

  const int n_runnr = 4;
  TString runnr_v[n_runnr]={"B","C","DE","F"};
  TString legend_v[n_runnr] ={"Run2017B","Run2017C","Run2017DE","Run2017F"};

  TFile* f_Res_mpf_nom;
  TFile* f_Res_mpf_down;
  TFile* f_Res_mpf_up;

  TFile* f_Res_rel_nom;
  TFile* f_Res_rel_down;
  TFile* f_Res_rel_up;

  TH1D* res_mpf_kfsrfit_nom;
  TH1D* res_mpf_kfsrfit_down;
  TH1D* res_mpf_kfsrfit_up;

  TH1D* res_rel_kfsrfit_nom;
  TH1D* res_rel_kfsrfit_down;
  TH1D* res_rel_kfsrfit_up;
 bool kSquare = true;

  for(int j=0;j<n_runnr;j++){
    TString runnr=runnr_v[j];
    TString path;

    if(runnr=="D"||runnr=="E"||runnr=="DE"||runnr=="F"||runnr=="DEF") path=CorrectionObject::_input_pathDi+TString("/Run") + runnr +"_"+CorrectionObject::_dname;
    else if(runnr=="B"||runnr=="C"||runnr=="BC") path =CorrectionObject::_input_pathSi+"/Run" + runnr +"_"+CorrectionObject::_dnameSi;
    else throw runtime_error("In L2ResAllOverlay_JEC(): could not identify if it is a DiJet or single Jet run.");    
  
    // std::cout<< "looking in path: "<<path<<std::endl;
    
    f_Res_mpf_nom  = new TFile(path +"_JERnominalV2_RunBCDEF"+ "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");   
    f_Res_mpf_down = new TFile(path +"_JERupV2"+ "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ"); 
    f_Res_mpf_up     = new TFile(path +"_JERdownV2"+ "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");

     f_Res_rel_nom  = new TFile(path +"_JERnominalV2_RunBCDEF"+ "/" + "Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");   
    f_Res_rel_down = new TFile(path +"_JERupV2"+ "/" + "Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ"); 
    f_Res_rel_up     = new TFile(path +"_JERdownV2"+ "/" + "Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root","READ");   
    
 
    res_mpf_kfsrfit_nom  = (TH1D*)f_Res_mpf_nom->Get("res_logpt_mpf");
    res_mpf_kfsrfit_down = (TH1D*)f_Res_mpf_down->Get("res_logpt_mpf");
    res_mpf_kfsrfit_up      = (TH1D*)f_Res_mpf_up  ->Get("res_logpt_mpf");
 
    res_rel_kfsrfit_nom  = (TH1D*)f_Res_rel_nom->Get("res_logpt_dijet");
    res_rel_kfsrfit_down = (TH1D*)f_Res_rel_down->Get("res_logpt_dijet");
    res_rel_kfsrfit_up      = (TH1D*)f_Res_rel_up  ->Get("res_logpt_dijet");
 
    TString lumitag;
    //lumitags
    if(runnr == "BC") lumitag      = "RunBC 14.4 fb^{-1}";   //2017! 
    //FIXME differentiate between 2016 (5.8 fb^{-1}) and 2017
    else if(runnr == "B") lumitag      = "RunB  4.8 fb^{-1}";
    else if(runnr == "C") lumitag      = "RunC  9.6 fb^{-1}";
    else if(runnr == "D") lumitag      = "RunD  4.2 fb^{-1}";
    else if(runnr == "E") lumitag      = "RunE  9.3 fb^{-1}";
    else if(runnr == "F") lumitag      = "RunF  13.5 fb^{-1}";
    else if(runnr == "DE") lumitag      = "RunDE  13.5 fb^{-1}"; 
    else if(runnr == "DEF") lumitag      = "RunDEF  26.9 fb^{-1}"; 
    else if(runnr == "BCD") lumitag      = "RunBCD  18.6 fb^{-1}";
    else if(runnr == "BCDEF") lumitag      = "RunBCDEF  41.3 fb^{-1}";
    else if(runnr == "CDEF") lumitag      = "RunCDEF 36.5 fb^{-1}";
    else throw runtime_error("In constructor: Invalid RunNr. specified.");
    
      
    TCanvas* c3 = new TCanvas("c1","c1",600,600);
    TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
    h->SetMaximum(1.2); //1.2
    h->SetMinimum(0.8); //0.8
    tdrCanvas(c3,"c3",h,4,10,kSquare,lumitag);
     
    TLegend leg1 = tdrLeg(0.20,0.18,0.35,0.40);
    TLine *line = new TLine(0.,1,5.191,1);              
    
    //style for mpf
    res_mpf_kfsrfit_nom  ->SetMarkerStyle(1);
    res_mpf_kfsrfit_down ->SetMarkerStyle(1); 
    res_mpf_kfsrfit_up     ->SetMarkerStyle(1); 

    res_mpf_kfsrfit_nom->SetLineWidth(2); 
    res_mpf_kfsrfit_down->SetLineWidth(2);
    res_mpf_kfsrfit_up   ->SetLineWidth(2);

    res_mpf_kfsrfit_nom->SetLineColor(kBlack); 
    res_mpf_kfsrfit_down->SetLineColor(kRed); 
    res_mpf_kfsrfit_up->SetLineColor(kBlue); 

    leg1.SetHeader("MPF");
    leg1.AddEntry(res_mpf_kfsrfit_nom, "nominal");
    leg1.AddEntry(res_mpf_kfsrfit_down, "down");
    leg1.AddEntry(res_mpf_kfsrfit_up, "up");
    leg1.Draw();

    res_mpf_kfsrfit_nom->Draw("E1 SAME"); 
    res_mpf_kfsrfit_down->Draw("E1 SAME"); 
    res_mpf_kfsrfit_up->Draw("E1 SAME"); 

    line->Draw("E1 SAME");
    
    // std::cout<< "saving in path: "<<path<<std::endl;

    c3->SaveAs(path+"_JERnominalV2_RunBCDEF/plots/L2Res_MPF_JER_summary.pdf");
   delete c3;
    delete h;
 

    TLegend leg2 = tdrLeg(0.20,0.18,0.35,0.40);

    TCanvas* c2 = new TCanvas();
    TH1D *t = new TH1D("t",";|#eta|;Relative correction",41,0,5.191);
    t->SetMaximum(1.2); //1.2
    t->SetMinimum(0.8); //0.8
     tdrCanvas(c2,"c2",t,4,10,kSquare,lumitag);

    res_rel_kfsrfit_nom  ->SetMarkerStyle(1);
    res_rel_kfsrfit_down ->SetMarkerStyle(1); 
    res_rel_kfsrfit_up     ->SetMarkerStyle(1); 

    res_rel_kfsrfit_nom->SetLineWidth(2); 
    res_rel_kfsrfit_down->SetLineWidth(2);
    res_rel_kfsrfit_up   ->SetLineWidth(2);

    res_rel_kfsrfit_nom->SetLineColor(kBlack); 
    res_rel_kfsrfit_down->SetLineColor(kRed); 
    res_rel_kfsrfit_up->SetLineColor(kBlue); 
    
    leg2.SetHeader("p_{T}-bal");
    leg2.AddEntry(res_rel_kfsrfit_nom, "nominal");
    leg2.AddEntry(res_rel_kfsrfit_down, "down");
    leg2.AddEntry(res_rel_kfsrfit_up, "up");
    leg2.Draw();

    res_rel_kfsrfit_nom->Draw("E1 SAME"); 
    res_rel_kfsrfit_down->Draw("E1 SAME"); 
    res_rel_kfsrfit_up->Draw("E1 SAME"); 
    line->Draw("E1 SAME");

    c2->SaveAs(path+"_JERnominalV2_RunBCDEF/plots/L2Res_DiJet_JER_summary.pdf");

    delete c2;
    delete t;
  }
}
