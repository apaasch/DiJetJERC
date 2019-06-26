#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TROOT.h>


using namespace std;

void CorrectionObject::kFSR_CorrectFormulae(float pt_min){
  cout << "--------------- Starting kFSR() ---------------" << endl << endl;
  TStyle* m_gStyle = new TStyle("beautiful","style to make plots more beautiful");
  m_gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2(kTRUE);
  m_gStyle->SetOptStat(kFALSE);
  m_gStyle->SetPalette(1);
  m_gStyle->SetTitleSize(0.05,"x");
  m_gStyle->SetTitleSize(0.05,"y");
  m_gStyle->SetTitleSize(0.05,"z");
  m_gStyle->SetLabelSize(0.05,"x");
  m_gStyle->SetLabelSize(0.05,"y");
  m_gStyle->SetLabelSize(0.05,"z");
  m_gStyle->SetTitleYOffset(1.12);
  m_gStyle->SetPadLeftMargin(0.15);
  m_gStyle->SetPadRightMargin(0.15);
  m_gStyle->SetPadColor(0);
  m_gStyle->SetCanvasColor(0);
  gROOT->SetStyle("beautiful");
  int n_pt_ = max(n_pt,n_pt_HF);
  bool eta_cut_bool;
  int n_pt_cutted;

 
  //get MC response
  double rel_r_mc[n_pt_-1][n_eta-1][n_alpha];
  double err_rel_r_mc[n_pt_-1][n_eta-1][n_alpha];
 
  double mpf_r_mc[n_pt_-1][n_eta-1][n_alpha];
  double err_mpf_r_mc[n_pt_-1][n_eta-1][n_alpha];
 
  //get DATA response
  double rel_r_data[n_pt_-1][n_eta-1][n_alpha];
  double err_rel_r_data[n_pt_-1][n_eta-1][n_alpha];
 
  double mpf_r_data[n_pt_-1][n_eta-1][n_alpha];
  double err_mpf_r_data[n_pt_-1][n_eta-1][n_alpha];
 
  
// get ratio for MC to DATA responses
  double ratio_al_rel_r[n_pt_-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al_rel_r[n_pt_-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins
  double ratio_al_mpf_r[n_pt_-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf_r[n_pt_-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins

  TProfile *pr_data_asymmetry[n_eta-1][n_alpha];// pT-balance response for data  
  TProfile *pr_data_B[n_eta-1][n_alpha];//MPF response for data
  TProfile *pr_mc_asymmetry[n_eta-1][n_alpha];// pT-balanse responce for MC  
  TProfile *pr_mc_B[n_eta-1][n_alpha];//MPF response for MC

  TH2D *hdata_asymmetry[n_eta-1][n_alpha];
  TH2D *hdata_B[n_eta-1][n_alpha];
  TH2D *hmc_asymmetry[n_eta-1][n_alpha];
  TH2D *hmc_B[n_eta-1][n_alpha];

  int n_entries_mc[n_eta-1][n_alpha][n_pt_-1];
  int n_entries_data[n_eta-1][n_alpha][n_pt_-1];
  int count = 0;

  TString name1 = "hist_data_asymmetry_";
  TString name2 = "hist_data_B_";
  TString name3 = "hist_mc_asymmetry_";
  TString name4 = "hist_mc_B_";

  for(int i=0; i<n_alpha; i++){
    for(int j=0; j<n_eta-1; j++){
      eta_cut_bool = fabs(eta_bins[j])>eta_cut;
      for(int k= 0 ; k < n_pt_-1  ; k++ ){
	rel_r_mc[k][j][i] = 0;      
	err_rel_r_mc[k][j][i] = 0;  
	mpf_r_mc[k][j][i] = 0;      
	err_mpf_r_mc[k][j][i] = 0;  
                                  
	rel_r_data[k][j][i] = 0;    
	err_rel_r_data[k][j][i] = 0;
	mpf_r_data[k][j][i] = 0;    
	err_mpf_r_data[k][j][i] = 0;
	
	ratio_al_rel_r[k][j][i] = 0;
	err_ratio_al_rel_r[k][j][i] = 0;
	ratio_al_mpf_r[k][j][i] = 0;
	err_ratio_al_mpf_r[k][j][i] = 0;
      }
      
      n_pt_cutted = ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 );
      TString name = name1; name+=count;
      hdata_asymmetry[j][i] = new TH2D(name,"A in DATA; p_{T}^{ave} [GeV]; A",n_pt_cutted , (eta_cut_bool?pt_bins_HF:pt_bins),nResponseBins, Response_min, Response_max);
      name = name2;name+=count;
      hdata_B[j][i]         = new TH2D(name,"B in DATA;p_{T}^{ave} [GeV];B",n_pt_cutted , (eta_cut_bool?pt_bins_HF:pt_bins),nResponseBins, Response_min, Response_max);
      name = name3; name+=count;
      hmc_asymmetry[j][i]   = new TH2D(name,"A in MC;p_{T}^{ave} [GeV];A",n_pt_cutted , (eta_cut_bool?pt_bins_HF:pt_bins),nResponseBins, Response_min, Response_max);
      name = name4; name+=count;
      hmc_B[j][i]           = new TH2D(name,"B in MC;p_{T}^{ave} [GeV];B",n_pt_cutted , (eta_cut_bool?pt_bins_HF:pt_bins),nResponseBins, Response_min, Response_max);
          
      count++;
    }
  }
  
  cout << "Set up a total of " << count << " histograms." << endl;

    for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
      for(int k= 0 ; k < n_pt_-1 ; k++ ){
	n_entries_mc[j][i][k] = 0;
	n_entries_data[j][i][k] = 0;
      }
     }
   }

  cout<<"Define Tree Readers.\n";
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
  TTreeReaderValue<Float_t> asymmetry_data(myReader_DATA, "asymmetry");
  TTreeReaderValue<Float_t> B_data(myReader_DATA, "B");   
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  TTreeReaderValue<Float_t> jet3_pt_data(myReader_DATA, "jet3_pt");
  int idx = 0;
  
  cout << "starting to loop over DATA events." << endl;

  while (myReader_DATA.Next()) {

    if(*alpha_data<0.05) continue; 

    for(int j=0; j<n_eta-1; j++){
      if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
      for(int i=0; i<n_alpha; i++){
	if(*alpha_data>alpha_bins[i]) continue;
	hdata_asymmetry[j][i]->Fill(*pt_ave_data,*asymmetry_data,*weight_data);
	hdata_B[j][i]->Fill(*pt_ave_data,*B_data,*weight_data);
	eta_cut_bool = fabs(eta_bins[j])>eta_cut;
	for(int k= 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){
	      if((*pt_ave_data < (eta_cut_bool?pt_bins_HF:pt_bins)[k]) || (*pt_ave_data >= (eta_cut_bool?pt_bins_HF:pt_bins)[k+1])) continue;
	      n_entries_data[j][i][k]++;
	}
	idx++;
	if(idx%5000000==0) cout << "looping over data-TTree: Idx = " << idx << endl;
      }
    }
  }
  cout << "Finished running over DATA events. Read in a total of " << idx << " events." << endl;

   // Get relevant quantities from MC, loop over events
   TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
   TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
   TTreeReaderValue<Float_t> asymmetry_mc(myReader_MC, "asymmetry");
   TTreeReaderValue<Float_t> B_mc(myReader_MC, "B");
   TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
  TTreeReaderValue<Float_t> jet3_pt_mc(myReader_MC, "jet3_pt");
   idx = 0;

   //TODO make a option to load MC rel responses instead of calculating them every time

   while (myReader_MC.Next()) {

    if(*alpha_mc<0.05) continue; 
     //     if(idx>1000) break;//for tetsts

     for(int j=0; j<n_eta-1; j++){
       if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
       for(int i=0; i<n_alpha; i++){
	 if(*alpha_mc>alpha_bins[i]) continue;
	 hmc_asymmetry[j][i]->Fill(*pt_ave_mc,*asymmetry_mc,*weight_mc);
	 hmc_B[j][i]->Fill(*pt_ave_mc,*B_mc,*weight_mc);
	 eta_cut_bool = fabs(eta_bins[j])>eta_cut;
	 for(int k= 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){  
	   if((*pt_ave_mc < (eta_cut_bool?pt_bins_HF:pt_bins)[k]) || (*pt_ave_mc >= (eta_cut_bool?pt_bins_HF:pt_bins)[k+1])) continue;            
	   n_entries_mc[j][i][k]++;
	 }
	 idx++;
	   if(idx%10000000==0) cout << "looping over MC-TTree: Idx = " << idx << endl;
       }
     }
   }

   //Check number of entries in each bin
   bool enough_entries[n_alpha][n_eta-1][n_pt-1];
   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
	 eta_cut_bool = fabs(eta_bins[j])>eta_cut;
	 for(int k= 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){  
	   enough_entries[i][j][k] = false;
	   if(n_entries_mc[j][i][k] > 50 && n_entries_data[j][i][k] > 50) enough_entries[i][j][k] = true;
	}
     }
   }
   //build profiles out of asymmetry and B 2d-histos to get <A> and <B> as a function of pT in bins of eta,alpha
   for(int j=0; j<n_eta-1; j++){
     for(int i=0; i<n_alpha; i++){
       //print TH2D histos
       TCanvas* c_dummy1 = new TCanvas();
       hdata_asymmetry[j][i]->Draw("COLZ");
       delete c_dummy1;
       TCanvas* c_dummy2 = new TCanvas();
       hmc_asymmetry[j][i]->Draw("COLZ");
       delete c_dummy2;
       TCanvas* c_dummy3 = new TCanvas();
       hdata_B[j][i]->Draw("COLZ");      
       delete c_dummy3;
       TCanvas* c_dummy4 = new TCanvas();
       hmc_B[j][i]->Draw("COLZ");
       delete c_dummy4;

       //build profiles
       pr_data_asymmetry[j][i] = (TProfile*)hdata_asymmetry[j][i]->ProfileX();
       pr_data_B[j][i]         = (TProfile*)hdata_B[j][i]->ProfileX();
       pr_mc_asymmetry[j][i]   = (TProfile*)hmc_asymmetry[j][i]->ProfileX();
       pr_mc_B[j][i]           = (TProfile*)hmc_B[j][i]->ProfileX();

       //print profiles
       TCanvas* c_dummy5 = new TCanvas();
       pr_data_asymmetry[j][i]->Draw();
       c_dummy5->SetLogx();
       pr_data_asymmetry[j][i]->GetYaxis()->SetTitle("<A>");
       pr_data_asymmetry[j][i]->SetLineWidth(2);
       pr_data_asymmetry[j][i]->SetMinimum(-0.3);
       pr_data_asymmetry[j][i]->SetMaximum(0.3);
       delete c_dummy5;

       TCanvas* c_dummy6 = new TCanvas();
       pr_mc_asymmetry[j][i]->Draw();
       c_dummy6->SetLogx();
       pr_mc_asymmetry[j][i]->GetYaxis()->SetTitle("<A>");
       pr_mc_asymmetry[j][i]->SetLineWidth(2);
       pr_mc_asymmetry[j][i]->SetMinimum(-0.3);
       pr_mc_asymmetry[j][i]->SetMaximum(0.3);
       delete c_dummy6;

       TCanvas* c_dummy7 = new TCanvas();
       pr_data_B[j][i]->Draw();
       c_dummy7->SetLogx();
       pr_data_B[j][i]->GetYaxis()->SetTitle("<B>");
       pr_data_B[j][i]->SetLineWidth(2);
       pr_data_B[j][i]->SetMinimum(-0.3);
       pr_data_B[j][i]->SetMaximum(0.3);
       delete c_dummy7;

       TCanvas* c_dummy8 = new TCanvas();
       pr_mc_B[j][i]->Draw();
       c_dummy8->SetLogx();
       pr_mc_B[j][i]->GetYaxis()->SetTitle("<B>");
       pr_mc_B[j][i]->SetLineWidth(2);
       pr_mc_B[j][i]->SetMinimum(-0.3);
       pr_mc_B[j][i]->SetMaximum(0.3);
       delete c_dummy8;
     }
   }
   
   //calculate response from <A> and <B> in bins of pt,eta,alpha
   //gaussian error propagation from errors on <A> and <B>
   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
       eta_cut_bool = fabs(eta_bins[j])>eta_cut;
       for(int k= 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){ 
	 //responses for data, MC separately. Only for bins with >= 100 entries
	 double mpf_mc = (1+pr_mc_B[j][i]->GetBinContent(k+1))/(1-pr_mc_B[j][i]->GetBinContent(k+1));
	 if(!enough_entries[i][j][k]) mpf_mc = 0;
	 double mpf_data = (1+pr_data_B[j][i]->GetBinContent(k+1))/(1-pr_data_B[j][i]->GetBinContent(k+1));

	 if(!enough_entries[i][j][k]) mpf_data = 0;
	 double rel_mc = (1+pr_mc_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_mc_asymmetry[j][i]->GetBinContent(k+1));

	 if(!enough_entries[i][j][k]) rel_mc = 0;
	 double rel_data = (1+pr_data_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_data_asymmetry[j][i]->GetBinContent(k+1));

	 if(!enough_entries[i][j][k]) rel_data = 0;
	 double err_mpf_mc = 2/(pow((1-pr_mc_B[j][i]->GetBinContent(k+1)),2)) * pr_mc_B[j][i]->GetBinError(k+1);
	 double err_mpf_data = 2/(pow((1-pr_data_B[j][i]->GetBinContent(k+1)),2)) * pr_data_B[j][i]->GetBinError(k+1);
	 double err_rel_mc = 2/(pow((1-pr_mc_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_mc_asymmetry[j][i]->GetBinError(k+1);
	 double err_rel_data = 2/(pow((1-pr_data_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_data_asymmetry[j][i]->GetBinError(k+1);


	 //MC responses
	 rel_r_mc[k][j][i] = rel_mc;
	 err_rel_r_mc[k][j][i] = err_rel_mc;
	 mpf_r_mc[k][j][i] = mpf_mc;
	 err_mpf_r_mc[k][j][i] = err_mpf_mc;
	 
	 //DATA responses
	 rel_r_data[k][j][i] = rel_data;
	 err_rel_r_data[k][j][i] = err_rel_data;
	 mpf_r_data[k][j][i] = mpf_data;
	 err_mpf_r_data[k][j][i] = err_mpf_data;

	 
	 //ratio of responses, again gaussian error propagation
	 if(rel_data > 0){
	   ratio_al_rel_r[k][j][i] = rel_mc/rel_data;
	     }
	 else ratio_al_rel_r[k][j][i] = 0;
	 err_ratio_al_rel_r[k][j][i] = sqrt(pow(1/rel_data*err_rel_mc,2) + pow(rel_mc/(rel_data*rel_data)*err_rel_data,2));
	 if(mpf_data > 0) ratio_al_mpf_r[k][j][i] = mpf_mc/mpf_data;
	 else ratio_al_mpf_r[k][j][i] = 0;
	 err_ratio_al_mpf_r[k][j][i] = sqrt(pow(1/mpf_data*err_mpf_mc,2) + pow(mpf_mc/(mpf_data*mpf_data)*err_mpf_data,2));
       }
     }
   }
 
   //Normalization of hists to value at alpha = alpha_cut
 
  //1) find bin with alpha = alpha_cut: bin no. al_ref
   int al_ref=0;
   for(int i=0; i<n_alpha; i++){
     if(fabs(alpha_bins[i]-alpha_cut)<1e-4) 
       al_ref=i;
   }

   //2) Normalize values and errors of responses to value at alpha = alpha_cut
   
   for(int j=0; j<n_eta-1; j++){
     eta_cut_bool = fabs(eta_bins[j])>eta_cut;
     for(int k= 0 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){ 
       double norm_alref_rel_r = ratio_al_rel_r[k][j][al_ref];
       double err_norm_alref_rel_r = err_ratio_al_rel_r[k][j][al_ref];
       double norm_alref_mpf_r = ratio_al_mpf_r[k][j][al_ref];
       double err_norm_alref_mpf_r = err_ratio_al_mpf_r[k][j][al_ref];
       for(int i=0; i<n_alpha; i++){

	 if(norm_alref_rel_r>0){ //FIXME WHAT IS HAPPENING HERE? NO PROPER ERROR PROPAGATION !?! Ask other group that does kFSR-Extrapolation about error propagation
	   cout<<"--- HEY ref eta="<<eta_bins[j]<<" pt="<<pt_bins[k]<<" Rmc/Rdata = "<<norm_alref_rel_r<<" ratio_al_rel_r(alpha="<<alpha_bins[i]<<") = "<<ratio_al_rel_r[k][j][i]<<endl;
	   ratio_al_rel_r[k][j][i] =   ratio_al_rel_r[k][j][i]/norm_alref_rel_r; //original
	   err_ratio_al_rel_r[k][j][i] = sqrt(abs(pow(err_ratio_al_rel_r[k][j][i],2)-pow(err_norm_alref_rel_r,2)));

	 }
	 if(norm_alref_mpf_r>0){
	   ratio_al_mpf_r[k][j][i] =   ratio_al_mpf_r[k][j][i]/norm_alref_mpf_r;
	   err_ratio_al_mpf_r[k][j][i] = sqrt(abs(pow(err_ratio_al_mpf_r[k][j][i],2)-pow(err_norm_alref_mpf_r,2)));
	 }
       }
     }
   }
   
   // Build the Multigraphs containing the responses (MC and DATA) as a function of alpha
   TGraphErrors *graph_rel_r_mc[n_pt_-1][n_eta-1];  //set of points vs alpha
   TGraphErrors *graph_mpf_r_mc[n_pt_-1][n_eta-1];  //set of points vs alpha

   TGraphErrors *graph_rel_r_data[n_pt_-1][n_eta-1];  //set of points vs alpha
   TGraphErrors *graph_mpf_r_data[n_pt_-1][n_eta-1];  //set of points vs alpha
 
   // Build the Multigraphs containing the ratio of responses (MC/DATA) as a function of alpha
   TGraphErrors *graph_rel_r[n_pt_-1][n_eta-1];  //set of points vs alpha
   TMultiGraph *pTgraph_rel_r[n_eta-1];         //set of different pT bins in on eta bin
   TGraphErrors *graph_mpf_r[n_pt_-1][n_eta-1];  //set of points vs alpha
   TMultiGraph *pTgraph_mpf_r[n_eta-1];         //set of different pT bins in on eta bin


   //Define legend
   TLegend *leg1;
   leg1 = new TLegend(0.17,0.68,0.65,0.89,"","brNDC");//x+0.1
   leg1->SetBorderSize(0);
   leg1->SetTextSize(0.03);
   leg1->SetFillColor(10);
   leg1->SetLineColor(1);
   leg1->SetTextFont(42);
   leg1->SetNColumns(2);

   TLegend *leg2;//legend without points
   leg2 = new TLegend(0.37,0.78,0.70,0.89,"","brNDC");//x+0.1
   leg2->SetBorderSize(0);
   leg2->SetTextSize(0.045);
   leg2->SetFillColor(10);
   leg2->SetLineColor(1);
   leg2->SetTextFont(42);
   leg2->SetNColumns(2);

   //Define legend HF
   TLegend *leg1_HF;
   leg1_HF = new TLegend(0.17,0.68,0.65,0.89,"","brNDC");//x+0.1
   leg1_HF->SetBorderSize(0);
   leg1_HF->SetTextSize(0.03);
   leg1_HF->SetFillColor(10);
   leg1_HF->SetLineColor(1);
   leg1_HF->SetTextFont(42);
   leg1_HF->SetNColumns(2);

    
   //dummy for tdrCanvas
   TH1D *h = new TH1D("h",";dummy;",41,0,5.191);
   h->SetMaximum(1.2);
   h->SetMinimum(0.8);  

   
   double xbin_tgraph[n_alpha],zero[n_alpha];
   for(int i=0;i<n_alpha;i++){
     xbin_tgraph[i] = alpha_bins[i];
     zero[i] = 0;
   }

    TCanvas* c_rel_r[n_eta-1][n_pt_-1];
    TCanvas* c_mpf_r[n_eta-1][n_pt_-1];
    TString name_rel_r[n_eta-1][n_pt_-1];
    TString name_mpf_r[n_eta-1][n_pt_-1];
 
    for(int j=0; j<n_eta-1; j++){
      eta_cut_bool = fabs(eta_bins[j])>eta_cut;
      for(int k = 1 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){ 
        graph_rel_r_mc[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,rel_r_mc[k][j],zero,err_rel_r_mc[k][j]);
        graph_rel_r_mc[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_rel_r_mc[k][j]);
        graph_mpf_r_mc[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,mpf_r_mc[k][j],zero,err_mpf_r_mc[k][j]);
        graph_mpf_r_mc[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_r_mc[k][j]);
 
        graph_rel_r_data[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,rel_r_data[k][j],zero,err_rel_r_data[k][j]);
        graph_rel_r_data[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_rel_r_data[k][j]);
        graph_mpf_r_data[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,mpf_r_data[k][j],zero,err_mpf_r_data[k][j]);
        graph_mpf_r_data[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_r_data[k][j]);
 
        graph_rel_r_mc[k][j]->SetMarkerSize(1.5);
        graph_rel_r_mc[k][j]->SetMarkerColor(kRed);
        graph_rel_r_mc[k][j]->SetLineColor(kRed);
        graph_rel_r_mc[k][j]->SetMarkerStyle(20);
        graph_rel_r_mc[k][j]->SetTitle("");
 
        graph_mpf_r_mc[k][j]->SetMarkerSize(1.5);
        graph_mpf_r_mc[k][j]->SetMarkerColor(kRed);
        graph_mpf_r_mc[k][j]->SetLineColor(kRed);
        graph_mpf_r_mc[k][j]->SetMarkerStyle(20);
        graph_mpf_r_mc[k][j]->SetTitle("");
 
        graph_rel_r_data[k][j]->SetMarkerSize(1.5);
        graph_rel_r_data[k][j]->SetMarkerColor(kBlack);
        graph_rel_r_data[k][j]->SetLineColor(kBlack);
        graph_rel_r_data[k][j]->SetMarkerStyle(20);
        graph_rel_r_data[k][j]->SetTitle("");
 
        graph_mpf_r_data[k][j]->SetMarkerSize(1.5);
        graph_mpf_r_data[k][j]->SetMarkerColor(kBlack);
        graph_mpf_r_data[k][j]->SetLineColor(kBlack);
        graph_mpf_r_data[k][j]->SetMarkerStyle(20);
        graph_mpf_r_data[k][j]->SetTitle("");
 
        name_rel_r[j][k]="rel_r_alpha_"+eta_range[j]+"_"+eta_range[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1];
        c_rel_r[j][k] = new TCanvas(name_rel_r[j][k], name_rel_r[j][k], 850,700);
 
        m_gStyle->SetOptTitle(0);
  
        graph_rel_r_mc[k][j]->Draw("AP");
        graph_rel_r_data[k][j]->Draw("P SAME");
 
        graph_rel_r_mc[k][j]->GetYaxis()->SetRangeUser(0.9,1.2);
        graph_rel_r_mc[k][j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
        graph_rel_r_mc[k][j]->GetYaxis()->SetTitle("R_{p_{T}-balance}");
        graph_rel_r_mc[k][j]->GetYaxis()->SetTitleSize(0.045);
        graph_rel_r_mc[k][j]->GetYaxis()->SetTitleOffset(1.);
        graph_rel_r_mc[k][j]->GetXaxis()->SetTitle("cut on #alpha");
        graph_rel_r_mc[k][j]->GetXaxis()->SetTitleSize(0.045);
 
        graph_rel_r_data[k][j]->GetYaxis()->SetRangeUser(0.9,1.2);
        graph_rel_r_data[k][j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
        graph_rel_r_data[k][j]->GetYaxis()->SetTitle("R_{p_{T}-balance}");
        graph_rel_r_data[k][j]->GetYaxis()->SetTitleSize(0.045);
        graph_rel_r_data[k][j]->GetYaxis()->SetTitleOffset(1.);
        graph_rel_r_data[k][j]->GetXaxis()->SetTitle("cut on #alpha");
        graph_rel_r_data[k][j]->GetXaxis()->SetTitleSize(0.045);
 
        
        TLegend* leg_rel = new TLegend(0.20,0.15,0.31,0.35,"","brNDC");//x+0.1
        leg_rel->SetBorderSize(0);
        leg_rel->SetTextSize(0.038);
        leg_rel->SetFillColor(10);
        leg_rel->SetFillStyle(0);
        leg_rel->SetLineColor(1);
        leg_rel->SetTextFont(42);
        leg_rel->SetHeader("Rel response, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]+", "+(eta_cut_bool?pt_range_HF:pt_range)[k]+"#leq p_{T}<"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]);
        leg_rel->AddEntry(graph_rel_r_mc[k][j], "R^{MC}","P");
        leg_rel->AddEntry(graph_rel_r_data[k][j], "R^{DATA}","P");
        leg_rel->Draw("SAME");
        
        //tex->DrawLatex(0.53,0.91,CorrectionObject::_lumitag+"(13TeV)");
        
        c_rel_r[j][k]->SaveAs(CorrectionObject::_outpath+"plots/control/Rel_ResponseVsAlpha_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]+".pdf");
 
        name_mpf_r[j][k]="mpf_r_alpha_"+eta_range[j]+"_"+eta_range[j+1]+"_pT_"+pt_range[k]+"_"+pt_range[k+1];
        c_mpf_r[j][k] = new TCanvas(name_mpf_r[j][k], name_mpf_r[j][k], 850,700);
 
        m_gStyle->SetOptTitle(0);
  
        graph_mpf_r_mc[k][j]->Draw("AP");
        graph_mpf_r_data[k][j]->Draw("P SAME");
 
        graph_mpf_r_mc[k][j]->GetYaxis()->SetRangeUser(0.9,1.2);
        graph_mpf_r_mc[k][j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
        graph_mpf_r_mc[k][j]->GetYaxis()->SetTitle("R_{MPF}");
        graph_mpf_r_mc[k][j]->GetYaxis()->SetTitleSize(0.045);
        graph_mpf_r_mc[k][j]->GetYaxis()->SetTitleOffset(1.);
        graph_mpf_r_mc[k][j]->GetXaxis()->SetTitle("cut on #alpha");
        graph_mpf_r_mc[k][j]->GetXaxis()->SetTitleSize(0.045);
 
        graph_mpf_r_data[k][j]->GetYaxis()->SetRangeUser(0.9,1.2);
        graph_mpf_r_data[k][j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
        graph_mpf_r_data[k][j]->GetYaxis()->SetTitle("R_{MPF}");
        graph_mpf_r_data[k][j]->GetYaxis()->SetTitleSize(0.045);
        graph_mpf_r_data[k][j]->GetYaxis()->SetTitleOffset(1.);
        graph_mpf_r_data[k][j]->GetXaxis()->SetTitle("cut on #alpha");
        graph_mpf_r_data[k][j]->GetXaxis()->SetTitleSize(0.045);
 
        
        TLegend* leg_mpf = new TLegend(0.20,0.15,0.31,0.35,"","brNDC");//x+0.1
        leg_mpf->SetBorderSize(0);
        leg_mpf->SetTextSize(0.038);
        leg_mpf->SetFillColor(10);
        leg_mpf->SetFillStyle(0);
        leg_mpf->SetLineColor(1);
        leg_mpf->SetTextFont(42);
        leg_mpf->SetHeader("MPF response, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]+", "+(eta_cut_bool?pt_range_HF:pt_range)[k]+"#leq p_{T}<"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]);
        leg_mpf->AddEntry(graph_mpf_r_mc[k][j], "R^{MC}","P");
        leg_mpf->AddEntry(graph_mpf_r_data[k][j], "R^{DATA}","P");
        leg_mpf->Draw("SAME");
        
        //tex->DrawLatex(0.53,0.91,CorrectionObject::_lumitag+"(13TeV)");
        
        c_mpf_r[j][k]->SaveAs(CorrectionObject::_outpath+"plots/control/MPF_ResponseVsAlpha_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]+".pdf");
      } 
    }
   
   bool multigraph_rel_empty[n_eta-1];
   bool multigraph_mpf_empty[n_eta-1];

   bool firstHFeta = true;
   
   for(int j=0; j<n_eta-1; j++){
     multigraph_rel_empty[j] = true;
     multigraph_mpf_empty[j] = true;
   }
   for(int j=0; j<n_eta-1; j++){
     eta_cut_bool = fabs(eta_bins[j])>eta_cut;
     pTgraph_rel_r[j] = new TMultiGraph();
     pTgraph_mpf_r[j] = new TMultiGraph();
     eta_cut_bool = fabs(eta_bins[j])>eta_cut;     
     for(int k = 1 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){ 
       graph_rel_r[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_rel_r[k][j],zero,err_ratio_al_rel_r[k][j]);
       graph_rel_r[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_rel_r[k][j]);

       graph_rel_r[k][j]->SetMarkerSize(1.5);
       graph_rel_r[k][j]->SetMarkerStyle(20);
       if(k<9){
	 graph_rel_r[k][j]->SetMarkerColor(k+1);
	 graph_rel_r[k][j]->SetLineColor(k+1);
       }
       else{
	 graph_rel_r[k][j]->SetMarkerColor(k+19);
	 graph_rel_r[k][j]->SetLineColor(k+19);
       }
       if(graph_rel_r[k][j]->GetN()>0){
	 pTgraph_rel_r[j]->Add(graph_rel_r[k][j]); //one multigraph consisting of several TGraphErrors. One Multigraph for each eta bin
	 multigraph_rel_empty[j] = false;
       }

       TString pTbin_label = "";
       pTbin_label+=(eta_cut_bool?pt_bins_HF:pt_bins)[k];
       pTbin_label+=" < p_T < ";
       pTbin_label+=(eta_cut_bool?pt_bins_HF:pt_bins)[k+1];
       if(j==0) leg1->AddEntry(graph_rel_r[k][j],pTbin_label,"epl");     
       if(eta_cut_bool && firstHFeta) leg1_HF->AddEntry(graph_rel_r[k][j],pTbin_label,"epl");       
      
       graph_mpf_r[k][j] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_mpf_r[k][j],zero,err_ratio_al_mpf_r[k][j]);
       graph_mpf_r[k][j] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_r[k][j]);
       graph_mpf_r[k][j]->SetMarkerSize(1.3);
       graph_mpf_r[k][j]->SetMarkerStyle(20);
       if(k<9){
	 graph_mpf_r[k][j]->SetMarkerColor(k+1);
	 graph_mpf_r[k][j]->SetLineColor(k+1);
       }
       else{
	 graph_mpf_r[k][j]->SetMarkerColor(k+19);
	 graph_mpf_r[k][j]->SetLineColor(k+19);
       }
       if(graph_mpf_r[k][j]->GetN()>0) {
	 pTgraph_mpf_r[j]->Add(graph_mpf_r[k][j]);
	 multigraph_mpf_empty[j] = false;
       }
     }
     firstHFeta = firstHFeta && !eta_cut_bool;
   }

  //************************************* Test 2D plot kFSR  *************************************


   //Create horizontal line for plotting ("ideal value")
   TLine *line = new TLine(alpha_bins[0],1,alpha_bins[n_alpha-1]+0.01,1);

  TH2Poly* h_kFSR_pt_eta_rel = new TH2Poly();
  h_kFSR_pt_eta_rel->SetName("kFSR_pt_eta_rel");  
  h_kFSR_pt_eta_rel->SetTitle("kFSR");
  h_kFSR_pt_eta_rel->GetXaxis()->SetTitle("|#eta|");
  h_kFSR_pt_eta_rel->GetYaxis()->SetTitle("p_{T}^{ave}"); 
  h_kFSR_pt_eta_rel->SetStats(kFALSE);
  int bincount=0;
for(int i=0; i<n_eta-1; i++){   
   eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
     int realbinnumber = h_kFSR_pt_eta_rel->AddBin(eta_bins[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j], eta_bins[i+1], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1] );
     bincount++;
   }
  }

  TH2Poly* h_chi2_kFSR_rel = new TH2Poly();
  h_chi2_kFSR_rel->SetName("chi2_kFSR_rel");  
  h_chi2_kFSR_rel->SetTitle("#chi^{2} kFSR");
  h_chi2_kFSR_rel->GetXaxis()->SetTitle("|#eta|");
  h_chi2_kFSR_rel->GetYaxis()->SetTitle("p_{T}^{ave}"); 
  h_chi2_kFSR_rel->SetStats(kFALSE);
  for(int i=0; i<n_eta-1; i++){   
   eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
     h_chi2_kFSR_rel->AddBin(eta_bins[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j], eta_bins[i+1], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1] );
   }
  }
  
   TCanvas* Rel[n_eta-1][n_pt_-1];
   TString plotname_rel[n_eta-1][n_pt_-1];

   int bincounter = 1;
   TF1 *pol_rel[n_eta-1][n_pt_-1];
    for(int i=0; i<n_eta-1; i++){
     eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
     for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){ 
	plotname_rel[i][j]="dijet_kfsr_eta_"+eta_range[i]+"_"+eta_range[i+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
	Rel[i][j] = new TCanvas(plotname_rel[i][j], plotname_rel[i][j], 800,700);
	m_gStyle->SetOptTitle(0);

	pol_rel[i][j] = new TF1("pol_rel","pol1",jet3pt_min/(eta_cut_bool?pt_bins_HF:pt_bins)[j],0.4);//,0.14,0.36);
       if(j==0) continue;
       graph_rel_r[j][i] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_rel_r[j][i],zero,err_ratio_al_rel_r[j][i]);
       graph_rel_r[j][i] = (TGraphErrors*)CleanEmptyPoints(graph_rel_r[j][i]);


       if(graph_rel_r[j][i]->GetN()>0) {
       pol_rel[i][j]->SetParameters(1.1,0.1);
       graph_rel_r[j][i]->Fit(pol_rel[i][j],"RM");
       }
       else {
       pol_rel[i][j]->SetParameters(-1,-1);
       pol_rel[i][j]->SetParError(0,1);
       pol_rel[i][j]->SetParError(1,1);
       }
     TLatex *tex_rel = new TLatex();
     tex_rel->SetNDC();
     tex_rel->SetTextSize(0.045); 
     tex_rel->DrawLatex(0.64,0.91,CorrectionObject::_lumitag);

     TString chi2_rel = "#chi^{2}/n.d.f = ";
     TLatex *tex2_rel = new TLatex();
     double chi2ndf_kFSR_rel = -1;
     TString bin_text="p_{T} balance, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+", "+(eta_cut_bool?pt_range_HF:pt_range)[j]+"#leq p_{T}<"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
      if(graph_rel_r[j][i]->GetN()>0){
       chi2_rel += trunc(pol_rel[i][j]->GetChisquare());
       chi2_rel +="/";
       chi2_rel +=trunc(pol_rel[i][j]->GetNDF());

       tex2_rel->SetNDC();
       tex2_rel->SetTextSize(0.035); 
       tex2_rel->DrawLatex(0.24,0.25,chi2_rel);
       tex2_rel->DrawLatex(0.24,0.35,bin_text);
       h_kFSR_pt_eta_rel->Fill(eta_bins[i]+0.01, (eta_cut_bool?pt_bins_HF:pt_bins)[j]+0.01, pol_rel[i][j]->GetParameter(0));
       chi2ndf_kFSR_rel = pol_rel[i][j]->GetChisquare() / pol_rel[i][j]->GetNDF();
       h_chi2_kFSR_rel->Fill(eta_bins[i]+0.01, (eta_cut_bool?pt_bins_HF:pt_bins)[j]+0.01, chi2ndf_kFSR_rel); 
      }

      //Cosmetics 
      graph_rel_r[j][i]->SetMarkerSize(1.3);
      graph_rel_r[j][i]->SetMarkerStyle(20);
      graph_rel_r[j][i]->GetYaxis()->SetRangeUser(0.92,1.08);
      graph_rel_r[j][i]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.1);
      graph_rel_r[j][i]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})_{#alpha<X}/(R_{MC}/R_{DATA})_{#alpha<0.3}");
      graph_rel_r[j][i]->GetXaxis()->SetTitle("cut on #alpha");
      graph_rel_r[j][i]->Draw("AP");//to change range draw graph first, change range, draw it again and than update canvas
      graph_rel_r[j][i]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.1);
      graph_rel_r[j][i]->Draw("AP");
      Rel[i][j]->Update();
      line->SetLineStyle(2);
      line ->Draw("SAME");
      tex2_rel->DrawLatex(0.24,0.85,bin_text);
      Rel[i][j]->Print(CorrectionObject::_outpath+"plots/control/kFSR_Pt_eta_"+eta_range2[i]+"_"+eta_range2[i+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1]+".pdf");
     delete tex2_rel;
     delete tex_rel;
     bincounter++;
     }
    }


 TCanvas* Pt_dep_rel[n_eta-1];
 TString name_rel[n_eta-1];
 TH1D* h_kFSR_pt_rel[n_eta-1];
 for(int i=0; i<n_eta-1; i++){
   eta_cut_bool = fabs(eta_bins[i])>eta_cut;
   name_rel[i]="dijet_kfsr_pt_dep_"+eta_range[i]+"_"+eta_range[i+1];
   h_kFSR_pt_rel[i] = new TH1D(name_rel[i],"kFSR pt dependence", n_pt_-1, pt_bins);
   Pt_dep_rel[i] = new TCanvas(name_rel[i], name_rel[i], 800,700);
   
   eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){ 
     h_kFSR_pt_rel[i]->SetBinContent(j+1, pol_rel[i][j]->GetParameter(0));
     h_kFSR_pt_rel[i]->SetBinError(j+1, pol_rel[i][j]->GetParError(0));
   }
   h_kFSR_pt_rel[i]->GetYaxis()->SetRangeUser(0.92,1.08);
   h_kFSR_pt_rel[i]->SetMarkerStyle(20);
   h_kFSR_pt_rel[i]->GetYaxis()->SetTitle("kFSR");
   h_kFSR_pt_rel[i]->GetXaxis()->SetTitle("p_{T}");
   h_kFSR_pt_rel[i]->Draw("E");
   leg2->SetHeader("p_{T} balance, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]);
   leg2->Draw();
   Pt_dep_rel[i]->SetLogx();
   //   Float_t value_multifit_rel = pol1[i]->GetParameter(0);
   // TLine *line_fit = new TLine(pt_bins[0],value_multifit_rel,pt_bins[n_pt_-1],value_multifit_rel);
   // line_fit->SetLineStyle(2);
   // line_fit->Draw("SAME");
   Pt_dep_rel[i]->Print(CorrectionObject::_outpath+"plots/control/kFSR_Pt_eta_"+eta_range2[i]+"_"+eta_range2[i+1]+".pdf");
 }

 
  TCanvas* c1 = new TCanvas();
  c1->SetLogy();
  h_kFSR_pt_eta_rel->GetZaxis()->SetRangeUser(0.9,1.1);
  h_kFSR_pt_eta_rel->GetZaxis()->SetTitle("kFSR");
  h_kFSR_pt_eta_rel->Draw("COLZ");
  c1->Print(CorrectionObject::_outpath+"plots/kFSR_Pt_"+CorrectionObject::_generator_tag+"_2DPlot.pdf");

  TCanvas* c3 = new TCanvas();
  c3->SetLogy();
  h_chi2_kFSR_rel->GetZaxis()->SetRangeUser(0,20);
  h_chi2_kFSR_rel->GetZaxis()->SetTitle("#chi^{2}");  
  h_chi2_kFSR_rel->Draw("COLZ");

  //
  c3->Print(CorrectionObject::_outpath+"plots/Chi2_kFSR_pT_"+CorrectionObject::_generator_tag+"_2DPlot.pdf");
  //*****************************************************************************************************


  //  return;//stop here for tests

  TH2Poly* h_kFSR_pt_eta_mpf = new TH2Poly();
  h_kFSR_pt_eta_mpf->SetName("kFSR_pt_eta_mpf");  
  h_kFSR_pt_eta_mpf->SetTitle("kFSR");
  h_kFSR_pt_eta_mpf->GetXaxis()->SetTitle("|#eta|");
  h_kFSR_pt_eta_mpf->GetYaxis()->SetTitle("p_{T}^{ave}"); 
  h_kFSR_pt_eta_mpf->SetStats(kFALSE);
  for(int i=0; i<n_eta-1; i++){   
   eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
     h_kFSR_pt_eta_mpf->AddBin(eta_bins[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j], eta_bins[i+1], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1] );
   }
  }

  TH2Poly* h_chi2_kFSR_mpf = new TH2Poly();
  h_chi2_kFSR_mpf->SetName("chi2_kFSR_rel");  
  h_chi2_kFSR_mpf->SetTitle("#chi^{2} kFSR");
  h_chi2_kFSR_mpf->GetXaxis()->SetTitle("|#eta|");
  h_chi2_kFSR_mpf->GetYaxis()->SetTitle("p_{T}^{ave}"); 
  h_chi2_kFSR_mpf->SetStats(kFALSE);
  for(int i=0; i<n_eta-1; i++){   
   eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
     h_chi2_kFSR_mpf->AddBin(eta_bins[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j], eta_bins[i+1], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1] );
   }
  }


  TString plotname2[n_eta-1][n_pt_-1];
  TCanvas* MPF[n_eta-1][n_pt_-1];
   TF1 *pol_mpf[n_eta-1][n_pt_-1];
   bincounter = 1;
    for(int i=0; i<n_eta-1; i++){
   eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
	plotname2[i][j]="MPF_kfsr_eta_"+eta_range[i]+"_"+eta_range[i+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
	MPF[i][j] = new TCanvas(plotname2[i][j], plotname2[i][j], 800,700);
	m_gStyle->SetOptTitle(0);

	pol_mpf[i][j] = new TF1("pol_mpf","pol1",jet3pt_min/(eta_cut_bool?pt_bins_HF:pt_bins)[j],0.4);//,0.14,0.36);
       if(j==0) continue;
       graph_mpf_r[j][i] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_mpf_r[j][i],zero,err_ratio_al_mpf_r[j][i]);
       graph_mpf_r[j][i] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_r[j][i]);

      //Cosmetics 
       graph_mpf_r[j][i]->SetMarkerSize(1.3);
       graph_mpf_r[j][i]->SetMarkerStyle(21);
       graph_mpf_r[j][i]->GetYaxis()->SetRangeUser(0.92,1.08);
       graph_mpf_r[j][i]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.1);
       graph_mpf_r[j][i]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})_{#alpha<X}/(R_{MC}/R_{DATA})_{#alpha<0.3}");
       graph_mpf_r[j][i]->GetXaxis()->SetTitle("cut on #alpha");

       graph_mpf_r[j][i]->Draw("AP");
       graph_mpf_r[j][i]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.1);
       graph_mpf_r[j][i]->Draw("AP");
       MPF[i][j]->Update();
       line->SetLineStyle(2);
       line ->Draw("SAME");

       if(graph_mpf_r[j][i]->GetN()>0) {
       pol_mpf[i][j]->SetParameters(1.0,0.0);
       graph_mpf_r[j][i]->Fit(pol_mpf[i][j],"RM");
       }
       else {
       pol_mpf[i][j]->SetParameters(-1,-1);
       pol_mpf[i][j]->SetParError(0,1);
       pol_mpf[i][j]->SetParError(1,1);
     }

    TLatex *tex_mpf = new TLatex();
     tex_mpf->SetNDC();
     tex_mpf->SetTextSize(0.045); 
     tex_mpf->DrawLatex(0.64,0.91,CorrectionObject::_lumitag);

     TString chi2_mpf = "#chi^{2}/n.d.f = ";
     TLatex *tex2_mpf = new TLatex();
     double chi2ndf_kFSR_mpf = -1;
      if(graph_rel_r[j][i]->GetN()>0){
       chi2_mpf += trunc(pol_mpf[i][j]->GetChisquare());
       chi2_mpf +="/";
       chi2_mpf +=trunc(pol_mpf[i][j]->GetNDF());

       tex2_mpf->SetNDC();
       tex2_mpf->SetTextSize(0.035); 
       tex2_mpf->DrawLatex(0.24,0.25,chi2_mpf);
       TString bin_text="MPF, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+", "+(eta_cut_bool?pt_range_HF:pt_range)[j]+"#leq p_{T}<"+(eta_cut_bool?pt_range_HF:pt_range)[j+1];
       tex2_mpf->DrawLatex(0.24,0.35,bin_text);
       h_kFSR_pt_eta_mpf->Fill(eta_bins[i]+0.01, (eta_cut_bool?pt_bins_HF:pt_bins)[j]+0.01, pol_rel[i][j]->GetParameter(0)); 
       chi2ndf_kFSR_mpf = pol_mpf[i][j]->GetChisquare() / pol_mpf[i][j]->GetNDF();
       h_chi2_kFSR_mpf->Fill(eta_bins[i]+0.01, (eta_cut_bool?pt_bins_HF:pt_bins)[j]+0.01,chi2ndf_kFSR_mpf); 
      }

      MPF[i][j]->Print(CorrectionObject::_outpath+"plots/control/kFSR_MPF_eta_"+eta_range2[i]+"_"+eta_range2[i+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[j]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[j+1]+".pdf");

     delete tex2_mpf;
     delete tex_mpf;
     bincounter++;
     }
   }


 TCanvas* Pt_dep_mpf[n_eta-1];
 TString name_mpf[n_eta-1];
 TH1D* h_kFSR_pt_mpf[n_eta-1]; 
 for(int i=0; i<n_eta-1; i++){
   name_mpf[i]="MPF_kfsr_pt_dep_"+eta_range[i]+"_"+eta_range[i+1];
   eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   h_kFSR_pt_mpf[i]  = new TH1D(name_mpf[i],"kFSR pt dependence", ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) , (eta_cut_bool ? pt_bins_HF : pt_bins));
   Pt_dep_mpf[i] = new TCanvas(name_mpf[i], name_mpf[i], 800,700);
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
     h_kFSR_pt_mpf[i]->SetBinContent(j+1, pol_mpf[i][j]->GetParameter(0));
     h_kFSR_pt_mpf[i]->SetBinError(j+1, pol_rel[i][j]->GetParError(0));
   }
   h_kFSR_pt_mpf[i]->SetMarkerSize(1.5);
   h_kFSR_pt_mpf[i]->GetYaxis()->SetRangeUser(0.92,1.08);
   h_kFSR_pt_mpf[i]->GetYaxis()->SetTitle("kFSR");
   h_kFSR_pt_mpf[i]->GetXaxis()->SetTitle("p_{T}");
   h_kFSR_pt_mpf[i]->SetMarkerStyle(21);
   h_kFSR_pt_mpf[i]->Draw("");
   leg2->SetHeader("MPF, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]);
   leg2->Draw();

   Pt_dep_mpf[i]->SetLogx();
   Pt_dep_mpf[i]->Print(CorrectionObject::_outpath+"plots/control/kFSR_MPF_eta_"+eta_range2[i]+"_"+eta_range2[i+1]+".pdf");
 }

 
  TCanvas* c2 = new TCanvas();
  c2->SetLogy();
  h_kFSR_pt_eta_mpf->GetZaxis()->SetRangeUser(0.9,1.1);
  h_kFSR_pt_eta_mpf->GetZaxis()->SetTitle("kFSR");  
  h_kFSR_pt_eta_mpf->Draw("COLZ");
  
  c2->Print(CorrectionObject::_outpath+"plots/kFSR_MPF_"+CorrectionObject::_generator_tag+"_2DPlot.pdf");
  
  TCanvas* c4 = new TCanvas();
  c4->SetLogy();
  h_chi2_kFSR_mpf->GetZaxis()->SetRangeUser(0,20);  
  h_chi2_kFSR_mpf->GetZaxis()->SetTitle("#chi^{2}");  
  h_chi2_kFSR_mpf->Draw("COLZ");
  m_gStyle->SetOptStat(kFALSE);  

  c4->Print(CorrectionObject::_outpath+"plots/Chi2_kFSR_MPF_"+CorrectionObject::_generator_tag+"_2DPlot.pdf");
 //**********************************************************************************************

   /* +++++++++++++++++ PLOTS ++++++++++++++++++ */
   //Do the well-known, colorful kFSR linear fit plots
  
   //First set up the output files   
   //Create output .dat file, including the kFSR extrapolation (alpha->0)
   FILE *fp_rel_r, *fp_mpf_r; 
   TH1D *kFSR_MPF, *kFSR_DiJet, *plotkfsr;

   CorrectionObject::make_path( _outpath+"output/");
   cout << "Opening .dat files at: " << CorrectionObject::_outpath+"output/KFSR_MPF_extrapolation.dat" << endl;
   fp_mpf_r = fopen(CorrectionObject::_outpath+"output/KFSR_MPF_extrapolation.dat","w");
   fp_rel_r = fopen(CorrectionObject::_outpath+"output/KFSR_DiJet_extrapolation.dat","w");

   kFSR_MPF = new TH1D("kfsr_mpf","kfsr_mpf", n_eta-1,eta_bins);
   kFSR_DiJet = new TH1D("kfsr_dijet","kfsr_dijet", n_eta-1,eta_bins);
   plotkfsr = new TH1D("kfsr","kfsr", n_eta-1,eta_bins);

  // Start with pT-balance
  //*****************************************************************************************************************************************************
   //create plots
   TCanvas* a[n_eta-1];
   TString plotname[n_eta-1];
   TF1 *pol1[n_eta-1];
   for (int j=0; j<n_eta-1; j++){ 
     eta_cut_bool = fabs(eta_bins[j])>eta_cut;
     
     plotname[j]="dijet_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
     a[j] = new TCanvas(plotname[j], plotname[j], 800,700);
     m_gStyle->SetOptTitle(0);
     pTgraph_rel_r[j]->Draw("AP");

     if(!multigraph_rel_empty[j]){
       pTgraph_rel_r[j]->GetYaxis()->SetRangeUser(0.92,1.08);
       pTgraph_rel_r[j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
       pTgraph_rel_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})");
       pTgraph_rel_r[j]->GetXaxis()->SetTitle("cut on #alpha");

     }


     pol1[j] = new TF1("pol1","pol1",0.1,0.4);//,0.14,0.36);
     pol1[j]->SetParameters(1.0,0.1);
     if(j == 13){
       pol1[j]->SetParameters(0.985,0.05);  
     }
 
     if (multigraph_rel_empty[j]) cout << "Eta bin no. " << j << ", multigraph empty!" << endl;
     else cout << "Eta bin no. " << j << ", multigraph filled!" << endl;
     if(!multigraph_rel_empty[j]){
       pol1[j]->SetParameters(1.0,0.1);
       pTgraph_rel_r[j]->Fit(pol1[j],"RM");
     }
     else {
       pol1[j]->SetParameters(-1,-1);
       pol1[j]->SetParError(0,1);
       pol1[j]->SetParError(1,1);
     }
     line->SetLineStyle(2);
     line->Draw("SAME");

     // fill the output.dat file
     if (fp_rel_r!=NULL) {
       Float_t value = pol1[j]->GetParameter(0);
       Float_t uncert = pol1[j]->GetParError(0);
       fprintf(fp_rel_r, "%f %f\n",value,uncert);
     }
     if(!multigraph_rel_empty[j]){
       plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0));
       plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
       kFSR_DiJet->SetBinContent(j+1,pol1[j]->GetParameter(0));
       kFSR_DiJet->SetBinError(j+1,pol1[j]->GetParError(0));
     }

     if(eta_cut_bool){
       leg1_HF->SetHeader("p_{T} balance, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
       leg1_HF->Draw();
     }
     else{
       leg1->SetHeader("p_{T} balance, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
       leg1->Draw();
     }
     
     TLatex *tex = new TLatex();
     tex->SetNDC();
     tex->SetTextSize(0.045); 
     tex->DrawLatex(0.64,0.91,CorrectionObject::_lumitag);

     TString chi2_loglin = "#chi^{2}/n.d.f = ";
     TLatex *tex2 = new TLatex();
     if(!multigraph_rel_empty[j]){

       chi2_loglin += trunc(pol1[j]->GetChisquare());
       chi2_loglin +="/";
       chi2_loglin +=trunc(pol1[j]->GetNDF());


       tex2->SetNDC();
       tex2->SetTextSize(0.035); 
       tex2->DrawLatex(0.64,0.35,chi2_loglin);
     }
     cout << "Printing kFSR plots to " << CorrectionObject::_outpath+"plots/control/kFSR_Pt_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf" << endl;
     a[j]->Print(CorrectionObject::_outpath+"plots/control/kFSR_Pt_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");


     //delete stuff

     delete tex2;
     delete tex;
   }
   
   
   cout << endl << endl << "finished all fits for rel" << endl  << endl;
   fclose(fp_rel_r);
   cout << "closed some file, now opening output file" << endl;
   
   // create output file including the kFSR plot
   TFile* outputfile_rel_r;
   cout << "Creating output-rootfile:" << CorrectionObject::_outpath+"Histo_KFSR_DiJet_"+CorrectionObject::_generator_tag+"_L1.root" << endl;
   outputfile_rel_r = new TFile(CorrectionObject::_outpath+"Histo_KFSR_DiJet_"+CorrectionObject::_generator_tag+"_L1.root","RECREATE");
   
   cout << "now writing kFSR values" << endl;
   kFSR_DiJet->Write();
   h_kFSR_pt_eta_rel->Write();
   outputfile_rel_r->Write();
   cout << "closing output file" << endl;
   outputfile_rel_r->Close();
   cout << "closed output file" << endl;



   // And now MPF results
   //create plots
   TCanvas* b[n_eta-1];
   TString plotname1[n_eta-1];
    for (int j=0; j<n_eta-1; j++){
      eta_cut_bool = fabs(eta_bins[j])>eta_cut;
      
     plotname1[j]="mpf_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
     b[j] = new TCanvas(plotname1[j], plotname1[j], 800,700);
     m_gStyle->SetOptTitle(0);

     pTgraph_mpf_r[j]->Draw("AP");
     if(!multigraph_mpf_empty[j]){
       pTgraph_mpf_r[j]->GetYaxis()->SetRangeUser(0.92,1.08);
       pTgraph_mpf_r[j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
       pTgraph_mpf_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})");
       pTgraph_mpf_r[j]->GetXaxis()->SetTitle("cut on #alpha");
     }

     pol1[j] = new TF1("pol1","pol1",0.1,0.4);//,0.14,0.36);  //TEST AK4
     pol1[j]->SetParameters(0,0);

     if(!multigraph_mpf_empty[j]){
       pTgraph_mpf_r[j]->Fit(pol1[j],"R");
     }
     else{
       pol1[j]->SetParameters(1.03,-0.1);
       pol1[j]->SetParError(0,1);
       pol1[j]->SetParError(1,1);
     }
     line->SetLineStyle(2);
     line->Draw("SAME");


     // fill the output.dat file
     if (fp_mpf_r!=NULL) {
       Float_t value = pol1[j]->GetParameter(0);
       Float_t uncert = pol1[j]->GetParError(0);
       fprintf(fp_mpf_r, "%f %f\n",value,uncert);
     }
     if(!multigraph_mpf_empty[j]) {
       plotkfsr->SetBinContent(j+1,pol1[j]->GetParameter(0));
       plotkfsr->SetBinError(j+1,pol1[j]->GetParError(0));
       kFSR_MPF->SetBinContent(j+1,pol1[j]->GetParameter(0));
       kFSR_MPF->SetBinError(j+1,pol1[j]->GetParError(0));
     }

     if(eta_cut_bool){
       leg1_HF->SetHeader("MPF, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
       leg1_HF->Draw();
     }
     else{
       leg1->SetHeader("MPF, "+eta_range[j]+"#leq|#eta|<"+eta_range[j+1]);
       leg1->Draw();
     }
     
     TLatex *tex = new TLatex();
     tex->SetNDC();
     tex->SetTextSize(0.045); 
     tex->DrawLatex(0.64,0.91,CorrectionObject::_lumitag);

     TString chi2_loglin = "#chi^{2}/n.d.f = ";
     TLatex *tex2 = new TLatex();
     if(!multigraph_mpf_empty[j]){

       chi2_loglin += trunc(pol1[j]->GetChisquare());
       chi2_loglin +="/";
       chi2_loglin +=trunc(pol1[j]->GetNDF());

 
       tex2->SetNDC();
       tex2->SetTextSize(0.035); 
       tex2->DrawLatex(0.64,0.35,chi2_loglin);
     }

     //save the plots
     cout << "Saving the MPF plots to " << CorrectionObject::_outpath + "plots/control/kFSR_MPF_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf" << endl;
     b[j]->Print(CorrectionObject::_outpath + "plots/control/kFSR_MPF_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");

     //delete stuff
     delete tex2;
     delete tex;
   }
   fclose(fp_mpf_r);

   // create output file including the kFSR plot
   TFile* outputfile_mpf_r;
   cout << "Creating output-rootfile: " << CorrectionObject::_outpath+"Histo_KFSR_MPF_"+CorrectionObject::_generator_tag+"_L1.root" << endl;
   outputfile_mpf_r = new TFile(CorrectionObject::_outpath+"Histo_KFSR_MPF_"+CorrectionObject::_generator_tag+"_L1.root","RECREATE");
   kFSR_MPF->Write();
   h_kFSR_pt_eta_mpf->Write();
   outputfile_mpf_r->Write();
   outputfile_mpf_r->Close();


   cout << "+++++++++++++++++ Finished kFSR() +++++++++++++++++++" << endl;

   //delete everything
   
   for (int j=0; j<n_eta-1; j++){ //n_eta-1
     delete pol1[j];
   }

   delete outputfile_mpf_r;
   delete outputfile_rel_r;
   
   delete plotkfsr;
   delete kFSR_DiJet;
   delete kFSR_MPF;
   delete line;
   
   for(int j=0; j<n_eta-1; j++){
     eta_cut_bool = fabs(eta_bins[j])>eta_cut;
     for(int k= 1 ; k <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; k++ ){
       delete graph_rel_r[k][j];
       delete graph_mpf_r[k][j];
     }
   }
   
   for(int j=0; j<n_eta-1; j++){
     delete pTgraph_rel_r[j];
     delete pTgraph_mpf_r[j];
   }
   delete leg1;
   delete leg1_HF;   
   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
       delete pr_data_asymmetry[j][i];
       delete pr_data_B[j][i];
       delete pr_mc_asymmetry[j][i];
       delete pr_mc_B[j][i];
       delete hdata_asymmetry[j][i];
       delete hdata_B[j][i];
       delete hmc_asymmetry[j][i];
       delete hmc_B[j][i];
       }
     }
	       
   delete  m_gStyle;    
   cout << "++++++++++++ Deleted everything in kFSR(), exiting ++++++++++++++" << endl;
}
