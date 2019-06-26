#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
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



using namespace std;

void CorrectionObject::kFSR_CorrectFormulae_eta(){
  cout << "--------------- Starting kFSR() ---------------" << endl << endl;
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2(kTRUE);

// get ratio for MC to DATA responses
  double ratio_al_rel_r[n_pt-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al_rel_r[n_pt-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins
  double ratio_al_mpf_r[n_pt-1][n_eta-1][n_alpha]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf_r[n_pt-1][n_eta-1][n_alpha]; //error of ratio at pt,eta,alpha bins
  TProfile *pr_data_asymmetry[n_eta-1][n_alpha];// pT-balance response for data  
  TProfile *pr_data_B[n_eta-1][n_alpha];//MPF response for data
  TProfile *pr_mc_asymmetry[n_eta-1][n_alpha];// pT-balanse responce for MC  
  TProfile *pr_mc_B[n_eta-1][n_alpha];//MPF response for MC
  TH2D *hdata_asymmetry[n_eta-1][n_alpha];
  TH2D *hdata_B[n_eta-1][n_alpha];
  TH2D *hmc_asymmetry[n_eta-1][n_alpha];
  TH2D *hmc_B[n_eta-1][n_alpha];
  int n_entries_mc[n_eta-1][n_alpha][n_pt-1];
  int n_entries_data[n_eta-1][n_alpha][n_pt-1];
  int count = 0;
  TString name1 = "hist_data_asymmetry_";
  TString name2 = "hist_data_B_";
  TString name3 = "hist_mc_asymmetry_";
  TString name4 = "hist_mc_B_";
  for(int i=0; i<n_alpha; i++){
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){
	ratio_al_rel_r[k][j][i] = 0;
	err_ratio_al_rel_r[k][j][i] = 0;
	ratio_al_mpf_r[k][j][i] = 0;
	err_ratio_al_mpf_r[k][j][i] = 0;
      }
      
      TString name = name1; name+=count;
      hdata_asymmetry[j][i] = new TH2D(name,"A in DATA; p_{T}^{ave} [GeV]; A",n_pt-1,pt_bins,nResponseBins, Response_min, Response_max);
      name = name2;name+=count;
      hdata_B[j][i] = new TH2D(name,"B in DATA;p_{T}^{ave} [GeV];B",n_pt-1,pt_bins,nResponseBins, Response_min, Response_max);
      name = name3; name+=count;
      hmc_asymmetry[j][i] = new TH2D(name,"A in MC;p_{T}^{ave} [GeV];A",n_pt-1,pt_bins,nResponseBins, Response_min, Response_max);
      name = name4; name+=count;
      hmc_B[j][i] = new TH2D(name,"B in MC;p_{T}^{ave} [GeV];B",n_pt-1,pt_bins,nResponseBins, Response_min, Response_max);
  
      count++;
    }
  }
  
  cout << "Set up a total of " << count << " histograms." << endl;

    for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
       for(int k=0; k<n_pt-1; k++){
	 n_entries_mc[j][i][k] = 0;
	 n_entries_data[j][i][k] = 0;
       }
     }
   }
  
  // Get relevant quantities from DATA, loop over events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
  TTreeReaderValue<Float_t> asymmetry_data(myReader_DATA, "asymmetry");
  TTreeReaderValue<Float_t> B_data(myReader_DATA, "B");   
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  int idx = 0;
  
  cout << "starting to loop over DATA events." << endl;

  while (myReader_DATA.Next()) {
    if(*alpha_data<0.05) continue; 
    for(int j=0; j<n_eta-1; j++){
      if(*probejet_eta_data>eta_bins[j+1] || *probejet_eta_data<eta_bins[j]) continue;
      for(int i=0; i<n_alpha; i++){
	if(*alpha_data>alpha_bins[i]) continue;
	else{
 	  hdata_asymmetry[j][i]->Fill(*pt_ave_data,*asymmetry_data,*weight_data);
	  hdata_B[j][i]->Fill(*pt_ave_data,*B_data,*weight_data);
	  for(int k=0; k<n_pt-1; k++){                                                 ///int k=0; k<n_pt-1; k++
	      if((*pt_ave_data < pt_bins[k]) || (*pt_ave_data >= pt_bins[k+1])) continue;
	       n_entries_data[j][i][k]++;
	    }
	  idx++;
	  if(idx%1000000==0) cout << "looping over data-TTree: Idx = " << idx << endl;
	}
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
   idx = 0;

   while (myReader_MC.Next()) {
     if(*alpha_mc<0.05) continue; 
     for(int j=0; j<n_eta-1; j++){
       if(*probejet_eta_mc>eta_bins[j+1] || *probejet_eta_mc<eta_bins[j]) continue;
       for(int i=0; i<n_alpha; i++){
	 if(*alpha_mc>alpha_bins[i]) continue;
	 else{
	   hmc_asymmetry[j][i]->Fill(*pt_ave_mc,*asymmetry_mc,*weight_mc);
	   hmc_B[j][i]->Fill(*pt_ave_mc,*B_mc,*weight_mc);
	   for(int k=0; k<n_pt-1; k++){                                                          ///int k=0; k<n_pt-1; k++
	     if((*pt_ave_mc < pt_bins[k]) || (*pt_ave_mc >= pt_bins[k+1])) continue;              //pt_bins[k]    pt_bins[k+1]
	       n_entries_mc[j][i][k]++;
	    }
	   idx++;
	   if(idx%1000000==0) cout << "looping over MC-TTree: Idx = " << idx << endl;
	 }
       }
     }
   }


   bool enough_entries[n_alpha][n_eta-1][n_pt-1];
   for(int i=0; i<n_alpha; i++){
     for(int j=0; j<n_eta-1; j++){
       for(int k=0; k<n_pt-1; k++){
	 enough_entries[i][j][k] = false;
	 if(n_entries_mc[j][i][k] > 100 && n_entries_data[j][i][k] > 100) enough_entries[i][j][k] = true;
	}
     }
   }
   //build profiles out of asymmetry and B 2d-histos to get <A> and <B> as a function of pT in bins of eta,alpha
   for(int j=0; j<n_eta-1; j++){
     for(int i=0; i<n_alpha; i++){
       //print TH2D histos
       TCanvas* c_dummy1 = new TCanvas();
       hdata_asymmetry[j][i]->Draw("COLZ");
       c_dummy1->SaveAs(CorrectionObject::_outpath+"plots/control/TH2_A_DATA_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_"+alpha_range[i]+".pdf");
       delete c_dummy1;
       TCanvas* c_dummy2 = new TCanvas();
       hmc_asymmetry[j][i]->Draw("COLZ");
       c_dummy2->SaveAs(CorrectionObject::_outpath+"plots/control/TH2_A_MC_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_"+alpha_range[i]+".pdf");
       delete c_dummy2;
       TCanvas* c_dummy3 = new TCanvas();
       hdata_B[j][i]->Draw("COLZ");
       c_dummy3->SaveAs(CorrectionObject::_outpath+"plots/control/TH2_B_DATA_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_"+alpha_range[i]+".pdf");
       delete c_dummy3;
       TCanvas* c_dummy4 = new TCanvas();
       hmc_B[j][i]->Draw("COLZ");
       c_dummy4->SaveAs(CorrectionObject::_outpath+"plots/control/TH2_B_MC_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_"+alpha_range[i]+".pdf");
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
       c_dummy5->SaveAs(CorrectionObject::_outpath+"plots/control/Profile_A_DATA_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_"+alpha_range[i]+".pdf");
       delete c_dummy5;
       TCanvas* c_dummy6 = new TCanvas();
       pr_mc_asymmetry[j][i]->Draw();
       c_dummy6->SetLogx();
       pr_mc_asymmetry[j][i]->GetYaxis()->SetTitle("<A>");
       pr_mc_asymmetry[j][i]->SetLineWidth(2);
       pr_mc_asymmetry[j][i]->SetMinimum(-0.3);
       pr_mc_asymmetry[j][i]->SetMaximum(0.3);
       c_dummy6->SaveAs(CorrectionObject::_outpath+"plots/control/Profile_A_MC_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_"+alpha_range[i]+".pdf");
       delete c_dummy6;
       TCanvas* c_dummy7 = new TCanvas();
       pr_data_B[j][i]->Draw();
       c_dummy7->SetLogx();
       pr_data_B[j][i]->GetYaxis()->SetTitle("<B>");
       pr_data_B[j][i]->SetLineWidth(2);
       pr_data_B[j][i]->SetMinimum(-0.3);
       pr_data_B[j][i]->SetMaximum(0.3);
       c_dummy7->SaveAs(CorrectionObject::_outpath+"plots/control/Profile_B_DATA_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_"+alpha_range[i]+".pdf");
       delete c_dummy7;
       TCanvas* c_dummy8 = new TCanvas();
       pr_mc_B[j][i]->Draw();
       c_dummy8->SetLogx();
       pr_mc_B[j][i]->GetYaxis()->SetTitle("<B>");
       pr_mc_B[j][i]->SetLineWidth(2);
       pr_mc_B[j][i]->SetMinimum(-0.3);
       pr_mc_B[j][i]->SetMaximum(0.3);
       c_dummy8->SaveAs(CorrectionObject::_outpath+"plots/control/Profile_B_MC_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+"_"+alpha_range[i]+".pdf");
       delete c_dummy8;
     }
   }
   
   //calculate response from <A> and <B> in bins of pt,eta,alpha
   //gaussian error propagation from errors on <A> and <B>

   for(int k=0; k<n_pt-1; k++){
     for(int j=0; j<n_eta-1; j++){
       for(int i=0; i<n_alpha; i++){
	 //responses for data, MC separately. Only for bins with >= 30 entries
	 double mpf_mc = (1+pr_mc_B[j][i]->GetBinContent(k+1))/(1-pr_mc_B[j][i]->GetBinContent(k+1));
	 //if(pr_mc_B[j][i]->GetBinEntries(k+1) < 30) mpf_mc = 0;
	 if(!enough_entries[i][j][k]) mpf_mc = 0;
	 double mpf_data = (1+pr_data_B[j][i]->GetBinContent(k+1))/(1-pr_data_B[j][i]->GetBinContent(k+1));
	 //if(pr_data_B[j][i]->GetBinEntries(k+1) < 30) mpf_data = 0;
	 if(!enough_entries[i][j][k]) mpf_data = 0;
	 double rel_mc = (1+pr_mc_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_mc_asymmetry[j][i]->GetBinContent(k+1));
	 //if(pr_mc_asymmetry[j][i]->GetBinEntries(k+1) < 30) rel_mc = 0;
	 if(!enough_entries[i][j][k]) rel_mc = 0;
	 double rel_data = (1+pr_data_asymmetry[j][i]->GetBinContent(k+1))/(1-pr_data_asymmetry[j][i]->GetBinContent(k+1));
	 //if(pr_data_asymmetry[j][i]->GetBinEntries(k+1) < 30) rel_data = 0;
	 if(!enough_entries[i][j][k]) rel_data = 0;
	 double err_mpf_mc = 2/(pow((1-pr_mc_B[j][i]->GetBinContent(k+1)),2)) * pr_mc_B[j][i]->GetBinError(k+1);
	 double err_mpf_data = 2/(pow((1-pr_data_B[j][i]->GetBinContent(k+1)),2)) * pr_data_B[j][i]->GetBinError(k+1);
	 double err_rel_mc = 2/(pow((1-pr_mc_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_mc_asymmetry[j][i]->GetBinError(k+1);
	 double err_rel_data = 2/(pow((1-pr_data_asymmetry[j][i]->GetBinContent(k+1)),2)) * pr_data_asymmetry[j][i]->GetBinError(k+1);

	 //ratio of responses, again gaussian error propagation
	 if(rel_data > 0) ratio_al_rel_r[k][j][i] = rel_mc/rel_data;
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
     for(int k=0; k<n_pt-1; k++){
       double norm_alref_rel_r = ratio_al_rel_r[k][j][al_ref];
       double err_norm_alref_rel_r = err_ratio_al_rel_r[k][j][al_ref];
       double norm_alref_mpf_r = ratio_al_mpf_r[k][j][al_ref];
       double err_norm_alref_mpf_r = err_ratio_al_mpf_r[k][j][al_ref];
       for(int i=0; i<n_alpha; i++){

	 if(norm_alref_rel_r>0){ //WHAT IS HAPPENING HERE? NO PROPER ERROR PROPAGATION !?! Ask other group that does kFSR-Extrapolation about error propagation
	   ratio_al_rel_r[k][j][i] =   ratio_al_rel_r[k][j][i]/norm_alref_rel_r; //original
	   err_ratio_al_rel_r[k][j][i] = sqrt(abs(pow(err_ratio_al_rel_r[k][j][i],2)-pow(err_norm_alref_rel_r,2)));//Original
	   //  err_ratio_al_rel_r[k][j][i] = sqrt(abs(pow(err_ratio_al_rel_r[k][j][i]/ratio_al_rel_r[k][j][i],2)-pow(err_norm_alref_rel_r/norm_alref_rel_r,2)));
	   //  err_ratio_al_rel_r[k][j][i] = sqrt(abs(pow(err_ratio_al_rel_r[k][j][i] / (ratio_al_rel_r[k][j][i]) ,2)+pow(err_norm_alref_rel_r / norm_alref_rel_r,2))) * ratio_al_rel_r[k][j][i] / norm_alref_rel_r ; //self
	  
	   if(i == al_ref) err_ratio_al_rel_r[k][j][i] = 0.;
	  

	 }
	 if(norm_alref_mpf_r>0){
	   ratio_al_mpf_r[k][j][i] =   ratio_al_mpf_r[k][j][i]/norm_alref_mpf_r;
	   err_ratio_al_mpf_r[k][j][i] = sqrt(abs(pow(err_ratio_al_mpf_r[k][j][i],2)-pow(err_norm_alref_mpf_r,2))); //Original
	   //  err_ratio_al_mpf_r[k][j][i] = sqrt(abs(pow(err_ratio_al_mpf_r[k][j][i]/ratio_al_mpf_r[k][j][i],2)-pow(err_norm_alref_mpf_r/norm_alref_mpf_r,2)));
	   //  err_ratio_al_mpf_r[k][j][i] = sqrt(abs(pow(err_ratio_al_mpf_r[k][j][i] / (ratio_al_mpf_r[k][j][i]) ,2)+pow(err_norm_alref_mpf_r / norm_alref_mpf_r,2))) * err_ratio_al_mpf_r[k][j][i] / norm_alref_mpf_r;
	  
	   if(i == al_ref) err_ratio_al_mpf_r[k][j][i] = 0.;
	 }
       }
     }
   }


   // Build the Multigraphs containing the ratio of responses (MC/DATA) as a function of alpha
   TGraphErrors *graph_rel_r[n_pt-1][n_eta-1];  //set of points vs alpha
   TMultiGraph *pTgraph_rel_r[n_eta-1];         //set of different pT bins in on eta bin
   TGraphErrors *graph_mpf_r[n_pt-1][n_eta-1];  //set of points vs alpha
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


   double xbin_tgraph[n_alpha],zero[n_alpha];
   for(int i=0;i<n_alpha;i++){
     xbin_tgraph[i] = alpha_bins[i];
     zero[i] = 0;
   }

   bool multigraph_rel_empty[n_eta-1];
   bool multigraph_mpf_empty[n_eta-1];
   for(int j=0; j<n_eta-1; j++){
     multigraph_rel_empty[j] = true;
     multigraph_mpf_empty[j] = true;
   }
   for(int j=0; j<n_eta-1; j++){
     pTgraph_rel_r[j] = new TMultiGraph();
     pTgraph_mpf_r[j] = new TMultiGraph();
     for(int k=0; k<n_pt-1; k++){
       if(pt_bins[k]<73) continue;
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
       pTbin_label+=pt_bins[k];
       pTbin_label+=" < p_{T} < ";
       pTbin_label+=pt_bins[k+1];
       if(j==0) leg1->AddEntry(graph_rel_r[k][j],pTbin_label,"epl");
      
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
   }

  //************************************* Test 2D plot kFSR  *************************************
   TH2D* h_kFSR_pt_eta_rel = new TH2D("kFSR_pt_eta_rel","kFSR;|#eta|;p_{T}^{ave}", n_eta-1, eta_bins, n_pt-1, pt_bins);
   TH2D* h_chi2_kFSR_rel = new TH2D("chi2_kFSR_rel","#chi^{2} kFSR;|#eta|;p_{T}^{ave}", n_eta-1, eta_bins, n_pt-1, pt_bins);
   TCanvas* Rel[n_eta-1][n_pt-1];
   TString plotname1[n_eta-1][n_pt-1];

   TF1 *pol_rel[n_eta-1][n_pt-1];
    for(int i=0; i<n_eta-1; i++){
      for(int j=0; j<n_pt-1; j++){
 
	plotname1[i][j]="dijet_kfsr_eta_"+eta_range[i]+"_"+eta_range[i+1]+"_pT_"+pt_range[j]+"_"+pt_range[j+1];
	Rel[i][j] = new TCanvas(plotname1[i][j], plotname1[i][j], 800,700);
	m_gStyle->SetOptTitle(0);
	

       pol_rel[i][j] = new TF1("pol_rel","pol1",0.14,0.36);
       if(pt_bins[j]<73) continue;
       graph_rel_r[j][i] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_rel_r[j][i],zero,err_ratio_al_rel_r[j][i]);
       graph_rel_r[j][i] = (TGraphErrors*)CleanEmptyPoints(graph_rel_r[j][i]);

       graph_rel_r[j][i]->Draw("AP");
       if(graph_rel_r[j][i]->GetN()>0) {

       pol_rel[i][j]->SetParameters(1.5,-0.5);
       graph_rel_r[j][i]->Fit(pol_rel[i][j],"RM");
       }
       else {
       pol_rel[i][j]->SetParameters(-1,-1);
       pol_rel[i][j]->SetParError(0,1);
       pol_rel[i][j]->SetParError(1,1);
     }
       h_kFSR_pt_eta_rel->SetBinContent(i+1, j+1,pol_rel[i][j]->GetParameter(0));
       h_kFSR_pt_eta_rel->SetBinError(i+1, j+1, pol_rel[i][j]->GetParError(0));

      double chi2ndf_kFSR_rel = pol_rel[i][j]->GetChisquare() / pol_rel[i][j]->GetNDF();
      h_chi2_kFSR_rel->SetBinContent(i+1, j+1,chi2ndf_kFSR_rel);

     Rel[i][j]->Print(CorrectionObject::_outpath+"plots/kFSR_Pt_eta_"+eta_range2[i]+"_"+eta_range2[i+1]+"_pT_"+pt_range[j]+"_"+pt_range[j+1]+".pdf");
     }
   }

  TCanvas* c1 = new TCanvas();
  m_gStyle->SetOptStat(0);
  h_kFSR_pt_eta_rel->Draw("COLZ");
  c1->Print(CorrectionObject::_outpath+"plots/kFSR_Pt_"+CorrectionObject::_generator_tag+"_2DPlot.pdf");

  TCanvas* c3 = new TCanvas();
  m_gStyle->SetOptStat(0);
  h_chi2_kFSR_rel->Draw("COLZ");
  c3->Print(CorrectionObject::_outpath+"plots/Chi2_kFSR_pT_"+CorrectionObject::_generator_tag+"_2DPlot.pdf");
  //*****************************************************************************************************
  
   TH2D* h_kFSR_pt_eta_mpf = new TH2D("kFSR_pt_eta_mpf","kFSR;|#eta|;p_{T}^{ave}", n_eta-1, eta_bins, n_pt-1, pt_bins);
   TH2D* h_chi2_kFSR_mpf = new TH2D("chi2_kFSR_rel","#chi^{2} kFSR;|#eta|;p_{T}^{ave}", n_eta-1, eta_bins, n_pt-1, pt_bins);

   TF1 *pol_mpf[n_eta-1][n_pt-1];
    for(int i=0; i<n_eta-1; i++){
      for(int j=0; j<n_pt-1; j++){
       pol_mpf[i][j] = new TF1("pol_mpf","pol1",0.14,0.36);
       if(pt_bins[j]<73) continue;
       graph_mpf_r[j][i] = new TGraphErrors(n_alpha,xbin_tgraph,ratio_al_mpf_r[j][i],zero,err_ratio_al_mpf_r[j][i]);
       graph_mpf_r[j][i] = (TGraphErrors*)CleanEmptyPoints(graph_mpf_r[j][i]);
       if(graph_mpf_r[j][i]->GetN()>0) {

       pol_mpf[i][j]->SetParameters(1.5,-0.5);
       graph_mpf_r[j][i]->Fit(pol_mpf[i][j],"RM");
       }
       else {
       pol_mpf[i][j]->SetParameters(-1,-1);
       pol_mpf[i][j]->SetParError(0,1);
       pol_mpf[i][j]->SetParError(1,1);
     }
       h_kFSR_pt_eta_mpf->SetBinContent(i+1, j+1,pol_mpf[i][j]->GetParameter(0));
       h_kFSR_pt_eta_mpf->SetBinError(i+1, j+1, pol_mpf[i][j]->GetParError(0));

      double chi2ndf_kFSR_mpf = pol_mpf[i][j]->GetChisquare() / pol_mpf[i][j]->GetNDF();
      h_chi2_kFSR_mpf->SetBinContent(i+1, j+1,chi2ndf_kFSR_mpf);
     }
   }

  TCanvas* c2 = new TCanvas();
    m_gStyle->SetOptTitle(0);
    h_kFSR_pt_eta_mpf->Draw("COLZ");
  c2->Print(CorrectionObject::_outpath+"plots/kFSR_MPF_"+CorrectionObject::_generator_tag+"_2DPlot.pdf");
  
  TCanvas* c4 = new TCanvas();
  m_gStyle->SetOptStat(0);
  h_chi2_kFSR_mpf->Draw("COLZ");
  c4->Print(CorrectionObject::_outpath+"plots/Chi2_kFSR_MPF_"+CorrectionObject::_generator_tag+"_2DPlot.pdf");
 //**********************************************************************************************






   /* +++++++++++++++++ PLOTS ++++++++++++++++++ */
   //Do the well-known, colorful kFSR linear fit plots
  

   //Create horizontal line for plotting ("ideal value")
   TLine *line = new TLine(alpha_bins[0]-0.01,1,alpha_bins[n_alpha-1]+0.01,1);

   //First set up the output files   
   //Create output .dat file, including the kFSR extrapolation (alpha->0)
   FILE *fp_rel_r, *fp_mpf_r; 
   TH1D *kFSR_MPF, *kFSR_DiJet, *plotkfsr;
   
   cout << "Opening .dat files at: " << CorrectionObject::_outpath+"output/KFSR_MPF_extrapolation.dat" << endl;
   fp_mpf_r = fopen(CorrectionObject::_outpath+"output/KFSR_MPF_extrapolation.dat","w");
   fp_rel_r = fopen(CorrectionObject::_outpath+"output/KFSR_DiJet_extrapolation.dat","w");

   kFSR_MPF = new TH1D("kfsr_mpf","kfsr_mpf", n_eta-1,eta_bins);
   kFSR_DiJet = new TH1D("kfsr_dijet","kfsr_dijet", n_eta-1,eta_bins);
   plotkfsr = new TH1D("kfsr","kfsr", n_eta-1,eta_bins);


   // Start with pT-balance
   //create plots
   TCanvas* a[n_eta-1];
   TString plotname[n_eta-1];
   TF1 *pol1[n_eta-1];
   for (int j=0; j<n_eta-1; j++){ //n_eta-1
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

     // pol1[j] = new TF1("pol1","pol1",0.09,0.38);
     pol1[j] = new TF1("pol1","pol1",0.14,0.36); //0.36
     //     pol1[j] = new TF1("pol1","pol1",0.09,0.42);
 
     pol1[j]->SetParameters(1.5,-0.5);
    
 
     if (multigraph_rel_empty[j]) cout << "Eta bin no. " << j << ", multigraph empty!" << endl;
     else cout << "Eta bin no. " << j << ", multigraph filled!" << endl;
     if(!multigraph_rel_empty[j]){
       cout<<"Set Parameters and Fit!"<<endl;
       pol1[j]->SetParameters(1.5,-0.5);
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

     leg1->SetHeader("p_{T} balance, "+eta_range[j]+"#leq #eta <"+eta_range[j+1]);
     leg1->Draw();

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
     cout << "Printing kFSR plots to " << CorrectionObject::_outpath+"plots/kFSR_Pt_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf" << endl;
     a[j]->Print(CorrectionObject::_outpath+"plots/kFSR_Pt_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");


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
   h_chi2_kFSR_rel->Write();
   outputfile_rel_r->Write();
   cout << "closing output file" << endl;
   outputfile_rel_r->Close();
   cout << "closed output file" << endl;



   // And now MPF results
   //create plots
   TCanvas* b[n_eta-1];
   for (int j=0; j<n_eta-1; j++){
     plotname[j]="mpf_kfsr_diffPt_eta_"+eta_range[j]+"_"+eta_range[j+1];
     b[j] = new TCanvas(plotname[j], plotname[j], 800,700);
     m_gStyle->SetOptTitle(0);

     pTgraph_mpf_r[j]->Draw("AP");
     if(!multigraph_mpf_empty[j]){
       pTgraph_mpf_r[j]->GetYaxis()->SetRangeUser(0.92,1.08);
       pTgraph_mpf_r[j]->GetXaxis()->SetRangeUser(0,alpha_bins[n_alpha-1]+0.01);
       pTgraph_mpf_r[j]->GetYaxis()->SetTitle("(R_{MC}/R_{DATA})");
       pTgraph_mpf_r[j]->GetXaxis()->SetTitle("cut on #alpha");
     }

     //pol1[j] = new TF1("pol1","pol1",0.09,0.36);  //TEST AK4
     pol1[j] = new TF1("pol1","pol1",0.14,0.36);  //TEST AK4
     //pol1[j] = new TF1("pol1","pol1",0.09,0.42);  //TEST AK4
     pol1[j]->SetParameters(0,0);

     if(!multigraph_mpf_empty[j]) pTgraph_mpf_r[j]->Fit(pol1[j],"R");
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

     leg1->SetHeader("MPF, "+eta_range[j]+"#leq #eta <"+eta_range[j+1]);
     leg1->Draw();

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
     cout << "Saving the MPF plots to " << CorrectionObject::_outpath + "plots/kFSR_MPF_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf" << endl;
     b[j]->Print(CorrectionObject::_outpath + "plots/kFSR_MPF_"+CorrectionObject::_generator_tag+"_eta_"+eta_range2[j]+"_"+eta_range2[j+1]+".pdf");

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
   h_chi2_kFSR_mpf->Write();
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
     for(int k=0; k<n_pt-1; k++){
       delete graph_rel_r[k][j];
       delete graph_mpf_r[k][j];
     }
   }
   for(int j=0; j<n_eta-1; j++){
     delete pTgraph_rel_r[j];
     delete pTgraph_mpf_r[j];
   }

   delete leg1;
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
