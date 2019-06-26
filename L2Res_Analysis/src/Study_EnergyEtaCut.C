//Translate kinatic edge cut to pt,ave
// author: A.Karavdina
// date: 13.12.2018
// Run it with following command:
// root -l -b -q Study_EnergyEtaCut.C

void Study_EnergyEtaCut(){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);

 Double_t w = 600;
 Double_t h = 600;
 TString gl_label = "HLT_DiJet";
 //2016
 const int n_pt = 11;
 const double pt_bins[n_pt] = {51, 74, 96, 165, 232, 300, 366, 456, 569, 1000, 2000};
 const int n_pt_HF = 10;
 const double pt_bins_HF[n_pt_HF] = {51, 72, 95, 118, 188, 257, 354, 450, 1000, 2000};

 /* //2017
 const int n_pt = 14;
 const double pt_bins[n_pt] = {51,73,85,97,140,179,240,307,370,434,520,649,1000,2000};
 // const int n_pt_HF = 9;
 // const double pt_bins_HF[n_pt_HF] = {51,73,93,113,176,239,318,1000,2000};
 const int n_pt_HF = 10;
 const double pt_bins_HF[n_pt_HF] = {51,73,93,113,176,239,318,370,1000,2000};
 */
// TString gl_label = "HLT_SingleJet";
// const int n_pt = 15;
// const double pt_bins[n_pt] = {40, 72, 95, 130, 160, 190, 226, 250, 283, 344, 443, 577, 606, 1000,2000};
// const int n_pt_HF = 15;
// const double pt_bins_HF[n_pt_HF] = {40, 72, 95, 130, 160, 190, 226, 250, 283, 344, 443, 577,606, 1000,2000};


const int n_eta = 19;
const double eta_bins[n_eta] ={0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
const double eta_cut = 2.853-1e-5;

 TH2Poly* kinematic_edge = new TH2Poly();
 kinematic_edge->SetName("kinematic_edge");  
 kinematic_edge->SetTitle("Kinematicly allowed bins");
 kinematic_edge->GetXaxis()->SetTitle("|#eta|");
 kinematic_edge->GetYaxis()->SetTitle("p_{T}^{ave}"); 
 kinematic_edge->SetStats(kFALSE);
 int bincount=0;
 for(int i=0; i<n_eta-1; i++){   
   bool eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
     int realbinnumber = kinematic_edge->AddBin(eta_bins[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j], eta_bins[i+1], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1] );
     bincount++;
   }
  }

 for(int i=0; i<n_eta-1; i++){   
   bool eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
   for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
     int fill_me = 1;
     double pt_ave = (eta_cut_bool?pt_bins_HF:pt_bins)[j]+0.01;
     if(pt_ave*cosh(eta_bins[i]+0.01)>3250) fill_me=0;
     kinematic_edge->Fill(eta_bins[i]+0.01, (eta_cut_bool?pt_bins_HF:pt_bins)[j]+0.01, fill_me);
   }
 }

 TCanvas * c1_kin_edge = new TCanvas("ckin_edge", "c", w, h);
  kinematic_edge->Draw("colz");
  //c1_kin_edge->SetLogy();
 c1_kin_edge->SaveAs("kin_edge_"+gl_label+".pdf");

 
TCanvas * c1_kin_edge1D = new TCanvas("ckin_edge1D", "c", w, h);
TF1 *fa1 = new TF1("fa1","3250/cosh(x)",0,7.0);
 fa1->Draw();
 fa1->GetYaxis()->SetTitle("max pt^{ave}, GeV");
 fa1->GetXaxis()->SetTitle("#eta");
 c1_kin_edge1D->SetLogy();
 c1_kin_edge1D->SetGrid();
 c1_kin_edge1D->SaveAs("kin_edge_PtvsEta_1D.pdf");
}
