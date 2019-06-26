// Plot results of compression test results
// author: A.Karavdina
// date: 28.05.2019
// Run it with following command:
// root -l -b -q Plot_compressionResults.C

void Plot_compressionResults(){
  gStyle->SetOptStat(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(1.25);
 Double_t w = 600;
 Double_t h = 600;
 TString gl_label = "Real time";
 //2016
 const int n_Ev_bins = 6;
 const double n_Ev[n_Ev_bins] = {1000,1e4,5e4,1e5,2.5e5,5e5};
 const double n_Ev_err[n_Ev_bins] = {0,0,0,0,0,0};

 const double real_time_ZLIB1[n_Ev_bins] = {15.926,60.466,264.992,535.652,1300.884,2557.804};
 const double real_time_err_ZLIB1[n_Ev_bins] = {3.634,1.213,5.153,36.932,55.901,46.275};
 const double cpu_time_ZLIB1[n_Ev_bins] = {13.462, 59.044,261.564,517.324,1269.540,2528.848};
 const double cpu_time_err_ZLIB1[n_Ev_bins] = {0.031,0.100,0.750,9.562,2.019,7.152};

 const double real_time_LZMA9[n_Ev_bins] = {17.996,66.982,282.918,549.98,1368.194,2718.480};
 const double real_time_err_LZMA9[n_Ev_bins] = {5.656,6.574,7.294,6.665,41.303,58.994};
 const double cpu_time_LZMA9[n_Ev_bins] = {14.628,62.874,278.404,544.528,1344.784,2683.846};
 const double cpu_time_err_LZMA9[n_Ev_bins] = {1e-4,0.089,3.425,2.481,2.602,6.941};

 double rel_diff_real_time[n_Ev_bins];
 double rel_diff_cpu_time[n_Ev_bins];

 for(int i=0;i<n_Ev_bins; i++){
   rel_diff_real_time[i] = 100.*(real_time_LZMA9[i]-real_time_ZLIB1[i])/real_time_ZLIB1[i];
   rel_diff_cpu_time[i] = 100.*(cpu_time_LZMA9[i]-cpu_time_ZLIB1[i])/cpu_time_ZLIB1[i];
 }

 TGraphErrors* gr_real_time_ZLIB1 = new TGraphErrors(n_Ev_bins,n_Ev,real_time_ZLIB1,n_Ev_err,real_time_err_ZLIB1);
 gr_real_time_ZLIB1->SetMarkerStyle(20);
 gr_real_time_ZLIB1->SetMarkerColor(kBlack);
 gr_real_time_ZLIB1->SetMarkerSize(1.2);

 TGraphErrors* gr_real_time_LZMA9 = new TGraphErrors(n_Ev_bins,n_Ev,real_time_LZMA9,n_Ev_err,real_time_err_LZMA9);
 gr_real_time_LZMA9->SetMarkerStyle(21);
 gr_real_time_LZMA9->SetMarkerColor(kRed);
 gr_real_time_LZMA9->SetMarkerSize(1.2);

 auto legend = new TLegend(0.20,0.8,0.32,0.9);
 legend->SetHeader("");
 legend->AddEntry(gr_real_time_ZLIB1,"ZLIB1","p");
 legend->AddEntry(gr_real_time_LZMA9,"LZMA9","p");


 TCanvas * c1_realTime = new TCanvas("crealTime", "c", w, h);
 c1_realTime->SetLogx(); 
 TMultiGraph *mg_real_time = new TMultiGraph();
 mg_real_time->Add(gr_real_time_LZMA9,"p");
 mg_real_time->Add(gr_real_time_ZLIB1,"p");
 mg_real_time->Draw("a");
 mg_real_time->GetXaxis()->SetTitle("N events");
 mg_real_time->GetYaxis()->SetTitle("Real time, s");
 legend->Draw();
 c1_realTime->SaveAs("RealTime_abs.pdf");

 TGraphErrors* gr_cpu_time_ZLIB1 = new TGraphErrors(n_Ev_bins,n_Ev,cpu_time_ZLIB1,n_Ev_err,cpu_time_err_ZLIB1);
 gr_cpu_time_ZLIB1->SetMarkerStyle(20);
 gr_cpu_time_ZLIB1->SetMarkerColor(kBlack);
 gr_cpu_time_ZLIB1->SetMarkerSize(1.2);

 TGraphErrors* gr_cpu_time_LZMA9 = new TGraphErrors(n_Ev_bins,n_Ev,cpu_time_LZMA9,n_Ev_err,cpu_time_err_LZMA9);
 gr_cpu_time_LZMA9->SetMarkerStyle(21);
 gr_cpu_time_LZMA9->SetMarkerColor(kRed);
 gr_cpu_time_LZMA9->SetMarkerSize(1.2);


 TCanvas * c1_cpuTime = new TCanvas("ccpuTime", "c", w, h);
 c1_cpuTime->SetLogx(); 
 TMultiGraph *mg_cpu_time = new TMultiGraph();
 mg_cpu_time->Add(gr_cpu_time_LZMA9,"p");
 mg_cpu_time->Add(gr_cpu_time_ZLIB1,"p");
 mg_cpu_time->Draw("a");
 mg_cpu_time->GetXaxis()->SetTitle("N events");
 mg_cpu_time->GetYaxis()->SetTitle("CPU time, s");
 legend->Draw();
 c1_cpuTime->SaveAs("CpuTime_abs.pdf");

 TCanvas * c1_realTime_diff = new TCanvas("crealTime_diff", "c", w, h);
 c1_realTime_diff->SetLogx();
 TGraphErrors* gr_real_time_diff = new TGraphErrors(n_Ev_bins,n_Ev,rel_diff_real_time,n_Ev_err,n_Ev_err);
 gr_real_time_diff->SetMarkerStyle(47);
 gr_real_time_diff->SetMarkerColor(kBlue);
 gr_real_time_diff->SetMarkerSize(1.5);
 gr_real_time_diff->Draw("AP");
 gr_real_time_diff->SetTitle("");
 gr_real_time_diff->GetXaxis()->SetTitle("N events");
 gr_real_time_diff->GetYaxis()->SetTitle("REAL (LZMA9-ZLIB1)/ZLIB1, %");
 gr_real_time_diff->GetYaxis()->SetRangeUser(0,20);
 c1_realTime_diff->SaveAs("RealTime_rel.pdf");

 TCanvas * c1_cpuTime_diff = new TCanvas("ccpuTime_diff", "c", w, h);
 c1_cpuTime_diff->SetLogx();
 TGraphErrors* gr_cpu_time_diff = new TGraphErrors(n_Ev_bins,n_Ev,rel_diff_cpu_time,n_Ev_err,n_Ev_err);
 gr_cpu_time_diff->SetMarkerStyle(47);
 gr_cpu_time_diff->SetMarkerColor(kBlue);
 gr_cpu_time_diff->SetMarkerSize(1.5);
 gr_cpu_time_diff->Draw("AP");
 gr_cpu_time_diff->SetTitle("");
 gr_cpu_time_diff->GetXaxis()->SetTitle("N events");
 gr_cpu_time_diff->GetYaxis()->SetTitle("CPU (LZMA9-ZLIB1)/ZLIB1, %");
 gr_cpu_time_diff->GetYaxis()->SetRangeUser(0,20);
 c1_cpuTime_diff->SaveAs("CpuTime_rel.pdf");



 const double RSS_ZLIB1[n_Ev_bins] = {1e-6*1101420,1e-6*1125376,1e-6*1155036,1e-6*1162452,1e-6*1262452,1e-6*1394364};
 const double RSS_LZMA9[n_Ev_bins] = {1e-6*1108656,1e-6*1135596,1e-6*1170440,1e-6*1181156,1e-6*1281156,1e-6*1414796};
 TGraphErrors* gr_rss_ZLIB1 = new TGraphErrors(n_Ev_bins,n_Ev,RSS_ZLIB1,n_Ev_err,n_Ev_err);
 TGraphErrors* gr_rss_LZMA9 = new TGraphErrors(n_Ev_bins,n_Ev,RSS_LZMA9,n_Ev_err,n_Ev_err);
 gr_rss_ZLIB1->SetMarkerStyle(20);
 gr_rss_ZLIB1->SetMarkerColor(kBlack);
 gr_rss_ZLIB1->SetMarkerSize(1.2);
 gr_rss_LZMA9->SetMarkerStyle(21);
 gr_rss_LZMA9->SetMarkerColor(kRed);
 gr_rss_LZMA9->SetMarkerSize(1.2);
 TCanvas * c1_rss = new TCanvas("crss", "c", w, h);
 c1_rss->SetLogx();
 TMultiGraph *mg_rss = new TMultiGraph();
 mg_rss->Add(gr_rss_LZMA9,"p");
 mg_rss->Add(gr_rss_ZLIB1,"p");
 mg_rss->Draw("a");
 mg_rss->GetXaxis()->SetTitle("N events");
 mg_rss->GetYaxis()->SetTitle("RSS, GB");
 legend->Draw();
 c1_rss->SaveAs("RSS.pdf");
}

// const int n_eta = 19;
// const double eta_bins[n_eta] ={0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
// const double eta_cut = 2.853-1e-5;

//  TH2Poly* kinematic_edge = new TH2Poly();
//  kinematic_edge->SetName("kinematic_edge");  
//  kinematic_edge->SetTitle("Kinematicly allowed bins");
//  kinematic_edge->GetXaxis()->SetTitle("|#eta|");
//  kinematic_edge->GetYaxis()->SetTitle("p_{T}^{ave}"); 
//  kinematic_edge->SetStats(kFALSE);
//  int bincount=0;
//  for(int i=0; i<n_eta-1; i++){   
//    bool eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
//    for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
//      int realbinnumber = kinematic_edge->AddBin(eta_bins[i], (eta_cut_bool?pt_bins_HF:pt_bins)[j], eta_bins[i+1], (eta_cut_bool?pt_bins_HF:pt_bins)[j+1] );
//      bincount++;
//    }
//   }

//  for(int i=0; i<n_eta-1; i++){   
//    bool eta_cut_bool = fabs(eta_bins[i])>eta_cut;     
//    for(int j= 0 ; j <  ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ) ; j++ ){
//      int fill_me = 1;
//      double pt_ave = (eta_cut_bool?pt_bins_HF:pt_bins)[j]+0.01;
//      if(pt_ave*cosh(eta_bins[i]+0.01)>3250) fill_me=0;
//      kinematic_edge->Fill(eta_bins[i]+0.01, (eta_cut_bool?pt_bins_HF:pt_bins)[j]+0.01, fill_me);
//    }
//  }

//  TCanvas * c1_kin_edge = new TCanvas("ckin_edge", "c", w, h);
//   kinematic_edge->Draw("colz");
//   //c1_kin_edge->SetLogy();
//  c1_kin_edge->SaveAs("kin_edge_"+gl_label+".pdf");

 
// TCanvas * c1_kin_edge1D = new TCanvas("ckin_edge1D", "c", w, h);
// TF1 *fa1 = new TF1("fa1","3250/cosh(x)",0,7.0);
//  fa1->Draw();
//  fa1->GetYaxis()->SetTitle("max pt^{ave}, GeV");
//  fa1->GetXaxis()->SetTitle("#eta");
//  c1_kin_edge1D->SetLogy();
//  c1_kin_edge1D->SetGrid();
//  c1_kin_edge1D->SaveAs("kin_edge_PtvsEta_1D.pdf");
// }
