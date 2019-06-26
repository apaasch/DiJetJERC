//investigation of difference of response value obtained with hist of different width
// author: A.Karavdina
// date: 18.11.2018
// Run it with following command:
// root -l -b -q Compare_response_hist_settings.C

void Compare_response_hist_settings(){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);

 Double_t w = 600;
 Double_t h = 600;

 // bool isWeights = true;
 bool isWeights = false;

 //Files after UHH2/BaconJet (dijet back to back) selection
 TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETs_CHSPFMETjetRaw50GeV_minGenJetPt15_minJetPt15/";

 /*  TString name = "uhh2.AnalysisModuleRunner.MC.QCDPt_pythia8_AK4CHS.root";
     TString gl_label = "QCD_Pythia8"; */
 TString name = "uhh2.AnalysisModuleRunner.MC.QCD_HT_MadGraph_AK4CHS.root";
 TString gl_label = "QCD_MadGraph";

 TFile *input = TFile::Open(path+name);
 TTree *TTreeAna = (TTree*)input->Get("AnalysisTree");
 TTreeAna->SetAlias("Ngenjet", "Ngenjet");
 TTreeAna->SetAlias("Njet", "Njet");
 TTreeAna->SetAlias("MET", "MET");
 TTreeAna->SetAlias("alpha", "alpha");
 TTreeAna->SetAlias("probejet_eta", "probejet_eta");
 TTreeAna->SetAlias("rel_r", "rel_r");
 TTreeAna->SetAlias("mpd_r", "mpf_r");
 TTreeAna->SetAlias("pt_ave", "pt_ave");


 TString condition="weight*(abs(probejet_eta)>2.5 && abs(probejet_eta)<2.65 && pt_ave>520 && pt_ave<648 && alpha<0.3)";

 //1D response hists
 TH1D *Rel_response_old = new TH1D("Rel_response_old",";pt-balance response",100,-1.2,1.2);
 TTreeAna->Project("Rel_response_old","rel_r",condition);
 Rel_response_old->SetMarkerColor(kRed);
 Rel_response_old->SetMarkerSize(1.0);
 Rel_response_old->SetMarkerStyle(21);
 double rel_l_old = Rel_response_old->GetMean();
 char rel_l_old_txt[100];
 sprintf (rel_l_old_txt, "%.3f", rel_l_old);
 TString old_txt="OLD, Mean=";
 old_txt+=rel_l_old_txt;

 TH1D *Rel_response_new = new TH1D("Rel_response_new",";pt-balance response",83,0,2.0);
 TTreeAna->Project("Rel_response_new","rel_r",condition);
 Rel_response_new->SetMarkerColor(kBlue);
 Rel_response_new->SetMarkerSize(1.1);
 Rel_response_new->SetMarkerStyle(20);
 double rel_l_new = Rel_response_new->GetMean();
 char rel_l_new_txt[100];
 sprintf (rel_l_new_txt, "%.3f", rel_l_new);
 TString new_txt="NEW, Mean=";
 new_txt+=rel_l_new_txt;


 TH1D *Rel_response_new2 = new TH1D("Rel_response_new2",";pt-balance response",62,0,1.5);
 TTreeAna->Project("Rel_response_new2","rel_r",condition);
 Rel_response_new2->SetMarkerColor(kGreen);
 Rel_response_new2->SetMarkerSize(1.1);
 Rel_response_new2->SetMarkerStyle(20);
 double rel_l_new2 = Rel_response_new2->GetMean();
 char rel_l_new2_txt[100];
 sprintf (rel_l_new2_txt, "%.3f", rel_l_new2);
 TString new2_txt="NEW2, Mean=";
 new2_txt+=rel_l_new2_txt;

 // TH1D *Rel_response_tmp = new TH1D("Rel_response_tmp",";pt-balance response",167,-2.0,2.0);
 // TTreeAna->Project("Rel_response_tmp","rel_r",condition);
 // Rel_response_tmp->SetMarkerColor(kBlue);
 // // Rel_response_tmp->SetMarkerSize(1.5);
 // // Rel_response_tmp->SetMarkerStyle(21);
 // double rel_l_tmp = Rel_response_tmp->GetMean();





 TCanvas * c1_RelResponse = new TCanvas("cRelResponse", "c", w, h);
 TLegend *leg_RelResponse = new TLegend(0.60,0.71,0.98,0.96);
 leg_RelResponse->SetHeader(gl_label);
 leg_RelResponse->AddEntry(Rel_response_old,old_txt,"lep");
 leg_RelResponse->AddEntry(Rel_response_new,new_txt,"lep");
 leg_RelResponse->AddEntry(Rel_response_new2,new2_txt,"lep");
 leg_RelResponse->SetTextSize(0.04);
 // Rel_response_tmp->Draw();
 Rel_response_new->Draw("P");
 Rel_response_new2->Draw("P SAME");
 Rel_response_old->Draw("P SAME");
 leg_RelResponse->Draw();
 c1_RelResponse->SaveAs("RelResponse_"+gl_label+".pdf");

}
