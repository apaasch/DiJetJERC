//investigation of difference between number of GEN and RECO jets
// author: A.Karavdina
// date: 16.11.2018
// Run it with following command:
// root -l -b -q Compare_jet_multiplicities.C

void Compare_jet_multiplicities(){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);

 Double_t w = 600;
 Double_t h = 600;

 // bool isWeights = true;
 bool isWeights = false;

 //Files after UHH2/BaconJet (dijet back to back) selection
 // TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETs_CHSPFMETjetRaw50GeV_minJetPt0/";
 // TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETs_CHSPFMETjetRaw50GeV_minGenJetPt15_minJetPt15/";
 TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETs_CHSPFMETjetRaw50GeV_minGenJetPt10_minJetPt10/";

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



 //vs Ngenjet
 TH2F *dNjets_Ngenjet_all = new TH2F("dNjets_Ngenjet_all",";N_{GEN jets};N_{GEN jets} - N_{RECO jets}",28,2,30,40,-20,20);
 TTreeAna->Project("dNjets_Ngenjet_all","(Ngenjet-Njet):Ngenjet","weight*(Ngenjet<30. && Njet<30)");
 dNjets_Ngenjet_all->FitSlicesY();
 TH1F *dNjets_Ngenjet_all_1 = (TH1F*)gDirectory->Get("dNjets_Ngenjet_all_1");
 dNjets_Ngenjet_all_1->SetMarkerColor(kBlack);
 dNjets_Ngenjet_all_1->SetMarkerSize(1.5);
 dNjets_Ngenjet_all_1->SetMarkerStyle(21);
 TH2F *dNjets_Ngenjet_EC = new TH2F("dNjets_Ngenjet_EC",";Ngenjets;Ngenjets - Nrecojets",28,2,30,40,-20,20);
 TTreeAna->Project("dNjets_Ngenjet_EC","(Ngenjet-Njet):Ngenjet","weight*(Ngenjet<30. && Njet<30 && fabs(probejet_eta)>2.5 && fabs(probejet_eta)<3.0 && alpha<0.3)");
 dNjets_Ngenjet_EC->FitSlicesY();
 TH1F *dNjets_Ngenjet_EC_1 = (TH1F*)gDirectory->Get("dNjets_Ngenjet_EC_1");
 dNjets_Ngenjet_EC_1->SetMarkerColor(kRed);
 dNjets_Ngenjet_EC_1->SetMarkerSize(1.5);
 dNjets_Ngenjet_EC_1->SetMarkerStyle(20);
 TH2F *dNjets_Ngenjet_HF = new TH2F("dNjets_Ngenjet_HF",";Ngenjets;Ngenjets - Nrecojets",28,2,30,40,-20,20);
 TTreeAna->Project("dNjets_Ngenjet_HF","(Ngenjet-Njet):Ngenjet","weight*(Ngenjet<30. && Njet<30 && fabs(probejet_eta)>3.2)");
 dNjets_Ngenjet_HF->FitSlicesY();
 TH1F *dNjets_Ngenjet_HF_1 = (TH1F*)gDirectory->Get("dNjets_Ngenjet_HF_1");
 dNjets_Ngenjet_HF_1->SetMarkerColor(kGreen);
 dNjets_Ngenjet_HF_1->SetMarkerSize(1.4);
 dNjets_Ngenjet_HF_1->SetMarkerStyle(22);
 TH2F *dNjets_Ngenjet_BB = new TH2F("dNjets_Ngenjet_BB",";Ngenjets;Ngenjets - Nrecojets",28,2,30,40,-20,20);
 TTreeAna->Project("dNjets_Ngenjet_BB","(Ngenjet-Njet):Ngenjet","weight*(Ngenjet<30. && Njet<30 && fabs(probejet_eta)<1.3)");
 dNjets_Ngenjet_BB->FitSlicesY();
 TH1F *dNjets_Ngenjet_BB_1 = (TH1F*)gDirectory->Get("dNjets_Ngenjet_BB_1");
 dNjets_Ngenjet_BB_1->SetMarkerColor(kBlue);
 dNjets_Ngenjet_BB_1->SetMarkerSize(1.8);
 dNjets_Ngenjet_BB_1->SetMarkerStyle(20);


 TCanvas * c1_Njets = new TCanvas("cNjets", "c", w, h);
 TLegend *leg_Njets = new TLegend(0.17,0.21,0.68,0.46);
 leg_Njets->SetHeader(gl_label);
 leg_Njets->AddEntry(dNjets_Ngenjet_BB_1,"probejet |#eta|<1.3","lep");
 leg_Njets->AddEntry(dNjets_Ngenjet_EC_1,"2.5<probejet |#eta|<3.0","lep");
 leg_Njets->AddEntry(dNjets_Ngenjet_HF_1,"probejet |#eta|>3.2","lep");
 leg_Njets->AddEntry(dNjets_Ngenjet_all_1,"0.0<probejet |#eta|<5.2","lep");
 leg_Njets->SetTextSize(0.045);
 dNjets_Ngenjet_BB_1->Draw();
 dNjets_Ngenjet_BB_1->SetTitle("");
 dNjets_Ngenjet_BB_1->GetYaxis()->SetTitle("N_{GEN jets} - N_{RECO jets}");
 dNjets_Ngenjet_BB_1->GetXaxis()->SetTitle("N_{GEN jets}");
 dNjets_Ngenjet_BB_1->GetYaxis()->SetRangeUser(-20, 20);
 dNjets_Ngenjet_BB_1->GetXaxis()->SetRangeUser(2, 10);
 dNjets_Ngenjet_EC_1->Draw("same");
 dNjets_Ngenjet_HF_1->Draw("same");
 dNjets_Ngenjet_all_1->Draw("same");
 leg_Njets->Draw();
 c1_Njets->SaveAs("dNjets_Ngenjets_"+gl_label+".pdf");

 TCanvas * c1_Njets_2D = new TCanvas("cNjets_2D", "c", w, h);
 dNjets_Ngenjet_all->SetTitle("all events");
 dNjets_Ngenjet_all->Draw("colz");
 c1_Njets_2D->SaveAs("2D_dNjets_Ngenjets_"+gl_label+".pdf");

 TCanvas * c1_Njets_2D_EC = new TCanvas("cNjets_2D_EC", "c", w, h);
 dNjets_Ngenjet_EC->SetTitle("2.5<probejet |#eta|<3.0");
 dNjets_Ngenjet_EC->Draw("colz");
 c1_Njets_2D_EC->SaveAs("2D_dNjets_Ngenjets_EC_"+gl_label+".pdf");
 TCanvas * c1_Njets_2D_HF = new TCanvas("cNjets_2D_HF", "c", w, h);
 dNjets_Ngenjet_HF->SetTitle("probejet |#eta|>3.2");
 dNjets_Ngenjet_HF->Draw("colz");
 c1_Njets_2D_HF->SaveAs("2D_dNjets_Ngenjets_HF_"+gl_label+".pdf");

}
