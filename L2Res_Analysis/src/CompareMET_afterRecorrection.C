//investigation of difference between Modified MET and Recorrected MET
// author: A.Karavdina
// date: 07.11.2018
// Run it with following command:
// root -l -b -q CompareMET_afterRecorrection.C

void CompareMET_afterRecorrection(){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);

 Double_t w = 600;
 Double_t h = 600;

 // bool isWeights = true;
 bool isWeights = false;

 //Files after UHH2/BaconJet (dijet back to back) selection
 TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_Summer16_23Sep2016_V4_noClosure_noJER_slimmedMETsModifiedMET_PFMET_jetRaw50GeV_v3Test_local/";
 TString name = "uhh2.AnalysisModuleRunner.MC.QCD_local_MadGraph_AK4CHS.root";
 TString gl_label = "QCDHT_slimmedMETsModifiedMET";

 /* // TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_Summer16_23Sep2016_V4_noClosure_noJER_slimmedMETsModifiedMET_PFMET_jetRaw50GeV_v3Test/";
 TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_Summer16_23Sep2016_V4_noClosure_noJER_slimmedMETsModifiedMET_PFMET_jetRaw50GeV_v3Test_minJetPt30/";
 TString name = "uhh2.AnalysisModuleRunner.MC.QCD_HT_MadGraph_AK4CHS.root";
 TString gl_label = "QCDHT_slimmedMETsModifiedMET";
 */
 // TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_Summer16_23Sep2016_V4_noClosure_noJER_slimmedMETsModifiedMET_PFMET_jetRaw50GeV/";
 // TString name = "uhh2.AnalysisModuleRunner.MC.QCD_HT_MadGraph_AK4CHS.root";
 // TString gl_label = "QCDHT_slimmedMETsModifiedMET";


 /* TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETs_PFMET_jetRaw15GeV/";
 TString name = "uhh2.AnalysisModuleRunner.MC.QCD_HT_MadGraph_AK4CHS.root";
 TString gl_label = "QCDHT_slimmedMETs"; */


 /* TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETs_PFMET_jetRaw50GeV/";
 TString name = "uhh2.AnalysisModuleRunner.MC.QCD_HT_MadGraph_AK4CHS.root";
 TString gl_label = "QCDHT_slimmedMETs";*/

 // TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETsModifiedMET_PFMET_jetRaw50GeV/";
 // TString name = "uhh2.AnalysisModuleRunner.MC.QCD_HT_MadGraph_AK4CHS.root";
 // TString gl_label = "QCDHT_slimmedMETsModifiedMET";

 /* TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETs_PFMET_jetRaw15GeV/";
 TString name = "uhh2.AnalysisModuleRunner.DATA.DATA_RunF_17Nov17_2017.root";
 TString gl_label = "slimmedMETs"; */

 /* TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETsModifiedMET_PFMET_jetRaw50GeV/";
 TString name = "uhh2.AnalysisModuleRunner.DATA.DATA_RunF_17Nov17_2017.root";
 TString gl_label = "slimmedMETsModifiedMET";*/

 TFile *input = TFile::Open(path+name);
 TTree *TTreeAna = (TTree*)input->Get("AnalysisTree");
 TTreeAna->SetAlias("dMET", "dMET");
 TTreeAna->SetAlias("MET", "MET");
 TTreeAna->SetAlias("alpha", "alpha");
 TTreeAna->SetAlias("probejet_eta", "probejet_eta");

 //vs alpha
 TH2F *dMetvsMET_alpha10 = new TH2F("dMetvsMET_alpha10",";MET;#Delta MET/MET",100,0,500,100,-2.,2.);
 TTreeAna->Project("dMetvsMET_alpha10","dMET/MET:MET","alpha<1.0");
 dMetvsMET_alpha10->FitSlicesY();
 TH1F *dMetvsMET_alpha10_1 = (TH1F*)gDirectory->Get("dMetvsMET_alpha10_1");
 dMetvsMET_alpha10_1->SetMarkerColor(kBlack);
 dMetvsMET_alpha10_1->SetMarkerSize(1.5);
 dMetvsMET_alpha10_1->SetMarkerStyle(21);
 TH2F *dMetvsMET_alpha03 = new TH2F("dMetvsMET_alpha03",";MET;#Delta MET/MET",100,0,500,100,-2.,2.);
 TTreeAna->Project("dMetvsMET_alpha03","dMET/MET:MET","alpha<0.3");
 dMetvsMET_alpha03->FitSlicesY();
 TH1F *dMetvsMET_alpha03_1 = (TH1F*)gDirectory->Get("dMetvsMET_alpha03_1");
 dMetvsMET_alpha03_1->SetMarkerColor(kRed);
 dMetvsMET_alpha03_1->SetMarkerSize(1.5);
 dMetvsMET_alpha03_1->SetMarkerStyle(20);
 TH2F *dMetvsMET_alpha01 = new TH2F("dMetvsMET_alpha01",";MET;#Delta MET/MET",100,0,500,100,-2.,2.);
 TTreeAna->Project("dMetvsMET_alpha01","dMET/MET:MET","alpha<0.1");
 dMetvsMET_alpha01->FitSlicesY();
 TH1F *dMetvsMET_alpha01_1 = (TH1F*)gDirectory->Get("dMetvsMET_alpha01_1");
 dMetvsMET_alpha01_1->SetMarkerColor(kGreen);
 dMetvsMET_alpha01_1->SetMarkerSize(1.8);
 dMetvsMET_alpha01_1->SetMarkerStyle(22);

 TCanvas * c1_alpha = new TCanvas("calpha", "c", w, h);
 TLegend *leg_alpha = new TLegend(0.17,0.67,0.68,0.87);
 leg_alpha->SetHeader(gl_label);
 leg_alpha->AddEntry(dMetvsMET_alpha01_1,"#alpha<0.1","lep");
 leg_alpha->AddEntry(dMetvsMET_alpha03_1,"#alpha<0.3","lep");
 leg_alpha->AddEntry(dMetvsMET_alpha10_1,"#alpha<1.0","lep");
 leg_alpha->SetTextSize(0.045);
 dMetvsMET_alpha01_1->Draw();
 dMetvsMET_alpha01_1->SetTitle("");
 dMetvsMET_alpha01_1->GetYaxis()->SetTitle("(MET_{before} - MET_{after})/MET");
 dMetvsMET_alpha01_1->GetYaxis()->SetRangeUser(-0.5, 0.5);
 dMetvsMET_alpha03_1->Draw("same");
 dMetvsMET_alpha10_1->Draw("same");
 leg_alpha->Draw();
 c1_alpha->SaveAs("dMET_alpha_"+gl_label+".pdf");
 TCanvas * c1_alpha10_2D = new TCanvas("calpha10_2D", "c", w, h);
 dMetvsMET_alpha10->SetTitle("#alpha<1.0");
 dMetvsMET_alpha10->Draw("colz");
 c1_alpha10_2D->SaveAs("2D_dMET_alpha10_"+gl_label+".pdf");
 TCanvas * c1_alpha03_2D = new TCanvas("calpha03_2D", "c", w, h);
 dMetvsMET_alpha03->SetTitle("#alpha<0.3");
 dMetvsMET_alpha03->Draw("colz");
 c1_alpha03_2D->SaveAs("2D_dMET_alpha03_"+gl_label+".pdf");
 TCanvas * c1_alpha01_2D = new TCanvas("calpha01_2D", "c", w, h);
 dMetvsMET_alpha01->SetTitle("#alpha<0.1");
 dMetvsMET_alpha01->Draw("colz");
 c1_alpha01_2D->SaveAs("2D_dMET_alpha01_"+gl_label+".pdf");

 //vs probe jet eta
TH2F *dMetvsMET_probejetEtaBB = new TH2F("dMetvsMET_probejetEtaBB",";MET;#Delta MET/MET",100,0,500,100,-2.,2.);
 TTreeAna->Project("dMetvsMET_probejetEtaBB","dMET/MET:MET","fabs(probejet_eta)<1.3 && alpha<0.3");
 dMetvsMET_probejetEtaBB->FitSlicesY();
 TH1F *dMetvsMET_probejetEtaBB_1 = (TH1F*)gDirectory->Get("dMetvsMET_probejetEtaBB_1");
 dMetvsMET_probejetEtaBB_1->SetMarkerColor(kBlack);
 dMetvsMET_probejetEtaBB_1->SetMarkerSize(1.5);
 dMetvsMET_probejetEtaBB_1->SetMarkerStyle(21);
 TH2F *dMetvsMET_probejetEtaEC = new TH2F("dMetvsMET_probejetEtaEC",";MET;#Delta MET/MET",100,0,500,100,-2.,2.);
 TTreeAna->Project("dMetvsMET_probejetEtaEC","dMET/MET:MET","fabs(probejet_eta)>2.5 && fabs(probejet_eta)<3.0 && alpha<0.3");
 dMetvsMET_probejetEtaEC->FitSlicesY();
 TH1F *dMetvsMET_probejetEtaEC_1 = (TH1F*)gDirectory->Get("dMetvsMET_probejetEtaEC_1");
 dMetvsMET_probejetEtaEC_1->SetMarkerColor(kRed);
 dMetvsMET_probejetEtaEC_1->SetMarkerSize(1.5);
 dMetvsMET_probejetEtaEC_1->SetMarkerStyle(20);
 TH2F *dMetvsMET_probejetEtaHF = new TH2F("dMetvsMET_probejetEtaHF",";MET;#Delta MET/MET",100,0,500,100,-2.,2.);
 TTreeAna->Project("dMetvsMET_probejetEtaHF","dMET/MET:MET","fabs(probejet_eta)>3.2 && alpha<0.3");
 dMetvsMET_probejetEtaHF->FitSlicesY();
 TH1F *dMetvsMET_probejetEtaHF_1 = (TH1F*)gDirectory->Get("dMetvsMET_probejetEtaHF_1");
 dMetvsMET_probejetEtaHF_1->SetMarkerColor(kGreen);
 dMetvsMET_probejetEtaHF_1->SetMarkerSize(1.8);
 dMetvsMET_probejetEtaHF_1->SetMarkerStyle(22);

 TCanvas * c1_probejetEta = new TCanvas("cprobejetEta", "c", w, h);
 TLegend *leg_probejetEta = new TLegend(0.17,0.67,0.68,0.87);
 leg_probejetEta->SetHeader(gl_label);
 leg_probejetEta->AddEntry(dMetvsMET_probejetEtaHF_1,"probejet |#eta|>3.2","lep");
 leg_probejetEta->AddEntry(dMetvsMET_probejetEtaEC_1,"2.5<probejet |#eta|<3.0","lep");
 leg_probejetEta->AddEntry(dMetvsMET_probejetEtaBB_1,"probejet |#eta|<1.3","lep");
 leg_probejetEta->SetTextSize(0.045);
 dMetvsMET_probejetEtaHF_1->Draw();
 dMetvsMET_probejetEtaHF_1->SetTitle("");
 dMetvsMET_probejetEtaHF_1->GetYaxis()->SetTitle("(MET_{before} - MET_{after})/MET");
 dMetvsMET_probejetEtaHF_1->GetYaxis()->SetRangeUser(-0.5, 0.5);
 dMetvsMET_probejetEtaEC_1->Draw("same");
 dMetvsMET_probejetEtaBB_1->Draw("same");
 leg_probejetEta->Draw();
 c1_probejetEta->SaveAs("dMET_ProbejetEta_"+gl_label+".pdf");
 TCanvas * c1_probejetEtaBB_2D = new TCanvas("cprobejetEtaBB_2D", "c", w, h);
 dMetvsMET_probejetEtaBB->SetTitle("probejet |#eta|<1.3");
 dMetvsMET_probejetEtaBB->Draw("colz");
 c1_probejetEtaBB_2D->SaveAs("2D_dMET_probejetEtaBB_"+gl_label+".pdf");
 TCanvas * c1_probejetEtaEC_2D = new TCanvas("cprobejetEtaEC_2D", "c", w, h);
 dMetvsMET_probejetEtaEC->SetTitle("2.5<probejet |#eta|<3.0");
 dMetvsMET_probejetEtaEC->Draw("colz");
 c1_probejetEtaEC_2D->SaveAs("2D_dMET_probejetEtaEC_"+gl_label+".pdf");
 TCanvas * c1_probejetEtaHF_2D = new TCanvas("cprobejetEtaHF_2D", "c", w, h);
 dMetvsMET_probejetEtaHF->SetTitle("probejet |#eta|>3.2");
 dMetvsMET_probejetEtaHF->Draw("colz");
 c1_probejetEtaHF_2D->SaveAs("2D_dMET_probejetEtaHF_"+gl_label+".pdf");

//Plots vs probe jet eta
 TH2F *dMetvsprobejetEta = new TH2F("dMetvsprobejetEta",";probe jet #eta;#Delta MET/MET",20,0,5.2,100,-5.,5.);
 TTreeAna->Project("dMetvsprobejetEta","dMET/MET:fabs(probejet_eta)","alpha<0.3 && MET>50");
 TCanvas * c1_probejetEta_2D = new TCanvas("cprobejetEta_2D", "c", w, h);
 dMetvsprobejetEta->Draw("colz");
 c1_probejetEta_2D->SaveAs("2D_dMET_vs_probejetEta_"+gl_label+".pdf");

 TH2F *dMetvsProbeJet_alpha10 = new TH2F("dMetvsProbeJet_alpha10",";MET;#Delta MET/MET",20,0,5.2,100,-5.,5.);
 TTreeAna->Project("dMetvsProbeJet_alpha10","dMET/MET:fabs(probejet_eta)","alpha<1.0 && MET>50");
 dMetvsProbeJet_alpha10->FitSlicesY();
 TH1F *dMetvsProbeJet_alpha10_1 = (TH1F*)gDirectory->Get("dMetvsProbeJet_alpha10_1");
 dMetvsProbeJet_alpha10_1->SetMarkerColor(kBlack);
 dMetvsProbeJet_alpha10_1->SetMarkerSize(1.5);
 dMetvsProbeJet_alpha10_1->SetMarkerStyle(21);
 TH2F *dMetvsProbeJet_alpha03 = new TH2F("dMetvsProbeJet_alpha03",";MET;#Delta MET/MET",20,0,5.2,100,-5.,5.);
 TTreeAna->Project("dMetvsProbeJet_alpha03","dMET/MET:fabs(probejet_eta)","alpha<0.3 && MET>50");
 dMetvsProbeJet_alpha03->FitSlicesY();
 TH1F *dMetvsProbeJet_alpha03_1 = (TH1F*)gDirectory->Get("dMetvsProbeJet_alpha03_1");
 dMetvsProbeJet_alpha03_1->SetMarkerColor(kRed);
 dMetvsProbeJet_alpha03_1->SetMarkerSize(1.5);
 dMetvsProbeJet_alpha03_1->SetMarkerStyle(20);
 TH2F *dMetvsProbeJet_alpha01 = new TH2F("dMetvsProbeJet_alpha01",";MET;#Delta MET/MET",20,0,5.2,100,-5.,5.);
 TTreeAna->Project("dMetvsProbeJet_alpha01","dMET/MET:fabs(probejet_eta)","alpha<0.1 && MET>50");
 dMetvsProbeJet_alpha01->FitSlicesY();
 TH1F *dMetvsProbeJet_alpha01_1 = (TH1F*)gDirectory->Get("dMetvsProbeJet_alpha01_1");
 dMetvsProbeJet_alpha01_1->SetMarkerColor(kGreen);
 dMetvsProbeJet_alpha01_1->SetMarkerSize(1.8);
 dMetvsProbeJet_alpha01_1->SetMarkerStyle(22);

 TCanvas * c1_eta_alpha = new TCanvas("calpha", "c", w, h);
 TLegend *leg_eta_alpha = new TLegend(0.17,0.67,0.68,0.87);
 leg_eta_alpha->SetHeader(gl_label);
 leg_eta_alpha->AddEntry(dMetvsProbeJet_alpha01_1,"#alpha<0.1","lep");
 leg_eta_alpha->AddEntry(dMetvsProbeJet_alpha03_1,"#alpha<0.3","lep");
 leg_eta_alpha->AddEntry(dMetvsProbeJet_alpha10_1,"#alpha<1.0","lep");
 leg_eta_alpha->SetTextSize(0.045);
 dMetvsProbeJet_alpha01_1->Draw();
 dMetvsProbeJet_alpha01_1->SetTitle("");
 dMetvsProbeJet_alpha01_1->GetYaxis()->SetTitle("(MET_{before} - MET_{after})/MET");
 dMetvsProbeJet_alpha01_1->GetYaxis()->SetRangeUser(-0.5, 0.5);
 dMetvsProbeJet_alpha03_1->Draw("same");
 dMetvsProbeJet_alpha10_1->Draw("same");
 leg_eta_alpha->Draw();
 c1_eta_alpha->SaveAs("dMETvsprobejetEta_alpha_"+gl_label+".pdf");

}
