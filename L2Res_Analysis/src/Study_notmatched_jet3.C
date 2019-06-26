//investigation of 3rd jet, which was not matched to parton in MadGraph simulation
// author: A.Karavdina
// date: 16.11.2018
// Run it with following command:
// root -l -b -q Study_notmatched_jet3.C

void Study_notmatched_jet3(){
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
 TString path = "/nfs/dust/cms/user/karavdia/JERC/17Nov2017/Residuals/Run17DEF_Data_JEC_V23_noClosure_noJER_slimmedMETs_CHSPFMETjetRaw50GeV_minGenJetPt15_minJetPt15/";

 /* TString name = "uhh2.AnalysisModuleRunner.MC.QCDPt_pythia8_AK4CHS.root";
 // TString gl_label = "QCD_Pythia8_ptAveAbove100";
 TString gl_label = "QCD_Pythia8";
 */

 TString name = "uhh2.AnalysisModuleRunner.MC.QCD_HT_MadGraph_AK4CHS.root";
 //TString gl_label = "QCD_MadGraph";
 //TString gl_label = "QCD_MadGraph_3rdJetnotEC";
 // TString gl_label = "QCD_MadGraph_ptAveAbove200";


 TFile *input = TFile::Open(path+name);
 TTree *TTreeAna = (TTree*)input->Get("AnalysisTree");
 TTreeAna->SetAlias("Ngenjet", "Ngenjet");
 TTreeAna->SetAlias("Njet", "Njet");
 TTreeAna->SetAlias("MET", "MET");
 TTreeAna->SetAlias("alpha", "alpha");
 TTreeAna->SetAlias("probejet_eta", "probejet_eta");
 TTreeAna->SetAlias("flavor3rdjet", "flavor3rdjet");

 // TString condition="weight*(Ngenjet<30. && Njet<30 && flavor3rdjet<0 && pt_ave>0 && (fabs(jet3_eta)<2.5 || fabs(jet3_eta)>3.0))";

 /* // TString condition="weight*(Ngenjet<30. && Njet<30 && flavor3rdjet>0 && pt_ave>0)";
 TString condition="weight*(Ngenjet<30. && Njet<30 && flavor3rdjet>0 && pt_ave>0 && alpha<0.4)";
 TString plotTitel = "events with matched 3rd jet";//Title of 2D plot
 // TString gl_label = "QCD_MadGraph_3rdJetMATCHED"; //name of pdf file
 TString gl_label = "QCD_MadGraph_3rdJetMATCHED_alpha04"; //name of pdf file
 */
 // TString condition="weight*(Ngenjet<30. && Njet<30 && flavor3rdjet<0 && pt_ave>0)";
 // TString gl_label = "QCD_MadGraph"; //name of pdf file
 // TString condition="weight*(Ngenjet<30. && Njet<30 && flavor3rdjet<0 && pt_ave>0 && alpha<0.4)";
 // TString gl_label = "QCD_MadGraph_alpha04"; //name of pdf file
 TString plotTitel = "events with not matched 3rd jet";//Title of 2D plot
 TString condition="weight*(Ngenjet<30. && Njet<30 && flavor3rdjet<0 && pt_ave>300)";
 TString gl_label = "QCD_MadGraph_ptAveAbove300";


 //vs Ngenjet
 TH2F *dNjets_Ngenjet_all = new TH2F("dNjets_Ngenjet_all",";N_{GEN jets};N_{GEN jets} - N_{RECO jets}",28,2,30,40,-20,20);
 TTreeAna->Project("dNjets_Ngenjet_all","Ngenjet-Njet:Ngenjet",condition);

 TCanvas * c1_Njets_2D = new TCanvas("cNjets_2D", "c", w, h);
 dNjets_Ngenjet_all->SetTitle(plotTitel);
 dNjets_Ngenjet_all->Draw("colz");
 c1_Njets_2D->SaveAs("2D_dNjets_Ngenjets_notMatched3rdJet_"+gl_label+".pdf");

 TH2F *Jet3eta_Jet3pt_all = new TH2F("Jet3eta_Jet3pt_all",";jet3 #eta;jet3 pt, GeV",20,-5.2,5.2,20,0,100);
 TTreeAna->Project("Jet3eta_Jet3pt_all","jet3_pt:jet3_eta",condition);
 TCanvas * c1_Jet3eta_Jet3pt_2D = new TCanvas("cJet3eta_Jet3pt_2D", "c", w, h);
 Jet3eta_Jet3pt_all->SetTitle(plotTitel);
 Jet3eta_Jet3pt_all->Draw("colz");
 c1_Jet3eta_Jet3pt_2D->SaveAs("2D_Jet3eta_Jet3pt_notMatched3rdJet_"+gl_label+".pdf");

 TH2F *JetProbeeta_JetProbept_all = new TH2F("JetProbeeta_JetProbept_all",";probe jet #eta;probe jet pt, GeV",20,-5.2,5.2,150,0,600);
 TTreeAna->Project("JetProbeeta_JetProbept_all","probejet_pt:probejet_eta",condition);
 TCanvas * c1_JetProbeeta_JetProbept_2D = new TCanvas("cJetProbeeta_JetProbept_2D", "c", w, h);
 JetProbeeta_JetProbept_all->SetTitle(plotTitel);
 JetProbeeta_JetProbept_all->Draw("colz");
 c1_JetProbeeta_JetProbept_2D->SaveAs("2D_JetProbeeta_JetProbept_notMatched3rdJet_"+gl_label+".pdf");

 TH2F *JetProbeeta_JetAvept_all = new TH2F("JetProbeeta_JetAvept_all",";probe jet #eta;pt_{ave}, GeV",20,-5.2,5.2,50,0,500);
 TTreeAna->Project("JetProbeeta_JetAvept_all","pt_ave:probejet_eta",condition);
 TCanvas * c1_JetProbeeta_JetAvept_2D = new TCanvas("cJetProbeeta_JetAvept_2D", "c", w, h);
 JetProbeeta_JetAvept_all->SetTitle(plotTitel);
 JetProbeeta_JetAvept_all->Draw("colz");
 c1_JetProbeeta_JetAvept_2D->SaveAs("2D_JetProbeeta_JetAvept_notMatched3rdJet_"+gl_label+".pdf");

 TH2F *JetProbeeta_Jet3eta_all = new TH2F("JetProbeeta_Jet3eta_all",";probe jet #eta;jet3 #eta",20,-5.2,5.2,20,-5.2,5.2);
 TTreeAna->Project("JetProbeeta_Jet3eta_all","jet3_eta:probejet_eta",condition);
 TCanvas * c1_JetProbeeta_Jet3eta_2D = new TCanvas("cJetProbeeta_Jet3eta_2D", "c", w, h);
 JetProbeeta_Jet3eta_all->SetTitle(plotTitel);
 JetProbeeta_Jet3eta_all->Draw("colz");
 c1_JetProbeeta_Jet3eta_2D->SaveAs("2D_JetProbeeta_Jet3eta_notMatched3rdJet_"+gl_label+".pdf");

 TH2F *JetAvept_alpha_all = new TH2F("JetAvept_alpha_all",";pt_{ave}, GeV;#alpha",50,0,500,10,0,1.0);
 TTreeAna->Project("JetAvept_alpha_all","alpha:pt_ave",condition);
 TCanvas * c1_JetAvept_alpha_2D = new TCanvas("cJetAvept_alpha_2D", "c", w, h);
 JetAvept_alpha_all->SetTitle(plotTitel);
 JetAvept_alpha_all->Draw("colz");
 c1_JetAvept_alpha_2D->SaveAs("2D_JetAvept_alpha_notMatched3rdJet_"+gl_label+".pdf");


}
