// investigation of difference between output trees after ntuplewriter or UHH2 AnalysisModule
// At the moment works only for variables in single branches, i.e does not work for variables from object collections/vectors
// author: A.Karavdina
// date: 03.07.2019
// Run it with following command:
// root -l -b -q CompareTrees_80Xvs102X.C\(\"MET\",\"met_pt\",\"slimmedMETs.m_pt\",true,90,0,900,10\)


void CompareTrees_80Xvs102X(TString fileName="QCD", TString var_name="MET",TString branch_name_old="met_pt",TString branch_name_new="met_pt", bool isFloat=false, int nbins=10, double minx=0,double maxx=10, int maxevents=100){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);

  double w = 600;
  double h = 600;

  // bool isWeights = true;
  bool isWeights = false;

  int nunmatched = 0;

  //Files after UHH2
  // TString path_old = "/nfs/dust/cms/user/karavdia/sframe_all/JER2018_addEvIDinfo_woPUid_addEnergyEtaCut_fixSortJets_DATA/Autumn18_V15/AK4CHS/";
  TString path_old = "/nfs/dust/cms/user/karavdia/sframe_all/JER2018_addEvIDinfo_woPUid_addEnergyEtaCut_fixSortJets_"; path_old+=(fileName.Contains("Run"))? "DATA" :"QCD"; path_old+="/Autumn18_V13h/AK4CHS/";
  // TString path_old = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/DiJetJERC/";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT100to200.root";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT1000to1500.root";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT1500to2000.root";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT2000toInf.root";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT200to300.root";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT300to500.root";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT500to700.root";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT50to100.root";
  // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCDHT700to1000.root";

  // TString name_old = "uhh2.AnalysisModuleRunner.DATA.DATA_RunA.root";
  // TString name_old = "test_QCD.root";
  TString name_old = "uhh2.AnalysisModuleRunner"; name_old+= (fileName.Contains("Run"))? ".DATA.DATA_" : ".MC."; name_old+= fileName; name_old+= ".root";


  TString path_new = "/nfs/dust/cms/user/amalara//sframe_all//DiJetJERC_DiJetHLT_"; path_new+=(fileName.Contains("Run"))? "DATA" :"QCD"; path_new+="/Autumn18_V15/AK4CHS/";
  // TString path_new = "/nfs/dust/cms/user/amalara//sframe_all//DiJetJERC_DiJetHLT_DATA/Autumn18_V15/AK4CHS/";
  // TString path_new = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/DiJetJERC/";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT100to200.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT1000to1500.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT1500to2000.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT2000toInf.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT200to300.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT300to500.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT500to700.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT50to100.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCDHT700to1000.root";
  // TString name_new = "uhh2.AnalysisModuleRunner.DATA.DATA_RunA_2018.root";
  // TString name_new = "testA_QCD.root";
  TString name_new = "uhh2.AnalysisModuleRunner"; name_new+= (fileName.Contains("Run"))? ".DATA.DATA_" : ".MC."; name_new+= fileName; name_new+= (fileName.Contains("Run"))? "_2018.root" : ".root";


  TString gl_label_old = "JERSF";
  TString gl_label_new = "DiJet";

  cout << path_old+name_old << endl;
  cout << path_new+name_new << endl;

  TFile *input_old = TFile::Open(path_old+name_old);
  TTree *TTreeAna_old = (TTree*)input_old->Get("AnalysisTree");
  TFile *input_new = TFile::Open(path_new+name_new);
  TTree *TTreeAna_new = (TTree*)input_new->Get("AnalysisTree");

  // int eventID,feventID;
  long long eventID,feventID;
  int run,frun;
  int lumi,flumi;
  Float_t var_val_f,fval_val_f;
  Int_t var_val_i,fval_val_i;

  Float_t var_asy,fval_asy;
  float barreljet_eta, probejet_eta;

  TTreeAna_old->SetBranchAddress("eventID",&eventID);
  TTreeAna_new->SetBranchAddress("eventID",&feventID);
  TTreeAna_old->SetBranchAddress("run",&run);
  TTreeAna_new->SetBranchAddress("run",&frun);
  TTreeAna_old->SetBranchAddress("lumi_sec",&lumi);
  TTreeAna_new->SetBranchAddress("lumi_sec",&flumi);
  TTreeAna_old->SetBranchAddress("asymmetry",&var_asy);
  TTreeAna_new->SetBranchAddress("probejet_eta", &barreljet_eta);
  TTreeAna_new->SetBranchAddress("barreljet_eta",&probejet_eta);

  TTreeAna_new->BuildIndex("eventID");
  TTreeAna_old->BuildIndex("eventID");

  if(isFloat){
    cout<<"branch_name = "<<branch_name_old.Data()<<" ";
    cout<<TTreeAna_old->GetListOfBranches()->FindObject(branch_name_old)<<endl;
    cout<<"branch_name = "<<branch_name_new.Data()<<" ";
    cout<<TTreeAna_new->GetListOfBranches()->FindObject(branch_name_new)<<endl;
    TTreeAna_old->SetBranchAddress(branch_name_old,&var_val_f);
    TTreeAna_new->SetBranchAddress(branch_name_new,&fval_val_f);
  }
  else{
    cout<<"branch_name = "<<branch_name_old.Data()<<" ";
    cout<<TTreeAna_old->GetListOfBranches()->FindObject(branch_name_old)<<endl;
    TTreeAna_old->SetBranchAddress(branch_name_old,&var_val_i);
    TTreeAna_new->SetBranchAddress(branch_name_new,&fval_val_i);
  }

  TTreeAna_old->AddFriend(TTreeAna_new);


  TH2D *hist_comp = new TH2D("h"+var_name,var_name+";"+gl_label_old+";"+gl_label_new,nbins,minx,maxx,nbins,minx,maxx);

  int nentries = TTreeAna_old->GetEntries();
  int nok = 0;
  if(maxevents>0) nentries = maxevents;
  for (Long64_t i=0;i<nentries;i++) {
    TTreeAna_old->GetEntry(i);
    // cout << TMath::Abs(barreljet_eta) << " " << TMath::Abs(probejet_eta) << endl;
    if (!(TMath::Abs(probejet_eta)>2.85 && TMath::Abs(barreljet_eta)<2.85)) continue;
    // cout << TMath::Abs(barreljet_eta) << " " << TMath::Abs(probejet_eta) << endl;
    // cout<<"feventID = "<<feventID<<end;
    // if (var_val_f<0.6) continue;
    if(feventID==eventID && lumi!=flumi) {nunmatched++; cout<<"lumi = "<<lumi<<" flumi ="<<flumi<<endl;}
    if (feventID==eventID && run==frun && lumi==flumi) {
      //     cout<<"lumi = "<<lumi<<" flumi ="<<flumi<<endl;
      // var_val_i = TMath::Abs(var_val_i);
      // fval_val_i = TMath::Abs(fval_val_i);
      // var_val_f = TMath::Abs(var_val_f);
      // fval_val_f = TMath::Abs(fval_val_f);
      nok++;
      if(isFloat){
        hist_comp->Fill(var_val_f,fval_val_f);
        // if(fabs(var_val_f-fval_val_f)/var_val_f>0.5) cout<<"Difference above 50% for eventID="<<eventID<<", run="<<run<<" lumi="<<lumi<<endl;
      }
      else{
        hist_comp->Fill(var_val_i,fval_val_i);
        if(double(fabs(var_val_i-fval_val_i))/var_val_i>0.5) cout<<"Difference above 50% for eventID="<<eventID<<", run="<<run<<" lumi="<<lumi<<endl;
      }
    }
  }
  printf("nok = %d, fentries=%lld\n",nok,TTreeAna_new->GetEntries());

  TCanvas * c1_canvas = new TCanvas("canvas", "c", w, h);
  c1_canvas->SetLogz(1);
  hist_comp->Draw("colz");
  TLine *line = new TLine(minx,minx,maxx,maxx);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->Draw();
  // c1_canvas->SaveAs("CompareTrees_PUid/Compare_"+branch_name+"__"+gl_label_old+"_vs_"+gl_label_new+".pdf");
  c1_canvas->SaveAs("./Compare_"+var_name+"_"+fileName+"_"+gl_label_old+"_vs_"+gl_label_new+".pdf");
  cout<< "nunmatched: " << nunmatched << endl;

}
