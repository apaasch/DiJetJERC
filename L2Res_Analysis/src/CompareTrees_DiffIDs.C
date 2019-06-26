// investigation of difference between output trees after ntuplewriter or UHH2 AnalysisModule (JERSF)
// implemented for PUid checks,  but might be useful in general
// At the moment works only for variables in single branches, i.e does not work for variables from object collections/vectors
// author: A.Karavdina
// date: 27.05.2019
// Run it with following command:
// root -l -b -q CompareTrees_DiffIDs.C\(\"Njet\",\"Njet\",false,50,0,50,10\)


void CompareTrees_DiffIDs(TString var_name="Ngenjet",TString branch_name="Ngenjet", bool isFloat=false, int nbins=10, double minx=0,double maxx=10, int maxevents=100){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);

 double w = 600;
 double h = 600;

 // bool isWeights = true;
 bool isWeights = false;

 //Files after UHH2
 // TString path_old = "/nfs/dust/cms/user/karavdia/sframe_all/";
 // TString name_old = "JER2018_woPUid_QCD/Autumn18_V13h/AK4CHS/uhh2.AnalysisModuleRunner.MC.QCDHT.root";
 // TString path_new = "/nfs/dust/cms/user/karavdia/sframe_all/";
 // TString name_new = "JER2018_wPUid_QCD/Autumn18_V13h/AK4CHS/uhh2.AnalysisModuleRunner.MC.QCDHT.root";

 // TString path_old = "/nfs/dust/cms/user/karavdia/sframe_all/JER2018_addEvIDinfo_woPUid_addEnergyEtaCut_DATA/Autumn18_V13h/AK4CHS/";
 // TString name_old = "uhh2.AnalysisModuleRunner.DATA.DATA_RunABC.root";
 // TString path_new = "/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V13h_L1L2L3_L2L3Res_NoJERsmear_NoPUid_Lead2jetsPassJetID_addJERSMbool/";
 // TString name_new = "uhh2.AnalysisModuleRunner.DATA.DATA_RunABC_2018_AK4CHS.root";

 TString path_old = "/nfs/dust/cms/user/karavdia/sframe_all/JER2018_addEvIDinfo_woPUid_addEnergyEtaCut_fixSortJets_DATA/Autumn18_V13h/AK4CHS/";
 TString name_old = "uhh2.AnalysisModuleRunner.DATA.DATA_RunABC.root";
 // TString name_old = "uhh2.AnalysisModuleRunner.DATA.DATA_RunC.root";
 TString path_new = "/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V13h_L1L2L3_L2L3Res_NoJERsmear_NoPUid_Lead2jetsPassJetID_addJERSMbool_fixJetSorting/";
 TString name_new = "uhh2.AnalysisModuleRunner.DATA.DATA_RunABC_2018_AK4CHS.root";
 //TString name_new = "uhh2.AnalysisModuleRunner.DATA.DATA_RunC_2018_AK4CHS.root";

 // TString gl_label_old = "woPUid";
 // TString gl_label_new = "withPUid";

 TString gl_label_old = "JERSF";
 TString gl_label_new = "L2Res";

 TFile *input_old = TFile::Open(path_old+name_old);
 TTree *TTreeAna_old = (TTree*)input_old->Get("AnalysisTree");
 TFile *input_new = TFile::Open(path_new+name_new);
 TTree *TTreeAna_new = (TTree*)input_new->Get("AnalysisTree");

 // int eventID,feventID;
 long long eventID,feventID;
 int run,frun;
 Float_t var_val_f,fval_val_f;
 Int_t var_val_i,fval_val_i;

 TTreeAna_old->SetBranchAddress("eventID",&eventID);
 TTreeAna_new->SetBranchAddress("eventID",&feventID);
 TTreeAna_old->SetBranchAddress("run",&run);
 TTreeAna_new->SetBranchAddress("run",&frun);
 TTreeAna_new->BuildIndex("eventID");
 TTreeAna_old->BuildIndex("eventID");

 if(isFloat){
   cout<<"branch_name = "<<branch_name.Data()<<" ";
   cout<<TTreeAna_old->GetListOfBranches()->FindObject(branch_name)<<endl;
   TTreeAna_old->SetBranchAddress(branch_name,&var_val_f);
   TTreeAna_new->SetBranchAddress(branch_name,&fval_val_f);
 }
 else{
   cout<<"branch_name = "<<branch_name.Data()<<" ";
   cout<<TTreeAna_old->GetListOfBranches()->FindObject(branch_name)<<endl;
   TTreeAna_old->SetBranchAddress(branch_name,&var_val_i);
   TTreeAna_new->SetBranchAddress(branch_name,&fval_val_i);
 }

 TTreeAna_old->AddFriend(TTreeAna_new);


 TH2D *hist_comp = new TH2D("h"+branch_name,var_name+";"+gl_label_old+";"+gl_label_new,nbins,minx,maxx,nbins,minx,maxx);

 int nentries = TTreeAna_old->GetEntries();
 int nok = 0;
 if(maxevents>0) nentries = maxevents;
 for (Long64_t i=0;i<nentries;i++) {
   TTreeAna_old->GetEntry(i);
   if (feventID==eventID && run==frun) {
     nok++;
     if(isFloat){
       hist_comp->Fill(var_val_f,fval_val_f);
       if(fabs(var_val_f-fval_val_f)/var_val_f>0.1) cout<<"Difference above 10% for eventID="<<eventID<<", run="<<run<<endl;
     }
     else{
       hist_comp->Fill(var_val_i,fval_val_i);
       if(double(fabs(var_val_i-fval_val_i))/var_val_i>0.1) cout<<"Difference above 10% for eventID="<<eventID<<", run="<<run<<endl;
     }
   }
 }
 printf("nok = %d, fentries=%lld\n",nok,TTreeAna_new->GetEntries());

 TCanvas * c1_canvas = new TCanvas("canvas", "c", w, h);
 hist_comp->Draw("colz");
 // c1_canvas->SaveAs("CompareTrees_PUid/Compare_"+branch_name+"__"+gl_label_old+"_vs_"+gl_label_new+".pdf");
 c1_canvas->SaveAs("CompareTrees_DiJet_JER_L2Res/Compare_"+branch_name+"__"+gl_label_old+"_vs_"+gl_label_new+".pdf");

}
