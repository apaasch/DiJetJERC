// investigation of difference between output trees after ntuplewriter or UHH2 AnalysisModule
// implemented for RelVals checks,  but might be useful in general
// At the moment works only for variables in single branches, i.e does not work for variables from object collections/vectors
// author: A.Karavdina
// date: 15.03.2019
// Run it with following command:
// root -l -b -q CompareTrees_RelVals.C\(\"Njet\",\"Njet\",false,50,0,50,10\)


void CompareTrees_RelVals(TString var_name="Ngenjet",TString branch_name="Ngenjet", bool isFloat=false, int nbins=10, double minx=0,double maxx=10, int maxevents=100){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);

 double w = 600;
 double h = 600;

 // bool isWeights = true;
 bool isWeights = false;

 //Files after UHH2
 TString path = "/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1__RelVals_104/AK4CHS_minJet10/WithJEC_2017V32_2018V7_2016V11/";
 TString name_old = "uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_refHCAL_AK4CHS.root";
 TString name_new = "uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_newHCAL_AK4CHS.root";


 // TString path = "/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1__RelVals_105rel2/AK4CHS_minJet10/WithJEC_2017V32_2018V7_2016V11/";
 // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_refHCAL_AK4CHS.root";
 // // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_newHCAL_AK4CHS.root";
 // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_newECALHCAL_AK4CHS.root";
 // //TEST ---
 // TString path = "/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1__RelVals_105rel2/AK4CHS/WithJEC_2017V32_2018V7_2016V11/workdir_DiJet_RunA/";
 // TString name_old = "uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_refHCAL_AK4CHS_9.root";
 // TString name_new = "uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_newHCAL_AK4CHS_9.root";

 TString gl_label_old = "refHCAL";
 TString gl_label_new = "newHCAL";

 TFile *input_old = TFile::Open(path+name_old);
 TTree *TTreeAna_old = (TTree*)input_old->Get("AnalysisTree");
 TFile *input_new = TFile::Open(path+name_new);
 TTree *TTreeAna_new = (TTree*)input_new->Get("AnalysisTree");

 int eventID,feventID;
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
   if (feventID==eventID) {
     nok++;
     if(isFloat){
       hist_comp->Fill(var_val_f,fval_val_f);
     }
     else{
       hist_comp->Fill(var_val_i,fval_val_i);
     }
   }
 }
 printf("nok = %d, fentries=%lld\n",nok,TTreeAna_new->GetEntries());

 TCanvas * c1_canvas = new TCanvas("canvas", "c", w, h);
 hist_comp->Draw("colz");
 c1_canvas->SaveAs("CompareTrees_plots_104_minJet10/Compare_"+branch_name+"__"+gl_label_old+"_vs_"+gl_label_new+".pdf");

}
