// #include "../../../include/constants.h"
#include "TSystem.h"

using namespace std;

TString study = "-STUDY-";
TString label = "-LABEL-";
TString corr = "-CORR-";
TString qcd = "-SAMPLE-";


bool debug = false;

void fillHists(TFile* f);
// void FillAllHists(vector<TString> window, vector<vector<double>> cutvalues, TString var);
// TH1F* get_hist(TFile* file, vector<double> cut, TString weightname, TString var, bool weight_cut=false);
inline double Round(double wert, int nachkommastellen) {return ((int) ((wert * pow(10.0, nachkommastellen)) + ((wert < 0)? - 0.5 : 0.5))) / pow(10.0, nachkommastellen);}

map<TString, TFile*> files;
//  [eta]        [ptgen]      [hist]
map<TString, map<TString, map<TString, TH1F*>>> hists;

map<TString, map<TString, double>> eta_bins = {
  {  "1", { {"min", 0.000}, {"max", 0.261} } },
  {  "2", { {"min", 0.261}, {"max", 0.522} } },
  {  "3", { {"min", 0.522}, {"max", 0.783} } },
  {  "4", { {"min", 0.783}, {"max", 1.044} } },
  {  "5", { {"min", 1.044}, {"max", 1.305} } },
  {  "6", { {"min", 1.305}, {"max", 1.566} } },
  {  "7", { {"min", 1.566}, {"max", 1.740} } },
  {  "8", { {"min", 1.740}, {"max", 1.930} } },
  {  "9", { {"min", 1.930}, {"max", 2.043} } },
  { "10", { {"min", 2.043}, {"max", 2.172} } },
  { "11", { {"min", 2.172}, {"max", 2.322} } },
  { "12", { {"min", 2.322}, {"max", 2.500} } },
  { "13", { {"min", 2.500}, {"max", 2.650} } },
  { "14", { {"min", 2.650}, {"max", 2.853} } },
  { "15", { {"min", 2.853}, {"max", 2.964} } },
  { "16", { {"min", 2.964}, {"max", 3.139} } },
  { "17", { {"min", 3.139}, {"max", 3.489} } },
  { "18", { {"min", 3.489}, {"max", 3.839} } },
  { "19", { {"min", 3.839}, {"max", 5.191} } },
};

map<TString, map<TString, double>> pt_bins_forward = {
  {  "1", { {"min",   0}, {"max",   93} } },
  {  "2", { {"min",  93}, {"max",   99} } },
  {  "3", { {"min",  99}, {"max",  106} } },
  {  "4", { {"min", 106}, {"max",  116} } },
  {  "5", { {"min", 116}, {"max",  122} } },
  {  "6", { {"min", 122}, {"max",  130} } },
  {  "7", { {"min", 130}, {"max",  142} } },
  {  "8", { {"min", 142}, {"max",  154} } },
  {  "9", { {"min", 154}, {"max",  172} } },
  { "10", { {"min", 172}, {"max",  210} } },
  { "11", { {"min", 210}, {"max",  220} } },
  { "12", { {"min", 220}, {"max",  240} } },
  { "13", { {"min", 240}, {"max",  279} } },
  { "14", { {"min", 279}, {"max",  379} } },
  { "15", { {"min", 379}, {"max", 1500} } },
};

map<TString, map<TString, double>> pt_bins_central = {
  {  "1", { {"min",    0}, {"max",   66} } },
  {  "2", { {"min",   66}, {"max",   71} } },
  {  "3", { {"min",   71}, {"max",   77} } },
  {  "4", { {"min",   77}, {"max",   93} } },
  {  "5", { {"min",   93}, {"max",   98} } },
  {  "6", { {"min",   98}, {"max",  106} } },
  {  "7", { {"min",  106}, {"max",  118} } },
  {  "8", { {"min",  118}, {"max",  128} } },
  {  "9", { {"min",  128}, {"max",  145} } },
  { "10", { {"min",  145}, {"max",  189} } },
  { "11", { {"min",  189}, {"max",  203} } },
  { "12", { {"min",  203}, {"max",  223} } },
  { "13", { {"min",  223}, {"max",  257} } },
  { "14", { {"min",  257}, {"max",  291} } },
  { "15", { {"min",  291}, {"max",  325} } },
  { "16", { {"min",  325}, {"max",  358} } },
  { "17", { {"min",  358}, {"max",  391} } },
  { "18", { {"min",  391}, {"max",  434} } },
  { "19", { {"min",  434}, {"max",  478} } },
  { "20", { {"min",  478}, {"max",  531} } },
  { "21", { {"min",  531}, {"max",  546} } },
  { "22", { {"min",  546}, {"max",  563} } },
  { "23", { {"min",  563}, {"max",  585} } },
  { "24", { {"min",  585}, {"max",  644} } },
  { "25", { {"min",  644}, {"max",  730} } },
  { "26", { {"min",  730}, {"max",  790} } },
  { "27", { {"min",  790}, {"max",  840} } },
  { "28", { {"min",  840}, {"max",  920} } },
  { "29", { {"min",  920}, {"max", 1020} } },
  { "30", { {"min", 1020}, {"max", 1120} } },
  { "31", { {"min", 1120}, {"max", 1500} } },
};

vector<double> vec_eta_bins = { 0.000,  0.261,  0.522,  0.783,  1.044,  1.305,  1.566,  1.740,  1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,  2.964,  3.139,  3.489,  3.839,  5.191 };

vector<double> ptbins_central = {0, 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120 };
vector<double> ptbins_forward = {0, 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 };




int Response_MCTruth_Parallel(){

  cout << "-------------------------------------------" << endl;
  cout << "Usage: ./Response_MCTruth" << endl;
  cout << "-------------------------------------------" << endl;

  TH1::AddDirectory(kFALSE);

  TFile *file = new TFile("/nfs/dust/cms/user/paaschal/sframe_all/DiJetJERC_DiJetHLT/UL18/"+study+"/"+corr+"/"+label+"/uhh2.AnalysisModuleRunner.MC."+qcd+"_UL18.root");

  for(unsigned int e=0; e<vec_eta_bins.size(); e++){
    vector<double> ptbins = vec_eta_bins[e+1]<2.853?ptbins_central:ptbins_forward;
    // if(vec_eta_bins[e]<2.853) cout << e << " " << endl;
    for(unsigned int p=0; p<ptbins.size(); p++){
      TString ebin = to_string(e+1);
      TString pbin = to_string(p+1);
      // cout << ebin << " " << pbin << endl;
      hists["response"][ebin][pbin] = new TH1F("response_eta"+ebin+"_pt"+pbin, "response_eta"+ebin+"_pt"+pbin, 200, 0, 2.0);
      hists["response"][ebin][pbin]->Sumw2();
      hists["ptgen"][ebin][pbin] = new TH1F("ptgen_eta"+ebin+"_pt"+pbin, "ptgen_eta"+ebin+"_pt"+pbin, 150, 0, 1500);
      hists["ptgen"][ebin][pbin]->Sumw2();
    }
  }

  cout << "Start processing " << qcd << " ... " << endl;
  auto start = chrono::high_resolution_clock::now(); // Calculation time - start
  fillHists(file);
  auto stop = chrono::high_resolution_clock::now();  // Calculation time - stop
  auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
  cout << "\t ... End processing after " << duration.count() << "s" << endl;

  // study.ReplaceAll("eta_", "");
  TString outpath = std::getenv("CMSSW_BASE");
  outpath += "/src/UHH2/DiJetJERC/JERSF_Analysis/MCTruth/wide_eta_bin/files/parallel/";
  outpath += study+"/";
  outpath += corr+"/";
  outpath += label+"/";
  TFile *fout = new TFile(outpath+"/"+qcd+".root", "recreate");
  fout->cd();
  for(auto h: hists){
    for(auto e: hists[h.first]){
      for(auto p: hists[h.first][e.first]){
        hists[h.first][e.first][p.first]->Write(hists[h.first][e.first][p.first]->GetTitle());
      }
    }
  }

  return 1;
}

// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------

void fillHists(TFile* f){

  TTree* tree = (TTree *) f->Get("AnalysisTree");

  Float_t barrelgenjet_pt, barrelgenjet_eta, barrelgenjet_phi;
  Float_t barreljet_pt, barreljet_eta, barreljet_phi;
  Float_t probegenjet_pt, probegenjet_eta, probegenjet_phi;
  Float_t probejet_pt, probejet_eta, probejet_phi;
  Float_t weight;
  Int_t njet, ngenjet;

  double gen_thr = 10;
  double jet_thr = 15;
  double s_delta_phi = 2.7;
  double s_delta_R = 0.3;

  TString ietaBarrel = "0"; TString iptBarrel = "0";
  TString ietaProbe = "0"; TString iptProbe = "0";

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("barrelgenjet_pt", &barrelgenjet_pt);
  tree->SetBranchAddress("barrelgenjet_eta", &barrelgenjet_eta);
  tree->SetBranchAddress("barrelgenjet_phi", &barrelgenjet_phi);
  tree->SetBranchAddress("barreljet_pt", &barreljet_pt);
  tree->SetBranchAddress("barreljet_eta", &barreljet_eta);
  tree->SetBranchAddress("barreljet_phi", &barreljet_phi);
  tree->SetBranchAddress("probegenjet_pt", &probegenjet_pt);
  tree->SetBranchAddress("probegenjet_eta", &probegenjet_eta);
  tree->SetBranchAddress("probegenjet_phi", &probegenjet_phi);
  tree->SetBranchAddress("probejet_pt", &probejet_pt);
  tree->SetBranchAddress("probejet_eta", &probejet_eta);
  tree->SetBranchAddress("probejet_phi", &probejet_phi);
  tree->SetBranchAddress("weight", &weight);

  tree->SetBranchAddress("Njet", &njet);
  tree->SetBranchAddress("Ngenjet", &ngenjet);

  tree->SetBranchStatus("*",1);

  int Events = 0;
  double TotalEvents = tree->GetEntriesFast();
  cout << "\t ... Total Entries " << TotalEvents << endl;
  auto start = chrono::high_resolution_clock::now(); // Calculation time - start

  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(debug) cout << "\n +++ NEW EVENT +++ \n" << endl;
    Events++;
    if ( Events%1000000 == 0 ) {
      auto stop  = chrono::high_resolution_clock::now();  // Calculation time - stop
      auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
      auto time = duration.count();
      double left = Round( (TotalEvents/Events-1)*time / 60, 1);
      std::cout << "\t\t ... Analyzing event #" << Events << " (" << time << "s) - " << left << " min left" << std::endl;
    }
    if(tree->GetEntry(ievent)<=0) break;

    // -----------------------------------------------------------------------
    // --- Selection

    if(njet<2) continue;
    if(ngenjet<2) continue;

    if(barreljet_pt<jet_thr && barrelgenjet_pt<gen_thr) continue;
    if ( TMath::Abs(TVector2::Phi_mpi_pi((probejet_phi - barreljet_phi))) < s_delta_phi-0.1 ) { /*std::cout << "Jets are not back to back: " << TMath::Abs(TVector2::Phi_mpi_pi((probejet_phi - barreljet_phi))) << std::endl;*/ continue;}
    if ( TMath::Abs(TVector2::Phi_mpi_pi((barrelgenjet_phi - probegenjet_phi))) < s_delta_phi-0.1 ) { /*std::cout << "GenJets are not back to back: " << TMath::Abs(TVector2::Phi_mpi_pi((barrelgenjet_phi - probegenjet_phi))) << " " << event << " " << run << std::endl;*/ continue;}


    // -----------------------------------------------------------------------
    // --- Collect eta and pt bin

    TString ebin = "EMPTY";
    TString pbin = "EMPTY";

    Float_t barrel_eta = fabs(barreljet_eta);
    Float_t barrel_pt = fabs(barrelgenjet_pt);
    Float_t probe_eta = fabs(probejet_eta);
    Float_t probe_pt = fabs(probegenjet_pt);
    bool probeassigned = false;
    bool barrelassigned = false;

    for(auto e: eta_bins){
      if(debug) cout << "In eta bin " << e.first << " " << eta_bins[e.first]["min"] << " " << eta_bins[e.first]["max"] << " | " << barreljet_eta << " " << probejet_eta << endl;
      if(eta_bins[e.first]["min"]<barrel_eta && barrel_eta<eta_bins[e.first]["max"] && !barrelassigned){
        ebin = e.first;
        map<TString, map<TString, double>> ptbins = vec_eta_bins[e.first.Atoi()]<2.853?pt_bins_central:pt_bins_forward;
        for(auto p: ptbins){
          if(barrelassigned) continue;
          bool greater = ptbins[p.first]["min"]<barrel_pt;
          bool smaller = barrel_pt<ptbins[p.first]["max"];
          if(debug) cout << "\t\t ... PT BARREL " << p.first << " " << ptbins[p.first]["min"] << " " << barrel_pt << " " << ptbins[p.first]["max"] <<endl;
          if(ptbins[p.first]["min"]<barrel_pt && barrel_pt<ptbins[p.first]["max"]){
            pbin = p.first;
            ietaBarrel = ebin; iptBarrel = pbin;
            barrelassigned = true;
          }
        }
      }

      if(eta_bins[e.first]["min"]<probe_eta && probe_eta<eta_bins[e.first]["max"] && !probeassigned){
        ebin = e.first;
        map<TString, map<TString, double>> ptbins = vec_eta_bins[e.first.Atoi()]<2.853?pt_bins_central:pt_bins_forward;
        if(debug) cout << "Size " << ptbins.size() << endl;
        for(auto p: ptbins){
          if(debug) cout << "\t\t ... PT Probe " << p.first << " " << ptbins[p.first]["min"] << " " << probe_pt << " " << ptbins[p.first]["max"] <<endl;
          if(ptbins[p.first]["min"]<probe_pt && probe_pt<ptbins[p.first]["max"]){
            if(probeassigned) continue;
            pbin = p.first;
            ietaProbe = ebin; iptProbe = pbin;
            probeassigned = true;
          }
        }
      }

      if(barrelassigned && probeassigned){
        if(debug) cout << "Passed Event " << Events << endl;
      }
    }

    // -----------------------------------------------------------------------
    // --- Fill Hists

    double DPhiP = TMath::Abs(TVector2::Phi_mpi_pi(probejet_phi - probegenjet_phi));
    double DPhiB = TMath::Abs(TVector2::Phi_mpi_pi(barreljet_phi- barrelgenjet_phi));
    double DEtaP = TMath::Abs(probejet_eta - probegenjet_eta);
    double DEtaB = TMath::Abs(barreljet_eta- barrelgenjet_eta);
    double DRP = TMath::Sqrt( TMath::Power(DPhiB, 2 ) + TMath::Power(DEtaB, 2));
    double DRB = TMath::Sqrt( TMath::Power(DPhiP, 2 ) + TMath::Power(DEtaP, 2));

    double response_barrel = barreljet_pt/barrelgenjet_pt;
    double response_probe = probejet_pt/probegenjet_pt;

    if(debug) cout << ietaBarrel << " " << iptBarrel << " " << ietaProbe << " " << iptProbe << endl;
    if(DRB < s_delta_R){
      hists["response"][ietaBarrel][iptBarrel]->Fill(response_barrel, weight);
      hists["ptgen"][ietaBarrel][iptBarrel]->Fill(barrelgenjet_pt, weight);
    }
    if(DRP < s_delta_R){
      hists["response"][ietaProbe][iptProbe]->Fill(response_probe, weight);
      hists["ptgen"][ietaProbe][iptProbe]->Fill(probegenjet_pt, weight);
    }


  }
}
