#include "UHH2/BaconJets/include/JECAnalysisHists.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"

// #include "UHH2/BaconTrans/baconheaders/TJet.hh"
// #include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"
// #include "UHH2/BaconTrans/baconheaders/BaconAnaDefs.hh"
// #include "UHH2/BaconTrans/baconheaders/TVertex.hh"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <getopt.h>
using namespace std;
using namespace uhh2;
//using namespace baconhep;
//    uhh2::Event::Handle<TClonesArray> h_pv;
JECAnalysisHists::JECAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    book<TH1F>("N_jets", "N_{jets}", 50, -0.5, 49.5);
    book<TH1F>("pt_hat", "p_{T} hat", 150, 0, 6000);
    book<TH1F>("pt","p_{T} all jets; p_{T} (GeV)",100,0,1500);
    book<TH1F>("pt_modT1METjets","p_{T} all jets in 2.650<#eta<3.139; p_{T} (GeV)",100,0,1500);
    book<TH1F>("eta","#eta all jets; #eta",100,-5,5);
    double eta_bins[n_eta];
    for(int i=0; i<n_eta; i++) eta_bins[i] = eta_range[i];
    book<TH1F>("eta_binned","|#eta| all jets; |#eta|",n_eta-1, eta_bins);
    book<TH1F>("phi","#phi all jets; #phi",50,-M_PI,M_PI);
    book<TH1F>("MET","PF MET; PF MET",1000,0,1000);
    book<TH1F>("rawchsMET","raw CHS MET; raw CHS MET",1000,0,1000);

    book<TH1F>("nPu","Number of PU events",60,0,60);
    book<TH1F>("matchJetId_0","Matching Result for first TrgObj Jet",13,-3,10);
    book<TH1F>("matchJetId_1","Matching Result for second TrgObj Jet",13,-3,10);  
    book<TH1F>("N_PV","Number of PVtx",60,0,60);
    book<TH1F>("weight_histo","weight_histo ",20,0,2);

    book<TH1F>("pt_1","p_{T} jet 1",100,0,600);
    book<TH1F>("eta_1","#eta jet 1",100,-5,5);

    book<TH1F>("pt_2","p_{T} jet 2",100,0,600);
    book<TH1F>("eta_2","#eta jet 2",100,-5,5);

    book<TH1F>("pt_3","p_{T} jet 3",100,0,600);
    book<TH1F>("eta_3","#eta jet 3",100,-5,5);
    
    book<TH1F>("pt_1_onoff_Resp","p_{T,offline} / p_{T,online} jet 1",100,0,2);
    book<TH1F>("pt_2_onoff_Resp","p_{T,offline} / p_{T,online} jet 2",100,0,2);
     
    book<TH1F>("ptRaw_barrel","p^{Raw}_{T} barrel jet; p_{T}^{Raw,barrel} (GeV)",100,0,600);
    book<TH1F>("ptRaw_probe","p^{Raw}_{T} probe jet; p_{T}^{Raw,probe} (GeV)",100,0,600);
    book<TH1F>("pt_barrel","p_{T} barrel jet; p_{T}^{barrel} (GeV)",100,0,600);
    book<TH1F>("pt_probe","p_{T} probe jet; p_{T}^{probe} (GeV)",100,0,600);
    book<TH1F>("eta_barrel","#eta barrel jet; #eta^{barrel}",100,-5,5);
    book<TH1F>("eta_probe","#eta probe jet; #eta^{probe}",100,-5,5);
    book<TH1F>("eta_probe_pos","#eta probe jet >=0; #eta^{probe}",50,0,5);
    book<TH1F>("eta_probe_neg","#eta probe jet <0; |#eta^{probe}|",50,0,5);
    book<TH1F>("pt_ave","p_{T} ave jet; p_{T}^{ave} (GeV)",600,0,600);
    book<TH1F>("pt_ave_pthat","p_{T} ave jet; p_{T}^{ave} - p_{T}^{hat})/p_{T}^{hat}(GeV)",100,0,600);
    book<TH1F>("pt_ave_rebin","p_{T} ave jet; p_{T}^{ave} (GeV)",300,0,3000);

    book<TH1F>("pt_rel","p_{T}^{jet3} / p_{T}^{ave}; #alpha ", 50, 0, 1);

    book<TH1F>("asym","asymmetry jet 1 and jet 2; Asymmetry",150,-1.5,1.5);
    book<TH1F>("mpf","MPF response; MPF response",250,0.,2.5);
    book<TH1F>("r_rel","R_{rel}; R_{rel}; Relative response",250,0.,2.5);
 

    book<TH2D>("mpf_vs_etaProbe","MPF response vs. #eta probe jet; #eta probe; MPF response",100,-5,5,100,0.,2.);
    book<TH2D>("r_rel_vs_etaProbe","Relative response vs. #eta probe jet; #eta probe; R_{rel}",100,-5,5,100,0.,2.);
    book<TH2D>("pt_ave_vs_etaProbe","pt ave vs #eta probe jet; #eta probe; pT_{ave}, GeV",100,-5.2,5.2,200,0,1000);
    book<TH2D>("dPhi_vs_alpha","#Delta #Phi vs #alpha; #alpha; #Delta #Phi",150,0,3,200,-M_PI,4);

    book<TH2D>("phi_vs_etaProbe","#phi vs. #eta probe jet; #eta probe; #phi probe",100,-5,5,100,-M_PI,M_PI);
    book<TH2D>("phi_vs_etaAll","#phi vs. #eta all jet; #eta; #phi",100,-5,5,100,-M_PI,M_PI);

    book<TH2D>("phi_vs_etaAllCut30","#phi vs. #eta all jet, cut at Pt > 30GeV; #eta; #phi",100,-6,6,100,-M_PI,M_PI);   

    tt_gen_pthat  = ctx.get_handle<float>("gen_pthat");
    tt_gen_weight = ctx.get_handle<float>("gen_weight");
    tt_jet1_pt = ctx.get_handle<float>("jet1_pt");
    tt_jet2_pt = ctx.get_handle<float>("jet2_pt");
    tt_jet3_pt = ctx.get_handle<float>("jet3_pt");
    tt_jet1_ptRaw = ctx.get_handle<float>("jet1_ptRaw");
    tt_jet2_ptRaw = ctx.get_handle<float>("jet2_ptRaw");
    tt_jet3_ptRaw = ctx.get_handle<float>("jet3_ptRaw");
    tt_jet1_pt_onoff_Resp = ctx.get_handle<float>("jet1_pt_onoff_Resp");
    tt_jet2_pt_onoff_Resp = ctx.get_handle<float>("jet2_pt_onoff_Resp");    
    tt_nvertices = ctx.get_handle<int>("nvertices");
    tt_probejet_eta = ctx.get_handle<float>("probejet_eta");
    tt_probejet_phi = ctx.get_handle<float>("probejet_phi");
    tt_probejet_pt = ctx.get_handle<float>("probejet_pt");
    tt_probejet_ptRaw = ctx.get_handle<float>("probejet_ptRaw");
    tt_barreljet_eta = ctx.get_handle<float>("barreljet_eta");
    tt_barreljet_phi = ctx.get_handle<float>("barreljet_phi");
    tt_barreljet_pt = ctx.get_handle<float>("barreljet_pt");
    tt_barreljet_ptRaw = ctx.get_handle<float>("barreljet_ptRaw");
    tt_pt_ave = ctx.get_handle<float>("pt_ave");
    tt_alpha = ctx.get_handle<float>("alpha");
    tt_alpha_sum = ctx.get_handle<float>("alpha_sum");
    tt_rel_r = ctx.get_handle<float>("rel_r");
    tt_mpf_r = ctx.get_handle<float>("mpf_r");
    tt_asymmetry = ctx.get_handle<float>("asymmetry");
    tt_nPU = ctx.get_handle<int>("nPU");
    tt_matchJetId_0 =ctx.get_handle<int>("matchJetId_0");
    tt_matchJetId_1 =ctx.get_handle<int>("matchJetId_1");    
}

void JECAnalysisHists::fill(const uhh2::Event & ev){
  fill(ev, 0);
}
void JECAnalysisHists::fill(const uhh2::Event & ev, const int rand){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  // Don't forget to always use the weight when filling.


  double weight = ev.weight;
  const int njets = ev.jets->size();
  
  hist("N_jets")->Fill(njets, weight);

  try{
  if(!ev.isRealData)hist("pt_hat")->Fill(ev.get(tt_gen_pthat),weight);
    }
  catch(const std::runtime_error& error){
	  std::cout<<"Problem with gen in JECAnalysishists.cxx"<<std::endl;
	  std::cout<<error.what();
  }

  for (int i=0; i<njets; i++){
    Jet* jets = &ev.jets->at(i);
    hist("pt")->Fill(jets->pt(), weight);
    if(jets->eta()>2.650&&jets->eta()<3.139)
      hist("pt_modT1METjets")->Fill(jets->pt(), weight);    
    hist("eta")->Fill(jets->eta(), weight);
    hist("eta_binned")->Fill(jets->eta(), weight);
    hist("phi")->Fill(jets->phi(), weight);
    
    hist("pt_1_onoff_Resp")->Fill(ev.get(tt_jet1_pt_onoff_Resp), weight);
    hist("pt_2_onoff_Resp")->Fill(ev.get(tt_jet2_pt_onoff_Resp), weight);   
    
    hist("MET")->Fill(ev.met->pt(), weight);
    hist("rawchsMET")->Fill(hypot(ev.met->rawCHS_px(),ev.met->rawCHS_py()), weight);
    hist("nPu")->Fill(ev.get(tt_nPU), weight);
    hist("matchJetId_0")->Fill(ev.get(tt_matchJetId_0), weight);
    hist("matchJetId_1")->Fill(ev.get(tt_matchJetId_1), weight);			       
    hist("weight_histo")->Fill(weight, 1);
    
    ((TH2D*)hist("phi_vs_etaAll"))->Fill(jets->eta(),jets->phi(),weight);

    if(jets->pt()>30) ((TH2D*)hist("phi_vs_etaAllCut30"))->Fill(jets->eta(),jets->phi(),weight);

  }
  
  hist("N_PV")->Fill(ev.get(tt_nvertices), weight);
    
  

  Jet* jet1 = &ev.jets->at(0);
  hist("pt_1")->Fill(jet1->pt(), weight);
  hist("eta_1")->Fill(jet1->eta(), weight);
  Jet* jet2 = &ev.jets->at(1);
  hist("pt_2")->Fill(jet2->pt(), weight);
  hist("eta_2")->Fill(jet2->eta(), weight);
  // float ratio_pt = 0.;

  try{ 
     float ratio_pt = (ev.get(tt_pt_ave) - ev.get(tt_gen_pthat))/ev.get(tt_gen_pthat);
  }
  catch(const std::runtime_error& error){
	  std::cout<<"Problem with gen in JECAnalysishists.cxx"<<std::endl;
	  std::cout<<error.what();
  }  

  hist("pt_ave")          ->Fill(ev.get(tt_pt_ave), weight);
  hist("pt_ave_pthat")   ->Fill(ev.get(tt_pt_ave), weight); //FIXME why is this filled with pt_ave?
  hist("pt_ave_rebin") ->Fill(ev.get(tt_pt_ave), weight);  //FIXME why is this filled with pt_ave?
  hist("ptRaw_barrel")    ->Fill(ev.get(tt_barreljet_ptRaw), weight);
  hist("ptRaw_probe")     ->Fill(ev.get(tt_probejet_ptRaw) , weight);
  hist("pt_barrel")   ->Fill(ev.get(tt_barreljet_pt), weight);
  hist("pt_probe")    ->Fill(ev.get(tt_probejet_pt) , weight);
  hist("eta_barrel")  ->Fill(ev.get(tt_barreljet_eta), weight);
  hist("eta_probe")   ->Fill(ev.get(tt_probejet_eta) , weight);
  if(ev.get(tt_probejet_eta)>=0) hist("eta_probe_pos")->Fill(ev.get(tt_probejet_eta) , weight);
  else hist("eta_probe_neg")->Fill(fabs(ev.get(tt_probejet_eta)) , weight);
  hist("mpf")         ->Fill(ev.get(tt_mpf_r), weight);
  hist("asym")        ->Fill(ev.get(tt_asymmetry), weight);
  hist("r_rel")       ->Fill(ev.get(tt_rel_r), weight);

  ((TH2D*)hist("mpf_vs_etaProbe"))->Fill(ev.get(tt_probejet_eta),ev.get(tt_mpf_r),weight);
  ((TH2D*)hist("r_rel_vs_etaProbe"))->Fill(ev.get(tt_probejet_eta),ev.get(tt_rel_r),weight);
  float pt_ave = (0.5*(ev.get(tt_jet1_pt) + ev.get(tt_jet2_pt)));
  ((TH2D*)hist("pt_ave_vs_etaProbe"))->Fill(ev.get(tt_probejet_eta),pt_ave,weight);

  ((TH2D*)hist("phi_vs_etaProbe"))->Fill(ev.get(tt_probejet_eta),ev.get(tt_probejet_phi),weight);

  double deltaPhi = std::abs(TVector2::Phi_mpi_pi(jet1->phi() - jet2->phi()));
  ((TH2D*)hist("dPhi_vs_alpha"))->Fill(ev.get(tt_alpha),deltaPhi, weight);

  if (njets > 2){
    Jet* jet3 = &ev.jets->at(2);
    hist("pt_3")->Fill(ev.get(tt_jet3_pt), weight);
    hist("eta_3")->Fill(jet3->eta(), weight);
    hist("pt_rel")->Fill(ev.get(tt_jet3_pt)/(0.5*(ev.get(tt_barreljet_pt) + ev.get(tt_probejet_pt) )),weight);

  }

}

JECAnalysisHists::~JECAnalysisHists(){}
