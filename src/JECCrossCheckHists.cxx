#include "UHH2/DiJetJERC/include/JECCrossCheckHists.h"
#include "UHH2/DiJetJERC/include/constants.h"
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
JECCrossCheckHists::JECCrossCheckHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    book<TH1F>("sumweights", "Sum of event weights", 1, 0.5, 1.5);
    book<TH1F>("N_jets", "N_{jets}", 50, -0.5, 49.5);
    book<TH1F>("pt_hat", "p_{T} hat", 150, 0, 6000);
    book<TH1F>("PU_pt_hat","PU p_{T} hat", 150,0,6000);
    book<TH1F>("PU_pt_hat_Ratio","PU p_{T} hat  Ratio", 50,0,5);
    PU_vs_pt_hat = book<TH2F>("PU_pt_hat_vs_pt_hat","x=p_{T}hat y=PUp_{T}hat",150,0,6000 ,150,0,6000);

    book<TH1F>("pt","p_{T} all jets; p_{T} (GeV)",100,0,1500);
    book<TH1F>("pt_modT1METjets","p_{T} all jets in 2.650<#eta<3.139; p_{T} (GeV)",100,0,1500);
    book<TH1F>("eta","#eta all jets; #eta",100,-5,5);
    book<TH1F>("eta_binned","|#eta| all jets; |#eta|",n_eta-1, eta_bins);
    book<TH1F>("phi","#phi all jets; #phi",50,-M_PI,M_PI);

    book<TH1F>("MET","PF MET; PF MET",1000,0,1000);
    book<TH1F>("rawchsMET","raw CHS MET; raw CHS MET",1000,0,1000);

    book<TH1F>("nPu","Number of PU events",60,0,60);
    book<TH1F>("N_PV","Number of PVtx",60,0,60);
    book<TH1F>("weight_histo","weight_histo ",20,0,2);
    Weight_vs_pt_hat = book<TH2F>("Weight_vs_pt_hat","x=p_{T}hat y=Weight",150,0,3000 ,150,0,3000);
    Weight_vs_pt_hat_zoom = book<TH2F>("Weight_vs_pt_hat_zoom","x=p_{T}hat y=Weight",100,0,500 ,100,0,500);

    book<TH1F>("pt_1","p_{T} jet 1",100,0,600);
    book<TH1F>("eta_1","#eta jet 1",100,-5,5);

    book<TH1F>("pt_2","p_{T} jet 2",100,0,600);
    book<TH1F>("eta_2","#eta jet 2",100,-5,5);

    book<TH1F>("pt_3","p_{T} jet 3",100,0,600);
    book<TH1F>("eta_3","#eta jet 3",100,-5,5);
    book<TH1F>("pt_ave","p_{T} ave jet; p_{T}^{ave} (GeV)",600,0,600);
    book<TH1F>("pt_ave_rebin","p_{T} ave jet; p_{T}^{ave} (GeV)",300,0,3000);

    book<TH1F>("asym","asymmetry jet 1 and jet 2; Asymmetry",150,-1.5,1.5);

    //FIXME: commented in for madgraph samples because the genInfo seems to be missing!
       no_genp = true;
       if(no_genp) cout<<"!!! WARNING, no genparticle are used! !!!"<<endl;
}

void JECCrossCheckHists::fill(const uhh2::Event & ev){
  fill(ev, 0);
}
void JECCrossCheckHists::fill(const uhh2::Event & ev, const int rand){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  // Don't forget to always use the weight when filling.

  double weight = ev.weight;
  const int njets = ev.jets->size();

  hist("sumweights")->Fill(1, weight);
  hist("N_jets")->Fill(njets, weight);
   // if(!ev.isRealData && !no_genp){
  if(!ev.isRealData){
    // double pt_hat = ev.genInfo->binningValues()[0];
    double pt_hat = ev.genInfo->qScale();
    double PU_pt_hat = ev.genInfo->PU_pT_hat_max();
    hist("pt_hat")->Fill(pt_hat,weight);
    hist("PU_pt_hat")->Fill(PU_pt_hat, weight);
    hist("PU_pt_hat_Ratio")->Fill(PU_pt_hat/pt_hat, weight);
    PU_vs_pt_hat->Fill(pt_hat, PU_pt_hat,weight);
    Weight_vs_pt_hat->Fill(pt_hat, weight, weight);
  }

  for (int i=0; i<njets; i++){
    Jet* jets = &ev.jets->at(i);
    hist("pt")->Fill(jets->pt(), weight);
    if(jets->eta()>2.650&&jets->eta()<3.139)
      hist("pt_modT1METjets")->Fill(jets->pt(), weight);
    hist("eta")->Fill(jets->eta(), weight);
    hist("eta_binned")->Fill(jets->eta(), weight);
    hist("phi")->Fill(jets->phi(), weight);
  }
  hist("MET")->Fill(ev.met->pt(), weight);
  hist("rawchsMET")->Fill(hypot(ev.met->rawCHS_px(),ev.met->rawCHS_py()), weight);
  double nPU = 0;
   if(!ev.isRealData && !no_genp) nPU = ev.genInfo->pileup_TrueNumInteractions();
  // if(!ev.isRealData) nPU = ev.genInfo->pileup_TrueNumInteractions();
  hist("nPu")->Fill(nPU, weight);
  hist("weight_histo")->Fill(weight, 1);


  hist("N_PV")->Fill(ev.pvs->size(), weight);


  if(njets > 0){
    Jet* jet1 = &ev.jets->at(0);
    hist("pt_1")->Fill(jet1->pt(), weight);
    hist("eta_1")->Fill(jet1->eta(), weight);
  }
  if(njets >1){
    Jet* jet1 = &ev.jets->at(0);
    Jet* jet2 = &ev.jets->at(1);
    hist("pt_2")->Fill(jet2->pt(), weight);
    hist("eta_2")->Fill(jet2->eta(), weight);
    double pt_ave = (jet1->pt()+jet2->pt())/2;
    hist("pt_ave")       ->Fill(pt_ave, weight);
    hist("pt_ave_rebin") ->Fill(pt_ave, weight);
    double A = (jet2->pt() - jet1->pt()) / (jet2->pt() + jet1->pt());
    hist("asym")         ->Fill(A, weight);
  }
  if (njets > 2){
    Jet* jet3 = &ev.jets->at(2);
    hist("pt_3")->Fill(jet3->pt(), weight);
    hist("eta_3")->Fill(jet3->eta(), weight);

  }
}

JECCrossCheckHists::~JECCrossCheckHists(){}
