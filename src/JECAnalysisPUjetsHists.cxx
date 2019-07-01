#include "UHH2/DiJetJERC/include/JECAnalysisPUjetsHists.h"
#include "UHH2/DiJetJERC/include/constants.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/common/include/Utils.h"

// #include "UHH2/BaconTrans/baconheaders/TJet.hh"
// #include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"
// #include "UHH2/BaconTrans/baconheaders/BaconAnaDefs.hh"
// #include "UHH2/BaconTrans/baconheaders/TVertex.hh"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <getopt.h>
using namespace std;
using namespace uhh2;
//using namespace baconhep;
//    uhh2::Event::Handle<TClonesArray> h_pv;
JECAnalysisPUjetsHists::JECAnalysisPUjetsHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    double eta_bins[n_eta_RelVals];
    for(int i=0; i<n_eta_RelVals; i++) eta_bins[i] = eta_range_RelVals[i];
    book<TH2D>("Jet_eta_pt", "Jet_eta_pt;RECO jet #eta; jet pt", n_eta_RelVals-1, eta_bins, 100,0,1000);

    book<TH1F>("nPu","Number of PU events",100,0,100);
    book<TH1F>("N_PV","Number of PVtx",100,0,100);

    tt_nPU = ctx.get_handle<int>("nPU");
    tt_nvertices = ctx.get_handle<int>("nvertices");

    book<TH1F>("PU_id","PU id discrimintant",120,-1.1,1.1);
    book<TH2D>("Jet_eta_PUid", "Jet_eta_PUid;RECO jet #eta; PU id", n_eta_RelVals-1, eta_bins, 120,-1.1,1.1);

    book<TH2D>("Jet_eta_chEM", "Jet_eta_chEM;RECO jet #eta; Ch.EM", n_eta_RelVals-1, eta_bins, 110,0.0,1.1);
    book<TH2D>("Jet_eta_neuEM", "Jet_eta_neuEM;RECO jet #eta; Neut.EM", n_eta_RelVals-1, eta_bins, 110,0.0,1.1);
    book<TH2D>("Jet_eta_chHAD", "Jet_eta_chEM;RECO jet #eta; Ch.HAD", n_eta_RelVals-1, eta_bins, 110,0.0,1.1);
    book<TH2D>("Jet_eta_neuHAD", "Jet_eta_neuEM;RECO jet #eta; Neut.HAD", n_eta_RelVals-1, eta_bins, 110,0.0,1.1);

    book<TH1F>("PU_id_BB1","PU id discrimintant",120,-1.1,1.1);
    book<TH1F>("PU_id_BB2","PU id discrimintant",120,-1.1,1.1);
    book<TH1F>("PU_id_EC1","PU id discrimintant",120,-1.1,1.1);
    book<TH1F>("PU_id_EC2","PU id discrimintant",120,-1.1,1.1);
    book<TH1F>("PU_id_HF","PU id discrimintant",120,-1.1,1.1);

    book<TH2D>("Jet_PU_id_chEM", "Jet_PU_id_chEM;RECO jet PU id; Ch.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuEM", "Jet_PU_id_neuEM;RECO jet PU id; Neut.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_chHAD", "Jet_PU_id_chEM;RECO jet PU id; Ch.HAD", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuHAD", "Jet_PU_id_neuEM;RECO jet PU id; Neut.HAD", 44,-1.1,1.1, 22,0.0,1.1);

    book<TH2D>("Jet_PU_id_chEM_BB1", "Jet_PU_id_chEM_BB1;RECO jet PU id; Ch.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuEM_BB1", "Jet_PU_id_neuEM_BB1;RECO jet PU id; Neut.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_chHAD_BB1", "Jet_PU_id_chEM_BB1;RECO jet PU id; Ch.HAD", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuHAD_BB1", "Jet_PU_id_neuEM_BB1;RECO jet PU id; Neut.HAD", 44,-1.1,1.1, 22,0.0,1.1);

    book<TH2D>("Jet_PU_id_chEM_BB2", "Jet_PU_id_chEM_BB2;RECO jet PU id; Ch.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuEM_BB2", "Jet_PU_id_neuEM_BB2;RECO jet PU id; Neut.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_chHAD_BB2", "Jet_PU_id_chEM_BB2;RECO jet PU id; Ch.HAD", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuHAD_BB2", "Jet_PU_id_neuEM_BB2;RECO jet PU id; Neut.HAD", 44,-1.1,1.1, 22,0.0,1.1);

    book<TH2D>("Jet_PU_id_chEM_EC1", "Jet_PU_id_chEM_EC1;RECO jet PU id; Ch.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuEM_EC1", "Jet_PU_id_neuEM_EC1;RECO jet PU id; Neut.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_chHAD_EC1", "Jet_PU_id_chEM_EC1;RECO jet PU id; Ch.HAD", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuHAD_EC1", "Jet_PU_id_neuEM_EC1;RECO jet PU id; Neut.HAD", 44,-1.1,1.1, 22,0.0,1.1);

    book<TH2D>("Jet_PU_id_chEM_EC2", "Jet_PU_id_chEM_EC2;RECO jet PU id; Ch.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuEM_EC2", "Jet_PU_id_neuEM_EC2;RECO jet PU id; Neut.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_chHAD_EC2", "Jet_PU_id_chEM_EC2;RECO jet PU id; Ch.HAD", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuHAD_EC2", "Jet_PU_id_neuEM_EC2;RECO jet PU id; Neut.HAD", 44,-1.1,1.1, 22,0.0,1.1);

    book<TH2D>("Jet_PU_id_chEM_HF", "Jet_PU_id_chEM_HF;RECO jet PU id; Ch.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuEM_HF", "Jet_PU_id_neuEM_HF;RECO jet PU id; Neut.EM", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_chHAD_HF", "Jet_PU_id_chEM_HF;RECO jet PU id; Ch.HAD", 44,-1.1,1.1, 22,0.0,1.1);
    book<TH2D>("Jet_PU_id_neuHAD_HF", "Jet_PU_id_neuEM_HF;RECO jet PU id; Neut.HAD", 44,-1.1,1.1, 22,0.0,1.1);

    //    book<TH2D>("Jet_pt_PUid", "Jet_pt_PUid;RECO jet pt; PU is", 100,0,1000, 120,-1.1,1.1);



}

void JECAnalysisPUjetsHists::fill(const uhh2::Event & ev){
  fill(ev, 0);
}
void JECAnalysisPUjetsHists::fill(const uhh2::Event & ev, const int reco_jet_id){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  // Don't forget to always use the weight when filling.

  //if(ev.isRealData) return;//fill only for MC
  Jet & probe_jet = ev.jets->at(reco_jet_id);// probe RECO jet
  double weight = ev.weight;
  double probe_jet_eta = fabs(probe_jet.eta());
  double probe_jet_pileupID = probe_jet.pileupID();
  double probe_jet_chEM = probe_jet.chargedEmEnergyFraction();
  double probe_jet_neutEM = probe_jet.neutralEmEnergyFraction();
  double probe_jet_chHAD = probe_jet.chargedHadronEnergyFraction();
  double probe_jet_neutHAD = probe_jet.neutralHadronEnergyFraction();
  //  if(probe_jet_chEM<1e-4 && probe_jet_eta<3.0) return;//Skip buggy events. Noise???
  //  if(probe_jet_chEM<1e-4) probe_jet_chEM = 1e-4;//Set small non-zero value for hists
  ((TH2D*)hist("Jet_eta_pt"))->Fill(probe_jet_eta,probe_jet.pt(),weight);
  ((TH2D*)hist("Jet_eta_PUid"))->Fill(probe_jet_eta,probe_jet_pileupID,weight);
  //  ((TH2D*)hist("Jet_pt_PUid"))->Fill(probe_jet.pt(),probe_jet_pileupID,weight);

  ((TH2D*)hist("Jet_eta_chEM"))->Fill(probe_jet_eta,probe_jet_chEM,weight);
  ((TH2D*)hist("Jet_eta_neuEM"))->Fill(probe_jet_eta,probe_jet_neutEM,weight);
  ((TH2D*)hist("Jet_eta_chHAD"))->Fill(probe_jet_eta,probe_jet_chHAD,weight);
  ((TH2D*)hist("Jet_eta_neuHAD"))->Fill(probe_jet_eta,probe_jet_neutHAD,weight);

  ((TH2D*)hist("Jet_PU_id_chEM"))->Fill(probe_jet_pileupID,probe_jet_chEM,weight);
  ((TH2D*)hist("Jet_PU_id_neuEM"))->Fill(probe_jet_pileupID,probe_jet_neutEM,weight);
  ((TH2D*)hist("Jet_PU_id_chHAD"))->Fill(probe_jet_pileupID,probe_jet_chHAD,weight);
  ((TH2D*)hist("Jet_PU_id_neuHAD"))->Fill(probe_jet_pileupID,probe_jet_neutHAD,weight);

  hist("PU_id")->Fill(probe_jet_pileupID, weight);
  hist("nPu")->Fill(ev.get(tt_nPU), weight);
  hist("N_PV")->Fill(ev.get(tt_nvertices), weight);

  if(probe_jet_eta<1.3){
    hist("PU_id_BB1")->Fill(probe_jet_pileupID, weight);
    ((TH2D*)hist("Jet_PU_id_chEM_BB1"))->Fill(probe_jet_pileupID,probe_jet_chEM,weight);
    ((TH2D*)hist("Jet_PU_id_neuEM_BB1"))->Fill(probe_jet_pileupID,probe_jet_neutEM,weight);
    ((TH2D*)hist("Jet_PU_id_chHAD_BB1"))->Fill(probe_jet_pileupID,probe_jet_chHAD,weight);
    ((TH2D*)hist("Jet_PU_id_neuHAD_BB1"))->Fill(probe_jet_pileupID,probe_jet_neutHAD,weight);
    //    if(fabs(probe_jet_pileupID)>1)
    // if(probe_jet_neutEM>0.8)
    //   cout<<"probe_jet_pileupID, chEm, neuEm = "<<probe_jet_pileupID<<", "<<probe_jet_chEM<<" "<<probe_jet_neutEM<<" chHad, neutHad = "<<probe_jet_chHAD<<", "<<probe_jet_neutHAD<<" weight = "<<weight<<" chargedMultiplicity() = "<<probe_jet.chargedMultiplicity()<<", PF sum = "<<probe_jet_chEM+probe_jet_neutEM+probe_jet_chHAD+probe_jet_neutHAD<<endl;
  }
  if(probe_jet_eta>1.3 && probe_jet_eta<2.1){
    hist("PU_id_BB2")->Fill(probe_jet_pileupID, weight);
    ((TH2D*)hist("Jet_PU_id_chEM_BB2"))->Fill(probe_jet_pileupID,probe_jet_chEM,weight);
    ((TH2D*)hist("Jet_PU_id_neuEM_BB2"))->Fill(probe_jet_pileupID,probe_jet_neutEM,weight);
    ((TH2D*)hist("Jet_PU_id_chHAD_BB2"))->Fill(probe_jet_pileupID,probe_jet_chHAD,weight);
    ((TH2D*)hist("Jet_PU_id_neuHAD_BB2"))->Fill(probe_jet_pileupID,probe_jet_neutHAD,weight);
  }
  if(probe_jet_eta>2.1 && probe_jet_eta<2.5){
    hist("PU_id_EC1")->Fill(probe_jet_pileupID, weight);
    ((TH2D*)hist("Jet_PU_id_chEM_EC1"))->Fill(probe_jet_pileupID,probe_jet_chEM,weight);
    ((TH2D*)hist("Jet_PU_id_neuEM_EC1"))->Fill(probe_jet_pileupID,probe_jet_neutEM,weight);
    ((TH2D*)hist("Jet_PU_id_chHAD_EC1"))->Fill(probe_jet_pileupID,probe_jet_chHAD,weight);
    ((TH2D*)hist("Jet_PU_id_neuHAD_EC1"))->Fill(probe_jet_pileupID,probe_jet_neutHAD,weight);
  }
  if(probe_jet_eta>2.5 && probe_jet_eta<3.0){
    hist("PU_id_EC2")->Fill(probe_jet_pileupID, weight);
    ((TH2D*)hist("Jet_PU_id_chEM_EC2"))->Fill(probe_jet_pileupID,probe_jet_chEM,weight);
    ((TH2D*)hist("Jet_PU_id_neuEM_EC2"))->Fill(probe_jet_pileupID,probe_jet_neutEM,weight);
    ((TH2D*)hist("Jet_PU_id_chHAD_EC2"))->Fill(probe_jet_pileupID,probe_jet_chHAD,weight);
    ((TH2D*)hist("Jet_PU_id_neuHAD_EC2"))->Fill(probe_jet_pileupID,probe_jet_neutHAD,weight);
  }
  if(probe_jet_eta>3.0){
    hist("PU_id_HF")->Fill(probe_jet_pileupID, weight);
    ((TH2D*)hist("Jet_PU_id_chEM_HF"))->Fill(probe_jet_pileupID,probe_jet_chEM,weight);
    ((TH2D*)hist("Jet_PU_id_neuEM_HF"))->Fill(probe_jet_pileupID,probe_jet_neutEM,weight);
    ((TH2D*)hist("Jet_PU_id_chHAD_HF"))->Fill(probe_jet_pileupID,probe_jet_chHAD,weight);
    ((TH2D*)hist("Jet_PU_id_neuHAD_HF"))->Fill(probe_jet_pileupID,probe_jet_neutHAD,weight);
  }


}

JECAnalysisPUjetsHists::~JECAnalysisPUjetsHists(){}
