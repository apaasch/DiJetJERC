
#include "UHH2/DiJetJERC/include/JECAnalysisRecoGenMatchedHistsFractions.h"
#include "UHH2/DiJetJERC/include/constants.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/common/include/Utils.h"

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
JECAnalysisRecoGenMatchedHistsFractions::JECAnalysisRecoGenMatchedHistsFractions(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    double eta_bins[n_eta_RelVals];
    for(int i=0; i<n_eta_RelVals; i++) eta_bins[i] = eta_range_RelVals[i];
    const int fr_vec_nr_pf = 110;
    double frac_vec_pf[fr_vec_nr_pf];
    for(int i=0;i<fr_vec_nr_pf;i++){
      frac_vec_pf[i] = i*1.09/fr_vec_nr_pf;
    }
    const int pf_label_nr = 6;
    double pf_label[pf_label_nr];
    for(int i=0;i<pf_label_nr;i++) pf_label[i] = i;
    book<TH3D>("PF_frac_event_eta", "PF_frac_event_eta;RECO jet #eta;PF fraction;", n_eta_RelVals-1, eta_bins, fr_vec_nr_pf-1, frac_vec_pf, pf_label_nr-1, pf_label);
    ((TH3D*)hist("PF_frac_event_eta"))->GetZaxis()->SetBinLabel(1,"Neut. Had");
    ((TH3D*)hist("PF_frac_event_eta"))->GetZaxis()->SetBinLabel(2,"Ch. Had");
    ((TH3D*)hist("PF_frac_event_eta"))->GetZaxis()->SetBinLabel(3,"Neut. EM");
    ((TH3D*)hist("PF_frac_event_eta"))->GetZaxis()->SetBinLabel(4,"Ch. EM");
    ((TH3D*)hist("PF_frac_event_eta"))->GetZaxis()->SetBinLabel(5,"#mu");

    book<TH2D>("PF_frac_event", "PF_frac_event;;PF fraction",pf_label_nr-1, pf_label, fr_vec_nr_pf-1, frac_vec_pf);
    ((TH2D*)hist("PF_frac_event"))->GetXaxis()->SetBinLabel(1,"Neut. Had");
    ((TH2D*)hist("PF_frac_event"))->GetXaxis()->SetBinLabel(2,"Ch. Had");
    ((TH2D*)hist("PF_frac_event"))->GetXaxis()->SetBinLabel(3,"Neut. EM");
    ((TH2D*)hist("PF_frac_event"))->GetXaxis()->SetBinLabel(4,"Ch. EM");
    ((TH2D*)hist("PF_frac_event"))->GetXaxis()->SetBinLabel(5,"#mu");



    const int fr_vec_nr = 202;
    //    const int fr_vec_nr = 109;
    double frac_vec[fr_vec_nr];
    for(int i=0;i<fr_vec_nr;i++){
      frac_vec[i] = i*4.04/fr_vec_nr;
    }

    const int pf_label_nr_had = 7;
    double pf_label_had[pf_label_nr_had];
    for(int i=0;i<pf_label_nr_had;i++) pf_label_had[i] = i;
    book<TH3D>("PF_to_HAD_event_eta", "PF_to_HAD_event_eta;RECO jet #eta;PF jet fraction/Hadron jet fraction;", n_eta_RelVals-1, eta_bins, fr_vec_nr-1, frac_vec, pf_label_nr_had-1,pf_label_had);

    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(1,"(chEM+chHAD)/(GEN charged)");
    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(2,"(neutEM+neutHAD)/(GEN neutral)");
    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(3,"chEM/(GEN charged EM)");
    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(4,"neutEM/(GEN neutral EM)");
    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(5,"chHAD/(GEN charged HAD)");
    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(6,"neutHAD/(GEN neutral HAD)");


    book<TH2D>("Jets_eta_GenVsReco", "Jets_eta_GenVsReco;GEN jet #eta;RECO jet #eta", n_eta_RelVals-1, eta_bins, n_eta_RelVals-1, eta_bins);
    book<TH2D>("Jets_eta_diffGenVsReco", "Jets_eta_diffGenVsReco;GEN jet #eta;(GEN jet #eta-RECO jet #eta)/GEN jet #eta", n_eta_RelVals-1, eta_bins, 120, 0, 1.2);
    // const int dR_vec_nr = 5*24;
    // double dR_vec[dR_vec_nr];
    // for(int i=0;i<dR_vec_nr;i++){
    //   dR_vec[i] = i*1.2/dR_vec_nr;
    // }

    book<TH1D>("Response", "Response;pt^{RECO}/pt^{GEN}", 100, 0, 5.);
    book<TH2D>("Response_eta", "Response; jet #eta; pt^{RECO}/pt^{GEN}", n_eta_RelVals-1, eta_bins,100,0,5.);

}

void JECAnalysisRecoGenMatchedHistsFractions::fill(const uhh2::Event & ev){
  fill(ev, 0, 0);
}
void JECAnalysisRecoGenMatchedHistsFractions::fill(const uhh2::Event & ev, const int gen_jet_id, const int reco_jet_id){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  // Don't forget to always use the weight when filling.

  if(ev.isRealData) return;//fill only for MC
  Jet & probe_jet = ev.jets->at(reco_jet_id);// probe RECO jet
  double dR_GenReco = deltaR(ev.jets->at(reco_jet_id), ev.genjets->at(gen_jet_id));
  // if(dR_GenReco>0.4)  cout<<"Attention., dR is "<<dR_GenReco<<endl;
  GenJet & genj = ev.genjets->at(gen_jet_id);
  double genjet_eta=fabs(genj.eta());
  double dEta = fabs(genjet_eta-probe_jet.eta());
  //  if(genj.genparticles_indices().size()<1) return;//FixME: make the code less dependent on constituents stored, because the most important vars are filled during ntuplization (e.g GEN fractions)
  double weight = ev.weight;
  ((TH2D*)hist("Response_eta"))->Fill(fabs(probe_jet.eta()),probe_jet.pt()/genj.pt(),weight);
  ((TH1D*)hist("Response"))->Fill(probe_jet.pt()/genj.pt(),weight);


  ((TH3D*)hist("PF_frac_event_eta"))->Fill(fabs(probe_jet.eta()),probe_jet.neutralHadronEnergyFraction(),"Neut. Had",weight);
  ((TH3D*)hist("PF_frac_event_eta"))->Fill(fabs(probe_jet.eta()),probe_jet.chargedHadronEnergyFraction(),"Ch. Had",weight);
  ((TH3D*)hist("PF_frac_event_eta"))->Fill(fabs(probe_jet.eta()),probe_jet.neutralEmEnergyFraction(),"Neut. EM",weight);
  ((TH3D*)hist("PF_frac_event_eta"))->Fill(fabs(probe_jet.eta()),probe_jet.chargedEmEnergyFraction(),"Ch. EM",weight);
  ((TH3D*)hist("PF_frac_event_eta"))->Fill(fabs(probe_jet.eta()),probe_jet.muonEnergyFraction(),"#mu",weight);

  ((TH2D*)hist("PF_frac_event"))->Fill("Neut. Had",probe_jet.neutralHadronEnergyFraction(),weight);
  ((TH2D*)hist("PF_frac_event"))->Fill("Ch. Had",probe_jet.chargedHadronEnergyFraction(),weight);
  ((TH2D*)hist("PF_frac_event"))->Fill("Neut. EM",probe_jet.neutralEmEnergyFraction(),weight);
  ((TH2D*)hist("PF_frac_event"))->Fill("Ch. EM",probe_jet.chargedEmEnergyFraction(),weight);
  ((TH2D*)hist("PF_frac_event"))->Fill("#mu",probe_jet.muonEnergyFraction(),weight);


  double chEMpf = probe_jet.chargedEmEnergyFraction();
  double neutEMpf = probe_jet.neutralEmEnergyFraction();
  double chHADpf = probe_jet.chargedHadronEnergyFraction();
  double neutHADpf = probe_jet.neutralHadronEnergyFraction();
  double neutEnergyPF=neutEMpf+neutHADpf;
  double chEnergyPF=chEMpf+chHADpf;


  double chEMgen = genj.cef();
  double chHADgen = genj.chf();
  double neutEMgen = genj.nef();
  double neutHADgen = genj.nhf();
  double neutEnergyGEN=genj.nhf()+genj.nef();
  double chEnergyGEN=genj.chf()+genj.cef();

  ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),chEnergyPF/chEnergyGEN,"(chEM+chHAD)/(GEN charged)",weight);
  ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),neutEnergyPF/neutEnergyGEN,"(neutEM+neutHAD)/(GEN neutral)",weight);

  ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),chEMpf/chEMgen,"chEM/(GEN charged EM)",weight);
  ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),neutEMpf/neutEMgen,"neutEM/(GEN neutral EM)",weight);
  ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),chHADpf/chHADgen,"chHAD/(GEN charged HAD)",weight);
  ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),neutHADpf/neutHADgen,"neutHAD/(GEN neutral HAD)",weight);

  ((TH2D*)hist("Jets_eta_diffGenVsReco"))->Fill(genjet_eta,fabs(genjet_eta-fabs(probe_jet.eta()))/genjet_eta,weight);
  ((TH2D*)hist("Jets_eta_GenVsReco"))->Fill(genjet_eta,fabs(probe_jet.eta()),weight);
}

JECAnalysisRecoGenMatchedHistsFractions::~JECAnalysisRecoGenMatchedHistsFractions(){}
