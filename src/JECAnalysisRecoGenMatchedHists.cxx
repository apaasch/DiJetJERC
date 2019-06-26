#include "UHH2/BaconJets/include/JECAnalysisRecoGenMatchedHists.h"
#include "UHH2/BaconJets/include/constants.h"
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
JECAnalysisRecoGenMatchedHists::JECAnalysisRecoGenMatchedHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    book<TH1D>("Hadrons", "Hadrons", 7000, -3500, 3500);
    book<TH1D>("Hadrons_total_energy_rel_event", "Hadrons_total_energy_rel_event; #sum hadrons energy/GEN jet energy", 100,0,2.0);
    book<TH1D>("Hadrons_total_count_event", "Hadrons_total_count_event; Number of hadrons per jet", 100,0,100);


    //    book<TH2D>("Hadrons_energy", "Hadrons_energy;;hadron energy, GeV", 8, 0, 8,1000,1e-6,1000);
    book<TH2D>("Hadrons_energy", "Hadrons_energy;;hadron energy, GeV", 8, 0, 8,2000,1e-6,1000);
    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(1,"#gamma");
    //    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(2,"#pi^{0}");
    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(2,"#pi^{#pm}");
    //    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(4,"K^{0}");
    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(3,"K^{0}_{L}");
    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(4,"K^{0}_{S}");
    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(5,"K^{#pm}");
    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(6,"p");
    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(7,"n");
    ((TH2D*)hist("Hadrons_energy"))->GetXaxis()->SetBinLabel(8,"other");

    book<TH2D>("Hadrons_energy_rel", "Hadrons_energy_rel;;hadron energy/GEN jet energy", 8, 0, 8,101,0,1.01);
    ((TH2D*)hist("Hadrons_energy_rel"))->GetXaxis()->SetBinLabel(1,"#gamma");
    ((TH2D*)hist("Hadrons_energy_rel"))->GetXaxis()->SetBinLabel(2,"#pi^{#pm}");
    ((TH2D*)hist("Hadrons_energy_rel"))->GetXaxis()->SetBinLabel(3,"K^{0}_{L}");
    ((TH2D*)hist("Hadrons_energy_rel"))->GetXaxis()->SetBinLabel(4,"K^{0}_{S}");
    ((TH2D*)hist("Hadrons_energy_rel"))->GetXaxis()->SetBinLabel(5,"K^{#pm}");
    ((TH2D*)hist("Hadrons_energy_rel"))->GetXaxis()->SetBinLabel(6,"p");
    ((TH2D*)hist("Hadrons_energy_rel"))->GetXaxis()->SetBinLabel(7,"n");
    ((TH2D*)hist("Hadrons_energy_rel"))->GetXaxis()->SetBinLabel(8,"other");

    double eta_bins[n_eta_RelVals];                                                                                                                                           
    for(int i=0; i<n_eta_RelVals; i++) eta_bins[i] = eta_range_RelVals[i];  
    book<TH2D>("Hadrons_genjeteta", "Hadrons_genjeteta;;GEN jet #eta (matched to probe RECO jet)", 8,0,8,n_eta_RelVals-1, eta_bins);
    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(1,"#gamma");
    //    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(2,"#pi^{0}");
    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(2,"#pi^{#pm}");
    //    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(4,"K^{0}");
    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(3,"K^{0}_{L}");
    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(4,"K^{0}_{S}");
    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(5,"K^{#pm}");
    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(6,"p");
    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(7,"n");
    ((TH2D*)hist("Hadrons_genjeteta"))->GetXaxis()->SetBinLabel(8,"other");

    // For connection between particles and PF candidates see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlow
    //    book<TH2D>("Hadrons_energy_rel_event", "Hadrons_energy_rel_event;;hadron energy/GEN jet energy (event)", 8, 0, 8,100,1e-6,1.01);
    book<TH2D>("Hadrons_energy_rel_event", "Hadrons_energy_rel_event;;#sum(hadron E)/GEN jet E", 8, 0, 8,101,0,1.01);
    ((TH2D*)hist("Hadrons_energy_rel_event"))->GetXaxis()->SetBinLabel(1,"#gamma");
    ((TH2D*)hist("Hadrons_energy_rel_event"))->GetXaxis()->SetBinLabel(2,"#pi^{#pm}");
    ((TH2D*)hist("Hadrons_energy_rel_event"))->GetXaxis()->SetBinLabel(3,"K^{0}_{L}");
    ((TH2D*)hist("Hadrons_energy_rel_event"))->GetXaxis()->SetBinLabel(4,"K^{0}_{S}");
    ((TH2D*)hist("Hadrons_energy_rel_event"))->GetXaxis()->SetBinLabel(5,"K^{#pm}");
    ((TH2D*)hist("Hadrons_energy_rel_event"))->GetXaxis()->SetBinLabel(6,"p");
    ((TH2D*)hist("Hadrons_energy_rel_event"))->GetXaxis()->SetBinLabel(7,"n");
    ((TH2D*)hist("Hadrons_energy_rel_event"))->GetXaxis()->SetBinLabel(8,"other");

    const int fr_vec_nr_pf = 110;
    double frac_vec_pf[fr_vec_nr_pf];
    for(int i=0;i<fr_vec_nr_pf;i++){
      //      frac_vec_pf[i] = -0.01+i*1.09/fr_vec_nr_pf;
      //      frac_vec_pf[i] = 1e-6+i*1.09/fr_vec_nr_pf;
      frac_vec_pf[i] = i*1.09/fr_vec_nr_pf;
    }

    const int hadron_label_nr = 9;
    //    const int hadron_label_nr = 7;
    double hadron_label[hadron_label_nr];
    for(int i=0;i<hadron_label_nr;i++) hadron_label[i] = i;

    book<TH3D>("Hadrons_energy_rel_event_eta", "Hadrons_energy_rel_event_eta;GEN jet #eta;#sum(hadron E)/GEN jet E;", n_eta_RelVals-1, eta_bins, fr_vec_nr_pf-1, frac_vec_pf, hadron_label_nr-1, hadron_label);
    ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->GetZaxis()->SetBinLabel(1,"#gamma");
    ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->GetZaxis()->SetBinLabel(2,"#pi^{#pm}");
    ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->GetZaxis()->SetBinLabel(3,"K^{0}_{L}");
    ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->GetZaxis()->SetBinLabel(4,"K^{0}_{S}");
    ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->GetZaxis()->SetBinLabel(5,"K^{#pm}");
    ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->GetZaxis()->SetBinLabel(6,"p");
    ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->GetZaxis()->SetBinLabel(7,"n");
    ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->GetZaxis()->SetBinLabel(8,"other");


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

    const int pf_label_nr_had = 3;
    double pf_label_had[pf_label_nr_had];
    for(int i=0;i<pf_label_nr_had;i++) pf_label_had[i] = i;   
    book<TH3D>("PF_to_HAD_event_eta", "PF_to_HAD_event_eta;RECO jet #eta;PF jet fraction/Hadron jet fraction;", n_eta_RelVals-1, eta_bins, fr_vec_nr-1, frac_vec, pf_label_nr_had-1,pf_label_had);

    //    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(1,"(chEM+chHAD)/(#pi^{#pm}+K^{#pm}+p+X^{#pm})");
    // ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(2,"(neutEM+neutHAD)/(#gamma+K^{0}_{L}+K^{0}_{S}+n+X^{0})");

    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(1,"(chEM+chHAD)/(GEN charged)");
    ((TH3D*)hist("PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(2,"(neutEM+neutHAD)/(GEN neutral)");

    const int fr_vec_nr_N = 40;
    double frac_vec_N[fr_vec_nr_N];
    for(int i=0;i<fr_vec_nr_N;i++){
      frac_vec_N[i] = i*10.0/fr_vec_nr_N;
    }
    book<TH3D>("N_PF_to_HAD_event_eta", "N_PF_to_HAD_event_eta;RECO jet #eta;N PF/N hadrons;", n_eta_RelVals-1, eta_bins, fr_vec_nr_N-1, frac_vec_N, pf_label_nr_had-1,pf_label_had);
    ((TH3D*)hist("N_PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(1,"N^{PF}_{charged}/N^{GEN}_{charged}");
    ((TH3D*)hist("N_PF_to_HAD_event_eta"))->GetZaxis()->SetBinLabel(2,"N^{PF}_{neutral}/N^{GEN}_{neutral}");


    book<TH2D>("Hadrons_count_event", "Hadrons_count_event;;number of hadrons (event) ", 8, 0, 8,40,0,40);
    ((TH2D*)hist("Hadrons_count_event"))->GetXaxis()->SetBinLabel(1,"#gamma");
    ((TH2D*)hist("Hadrons_count_event"))->GetXaxis()->SetBinLabel(2,"#pi^{#pm}");
    ((TH2D*)hist("Hadrons_count_event"))->GetXaxis()->SetBinLabel(3,"K^{0}_{L}");
    ((TH2D*)hist("Hadrons_count_event"))->GetXaxis()->SetBinLabel(4,"K^{0}_{S}");
    ((TH2D*)hist("Hadrons_count_event"))->GetXaxis()->SetBinLabel(5,"K^{#pm}");
    ((TH2D*)hist("Hadrons_count_event"))->GetXaxis()->SetBinLabel(6,"p");
    ((TH2D*)hist("Hadrons_count_event"))->GetXaxis()->SetBinLabel(7,"n");
    ((TH2D*)hist("Hadrons_count_event"))->GetXaxis()->SetBinLabel(8,"other");

    book<TH2D>("Jets_eta_GenVsReco", "Jets_eta_GenVsReco;GEN jet #eta;RECO jet #eta", n_eta_RelVals-1, eta_bins, n_eta_RelVals-1, eta_bins);
    book<TH2D>("Jets_eta_diffGenVsReco", "Jets_eta_diffGenVsReco;GEN jet #eta;(GEN jet #eta-RECO jet #eta)/GEN jet #eta", n_eta_RelVals-1, eta_bins, 120, 0, 1.2);
    const int dR_vec_nr = 5*24;
    double dR_vec[dR_vec_nr];
    for(int i=0;i<dR_vec_nr;i++){
      dR_vec[i] = i*1.2/dR_vec_nr;
    }
    book<TH2D>("Hadrons_energy_rel_dRmax", "Hadrons_energy_rel_dRmax;Max allowed #DeltaR(gen.particle,gen.jet);#sum (hadrons energy)/GEN jet energy", dR_vec_nr-1,dR_vec,100,1e-2,2.01);
    book<TH2D>("Hadron_energy_rel_dR", "Hadron_energy_rel_dR;#DeltaR(gen.particle,gen.jet);single hadron energy/GEN jet energy", dR_vec_nr-1,dR_vec,100,1e-2,2.01);

    book<TH1D>("Response", "Response;pt^{RECO}/pt^{GEN}", 100, 0, 5.);
    book<TH2D>("Response_eta", "Response; jet #eta; pt^{RECO}/pt^{GEN}", n_eta_RelVals-1, eta_bins,100,0,5.);


}

void JECAnalysisRecoGenMatchedHists::fill(const uhh2::Event & ev){
  fill(ev, 0, 0);
}
void JECAnalysisRecoGenMatchedHists::fill(const uhh2::Event & ev, const int gen_jet_id, const int reco_jet_id){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  // Don't forget to always use the weight when filling.

  if(ev.isRealData) return;//fill only for MC
  Jet & probe_jet = ev.jets->at(reco_jet_id);// probe RECO jet
  double dR_GenReco = deltaR(ev.jets->at(reco_jet_id), ev.genjets->at(gen_jet_id));
  GenJet & genj = ev.genjets->at(gen_jet_id);
  //  if(genj.genparticles_indices().size()<1) return;//FixME: make the code less dependent on constituents stored, because the most important vars are filled during ntuplization (e.g GEN fractions)
  double weight = ev.weight;
  ((TH2D*)hist("Response_eta"))->Fill(fabs(probe_jet.eta()),probe_jet.pt()/genj.pt(),weight);
  ((TH1D*)hist("Response"))->Fill(probe_jet.pt()/genj.pt(),weight);


  //  if(dR_GenReco>0.2)  cout<<"Attention., dR is "<<dR_GenReco<<endl;
  double genjet_eta=fabs(genj.eta());
  double dEta = fabs(genjet_eta-probe_jet.eta());
  //  if(dEta>1.0) cout<<"Attention., dETA is "<<dEta<<endl;
  double genjet_energy = genj.energy();
  double recjet_energy = probe_jet.energy();

  //To select max_dR ---------
  const int dR_vec_nr = 5*24;
  double dR_vec[dR_vec_nr];
  for(int i=0;i<dR_vec_nr;i++){
    dR_vec[i] = i*1.2/dR_vec_nr;
  }

  for(int i=0;i<dR_vec_nr;i++){
    double genjet_energy_recalcTMP=0;
    for (const auto candInd : genj.genparticles_indices()) {
      GenParticle & genp = ev.genparticles->at(candInd);
      double dr = deltaR(genj,genp); 
      if(dr<dR_vec[i]){
	genjet_energy_recalcTMP += genp.energy();
      }
    }
    //cout<<"dR = "<<dR_vec[i]<<" energy Response = "<<genjet_energy_recalcTMP/genjet_energy<<endl;
    ((TH2D*)hist("Hadrons_energy_rel_dRmax"))->Fill(dR_vec[i],genjet_energy_recalcTMP/genjet_energy,weight);
  }
  //[END] to select max_dR ---
  //  double dr_max=0.2;
  //  double dr_max=0.4;
  //  double dr_max=0.5;

  double dr_max=1000.;
 
  
  double gamma_energy=0; double pipm_energy=0; double K0L_energy=0; 
  double K0S_energy=0; double Kpm_energy=0; double p_energy=0; double n_energy=0;
  double other_energy=0;
  double ch_energy=0; double neut_energy=0;
  double ch_N=0; double neut_N=0;
  int gamma_count=0; int pipm_count=0; int K0L_count=0; 
  int K0S_count=0; int Kpm_count=0; int p_count=0; int n_count=0;
  int other_count=0;
  double genjet_energy_recalc=0;
  //  cout<<"genj.genparticles_indices().size() = "<<genj.genparticles_indices().size()<<endl;
  //  sort_by_pt<GenParticle>(*ev.genparticles);
  for (const auto candInd : genj.genparticles_indices()) {
    GenParticle & genp = ev.genparticles->at(candInd);
    double dr = deltaR(genj,genp); 
    //    cout<<"deltaR(genj,genp) = "<<deltaR(genj,genp)<<endl;
    //    double dr=-1;
    if(dr<dr_max){
      //      int energycl_id = genp.partonFlavour() ;
      int energycl_id = genp.pdgId();
      double energycl_energy = genp.energy();
      ((TH2D*)hist("Hadron_energy_rel_dR"))->Fill(dr,energycl_energy/genjet_energy,weight);
    
      genjet_energy_recalc +=energycl_energy;
      //      cout<<"charge="<<genp.charge()<<" energy="<<energycl_energy<<" energycl_id = "<<energycl_id<<endl;
      if(genp.charge()==0){
	neut_energy+=energycl_energy;
	neut_N++;
      }
      else{ 
	ch_energy+=energycl_energy;
	ch_N++;
      }
      hist("Hadrons")->Fill(energycl_id,weight);
      if(fabs(energycl_id)==22){
	((TH2D*)hist("Hadrons_energy"))->Fill("#gamma",energycl_energy,weight);
	((TH2D*)hist("Hadrons_energy_rel"))->Fill("#gamma",energycl_energy/genjet_energy,weight);
	((TH2D*)hist("Hadrons_genjeteta"))->Fill("#gamma",genjet_eta,weight);
	gamma_energy += energycl_energy;
	gamma_count++;
      }
      if(fabs(energycl_id)==211){
	((TH2D*)hist("Hadrons_energy"))->Fill("#pi^{#pm}",energycl_energy,weight);
	((TH2D*)hist("Hadrons_energy_rel"))->Fill("#pi^{#pm}",energycl_energy/genjet_energy,weight);
	((TH2D*)hist("Hadrons_genjeteta"))->Fill("#pi^{#pm}",genjet_eta,weight);
	pipm_energy += energycl_energy;
	pipm_count++;
      }
      if(fabs(energycl_id)==130){
	((TH2D*)hist("Hadrons_energy"))->Fill("K^{0}_{L}",energycl_energy,weight);
	((TH2D*)hist("Hadrons_energy_rel"))->Fill("K^{0}_{L}",energycl_energy/genjet_energy,weight);
	((TH2D*)hist("Hadrons_genjeteta"))->Fill("K^{0}_{L}",genjet_eta,weight);
	K0L_energy += energycl_energy;
	K0L_count++;
      }
      if(fabs(energycl_id)==310){
	((TH2D*)hist("Hadrons_energy"))->Fill("K^{0}_{S}",energycl_energy,weight);
	((TH2D*)hist("Hadrons_energy_rel"))->Fill("K^{0}_{S}",energycl_energy/genjet_energy,weight);
	((TH2D*)hist("Hadrons_genjeteta"))->Fill("K^{0}_{S}",genjet_eta,weight);
	K0S_energy += energycl_energy;
	K0S_count++;
      }
      if(fabs(energycl_id)==321){
	((TH2D*)hist("Hadrons_energy"))->Fill("K^{#pm}",energycl_energy,weight);
	((TH2D*)hist("Hadrons_energy_rel"))->Fill("K^{#pm}",energycl_energy/genjet_energy,weight);
	((TH2D*)hist("Hadrons_genjeteta"))->Fill("K^{#pm}",genjet_eta,weight);
	Kpm_energy += energycl_energy;
	Kpm_count++;
      }
      if(fabs(energycl_id)==2212){
	((TH2D*)hist("Hadrons_energy"))->Fill("p",energycl_energy,weight);
	((TH2D*)hist("Hadrons_energy_rel"))->Fill("p",energycl_energy/genjet_energy,weight);
	((TH2D*)hist("Hadrons_genjeteta"))->Fill("p",genjet_eta,weight);
	p_energy += energycl_energy;
	p_count++;
      }
      if(fabs(energycl_id)==2112){
	((TH2D*)hist("Hadrons_energy"))->Fill("n",energycl_energy,weight);
	((TH2D*)hist("Hadrons_energy_rel"))->Fill("n",energycl_energy/genjet_energy,weight);
	((TH2D*)hist("Hadrons_genjeteta"))->Fill("n",genjet_eta,weight);
	n_energy += energycl_energy;
	n_count++;
      }
      if(fabs(energycl_id)!=2112 && fabs(energycl_id)!=2212 && fabs(energycl_id)!=321 && fabs(energycl_id)!=310 && fabs(energycl_id)!=130 && fabs(energycl_id)!=211 && fabs(energycl_id)!=111 && fabs(energycl_id)!=22){
	//	cout<<"fabs(energycl_id)="<<fabs(energycl_id)<<endl;
	((TH2D*)hist("Hadrons_energy"))->Fill("other",energycl_energy,weight);
	((TH2D*)hist("Hadrons_energy_rel"))->Fill("other",energycl_energy/genjet_energy,weight);
	((TH2D*)hist("Hadrons_genjeteta"))->Fill("other",genjet_eta,weight);
	other_energy += energycl_energy;
	other_count++;
      }
    }
  }
  double tot_hadron_energy = gamma_energy+pipm_energy+K0L_energy+K0S_energy+Kpm_energy+p_energy+n_energy+other_energy;
  // cout<<" GEN energy = "<<genjet_energy<<", from "<<gamma_count<<" photons:"<<gamma_energy/genjet_energy<<", "<<pipm_count<<" pi^{+/-}: "<<pipm_energy/genjet_energy
  //     <<", "<<K0L_count<<" K0L: "<<K0L_energy/genjet_energy<<", "<<K0S_count<<" K0S_energy:"<<K0S_energy/genjet_energy
  //     <<", "<<Kpm_count<<" K^{+/-}:"<<Kpm_energy/genjet_energy<<", "<<p_count<<" protons: "<<p_energy/genjet_energy
  //     <<", "<<n_count<<" neutrons:"<<n_energy/genjet_energy<<","<<other_count<<" other:"<<other_energy/genjet_energy<<", in total:"<<tot_hadron_energy/genjet_energy<<endl;
  int tot_hadron_count = gamma_count+pipm_count+K0L_count+K0S_count+Kpm_count+p_count+n_count+other_count;
  ((TH2D*)hist("Hadrons_energy_rel_event"))->Fill("#gamma",gamma_energy/genjet_energy,weight);
  ((TH2D*)hist("Hadrons_energy_rel_event"))->Fill("#pi^{#pm}",pipm_energy/genjet_energy,weight);
  ((TH2D*)hist("Hadrons_energy_rel_event"))->Fill("K^{0}_{L}",K0L_energy/genjet_energy,weight);
  ((TH2D*)hist("Hadrons_energy_rel_event"))->Fill("K^{0}_{S}",K0S_energy/genjet_energy,weight);
  ((TH2D*)hist("Hadrons_energy_rel_event"))->Fill("K^{#pm}",Kpm_energy/genjet_energy,weight);
  ((TH2D*)hist("Hadrons_energy_rel_event"))->Fill("p",p_energy/genjet_energy,weight);
  ((TH2D*)hist("Hadrons_energy_rel_event"))->Fill("n",n_energy/genjet_energy,weight);
  ((TH2D*)hist("Hadrons_energy_rel_event"))->Fill("other",other_energy/genjet_energy,weight);

  ///
  ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->Fill(genjet_eta,gamma_energy/genjet_energy,"#gamma",weight);
  ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->Fill(genjet_eta,pipm_energy/genjet_energy,"#pi^{#pm}",weight);
  ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->Fill(genjet_eta,K0L_energy/genjet_energy,"K^{0}_{L}",weight);
  ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->Fill(genjet_eta,K0S_energy/genjet_energy,"K^{0}_{S}",weight);
  ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->Fill(genjet_eta,Kpm_energy/genjet_energy,"K^{#pm}",weight);
  ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->Fill(genjet_eta,p_energy/genjet_energy,"p",weight);
  ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->Fill(genjet_eta,n_energy/genjet_energy,"n",weight);
  ((TH3D*)hist("Hadrons_energy_rel_event_eta"))->Fill(genjet_eta,other_energy/genjet_energy,"other",weight);

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



  double neutEnergyPF=probe_jet.neutralHadronEnergyFraction()+probe_jet.neutralEmEnergyFraction();
  double chEnergyPF=probe_jet.chargedHadronEnergyFraction()+probe_jet.chargedEmEnergyFraction();
  double neutEnergyGEN=neut_energy/genjet_energy;
  double chEnergyGEN=ch_energy/genjet_energy;
  // ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),chEnergyPF/chEnergyGEN,"(chEM+chHAD)/(#pi^{#pm}+K^{#pm}+p+X^{#pm})",weight);
  // ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),neutEnergyPF/neutEnergyGEN,"(neutEM+neutHAD)/(#gamma+K^{0}_{L}+K^{0}_{S}+n+X^{0})",weight);

  ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),chEnergyPF/chEnergyGEN,"(chEM+chHAD)/(GEN charged)",weight);
  ((TH3D*)hist("PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),neutEnergyPF/neutEnergyGEN,"(neutEM+neutHAD)/(GEN neutral)",weight);

  // ((TH3D*)hist("N_PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),double(probe_jet.chargedMultiplicity())/ch_N,"N^{PF}_{charged}/(N_{#pi^{#pm}}+N_{K^{#pm}}+N_{p}+N_{X^{#pm}})",weight);
  // ((TH3D*)hist("N_PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),double(probe_jet.neutralMultiplicity())/neut_N,"N^{PF}_{neutral}/(N_{#gamma}+N_{K^{0}_{L}}+N_{K^{0}_{S}}+N_{n}+N_{X^{0}})",weight);
  //  cout<<"probe_jet.chargedMultiplicity() = "<<probe_jet.chargedMultiplicity()<<" ch_N = "<<ch_N<<" ratio = "<<double(probe_jet.chargedMultiplicity())/ch_N<<endl;
  ((TH3D*)hist("N_PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),double(probe_jet.chargedMultiplicity())/ch_N,"N^{PF}_{charged}/N^{GEN}_{charged}",weight);
  ((TH3D*)hist("N_PF_to_HAD_event_eta"))->Fill(fabs(probe_jet.eta()),double(probe_jet.neutralMultiplicity())/neut_N,"N^{PF}_{neutral}/N^{GEN}_{neutral}",weight);
 
  ((TH2D*)hist("Hadrons_count_event"))->Fill("#gamma",gamma_count,weight);
  ((TH2D*)hist("Hadrons_count_event"))->Fill("#pi^{#pm}",pipm_count,weight);
  ((TH2D*)hist("Hadrons_count_event"))->Fill("K^{0}_{L}",K0L_count,weight);
  ((TH2D*)hist("Hadrons_count_event"))->Fill("K^{0}_{S}",K0S_count,weight);
  ((TH2D*)hist("Hadrons_count_event"))->Fill("K^{#pm}",Kpm_count,weight);
  ((TH2D*)hist("Hadrons_count_event"))->Fill("p",p_count,weight);
  ((TH2D*)hist("Hadrons_count_event"))->Fill("n",n_count,weight);
  ((TH2D*)hist("Hadrons_count_event"))->Fill("other",other_count,weight);

  hist("Hadrons_total_energy_rel_event")->Fill(tot_hadron_energy/genjet_energy,weight);
  hist("Hadrons_total_count_event")->Fill(tot_hadron_count,weight);

  ((TH2D*)hist("Jets_eta_diffGenVsReco"))->Fill(genjet_eta,fabs(genjet_eta-fabs(probe_jet.eta()))/genjet_eta,weight);
  ((TH2D*)hist("Jets_eta_GenVsReco"))->Fill(genjet_eta,fabs(probe_jet.eta()),weight);
}

JECAnalysisRecoGenMatchedHists::~JECAnalysisRecoGenMatchedHists(){}
