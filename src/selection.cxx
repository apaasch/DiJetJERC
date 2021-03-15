#include "UHH2/DiJetJERC/include/selection.h"

#include <iostream>
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/L1Jet.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/GenericEvent.h"
#include "UHH2/core/include/PrimaryVertex.h"
// #include "UHH2/BaconTrans/baconheaders/TJet.hh"
//#include "UHH2/BaconTrans/baconheaders/TVertex.hh"
#include "UHH2/DiJetJERC/include/constants.h"

#include "TVector2.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <vector>
#include <stdexcept>

using namespace std;
namespace uhh2DiJetJERC {

  Selection::Selection(uhh2::Context & ctx) : context(ctx), event(0) {

    //FixME: why these variables are declared here?!
    tt_gen_pthat = ctx.declare_event_output<float>("gen_pthat");
    tt_gen_weight = ctx.declare_event_output<float>("gen_weight");
    tt_jet1_pt = ctx.declare_event_output<float>("jet1_pt");
    tt_jet2_pt = ctx.declare_event_output<float>("jet2_pt");
    tt_jet3_pt = ctx.declare_event_output<float>("jet3_pt");
    tt_jet1_ptRaw = ctx.declare_event_output<float>("jet1_ptRaw");
    tt_jet2_ptRaw = ctx.declare_event_output<float>("jet2_ptRaw");
    tt_jet3_ptRaw = ctx.declare_event_output<float>("jet3_ptRaw");
    tt_nvertices = ctx.declare_event_output<int>("nvertices");
    tt_probejet_eta = ctx.declare_event_output<float>("probejet_eta");
    tt_probejet_phi = ctx.declare_event_output<float>("probejet_phi");
    tt_probejet_pt = ctx.declare_event_output<float>("probejet_pt");
    tt_probejet_ptRaw = ctx.declare_event_output<float>("probejet_ptRaw");
    tt_barreljet_eta = ctx.declare_event_output<float>("barreljet_eta");
    tt_barreljet_phi = ctx.declare_event_output<float>("barreljet_phi");
    tt_barreljet_pt = ctx.declare_event_output<float>("barreljet_pt");
    tt_barreljet_ptRaw = ctx.declare_event_output<float>("barreljet_ptRaw");
    tt_pt_ave = ctx.declare_event_output<float>("pt_ave");
    tt_alpha = ctx.declare_event_output<float>("alpha");
    tt_rel_r = ctx.declare_event_output<float>("rel_r");
    tt_mpf_r = ctx.declare_event_output<float>("mpf_r");
    tt_asymmetry = ctx.declare_event_output<float>("asymmetry");
    tt_nPU = ctx.declare_event_output<float>("nPU");
    tt_probejet_chEmEF = ctx.declare_event_output<float>("probejet_chEmEF");

    // dataset_version = ctx.get("dataset_version");
    //
    // if (dataset_version.Contains("RunA")) Cut_Dir = ctx.get("Cut_dir") + "hotjets-18runA.root";
    // else if (dataset_version.Contains("RunB")) Cut_Dir = ctx.get("Cut_dir") + "hotjets-18runB.root";
    // else if (dataset_version.Contains("RunC")) Cut_Dir = ctx.get("Cut_dir") + "hotjets-18runC-ver3.root";
    // else if (dataset_version.Contains("RunD")) Cut_Dir = ctx.get("Cut_dir") + "hotjets-18runD.root";
    // else Cut_Dir = ctx.get("Cut_dir") + "hotjets-18MC.root";

    // cut_map = new TFile(Cut_Dir,"READ");
    // h_map = (TH2D*) cut_map->Get("h2hotfilter");
    // h_map->SetDirectory(0);
    // cut_map->Close();

    try{
      // diJetTrg  = (ctx.get("Trigger_Single") == "false");
      central = (ctx.get("Trigger_Central") == "true");
      fwd     = (ctx.get("Trigger_FWD") == "true");
    }
    catch(const runtime_error& error){
      cout<<"Got runtime error while looking for setting Trigger_Single"<<endl;
      cout << error.what() << "\n";
      cout<<"continue with diJetTrg settings, only relevant if jet trgObj matching is used"<<endl;
      // diJetTrg = true;
      central = true;
      fwd = true;
    }
    // //DEBUG
    // cout<<"\n!!!!!!!! selection diJetTrg "<<diJetTrg<<endl<<endl;
    bool isMC = (ctx.get("dataset_type") == "MC");
    // if (!isMC){
    //   if(!diJetTrg){
    //     handle_trigger40 = ctx.declare_event_input< vector< FlavorParticle > >(  "triggerObjects_hltSinglePFJet40" );
    //     handle_trigger60 = ctx.declare_event_input< vector< FlavorParticle > >(  "triggerObjects_hltSinglePFJet60" );
    //     handle_trigger80 = ctx.declare_event_input< vector< FlavorParticle > >(  "triggerObjects_hltSinglePFJet80" );
    //     handle_trigger140 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet140" );
    //     handle_trigger200 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet200" );
    //     handle_trigger260 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet260" );
    //     handle_trigger320 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet320" );
    //     handle_trigger400 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet400" );
    //     handle_trigger450 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet450" );
    //     handle_trigger500 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet500" );
    //     //dummies
    //     handle_trigger60_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet500" );
    //     handle_trigger80_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet500" );
    //     handle_trigger100_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet500" );
    //     handle_trigger160_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet500" );
    //     handle_trigger220_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet500" );
    //     handle_trigger300_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet500" );
    //   }
    //   else{
    //     if(central){
    //       handle_trigger40 = ctx.declare_event_input< vector< FlavorParticle > >(  "triggerObjects_hltDiPFJetAve40" );
    //       handle_trigger60 = ctx.declare_event_input< vector< FlavorParticle > >(  "triggerObjects_hltDiPFJetAve60" );
    //       handle_trigger80 = ctx.declare_event_input< vector< FlavorParticle > >(  "triggerObjects_hltDiPFJetAve80" );
    //       handle_trigger140 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve140" );
    //       handle_trigger200 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve200" );
    //       handle_trigger260 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve260" );
    //       handle_trigger320 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve320" );
    //       handle_trigger400 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve400" );
    //       handle_trigger500 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve500" );
    //       //as dummy
    //       handle_trigger450 = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltSinglePFJet400" );
    //     }
    //     if(fwd){
    //       handle_trigger60_HF = ctx.declare_event_input< vector< FlavorParticle > >(  "triggerObjects_hltDiPFJetAve60ForHFJEC" );
    //       handle_trigger80_HF = ctx.declare_event_input< vector< FlavorParticle > >(  "triggerObjects_hltDiPFJetAve80ForHFJEC" );
    //       handle_trigger100_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve100ForHFJEC" );
    //       handle_trigger160_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve160ForHFJEC" );
    //       handle_trigger220_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve220ForHFJEC" );
    //       handle_trigger300_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve300ForHFJEC" );
    //     }
    //   }
    // }
    bool no_genp = true;
    if(no_genp) cout<<"\n\n!!! WARNING, no genparticle are used! !!!\n\n"<<endl;

    handle_l1jet_seeds = ctx.declare_event_input< vector< L1Jet>>("L1Jet_seeds");

  }

  void Selection::SetEvent(uhh2::Event& evt) {
    event = &evt;
    assert(event);
  }

  bool Selection::PtMC() {
    assert(event);
    return (event->get(tt_pt_ave) > Pt_AveMC_cut);
  }

  bool Selection::DiJet() {
    assert(event);
    const int njets = event->jets->size();
    return njets >= 2;
  }

  bool Selection::DiJetAdvanced() {
    assert(event);

    const int njets = event->jets->size();
    if (njets < 2) return false;
    Jet* jet1 = &event->jets->at(0);// leading jet
    Jet* jet2 = &event->jets->at(1);// sub-leading jet

    double deltaPhi = std::abs(TVector2::Phi_mpi_pi(jet1->phi() - jet2->phi()));
    if (deltaPhi < s_delta_phi) return false;

    // |asymm| < 1
    if (fabs((event->get(tt_jet2_pt) - event->get(tt_jet1_pt)) / (event->get(tt_jet2_pt) + event->get(tt_jet1_pt))) > s_asymm) return false;


    no_genp=true; //FIXME needed for madgraph, should already be set globally, no idea why it does not work...

    //(pTgen1 < 1.5*pThat || pTreco1 < 1.5* pTgen1)
    // if(!event->isRealData && !no_genp){
    if(!event->isRealData){
      if(event->genjets->size() < 1) return false;
      // if(!(event->genjets->at(0).pt() < 1.5*event->genInfo->binningValues()[0] || event->jets->at(0).pt() < 1.5*event->genjets->at(0).pt())) return false; //FIXME
      if(!(event->genjets->at(0).pt() < 1.5*event->genInfo->qScale() || event->jets->at(0).pt() < 1.5*event->genjets->at(0).pt())) return false;
    }
    return true;
  }

  bool Selection::PUpthat() {
    assert(event);

    // if(no_genp) return true;

    double  gen_pthat = event->genInfo->qScale();
    double  PU_pt_hat = event->genInfo->PU_pT_hat_max();
    return (PU_pt_hat/gen_pthat)<1;

  }

  bool Selection::PtaveVsQScale(double cutValue) {
    assert(event);

    // if(no_genp) return true;

    double  pt_hat = event->genInfo->qScale();
    double ptave = event->get(tt_pt_ave);
    return (ptave/pt_hat)<cutValue;

    return false;
  }

  bool Selection::EnergyEtaCut() {
    // cut away events with jets containing energy more than sqrt(s)/2
    assert(event);
    double probejet_eta = event->get(tt_probejet_eta);
    double ptave = event->get(tt_pt_ave);

    if(ptave*cosh(probejet_eta)>3250) return false; //3250 GeV =sqrt(s)/2 with s=13 TeV
    return true;
  }

  bool Selection::EtaPhiCleaning() {
    assert(event);
    if (((dataset_version.Contains("RunC") || dataset_version.Contains("RunD")) && dataset_version.Contains("18") )) {//TODO for year dependance

      // None of the first 3 jets should be in that region.
      for (size_t i = 0; i < 3; i++) {
        if (event->jets->size()<i+1) continue;
        if (-2.964 < event->jets->at(i).eta() && event->jets->at(i).eta() < -1.305 && -1.6 < event->jets->at(i).phi() && event->jets->at(i).phi() < -0.7) return false;
      }
    }
    return true;

  }


  // bool Selection::ApplyHotMap() {
  //   assert(event);
  //   for(auto jet: *event->jets) {
  //     if(h_map->GetBinContent(h_map->GetXaxis()->FindBin(jet.eta()), h_map->GetYaxis()->FindBin(jet.phi())) > 0) return false;
  //   }
  //   return true;
  // }


  bool Selection::L1JetBXclean(Jet& jet, bool usePtRatioFilter){
    assert(event);

    std::vector< L1Jet>* l1jets = &event->get(handle_l1jet_seeds);
    bool _return = true;

    unsigned int n_l1jets = l1jets->size();

    if(n_l1jets<2) _return = false;

    if(_return){
      double dRmin = 100.;
      int dRmin_seed_idx = -1;
      float dR;

      for(unsigned int i = 0; i<n_l1jets; i++){
        dR=uhh2::deltaR(l1jets->at(i),jet);

        if(dR<0.4 && dR < dRmin){
          dRmin=dR;
          dRmin_seed_idx = i;
        }
      }
      if(dRmin_seed_idx>0){
        if(l1jets->at(dRmin_seed_idx).bx() == -1){
          if(usePtRatioFilter){
            _return = ( l1jets->at(dRmin_seed_idx).pt() / jet.pt() ) < 0.2;
          }
          else _return = false;
        }
      }
    }

    return _return;
  }

  bool Selection::L1JetBXcleanSmart(){
    assert(event);
    bool _return = true;

    std::vector< Jet>* jets = event->jets;
    unsigned int n_jets = jets->size();

    n_jets = std::min(int(n_jets),3);

    for(unsigned int j = 0; j<n_jets && _return ; j++){
      _return *=L1JetBXclean(jets->at(j), true);
    }

    return _return;
  }

  bool Selection::L1JetBXcleanFull(){
    assert(event);
    bool _return = true;

    std::vector< Jet>* jets = event->jets;
    unsigned int n_jets = jets->size();

    for(unsigned int j = 0; j<n_jets && _return ; j++){
      _return *=L1JetBXclean(jets->at(j));
    }

    return _return;
  }

  bool Selection::Unprefirable(std::vector<run_lumi_ev> rlsev){
    assert(event);
    bool _return = false;

    for(size_t i=0; i< rlsev.size(); i++){
      _return = rlsev[i].run == event->run;
      _return *= rlsev[i].lumiblock == event->luminosityBlock;
      _return *= rlsev[i].event == event->event;
      if(_return) break;
    }

    return _return;
  }

  Selection::~Selection() { }


  int Selection::goodPVertex() {
    assert(event);
    Int_t nvertices = event->pvs->size();
    // require in the event that there is at least one reconstructed vertex
    if(nvertices<=0) return 0;//false;
    float nPrVer = 0;
    // pick the first (i.e. highest sum pt) verte
    for (int i=0; i<nvertices; i++){
      PrimaryVertex* vertices = &event->pvs->at(i);
      // require that the vertex meets certain criteria
      //if(vertices->nTracksFit)
      //std::cout<<" vertices->nTracksFit = "<<vertices->nTracksFit<<" "<<fabs(vertices->z)<<" cut at "<<s_n_PvTracks<<std::endl;
      //        if((vertices->nTracksFit > s_n_PvTracks) && (fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy)){
      //SHOULD BE LIKE
      //cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2.")
      // std::cout<<fabs(vertices->z())<<" "<<fabs(vertices->rho())<<" "<<vertices->ndof()<<" "<<vertices->chi2()<<" "<<std::endl;
      if((fabs(vertices->z()) < 24.) && (fabs(vertices->rho()) < 2.) && (vertices->ndof() >= 4)
      && (vertices->chi2()) > 0){
        nPrVer++;
      }
    }
    return nPrVer;
    //    std::cout<<" nPrVer = "<<nPrVer<<" all vtxs = "<<nvertices<<std::endl;
    //    event->set(tt_nGoodvertices,nPrVer);
    // goodVtx = nPrVer;
    // if(nPrVer<=0) return false;
    // else
    //   return true;
  }

  bool Selection::EtaPtCut()
  {
    assert(event);

    double probejet_eta = event->get(tt_probejet_eta);
    double ptave = event->get(tt_pt_ave);

    if( fabs(probejet_eta)>2.500 && fabs(probejet_eta)<3.139 && ptave > 370 ){
      return false;
    }

    return true;
  }

  bool Selection::ChEMFrakCut() {
    assert(event);

    double probejet_eta = event->get(tt_probejet_eta);
    double chEM = event->get(tt_probejet_chEmEF);

    if( fabs(probejet_eta)>2.650 && fabs(probejet_eta)<2.853 && chEM > 0.1 ){
      return false;
    }

    return true;
  }

  int Selection::FindMatchingJet(unsigned int jetid, unsigned int trigger_th, bool use_fwd){
    assert(event);

    float dR_min=0.2;

    float eta = 1000;
    float phi = 1000;

    //DEBUG
    // cout<<"DEBUG in Selection::FindMatchingJet"<<endl;

    uhh2::GenericEvent::Handle<std::vector<FlavorParticle>> handle_sw;

    // if(use_fwd){
    //   switch(trigger_th){
    //   case 60: handle_sw = handle_trigger60_HF;
    // 	break;
    //   case 80: handle_sw = handle_trigger80_HF;
    // 	break;
    //   case 100: handle_sw = handle_trigger100_HF;
    // 	break;
    //   case 160: handle_sw = handle_trigger160_HF;
    // 	break;
    //   case 220: handle_sw = handle_trigger220_HF;
    // 	break;
    //   case 300: handle_sw = handle_trigger300_HF;
    // 	break;
    //   }
    // }
    // else{
    // switch(trigger_th){
    // case 40: handle_sw = handle_trigger40;
    //   break;
    // case 60: handle_sw = handle_trigger60;
    //   break;
    // case 80: handle_sw = handle_trigger80;
    //   break;
    // case 140: handle_sw = handle_trigger140;
    //   break;
    // case 200: handle_sw = handle_trigger200;
    //   break;
    // case 260: handle_sw = handle_trigger260;
    //   break;
    // case 320: handle_sw = handle_trigger320;
    //   break;
    // case 450: handle_sw = handle_trigger450;
    //   break;
    // case 400: handle_sw = handle_trigger400;
    //   break;
    // case 500: handle_sw = handle_trigger500;
    //   break;
    // }
    // }
    //DEBUG
    // cout<<"Selection::FindMatchingJet after switch, jet id: "<<jetid<<" th: "<<trigger_th<<endl;

    if(jetid >= event->get(handle_sw).size()){
      return -1;
    }

    const unsigned int njets = event->jets->size();
    unsigned int jetid_new = -2;

    float dR = 1000.;
    float dR_min_ = dR_min;
    for(unsigned int i = 0 ; i < njets ; i++){
      dR = uhh2::deltaR(event->jets->at(i), event->get(handle_sw).at(jetid));
      if(dR < dR_min_){
        dR_min_ = dR;
        jetid_new = i;
      }
    }

    //DEBUG
    // cout<<jetid_new<<"  "<<dR<<endl;

    return jetid_new;
  }

  // bool Selection::EtaPhiCleaning() {
  //   assert(event);
  //
  //   int n_bins_x = h_map->GetNbinsX();
  //   int n_bins_y = h_map->GetNbinsY();
  //
  //
  //   double xMin = h_map->GetXaxis()->GetXmin();
  //   double xWidth = h_map->GetXaxis()->GetBinWidth(1);
  //
  //
  //   double yMin = h_map->GetYaxis()->GetXmin();
  //   double yWidth = h_map->GetYaxis()->GetBinWidth(1);
  //   double cutValue=0;
  //
  //   const int njets = event->jets->size();
  //
  //   for(int i=0; i < njets; i++){
  //     int idx_x = 0;
  //     int idx_y = 0;
  //     Jet* jet = &event->jets->at(i);// loop over all jets in event
  //
  //     while(jet->eta() > xMin+xWidth + idx_x * xWidth) idx_x++;
  //     while(jet->phi() > yMin+yWidth + idx_y * yWidth) idx_y++;
  //
  //     cutValue = h_map->GetBinContent(idx_x+1, idx_y+1);
  //
  //     if(cutValue > 0) break;
  //   }
  //
  //   if(cutValue > 0) return false;
  //   return true;
  // }

  // bool Selection::EtaPhi() {
  //   assert(event);
  //
  //   double EtaPhi_regions[8][4]={{2.853, 2.964, 0.6, 1.},
  //   {-2.964,-2.853, 0.6, 1.},
  //   {2.853, 2.964, 2.2, 2.6},
  //   {-2.964,-2.853, 2.2, 2.6},
  //   {2.853, 2.964, -2.8, -2.2},
  //   {-2.964,-2.853, -2.8, -2.2},
  //   {-2.964,-2.853, 2.9, 3.1},
  //   {2.853, 2.964, 2.9, 3.1}};
  //
  //   double probejet_eta = event->get(tt_probejet_eta);
  //   double probejet_phi = event->get(tt_probejet_phi);
  //
  //   for(int i=0; i<8; i++){
  //     if(probejet_eta > EtaPhi_regions[i][0] && probejet_eta < EtaPhi_regions[i][1] && probejet_phi > EtaPhi_regions[i][2] && probejet_phi < EtaPhi_regions[i][3]){
  //       //	cout<<"Event rejected!"<<endl<<endl;
  //       return false;
  //     }
  //   }
  //   return true;
  // }

}


// Haikus are easy
// But sometimes they don't make sense
// Refrigerator
