#pragma once
//#include <UHH2/BaconJets/include/HistsBASE.h>
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "TClonesArray.h"
//#include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"
//#include "UHH2/BaconTrans/baconheaders/BaconAnaDefs.hh"
/**  \brief Example class for booking and filling histograms
*
* NOTE: This class uses the 'hist' method to retrieve histograms.
* This requires a string lookup and is therefore slow if you have
* many histograms. Therefore, it is recommended to use histogram
* pointers as member data instead, like in 'common/include/ElectronHists.h'.
*/
class JECAnalysisHists: public uhh2::Hists {
  //class JECAnalysisHists: public HistsBASE {
public:
  // use the same constructor arguments as Hists for forwarding:
  JECAnalysisHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  void fill(const uhh2::Event & ev, const int rand);
  virtual ~JECAnalysisHists();

  // TH2D *ptjet3_vs_alpha,*pt_ave_vs_alpha;
private:
  //    uhh2::Event::Handle<TClonesArray> h_jets;
  //    uhh2::Event::Handle<baconhep::TEventInfo> h_eventInfo;
  //    uhh2::Event::Handle<TClonesArray> h_pv;
  //Variables specific for dijet events
  uhh2::Event::Handle<float> tt_gen_pthat; uhh2::Event::Handle<float> tt_gen_weight;
  uhh2::Event::Handle<float> tt_jet1_pt;     uhh2::Event::Handle<float> tt_jet2_pt;     uhh2::Event::Handle<float> tt_jet3_pt;
  uhh2::Event::Handle<float> tt_jet1_ptRaw;  uhh2::Event::Handle<float> tt_jet2_ptRaw;  uhh2::Event::Handle<float> tt_jet3_ptRaw;
  uhh2::Event::Handle<int> tt_nvertices;
  uhh2::Event::Handle<float> tt_probejet_eta;  uhh2::Event::Handle<float> tt_probejet_phi; uhh2::Event::Handle<float> tt_probejet_pt; uhh2::Event::Handle<float> tt_probejet_ptRaw;
  uhh2::Event::Handle<float> tt_barreljet_eta;  uhh2::Event::Handle<float> tt_barreljet_phi; uhh2::Event::Handle<float> tt_barreljet_pt; uhh2::Event::Handle<float> tt_barreljet_ptRaw;
  uhh2::Event::Handle<float> tt_pt_ave;
  uhh2::Event::Handle<float> tt_alpha; uhh2::Event::Handle<float> tt_alpha_sum;
  uhh2::Event::Handle<float> tt_rel_r; uhh2::Event::Handle<float> tt_mpf_r; uhh2::Event::Handle<float> tt_asymmetry; uhh2::Event::Handle<float> tt_nPU;
  uhh2::Event::Handle<int> tt_matchJetId_0; uhh2::Event::Handle<int> tt_matchJetId_1;

};
