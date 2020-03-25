#include <iostream>
#include <memory>
#include <stdlib.h>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/EventHelper.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/L1Jet.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PrintingModules.h"

#include "UHH2/DiJetJERC/include/JECAnalysisHists.h"
#include "UHH2/DiJetJERC/include/JECCrossCheckHists.h"
#include "UHH2/DiJetJERC/include/JECRunnumberHists.h"
#include "UHH2/DiJetJERC/include/JECAnalysisRecoGenMatchedHistsFractions.h"
#include "UHH2/DiJetJERC/include/JECAnalysisPUjetsHists.h"
#include "UHH2/DiJetJERC/include/LumiHists.h"
#include "UHH2/DiJetJERC/include/selection.h"
#include "UHH2/DiJetJERC/include/constants.h"
#include "TClonesArray.h"
#include "TString.h"
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2Poly.h"
#include "UHH2/common/include/CollectionProducer.h"

using namespace std;
using namespace uhh2;

class testJER: public uhh2::AnalysisModule {

public:
  explicit testJER(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  ~testJER();

protected:

  // correctors
  std::unordered_map<std::string, std::unique_ptr<GenericJetCorrector> > JetCorr;
  std::unique_ptr<GenericJetResolutionSmearer> jet_resolution_smearer_nominal, jet_resolution_smearer_down, jet_resolution_smearer_up;

  std::unique_ptr<AnalysisModule> CollectionProducer_module_nominal;
  std::unique_ptr<AnalysisModule> CollectionProducer_module_down;
  std::unique_ptr<AnalysisModule> CollectionProducer_module_up;
  // cleaners
  std::unique_ptr<JetCleaner> jetID;

  // selections
  std::unique_ptr<uhh2::AnalysisModule> PVCleaner;
  std::unique_ptr<uhh2::AndSelection> metfilters_sel;

  //// Data/MC scale factors
  std::unique_ptr<uhh2::AnalysisModule> pileupSF;

  Event::Handle<double> jet_pt;
  Event::Handle<double> jet_pt_nominal;
  Event::Handle<double> jet_pt_down;
  Event::Handle<double> jet_pt_up;

  Event::Handle<double> jet_eta;
  Event::Handle<double> jet_eta_nominal;
  Event::Handle<double> jet_eta_down;
  Event::Handle<double> jet_eta_up;

  Event::Handle<vector<Jet>> jets_nominal;
  Event::Handle<vector<Jet>> jets_down;
  Event::Handle<vector<Jet>> jets_up;


  //useful booleans
  bool debug, no_genp;
  bool isMC, JECClosureTest, JERClosureTest, apply_EtaPhi_cut, apply_EtaPhi_HCAL, trigger_central, trigger_fwd, DO_Pu_ReWeighting, apply_lumiweights, apply_L1seed_from_bx1_filter, apply_PUid;
  bool is2016v2, is2016v3, is2017, is2018;
  std::unordered_map<std::string, std::vector<std::string>> runs = { {"2016", runPeriods2016}, {"2017", runPeriods2017}, {"UL17", runPeriods2017}, {"2018", runPeriods2018}};
  std::string year;
  bool isAK8, ispuppi;
  string SysType_PU;
  TString dataset_version, jetLabel;
  string JEC_Version, jecTag, jecVer, JEC_Level, jet_coll;


};

testJER::testJER(uhh2::Context & ctx) {

  no_genp = true;
  dataset_version = ctx.get("dataset_version");
  isMC = (ctx.get("dataset_type") == "MC");
  jetLabel = ctx.get("JetLabel");
  isAK8 = (jetLabel == "AK8CHS" || jetLabel == "AK8Puppi");
  ispuppi = (jetLabel == "AK4Puppi" || jetLabel == "AK8Puppi");
  jet_coll = isAK8? "AK8" : "AK4"; jet_coll += ispuppi? "PFPuppi" : "PFchs";
  JEC_Level = ctx.get("JEC_Level", "L1L2L3Residual");

  JECClosureTest = true;
  JERClosureTest = string2bool(ctx.get("JERClosureTest","false"));
  apply_EtaPhi_cut = string2bool(ctx.get("EtaPhi_cut", "true"));
  apply_EtaPhi_HCAL = string2bool(ctx.get("EtaPhi_HCAL", "true"));
  apply_PUid = string2bool(ctx.get("Apply_PUid_3rdjet", "false"));
  cout << "Dataset is " << ((isMC) ? " mc " : " data") << endl;

  apply_L1seed_from_bx1_filter =  (ctx.get("Apply_L1Seed_From_BX1_Filter","false") == "true" && !isMC);
  is2016v2 = (ctx.get("dataset_version").find("2016v2") != std::string::npos);
  is2016v3 = (ctx.get("dataset_version").find("2016v3") != std::string::npos);
  is2017 = (ctx.get("dataset_version").find("2017") != std::string::npos);
  is2018 = (ctx.get("dataset_version").find("2018") != std::string::npos);
  year = (is2016v2 || is2016v3)? "2016": (is2017?"2017" : (is2018?"2018" :"" ));

  JEC_Version = is2017?"Fall17_17Nov2017_V32" : (is2018?"Autumn18_V19" :"" );
  std::cout << "JEC_Version " << JEC_Version << '\n';
  jecTag = JEC_Version.substr(0,JEC_Version.find("_V"));
  jecVer = JEC_Version.substr(JEC_Version.find("_V")+2,JEC_Version.size()-JEC_Version.find("_V")-2);

  runs[year].push_back("MC");

  debug = string2bool(ctx.get("Debug","false"));
  // debug= true;

  PVCleaner.reset(new PrimaryVertexCleaner());
  /* MET filters */
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
  PrimaryVertexId pvid = StandardPrimaryVertexId();
  metfilters_sel.reset(new AndSelection(ctx, "metfilters"));
  metfilters_sel->add<TriggerSelection>("goodVertices", "Flag_goodVertices");
  metfilters_sel->add<NPVSelection>("1 good PV",1,-1,pvid); /* Not a metfilter. Used to select 1 good PV */
  metfilters_sel->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter");
  metfilters_sel->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_sel->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_sel->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  if (year != "2016") metfilters_sel->add<EcalBadCalibSelection>("EcalBadCalibSelection"); /*TODO check 2016*/ // Use this instead of Flag_ecalBadCalibFilter, uses ecalBadCalibReducedMINIAODFilter in ntuple_generator
  if (year != "2016") metfilters_sel->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter"); /*TODO check 2016, maybe Extra_BadPFMuonFilter */
  if (!isMC) metfilters_sel->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter"); /* TODO Not recommended for MC, but do check */
  /* metfilters_sel->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter"); TODO Not recommended, under review.Separate module in ntuple_generator for 2016v2*/

  for (const std::string & run : runs[year]) {
    //std::vector<std::string> JEC_corr = (run=="MC")? JERFiles::JECFilesMC(jecTag, jecVer, jet_coll) : JERFiles::JECFilesDATA(jecTag, jecVer, jet_coll, run,JECClosureTest? JERFiles::L1L2L3Residual : JERFiles::L1L2L3);//TODO
    // JetCorr[run].reset(new GenericJetCorrector(ctx, JEC_corr,"jets"));
    if (run=="MC") JetCorr[run].reset(new GenericJetCorrector(ctx, JERFiles::JECFilesMC(jecTag, jecVer, jet_coll),"jets"));
    else {
      if (JEC_Level=="L1L2") JetCorr[run].reset(new GenericJetCorrector(ctx, JERFiles::JECFilesDATA(jecTag, jecVer, jet_coll, run, JECClosureTest? JERFiles::L1L2 : JERFiles::L1L2L3),"jets"));
      else if (JEC_Level=="L1L2Residual") JetCorr[run].reset(new GenericJetCorrector(ctx, JERFiles::JECFilesDATA(jecTag, jecVer, jet_coll, run, JECClosureTest? JERFiles::L1L2Residual : JERFiles::L1L2L3),"jets"));
      else if (JEC_Level=="L1L2L3Residual") JetCorr[run].reset(new GenericJetCorrector(ctx, JERFiles::JECFilesDATA(jecTag, jecVer, jet_coll, run, JECClosureTest? JERFiles::L1L2L3Residual : JERFiles::L1L2L3),"jets"));
      else throw std::invalid_argument(JEC_Level+" is not implemented");
    }
  }

  const JetId jetId(AndId<Jet> (JetPUid(JetPUid::WP_TIGHT), ispuppi?JetPFID(JetPFID::WP_TIGHT_PUPPI):JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30, 6)));
  jetID.reset(new JetCleaner(ctx, jetId));


  //JER Smearing for corresponding JEC-Version

  if (is2018) {
    jet_resolution_smearer_nominal.reset(new GenericJetResolutionSmearer(ctx, 0, "nominal", "genjets", JERSmearing::SF_13TeV_Autumn18_RunABCD_V4, "2018/Autumn18_V4_MC_PtResolution_AK4PFchs.txt"));
    jet_resolution_smearer_down.reset(new GenericJetResolutionSmearer(ctx, -1, "down", "genjets", JERSmearing::SF_13TeV_Autumn18_RunABCD_V4, "2018/Autumn18_V4_MC_PtResolution_AK4PFchs.txt"));
    jet_resolution_smearer_up.reset(new GenericJetResolutionSmearer(ctx, 1, "up", "genjets", JERSmearing::SF_13TeV_Autumn18_RunABCD_V4, "2018/Autumn18_V4_MC_PtResolution_AK4PFchs.txt"));
  }

  if (is2017) {
    jet_resolution_smearer_nominal.reset(new GenericJetResolutionSmearer(ctx, 0, "nominal", "genjets", JERSmearing::SF_13TeV_Fall17_V3_RunBCDEF_Madgraph,"2017/Fall17_V3_MC_PtResolution_AK4PFchs.txt"));
    jet_resolution_smearer_down.reset(new GenericJetResolutionSmearer(ctx, -1, "down", "genjets", JERSmearing::SF_13TeV_Fall17_V3_RunBCDEF_Madgraph,"2017/Fall17_V3_MC_PtResolution_AK4PFchs.txt"));
    jet_resolution_smearer_up.reset(new GenericJetResolutionSmearer(ctx, 1, "up", "genjets", JERSmearing::SF_13TeV_Fall17_V3_RunBCDEF_Madgraph,"2017/Fall17_V3_MC_PtResolution_AK4PFchs.txt"));
  }

  CollectionProducer_module_nominal.reset(new CollectionProducer<Jet>( ctx, "jets", "nominal"));
  CollectionProducer_module_down.reset(new CollectionProducer<Jet>( ctx, "jets", "down"));
  CollectionProducer_module_up.reset(new CollectionProducer<Jet>( ctx, "jets", "up"));


  jets_nominal = ctx.get_handle<vector<Jet>>("nominal");
  jets_down = ctx.get_handle<vector<Jet>>("down");
  jets_up = ctx.get_handle<vector<Jet>>("up");

  //output
  ctx.undeclare_all_event_output();
  jet_pt = ctx.declare_event_output<double>("jet_pt");
  jet_pt_nominal = ctx.declare_event_output<double>("jet_pt_nominal");
  jet_pt_down = ctx.declare_event_output<double>("jet_pt_down");
  jet_pt_up = ctx.declare_event_output<double>("jet_pt_up");

  jet_eta = ctx.declare_event_output<double>("jet_eta");
  jet_eta_nominal = ctx.declare_event_output<double>("jet_eta_nominal");
  jet_eta_down = ctx.declare_event_output<double>("jet_eta_down");
  jet_eta_up = ctx.declare_event_output<double>("jet_eta_up");

  // Do pileup reweighting (define it after undeclaring all other variables to keep the weights in the output)
  apply_lumiweights = string2bool(ctx.get("apply_lumiweights","true"));
  apply_lumiweights = apply_lumiweights && isMC;

  SysType_PU = ctx.get("SysType_PU");
  pileupSF.reset(new MCPileupReweight(ctx,SysType_PU));

  cout<<"end of AnalyseModule Constructor" << '\n';

};

testJER::~testJER() {
}

bool testJER::process(Event & event) {

  //###############################################################
  //
  //Selection Module for DiJetJERC Calculation
  //
  //Select Di-Jet Events
  //Define Barrel and Probe Jet
  //Use possible third Jet to estimate alpha
  //Apply MC-Weights for Reweighting (Second Iteration)
  //
  //###############################################################

  #define ak4jets event.jets

  if(!pileupSF->process(event)) return false;

  if (!PVCleaner->process(event)) return false;


  // MET filters
  if(!metfilters_sel->passes(event)) return false;

  jetID->process(event);
  if (ak4jets->size()<1) return false;

  sort_by_pt<Jet>(*ak4jets);


  std::unordered_map<std::string, bool > apply_run;
  for (const std::string & run : runs[year]) apply_run[run] = false;

  bool apply_all = false;
  if (!isMC) {
    for (const std::string & run : runs[year]) {
      if (run=="MC") continue;
      if (year=="UL17") {
        if (run_number_map.at("2017").at(run).first <= event.run && event.run <= run_number_map.at("2017").at(run).second) apply_run[run] = true;
      } else {
        if (run_number_map.at(year).at(run).first <= event.run && event.run <= run_number_map.at(year).at(run).second) apply_run[run] = true;
      }
      apply_all+=apply_run[run];
    }
  } else {apply_run["MC"] = true; apply_all+=apply_run["MC"];}
  if (apply_all != 1) throw std::runtime_error("In testJER.cxx: Sum of apply_all when applying JECs is not == 1. Fix this.");


  for (const std::string run : runs[year]) {
    if (apply_run[run]){
      JetCorr[run]->process(event);
    }
  }

  //Apply JER to all jet collections
  CollectionProducer_module_nominal->process(event);
  CollectionProducer_module_down->process(event);
  CollectionProducer_module_up->process(event);

  jet_resolution_smearer_nominal->process(event);
  jet_resolution_smearer_down->process(event);
  jet_resolution_smearer_up->process(event);


  for (size_t i = 0; i < ak4jets->size(); i++) {
    auto jet = ak4jets->at(i);
    auto jet_nominal = event.get(jets_nominal).at(i);
    auto jet_down = event.get(jets_down).at(i);
    auto jet_up = event.get(jets_up).at(i);

    event.set(jet_pt,jet.pt());
    event.set(jet_eta,jet.eta());
    event.set(jet_pt_nominal,jet_nominal.pt());
    event.set(jet_eta_nominal,jet_nominal.eta());
    event.set(jet_pt_down,jet_down.pt());
    event.set(jet_eta_down,jet_down.eta());
    event.set(jet_pt_up,jet_up.pt());
    event.set(jet_eta_up,jet_up.pt());
  }


  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the ExampleModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(testJER)
