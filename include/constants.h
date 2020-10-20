#pragma once

#include <TString.h>
#include <map>
#include <array>

// barrel region (|eta| < 1.131)
constexpr static double s_eta_barr = 1.131;
// HF region (|eta| > 2.853)
const double eta_cut = 2.853;
// two back-to-back leading jets (delta_phi(j1,j2) = min(|phi1 - phi2|, 2PI - |phi2 - phi1|) > 2.7)
constexpr static double s_delta_phi = 2.7;

constexpr static double jet_threshold      = 15;
constexpr static double jet_threshold_min  = 10;

const int n_eta_bins_JER = 15;
const double eta_bins_JER[n_eta_bins_JER] = { 0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.650, 2.853, 2.964, 3.139, 5.191 };

constexpr static double s_delta_R   = 0.3;

//Eta bins:
//Abs eta range:
const int n_eta = 19;
const TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};
const TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};
const double eta_bins[n_eta]     = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

//Negative Eta Range
const int n_eta_full = 37;
const TString eta_range_full[n_eta_full] = {"-5.191","-3.839","-3.489","-3.139","-2.964","-2.853", "-2.65", "-2.5", "-2.322", "-2.172", "-1.93", "-1.653", "-1.479", "-1.305", "-1.044", "-0.783", "-0.522", "-0.261"," 0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};
const TString eta_range2_full[n_eta_full] = {"-5191","-3839","-3489","-3139","-2964","-2853", "-265", "-25", "-2322", "-2172", "-193", "-1653", "-1479", "-1305", "-1044", "-0783", "-0522", "-0261","00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};
const double eta_bins_full[n_eta_full]     = {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -1.93, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

const int n_pt_Si = 12;

/* //SingleJet triggers highest double checked with combined BCDEF------------- */
const double pt_bins_Si[n_pt_Si] = {
  40 ,
  72 ,
  95 ,
  160,
  226,
  283,
  344,
  443,
  577,
  606,1000,2000};

  const TString pt_range_Si[n_pt_Si]={
    "40",
    "72",
    "95",
    "160",
    "226",
    "283",
    "344",
    "443",
    "577",
    "606","1000","2000"};

    /* //SingleJet triggers highest End------------- */

    //19 bin edges, 18 actual bins
    //constexpr static int n_eta = 19;
    //static std::vector<double>   eta_range  =  {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
    static std::vector<double>   eta_range_mikko  = {0, 0.783, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
    //Eta bins for GlobalFit files
    const int n_eta_common = 19;
    const double eta_common_bins[n_eta_common] ={0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
    const TString eta_common_range[n_eta_common] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};

    const TString eta_output[n_eta_common-1] = {"eta_00_03", "eta_03_05","eta_05_08","eta_08_10","eta_10_13","eta_13_15","eta_15_17", "eta_17_19", "eta_19_22", "eta_22_23", "eta_23_25", "eta_25_27", "eta_27_29", "eta_29_30", "eta_30_31", "eta_31_35", "eta_35_38", "eta_38_52"};

    const int n_eta_common_2 = 7;
    const double eta_common_bins_2[n_eta_common_2] ={
      0.,
      1.305,
      1.93,
      2.5,
      2.964,
      3.2,
      5.191
    };
    const TString eta_common_range_2[n_eta_common_2] = {
      "0.000",
      "1.305",
      "1.930",
      "2.500",
      "2.964",
      "3.200",
      "5.191"};

      const TString eta_output_2[n_eta_common_2-1] = {"eta_00_13", "eta_13_19","eta_19_25","eta_25_30","eta_30_32","eta_32_52"};


      //eta binning for checks with RelVals
      constexpr static int n_eta_RelVals = 19;
      static std::vector<double>   eta_range_RelVals  =  {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};


      //alpha
      //static std::vector<double>   alpha_range= {0., 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};
      const double alpha_cut = 0.3;
      const TString s_alpha_cut = "0.3";
      const TString s_alpha_cut_name = "03";
      const int n_alpha = 11;
      const TString alpha_range[n_alpha] = {"a005", "a010", "a015", "a020", "a025", "a030", "a035", "a040", "a045","a050","a055"};
      const double alpha_bins[n_alpha] = {0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,  0.400, 0.450,0.50,0.55};

      const int n_alpha_common = 4;
      const TString alpha_range_common[n_alpha_common] = {"a10","a15", "a20", "a30"};//for Global fit we produce result only in few alpha values
      const double alpha_bins_common[n_alpha_common] = {0.100, 0.150, 0.200, 0.300};

      /** \brief Dijet event selection **/

      // cut on the extreme asymmetry for events with two jets  |(j2->pt - j1->pt /(j2->pt + j1->pt)| < 0.70
      constexpr static float s_asymm = 1.;

      /** \brief good Primary Vertex reconstruction **/
      // more than four tracks
      constexpr static float s_n_PvTracks = 4;
      // PV is located within 24cm in z vertex
      constexpr static float s_n_Pv_z = 24.0;
      // PV is located within 2cm in xy direction from the nominal interaction point
      constexpr static float s_n_Pv_xy = 2.0;


      /** \brief The trigger thresholds of pt_ave **/

      constexpr static float Pt_AveMC_cut = 51; // TODO to generalise for years
      constexpr static float PT_trigger_max = 7000; // TODO to generalise for years


      const std::map<std::string, std::vector<std::string> > pt_indexes = {
        {"DiJet_central",     {"trigger40", "trigger60", "trigger80", "trigger140", "trigger200", "trigger260", "trigger320", "trigger400", "trigger500"}},
        {"SingleJet_central", {"trigger40", "trigger60", "trigger80", "trigger140", "trigger200", "trigger260", "trigger320", "trigger400", "trigger450", "trigger500"}}, // TODO add trigger550 (may vary for each year)
        {"DiJet_forward",     {"trigger60_HFJEC", "trigger80_HFJEC",  "trigger100_HFJEC", "trigger160_HFJEC", "trigger220_HFJEC", "trigger300_HFJEC"}},
        {"SingleJet_forward", {"trigger40_HFJEC", "trigger60_HFJEC",  "trigger80_HFJEC",  "trigger140_HFJEC", "trigger200_HFJEC", "trigger260_HFJEC", "trigger320_HFJEC", "trigger400_HFJEC"}},
      };

      const std::map<std::string, std::vector<double> > pt_trigger_thr = {
        // 2017 AK4
        {"DiJet_central_2017",                { 73, 85, 97, 179, 307, 370, 434, 520, 649 }}, //for Di triggers 94X 17Nov2017
        {"DiJet_central_UL17",                { 73, 85, 97, 179, 307, 370, 434, 520, 649 }}, //for Di triggers 94X 17Nov2017
        {"DiJet_forward_2017",                { 73, 93, 113, 176, 239, 318 }}, // for Dijet_HFJEC 2017 94X 17Nov2017
        {"DiJet_forward_UL17",                { 73, 93, 113, 176, 239, 318 }}, // for Dijet_HFJEC 2017 94X 17Nov2017
        {"SingleJet_central_2017",            { 40, 72,  95, 160, 226, 283, 344, 443, 577, 606 }}, //for Single triggers 94X 17Nov2017
        {"SingleJet_central_UL17",            { 70, 87, 111, 180, 247, 310, 373, 457, 510, 562 }}, //from https://indico.cern.ch/event/887399/contributions/3741486/attachments/1984788/3306787/triggerUL17RunBCDEF.pdf
        {"SingleJet_central_UL17_ptbins",     { 70, 87, 111, 180, 247, 310, 373, 457, 510, 562 }}, //from https://indico.cern.ch/event/887399/contributions/3741486/attachments/1984788/3306787/triggerUL17RunBCDEF.pdf
        {"SingleJet_forward_2017",            { 73, 93, 113, 176, 239, 318 }}, // for Singlejet_HFJEC 2017 94X 17Nov2017 //TODO not used
        {"SingleJet_forward_UL17",            { 73, 93, 113, 176, 239, 318 }}, // for Singlejet_HFJEC 2017 94X 17Nov2017 //TODO not used
        {"SingleJet_forward_UL17_ptbins",     { 40, 72,  95, 160, 226, 283, 344, 443, 577, 606 }}, //for Single triggers 94X 17Nov2017
        //{"SingleJet_forward_UL17",            { 70, 87, 111, 180, 247, 310, 373, 457, 510, 562 }}, //from https://indico.cern.ch/event/887399/contributions/3741486/attachments/1984788/3306787/triggerUL17RunBCDEF.pdf
        //{"SingleJet_forward_UL17_ptbins",     { 70, 87, 111, 180, 247, 310, 373, 457, 510, 562 }}, //from https://indico.cern.ch/event/887399/contributions/3741486/attachments/1984788/3306787/triggerUL17RunBCDEF.pdf
        // 2017 AK8
        {"SingleJet_central_AK8_2017",        { 73, 90, 115, 181, 251, 312, 378, 457, 519, 566 }}, //for Single triggers 94X 17Nov2017
        {"SingleJet_central_AK8_UL17",        { 73, 90, 115, 181, 251, 312, 378, 457, 519, 566 }}, //for Single triggers 94X 17Nov2017
        // 2018 AK4
        {"DiJet_central_2018",                { 66, 93, 118, 189, 257, 325, 391, 478, 585 }}, //for Di triggers 2018, RunABC, ReReco https://indico.cern.ch/event/801509/contributions/3331436/attachments/1801472/2938522/L2Res-Triggers-25Feb2019.pdf
        {"DiJet_central_2018_ptbins",         { 66, 93, 118, 189, 257, 291, 325, 358, 391, 434, 478, 531, 585}},
        {"DiJet_forward_2018",                { 93, 116, 142, 210, 279, 379 }},
        {"DiJet_forward_2018_ptbins",         { 93, 116, 142, 210, 279, 379 }},
        {"DiJet_central_UL18",                { 66, 93, 118, 189, 257, 325, 391, 478, 585 }}, //for Di triggers 2018, RunABC, ReReco https://indico.cern.ch/event/801509/contributions/3331436/attachments/1801472/2938522/L2Res-Triggers-25Feb2019.pdf
        {"DiJet_central_UL18_ptbins",         { 66, 93, 118, 189, 257, 291, 325, 358, 391, 434, 478, 531, 585}},
        {"DiJet_forward_UL18",                { 93, 116, 142, 210, 279, 379 }},
        {"DiJet_forward_UL18_ptbins",         { 93, 116, 142, 210, 279, 379 }},
        // 2018 AK8
        {"SingleJet_central_AK8_2018",        { 78, 96, 119, 193, 262, 328, 393, 481, 534, 588 }}, //HLT AK8PFJet* //https://indico.desy.de/indico/event/24350/contribution/0/material/slides/0.pdf
        {"SingleJet_central_AK8_2018_ptbins", { 78, 96, 119, 193, 262, 328, 393, 481, 534, 588 }},
        {"SingleJet_forward_AK8_2018",        { 62, 95, 110, 182, 260, 339, 420, 508 }}, // HLT AK8PFJetFwd* //https://indico.desy.de/indico/event/24423/contribution/1/material/slides/1.pdf
        {"SingleJet_forward_AK8_2018_ptbins", { 62, 95, 110, 182, 260, 339, 420, 508 }},
        {"SingleJet_central_AK8_UL18",        { 78, 96, 119, 193, 262, 328, 393, 481, 534, 588 }}, //HLT AK8PFJet* //https://indico.desy.de/indico/event/24350/contribution/0/material/slides/0.pdf
        {"SingleJet_central_AK8_UL18_ptbins", { 78, 96, 119, 193, 262, 328, 393, 481, 534, 588 }},
        {"SingleJet_forward_AK8_UL18",        { 62, 95, 110, 182, 260, 339, 420, 508 }}, // HLT AK8PFJetFwd* //https://indico.desy.de/indico/event/24423/contribution/1/material/slides/1.pdf
        {"SingleJet_forward_AK8_UL18_ptbins", { 62, 95, 110, 182, 260, 339, 420, 508 }},
      };


      // RunII pt-bins used in LumiHist, Reco-GEN matched plots and L2Res analysis (2nd step)
      //2018
      const int n_pt = 14;
      const double pt_bins[n_pt] = {//TODO check which code assumed the "min Bias" bin
      66,
      93,
      118,
      189,
      257,
      291,
      325,
      358,
      391,
      434,
      478,
      531,585,1000};

      const TString pt_range[n_pt]={
        "66",
        "93",
        "118",
        "189",
        "257",
        "291",
        "325",
        "358",
        "391",
        "434",
        "478",
        "531","585","1000"};

        const int n_pt_HF = 10;
        const double pt_bins_HF[n_pt_HF] = {
          66,
          73,
          93,
          116,
          142,
          210,
          279,
          379,
          1000,
          2000};

          const TString pt_range_HF[n_pt_HF]={"51","73","93","113","176","239","318","370","1000","2000"};

          //LumiBins for the Runs
          const int lumibins_BC[3] = {0, 1, 2};
          const int lumibins_B[2] = {0, 1};
          const int lumibins_C[1] = {2};


          //Other consts used in L2Res analysis
          const int nResponseBins = 100;
          const double Response_min = -1.2; //min of asymmetry hist
          const double Response_max = 1.2; //max of asymmetry hist
          const int n_etabarr=5; // needed for the normalization to 1 in the barrel.
          const int n_enough_entries = 100; //min N events per eta-pt bin to consider it for the L2Res derivation
          const double jet3pt_min=15;//jets with pt below were disregarded during alpha calculation for kFSR
