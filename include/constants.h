#pragma once

#include <TString.h>
#include <map>
#include <array>
#include <TMath.h>

// barrel region (|eta| < 1.131)
constexpr static double s_eta_barr = 1.131;
// constexpr static double s_eta_barr = 1.044; // changed for shift 4
// HF region (|eta| > 2.853)
const double pi = TMath::Pi();
const double eta_cut = 2.853;
// const double eta_cut_barrel = 1.566;
// two back-to-back leading jets (delta_phi(j1,j2) = min(|phi1 - phi2|, 2PI - |phi2 - phi1|) > 2.7)
constexpr static double s_delta_phi = 2.7;

constexpr static double jet_threshold      = 15;
constexpr static double jet_threshold_min  = 10;

const double phi_bins[13] = { 0*pi/4, 1*pi/4, 2*pi/4, 3*pi/4, 4*pi/4, 5*pi/4, 6*pi/4, 7*pi/4, 8*pi/4, 9*pi/4, 10*pi/4, 11*pi/4, 12*pi/4};

const int n_eta_bins_JER = 15;
const int n_eta_bins_L2R = 19;
const int n_eta_bins_JERC = 20;
const int n_eta_bins_common = 18;
const int n_eta_bins_narrow = 22;
const int n_eta_bins_simple = 10;
const int n_eta_bins_calo = 42;
const double eta_bins_JER[n_eta_bins_JER]       = { 0.000,        0.522, 0.783,        1.131, 1.305,               1.740, 1.930, 2.043,        2.322, 2.500, 2.650, 2.853, 2.964, 3.139,               5.191};
const double eta_bins_L2R[n_eta_bins_L2R]       = { 0.000, 0.261, 0.522, 0.783, 1.044,        1.305, 1.479, 1.653,        1.930,        2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
const double eta_bins_JERC[n_eta_bins_JERC]     = { 0.000, 0.261, 0.522, 0.783, 1.044,        1.305,        1.566, 1.740, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191}; // old common binning
const double eta_bins_common[n_eta_bins_common] = { 0.000, 0.261, 0.522, 0.783, 1.044,        1.305,        1.566, 1.740, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139,               5.191}; // Combine last three bins for statistics
const double eta_bins_narrow[n_eta_bins_narrow] = { 0.000, 0.261, 0.522, 0.783, 1.044, 1.131, 1.305, 1.479, 1.653, 1.740, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
const double eta_bins_simple[n_eta_bins_simple] = { 0.000,                                    1.305,                                                  2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
const double eta_bins_calo[n_eta_bins_calo]     = { 0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

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
  606,1000,2000
};

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
  "606","1000","2000"
};

/* //SingleJet triggers highest End------------- */

//19 bin edges, 18 actual bins
//constexpr static int n_eta = 19;
//static std::vector<double>   eta_range  =  {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
static std::vector<double>   eta_range_mikko  = {0, 0.783, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
//Eta bins for GlobalFit files
// const int n_eta_bins_common = 20;
// const double eta_bins_common[n_eta_bins_common] = { 0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.566, 1.740, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
const int n_eta_common = 19;
const double eta_common_bins[n_eta_common] = {0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
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
  "5.191"
};

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
  {"SingleJet_central", {"trigger40", "trigger60", "trigger80", "trigger140", "trigger200", "trigger260", "trigger320", "trigger400", "trigger450", "trigger500"}},
  // {"SingleJet_central", {"trigger40", "trigger60", "trigger80", "trigger140", "trigger200", "trigger260", "trigger320", "trigger400", "trigger450", "trigger500", "trigger550"}},// 550 not present in 2016. Hack in the code.
  {"DiJet_forward",     {"trigger60_HFJEC", "trigger80_HFJEC",  "trigger100_HFJEC", "trigger160_HFJEC", "trigger220_HFJEC", "trigger300_HFJEC"}},
  {"SingleJet_forward", {"trigger40_HFJEC", "trigger60_HFJEC",  "trigger80_HFJEC",  "trigger140_HFJEC", "trigger200_HFJEC", "trigger260_HFJEC", "trigger320_HFJEC", "trigger400_HFJEC"}}, // last one added do to Config , "trigger500_HFJEC"
};


const std::map<std::string, std::vector<double> > pt_trigger_thr = {
  // Mix between 2017 and 2018 to be fixed!
  {"DiJet_central_Legacy",                    { 73, 85, 97, 179, 307, 370, 434, 520, 649 }},
  {"DiJet_forward_Legacy",                    { 93, 116, 142, 210, 279, 379 }},
  {"DiJet_central_Legacy_ptbins",             { 73, 85, 97, 179, 307, 370, 434, 520, 649 }},
  {"DiJet_forward_Legacy_ptbins",             { 93, 116, 142, 210, 279, 379 }},

  // UL16preVFP AK4
  {"DiJet_central_UL16preVFP",                { 59, 85, 104, 170, 236, 302, 370, 460, 575 }}, // from Jindrich
  {"DiJet_central_UL16preVFP_ptbins",         { 59, 75,  85,  98, 104, 128, 145, 170, 190, 220, 236, 280, 302, 325, 350, 370, 400, 430, 460, 500, 530, 555, 575, 644, 730, 790, 840, 920, 1020, 1120, 1400}}, // 30+1 bins
  {"DiJet_central_UL16preVFP_ptbins_default", { 59, 85, 104, 170, 236, 302, 370, 460, 575 }},
  {"DiJet_central_UL16preVFP_ptbins_fine_v1", { 59, 75,  85,  98, 104, 128, 145, 170, 190, 220, 236, 280, 302, 325, 350, 370, 400, 430, 460, 500, 530, 555, 575, 644, 730, 790, 840, 920, 1020, 1120, 1400}}, // 30+1 bins
  {"DiJet_forward_UL16preVFP",                { 86, 110, 132, 204, 279, 373 }}, // from Jindrich
  {"DiJet_forward_UL16preVFP_ptbins",         { 86, 93, 100, 110, 122, 132, 150, 176, 180, 204, 220, 239, 279, 373 }}, // 14+1 bins
  {"DiJet_forward_UL16preVFP_ptbins_default", { 86, 110, 132, 204, 279, 373 }}, // from Jindrich
  {"DiJet_forward_UL16preVFP_ptbins_fine_v1", { 86, 93, 100, 110, 122, 132, 150, 176, 180, 204, 220, 239, 279, 373 }}, // 14+1 bins
  // UL16preVFP AK8
  {"SingleJet_central_AK8_UL16preVFP",        { 77, 96, 117, 190, 256, 321, 386, 473, 526, 581 }}, // Copied from UL18; Also for SJ AK4 for RunB
  {"SingleJet_central_AK8_UL16preVFP_ptbins", { 77, 96, 117, 190, 256, 321, 386, 473, 526, 581 }}, // Copied from UL18; Also for SJ AK4 for RunB
  {"SingleJet_forward_AK8_UL16preVFP",        { 65, 103, 115, 179, 252, 317, 410, 519 }}, // Copied from UL18; Also for SJ AK4 for RunB
  {"SingleJet_forward_AK8_UL16preVFP_ptbins", { 65, 103, 115, 179, 252, 317, 410, 519 }}, // Copied from UL18; Also for SJ AK4 for RunB

  // UL16postVFP AK4
  {"DiJet_central_UL16postVFP",                { 59, 85, 104, 170, 236, 302, 370, 460, 575 }}, // from Jindrich
  {"DiJet_central_UL16postVFP_ptbins",         { 59, 75,  85,  98, 104, 128, 145, 170, 190, 220, 236, 280, 302, 325, 350, 370, 400, 430, 460, 500, 530, 555, 575, 644, 730, 790, 840, 920, 1020, 1120, 1400}}, // 30+1 bins
  {"DiJet_central_UL16postVFP_ptbins_default", { 59, 85, 104, 170, 236, 302, 370, 460, 575 }}, // from Jindrich
  {"DiJet_central_UL16postVFP_ptbins_fine_v1", { 59, 75,  85,  98, 104, 128, 145, 170, 190, 220, 236, 280, 302, 325, 350, 370, 400, 430, 460, 500, 530, 555, 575, 644, 730, 790, 840, 920, 1020, 1120, 1400}}, // 30+1 bins
  {"DiJet_forward_UL16postVFP",                { 86, 110, 132, 204, 279, 373 }}, // from Jindrich
  {"DiJet_forward_UL16postVFP_ptbins",         { 86, 93, 100, 110, 122, 132, 150, 176, 180, 204, 220, 239, 279, 373 }}, // 14+1 bins
  {"DiJet_forward_UL16postVFP_ptbins_default", { 86, 110, 132, 204, 279, 373 }}, // from Jindrich
  {"DiJet_forward_UL16postVFP_ptbins_fine_v1", { 86, 93, 100, 110, 122, 132, 150, 176, 180, 204, 220, 239, 279, 373 }}, // 14+1 bins
  // UL16 AK8 COPY FROM 2018
  {"SingleJet_central_AK8_UL16postVFP",       { 77, 96, 117, 190, 256, 321, 386, 473, 526, 581 }}, // Copied from UL18
  {"SingleJet_central_AK8_UL16postVFP_ptbins",{ 77, 96, 117, 190, 256, 321, 386, 473, 526, 581 }}, // Copied from UL18
  {"SingleJet_forward_AK8_UL16postVFP",       { 65, 103, 115, 179, 252, 317, 410, 519 }}, // Copied from UL18
  {"SingleJet_forward_AK8_UL16postVFP_ptbins",{ 65, 103, 115, 179, 252, 317, 410, 519 }}, // Copied from UL18

  // 2017 AK4
  {"DiJet_central_2017",                { 73, 85, 97, 179, 307, 370, 434, 520, 649 }}, //for Di triggers 94X 17Nov2017
  {"DiJet_forward_2017",                { 73, 93, 113, 176, 239, 318 }}, // for Dijet_HFJEC 2017 94X 17Nov2017
  {"SingleJet_central_2017",            { 40, 72,  95, 160, 226, 283, 344, 443, 577, 606 }}, //for Single triggers 94X 17Nov2017
  {"SingleJet_forward_2017",            { 73, 93, 113, 176, 239, 318 }}, // for Singlejet_HFJEC 2017 94X 17Nov2017 //TODO not used
  {"SingleJet_central_AK8_2017",        { 73, 90, 115, 181, 251, 312, 378, 457, 519, 566 }}, //for Single triggers 94X 17Nov2017
  // UL17 AK4
  {"DiJet_central_UL17",                { 70, 87, 111, 180, 247, 310, 373, 457, 562 }}, // SJ AK4 UL17 (Jindrich); No DJ for B,C
  {"DiJet_central_UL17_ptbins",         { 70, 73, 87, 93, 111, 113, 145, 176, 180, 239, 247, 280, 310, 318, 350, 370, 400, 430, 460, 500, 530, 555, 575, 644, 730, 790, 840, 920, 1020, 1120, 1400}}, // 30+1 bins
  {"DiJet_central_UL17_ptbins_default", { 70, 73, 87, 93, 111, 113, 176, 180, 239, 247, 310, 318, 373, 457, 562 }},// Combined DJ fwd and SJ central AK4
  {"DiJet_central_UL17_ptbins_fine_v1", { 70, 73, 87, 93, 111, 113, 145, 176, 180, 239, 247, 280, 310, 318, 350, 370, 400, 430, 460, 500, 530, 555, 575, 644, 730, 790, 840, 920, 1020, 1120, 1400}}, // 30+1 bins
  {"DiJet_forward_UL17",                { 73, 93, 113, 176, 239, 318 }}, // Copied from 2017
  {"DiJet_forward_UL17_ptbins",         { 70, 73, 87, 93, 111, 113, 140, 176, 180, 205, 239, 247, 310, 318}}, // 14+1 bins
  {"DiJet_forward_UL17_ptbins_fine_v1", { 70, 73, 87, 93, 111, 113, 140, 176, 180, 205, 239, 247, 310, 318}}, // 14+1 bins
  {"DiJet_forward_UL17_ptbins_default", { 70, 73, 87, 93, 111, 113, 176, 180, 239, 247, 310, 318}}, // Combined DJ fwd and SJ central AK4
  {"SingleJet_central_UL17",            { 70, 87, 111, 180, 247, 310, 373, 457, 562 }}, // From Jindrich / Skip trigger 450 since no DiJet trigger
  {"SingleJet_central_UL17_ptbins",     { 70, 73, 87, 93, 111, 113, 176, 180, 239, 247, 310, 318, 373, 457, 562 }},// Combined DJ fwd and SJ central AK4
  // UL17 AK8
  {"SingleJet_central_AK8_UL17",        { 77, 96, 117, 190, 256, 321, 386, 473, 526, 581 }}, // Copied from UL18
  {"SingleJet_central_AK8_UL17_ptbins", { 65, 77, 96, 103, 115, 117, 179, 190, 252, 256, 317, 321, 410, 519, 526, 581 }},
  {"SingleJet_forward_AK8_UL17",        { 65, 103, 115, 179, 252, 317, 410, 519 }}, // Copied from UL18
  {"SingleJet_forward_AK8_UL17_ptbins", { 65, 77, 96, 103, 115, 117, 179, 190, 252, 256, 317, 321, 410, 519 }},

  // 2018 AK4
  {"DiJet_central_2018",                { 66, 93, 118, 189, 257, 325, 391, 478, 585 }}, //for Di triggers 2018, RunABC, ReReco https://indico.cern.ch/event/801509/contributions/3331436/attachments/1801472/2938522/L2Res-Triggers-25Feb2019.pdf
  {"DiJet_central_2018_ptbins",         { 66, 93, 118, 189, 257, 291, 325, 358, 391, 434, 478, 531, 585}},
  {"DiJet_forward_2018",                { 93, 116, 142, 210, 279, 379 }},
  {"DiJet_forward_2018_ptbins",         { 93, 116, 142, 210, 279, 379 }},
  // 2018 AK8
  {"SingleJet_central_AK8_2018",        { 78, 96, 119, 193, 262, 328, 393, 481, 534, 588 }}, //HLT AK8PFJet* //https://indico.desy.de/indico/event/24350/contribution/0/material/slides/0.pdf
  {"SingleJet_central_AK8_2018_ptbins", { 78, 96, 119, 193, 262, 328, 393, 481, 534, 588 }},
  {"SingleJet_forward_AK8_2018",        { 62, 95, 110, 182, 260, 339, 420, 508 }}, // HLT AK8PFJetFwd* //https://indico.desy.de/indico/event/24423/contribution/1/material/slides/1.pdf
  {"SingleJet_forward_AK8_2018_ptbins", { 62, 95, 110, 182, 260, 339, 420, 508 }},
  // UL18 AK4
  {"DiJet_central_UL18",                { 66, 93, 118, 189, 257, 325, 391, 478, 585 }}, //for Di triggers 2018, RunABC, ReReco https://indico.cern.ch/event/801509/contributions/3331436/attachments/1801472/2938522/L2Res-Triggers-25Feb2019.pdf
  {"DiJet_central_UL18_ptbins",         { 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1400, 1500, 1600, 2000, 2500}}, // 30+1 bins
  {"DiJet_central_UL18_ptbins_default", { 66, 93, 118, 189, 257, 325, 391, 478, 585 }}, //for Di triggers 2018, RunABC, ReReco https://indico.cern.ch/event/801509/contributions/3331436/attachments/1801472/2938522/L2Res-Triggers-25Feb2019.pdf
  {"DiJet_central_UL18_ptbins_quick",   { 66 }}, //for Di triggers 2018, RunABC, ReReco https://indico.cern.ch/event/801509/contributions/3331436/attachments/1801472/2938522/L2Res-Triggers-25Feb2019.pdf
  {"DiJet_central_UL18_ptbins_fine",    { 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120}}, // 30 bins
  {"DiJet_central_UL18_ptbins_fine_v1", { 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120, 1400}}, // 30+1 bins
  {"DiJet_central_UL18_ptbins_fine_v2", { 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120, 1390, 1670}}, // 31+1 bins
  {"DiJet_central_UL18_ptbins_fine_v3", { 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1400, 1500, 1600, 2000, 2500}}, // 33+1 bins
  {"DiJet_central_UL18_ptbins_fine_v4", { 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120, 1240, 1390, 1670}}, // 32+1 bins
  {"DiJet_central_UL18_ptbins_fine_v5", { 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120, 1240, 1400, 1500, 1600, 2000, 2500}}, // 35+1 bins
  {"DiJet_central_UL18_ptbins_fine_v6", { 66, 71,  77,  93,  98, 106, 118, 128, 145, 160, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478,  531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120, 1240, 1400, 1500, 1600, 1750, 2000, 2500}}, // 37+1 bins
  // {"DiJet_central_UL18_ptbins_alpha",   { 15, 22, 28, 32, 38, 42, 48, 66, 93, 118, 189, 257, 325, 391, 478, 585, 644, 730, 790, 840, 920, 1020, 1120, 1120, 1400}}, // 24+1 bins
  {"DiJet_central_UL18_ptbins_alpha",   { 15, 66, 71,  77,  93,  98, 106, 118, 128, 145, 189, 203, 223, 257, 291, 325, 358, 391, 434, 478, 531, 546, 563, 585, 644, 730, 790, 840, 920, 1020, 1120, 1400}}, // 31+1 bins
  {"DiJet_forward_UL18",                { 93, 116, 142, 210, 279, 379 }},
  {"DiJet_forward_UL18_ptbins",         { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }}, // 14+1 bins
  {"DiJet_forward_UL18_ptbins_default", { 93, 116, 142, 210, 279, 379 }},
  {"DiJet_forward_UL18_ptbins_quick",   { 93 }},
  {"DiJet_forward_UL18_ptbins_fine",    { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }},
  {"DiJet_forward_UL18_ptbins_fine_v1", { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }},
  {"DiJet_forward_UL18_ptbins_fine_v2", { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }}, // only change in central region
  {"DiJet_forward_UL18_ptbins_fine_v3", { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }}, // only change in central region
  {"DiJet_forward_UL18_ptbins_fine_v4", { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }}, // only change in central region
  {"DiJet_forward_UL18_ptbins_fine_v5", { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }}, // only change in central region
  {"DiJet_forward_UL18_ptbins_fine_v6", { 93, 96,  99, 102, 106, 111, 116, 119, 122, 125, 130, 135, 142, 146, 154, 160, 172, 185, 210, 220, 240, 279, 379 }}, // only change in central region
  {"DiJet_forward_UL18_ptbins_alpha",   { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }}, // only change in forward region
  {"DiJet_forward_UL18_ptbins_common",  { 93, 99, 106, 116, 122, 130, 142, 154, 172, 210, 220, 240, 279, 379 }}, // same as fine_v1, but keep for consitency
  {"DiJet_central_UL18_Jindrich",       { 69, 95, 120, 193, 264, 326, 399, 474, 595 }}, // https://indico.cern.ch/event/1165351/contributions/4893655/attachments/2451202/4200382/L2Res_05_2022-2.pdf
  {"DiJet_forward_UL18_Jindrich",       { 95, 120, 144, 215, 282, 388 }}, // https://indico.cern.ch/event/1165351/contributions/4893655/attachments/2451202/4200382/L2Res_05_2022-2.pdf
  // UL18 AK8
  {"SingleJet_central_AK8_UL18",                { 77, 96, 117, 190, 256, 321, 386, 473, 526, 581 }}, // https://indico.cern.ch/event/983310/contributions/4144228/attachments/2159426/3643048/L2Res_09_12_2020.pdf
  {"SingleJet_central_AK8_UL18_ptbins",         { 77, 96, 117, 190, 256, 321, 386, 473, 526, 581 }},
  {"SingleJet_central_AK8_UL18_ptbins_default", { 77, 96, 117, 190, 256, 321, 386, 473, 526, 581 }},
  {"SingleJet_central_AK8_UL18_ptbins_fine_v1", { 77, 80,  96, 117, 150, 190, 220, 256, 280, 321, 386, 420, 473, 526, 581, 1000 }},
  {"SingleJet_forward_AK8_UL18",                { 65, 103, 115, 179, 252, 317, 410, 519 }}, // https://indico.cern.ch/event/983310/contributions/4144228/attachments/2159426/3643048/L2Res_09_12_2020.pdf
  {"SingleJet_forward_AK8_UL18_ptbins",         { 65, 103, 115, 179, 252, 317, 410, 519 }},
  {"SingleJet_forward_AK8_UL18_ptbins_default", { 65, 103, 115, 179, 252, 317, 410, 519 }},
  {"SingleJet_forward_AK8_UL18_ptbins_fine_v1", { 65,  80, 103, 115, 140, 179, 210, 252, 317, 410, 519, 1000 }},
  // 2022
  {"DiJet_central_2022preEE",                { 66, 83, 108, 177, 246, 313, 378, 464, 572 }},
  {"DiJet_central_2022preEE_ptbins_default", { 66, 75, 83, 95, 108, 128, 145, 177, 195, 225, 246, 285, 313, 340, 378, 420, 464, 490, 530, 572, 644, 730, 790, 840, 920, 1020, 1120, 1400 }},
  {"DiJet_forward_2022preEE",                { 79, 99, 122, 194, 273, 365 }},
  {"DiJet_forward_2022preEE_ptbins_default", { 79, 99, 109, 122, 140, 150, 194, 220, 250, 273, 365 }},

  {"DiJet_central_2022postEE",                { 66, 83, 108, 177, 246, 313, 378, 464, 572 }},
  {"DiJet_central_2022postEE_ptbins_default", { 66, 75, 83, 95, 108, 128, 145, 177, 195, 225, 246, 285, 313, 340, 378, 420, 464, 490, 530, 572, 644, 730, 790, 840, 920, 1020, 1120, 1400 }},
  {"DiJet_central_2022postEE_ptbins_fine_v1", { 66, 75, 83, 95, 108, 135, 177, 195, 225, 246, 285, 313, 350, 378, 420, 464, 490, 530, 572, 620, 700, 750, 800, 850, 900, 1000 }}, // 25 bins
  {"DiJet_central_2022postEE_ptbins_fine_v2", { 66, 75, 83, 95, 108, 135, 177, 195, 225, 246, 285, 313, 350, 378, 420, 464, 490, 530, 572, 620, 700, 800, 900, 1000 }}, // 25 bins
  {"DiJet_central_2022postEE_ptbins_fine_v3", { 66, 83, 108, 140, 177, 210, 246, 275, 313, 378, 434, 464, 500, 525, 550, 572, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900 }}, // 25 bins
  {"DiJet_central_2022postEE_ptbins_fine_v4", { 66, 83, 108, 140, 177, 210, 246, 275, 313, 378, 434, 464, 500, 525, 550, 572, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 1000, 1200 }}, // 25 bins
  {"DiJet_forward_2022postEE",                { 79, 99, 122, 194, 273, 365 }},
  {"DiJet_forward_2022postEE_ptbins_default", { 79, 99, 109, 122, 140, 150, 194, 220, 250, 273, 365 }},
  {"DiJet_forward_2022postEE_ptbins_fine_v1", { 79, 99, 122, 194, 273, 365 }},
  {"DiJet_forward_2022postEE_ptbins_fine_v2", { 79, 99, 122, 194, 273, 365 }},
  {"DiJet_forward_2022postEE_ptbins_fine_v3", { 79, 99, 122, 194, 240, 273, 300, 365, 400, 500 }},
  {"DiJet_forward_2022postEE_ptbins_fine_v4", { 79, 99, 122, 194, 240, 273, 300, 365, 400, 500 }},
  // 2023 pre BPix
  {"DiJet_central_2023preBPix",                { 66, 83, 108, 177, 246, 313, 378, 464, 572 }},
  {"DiJet_central_2023preBPix_ptbins_default", { 66, 75, 83, 95, 108, 128, 145, 177, 195, 225, 246, 285, 313, 340, 378, 420, 464, 490, 530, 572, 644, 730, 790, 840, 920, 1020, 1120, 1400 }},
                                              //  { 10., 15., 21., 28., 37., 49., 59., 86.,110., 132.,170.,204.,236.,279.,302.,373.,460.,575., 638.,737.,846.,  967., 1101., 1248., 1410., 1588., 1784., 2000., 2238., 2500., 2787., 3103., 3450., 4037., 5220.]
  {"DiJet_central_2023preBPix_ptbins_finev1",  { 66, 75, 83, 95, 108, 140, 177, 210, 246, 285, 313, 340, 378, 420, 464, 500, 572, 600, 700, 800, 900, 1050, 1150, 1500 }},
  {"DiJet_forward_2023preBPix",                { 79, 99, 122, 194, 273, 365 }},
  {"DiJet_forward_2023preBPix_ptbins_default", { 79, 99, 109, 122, 140, 150, 194, 220, 250, 273, 365 }},
  {"DiJet_forward_2023preBPix_ptbins_finev1",  { 79, 99, 109, 122, 140, 150, 194, 220, 250, 273, 365 }},
  // 2023 post BPix
  {"DiJet_central_2023postBPix",                { 66, 83, 108, 177, 246, 313, 378, 464, 572 }},
  {"DiJet_central_2023postBPix_ptbins_default", { 66, 75, 83, 95, 108, 128, 145, 177, 195, 225, 246, 285, 313, 340, 378, 420, 464, 490, 530, 572, 644, 730, 790, 840, 920, 1020, 1120, 1400 }},
  {"DiJet_central_2023postBPix_ptbins_finev1",  { 66, 75, 83, 95, 108, 140, 177, 210, 246, 285, 313, 340, 378, 420, 464, 500, 572, 600, 700, 800, 900, 1050, 1150, 1500 }},
  {"DiJet_forward_2023postBPix",                { 79, 99, 122, 194, 273, 365 }},
  {"DiJet_forward_2023postBPix_ptbins_default", { 79, 99, 109, 122, 140, 150, 194, 220, 250, 273, 365 }},
  {"DiJet_forward_2023postBPix_ptbins_finev1",  { 79, 99, 109, 122, 140, 150, 194, 220, 250, 273, 365 }},
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
