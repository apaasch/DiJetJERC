#pragma once

#ifndef  CONSTANTS_H
#define  CONSTANTS_H

const double eta_cut = 2.853-1e-5;
  
const int n_eta_full = 37;

const TString eta_range_full[n_eta_full] = {"-5.191","-3.839","-3.489","-3.139","-2.964","-2.853", "-2.65", "-2.5", "-2.322", "-2.172", "-1.93", "-1.653", "-1.479", "-1.305", "-1.044", "-0.783", "-0.522", "-0.261"," 0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};

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
constexpr static int n_eta = 19;
static std::vector<double>   eta_range  =  {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
static std::vector<double>   eta_range_mikko  = {0, 0.783, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
static std::vector<double>   alpha_range= {0., 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};

/* constexpr static int n_eta_RelVals = 14; */
/* static std::vector<double>   eta_range_RelVals  =  {0, 0.522, 1.044, 1.305, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 3.139, 3.489, 3.839, 5.191}; */
constexpr static int n_eta_RelVals = 19;
static std::vector<double>   eta_range_RelVals  =  {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

/** \brief Dijet event selection **/
// barrel region (|eta| < 1.3)
constexpr static float s_eta_barr = 1.3;
// two back-to-back leading jets (delta_phi(j1,j2) = min(|phi1 - phi2|, 2PI - |phi2 - phi1|) > 2.9)
constexpr static float s_delta_phi = 2.7;
// cut on the asymmetry for events with two jets  |(j2->pt - j1->pt /(j2->pt + j1->pt)| < 0.70
// constexpr static float s_asymm = 0.7;

constexpr static float s_asymm = 1.; //FIXME change back to 0.7?

 // cut on the asymmetry for events with two jets  |(j2->pt - j1->pt /(j2->pt + j1->pt)| < 0.2
/* constexpr static float s_asymm = 0.2; */
// relative third jet fraction pt_rel = 2*j3_pt/(j1_pt + j2_pt) < 0.2
constexpr static float s_pt_rel = 0.4;

/** \brief good Primary Vertex reconstruction **/
// more than four tracks
constexpr static float s_n_PvTracks = 4;
// PV is located within 24cm in z vertex
constexpr static float s_n_Pv_z = 24.0;
// PV is located within 2cm in xy direction from the nominal interaction point
constexpr static float s_n_Pv_xy = 2.0;


/** \brief The trigger thresholds of pt_ave **/


constexpr static int n_pt_bins = 9;


/* //from 94X max RunB to RunD */
 constexpr static float s_Pt_AveMC_cut_2017   =40;
 constexpr static float s_Pt_Ave40_cut_2017   =40 ;
 constexpr static float s_Pt_Ave60_cut_2017   =72 ;
 constexpr static float s_Pt_Ave80_cut_2017   =95 ;
 constexpr static float s_Pt_Ave140_cut_2017  =160;
 constexpr static float s_Pt_Ave200_cut_2017  =226;
 constexpr static float s_Pt_Ave260_cut_2017  =283;
 constexpr static float s_Pt_Ave320_cut_2017  =344;
 constexpr static float s_Pt_Ave400_cut_2017  =443;
 constexpr static float s_Pt_Ave450_cut_2017  =577;
 constexpr static float s_Pt_Ave500_cut_2017  =606;

//FIXME check the numbers below, are copy of central triggers atm, are they still?? does not look so...
  constexpr static float s_Pt_Ave40HF_cut_2017  =40 ;
  constexpr static float s_Pt_Ave60HF_cut_2017  =60 ;
  constexpr static float s_Pt_Ave80HF_cut_2017  =72 ;
 constexpr static float s_Pt_Ave140HF_cut_2017  =160;
 constexpr static float s_Pt_Ave200HF_cut_2017  =190;
 constexpr static float s_Pt_Ave260HF_cut_2017  =290;
 constexpr static float s_Pt_Ave320HF_cut_2017  =310;
 constexpr static float s_Pt_Ave400HF_cut_2017  =390;
 constexpr static float s_Pt_Ave450HF_cut_2017  =444;
 constexpr static float s_Pt_Ave500HF_cut_2017  =480;

//FixME: 2018 is copy of 2017 at the moment
constexpr static float s_Pt_AveMC_cut_2018   =40;
 constexpr static float s_Pt_Ave40_cut_2018   =40 ;
 constexpr static float s_Pt_Ave60_cut_2018   =72 ;
 constexpr static float s_Pt_Ave80_cut_2018   =95 ;
 constexpr static float s_Pt_Ave140_cut_2018  =160;
 constexpr static float s_Pt_Ave200_cut_2018  =226;
 constexpr static float s_Pt_Ave260_cut_2018  =283;
 constexpr static float s_Pt_Ave320_cut_2018  =344;
 constexpr static float s_Pt_Ave400_cut_2018  =443;
 constexpr static float s_Pt_Ave450_cut_2018  =577;
 constexpr static float s_Pt_Ave500_cut_2018  =606;

//FIXME check the numbers below, are copy of central triggers atm, are they still?? does not look so...
  constexpr static float s_Pt_Ave40HF_cut_2018  =40 ;
  constexpr static float s_Pt_Ave60HF_cut_2018  =60 ;
  constexpr static float s_Pt_Ave80HF_cut_2018  =72 ;
 constexpr static float s_Pt_Ave140HF_cut_2018  =160;
 constexpr static float s_Pt_Ave200HF_cut_2018  =190;
 constexpr static float s_Pt_Ave260HF_cut_2018  =290;
 constexpr static float s_Pt_Ave320HF_cut_2018  =310;
 constexpr static float s_Pt_Ave400HF_cut_2018  =390;
 constexpr static float s_Pt_Ave450HF_cut_2018  =444;
 constexpr static float s_Pt_Ave500HF_cut_2018  =480;


 constexpr static float Pt_AveMC_cut   =  51;
//from Di triggers 94X 17Nov2017
 constexpr static float d_Pt_Ave40_cut_2017   =  73;
 constexpr static float d_Pt_Ave60_cut_2017   =  85;
 constexpr static float d_Pt_Ave80_cut_2017   =  97;
 constexpr static float d_Pt_Ave140_cut_2017  = 179;
 constexpr static float d_Pt_Ave200_cut_2017  = 307;
 constexpr static float d_Pt_Ave260_cut_2017  = 370;
 constexpr static float d_Pt_Ave320_cut_2017  = 434;
 constexpr static float d_Pt_Ave400_cut_2017  = 520;
 constexpr static float d_Pt_Ave500_cut_2017  = 649;
// Dijet_HFJEC 2017 94X 17Nov2017
constexpr static float d_Pt_Ave60HF_cut_2017   = 73 ;
constexpr static float d_Pt_Ave80HF_cut_2017   = 93 ;
constexpr static float d_Pt_Ave100HF_cut_2017  = 113;
constexpr static float d_Pt_Ave160HF_cut_2017  = 176;
constexpr static float d_Pt_Ave220HF_cut_2017  = 239;
constexpr static float d_Pt_Ave300HF_cut_2017  = 318;



//from Di triggers 2018, RunABC, ReReco
//https://indico.cern.ch/event/801509/contributions/3331436/attachments/1801472/2938522/L2Res-Triggers-25Feb2019.pdf
 constexpr static float d_Pt_Ave40_cut_2018   =  66;
 constexpr static float d_Pt_Ave60_cut_2018   =  93;
 constexpr static float d_Pt_Ave80_cut_2018   =  118;
 constexpr static float d_Pt_Ave140_cut_2018  = 189;
 constexpr static float d_Pt_Ave200_cut_2018  = 257;
 constexpr static float d_Pt_Ave260_cut_2018  = 325;
 constexpr static float d_Pt_Ave320_cut_2018  = 391;
 constexpr static float d_Pt_Ave400_cut_2018  = 478;
 constexpr static float d_Pt_Ave500_cut_2018  = 585;
//Dijet_HFJEC 2018, RunABC, ReReco
constexpr static float d_Pt_Ave60HF_cut_2018   = 93 ;
constexpr static float d_Pt_Ave80HF_cut_2018   = 116 ;
constexpr static float d_Pt_Ave100HF_cut_2018  = 142;
constexpr static float d_Pt_Ave160HF_cut_2018  = 210;
constexpr static float d_Pt_Ave220HF_cut_2018  = 279;
constexpr static float d_Pt_Ave300HF_cut_2018  = 379;

// RunII rough bins, used in LumiHist only and Reco-GEN matched plots for pt binning
const int n_pt = 15;
const double pt_bins[n_pt] = {
  10,
  15,
  20,
  25,
  30,
  40,
  50,
  70,
  120,
  250,
  400,
  650,
  1000,1500,2000};

const TString pt_range[n_pt]={
  "10",
  "15",
  "20",
  "25",
  "30",
  "40",
  "50",
  "70",
  "120",
  "250",
  "400",
  "650",
  "1000","1500","2000"};


const int n_pt_HF = 8;
const double pt_bins_HF[n_pt_HF] = {
  70 ,
  95 ,
  120,
  150,
  210,
  280,
  1000,2000};

const TString pt_range_HF[n_pt_HF]={
  "70",
  "95",
  "120",
  "150",
  "210",
  "280",
  "1000","2000"};

/* 2017 */
/* const int n_pt = 12; */
/* const double pt_bins[n_pt] = { */
/*   51, */
/*   73, */
/*   85, */
/*   97, */
/*   179, */
/*   307, */
/*   370, */
/*   434, */
/*   520, */
/*   649,1000,2000}; */

/* const TString pt_range[n_pt]={ */
/*   "51", */
/*   "73", */
/*   "85", */
/*   "97", */
/*   "179", */
/*   "307", */
/*   "370", */
/*   "434", */
/*   "520", */
/*   "648","1000","2000"}; */


/* const int n_pt_HF = 8; */
/* const double pt_bins_HF[n_pt_HF] = { */
/*   73 , */
/*   93 , */
/*   113, */
/*   176, */
/*   239, */
/*   318, */
/*   1000,2000}; */

/* const TString pt_range_HF[n_pt_HF]={ */
/*   "73", */
/*   "93", */
/*   "113", */
/*   "176", */
/*   "239", */
/*   "318", */
/*   "1000","2000"}; */



/* //2017 */
/* //from Si triggers */
/* constexpr static int trg_vals_Si[10]  = {40,60,80,140,200,260,320,400,450,500}; */
/* constexpr static int trg_vals_Sifwd[7]  = {60,80,140,200,260,320,400}; */
/* constexpr static int trg_vals_Si_plusfwd[17]  = {40,60,80,140,200,260,320,400,450,500,60,80,140,200,260,320,400}; */
/* //from Di triggersint trg_val */
/* constexpr static int trg_vals_Di[9]  = {40,60,80,140,200,260,320,400,500}; */
/* constexpr static int trg_vals_HF[6]  = {60,80,100,160,220,300}; */
/* constexpr static int trg_vals_Di_plusHF[15]  = {40,60,80,140,200,260,320,400,500,60,80,100,160,220,300}; */


// Runnumbers for applying different corrections
// taken from PdmV, i.e 
// https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis
  constexpr static int s_runnr_B_2016  = 275376; //up to this one, including this one
  constexpr static int s_runnr_C_2016  = 276283; //up to this one, including this one
  constexpr static int s_runnr_D_2016 =  276811; //up to this one, including this one
  constexpr static int s_runnr_E_2016 =  277420; //up to this one, including this one
  constexpr static int s_runnr_F_2016 =  278801; //up to this one, including this one = Fearly
  constexpr static int s_runnr_G_2016 =  280385; //up to this one, including this one
  constexpr static int s_runnr_H_2016 =  284044; //up to this one, including this one

  constexpr static int s_runnr_B_2017  = 299329; //up to this one, including this one
  constexpr static int s_runnr_C_2017  = 302029; //up to this one, including this one
  constexpr static int s_runnr_D_2017 =  303434; //up to this one, including this one
  constexpr static int s_runnr_E_2017 =  304826; //up to this one, including this one
  constexpr static int s_runnr_F_2017  = 306462; //up to this one, including this one

  constexpr static int s_runnr_A_2018  = 316995; //up to this one, including this one
  constexpr static int s_runnr_B_2018  = 319310; //up to this one, including this one
  constexpr static int s_runnr_C_2018  = 320065; //up to this one, including this one
  constexpr static int s_runnr_D_2018 =  325175; //up to this one, including this one

/* //FixME: is it correct? */
/* constexpr static float s_lumi_HF_60 = 19.605; */
/* constexpr static float s_lumi_HF_80 = 28.773; */
/* constexpr static float s_lumi_HF_100 = 98.341; */
/* constexpr static float s_lumi_HF_160 = 561.485; */
/* constexpr static float s_lumi_HF_220 = 28983.619; */
/* constexpr static float s_lumi_HF_300 = 35672.846; */

/* // 2017 RunB PromptReco, not completly filled, not used if apply lumi weights is false */
/* //Lumi recorded by different triggers in pb-1 */
/* constexpr static float s_lumi_cent_40 = 0.09; */
/* constexpr static float s_lumi_cent_60 = 0.024; */
/* constexpr static float s_lumi_cent_80 = 0.; */
/* constexpr static float s_lumi_cent_140 = 27.769; */
/* constexpr static float s_lumi_cent_200 = 139.732; */
/* constexpr static float s_lumi_cent_260 = 521.863; */
/* constexpr static float s_lumi_cent_320 = 2965.774; */
/* constexpr static float s_lumi_cent_400 = 9021.331; */
/* constexpr static float s_lumi_cent_450 = 29280.311; */
/* constexpr static float s_lumi_cent_500 = 29280.311; */

#endif
