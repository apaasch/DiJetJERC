#pragma once

#include <TString.h>
#include <map>
#include <array>

//TODO get everything in a namespace, so that it can be checked against the constants from constants.h

using namespace std;

const double eta_cut = 2.853-1e-5;
//const double eta_cut = 100.; //central only

const int nResponseBins = 100;
const double Response_min = -1.2; //min of asymmetry hist
const double Response_max = 1.2; //max of asymmetry hist

const int n_etabarr=5; // needed for the normalization to 1 in the barrel. ???

const int n_enough_entries = 100;

//Alpha: 
const double alpha_cut = 0.3;
const TString s_alpha_cut = "0.3";
const TString s_alpha_cut_name = "03";
const double jet3pt_min=15;//jets with pt below were disregarded during processing
/* const double alpha_cut = 1.0; */
/* const TString s_alpha_cut = "1.0"; */
/* const TString s_alpha_cut_name = "10"; */



const int n_alpha = 11;
const TString alpha_range[n_alpha] = {"a005", "a010", "a015", "a020", "a025", "a030", "a035", "a040", "a045","a050","a055"};
const double alpha_bins[n_alpha] = {0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,  0.400, 0.450,0.50,0.55};

//Pt bins:


const int n_triggerDi = 9;
const int n_triggerDi_HF = 6;
const int n_triggerSi = 10;
const int n_triggerSi_HF = 9;
const int triggerValSi[n_triggerSi] = {40, 60, 80, 140, 200, 260, 320, 400, 450, 500};
const int triggerValDi[n_triggerDi] = {40, 60, 80, 140, 200, 260, 320, 400, 500};
const int triggerValDi_HF[n_triggerDi_HF]  = {60,80,100,160,220,300};
const int triggerValSi_HF[n_triggerSi_HF]  = {60,80,140,200,260,320,400,450,500};

// const int n_pt = 11;

//const int n_pt = 14;
const int n_pt_Si = 15;


// // //for Si derived trigger
// const TString pt_range[n_pt]= {"40",
// 			       "62",
// 			       "85",
// 			       "153",
// 			       "222",
// 			       "290",
// 			       "359",
// 			       "443",
// 			       "486",
// 			       "536",
// 			       "1000",
// 			       "2000"};




//for Di derived triggers
// const double pt_bins[n_pt] = {51, 51, 74, 96, 148, 214, 285, 346, 426, 525, 1000, 2000};
// const TString pt_range[n_pt]= {"51", "51", "74", "96", "148", "214", "285", "346", "426", "525", "1000", "2000"};

/*
//for Di derived triggers 94X 17Nov2017
//DiJet triggers -------------
const double pt_bins[n_pt] = {//TODO check which code assumed the "min Bias" bin
  51,
  73,
  85,
  97,
  140,
  179,
  240,
  307,
  370,
  434,
  520,
  649,1000,2000};

const TString pt_range[n_pt]={
  "51",
  "73",
  "85",
  "97",
  "140",
  "179",
  "240",
  "307",
  "370",
  "434",
  "520",
  "648","1000","2000"};
// // // //for Di derived triggers 94X 17Nov2017

const int n_pt_HF = 10;
const double pt_bins_HF[n_pt_HF] = {
 51,
 73,
 93,
 113,
 176,
 239,
 318,
 370,
 1000,
 2000};

const TString pt_range_HF[n_pt_HF]={"51","73","93","113","176","239","318","370","1000","2000"};
*/

//Dijet triggers 2018
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
// // // //for Di derived triggers 94X 17Nov2017

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



//SingleJet triggers
/* const double pt_bins[n_pt_Si] = { */
/* 40 , */
/* 72  , */
/* 95 , */
/* 130, */
/* 160 , */
/* 190, */
/* 226 , */
/* 250, */
/* 283 , */
/* 344 , */
/* 443 , */
/* 577 , */
/* 606 ,1000,2000}; */

/* const TString pt_range[n_pt_Si]={ */
/*   "40" , */
/*   "72" , */
/*   "95" , */
/*   "130", */
/*   "160", */
/*   "190", */
/*   "226", */
/*   "250", */
/*   "283", */
/*   "344", */
/*   "443", */
/*   "577", */
/*   "606","1000","2000"}; */
/* // [END]SingleJet triggers? ------------- */


// const double pt_bins_Si[n_pt_Si] = {
//  40 ,
//  65 ,
//  74  ,
//  100 ,
//  165 ,
//  226 ,
//  285 ,
//  346 ,
//  445 ,
//  566 ,
//  608 ,1000,2000};





// const double pt_bins_HF[n_pt] = {60,80,100,160,220,300, 1000, 2000};
// const TString pt_range_HF[n_pt]= {"60","80","100","160","220","300", "1000", "2000"};

// di from back to back with trg matching
// 	 50.02 
// 	 70.91 
// 	 91.80 
// 	 161.66 
// 	  240
// 	 278.90 
// 	 343.21 
// 	 425.20 
// 	 523.62 

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


const double eta_bins2[19]     = { 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

const int n_alpha_common = 4;
const TString alpha_range_common[n_alpha_common] = {"a10","a15", "a20", "a30"};//for Global fit we produce result only in few alpha values
const double alpha_bins_common[n_alpha_common] = {0.100, 0.150, 0.200, 0.300};

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


//LumiBins for the Runs
const int lumibins_BC[3] = {0, 1, 2};
const int lumibins_B[2] = {0, 1};
const int lumibins_C[1] = {2};
