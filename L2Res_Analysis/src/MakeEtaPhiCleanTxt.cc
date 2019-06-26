//#include "../include/parameters.h"
#include "../../include/constants.h"
#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"
#include "../include/tdrstyle_mod15.h"

#include <TStyle.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TFile.h>
#include <TH2D.h>

#include <iostream>
#include <fstream>
using namespace std;

void CorrectionObject::MakeEtaPhiCleanTxt(){
  cout << "--------------- Starting MakeEtaPhiCleanTxT() ---------------" << endl << endl;
  cout <<"Asuming that the eventID, lumi and run numbers from events that would be cutted with the Eta-Phi cleaning are saved in the AnalysisTree while the eventID of all other events is saved as 0."<<endl<<endl;
  gStyle->SetOptStat(0);
  
  CorrectionObject::make_path(CorrectionObject::_outpath);
  ofstream _file;
  _file.open (CorrectionObject::_outpath+"EventsCuttedByEtaPhiCleaning.txt");
  _file << "run,\tlumiBlock,\tevent\n";

  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<int> eventID(myReader_DATA, "eventID");
  TTreeReaderValue<int> lumiSec(myReader_DATA, "lumi_sec");
  TTreeReaderValue<int> run(myReader_DATA, "run");
  while (myReader_DATA.Next()) {
    if(not *eventID) continue;
    _file << *run<<",\t"<<*lumiSec<<",\t\t"<<*eventID<<"\n";
  }

  cout<<"wrote file "<<CorrectionObject::_outpath<<"EventsCuttedByEtaPhiCleaning.txt"<<endl<<endl;
  _file.close();
}
    
    


