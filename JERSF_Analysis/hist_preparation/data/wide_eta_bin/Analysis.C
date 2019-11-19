#include "MySelector.h"
#include <iostream>
#include <fstream>
#include <TChain.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TLatex.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <string>
#include <TCanvas.h>
#include <TFileCollection.h>
#include <TLorentzVector.h>
using namespace std;

int main ( int argc, char *argv[] ) {
  TString list_files;
  if ( argc > 3 ) { // argc should be 3 for correct execution
    std::cout<<"usage: "<< argv[0] <<" <filename>\n";
  } else {
    list_files = argv[1];
    std::ifstream path_to_file_with_list_of_inputfiles(list_files);
    // Always check to see if file opening succeeded
    if (!path_to_file_with_list_of_inputfiles.is_open()) std::cout<<"Could not open file\n";
    // the_file is closed implicitly here
  }

  bool isAK8 = list_files.Contains("AK8");
  std::string year = list_files.Contains("Autumn18")? "2018": (list_files.Contains("Fall17")? "2017": (list_files.Contains("Summer16")? "2016": "" )) ;
  std::cout << list_files << "\t" << isAK8 << "\t" << year << "\n";
  TString outdirname = "./";
  if (argc >= 3) {
    outdirname = argv[2];
  }

  std::cout << "outdirname: " << outdirname <<std::endl;

  MySelector *A = new MySelector(outdirname, year, isAK8);

  TChain* ch = new TChain("AnalysisTree");
  //TChain* ch = new TChain("AK4PFCHS/t");
  //TChain* c1 = new TChain("hlt/t");
  FILE *input;
  char filename[300];
  double TotalEvents = 0;
  TString path_to_file_with_list_of_inputfiles ( argv[1] );
  input = fopen( path_to_file_with_list_of_inputfiles, "r" );
  if (input != NULL) { // lets read each line and get the filename from it
    while (fscanf(input,"%s\n",filename) != EOF) {
      TFile* currentFile = TFile::Open(filename);
      std::cout << currentFile->GetName() << " " << ((TTree*)currentFile->Get("AnalysisTree"))->GetEntriesFast() << '\n';
      TotalEvents += ((TTree*)currentFile->Get("AnalysisTree"))->GetEntriesFast();
      ch->Add(filename);
      //ch->AddFriend(c1,filename);
    }
  }
  std::cout << "Total Events: " << TotalEvents << "\n";
  ch->Process(A);

}
