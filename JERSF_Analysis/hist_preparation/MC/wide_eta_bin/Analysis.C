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

int main ( int argc, char *argv[] )
//int main()
{
  if ( argc > 3 ) // argc should be 2 for correct execution
  // We print argv[0] assuming it is the program name
  std::cout<<"usage: "<< argv[0] <<" <filename>\n";
  else {
    // We assume argv[1] is a filename to open
    std::ifstream path_to_file_with_list_of_inputfiles ( argv[1] );
    // Always check to see if file opening succeeded
    if ( !path_to_file_with_list_of_inputfiles.is_open() )
    std::cout<<"Could not open file\n";
    else {
      char x;
      // path_to_file_with_list_of_inputfiles.get ( x ) returns false if the end of the file
      //  is reached or an error occurs
      int i;
      while ( path_to_file_with_list_of_inputfiles.get ( x ) && i < 3 ){
        std::cout<< x;
        i++;
      }
      std::cout << "..." << std::endl;
    }
    // the_file is closed implicitly here
  }

  TString outdirname = "./";
  if (argc >= 3) {
    outdirname = argv[2];
  }

  std::cout << "outdirname: " << outdirname <<std::endl;

  MySelector *A = new MySelector(outdirname);

  TChain* ch = new TChain("AnalysisTree");
  //TChain* ch = new TChain("AK4PFCHS/t");
  //TChain* c1 = new TChain("hlt/t");
  FILE *input;
  char filename[300];



  TString path_to_file_with_list_of_inputfiles ( argv[1] );
  input = fopen( path_to_file_with_list_of_inputfiles, "r" );
  if (input != NULL)
  { // lets read each line and get the filename from it
    while (fscanf(input,"%s\n",filename) != EOF)
    {
      printf("%s\n",filename);
      ch->Add(filename);
      //ch->AddFriend(c1,filename);
    }
  }

  ch->Process(A);

}
