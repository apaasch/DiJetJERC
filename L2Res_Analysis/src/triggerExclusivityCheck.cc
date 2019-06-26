#include "../include/parameters.h"
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

using namespace std;

void CorrectionObject::triggerExclusivityCheck(){
  cout << "--------------- Starting triggerExclusivityCheck() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/fullAsym/");

  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<int> trg40(myReader_DATA, "trigger40");
  TTreeReaderValue<int> trg60(myReader_DATA, "trigger60");
  TTreeReaderValue<int> trg80(myReader_DATA, "trigger80");
  TTreeReaderValue<int> trg140(myReader_DATA, "trigger140");
  TTreeReaderValue<int> trg200(myReader_DATA, "trigger200");
  TTreeReaderValue<int> trg260(myReader_DATA, "trigger260");
  TTreeReaderValue<int> trg320(myReader_DATA, "trigger320");
  TTreeReaderValue<int> trg400(myReader_DATA, "trigger400");
  TTreeReaderValue<int> trg450(myReader_DATA, "trigger450");
  TTreeReaderValue<int> trg500(myReader_DATA, "trigger500");
    
  int myCount = 0;
  bool allExclusive = true;
  while (myReader_DATA.Next()) {
    bool exclusive = true;
    exclusive = (*trg40)^(*trg60)^(*trg80)^(*trg140)^(*trg200)^(*trg260)^(*trg320)^(*trg400)^(*trg450)^(*trg500);
    if(!exclusive){
      allExclusive=false;
      cout<<"event nr. "<< "is not exclusive"<<endl;
    }
    myCount++;   
  }

  //DEBUG
  std::cout<<"\ncount data "<<myCount<<std::endl;
  if(allExclusive) cout<<"All trigger are really exclusive."<<endl;
  

}
    
    


