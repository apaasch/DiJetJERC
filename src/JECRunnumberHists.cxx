#include "UHH2/BaconJets/include/JECRunnumberHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"

// #include "UHH2/BaconTrans/baconheaders/TJet.hh"
// #include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"
// #include "UHH2/BaconTrans/baconheaders/BaconAnaDefs.hh"
// #include "UHH2/BaconTrans/baconheaders/TVertex.hh"



#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <getopt.h>
using namespace std;
using namespace uhh2;
//using namespace baconhep;

JECRunnumberHists::JECRunnumberHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    book<TH1F>("runnr", "run number", 20001, 269999.5, 290000.5);
    
}

void JECRunnumberHists::fill(const uhh2::Event & ev){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  // Don't forget to always use the weight when filling.


  double weight = ev.weight;
  
  hist("runnr")->Fill(ev.run, weight);
    
  
}

JECRunnumberHists::~JECRunnumberHists(){}
