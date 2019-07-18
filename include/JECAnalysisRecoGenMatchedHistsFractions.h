#pragma once
//#include <UHH2/BaconJets/include/HistsBASE.h>
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "TClonesArray.h"
//#include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"
//#include "UHH2/BaconTrans/baconheaders/BaconAnaDefs.hh"
/**  \brief Example class for booking and filling histograms
*
* NOTE: This class uses the 'hist' method to retrieve histograms.
* This requires a string lookup and is therefore slow if you have
* many histograms. Therefore, it is recommended to use histogram
* pointers as member data instead, like in 'common/include/ElectronHists.h'.
*/
class JECAnalysisRecoGenMatchedHistsFractions: public uhh2::Hists {
  //class JECAnalysisRecoGenMatchedHistsFractions: public HistsBASE {
public:
  // use the same constructor arguments as Hists for forwarding:
  JECAnalysisRecoGenMatchedHistsFractions(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  void fill(const uhh2::Event & ev, const int genjet_id, const int recjet_id);
  virtual ~JECAnalysisRecoGenMatchedHistsFractions();

};
