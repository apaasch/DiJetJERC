#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "TClonesArray.h"

class JECRunnumberHists: public uhh2::Hists {

public:
  // use the same constructor arguments as Hists for forwarding:
  JECRunnumberHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~JECRunnumberHists();


};
