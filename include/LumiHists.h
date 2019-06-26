#include "UHH2/core/include/Hists.h"
#include <TString.h>
#include "../include/constants.h"
#include "UHH2/common/include/LuminosityHists.h"
#include <TProfile.h>
#pragma once

using namespace std;

/** \brief Create the "lumi plot", i.e. event yield vs. time in bins of equal integrated luminosity
 * 
 * Configuration from context:
 *  - "lumi_file": path to root file with luminosity information
 *  - "lumihists_lumi_per_bin": integrated luminosity per bin in the histogram (optional; default: 50.0)
 */


class LumiHists: public uhh2::Hists {
public:
    LumiHists(uhh2::Context & ctx,
                    const std::string & dirname,
                    const std::string & triggername = "",
                    bool do_inst_lumi_hist = false);

    virtual void fill(const uhh2::Event & ev) override;
    
private:
    
    // save the upper bin borders of those run/lumi numbers to
    // still include in the bin. Has size = nbins - 1, where
    // nbins is the number of bins in the lumi histogram
    std::vector<run_lumi> upper_binborders;

    std::map<run_lumi, double> rl2lumi;

    TH1D * hlumi;
    TH2D * hAsymLumi[37][13];
    TH2D * hBsymLumi[37][13];
    TProfile * pr_AsymLumi[37][13];
    TH1D * hinstlumi;
    TH1D * hinstlumi_ref;
    double lumi_per_bin;
    bool siTrg;
    bool trigger_fwd;
    const std::string triggername_;

  bool do_inst_lumi_hist_;
};
