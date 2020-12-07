#ifndef _LowMu2020_h_
#define _LowMu2020_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunLowMu2020(const TString filename,
                      TString outname,
                      TH1F *genPU,
                      TString era,
                      Bool_t debug=false);

#endif
