#ifndef _ExclusiveTop2020_h_
#define _ExclusiveTop2020_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunExclusiveTop2020(const TString filename,
                      TString outname,
                      TH1F *normH, 
                      TH1F *genPU,
                      TString era,
                      Bool_t debug=false,
                      Bool_t skimtree=false);

#endif
