#ifndef _ExclusiveDiTau_h_
#define _ExclusiveDiTau_h_

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunExclusiveDiTau(TString filename,
                     TString outname,
                     TH1F *genPU,
                     TString era,
                     Bool_t debug=false,
                     std::string systVar = "");



#endif
