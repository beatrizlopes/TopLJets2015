#ifndef _ExclusiveTop_h_
#define _ExclusiveTop_h_

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunExclusiveTop(TString filename,
                     TString outname,
                     Int_t channelSelection,
                     Int_t chargeSelection,
                     TH1F *normH,
                     TH1F *genPU,
                     TString era,
                     Bool_t debug=false,
                     std::string systVar = "");

void chistar(int& npar, double *deriv, double& f, double *par, int flag);
void analyticMattFit (TLorentzVector* bJet_Had, TLorentzVector* bJet_Lep, TLorentzVector* Lep, TLorentzVector* nu, TLorentzVector* lightJet0, TLorentzVector* lightJet1, TMatrixD covarianceMatrix, double* chiSquare);

#endif
