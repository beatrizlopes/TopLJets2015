#ifndef _RunFwdProtons_h_
#define _RunFwdProtons_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/PPSEff.h"

void RunFwdProtons(const TString filename,
                      TString outname,
                      TH1F *normH, 
                      TH1F *genPU, 					  
                      TString era,
                      Bool_t debug=false);

class PPSEfficiencies
{
  public:
  
  // Constructor
  PPSEfficiencies(){}
  
  // Destructor 
	~PPSEfficiencies() {}
 
  // class function:
  float GetPixelEff(float, float, int);
  float GetStripEff(float, float, int);
  void init(TString); 
  void init(); 
  private:
    TH2D * h_pixels45, * h_pixels56;
    TH2D * h_strips45, * h_strips56;
    
    TFile * f_strips;
    TFile * f_pixels;
    
    TString path = "/eos/project/c/ctpps/subsystems";
    int bin, isData=false;
};

#endif


