#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  MiniEvent_t()
  {
    g_nw=0; g_npsw=0; ng=0; ngtop=0; 
    ngamma=0; nl=0; nj=0; nbj=0;
  }

  static const int MAXWEIGHTS   =  10;  // was 1285 in ttbar MC
  static const int MAXPSWEIGHTS =  46;
  static const int MAXGENPAR    =  20;  // initialy was 500
  static const int MAXGENTOPAR  =  25;
  static const int MAXGAMMA     =  20;
  static const int MAXLEP       =  20;
  static const int MAXJET       =  100;
  static const int MAXMETSYS    =  14;
  static const int MAXJETSYS    =  29;
  static const int MAXRAWMU     =  50;
  static const int MAXPROTONS   =  50;
  
  Bool_t isData;
  UInt_t run,lumi;
  ULong64_t event;
  Float_t beamXangle, instLumi;
   
  //gen level event
  Int_t g_id1, g_id2;
  Float_t g_x1, g_x2, g_qscale;
  Int_t g_pu,g_putrue;
  Int_t g_nw, g_npsw;
  Float_t g_w[MAXWEIGHTS], g_psw[MAXPSWEIGHTS];
  Int_t ng,ngtop;
  Int_t g_id[MAXGENPAR],g_bid[MAXGENPAR],g_tagCtrs[MAXGENPAR];
  Bool_t g_isSemiLepBhad[MAXGENPAR];
  Float_t g_pt[MAXGENPAR],g_eta[MAXGENPAR],g_phi[MAXGENPAR],g_m[MAXGENPAR],g_xb[MAXGENPAR],g_xbp[MAXGENPAR];
  Int_t gtop_id[MAXGENTOPAR];
  Float_t gtop_pt[MAXGENTOPAR],gtop_pz[MAXGENTOPAR],gtop_eta[MAXGENTOPAR],gtop_phi[MAXGENTOPAR],gtop_m[MAXGENTOPAR];
  Int_t g_nchPV;
  Float_t g_sumPVChPt,g_sumPVChPz,g_sumPVChHt;

  //reco level event
  Int_t nvtx;
  Int_t triggerBits,addTriggerBits;
  Int_t zeroBiasPS;
  Float_t rho;

  //leptons
  Int_t nl;
  Bool_t l_isPromptFinalState[MAXLEP], l_isDirectPromptTauDecayProductFinalState[MAXLEP];
  Int_t l_id[MAXLEP],l_charge[MAXLEP],l_pid[MAXLEP],l_g[MAXLEP];
  Float_t l_pt[MAXLEP],l_eta[MAXLEP],l_phi[MAXLEP], l_mass[MAXLEP], l_highpt[MAXLEP],
    l_scaleUnc1[MAXLEP], l_scaleUnc2[MAXLEP], l_scaleUnc3[MAXLEP], l_scaleUnc4[MAXLEP], l_scaleUnc5[MAXLEP], l_scaleUnc6[MAXLEP], l_scaleUnc7[MAXLEP],
    l_miniIso[MAXLEP], l_chargedHadronIso[MAXLEP], l_relIso[MAXLEP], l_ip3d[MAXLEP], l_ip3dsig[MAXLEP],l_mva[MAXLEP],l_mvaCats[MAXLEP];

  Int_t ngamma;
  Bool_t gamma_isPromptFinalState[MAXGAMMA];
  Int_t gamma_pid[MAXGAMMA],gamma_idFlags[MAXGAMMA],gamma_g[MAXGAMMA];
  Float_t gamma_pt[MAXGAMMA],gamma_eta[MAXGAMMA],gamma_phi[MAXGAMMA],
    gamma_scaleUnc1[MAXGAMMA],gamma_scaleUnc2[MAXGAMMA],gamma_scaleUnc3[MAXGAMMA],gamma_scaleUnc4[MAXGAMMA],gamma_scaleUnc5[MAXGAMMA],gamma_scaleUnc6[MAXGAMMA],gamma_scaleUnc7[MAXGAMMA],
    gamma_mva[MAXGAMMA], gamma_mvaCats[MAXGAMMA],
    gamma_chargedHadronIso[MAXGAMMA],gamma_neutralHadronIso[MAXGAMMA],gamma_photonIso[MAXGAMMA],gamma_hoe[MAXGAMMA],gamma_r9[MAXGAMMA],gamma_sieie[MAXGAMMA];

  Int_t nj, nbj;
  Float_t j_pt[MAXJET],j_eta[MAXJET],j_phi[MAXJET],j_mass[MAXJET],j_area[MAXJET],j_rawsf[MAXJET];
  Float_t j_jerUp[MAXJET],j_jerDn[MAXJET],j_jecUp[MAXJETSYS][MAXJET],j_jecDn[MAXJETSYS][MAXJET];
  Float_t j_csv[MAXJET],j_deepcsv[MAXJET],j_pumva[MAXJET],j_emf[MAXJET],j_qg[MAXJET];
  Float_t j_c2_00[MAXJET],j_c2_02[MAXJET],j_c2_05[MAXJET],j_c2_10[MAXJET],j_c2_20[MAXJET];
  Float_t j_zg[MAXJET],j_mult[MAXJET],j_gaptd[MAXJET],j_gawidth[MAXJET],j_gathrust[MAXJET],j_tau32[MAXJET],j_tau21[MAXJET];
  Float_t j_vtxmass[MAXJET],j_vtx3DVal[MAXJET],j_vtx3DSig[MAXJET],j_vtxpx[MAXJET],j_vtxpy[MAXJET],j_vtxpz[MAXJET];
  Float_t e_j_px[MAXJET], e_j_py[MAXJET], e_j_pz[MAXJET]; // error on jet components
  
  Bool_t j_btag[MAXJET];
  Int_t j_vtxNtracks[MAXJET],j_flav[MAXJET],j_id[MAXJET],j_pid[MAXJET],j_hadflav[MAXJET],j_g[MAXJET];

  //met
  Float_t met_pt,met_phi,met_sig;
  Float_t e_met_px, e_met_py, e_met_pxpy;  // error on estimation of MET
  Float_t met_ptShifted[MAXMETSYS],met_phiShifted[MAXMETSYS];
  Int_t met_filterBits;

  //event energy fluxes (PF-based)
  Int_t nchPV;
  Float_t sumPVChPt,sumPVChPz,sumPVChHt;
  Int_t nPFCands[8],nPFChCands[8];
  Float_t sumPFHt[8],sumPFEn[8],sumPFPz[8],sumPFChHt[8],sumPFChEn[8],sumPFChPz[8];

  //CTPPS protons
  Short_t nfwdtrk,fwdtrk_pot[MAXPROTONS],fwdtrk_method[MAXPROTONS];
  Float_t fwdtrk_thetax[MAXPROTONS],fwdtrk_thetay[MAXPROTONS],
    fwdtrk_vx[MAXPROTONS],fwdtrk_vy[MAXPROTONS],fwdtrk_vz[MAXPROTONS],
    fwdtrk_time[MAXPROTONS],fwdtrk_timeError[MAXPROTONS],
    fwdtrk_chisqnorm[MAXPROTONS],fwdtrk_xi[MAXPROTONS],fwdtrk_xiError[MAXPROTONS],fwdtrk_xiSF[MAXPROTONS],fwdtrk_t[MAXPROTONS];
  
  //PPS tracks (needed for low PU run)
  Short_t nppstrk,ppstrk_pot[MAXPROTONS];
  Float_t ppstrk_x[MAXPROTONS],ppstrk_y[MAXPROTONS], ppstrk_xUnc[MAXPROTONS],ppstrk_yUnc[MAXPROTONS],
    ppstrk_tx[MAXPROTONS],ppstrk_ty[MAXPROTONS],ppstrk_txUnc[MAXPROTONS],ppstrk_tyUnc[MAXPROTONS],
    ppstrk_chisqnorm[MAXPROTONS];
  //ppstrk_t[MAXPROTONS],ppstrk_tUnc[MAXPROTONS]; // UFSD only (2018)

  //these are crazy variables for the cross check
  Int_t nrawmu;
  Short_t rawmu_pt[MAXRAWMU],rawmu_eta[MAXRAWMU],rawmu_phi[MAXRAWMU];
  Int_t rawmu_pid[MAXRAWMU];
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev,Int_t, std::vector<std::string>);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev);

#endif
