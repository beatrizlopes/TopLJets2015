#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/JSONWrapper.h"

#include "TFile.h"
#include "TSystem.h"

#include <iostream>


using namespace std;

//
EfficiencyScaleFactorsWrapper::EfficiencyScaleFactorsWrapper(bool isData,TString era)
{
  if(isData) return;
  init(era);
}

EfficiencyScaleFactorsWrapper::EfficiencyScaleFactorsWrapper(bool isData,TString era,std::map<TString,TString> cfgMap):
  cfgMap_(cfgMap)
{
  if(isData) return;
  init(era);
}

//
void EfficiencyScaleFactorsWrapper::init(TString era)
{
  if(era.Contains("era2017")) era_=2017;
  if(era.Contains("era2016")) era_=2016;

  cout << "[EfficiencyScaleFactorsWrapper]" << endl
       << "\tStarting efficiency scale factors for " << era << endl
       << "\tWarnings: no trigger SFs for any object" << endl
       << "\t          uncertainties returned are of statistical nature only" << endl
       << "\tDon't forget to fix these and update these items!" << endl;

  //ELECTRONS
  TString e_recoSF,e_idSF;
  if(era_==2016){
    e_recoSF=era+"/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root";
    e_idSF=era+"/2016LegacyReReco_ElectronTight_Fall17V2.root";
  }else if(era_==2017){
    e_recoSF=era+"/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root";
    e_idSF=era+"/2017_ElectronTight.root";
  }

  gSystem->ExpandPathName(e_recoSF);
  TFile *fIn=TFile::Open(e_recoSF);
  if(fIn && !fIn->IsZombie()) {
    cout << "electrons: reco SF from " << e_recoSF << endl;
    scaleFactorsH_["e_rec"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
    scaleFactorsH_["e_rec"]->SetDirectory(0);
    fIn->Close();
  }

  gSystem->ExpandPathName(e_idSF);
  fIn=TFile::Open(e_idSF);
  if(fIn && !fIn->IsZombie()) {
    cout << "electrons: id SF from " << e_idSF << endl;
    scaleFactorsH_["e_id"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
    scaleFactorsH_["e_id"]->SetDirectory(0);     
    fIn->Close();
  }

  
  //MUONS
  std::vector<float>   lumiWgts;
  std::vector<TString> m_tkSF, m_idSF,m_isoSF, m_trigSF;
  if(era_==2016){
    m_tkSF.push_back( era+"/MuonTracking_EfficienciesAndSF_BCDEF.root");
    m_tkSF.push_back( era+"/MuonTracking_EfficienciesAndSF_GH.root");
    m_idSF.push_back( era+"/RunBCDEF_SF_ID.root");
    m_idSF.push_back( era+"/RunGH_SF_ID.root");
    m_isoSF.push_back( era+"/RunBCDEF_SF_ISO.root");
    m_isoSF.push_back( era+"/RunGH_SF_ISO.root");
    lumiWgts.push_back(0.5);
    lumiWgts.push_back(0.5);
  }
  if(era_==2017){
    m_tkSF.push_back( era+"/Run2_UL_2017_Efficiency_muon_generalTracks_Run2017_UL_trackerMuon.root");
    m_idSF.push_back( era+"/Run2_UL_2017_2017_Z_Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root");
    m_isoSF.push_back( era+"/Run2_UL_2017_2017_Z_Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root");
    m_trigSF.push_back( era+"/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root");
    lumiWgts.push_back(1.0);
  }
  else{cout << "ERROR: Wrong era provided (era="<<era<<")"<<endl;}

  for(size_t i=0; i<m_tkSF.size(); i++) {
    gSystem->ExpandPathName(m_tkSF[i]);
    fIn=TFile::Open(m_tkSF[i]);
    if(fIn && !fIn->IsZombie()) {
      cout << "muons: tk SF from " << m_tkSF[i] << " with weight " << lumiWgts[i] << endl;
	  if(era_==2016){
        scaleFactorsGr_["m_tk"]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_aeta_dr030e030_corr");
	  }
	  if(era_==2017){
        scaleFactorsH_["m_tk"]=(TH2F *)fIn->Get("NUM_TrackerMuons_DEN_genTracks_abseta_pt")->Clone();
        scaleFactorsH_["m_tk"]->SetDirectory(0);
	  }
      fIn->Close();
    }
  }
  for(size_t i=0; i<m_idSF.size(); i++) {
    gSystem->ExpandPathName(m_idSF[i]);
    fIn=TFile::Open(m_idSF[i]);
    if(fIn && !fIn->IsZombie()) {
      cout << "muons: id SF from " << m_idSF[i] << " with weight " << lumiWgts[i] << endl;
      scaleFactorsH_["m_id"]=(TH2F *)fIn->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt")->Clone();
      scaleFactorsH_["m_id"]->SetDirectory(0);
      fIn->Close();
    }
  }
  for(size_t i=0; i<m_isoSF.size(); i++) {
    gSystem->ExpandPathName(m_isoSF[i]);
    fIn=TFile::Open(m_isoSF[i]);
    if(fIn && !fIn->IsZombie()) {
      cout << "muons: iso SF from " << m_isoSF[i] << " with weight " << lumiWgts[i] << endl;
      scaleFactorsH_["m_iso"]=(TH2F *)fIn->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt")->Clone();
      scaleFactorsH_["m_iso"]->SetDirectory(0);
      fIn->Close();
    }
  }
  for(size_t i=0; i<m_trigSF.size(); i++) {
    gSystem->ExpandPathName(m_trigSF[i]);
    fIn=TFile::Open(m_trigSF[i]);
    if(fIn && !fIn->IsZombie()) {
      cout << "muons: Trigger SF from " << m_trigSF[i] << " with weight " << lumiWgts[i] << endl;
	  if(era_==2017){
		//scaleFactorsH_["m_trig"]=(TH2F *)fIn->Get("IsoMu27_PtEtaBins/abseta_pt_ratio")->Clone();
		scaleFactorsH_["m_trig"]=(TH2F *)fIn->Get("NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt")->Clone();
		scaleFactorsH_["m_trig"]->SetDirectory(0);
	  }
	  else{ cout << "ERROR: no SF awailable for 2016 or 2018!"<<endl;}
      fIn->Close();
    }
  }
  
}

//
EffCorrection_t EfficiencyScaleFactorsWrapper::getDileptonTriggerCorrection(std::vector<Particle> &leptons){

  int dilCode=abs(leptons[0].id()*leptons[1].id());
  EffCorrection_t corr(1.0,0.0);
  if(era_==2016) {
    //values from AN 2016/392 (v3) 
    //as scale factors are approximately constant as function of number of jets take:
    //   - central value from events with 2 jets
    //   - uncertainty as 2 jets unc. + max. difference wrt to other jet mult
    //the offline ids to use must be tight ones
    float leta(fabs(leptons[0].Eta()));
    if(dilCode==11*13) {
      if(leta<0.09)      { corr.first=0.995; corr.second=TMath::Sqrt(0.002*0.002+0.005*0.005); }
      else if (leta<1.2) { corr.first=0.996; corr.second=TMath::Sqrt(0.003*0.003+0.011*0.011); }
      else if (leta<2.1) { corr.first=0.994; corr.second=TMath::Sqrt(0.003*0.003+0.019*0.019); }
      else               { corr.first=0.992; corr.second=TMath::Sqrt(0.011*0.011+0.016*0.016); } //0 jets has too large error
    }
    if(dilCode==11*11) {
      if(leta<0.09)      { corr.first=0.991; corr.second=TMath::Sqrt(0.002*0.002+0.027*0.027); }
      else if (leta<1.2) { corr.first=0.994; corr.second=TMath::Sqrt(0.004*0.004+0.018*0.018); }
      else if (leta<2.1) { corr.first=0.992; corr.second=TMath::Sqrt(0.004*0.004+0.025*0.025); }
      else               { corr.first=0.982; corr.second=TMath::Sqrt(0.014*0.014+0.029*0.029); }
    }
    if(dilCode==13*13) {
      if(leta<0.09)      { corr.first=0.994; corr.second=TMath::Sqrt(0.002*0.003+0.015*0.015); }
      else if (leta<1.2) { corr.first=0.994; corr.second=TMath::Sqrt(0.003*0.003+0.021*0.021); }
      else if (leta<2.1) { corr.first=0.989; corr.second=TMath::Sqrt(0.002*0.002+0.015*0.015); }
      else               { corr.first=0.977; corr.second=TMath::Sqrt(0.009*0.009+0.014*0.014); } //0 jets has too large error
    }
    
  }

  return corr;
}



//
EffCorrection_t EfficiencyScaleFactorsWrapper::getPhotonTrigCorrection(float apt,float mjj){
  
  EffCorrection_t corr(1.0,0.0);
  if(era_==2017){
    if(apt>=200){
      corr.first  = 0.5*0.993*(1.+TMath::Erf((apt-201.5)/(TMath::Sqrt(2.)*14.4)));
      corr.second = 0.03*corr.first;
    }else{
      corr.first  = 0.5*1.000*(1.+TMath::Erf((apt-71.73)/(TMath::Sqrt(2.)*1.005)));
      corr.first *= 0.5*1.000*(1.+TMath::Erf((mjj-264.906)/(TMath::Sqrt(2.)*100.)));
      float a_relUnc((0.2-0.03)/2.);
      float b_relUnc(0.2-a_relUnc*3.);
      float relUnc=min(max(a_relUnc*(mjj/1000.)+b_relUnc,0.03),0.2);
      corr.second=relUnc*corr.first;
    }    
  }else{
    if(apt>=200){
      corr.first  = 0.5*0.999*(1.+TMath::Erf((apt-162.62)/(TMath::Sqrt(2.)*4.391)));
      corr.second = 0.03*corr.first;
    }else{
      corr.first  = 0.5*1.000*(1.+TMath::Erf((apt-70.01)/(TMath::Sqrt(2.)*1.196)));
      corr.first *= 0.5*1.000*(1.+TMath::Erf((mjj-260.948)/(TMath::Sqrt(2.)*174.5)));
      float a_relUnc((0.05-0.03)/2.);
      float b_relUnc(0.05-a_relUnc*3.);
      float relUnc=min(max(a_relUnc*(mjj/1000.)+b_relUnc,0.03),0.05);
      corr.second = relUnc*corr.first;
    }
  }

  return corr;
}



//
EffCorrection_t EfficiencyScaleFactorsWrapper::getTriggerCorrection(std::vector<Particle> leptons, 
                                                                    std::vector<Particle> photons,
                                                                    std::vector<Particle> jets,
                                                                    TString period)
{
  EffCorrection_t corr(1.0,0.0);
  if(era_==2017)
    {
      if(leptons.size()==1)
        {
          TString hname(abs(leptons[0].id())==11 ? "e_trig" : "m_trig");          
          if( abs(leptons[0].id())==13 and scaleFactorsH_.find(hname)!=scaleFactorsH_.end() )
            {
              TH1 *h=scaleFactorsH_[hname];
              float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[0].eta())),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

              float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(leptons[0].pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);

              corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
              corr.second=h->GetBinError(etaBinForEff,ptBinForEff);
            }
          //electron histogram uses eta, not abs(eta)
          else if( abs(leptons[0].id())==11 and scaleFactorsH_.find(hname)!=scaleFactorsH_.end() )
            {
              TH1 *h=scaleFactorsH_[hname];
              float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(leptons[0].eta()),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

              float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(leptons[0].pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);

              corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
              corr.second=h->GetBinError(etaBinForEff,ptBinForEff);
            }
        }
    }

  return corr;
}

//
EffCorrection_t EfficiencyScaleFactorsWrapper::getOfflineCorrection(Particle p,TString period)
{
  return getOfflineCorrection(p.id(),p.pt(),p.eta(),period);
}

//
EffCorrection_t EfficiencyScaleFactorsWrapper::getOfflineCorrection(int pdgId,float pt,float eta,TString period)
{
  EffCorrection_t corr(1.0,0.00);

  TString corrSteps[]={"rec","tk","id","iso"};
  for(size_t icor=0; icor<sizeof(corrSteps)/sizeof(TString); icor++) {

    TString idstr("e");
    if(abs(pdgId)==13) idstr="m";
    if(abs(pdgId)==22) idstr="g";
    TString hname(idstr+"_"+corrSteps[icor]);
    if(abs(pdgId)==13) hname += period;
    
    Float_t iSF(1.0), iSFUnc(0.0);
    if(scaleFactorsGr_.find(hname)!=scaleFactorsGr_.end() ) {
      Double_t x(0.),xdiff(9999.),y(0.);
      for(Int_t ip=0; ip<scaleFactorsGr_[hname]->GetN(); ip++)
        {
          scaleFactorsGr_[hname]->GetPoint(ip,x,y);
          float ixdiff(TMath::Abs(fabs(eta)-x));
          if(ixdiff>xdiff) continue;
          xdiff=ixdiff;
          iSF=y;
          iSFUnc=scaleFactorsGr_[hname]->GetErrorY(ip);
        }
    }
    else if(scaleFactorsH_.find(hname)!=scaleFactorsH_.end() ) {
      TH2 *h=scaleFactorsH_[hname];
      Double_t xval(eta), yval(pt);
      if(idstr=="m") {
        if(era_==2017) {
          xval=pt;
          yval=fabs(eta);
        }
      }
      xval=TMath::Max( TMath::Min(xval,h->GetXaxis()->GetXmax()-0.001), h->GetXaxis()->GetXmin() );
      Int_t xbin=h->GetXaxis()->FindBin(xval);
      yval=TMath::Max( TMath::Min(yval,h->GetYaxis()->GetXmax()-0.001), h->GetYaxis()->GetXmin() );
      Int_t ybin=h->GetYaxis()->FindBin(yval);      

      iSF=h->GetBinContent(xbin,ybin);
      iSFUnc=h->GetBinError(xbin,ybin);
    }
   
    //corr.second = sqrt(pow(iSFUnc*corr.first,2)+pow(iSF*corr.second,2));
    corr.second = sqrt( iSFUnc*iSFUnc + corr.second*corr.second );
    corr.first  *= iSF;
  }
     
  //
  return corr;
}

//
EfficiencyScaleFactorsWrapper::~EfficiencyScaleFactorsWrapper()
{
}
