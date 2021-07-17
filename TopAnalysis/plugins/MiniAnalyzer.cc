//
// -*- C++ -*-
//
// Package:    TopLJets2015/TopAnalysis
// Class:      MiniAnalyzer
//
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Qamar Ul Hassan
//         Created:  Sun, 13 Jul 2014 06:22:18 GMT
//
//

// Main tools
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// TOTEM/PPS related
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLiteFwd.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"
#include "TopLJets2015/TopAnalysis/interface/PPSEff.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"

// Other tools
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/MyIPTools.h"
#include "TopLJets2015/TopAnalysis/interface/JetShapes.h"
#include "TopLJets2015/TopAnalysis/interface/RoccoR.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "FWCore/Framework/interface/Run.h"
#include "TopQuarkAnalysis/BFragmentationAnalyzer/interface/BFragmentationAnalyzerUtils.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom2.h"
#include "TTree.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;

typedef math::XYZTLorentzVector LorentzVector;

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual void endRun(const edm::Run&,const edm::EventSetup&);
private:
  virtual void beginJob() override;
  void genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  float getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			 const reco::Candidate* ptcl,
                         float r_iso_min, float r_iso_max, float kt_scale,
                         bool charged_only);

  bool isSoftMuon(const reco::Muon & recoMu,const reco::Vertex &vertex);
  bool isMediumMuon2016ReReco(const reco::Muon & recoMu);

  // member data
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorevtToken_;
  edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
  edm::EDGetTokenT<LHERunInfoProduct> generatorRunInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>  > genPhotonsToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>  > genLeptonsToken_, genJetsToken_;
  edm::EDGetTokenT<reco::METCollection> genMetsToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPUProtonsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> particleLevelToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_,metFilterBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_,l1triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_;
  edm::EDGetTokenT<edm::View<pat::Photon>  >  photonToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<edm::DetSetVector<TotemRPRecHit>> tokenStripHits_;
  edm::EDGetTokenT<edm::DetSetVector<TotemRPUVPattern>> tokenStripPatterns_;
  edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite> > ctppsToken_;
  std::vector< edm::EDGetTokenT<std::vector<reco::ForwardProton> > > tokenRecoProtons_;
  edm::EDGetTokenT<bool> BadChCandFilterToken_,BadPFMuonFilterToken_,BadPFMuonDzFilterToken_;

  std::unordered_map<std::string,TH1*> histContainer_;

  std::string jetIdToUse_, FilterType_, EGIDVersion_;
  std::vector<JetCorrectionUncertainty *> jecCorrectionUncs_;

  std::vector<std::string> triggersToUse_,metFiltersToUse_,ListVars_;

  bool saveTree_;
  TTree *tree_;
  MiniEvent_t ev_;

  RoccoR *muonRC_;
  
  PPSEff *PPS_eff_; int runNumber_;

  edm::Service<TFileService> fs;

  //counters
  int nrecleptons_, nrecphotons_, ngleptons_, ngphotons_, nmultiprotons_[2];
  int nrecjets_, nrecbjets_;

  //apply filter to save tree
  bool applyFilt_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig) :
  generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  generatorevtToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator",""))),
  generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
  generatorRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>({"externalLHEProducer"})),
  puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),
  genPhotonsToken_(consumes<std::vector<reco::GenParticle> >(edm::InputTag("particleLevel:photons"))),
  genLeptonsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:leptons"))),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:jets"))),
  genMetsToken_(consumes<reco::METCollection>(edm::InputTag("particleLevel:mets"))),
  genParticlesToken_(consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"))),
  genPUProtonsToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("PUprotons"))),
  prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  particleLevelToken_(consumes<reco::GenParticleCollection>(edm::InputTag("particleLevel"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  metFilterBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  l1triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("l1prescales"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  tokenStripHits_( mayConsume<edm::DetSetVector<TotemRPRecHit>>(iConfig.getParameter<edm::InputTag>("tagStripHits")) ),
  tokenStripPatterns_( consumes<edm::DetSetVector<TotemRPUVPattern>>(iConfig.getParameter<edm::InputTag>("tagStripPatterns")) ),
  ctppsToken_(consumes<std::vector<CTPPSLocalTrackLite> >(iConfig.getParameter<edm::InputTag>("ctppsLocalTracks"))),
  BadChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badChCandFilter"))),
  BadPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
  BadPFMuonDzFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuonDzFilter"))),
  saveTree_( iConfig.getParameter<bool>("saveTree") ),
  runNumber_( iConfig.getUntrackedParameter<int>("runNumber") ),
  applyFilt_( iConfig.getParameter<bool>("applyFilt") )
{
	
  tokenRecoProtons_.push_back( consumes<std::vector<reco::ForwardProton>>(iConfig.getParameter<InputTag>("tagRecoProtons")));
  tokenRecoProtons_.push_back( consumes<std::vector<reco::ForwardProton>>(iConfig.getParameter<InputTag>("tagMultiRecoProtons")));

  	
  //now do what ever initialization is needed
  electronToken_      = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  photonToken_        = mayConsume<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"));
  ListVars_           = iConfig.getParameter<std::vector<std::string> >("ListVars");
  FilterType_         = iConfig.getParameter<std::string>("FilterType");
  triggersToUse_      = iConfig.getParameter<std::vector<std::string> >("triggersToUse");
  metFiltersToUse_    = iConfig.getParameter<std::vector<std::string> >("metFiltersToUse");
  jetIdToUse_         = iConfig.getParameter<std::string>("jetIdToUse");
  EGIDVersion_        = iConfig.getParameter<std::string>("EGIDVersion");
  std::string jecUncFile(iConfig.getParameter<std::string>("jecUncFile"));
  //std::string jecUncFile(edm::FileInPath(iConfig.getParameter<std::string>("jecUncFile")).fullPath());
  for(auto name : iConfig.getParameter<std::vector<std::string> >("jecUncSources") ) {
    JetCorrectorParameters *p = new JetCorrectorParameters(jecUncFile,name.c_str());
    jecCorrectionUncs_.push_back(new JetCorrectionUncertainty(*p));
  }
  
  muonRC_ = new RoccoR();
  muonRC_->init(iConfig.getParameter<std::string>("RoccoR"));
  //muonRC_->init(edm::FileInPath(iConfig.getParameter<std::string>("RoccoR")).fullPath());
  
  PPS_eff_ = new PPSEff();
  //PPS_eff_ = new PPSEff(edm::FileInPath(iConfig.getParameter<std::string>("PPS_pixelEff")).fullPath());

  histContainer_["triggerList"] = fs->make<TH1F>("triggerList", ";Trigger bits;",triggersToUse_.size(),0,triggersToUse_.size());
  histContainer_["triggerPrescale"] = fs->make<TH1D>("triggerPrescale", ";Trigger prescale sum;",triggersToUse_.size(),0,triggersToUse_.size());
  for(size_t i=0; i<triggersToUse_.size(); i++) histContainer_["triggerList"] ->GetXaxis()->SetBinLabel(i+1,triggersToUse_[i].c_str());
  histContainer_["counter"]    = fs->make<TH1F>("counter", ";Counter;Events",4,0,4);
  histContainer_["RPcount"]    = fs->make<TH2F>("RPcount", ";Nhits (arm=0);NHits (arm=1)",3,0,3,3,0,3);
  histContainer_["fidcounter"] = (TH1 *)fs->make<TH2F>("fidcounter",    ";Variation;Events", 1500, 0., 1500.,11,0,11);
  histContainer_["pu"]         = fs->make<TH1F>("pu",      ";Pileup observed;Events / 1",100,0,100);
  histContainer_["putrue"]     = fs->make<TH1F>("putrue",  ";Pileup true;Events / 0.1",100,0,100);
  for(std::unordered_map<std::string,TH1*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2();

  //create a tree for the selected events
  if(saveTree_)
    {
      tree_ = fs->make<TTree>("tree","tree with selected events");
      createMiniEventTree(tree_,ev_,jecCorrectionUncs_.size(),ListVars_);
    }
}


//
MiniAnalyzer::~MiniAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
void MiniAnalyzer::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // PILEUP
  //
  edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(puToken_,PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator ipu;
  for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu)
    {
      if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
      ev_.g_pu=ipu->getPU_NumInteractions();
      ev_.g_putrue=ipu->getTrueNumInteractions();
    }
  histContainer_["pu"]->Fill(ev_.g_pu);
  histContainer_["putrue"]->Fill(ev_.g_putrue);

  //
  // GENERATOR WEIGHTS
  //
  ev_.g_nw=0; ev_.g_w[0]=1.0;
  edm::Handle<GenEventInfoProduct> evt;
  iEvent.getByToken( generatorToken_,evt);
  if(evt.isValid())
    {
      ev_.g_w[0] = evt->weight();
      ev_.g_nw++;

      //PDF info
      ev_.g_qscale = evt->pdf()->scalePDF;
      ev_.g_x1     = evt->pdf()->x.first;
      ev_.g_x2     = evt->pdf()->x.second;
      ev_.g_id1    = evt->pdf()->id.first;
      ev_.g_id2    = evt->pdf()->id.second;
	  
	  //parton shower weights:
	  ev_.g_npsw = evt->weights().size();
	  if(ev_.MAXPSWEIGHTS<ev_.g_npsw){
		  cout << "WARNING: expected MAXN PS weights ("<<ev_.MAXPSWEIGHTS<<") is smaller than the NPS weights in MC ("<<ev_.g_npsw<<")."<<endl;
		  cout <<"\t\t... will store only the first " << ev_.MAXPSWEIGHTS << "weights."<<endl;
		  ev_.g_npsw = ev_.MAXPSWEIGHTS;
	  }
	  for(int i=0; i< ev_.g_npsw; i++) ev_.g_psw[i]=evt->weights().at(i);
	}
  histContainer_["counter"]->Fill(1,ev_.g_w[0]);

  //alternative weights for systematics
  // see https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#Retrieving_the_weights for details
  // convention from https://twiki.cern.ch/twiki/pub/CMS/TopModGen/pwg-rwl.txt
  edm::Handle<LHEEventProduct> evet;
  iEvent.getByToken(generatorlheToken_, evet);
  if(evet.isValid()) 
    {

       // Store 5 scale weights and 103 PDF weights
	   ev_.g_nw+=(5+103);
	   //float origW = evet->originalXWGTUP();
	   for(int i=1;i<ev_.g_nw;i++) ev_.g_w[i]=0;   
       //cout << "Original wegiht = " << origW << endl;
	   for (unsigned int i=0; i<evet->weights().size(); i++) {
		   int id=atoi(evet->weights()[i].id.c_str());
		   // Scale uncertainties
		   if (id == 1001) ev_.g_w[1]=evet->weights()[i].wgt; // nominal
		   if (id == 1002) ev_.g_w[2]=evet->weights()[i].wgt; // muF*2.0
		   if (id == 1003) ev_.g_w[3]=evet->weights()[i].wgt; // muF*0.5
		   if (id == 1004) ev_.g_w[4]=evet->weights()[i].wgt; // muR*2.0
		   if (id == 1007) ev_.g_w[5]=evet->weights()[i].wgt; // muR*0.5
		   
		   // PDF uncertanties (NNPDF31_nnlo_hessian_pdfas, 306000)
		   if (id == 1010) ev_.g_w[6]=evet->weights()[i].wgt; // central value
		   if (id>=1011 && id <= 1110) ev_.g_w[id-1004]=evet->weights()[i].wgt; // PDF eig.
		   if (id == 1111) ev_.g_w[107]=evet->weights()[i].wgt; // aS=0.116
		   if (id == 1112) ev_.g_w[108]=evet->weights()[i].wgt; // aS=0.120
		   
		   //cout << i << ": "<< id << ", val = " << evet->weights()[i].wgt << endl;

		   
	   }

	   if(ev_.MAXWEIGHTS<ev_.g_nw){
		 cout << "WARNING: expected MAXN weights ("<<ev_.MAXWEIGHTS<<") is smaller than the N weights in MC ("<<ev_.g_nw<<")."<<endl;
		 cout <<"\t\t... will store only the first " << ev_.MAXWEIGHTS << "weights."<<endl;
		 ev_.g_nw = ev_.MAXWEIGHTS-1;
	   }	  
    }

  //
  // GENERATOR LEVEL EVENT
  //
  ev_.ng=0;
  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken_,genJets);
  std::map<const reco::Candidate *,int> jetConstsMap;
  //edm::Handle<edm::ValueMap<float> > petersonFrag;
  //iEvent.getByToken(petersonFragToken_,petersonFrag);
  int ngjets(0),ngbjets(0);
  if(genJets.isValid()){
    for(auto genJet=genJets->begin(); genJet!=genJets->end(); ++genJet)
      {
        edm::Ref<std::vector<reco::GenJet> > genJetRef(genJets,genJet-genJets->begin());

        //map the gen particles which are clustered in this jet
        JetFragInfo_t jinfo=analyzeJet(*genJet);

        std::vector< const reco::Candidate * > jconst=genJet->getJetConstituentsQuick();
        for(size_t ijc=0; ijc <jconst.size(); ijc++) jetConstsMap[ jconst[ijc] ] = ev_.ng;
        ev_.g_tagCtrs[ev_.ng]       = (jinfo.nbtags&0xf) | ((jinfo.nctags&0xf)<<4) | ((jinfo.ntautags&0xf)<<8);
        ev_.g_xb[ev_.ng]            = jinfo.xb;
        ev_.g_bid[ev_.ng]           = jinfo.leadTagId;
        ev_.g_isSemiLepBhad[ev_.ng] = jinfo.hasSemiLepDecay;
        ev_.g_id[ev_.ng]   = genJet->pdgId();
        ev_.g_pt[ev_.ng]   = genJet->pt();
        ev_.g_eta[ev_.ng]  = genJet->eta();
        ev_.g_phi[ev_.ng]  = genJet->phi();
        ev_.g_m[ev_.ng]    = genJet->mass();
        ev_.ng++;

        //gen level selection
        if(genJet->pt()>25 && fabs(genJet->eta())<2.5)
          {
            ngjets++;
            if(abs(genJet->pdgId())==5) ngbjets++;
          }
      }
  }

  //leptons
  edm::Handle<std::vector<reco::GenJet> > dressedLeptons;
  iEvent.getByToken(genLeptonsToken_,dressedLeptons);
  if(dressedLeptons.isValid()) {
    for(auto genLep = dressedLeptons->begin();  genLep != dressedLeptons->end(); ++genLep)
      {
        //map the gen particles which are clustered in this lepton
        std::vector< const reco::Candidate * > jconst=genLep->getJetConstituentsQuick();
        for(size_t ijc=0; ijc <jconst.size(); ijc++) jetConstsMap[ jconst[ijc] ] = ev_.ng;

        ev_.g_pt[ev_.ng]   = genLep->pt();
        ev_.g_id[ev_.ng]   = genLep->pdgId();
        ev_.g_eta[ev_.ng]  = genLep->eta();
        ev_.g_phi[ev_.ng]  = genLep->phi();
        ev_.g_m[ev_.ng]    = genLep->mass();
        ev_.ng++;

        //gen level selection
        if(genLep->pt()>25 && fabs(genLep->eta())<2.5) ngleptons_++;
      }
  }
    
  edm::Handle<std::vector<reco::GenParticle> > genPhotons;
  iEvent.getByToken(genPhotonsToken_,genPhotons);
  if(genPhotons.isValid()){
    for(auto genPhoton = genPhotons->begin();  genPhoton != genPhotons->end(); ++genPhoton)
      {
        if(genPhoton->pt()<15) continue;
        if(fabs(genPhoton->eta())>2.5) continue;

        ev_.g_pt[ev_.ng]   = genPhoton->pt();
        ev_.g_id[ev_.ng]   = genPhoton->pdgId();
        ev_.g_eta[ev_.ng]  = genPhoton->eta();
        ev_.g_phi[ev_.ng]  = genPhoton->phi();
        ev_.g_m[ev_.ng]    = genPhoton->mass();
        ev_.ng++;

        //gen level selection
        if(genPhoton->pt()>20 && fabs(genPhoton->eta())<2.5) ngphotons_++;
      }
  }

  if(ev_.MAXGENPAR<ev_.ng){
    cout << "ERROR: expected MAXN genpar ("<<ev_.MAXGENPAR<<") is smaller than the N genpar in MC ("<<ev_.ng<<")."<<endl;
	cout <<"\t\t... can cause memory leaks!!!"<<endl;
  }   

  //final state particles
  ev_.g_nchPV=0;
  ev_.g_sumPVChPt=0;
  ev_.g_sumPVChPz=0;
  ev_.g_sumPVChHt=0;
  edm::Handle<pat::PackedGenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_,genParticles);
  LorentzVector pvP4(0,0,0,0);
  if(genParticles.isValid()){
    for (size_t i = 0; i < genParticles->size(); ++i)
      {
        const pat::PackedGenParticle & genIt = (*genParticles)[i];
        if(genIt.pt()<0.9) continue;
        if(genIt.charge()==0) continue;
        ev_.g_nchPV++;
        pvP4+=genIt.p4();
        ev_.g_sumPVChPz+=fabs(genIt.pz());
        ev_.g_sumPVChHt+=genIt.pt();
      }
  }
  ev_.g_sumPVChPt=pvP4.Pt();

  // Save generator particles //
  ev_.ngtop=0; 
  //Save pileup protons (for PPS reconstruction)
  edm::Handle<reco::GenParticleCollection> genPUProtons;
  iEvent.getByToken(genPUProtonsToken_,genPUProtons);
  float XIMIN = 0.01, XIMAX = 0.2, EBEAM=6500;
  bool hasPUProtons (genPUProtons.isValid());
  if(hasPUProtons){
    for (size_t i = 0; i < genPUProtons->size(); ++i)
      {
        const reco::GenParticle & genIt = (*genPUProtons)[i];
        int absid=abs(genIt.pdgId());
        bool outGoingProton( absid==2212 && genIt.status()==1 && fabs(genIt.eta())>4.7 );
		bool tagProton( abs(genIt.pz())>(1. - XIMAX)*EBEAM && abs(genIt.pz())<(1. - XIMIN)*EBEAM);
        if(outGoingProton && tagProton)
          {
            ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId();
            ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
            ev_.gtop_pz[ ev_.ngtop ]  = genIt.pz();
            ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
            ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
            ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
            ev_.ngtop++;
          }
      }
  }
  
  //W,Z and top quarks (lastCopy)
  // don't keep protons if hasPUProtons and not using premix_stage2 modifier
  edm::Handle<reco::GenParticleCollection> prunedGenParticles;
  iEvent.getByToken(prunedGenParticlesToken_,prunedGenParticles);
  if(prunedGenParticles.isValid()){
    for (size_t i = 0; i < prunedGenParticles->size(); ++i)
      {
        const reco::GenParticle & genIt = (*prunedGenParticles)[i];
        int absid=abs(genIt.pdgId());
        bool outGoingProton( absid==2212 && genIt.status()==1 && fabs(genIt.eta())>4.7 && !hasPUProtons);
        bool topLastCopy(absid==6 && genIt.isLastCopy());
        bool wLastCopy(absid==24 && genIt.isLastCopy());
        bool zLastCopy(absid==23 && genIt.isLastCopy());
		bool toKeep = outGoingProton || topLastCopy || zLastCopy || wLastCopy;
        if(toKeep)
          {
            ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId();
            ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
            ev_.gtop_pz[ ev_.ngtop ]  = genIt.pz();
            ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
            ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
            ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
            ev_.ngtop++;
          }

        //save daughters
        if(topLastCopy || wLastCopy || zLastCopy) {
          size_t ndau=genIt.numberOfDaughters();
          for(size_t idau=0; idau<ndau; idau++){
            const reco::Candidate *d=genIt.daughter(idau);
			if(abs(d->pdgId())==24) continue; // keep all W's anyway
            ev_.gtop_id[ ev_.ngtop ]  = d->pdgId();
            ev_.gtop_pt[ ev_.ngtop ]  = d->pt();
            ev_.gtop_pz[ ev_.ngtop ]  = d->pz();
            ev_.gtop_eta[ ev_.ngtop ] = d->eta();
            ev_.gtop_phi[ ev_.ngtop ] = d->phi();
            ev_.gtop_m[ ev_.ngtop ]   = d->mass();
            ev_.ngtop++;
          }
        }
      }
  }
  
  if(ev_.MAXGENTOPAR<ev_.ngtop){
    cout << "ERROR: expected MAXN gentoppar ("<<ev_.MAXGENTOPAR<<") is smaller than the N gentoppar in MC ("<<ev_.ngtop<<")."<<endl;
	cout <<"\t\t... can cause memory leaks!!!"<<endl;
  }  

  //pseudo-tops
/*  edm::Handle<reco::GenParticleCollection> particleLevel;
  iEvent.getByToken(particleLevelToken_,particleLevel);
  for (size_t i = 0; i < particleLevel->size(); ++i)
    {
      const GenParticle & genIt = (*particleLevel)[i];
      ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId()*1000;
      ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
      ev_.gtop_pz[ ev_.ngtop ]  = genIt.pz();
      ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
      ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
      ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
      ev_.ngtop++;
    }
*/
  //gen met
  edm::Handle<reco::METCollection> genMet;
  iEvent.getByToken(genMetsToken_,genMet);
  if(genMet.isValid()){
    ev_.gtop_id[ ev_.ngtop ]  = 0;
    ev_.gtop_pt[ ev_.ngtop ]  = (*genMet)[0].pt();
    ev_.gtop_eta[ ev_.ngtop ] = 0;
    ev_.gtop_phi[ ev_.ngtop ] = (*genMet)[0].phi();
    ev_.gtop_m[ ev_.ngtop ]   = 0;
    ev_.ngtop++;
  }

  //fiducial counters
  for(Int_t iw=0; iw<ev_.g_nw; iw++)
    {
      Double_t x(iw);
      Double_t wgt(ev_.g_w[iw]);
      TH2F *fidCounter=(TH2F *)histContainer_["fidcounter"];
      fidCounter->Fill(x,0.,wgt);
      if(ngleptons_>0)               fidCounter->Fill(x, 1., wgt);
      if(ngleptons_>1)               fidCounter->Fill(x, 2., wgt);
      if(ngleptons_>0 && ngjets>0)   fidCounter->Fill(x, 3., wgt);
      if(ngleptons_>1 && ngjets>0)   fidCounter->Fill(x, 4., wgt);
      if(ngleptons_>0 && ngjets>1)   fidCounter->Fill(x, 5., wgt);
      if(ngleptons_>1 && ngjets>1)   fidCounter->Fill(x, 6., wgt);
      if(ngleptons_>0 && ngjets>2)   fidCounter->Fill(x, 7., wgt);
      if(ngleptons_>1 && ngjets>2)   fidCounter->Fill(x, 8., wgt);
      if(ngleptons_>0 && ngjets>3)   fidCounter->Fill(x, 9., wgt);
      if(ngleptons_>1 && ngjets>3)   fidCounter->Fill(x, 10.,wgt);
    }

}


//
void MiniAnalyzer::recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //VERTICES
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->size()<1) return; // skip the events no PV
  const reco::Vertex &primVtx = vertices->front();
  reco::VertexRef primVtxRef(vertices,0);
  ev_.nvtx=vertices->size();
  ev_.zPV = (ev_.nvtx>0) ? primVtx.z() : 999;
  ev_.zPV2 = (ev_.nvtx>1) ? fabs(primVtx.z() - vertices->at(1).z()): -1;
  unsigned int _second_vertex_index = 1;
  for (size_t ipv = 2; ipv < vertices->size(); ++ipv) {
	  float dz = fabs(primVtx.z() - vertices->at(ipv).z());
	  if (dz < ev_.zPV2) {
		  ev_.zPV2 = dz;
		  _second_vertex_index = ipv;
	  }
  }

  //RHO
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  float rho=*rhoH;
  ev_.rho=rho;

  //TRIGGER INFORMATION
  edm::Handle<edm::TriggerResults> h_trigRes;
  iEvent.getByToken(triggerBits_, h_trigRes);
  std::vector<string> triggerList;
  Service<service::TriggerNamesService> tns;
  tns->getTrigPaths(*h_trigRes,triggerList);
  edm::Handle<pat::PackedTriggerPrescales> h_trigPrescale;
  iEvent.getByToken(triggerPrescales_, h_trigPrescale);
  edm::Handle<pat::PackedTriggerPrescales> h_l1trigPrescale;
  iEvent.getByToken(l1triggerPrescales_, h_l1trigPrescale);
  ev_.triggerBits=0;
  ev_.addTriggerBits=0;
  for (unsigned int i=0; i< h_trigRes->size(); i++)
    {
      if( !(*h_trigRes)[i].accept() ) continue;
      for(size_t itrig=0; itrig<triggersToUse_.size(); itrig++)
	{
	  if (triggerList[i].find(triggersToUse_[itrig])==string::npos) continue;
          int prescale=h_trigPrescale->getPrescaleForIndex(i);
          int l1prescale=h_l1trigPrescale->getPrescaleForIndex(i);
          if(itrig<32)
            ev_.triggerBits |= (1 << itrig);
          else
            ev_.addTriggerBits |= (1 << (itrig-32));
	  histContainer_["triggerList"]->Fill(itrig);
          histContainer_["triggerPrescale"]->Fill(itrig,prescale);
          bool isZeroBias(triggerList[i].find("ZeroBias")!=string::npos);
          if(!isZeroBias) continue;
          ev_.zeroBiasPS=prescale*l1prescale;
	}
    }
  bool passTrigger((ev_.triggerBits + ev_.addTriggerBits)!=0);
  //if(!passTrigger) return; not obvious that triggers are simulated properly
  if(ev_.isData && !passTrigger) return;

  	    
  //
  //PPS local tracks (if present)
  //
  ev_.nppstrk=0;
  edm::Handle<CTPPSLocalTrackLiteCollection> recoPPSTracks;
  iEvent.getByToken(ctppsToken_, recoPPSTracks);
  if(recoPPSTracks.isValid()){
    for (const auto& trk : *recoPPSTracks)
      {
        CTPPSDetId detid(trk.getRPId());
        ev_.ppstrk_pot[ev_.nppstrk]       = 100*detid.arm()+10*detid.station()+detid.rp();
        ev_.ppstrk_x[ev_.nppstrk]         = trk.getX();
        ev_.ppstrk_y[ev_.nppstrk]         = trk.getY();
        ev_.ppstrk_xUnc[ev_.nppstrk]      = trk.getXUnc();
        ev_.ppstrk_yUnc[ev_.nppstrk]      = trk.getYUnc();
        ev_.ppstrk_tx[ev_.nppstrk]        = trk.getTx();
        ev_.ppstrk_ty[ev_.nppstrk]        = trk.getTy();
        ev_.ppstrk_txUnc[ev_.nppstrk]     = trk.getTxUnc();
        ev_.ppstrk_tyUnc[ev_.nppstrk]     = trk.getTyUnc();
        ev_.ppstrk_chisqnorm[ev_.nppstrk] = trk.getChiSquaredOverNDF();
		
		// https://github.com/cms-sw/cmssw/blob/master/DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrackRecoInfo.h
		CTPPSpixelLocalTrackReconstructionInfo  trackInfo = trk.getPixelTrackRecoInfo();
        ev_.ppstrk_RecoInfo[ev_.nppstrk]  = Short_t(trackInfo);

        /* UFSD only (2018)
        ev_.ppstrk_t[ev_.nppstrk] = trk.getTime();
        ev_.ppstrk_tUnc[ev_.nppstrk] = trk.getTimeUnc();
        */
        ev_.nppstrk++;
      }
  }

  //
  //PPS protons, loop over multi- and single-RP reco
  //
  ev_.nfwdtrk=0;
  for(size_t ip=0; ip<tokenRecoProtons_.size(); ip++){
    edm::Handle<vector<reco::ForwardProton>> recoProtons;
    iEvent.getByToken(tokenRecoProtons_[ip], recoProtons);
    if(recoProtons.isValid()){
      try{
        for (const auto & proton : *recoProtons)
          {
            if(!proton.validFit()) continue;

            CTPPSDetId detid( (*(proton.contributingLocalTracks().begin()))->getRPId() );

			//apply aperture cuts (Depending on the run number) https://twiki.cern.ch/twiki/bin/viewauth/CMS/TaggedProtonsGettingStarted#Fiducial_cuts
			if(proton.xi() > PPS_eff_->getXiHigh(detid.arm(),ev_.run,ev_.beamXangle)) continue;
			
			// keep only tracker information (no timing detectors)
			if(detid.rp()!=3) continue;
            ev_.fwdtrk_method[ev_.nfwdtrk]    = Short_t(proton.method());
			bool isPixel = (detid.station()==2), isMultiRP = (ev_.fwdtrk_method[ev_.nfwdtrk]==1);
			
			// Veto track: if not fiducial or is shifted (only for pixels)
			bool VetoTrack = false;
			if(isPixel || isMultiRP){ // check fiducial cut
				for (const auto &pr_tr : proton.contributingLocalTracks()) {
					CTPPSDetId rpId(pr_tr->getRPId());
					if(rpId.station()!=2) continue; // only pixels
					float x = pr_tr->getX(), y = pr_tr->getY();
					if(!PPS_eff_->InFiducial(detid.arm(), 2, ev_.run, x, y)) VetoTrack = true;
				}
			}
			if(VetoTrack) continue;
			
			// Reject shifted pixel tracks (https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsGettingStarted#Specific_features_and_warnings_f)
			CTPPSpixelLocalTrackReconstructionInfo pixtrackinfo = (*proton.contributingLocalTracks().begin())->getPixelTrackRecoInfo();
			bool _pixShift = !(pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::notShiftedRun || 
							pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::noShiftedPlanes ||
							pixtrackinfo == CTPPSpixelLocalTrackReconstructionInfo::invalid);
			//if(pixShift) skipTrack = true;
            ev_.fwdtrk_shifted[ev_.nfwdtrk]   = _pixShift ? 1 : 0;			


            ev_.fwdtrk_pot[ev_.nfwdtrk]       = 100*detid.arm()+10*detid.station()+detid.rp();
            ev_.fwdtrk_chisqnorm[ev_.nfwdtrk] = proton.normalizedChi2();
            ev_.fwdtrk_thetax[ev_.nfwdtrk]    = proton.thetaX();
            ev_.fwdtrk_thetay[ev_.nfwdtrk]    = proton.thetaY();
            ev_.fwdtrk_vx[ev_.nfwdtrk]        = proton.vx();
            ev_.fwdtrk_vy[ev_.nfwdtrk]        = proton.vy();
            ev_.fwdtrk_vz[ev_.nfwdtrk]        = proton.vz();
            ev_.fwdtrk_time[ev_.nfwdtrk]      = proton.time();
            ev_.fwdtrk_timeError[ev_.nfwdtrk] = proton.timeError();

            ev_.fwdtrk_xi[ev_.nfwdtrk]        = proton.xi();
            ev_.fwdtrk_xiSF[ev_.nfwdtrk]      = 1.;
            ev_.fwdtrk_xiError[ev_.nfwdtrk]   = proton.xiError();
            ev_.fwdtrk_t[ev_.nfwdtrk]         = proton.t();
			if(isMultiRP) nmultiprotons_[detid.arm()]++;
			
			// Get detector track x-y
			ev_.fwdtrk_FarX[ev_.nfwdtrk] = 0;
			ev_.fwdtrk_FarY[ev_.nfwdtrk] = 0;
			ev_.fwdtrk_NearX[ev_.nfwdtrk] = 0;
			ev_.fwdtrk_NearY[ev_.nfwdtrk] = 0;
			for (const auto &pr_tr : proton.contributingLocalTracks()) {
				CTPPSDetId rpId(pr_tr->getRPId());
				if(rpId.rp()==6) continue; // skip timing tracks
				if(rpId.station()){
					ev_.fwdtrk_FarX[ev_.nfwdtrk] = pr_tr->getX();
					ev_.fwdtrk_FarY[ev_.nfwdtrk] = pr_tr->getY();
				}
				else{
					  ev_.fwdtrk_NearX[ev_.nfwdtrk] = pr_tr->getX();
					  ev_.fwdtrk_NearY[ev_.nfwdtrk] = pr_tr->getY();
				}
			}			
            ev_.nfwdtrk++;
          }
      }
      catch(std::exception &e){
        cout << e.what() << endl;
      }
    }
  }// end loop over all protons (single-, multi-RP)
  
    
  // For 2017 check number of strip hits to count for truth/unsuff zero hits
  // Jan's presentation: https://indico.cern.ch/event/935869/
  // https://github.com/CTPPS/cmssw/blob/pps_event_cat/Validation/CTPPS/plugins/CTPPSEventCategoryPlotter.cc##L162
  edm::Handle<edm::DetSetVector<TotemRPRecHit>> hStripHits;
  iEvent.getByToken(tokenStripHits_,hStripHits);
  
  edm::Handle<edm::DetSetVector<TotemRPUVPattern>> hStripPatterns;  
  iEvent.getByToken(tokenStripPatterns_, hStripPatterns);
  
  if ( hStripHits.isValid() && hStripPatterns.isValid() ) { // true only for AOD
      	  
	  // Count number of planes with nHits>5:
	  unsigned int  n_too_full_u[2] = {0,0}, n_too_full_v[2] = {0,0};
	  for (const auto &dsRecHits : *hStripHits){
        TotemRPDetId planeId(dsRecHits.detId());
        unsigned int armId = planeId.arm();
        unsigned int planeIdx = planeId.plane();
		
		if(dsRecHits.size()>5){
		  if ((planeIdx % 2) == 1)
            n_too_full_u[armId]++;
          else
            n_too_full_v[armId]++;
	    }
      }
	  
	  // count number of fit patterns
	  unsigned int u_patterns[2] = {0,0}, v_patterns[2] = {0,0};
      for (const auto &dsPatterns : *hStripPatterns)
      {
        TotemRPDetId rpId(dsPatterns.detId());
        unsigned int armId = rpId.arm();

        for (const auto &pat : dsPatterns)
        {
          if (! pat.getFittable()) continue;

          if (pat.getProjection() == TotemRPUVPattern::projU)
            u_patterns[armId]++;
          if (pat.getProjection() == TotemRPUVPattern::projV)
            v_patterns[armId]++;
        }
      }
	  
	  // if nTracks>=1 and 0 multiRP assume events with nTracks>1 and increment the proton counter
	  for(int i=0;i<2;i++) {
		  bool is_suff = ((u_patterns[i]>= 1 || n_too_full_u[i]>=3) && (v_patterns[i]>= 1 || n_too_full_v[i]>=3));
		  if(nmultiprotons_[i]==0 &&  is_suff) nmultiprotons_[i]=2;
	  }
	  
	  //for(int i=0;i<2;i++){
	  //cout << "For rpId="<<i*100+3<<": n_too_full_u="<<n_too_full_u[i]<<" , u_patterns = " << u_patterns[i];
	  //cout << " n_too_full_v="<<n_too_full_v[i]<<" , v_patterns = " << v_patterns[i]<< ", suff  = " << is_suff[i] << endl;
      //}
	  
  } // end if(hStripHits.isValid() && hStripPatterns.isValid()))  
  
  if(ev_.MAXPROTONS<ev_.nfwdtrk || ev_.MAXPROTONS<ev_.nppstrk){
     cout << "ERROR: MAXPROTONS ("<<ev_.MAXPROTONS<<") is smaller than the N of RP hits/tracks ("<<ev_.nfwdtrk<<"/"<<ev_.nppstrk<<")."<<endl;
	 cout <<"\t\t... can cause memory leaks!!!"<<endl;
  }
	  
  //PF candidates (used for muon mini isolation) 
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfToken_,pfcands);
  	
  //
  //LEPTON SELECTION 
  ev_.nl=0; 
  
  //MUON SELECTION: cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  ev_.nrawmu=0;
  std::vector<Double_t> muStatUncReplicas(100,0);
  for (const pat::Muon &mu : *muons)
    {

      //raw muon info
      ev_.rawmu_pt[ev_.nrawmu]=(Short_t)mu.pt();
      ev_.rawmu_eta[ev_.nrawmu]=(Short_t)10*mu.eta();
      ev_.rawmu_phi[ev_.nrawmu]=(Short_t)10*mu.phi();
      ev_.rawmu_pid[ev_.nrawmu]= mu.selectors();
      ev_.nrawmu++;

      //apply correction
      float pt  = mu.pt();
      if(pt<2) continue; //no need to care about very low pt muons here... (corrections will tend to be meaningless)
      float eta = mu.eta();
      float phi = mu.phi();
      float q   = mu.charge();
      const reco::GenParticle * gen=mu.genLepton();

      //rochester corrections
      float sf(1.0),smearSeed(-1);
      float statUnc(0),zptUnc(0),ewkUnc(0),deltamUnc(0);
      if(iEvent.isRealData())
        {
          sf = muonRC_->kScaleDT(q, pt, eta, phi);
          for(int i=0; i<100; i++)
            muStatUncReplicas[i] = muonRC_->kScaleDT(q, pt, eta, phi, 1,i);
          statUnc = (1.0+TMath::StdDev(muStatUncReplicas.begin(),muStatUncReplicas.end()))/sf;
          zptUnc = muonRC_->kScaleDT(q, pt, eta, phi, 2)/sf;
          ewkUnc = muonRC_->kScaleDT(q, pt, eta, phi, 3)/sf;
          deltamUnc = muonRC_->kScaleDT(q, pt, eta, phi, 4)/sf;
        }
      else
        {
          smearSeed=gRandom->Rndm();
          int tlwm=(mu.innerTrack().isNonnull() ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0);
          sf = gen ?
            muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt()) :
            muonRC_->kSmearMC(q, pt, eta, phi, tlwm , smearSeed);
          for(int i=0; i<100; i++)
            muStatUncReplicas[i]=gen ?
              muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt(),1,i) :
              muonRC_->kSmearMC(q, pt, eta, phi, tlwm, smearSeed,1,i);
          statUnc=(1.0+TMath::StdDev(muStatUncReplicas.begin(),muStatUncReplicas.end()))/sf;
          zptUnc=(gen ?
                  muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt(),2) :
                  muonRC_->kSmearMC(q, pt, eta, phi, tlwm, smearSeed,2))/sf;
          ewkUnc=(gen ?
                  muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt(),3) :
                  muonRC_->kSmearMC(q, pt, eta, phi, tlwm, smearSeed,3))/sf;
          deltamUnc=(gen ?
                     muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt(),4) :
                     muonRC_->kSmearMC(q, pt, eta, phi, tlwm, smearSeed,4))/sf;
        }

      auto p4  = mu.p4() * sf;

      //kinematics
      bool passPt(p4.Pt() > 10);
      bool passEta(fabs(p4.Eta()) < 2.5);
      if(!passPt || !passEta) continue;

      //ID
      bool isLoose(muon::isLooseMuon(mu));
      if(!isLoose) continue;

      //save info
      ev_.l_isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.l_isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=13;
      ev_.l_g[ev_.nl]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=ev_.l_id[ev_.nl]) continue;
	  if(deltaR( mu.eta(),mu.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.l_g[ev_.nl]=ig;
	  break;
	}

      ev_.l_charge[ev_.nl]   = q;
      ev_.l_pt[ev_.nl]       = p4.Pt();
      ev_.l_eta[ev_.nl]      = p4.Eta();
      ev_.l_phi[ev_.nl]      = p4.Phi();
      ev_.l_mass[ev_.nl]     = p4.M();
      ev_.l_scaleUnc1[ev_.nl]= statUnc;
      ev_.l_scaleUnc2[ev_.nl]= zptUnc;
      ev_.l_scaleUnc3[ev_.nl]= ewkUnc;
      ev_.l_scaleUnc4[ev_.nl]= deltamUnc;
      ev_.l_highpt[ev_.nl]   = mu.tunePMuonBestTrack()->pt();
      ev_.l_scaleUnc5[ev_.nl]= mu.tunePMuonBestTrack()->ptError();
      ev_.l_mva[ev_.nl]      = 0;
      ev_.l_pid[ev_.nl]      = mu.selectors();
      ev_.l_chargedHadronIso[ev_.nl] = mu.pfIsolationR04().sumChargedHadronPt;
      ev_.l_miniIso[ev_.nl]  = getMiniIsolation(pfcands,&mu,0.05,0.2, 10., false);
      ev_.l_relIso[ev_.nl]   = (
				mu.pfIsolationR04().sumChargedHadronPt
				+ max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt)
				) / p4.Pt();
      ev_.l_ip3d[ev_.nl]    = -9999.;
      ev_.l_ip3dsig[ev_.nl] = -9999;
	  ev_.l_dz[ev_.nl] = -1;
	  // Trigger maching
	  //for(size_t i=0; i<triggersToUse_.size(); i++) 
	  //	  ev_.l_trigMatch[ev_.nl][i] = mu.triggered((triggersToUse_[i]+"*").c_str());
	  	  
      if(mu.innerTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::TrackRef>(mu.innerTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d[ev_.nl]    = ip3dRes.second.value();
	  ev_.l_ip3dsig[ev_.nl] = ip3dRes.second.significance();
	  
	  ev_.l_dz[ev_.nl] = fabs(mu.innerTrack()->dz(primVtx.position()));
	}
      ev_.nl++;

      if( p4.Pt()>25 && fabs(p4.Eta())<2.5 && isLoose) nrecleptons_++;
    }
	
  if(ev_.MAXRAWMU<ev_.nrawmu){
     cout << "ERROR: MAXRAWMU ("<<ev_.MAXRAWMU<<") is smaller than the N of raw muons ("<<ev_.nrawmu<<")."<<endl;
	 cout <<"\t\t... can cause memory leaks!!!"<<endl;
  }	

  // ELECTRON SELECTION: cf. https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electronToken_, electrons);
  for (const pat::Electron &e : *electrons)
    {

      float enSF(1.0);
      try{
        enSF=e.userFloat("ecalTrkEnergyPostCorr");
      }catch(...){
      }

      auto corrP4  = e.p4() * enSF / e.energy();

      //kinematics cuts
      bool passPt(corrP4.pt() > 15.0);
      bool passEta(fabs(corrP4.eta()) < 2.5 && (fabs(e.superCluster()->eta()) < 1.4442 || fabs(e.superCluster()->eta()) > 1.5660));
      if(!passPt || !passEta) continue;

      //full id+iso decisions
      bool isVeto( e.electronID("cutBasedElectronID-Fall17-94X-V"+EGIDVersion_+"-veto") );
      int vetoBits( e.userInt("cutBasedElectronID-Fall17-94X-V"+EGIDVersion_+"-veto")  );
      bool passVetoId( (vetoBits | 0xc0)== 0x3ff);  //mask isolation cuts and require all bits active
      bool isLoose( e.electronID("cutBasedElectronID-Fall17-94X-V"+EGIDVersion_+"-loose") );
      int looseBits( e.userInt("cutBasedElectronID-Fall17-94X-V"+EGIDVersion_+"-loose")  );
      bool passLooseId( (looseBits | 0xc0)== 0x3ff);  //mask isolation cuts and require all bits active
      bool isMedium( e.electronID("cutBasedElectronID-Fall17-94X-V"+EGIDVersion_+"-medium") );
      int mediumBits( e.userInt("cutBasedElectronID-Fall17-94X-V"+EGIDVersion_+"-medium")  );
      bool passMediumId( (mediumBits | 0xc0)== 0x3ff);  //mask isolation cuts and require all bits active
      bool isTight( e.electronID("cutBasedElectronID-Fall17-94X-V"+EGIDVersion_+"-tight") );
      int tightBits( e.userInt("cutBasedElectronID-Fall17-94X-V"+EGIDVersion_+"-tight") );
      bool passTightId( (tightBits | 0xc0)== 0x3ff);  //mask isolation cuts and require all bits active

      bool mvawp80(e.electronID("mvaEleID-Fall17-iso-V"+EGIDVersion_+"-wp80"));
      bool mvawp90(e.electronID("mvaEleID-Fall17-iso-V"+EGIDVersion_+"-wp90"));
      bool mvawploose(e.electronID("mvaEleID-Fall17-iso-V"+EGIDVersion_+"-wpLoose"));
      bool mvanonisowp80(e.electronID("mvaEleID-Fall17-noIso-V"+EGIDVersion_+"-wp80"));
      bool mvanonisowp90(e.electronID("mvaEleID-Fall17-noIso-V"+EGIDVersion_+"-wp90"));
      bool mvanonisowploose(e.electronID("mvaEleID-Fall17-noIso-V"+EGIDVersion_+"-wpLoose"));
      bool passHEEP(e.electronID("heepElectronID-HEEPV70"));

      //impact parameter cuts
      //see details in https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      bool passIpCuts(true); float dz = -1;
      if(e.gsfTrack().isNonnull())
	{
	  float dxy(fabs(e.gsfTrack()->dxy(primVtx.position())));
	  dz = (fabs(e.gsfTrack()->dz(primVtx.position())));
	  if(fabs(e.superCluster()->eta()) < 1.4442)
	    {
	      if(dxy>0.05 || dz>0.10) passIpCuts=false;
	    }
	  else
	    {
	      if(dxy>0.10 || dz>0.20) passIpCuts=false;
	    }
	}
      else
	{
	  passIpCuts=false;
	}

      //save the electron
      const reco::GenParticle * gen=e.genLepton();
      ev_.l_isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.l_isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=11;
      ev_.l_g[ev_.nl]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=ev_.l_id[ev_.nl]) continue;
	  if(deltaR( corrP4.eta(),corrP4.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.l_g[ev_.nl]=ig;
	  break;
	}
      ev_.l_mva[ev_.nl]=e.userFloat("ElectronMVAEstimatorRun2Fall17IsoV"+EGIDVersion_+"Values");
      ev_.l_mvaCats[ev_.nl]=e.userInt("ElectronMVAEstimatorRun2Fall17IsoV"+EGIDVersion_+"Categories");

      ev_.l_pid[ev_.nl]=0;
      ev_.l_pid[ev_.nl]= (passVetoId | (isVeto<<1)
			  | (passLooseId<<2) | (isLoose<<3)
			  | (passMediumId<<4) | (isMedium<<5)
			  | (passTightId<<6) | (isTight<<7)
			  | (passIpCuts<<8)
                          | (mvawp80<<9) | (mvawp90<<10) | (mvawploose<<11)
                          | (mvanonisowp80<<12) | (mvanonisowp90<<13) | (mvanonisowploose<<14)
                          | (passHEEP<<15)
			 );

      ev_.l_charge[ev_.nl]   = e.charge();
      ev_.l_dz[ev_.nl]   = dz;
      ev_.l_pt[ev_.nl]       = corrP4.pt();
      ev_.l_highpt[ev_.nl]   = corrP4.pt();
      ev_.l_eta[ev_.nl]      = corrP4.eta();
      ev_.l_phi[ev_.nl]      = corrP4.phi();
      ev_.l_mass[ev_.nl]     = corrP4.mass();
      ev_.l_scaleUnc1[ev_.nl] = 0.5*(e.userFloat("energyScaleStatUp")-e.userFloat("energyScaleStatDown"));
      ev_.l_scaleUnc2[ev_.nl] = 0.5*(e.userFloat("energyScaleGainUp")-e.userFloat("energyScaleGainDown"));
      ev_.l_scaleUnc3[ev_.nl] = 0.5*(e.userFloat("energyScaleSystUp")-e.userFloat("energyScaleSystDown"));
      ev_.l_scaleUnc4[ev_.nl] = 0.5*(e.userFloat("energySigmaUp")-e.userFloat("energySigmaDown"));
      ev_.l_scaleUnc5[ev_.nl] = 0.5*(e.userFloat("energySigmaPhiUp")-e.userFloat("energySigmaPhiDown"));
      ev_.l_scaleUnc6[ev_.nl] = 0.5*(e.userFloat("energySigmaRhoUp")-e.userFloat("energySigmaRhoDown"));
      ev_.l_scaleUnc7[ev_.nl] = 0.5*(e.userFloat("energyScaleUp")-e.userFloat("energyScaleDown"));
      ev_.l_miniIso[ev_.nl]  = getMiniIsolation(pfcands,&e,0.05, 0.2, 10., false);
      ev_.l_relIso[ev_.nl]   = (e.chargedHadronIso()+ max(0., e.neutralHadronIso() + e.photonIso()  - 0.5*e.puChargedHadronIso()))/corrP4.pt();
      ev_.l_chargedHadronIso[ev_.nl] = e.chargedHadronIso();
      ev_.l_ip3d[ev_.nl]     = -9999.;
      ev_.l_ip3dsig[ev_.nl]  = -9999;
      if(e.gsfTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::GsfTrackRef>(e.gsfTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d[ev_.nl]    = ip3dRes.second.value();
	  ev_.l_ip3dsig[ev_.nl] = ip3dRes.second.significance();
	}
      ev_.nl++;

      if( corrP4.pt()>25 && passEta && passLooseId ) nrecleptons_++;
    }
  if(ev_.MAXLEP<ev_.nl){
     cout << "ERROR: MAXLEP ("<<ev_.MAXLEP<<") is smaller than the N of leptons in the sample ("<<ev_.nl<<")."<<endl;
	 cout <<"\t\t... can cause memory leaks!!!"<<endl;
  }
  
  // PHOTON SELECTION: cf. https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
  ev_.ngamma=0;
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByToken(photonToken_, photons);
  for (const pat::Photon &g : *photons)
    {
      float enSF(1.0);
      try{
        enSF=g.userFloat("ecalEnergyPostCorr");
      }catch(...){
      }

      auto corrP4  = g.p4() * enSF / g.energy();

      //kinematics cuts
      bool passPt(corrP4.pt() > 30.0);
      float eta=corrP4.eta();
      bool passEta(fabs(eta) < 2.5 && (fabs(eta) < 1.4442 || fabs(eta) > 1.5660));
      if(!passPt || !passEta) continue;

      //full id+iso decisions
      int looseBits( g.photonID("cutBasedPhotonID-Fall17-94X-V"+EGIDVersion_+"-loose") );
      int mediumBits( g.photonID("cutBasedPhotonID-Fall17-94X-V"+EGIDVersion_+"-medium") );
      int tightBits( g.photonID("cutBasedPhotonID-Fall17-94X-V"+EGIDVersion_+"-tight") );
	  
      bool ismvawp80( g.photonID("mvaPhoID-RunIIFall17-v"+EGIDVersion_+"-wp80"));
      bool ismvawp90( g.photonID("mvaPhoID-RunIIFall17-v"+EGIDVersion_+"-wp90"));

      //save the photon
      const reco::GenParticle * gen=(const reco::GenParticle *)g.genPhoton();
      ev_.gamma_isPromptFinalState[ev_.ngamma] = gen ? gen->isPromptFinalState() : false;
      ev_.gamma_g[ev_.ngamma]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=22) continue;
	  if(deltaR( corrP4.eta(),corrP4.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.gamma_g[ev_.ngamma]=ig;
	  break;
	}

      ev_.gamma_mva[ev_.ngamma]=g.userFloat("PhotonMVAEstimatorRunIIFall17v"+EGIDVersion_+"Values");
      ev_.gamma_mvaCats[ev_.ngamma]=g.userInt("PhotonMVAEstimatorRunIIFall17v"+EGIDVersion_+"Categories");
      ev_.gamma_idFlags[ev_.ngamma]= g.passElectronVeto() | (g.hasPixelSeed()<<1) | (ismvawp80<<2) | (ismvawp90<<3);
      ev_.gamma_pid[ev_.ngamma]= ( (looseBits & 0x3ff)
                                   | ((mediumBits & 0x3ff)<<10)
                                   | ((tightBits & 0x3ff)<<20));
      ev_.gamma_pt[ev_.ngamma]  = corrP4.pt();
      ev_.gamma_eta[ev_.ngamma] = corrP4.eta();
      ev_.gamma_phi[ev_.ngamma] = corrP4.phi();
      ev_.gamma_scaleUnc1[ev_.ngamma] = 0.5*(g.userFloat("energyScaleStatUp")-g.userFloat("energyScaleStatDown"));
      ev_.gamma_scaleUnc2[ev_.ngamma] = 0.5*(g.userFloat("energyScaleGainUp")-g.userFloat("energyScaleGainDown"));
      ev_.gamma_scaleUnc3[ev_.ngamma] = 0.5*(g.userFloat("energyScaleSystUp")-g.userFloat("energyScaleSystDown"));
      ev_.gamma_scaleUnc4[ev_.ngamma] = 0.5*(g.userFloat("energySigmaUp")-g.userFloat("energySigmaDown"));
      ev_.gamma_scaleUnc5[ev_.ngamma] = 0.5*(g.userFloat("energySigmaPhiUp")-g.userFloat("energySigmaPhiDown"));
      ev_.gamma_scaleUnc6[ev_.ngamma] = 0.5*(g.userFloat("energySigmaRhoUp")-g.userFloat("energySigmaRhoDown"));
      ev_.gamma_scaleUnc7[ev_.ngamma] = 0.5*(g.userFloat("energyScaleUp")-g.userFloat("energyScaleDown"));
      ev_.gamma_chargedHadronIso[ev_.ngamma] = g.chargedHadronIso();
      ev_.gamma_neutralHadronIso[ev_.ngamma] = g.neutralHadronIso();
      ev_.gamma_photonIso[ev_.ngamma]        = g.photonIso();
      ev_.gamma_hoe[ev_.ngamma]              = g.hadTowOverEm();
      ev_.gamma_sieie[ev_.ngamma]            = g.full5x5_sigmaIetaIeta();
      ev_.gamma_r9[ev_.ngamma]               = g.full5x5_r9();
      ev_.ngamma++;
	  
      if(ev_.ngamma>ev_.MAXGAMMA) break;
      if( tightBits ) nrecphotons_++;
    }
  if(ev_.MAXGAMMA==ev_.ngamma){
     cout << "WARNING: MAXGAMMA ("<<ev_.MAXGAMMA<<") equal to stored N of photons in the sample ("<<ev_.ngamma<<")."<<endl;
	 cout <<"\t\t... check that the actuall number of photons is not larger!!!"<<endl;
  }
  
  // JETS
  ev_.nj=0;
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_,jets);
  JME::JetResolution jerResolution  = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
  JME::JetResolutionScaleFactor jerResolutionSF = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");
  std::vector< std::pair<const reco::Candidate *,int> > clustCands;
  for(auto j = jets->begin();  j != jets->end(); ++j)
    {
      //base kinematics
      if(j->pt()<10 || fabs(j->eta())>4.7) continue;

      //resolution corrections
      float jerSF[]={1.0,1.0,1.0};
      ev_.j_g[ev_.nj] = -1;
      if(!iEvent.isRealData())
        {
          //match to gen level jet/parton
          float genj_pt(0);
          const reco::Candidate *genParton = j->genParton();
          ev_.j_flav[ev_.nj]       = j->partonFlavour();
          ev_.j_hadflav[ev_.nj]    = j->hadronFlavour();
          ev_.j_pid[ev_.nj]        = genParton ? genParton->pdgId() : 0;
          for(int ig=0; ig<ev_.ng; ig++)
            {
              if(abs(ev_.g_id[ig])==11 || abs(ev_.g_id[ig])==13) continue;
              if(deltaR( j->eta(),j->phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
              genj_pt=ev_.g_pt[ig];
              ev_.j_g[ev_.nj]=ig;
              ev_.g_xbp[ig]  = genParton   ? ev_.g_xb[ig]*genj_pt/genParton->pt() : 0.;
              break;
            }

          //jet energy resolution
          JME::JetParameters jerParams = {{JME::Binning::JetPt, j->pt()},
                                          {JME::Binning::JetEta, j->eta()},
                                          {JME::Binning::Rho, rho}};
          float r = jerResolution.getResolution(jerParams);
          jerSF[0] = jerResolutionSF.getScaleFactor(jerParams);
          jerSF[1] = jerResolutionSF.getScaleFactor(jerParams, Variation::UP);
          jerSF[2] = jerResolutionSF.getScaleFactor(jerParams, Variation::DOWN);
          for(int i=0; i<3; i++) {
            //use stochasting smearing for unmatched jets
            if(genj_pt<=0)
              {
                float sigma = r * std::sqrt(std::max(float( pow(jerSF[i],2)- 1.0),float(0.)));
                jerSF[i] = std::max(float(1.0 + gRandom->Gaus(0, sigma)),float(0.));
              }
            else {
              float dPt = j->pt()-genj_pt;
              jerSF[i] = std::max(float(1.0 + (jerSF[i] - 1.) * dPt / j->pt()),float(0.));
            }
          }
          //make up/down variations relative
          if(jerSF[0]>0) { jerSF[1]/=jerSF[0]; jerSF[2]/=jerSF[0]; }		  
        }


      auto corrP4  = j->p4() * jerSF[0]; 
	  
	  // Skip jets with corrected PT < 10 GeV
      if(corrP4.pt()<10 ) continue;

      //jet id cf. for AK4CHS jets
      //2017 https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL#Preliminary_Recommendations_for
      float NHF  = j->neutralHadronEnergyFraction();
      float NEMF = j->neutralEmEnergyFraction();
      float CHF  = j->chargedHadronEnergyFraction();
      float MUF  = j->muonEnergyFraction();
      float CEMF = j->chargedEmEnergyFraction();
      float NumChargedParticles = j->chargedMultiplicity();
      float NumNeutralParticles = j->neutralMultiplicity();
      float NumConst = NumChargedParticles+NumNeutralParticles;
      float CHM = j->chargedMultiplicity();

      bool tightLepVeto(true),looseJetID(true);//,tightJetId(true);
      if(abs(j->eta())<2.6) {
        looseJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0);
        //tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99);
        tightLepVeto = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.80);
      }
      else if(abs(j->eta())<2.7) {
        looseJetID = (NHF<0.90 && NEMF<0.99 && CHM>0);
        //tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1);
        tightLepVeto = (NHF<0.90 && NEMF<0.99 && CHM>0 && MUF<0.8 && CEMF<0.80);
      }
      else if(abs(j->eta())<3.0) {
        looseJetID = (NEMF>0.01 && NEMF<0.99 && NumNeutralParticles>2);
        //tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
        tightLepVeto = (NEMF>0.01 && NEMF<0.99 && NumNeutralParticles>2);
      }
      else {
        looseJetID = (NHF>0.02 && NEMF<0.90 && NumNeutralParticles>10);
        tightLepVeto = (NHF>0.02 && NEMF<0.90 && NumNeutralParticles>10);
      }

      if(jetIdToUse_=="tightLepVeto") { if(!tightLepVeto) continue; }
      else { if(!looseJetID) continue; }

      //save jet
      ev_.j_area[ev_.nj]    = j->jetArea();
      ev_.j_jerUp[ev_.nj]   = jerSF[1];
      ev_.j_jerDn[ev_.nj]   = jerSF[2];
	  if(ev_.MAXJETSYS<jecCorrectionUncs_.size()){
         cout << "ERROR: MAXJETSYS ("<<ev_.MAXJETSYS<<") is smaller than the N jet syst. MC ("<<jecCorrectionUncs_.size()<<")."<<endl;
	     cout <<"\t\t... can cause memory leaks!!!"<<endl;
      }
      for(size_t iunc=0; iunc<jecCorrectionUncs_.size(); iunc++){
        jecCorrectionUncs_[iunc]->setJetPt(j->pt());
        jecCorrectionUncs_[iunc]->setJetEta(j->eta());
        ev_.j_jecUp[iunc][ev_.nj]=1.+jecCorrectionUncs_[iunc]->getUncertainty(true);
        jecCorrectionUncs_[iunc]->setJetPt(j->pt());
        jecCorrectionUncs_[iunc]->setJetEta(j->eta());
        ev_.j_jecDn[iunc][ev_.nj]=1.+jecCorrectionUncs_[iunc]->getUncertainty(false);
      }
      ev_.j_rawsf[ev_.nj]   = 1.0/jerSF[0];//j->pt()/corrP4.pt();
      //ev_.j_rawsf[ev_.nj]   = j->correctedJet("Uncorrected").pt()/j->pt();
      ev_.j_pt[ev_.nj]      = corrP4.pt();
      ev_.j_mass[ev_.nj]    = corrP4.mass();
      ev_.j_eta[ev_.nj]     = corrP4.eta();
      ev_.j_phi[ev_.nj]     = corrP4.phi();
      ev_.j_qg[ev_.nj]      = j->userFloat("QGTagger:qgLikelihood");
      ev_.j_pumva[ev_.nj]   = j->userFloat("pileupJetId:fullDiscriminant");
      ev_.j_id[ev_.nj]      = j->userInt("pileupJetId:fullId");
      ev_.j_csv[ev_.nj]     = j->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      ev_.j_deepcsv[ev_.nj] = j->bDiscriminator("pfDeepCSVJetTags:probb") + j->bDiscriminator("pfDeepCSVJetTags:probbb");
      //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
	  ev_.j_btag[ev_.nj]    = (ev_.j_deepcsv[ev_.nj]>0.4506);
      ev_.j_emf[ev_.nj]     = CEMF+NEMF;
	  
	  // jet momentum uncertainties (for exclusive ttbar analysis):
	  ev_.e_j_px[ev_.nj]  =  abs(cos(corrP4.phi()))*(abs(j->correctedJet("Uncorrected").pt()-j->pt()));
	  ev_.e_j_py[ev_.nj]  =  abs(sin(corrP4.phi()))*(abs(j->correctedJet("Uncorrected").pt()-j->pt()));
	  ev_.e_j_pz[ev_.nj]  =  (abs(j->correctedJet("Uncorrected").pt()-j->pt()))/(tan( 2.*atan(exp(-corrP4.eta())) ));

      //jet shape variables
      ev_.j_c2_00[ev_.nj]    = getC(2, 0.0, &(*j), true, 0.9);
      ev_.j_c2_02[ev_.nj]    = getC(2, 0.2, &(*j), true, 0.9);
      ev_.j_c2_05[ev_.nj]    = getC(2, 0.5, &(*j), true, 0.9);
      ev_.j_c2_10[ev_.nj]    = getC(2, 1.0, &(*j), true, 0.9);
      ev_.j_c2_20[ev_.nj]    = getC(2, 2.0, &(*j), true, 0.9);
      ev_.j_zg[ev_.nj]       = getZg(&(*j),true,0.9)[0];
      ev_.j_mult[ev_.nj]     = calcGA(0,0,&(*j),true,0.9);
      ev_.j_gaptd[ev_.nj]    = calcGA(0,2,&(*j),true,0.9);
      ev_.j_gawidth[ev_.nj]  = calcGA(1,1,&(*j),true,0.9);
      ev_.j_gathrust[ev_.nj] = calcGA(2,1,&(*j),true,0.9);
      //this function throws an exception in rare events
      try{
        ev_.j_tau32[ev_.nj]    = getTau(3,2,&(*j),true,0.9);
        ev_.j_tau21[ev_.nj]    = getTau(2,1,&(*j),true,0.9);
      }
      catch(...){
        ev_.j_tau32[ev_.nj]=-99;
        ev_.j_tau21[ev_.nj]=-99;
      }
      if( j->hasTagInfo("pfInclusiveSecondaryVertexFinder") )
	{
	  const reco::CandSecondaryVertexTagInfo *candSVTagInfo = j->tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
	  if( candSVTagInfo->nVertices() >= 1 )
	    {
	      math::XYZTLorentzVectorD vp4 = candSVTagInfo->secondaryVertex(0).p4();
	      ev_.j_vtxpx[ev_.nj]          = vp4.px();
	      ev_.j_vtxpy[ev_.nj]          = vp4.py();
	      ev_.j_vtxpz[ev_.nj]          = vp4.pz();
	      ev_.j_vtxmass[ev_.nj]        = vp4.mass();
	      ev_.j_vtxNtracks[ev_.nj]     = candSVTagInfo->nVertexTracks(0);
	      ev_.j_vtx3DVal[ev_.nj]       = candSVTagInfo->flightDistance(0).value();
	      ev_.j_vtx3DSig[ev_.nj]       = candSVTagInfo->flightDistance(0).significance();
	    }
	}
	  
	  // count reconstructed objects (used in the analysis)
	  if(ev_.j_pt[ev_.nj]>25 && abs(ev_.j_eta[ev_.nj])<2.5){
		  nrecjets_++;
		  if(ev_.j_btag[ev_.nj]) nrecbjets_++;
	  }
	  
	  // increment jet index
      ev_.nj++;	  

      //save all PF candidates central jet
      if(fabs(j->eta())>2.5) continue;
      for(size_t ipf=0; ipf<j->numberOfDaughters(); ipf++)
	{
	  const reco::Candidate *pf=j->daughter(ipf);
	  clustCands.push_back(std::pair<const reco::Candidate *,int>(pf,ev_.nj-1));
	}
    }
  if(ev_.MAXJET<ev_.nj){
     cout << "ERROR: MAXJET ("<<ev_.MAXJETSYS<<") is smaller than the N jets in the sample ("<<ev_.nj<<")."<<endl;
	 cout <<"\t\t... expect memory leaks!!!"<<endl;
  }
	  
  // MET
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  ev_.met_pt  = mets->at(0).pt();
  ev_.met_phi = mets->at(0).phi();
  ev_.met_sig = mets->at(0).significance();
  
  // Propogate JEC/JER corrections to MET: DONE AT THE ANALYSIS LEVEL
  //double px=mets->at(0).px(), py=mets->at(0).py();
  //for (int i=0;i<ev_.nj;i++){
//	  px+=ev_.j_pt[i]*cos(ev_.j_phi[i])*(1-ev_.j_rawsf[i]);
//	  py+=ev_.j_pt[i]*sin(ev_.j_phi[i])*(1-ev_.j_rawsf[i]);
  //}
  //ROOT::Math::PxPyPzM4D met(px,py,mets->at(0).pz(),0.0);
  //ev_.met_pt  = met.Pt();
  //ev_.met_phi = met.Phi();


  for(size_t i=0; i<ev_.MAXMETSYS; i++){
    ev_.met_ptShifted[i]  = mets->at(0).shiftedPt(pat::MET::METUncertainty(i));
    ev_.met_phiShifted[i] = mets->at(0).shiftedPhi(pat::MET::METUncertainty(i));
  }
  
  // MET errors (for exclusive ttbar analysis)
  ev_.e_met_px  = mets->at(0).getSignificanceMatrix()[0][0];
  ev_.e_met_py  = mets->at(0).getSignificanceMatrix()[1][1];
  ev_.e_met_pxpy= mets->at(0).getSignificanceMatrix()[0][1];
  
  //MET filter bits https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
  ev_.met_filterBits=0;
  edm::Handle<edm::TriggerResults> h_metFilters;
  iEvent.getByToken(metFilterBits_, h_metFilters);
  std::vector<string> metFilterNames;
  Service<service::TriggerNamesService> mfns;
  mfns->getTrigPaths(*h_metFilters,metFilterNames);
  for (unsigned int i=0; i< h_metFilters->size(); i++)
    {
      if( !(*h_metFilters)[i].accept() ) continue;
      for(size_t itrig=0; itrig<metFiltersToUse_.size(); itrig++)
	{
	  //cout << metFiltersToUse_[itrig] << " pass =  " << !(metFilterNames[i].find(metFiltersToUse_[itrig])==string::npos) << endl;
	  if (metFilterNames[i].find(metFiltersToUse_[itrig])==string::npos) continue;
	  ev_.met_filterBits |= (1<<itrig);
	  //cout << metFiltersToUse_[itrig] << " pass! " << endl;
	}
    }
  
  // Bad Charged Hadron and Bad Muon Filters cannot be accessed directly:
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters
  //try{
  //  edm::Handle<bool> ifilterbadPFMuon;
  //  iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
  //  bool filterbadPFMuon = *ifilterbadPFMuon;
  //	cout << "Check muon, answer = " << filterbadPFMuon << endl;
  //  ev_.met_filterBits |= (filterbadPFMuon<<(metFiltersToUse_.size()));
  //}
  //catch(...){
  //}

  try{
    edm::Handle<bool> ifilterbadPFMuonDz;
    iEvent.getByToken(BadPFMuonDzFilterToken_, ifilterbadPFMuonDz);
    bool filterbadPFMuonDz = *ifilterbadPFMuonDz;
    ev_.met_filterBits |= (filterbadPFMuonDz<<(metFiltersToUse_.size()));
  }
  catch(...){
  }
  
  //try{
  //  edm::Handle<bool> ifilterbadChCand;
  //  iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand);
  //  bool  filterbadChCandidate = *ifilterbadChCand;
  //  ev_.met_filterBits |= (filterbadChCandidate<<metFiltersToUse_.size()+2);
  //}
  //catch(...){
  //}	
  

  //PF candidates
  LorentzVector vtxPt[8]; 
  ev_.nchPV=0; ev_.sumPVChPt=0; ev_.sumPVChPz=0; ev_.sumPVChHt=0;
  ev_.ntrk=0;
  for(int i=0; i<8; i++){
	vtxPt[i].SetXYZT(0,0,0,0);
	ev_.nchPV_v[i]=0;
	ev_.sumPVChPt_v[i]=0;
	ev_.sumPVChPz_v[i]=0;
	ev_.sumPVChHt_v[i]=0;
	
    ev_.nPFCands[i]=0;
    ev_.sumPFHt[i]=0;
    ev_.sumPFEn[i]=0;
    ev_.sumPFPz[i]=0;
    ev_.nPFChCands[i]=0;
    ev_.sumPFChHt[i]=0;
    ev_.sumPFChEn[i]=0;
    ev_.sumPFChPz[i]=0;
  }
  for(auto pf = pfcands->begin();  pf != pfcands->end(); ++pf)
    {
      int ieta(-1);
      if(pf->eta()>-4.7) ieta=0;
      if(pf->eta()>-3)   ieta=1;
      if(pf->eta()>-2.5) ieta=2;
      if(pf->eta()>-1.5) ieta=3;
      if(pf->eta()>0)    ieta=4;
      if(pf->eta()>1.5)  ieta=5;
      if(pf->eta()>2.5)  ieta=6;
      if(pf->eta()>3.0)  ieta=7;
      if(pf->eta()>4.7)  ieta=-1;
      if(ieta<0) continue;
      ev_.nPFCands[ieta]++;
      ev_.sumPFHt[ieta] += pf->pt();
      ev_.sumPFEn[ieta] += pf->energy();
      ev_.sumPFPz[ieta] += pf->pz();
      if(pf->charge()!=0){
        ev_.nPFChCands[ieta]++;
        ev_.sumPFChHt[ieta] += pf->pt();
        ev_.sumPFChEn[ieta] += pf->energy();
        ev_.sumPFChPz[ieta] += (pf->pz());
        bool passChargeSel(pf->pt()>0.9 && fabs(pf->eta())<2.5); // split 2.1 and 2.5
        const pat::PackedCandidate::PVAssoc pvassoc=pf->fromPV(); 
		const pat::PackedCandidate::PVAssoc pvassoc2=pf->fromPV(_second_vertex_index); 
		int _bin;
        if(passChargeSel && pvassoc>=pat::PackedCandidate::PVTight){
          ev_.nchPV++;
          ev_.sumPVChPz+=(pf->pz());
          ev_.sumPVChHt+=pf->pt();
		  
		  // Add extra PV variables
		  _bin = 0;
		  ev_.nchPV_v[_bin]++;
		  ev_.sumPVChPz_v[_bin]+=(pf->pz());
		  ev_.sumPVChHt_v[_bin]+=pf->pt();
		  vtxPt[_bin]+=pf->p4();
		  if(fabs(pf->eta())<2.1) {_bin = 1;
			  ev_.nchPV_v[_bin]++;
			  ev_.sumPVChPz_v[_bin]+=(pf->pz());
			  ev_.sumPVChHt_v[_bin]+=pf->pt();
			  vtxPt[_bin]+=pf->p4();
			  if(ev_.ntrk<ev_.MAXTRACKS){
			    ev_.track_pt[ev_.ntrk] = pf->pt();
			    ev_.track_eta[ev_.ntrk] = pf->eta();
			    ev_.track_phi[ev_.ntrk] = pf->phi();
			    ev_.ntrk++;
			  }			  
		      if(pvassoc>pat::PackedCandidate::PVTight) {_bin=3;
			    ev_.nchPV_v[_bin]++;
			    ev_.sumPVChPz_v[_bin]+=(pf->pz());
			    ev_.sumPVChHt_v[_bin]+=pf->pt();
			    vtxPt[_bin]+=pf->p4();
		      }
		  }
		  if(pvassoc>pat::PackedCandidate::PVTight) {_bin=2;
			    ev_.nchPV_v[_bin]++;
			    ev_.sumPVChPz_v[_bin]+=(pf->pz());
			    ev_.sumPVChHt_v[_bin]+=pf->pt();
			    vtxPt[_bin]+=pf->p4();
		  }
		}
		if(passChargeSel && pvassoc2>=pat::PackedCandidate::PVTight){
		  _bin = 4;
		  ev_.nchPV_v[_bin]++;
		  ev_.sumPVChPz_v[_bin]+=(pf->pz());
		  ev_.sumPVChHt_v[_bin]+=pf->pt();
		  vtxPt[_bin]+=pf->p4();
		  if(fabs(pf->eta())<2.1) {_bin = 4+1;
			  ev_.nchPV_v[_bin]++;
			  ev_.sumPVChPz_v[_bin]+=(pf->pz());
			  ev_.sumPVChHt_v[_bin]+=pf->pt();
			  vtxPt[_bin]+=pf->p4();
		      if(pvassoc2>pat::PackedCandidate::PVTight) {_bin=4+3;
			    ev_.nchPV_v[_bin]++;
			    ev_.sumPVChPz_v[_bin]+=(pf->pz());
			    ev_.sumPVChHt_v[_bin]+=pf->pt();
			    vtxPt[_bin]+=pf->p4();
		      }
		  }
		  if(pvassoc2>pat::PackedCandidate::PVTight) {_bin=4+2;
			    ev_.nchPV_v[_bin]++;
			    ev_.sumPVChPz_v[_bin]+=(pf->pz());
			    ev_.sumPVChHt_v[_bin]+=pf->pt();
			    vtxPt[_bin]+=pf->p4();
		  }
		}		
      }
    }
  ev_.sumPVChPt=vtxPt[0].pt();
  for(int i=0; i<8; i++)
    ev_.sumPVChPt_v[i]=vtxPt[i].pt();
}

//cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
bool MiniAnalyzer::isSoftMuon(const reco::Muon & recoMu,const reco::Vertex &vertex)
{

  bool isGood(muon::isGoodMuon(recoMu, muon::TMOneStationTight));
  bool passLayersWithMeas(recoMu.innerTrack().isNonnull()
			  && recoMu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
			  && recoMu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 );
  bool matchesVertex(recoMu.innerTrack().isNonnull()
		     && fabs(recoMu.innerTrack()->dxy(vertex.position())) < 0.3
		     && fabs(recoMu.innerTrack()->dz(vertex.position())) < 20. );
  return (isGood && passLayersWithMeas && matchesVertex);
}

//cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Standard_MediumID_to_be_used_wit
bool MiniAnalyzer::isMediumMuon2016ReReco(const reco::Muon & recoMu)
{
  bool goodGlob = recoMu.isGlobalMuon() &&
    recoMu.globalTrack()->normalizedChi2() < 3 &&
    recoMu.combinedQuality().chi2LocalPosition < 12 &&
    recoMu.combinedQuality().trkKink < 20;
  bool isMedium = muon::isLooseMuon(recoMu) &&
    recoMu.innerTrack()->validFraction() > 0.8 &&
    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451);
  return isMedium;
}



// ------------ method called for each event  ------------
void MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get beam-crossing angle and LHC conditions
  try{
    edm::ESHandle<LHCInfo> hLHCInfo;
    std::string lhcInfoLabel("");
    iSetup.get<LHCInfoRcd>().get(lhcInfoLabel, hLHCInfo);
    if(hLHCInfo.isValid()){
      ev_.beamXangle=hLHCInfo->crossingAngle();
      ev_.instLumi=hLHCInfo->instLumi(); // instantaneous luminosity in ub-1
      ev_.betaStar=hLHCInfo->betaStar();
      ev_.fill=hLHCInfo->fillNumber();
    }
  }
  catch(...){
    ev_.beamXangle=0;
    ev_.instLumi=0;
    ev_.betaStar=0;
    ev_.fill=0;
  }

  histContainer_["counter"]->Fill(0);
  ngleptons_=0;   ngphotons_=0;
  nrecleptons_=0; nrecphotons_=0; nmultiprotons_[0]=nmultiprotons_[1]=0;
  nrecjets_=0; nrecbjets_=0;
  ev_.g_nw=0; ev_.ng=0; ev_.ngtop=0;
  ev_.nl=0; ev_.ngamma=0; ev_.nj=0; ev_.nfwdtrk=0; ev_.nrawmu=0;

  //analyze the event
  ev_.isData  = iEvent.isRealData();
  ev_.run     = ev_.isData ? iEvent.id().run() : runNumber_;
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event();

  if(!ev_.isData) genAnalysis(iEvent,iSetup);
  recAnalysis(iEvent,iSetup);
  

  //save event if at least one object at gen or reco level
  if(!saveTree_) return;
  // Define some filters, otherwise, you might store a lot of empty events (check returns in recAnalysis)
  if(applyFilt_){
	//if (FilterType_.find("ttbar")!=std::string::npos) if(nrecbjets_<2 || nrecjets_<4 || nrecleptons_==0) return;
      if (FilterType_.find("ttbar")!=std::string::npos) {
		  // skim (nJ>=4 and nL>0) OR (nL>1)
		  if((nrecjets_<4 || nrecleptons_==0) && (nrecleptons_<2)) return;
		  //if(!ev_.isData && nrecbjets_<2) return;
		  if(!ev_.isData && nrecbjets_<1) return;
		  if(!ev_.isData && nrecjets_>=4 && nrecbjets_<2) return;
	  }
	  if (FilterType_.find("dilep")!=std::string::npos) if(nrecleptons_<2) return;
	  if (FilterType_.find("lowmu")!=std::string::npos) {
		  if(ev_.nl==0 && ev_.nj==0 && nrecphotons_==0) return;
		  if(ev_.nl==0 && nrecphotons_==0 && ev_.nj>0 && ev_.j_pt[0]<100) return;
	  }
	// data - skim on event w/o forward protons but save the event count
	histContainer_["counter"]->Fill(2);
	if(ev_.isData) histContainer_["RPcount"]->Fill(nmultiprotons_[0],nmultiprotons_[1]);
	if(nrecbjets_!=0) histContainer_["counter"]->Fill(3);
	
    if (ev_.isData){ 
		if (FilterType_.find("ttbar")!=std::string::npos)
			if( nmultiprotons_[0]!=1 ||  nmultiprotons_[1]!=1 ) return;
		//if (FilterType_.find("dilep")!=std::string::npos)
		//	if( nmultiprotons_[0]!=1 &&  nmultiprotons_[1]!=1 ) return;
	}
  }
  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
MiniAnalyzer::beginJob(){
}

//
void
MiniAnalyzer::endRun(const edm::Run& iRun,
		     const EventSetup& iSetup)
{
  //try{

    //cout << "[MiniAnalyzer::endRun]" << endl;
	
	// Following lines list the generator weights as described in 
	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#Retrieving_the_weights
	
    //edm::Handle<LHERunInfoProduct> lheruninfo;
    //typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    //iRun.getByToken(generatorRunInfoToken_, lheruninfo );

    //LHERunInfoProduct myLHERunInfoProduct = *(lheruninfo.product());
	
	// Print all weights and corresponding integers
	//for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    //  std::cout << "tag="<<iter->tag() << std::endl;
    //  std::vector<std::string> lines = iter->lines();
    //  for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
	//	if(lines.at(iLine)=="") continue;
		//if(lines.at(iLine).find("weightgroup")==std::string::npos) continue;
    //    std::cout << lines.at(iLine) << std::endl;
    //  }
    //}
/*
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin();
	 iter!=myLHERunInfoProduct.headers_end();
	 iter++)
      {
	std::string tag("generator");
	if(iter->tag()!="") tag+="_"+iter->tag();

	std::vector<std::string> lines = iter->lines();
	std::vector<std::string> prunedLines;
	for (unsigned int iLine = 0; iLine<lines.size(); iLine++)
	  {
	    if(lines.at(iLine)=="") continue;
	    if(lines.at(iLine).find("weightgroup")!=std::string::npos) continue;
	    prunedLines.push_back( lines.at(iLine) );
	  }

	if(histContainer_.find(tag)==histContainer_.end())
	  {
	    std::cout << "Starting histo for " << tag << std::endl;
	    histContainer_[tag]=fs->make<TH1F>(tag.c_str(),tag.c_str(),prunedLines.size(),0,prunedLines.size());
	  }
	for (unsigned int iLine = 0; iLine<prunedLines.size(); iLine++)
	  histContainer_[tag]->GetXaxis()->SetBinLabel(iLine+1,prunedLines.at(iLine).c_str());
      }
	  */
  //}
  //catch(std::exception &e){
  //  std::cout << e.what() << endl
//	      << "Failed to retrieve LHERunInfoProduct" << std::endl;
  //}
  

}

//-------------
//cf. https://twiki.cern.ch/twiki/bin/view/CMS/MiniIsolationSUSY
float MiniAnalyzer::getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
				     const reco::Candidate* ptcl,
				     float r_iso_min, float r_iso_max, float kt_scale,
				     bool charged_only)
{

    if (ptcl->pt()<5.) return 99999.;

    float deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    float iso_nh(0.), iso_ch(0.), iso_ph(0.), iso_pu(0.);
    float ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    float r_iso = (float)TMath::Max((float)r_iso_min,
				    (float)TMath::Min((float)r_iso_max, (float)(kt_scale/ptcl->pt())));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;

      float dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;

      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    float iso(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      iso -= 0.5*iso_pu;
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }
    iso = iso/ptcl->pt();

    return iso;
}

// ------------ method called once each job just after ending the event loop  ------------
void
MiniAnalyzer::endJob()
{
  if(saveTree_) std::cout << "store tree with " << tree_->GetEntries() << " entries."<<endl<<"[MiniAnalyzer::endJob]"<<endl;
  else std::cout << "[MiniAnalyzer::endJob]" << endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
