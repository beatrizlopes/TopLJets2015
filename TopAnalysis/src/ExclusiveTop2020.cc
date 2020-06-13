#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveTop2020.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include "TMath.h"

using namespace std;

#define ADDVAR(x,name,t,tree) tree->Branch(name,x,TString(name)+TString(t))

void RunExclusiveTop2020(const TString in_fname,
                      TString outname,
                      TH1F *normH, 
                      TH1F *genPU, 
                      TString era,
                      bool debug,
					  bool skimtree) 
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  float EBEAM(6500);
  float XIMIN(0.03);
  float XIMAX(0.2);
  bool SKIMME(false);
  
  //////////////////////////////////
  // Constants used in ttbar reco //
  //////////////////////////////////
  float m_TOP = 173.1;
  float m_W   =  80.379;
  
  //const char* CMSSW_BASE = getenv("CMSSW_BASE");
  MiniEvent_t ev;  

  //preselection cuts to apply
  float minLeptonPt(30);
  float minJetPt(15);
  size_t minJetMultiplicity(4);
  size_t minBJetMultiplicity(2);

  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;
  
  //CORRECTIONS: LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  std::map<Int_t,Float_t> lumiPerRun=lumi.lumiPerRun();
  
  //CORRECTIONS: LEPTON EFFICIENCIES
  //EfficiencyScaleFactorsWrapper lepEffH(in_fname.Contains("Data13TeV"),era);

  //CORRECTIONS: L1-prefire 
  L1PrefireEfficiencyWrapper l1PrefireWR(in_fname.Contains("Data13TeV"),era);
  
  
  //CORRECTIONS: B-TAG CALIBRATION
  BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);

  //READ TREE FROM FILE
  TFile *f = TFile::Open(in_fname);
  if(f==NULL || f->IsZombie()) {
    cout << "Corrupted or missing file " << in_fname << endl;
    return;
  }
  TH1 *counter=(TH1 *)f->Get("analysis/counter");
  if(!counter) {cout << "Corrupted or missing counter: \"analysis/counter\" " << endl;return;}
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  if(!triggerList) {cout << "Corrupted or missing triggerList: \"analysis/triggerList\" " << endl;return;}
  TTree *t = (TTree*)f->Get("analysis/tree");
  if(!t) {cout << "Corrupted or missing tree: \"analysis/tree\" " << endl;return;}
  attachToMiniEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = min(1000,nentries); //restrict number of entries for testing
  t->GetEntry(0);
  std::cout << "...producing " << outname << " from " << nentries << " events" << std::endl;  
  bool isData = ev.isData;

  //PREPARE OUTPUT (BOOK SOME HISTOGRAMS)
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  
  //Variables (bools, ints, floats):
  TString bvars[]={"hasMTrigger","hasETrigger","hasEMTrigger","hasMETrigger","hasMMTrigger","hasEETrigger","hasSLT","hasDLT"};
  std::map<TString,bool> boutVars;
  for(size_t i=0; i<sizeof(bvars)/sizeof(TString); i++) boutVars[bvars[i]]=false;
  
  TString ivars[]={"nl","nbj","nlj","nxip","nxin"};
  std::map<TString,Int_t> ioutVars;
  for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++) ioutVars[ivars[i]]=0;
  
  TString fvars[]={"gen_mtt","gen_Ytt","mpp","Ypp","xip","xin","xip_truth","xin_truth",
                    "HT", "chi2min","mtt","Ytt","pTtt","dPhitt"};
  std::map<TString,Float_t> foutVars;
  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++) foutVars[fvars[i]]=0.;
  
  int pixel_pos_n(0), pixel_neg_n(0), pixel_n(0);
  int strip_pos_n(0), strip_neg_n(0), strip_n(0);
  int multi_n(0), mass_pixel_n(0);
  float pixel_pos_xi[ev.MAXPROTONS], pixel_neg_xi[ev.MAXPROTONS];
  float strip_pos_xi[ev.MAXPROTONS], strip_neg_xi[ev.MAXPROTONS]; 
  float mass_pixel_all[200];
  float multi_pos_xi,multi_neg_xi, multi_pos_xi_eff,multi_neg_xi_eff;
  int ttbar_assignment[6], ttbar_idx[4], top_type[2];
  float ttbar_assignment_res[6];
  float tops_jet_pt[6],tops_jet_eta[6],tops_jet_phi[6];
  
  int jet_n(0);
  float jet_pt[ev.MAXJET];

  TTree *outT=new TTree("tree","tree");
  if(skimtree) {
	  
	outT->Branch("run",&ev.run,"run/i");
	outT->Branch("event",&ev.event,"event/l");
	outT->Branch("lumi",&ev.lumi,"lumi/i");
	outT->Branch("nvtx",&ev.nvtx,"nvtx/I");
    if(t->FindBranch("nchPV")) outT->Branch("nchPV",&ev.nchPV,"nchPV/I");
	ADDVAR(&ev.beamXangle,"beamXangle","/F",outT);
	
	
	if(t->FindBranch("met_pt")){
	ADDVAR(&ev.met_pt,"met_pt","/F",outT);
	ADDVAR(&ev.met_phi,"met_phi","/F",outT);
	ADDVAR(&ev.met_sig,"met_sig","/F",outT);
	}

	outT->Branch("jet_n",&jet_n,"jet_n/I");
	outT->Branch("jet_pt",jet_pt,"jet_pt[jet_n]/F");
	
	// fill custom variables
	for(size_t i=0; i<sizeof(bvars)/sizeof(TString); i++){
	  ADDVAR(&(boutVars[bvars[i]]),bvars[i],"/O",outT);
	}

	for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++){
	  ADDVAR(&(ioutVars[ivars[i]]),ivars[i],"/I",outT);
	}
	
	for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++){
	  ADDVAR(&(foutVars[fvars[i]]),fvars[i],"/F",outT);
	}	
	
	//PPS info
    outT->Branch("multi_n",&multi_n,"multi_n/I");
    outT->Branch("multi_pos_xi",&multi_pos_xi,"multi_pos_xi/F");
    outT->Branch("multi_neg_xi",&multi_neg_xi,"multi_neg_xi/F");
    outT->Branch("multi_pos_xi_eff",&multi_pos_xi_eff,"multi_pos_xi_eff/F");
    outT->Branch("multi_neg_xi_eff",&multi_neg_xi_eff,"multi_neg_xi_eff/F");
    outT->Branch("pixel_pos_n",&pixel_pos_n,"pixel_pos_n/I");
    outT->Branch("pixel_neg_n",&pixel_neg_n,"pixel_neg_n/I");
    outT->Branch("pixel_n",&pixel_n,"pixel_n/i");
    outT->Branch("strip_pos_n",&strip_pos_n,"strip_pos_n/I");
    outT->Branch("strip_neg_n",&strip_neg_n,"strip_neg_n/I");
    outT->Branch("strip_n",&strip_n,"strip_n/I");	
	outT->Branch("pixel_pos_xi",pixel_pos_xi,"pixel_pos_xi[pixel_pos_n]/F");
	outT->Branch("pixel_neg_xi",pixel_neg_xi,"pixel_neg_xi[pixel_neg_n]/F");
	outT->Branch("strip_pos_xi",strip_pos_xi,"strip_pos_xi[strip_pos_n]/F");
	outT->Branch("strip_neg_xi",strip_neg_xi,"strip_neg_xi[strip_neg_n]/F");
    outT->Branch("mass_pixel_n",&mass_pixel_n,"mass_pixel_n/I");
	outT->Branch("mass_pixel_all",mass_pixel_all,"mass_pixel_all[mass_pixel_n]/F");
	
	// Jet truth assigment (blep, lep, bhad, jet1, jet2)
	if(t->FindBranch("ngtop")){
		outT->Branch("ttbar_assignment",ttbar_assignment,"ttbar_assignment[6]/I"); // index of jet or lepton
		outT->Branch("ttbar_assignment_res",ttbar_assignment_res,"ttbar_assignment_res[6]/F"); // index of jet or lepton
		outT->Branch("tops_jet_pt",tops_jet_pt,"tops_jet_pt[6]/F"); // index of jet or lepton
		outT->Branch("tops_jet_eta",tops_jet_eta,"tops_jet_eta[6]/F"); // index of jet or lepton
		outT->Branch("tops_jet_phi",tops_jet_phi,"tops_jet_phi[6]/F"); // index of jet or lepton
		outT->Branch("top_type",top_type,"top_type[2]/I");
	}
	// ttbar reconstrcted indices
	outT->Branch("ttbar_idx",ttbar_idx,"ttbar_idx[4]/I"); // index of jet or lepton
	
	outT->SetDirectory(fOut);
  } // end if SkimTree
  
  //BOOK HISTOGRAMS  
  HistTool ht;
  ht.setNsyst(0);
  //ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",50,0,100));
  //ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",10,0,10));
  //ht.addHist("mlb",          new TH1F("mlb",         ";m(l,b) [GeV];Events",20,0,250)); 
  //ht.addHist("nprotons",     new TH1F("nprotons",    ";Proton multiplicity; Events",6,0,6) );
  //ht.addHist("csi",          new TH1F("csi",         ";#xi = #deltap/p; Events",50,0,0.3) );
  //ht.addHist("x",            new TH1F("x",           ";x  [cm]; Events",50,0,25) );
  ht.addHist("ratevsrun",    new TH1F("ratevsrun",   ";Run number; #sigma [pb]",int(lumiPerRun.size()),0,float(lumiPerRun.size())));
  int i=0;
  for(auto key : lumiPerRun) {
    i++;
    ht.getPlots()["ratevsrun"]->GetXaxis()->SetBinLabel(i,Form("%d",key.first));
  }
  // normalization and event count
  ht.addHist("norm",     new TH1F("h_norm",    ";Category; Events",3,0,3) );
  ht.getPlots()["norm"]->GetXaxis()->SetBinLabel(1,"SumWeighted");
  ht.getPlots()["norm"]->GetXaxis()->SetBinLabel(2,"SumUnweighted");
  ht.getPlots()["norm"]->GetXaxis()->SetBinLabel(3,"inc");
  ht.getPlots()["norm"]->SetBinContent(1,counter->GetBinContent(1));
  ht.getPlots()["norm"]->SetBinContent(2,counter->GetBinContent(2));
  ht.getPlots()["norm"]->SetBinContent(3,nentries);
  
  std::cout << "initialization done" << std::endl;

  //EVENT SELECTION WRAPPER (GETS LISTS OF PHYSICS OBJECTS FROM THE INPUT)
  SelectionTool selector(in_fname, false, triggerList);
  
  //EVENT LOOP
  //select el/mu + >=4 jets events triggered by a single muon trigger
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)(iev+1)/(float)nentries); fflush(stdout); }
		
	  // Reset vars
      // ------------------------------------------------------------------------------------------------------- //
	  if(debug) cout << "Event #"<<iev<<" Reset vars" << endl;
	  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++) foutVars[fvars[i]]=0.;
	  for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++) ioutVars[ivars[i]]=0;
	  for(int i=0;i<4;i++) ttbar_idx[i]=-1;
	  for(int i=0;i<6;i++) {ttbar_assignment[i]=-1; ttbar_assignment_res[i]=tops_jet_pt[i]=tops_jet_eta[i]=tops_jet_phi[i]=0;}
      top_type[0] = top_type[1] = 0;
	  pixel_pos_n=0;pixel_neg_n=0;
	  strip_pos_n=0;strip_neg_n=0;
	  multi_n=0; multi_pos_xi=0; multi_neg_xi=0; multi_pos_xi_eff=0; multi_neg_xi_eff=0;
      for(int irp=0; irp<ev.MAXPROTONS; irp++) {
		pixel_pos_xi[irp]=0;pixel_neg_xi[irp]=0;
		strip_pos_xi[irp]=0;strip_neg_xi[irp]=0;
	  }
	  for(int irp2=0; irp2<200;irp2++) mass_pixel_all[irp2]=0;
      // ------------------------------------------------------------------------------------------------------- //

      //trigger
      boutVars["hasMTrigger"] = false;
      if(era.Contains("2016")) boutVars["hasMTrigger"]=(selector.hasTriggerBit("HLT_IsoMu24_v", ev.triggerBits) );     
      if(era.Contains("2017")) boutVars["hasMTrigger"]=(selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits) ); 
	  boutVars["hasETrigger"]=(selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v",                                  ev.triggerBits) || 
                               selector.hasTriggerBit("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v",                     ev.triggerBits) ||
							   selector.hasTriggerBit("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v", ev.triggerBits));
	  
      boutVars["hasMMTrigger"]=(selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev.triggerBits) ||
								selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                  ev.triggerBits) ||
								selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",          ev.triggerBits));
	  boutVars["hasEETrigger"]=(selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev.triggerBits) ||
								selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev.triggerBits) );
	  boutVars["hasEMTrigger"]=(selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits) );
	  boutVars["hasMETrigger"]=(selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
								selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits));
      
	  boutVars["hasSLT"] = boutVars["hasMTrigger"] || boutVars["hasETrigger"];
	  boutVars["hasDLT"] = (boutVars["hasMMTrigger"]||boutVars["hasEETrigger"]||boutVars["hasEMTrigger"]||boutVars["hasMETrigger"]);
      
	  if(debug) cout << "runNumber #"<<ev.run<<" pass trigger (hasSLT) = "  << boutVars["hasSLT"] <<endl;
	  if(!boutVars["hasSLT"] && SKIMME) continue;

      //lepton selection
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);     
      SelectionTool::QualityFlags muId(SelectionTool::TIGHT); //TIGHT
      leptons = selector.selLeptons(leptons,muId,SelectionTool::MVA90,minLeptonPt,2.1);
      bool passLeptons(leptons.size()==1);
	  if(debug) cout << "passLeptons = "<<passLeptons<< "(nl="<<leptons.size()<<")"<<endl;
      if(passLeptons && SKIMME) continue;
	  ioutVars["nl"] = leptons.size();
	  
      //select jets
      btvSF.addBTagDecisions(ev);
      if(!isData) btvSF.updateBTagDecisions(ev);      
      std::vector<Jet> allJets = selector.getGoodJets(ev,minJetPt,2.4,leptons,{});
      bool passJets(allJets.size()>=minJetMultiplicity);

	  if(debug) cout << "passJets = "<<passJets<< "(njets="<<allJets.size()<<")"<<endl;
      if(!passJets && SKIMME) continue;
	  foutVars["HT"] = 0;
	  jet_n = allJets.size();
	  std::vector<Jet>      bJets,lightJets;
      for(size_t ij=0; ij<allJets.size(); ij++) {
		  foutVars["HT"]+=allJets[ij].pt();
		  jet_pt[ij] = allJets[ij].pt();
          if(ev.j_btag[allJets[ij].getJetIndex()]>0) {
			  ioutVars["nbj"]++;
			  bJets.push_back(allJets[ij]);
		  }
		  else {
			  ioutVars["nlj"]++;
			  lightJets.push_back(allJets[ij]);
		  }	  
      }
	  bool passBJets(bJets.size()>=minBJetMultiplicity);
	  if(debug) cout << "passBJets = "<<passBJets<< "(nbjets="<<bJets.size()<<")"<<endl;
	  if(!passBJets && SKIMME) continue;
	  
      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
	  
	  // compute reconstructed tops (only for nJ>=4n nbJ>=2, nL>=1):
	  if(passJets && passBJets && passLeptons){
		neutrinoPzComputer.SetMET(met);
		neutrinoPzComputer.SetLepton(leptons[0].p4());
		float nupz=neutrinoPzComputer.Calculate();
		TLorentzVector neutrinoP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(met.Pt()*met.Pt()+nupz*nupz));
	  
		foutVars["chi2min"]=999; float invm_TOP = 1. / m_TOP; //, invm_W = 1. / m_W;
		for(int iblep=0;iblep<jet_n;iblep++){
		  if(ev.j_btag[allJets[iblep].getJetIndex()]==0) continue; // not a b jet
		  float _chi2 = TMath::Power( invm_TOP*(allJets[iblep].p4() + leptons[0].p4() + neutrinoP4).M() - 1, 2 );
		for(int ibhad=0;ibhad<jet_n;ibhad++){
		  if(ibhad==iblep) continue;
		  if(ev.j_btag[allJets[ibhad].getJetIndex()]==0) continue; // not a b jet
		for(int ij1=0;ij1<jet_n;ij1++){
		  if(ij1==ibhad || ij1==iblep) continue;
		for(int ij2=0;ij2<jet_n;ij2++){
		  if(ij2==ij1 || ij2==ibhad || ij2==iblep) continue;
		  _chi2+=  TMath::Power( invm_TOP*(allJets[ibhad].p4() + allJets[ij1].p4() + allJets[ij2].p4()).M() - 1, 2 );
		  //_chi2+=  TMath::Power( invm_W*(allJets[ij1].p4() + allJets[ij2].p4()).M() - 1, 2 );
		  if(_chi2 < foutVars["chi2min"]) {foutVars["chi2min"] = _chi2; ttbar_idx[0] = iblep; ttbar_idx[1] = ibhad; ttbar_idx[2] = ij1; ttbar_idx[3] = ij2;}
		}}}}
	  
		TLorentzVector tlep = (allJets.at(ttbar_idx[0]).p4() + leptons[0].p4() + neutrinoP4);
		TLorentzVector thad = (allJets.at(ttbar_idx[1]).p4() + allJets.at(ttbar_idx[2]).p4() + allJets.at(ttbar_idx[3]).p4());

		foutVars["pTtt"] = (tlep+thad).Pt();
		foutVars["mtt"] = (tlep+thad).M();
		foutVars["Ytt"] = (tlep+thad).Rapidity();
		foutVars["dPhitt"] = TMath::Abs(tlep.DeltaPhi(thad));
		if(debug) cout << "chi2min = "<<foutVars["chi2min"]<< "(ttbar_idx={"<<ttbar_idx[0]<<","<<ttbar_idx[1]<<","<<ttbar_idx[2]<<","<<ttbar_idx[3]<<"})"<<endl;
	  }
	  else{foutVars["pTtt"]=foutVars["mtt"]=-1; foutVars["Ytt"]=foutVars["dPhitt"]=-999;}
      //event weight
      float evWgt(1.0);
      
      //Data specific: check event rates after selection
      if(isData){
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
          ht.fill("ratevsrun",runBin,lumi,"inc");
        }else{
          //cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
		foutVars["gen_mtt"] = -1; foutVars["gen_Ytt"] = -999;
      }
	  else{ //MC specific: compute event weight and get truth info

        float normWgt(normH? normH->GetBinContent(1) : 1.0);        
        TString period = lumi.assignRunPeriod();
        float puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        EffCorrection_t selSF(1.0,0.0);// = lepEffH.getOfflineCorrection(leptons[0], period);
        EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});

        evWgt  = normWgt*puWgt*selSF.first*l1prefireProb.first;
        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);  
				
		// truth information (if available)
		TLorentzVector top1, top2, b1, b2, w1d1, w1d2, w2d1, w2d2;
	    if(debug) cout << "ngtop="<<ev.ngtop<<endl;
		for(int i=0;i<ev.ngtop;i++){
			if(ev.gtop_id[i]==6 && (i+7)<ev.ngtop)  {
				top1.SetPtEtaPhiM(ev.gtop_pt[i],ev.gtop_eta[i],ev.gtop_phi[i],ev.gtop_m[i]);
				b1.SetPtEtaPhiM(ev.gtop_pt[i+1],ev.gtop_eta[i+1],ev.gtop_phi[i+1],ev.gtop_m[i+1]);
				w1d1.SetPtEtaPhiM(ev.gtop_pt[i+6],ev.gtop_eta[i+6],ev.gtop_phi[i+6],ev.gtop_m[i+6]);
				w1d2.SetPtEtaPhiM(ev.gtop_pt[i+7],ev.gtop_eta[i+7],ev.gtop_phi[i+7],ev.gtop_m[i+7]);
				top_type[0]=100*fabs(ev.gtop_id[i+6])+fabs(ev.gtop_id[i+7]);
			}
			if(ev.gtop_id[i]==-6 && (i+6)<ev.ngtop) {
				top2.SetPtEtaPhiM(ev.gtop_pt[i],ev.gtop_eta[i],ev.gtop_phi[i],ev.gtop_m[i]);
				b2.SetPtEtaPhiM(ev.gtop_pt[i+1],ev.gtop_eta[i+1],ev.gtop_phi[i+1],ev.gtop_m[i+1]);
				w2d1.SetPtEtaPhiM(ev.gtop_pt[i+5],ev.gtop_eta[i+5],ev.gtop_phi[i+5],ev.gtop_m[i+5]);
				w2d2.SetPtEtaPhiM(ev.gtop_pt[i+6],ev.gtop_eta[i+6],ev.gtop_phi[i+6],ev.gtop_m[i+6]);
				top_type[1]=100*fabs(ev.gtop_id[i+5])+fabs(ev.gtop_id[i+6]);
			}
			if(ev.gtop_id[i]==2212){ // stable proton with xi < XIMAX (small momentum loss)
				if(ev.gtop_pz[i]>(1-XIMAX)*EBEAM) {foutVars["xip_truth"] = 1 - ev.gtop_pz[i]/EBEAM; ioutVars["nxip"]++;}
				else if(ev.gtop_pz[i]<(XIMAX-1)*EBEAM) {foutVars["xin_truth"] = 1 + ev.gtop_pz[i]/EBEAM; ioutVars["nxin"]++;}
			}
		}
		if(top1.Pt() && top2.Pt()){ // check if found pair of top quarks
			float mtt = (top1+top2).M();
			float Ytt = (top1+top2).Rapidity();
			foutVars["gen_mtt"] = mtt;
			foutVars["gen_Ytt"] = Ytt;

			float match_dR = 0.2;
			for(int ijet=0;ijet<jet_n;ijet++){
				if(b1.DeltaR(allJets[ijet].p4())<match_dR) {ttbar_assignment[0]=ijet; ttbar_assignment_res[0]=b1.Pt()/allJets[ijet].pt();
				tops_jet_pt[0]=b1.Pt();tops_jet_eta[0]=b1.Eta();tops_jet_phi[0]=b1.Phi();}
				if(w1d1.DeltaR(allJets[ijet].p4())<match_dR){ ttbar_assignment[1]=ijet; ttbar_assignment_res[1]=w1d1.Pt()/allJets[ijet].pt();
				tops_jet_pt[1]=w1d1.Pt();tops_jet_eta[1]=w1d1.Eta();tops_jet_phi[1]=w1d1.Phi();}
				if(w1d2.DeltaR(allJets[ijet].p4())<match_dR){ ttbar_assignment[2]=ijet; ttbar_assignment_res[2]=w1d2.Pt()/allJets[ijet].pt();
				tops_jet_pt[2]=w1d2.Pt();tops_jet_eta[2]=w1d2.Eta();tops_jet_phi[2]=w1d2.Phi();}
				if(b2.DeltaR(allJets[ijet].p4())<match_dR){ ttbar_assignment[3]=ijet; ttbar_assignment_res[3]=b2.Pt()/allJets[ijet].pt();
				tops_jet_pt[3]=b2.Pt();tops_jet_eta[3]=b2.Eta();tops_jet_phi[3]=b2.Phi();}
				if(w2d1.DeltaR(allJets[ijet].p4())<match_dR){ ttbar_assignment[4]=ijet; ttbar_assignment_res[4]=w2d1.Pt()/allJets[ijet].pt();
				tops_jet_pt[4]=w2d1.Pt();tops_jet_eta[4]=w2d1.Eta();tops_jet_phi[4]=w2d1.Phi();}
				if(w2d2.DeltaR(allJets[ijet].p4())<match_dR){ ttbar_assignment[5]=ijet; ttbar_assignment_res[5]=w2d2.Pt()/allJets[ijet].pt();
				tops_jet_pt[5]=w2d2.Pt();tops_jet_eta[5]=w2d2.Eta();tops_jet_phi[5]=w2d2.Phi();}
			}
			if(top_type[0]>1000 && passLeptons){
				if(w1d2.DeltaR(leptons[0].p4())<match_dR) {ttbar_assignment[1]=0;ttbar_assignment[2]=0; 
				ttbar_assignment_res[1]=w1d1.Pt()/leptons[0].pt();ttbar_assignment_res[2]=-1;}
				else {ttbar_assignment[1]=ttbar_assignment[2]=-1; ttbar_assignment_res[1]=ttbar_assignment_res[2]=0;}
			}
			if(top_type[1]>1000 && passLeptons){
				if(w2d1.DeltaR(leptons[0].p4())<match_dR) {ttbar_assignment[4]=0; ttbar_assignment[5]=0;
				ttbar_assignment_res[4]=w2d2.Pt()/leptons[0].pt();ttbar_assignment_res[5]=-1;}
				else {ttbar_assignment[4]=ttbar_assignment[5]=-1; ttbar_assignment_res[4]=ttbar_assignment_res[5]=0;}
			}		
		} // end if ttbar
	  } // end else isData (MC specific)
      

      //Fill data with roman pot information
      if (debug) cout << "Start PPS reco " << endl;
      for (int ift=0; ift<ev.nfwdtrk; ift++) { // nppstrk nfwdtrk

        const unsigned short pot_raw_id = ev.fwdtrk_pot[ift]; // ppstrk_pot fwdtrk_pot
		float xi = ev.fwdtrk_xi[ift]; 
		if(xi<XIMIN) continue; // cut on the event selection
//		float x  = ev.ppstrk_x[ift];

        //single pot reconstruction
        if(ev.fwdtrk_method[ift]==0){

		if(pot_raw_id==123){
			pixel_neg_xi[pixel_neg_n]=xi;
			pixel_neg_n++;}
		else if(pot_raw_id==23){
			pixel_pos_xi[pixel_pos_n]=xi;
			pixel_pos_n++;}
		else if(pot_raw_id==103){strip_neg_xi[strip_neg_n]=xi; strip_neg_n++;}
		else if(pot_raw_id==3){strip_pos_xi[strip_pos_n]=xi; strip_pos_n++;}
		else continue;
		}

        //multi reconstruction
        if(ev.fwdtrk_method[ift]==1){
		multi_n++;
		if (pot_raw_id<100){ multi_pos_xi=xi; multi_pos_xi_eff = ev.fwdtrk_xiSF[ift];}
		else { multi_neg_xi=xi;  multi_neg_xi_eff = ev.fwdtrk_xiSF[ift];}
		}
				
      } // loop over RP tracks
	  pixel_n=pixel_pos_n+pixel_neg_n;
	  strip_n=strip_pos_n+strip_neg_n;
	  
	  
	  // mass distribution for all:
	  mass_pixel_n=0;
	  for(int i=0;i<pixel_pos_n;i++){
		  for(int j=0;j<pixel_neg_n;j++){
			  mass_pixel_all[mass_pixel_n++]=13000.*sqrt(pixel_pos_xi[i]*pixel_neg_xi[j]);
		  }
	  }
	  
	  // RP variables:
	  if(multi_n==2){
	    foutVars["xip"] = multi_pos_xi;
		foutVars["xin"] = multi_neg_xi;
		foutVars["mpp"] = 13000.*sqrt(multi_pos_xi*multi_neg_xi);
		foutVars["Ypp"] = 0.5*TMath::Log(multi_pos_xi/multi_neg_xi);
	    if(debug) cout << "multi_n="<<multi_n<<", mpp = " << foutVars["mpp"] << endl;
	  }	  

	  // Save output tree
	  if(skimtree) outT->Fill();
 
    } // END EVENT LOOP
  
  //close input file
  f->Close();
  
  //save histos to file
  if(skimtree) cout << endl << "Writes " << fOut->GetName() << " with " << outT->GetEntries() <<  " events." << endl;  
  else cout << endl << "Writes " << fOut->GetName() << "." << endl;  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }  
  if(skimtree){outT->Write();}
  
  fOut->Close();
}
