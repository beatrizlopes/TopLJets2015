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

#include "TopLJets2015/TopAnalysis/src/GetPixEff.cc"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include "TMath.h"

using namespace std;

#define ADDVAR(x,name,t,tree) tree->Branch(name,x,TString(name)+TString(t))

double m_TOP = 173.1;
double m_W   =  80.379;
double m_NU  =  0.;

//
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
  //const char* CMSSW_BASE = getenv("CMSSW_BASE");
  MiniEvent_t ev;  

  //preselection cuts to apply
  float minLeptonPt(30);
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
  
  //Correctoin PPS
  PixEff PPS_eff;
  
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
  //t->GetEntry(0);
  std::cout << "...producing " << outname << " from " << nentries << " events" << std::endl;  
    

  //PREPARE OUTPUT (BOOK SOME HISTOGRAMS)
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  
  //Variables (bools, ints, floats):
  TString bvars[]={"hasMTrigger","hasETrigger","hasEMTrigger","hasMETrigger","hasMMTrigger","hasEETrigger","hasSLT","hasDLT"};
  std::map<TString,bool> boutVars;
  for(size_t i=0; i<sizeof(bvars)/sizeof(TString); i++) boutVars[bvars[i]]=false;
  
  TString ivars[]={"nl","nj","nbj","nlj","nxip","nxin"};
  std::map<TString,Int_t> ioutVars;
  for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++) ioutVars[ivars[i]]=0;
  
  TString fvars[]={"gen_mtt","gen_Ytt","mpp","Ypp","xip","xin","xip_truth","xin_truth",
                    "HT", "chi2min","mtt","Ytt","pTtt","dPhitt"};
  std::map<TString,Float_t> foutVars;
  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++) foutVars[fvars[i]]=0.;
  
  int pixel_pos_n(0), pixel_neg_n(0), pixel_n(0);
  int strip_pos_n(0), strip_neg_n(0), strip_n(0);
  int multi_n(0), mass_pixel_n(0);
  float pixel_pos_xi[50], pixel_neg_xi[50];
  float pixel_pos_xi_eff[50], pixel_neg_xi_eff[50];
  float strip_pos_xi[50], strip_neg_xi[50]; 
  float mass_pixel_all[200];
  float multi_pos_xi,multi_neg_xi;
	
  TTree *outT=new TTree("tree","tree");
  if(skimtree) {
	  
	outT->Branch("run",&ev.run,"run/i");
	outT->Branch("event",&ev.event,"event/l");
	outT->Branch("lumi",&ev.lumi,"lumi/i");
	outT->Branch("nvtx",&ev.nvtx,"nvtx/I");
    if(t->FindBranch("nchPV")) outT->Branch("nchPV",&ev.nchPV,"nchPV/I");
	ADDVAR(&ev.beamXangle,"beamXangle","/F",outT);
	
	
	ADDVAR(&ev.rho,"rho","/F",outT);
	if(t->FindBranch("met_pt")){
	ADDVAR(&ev.met_pt,"met_pt","/F",outT);
	ADDVAR(&ev.met_phi,"met_phi","/F",outT);
	ADDVAR(&ev.met_sig,"met_sig","/F",outT);
	}
	
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
    outT->Branch("pixel_pos_n",&pixel_pos_n,"pixel_pos_n/I");
    outT->Branch("pixel_neg_n",&pixel_neg_n,"pixel_neg_n/I");
    outT->Branch("pixel_n",&pixel_n,"pixel_n/i");
    outT->Branch("strip_pos_n",&strip_pos_n,"strip_pos_n/I");
    outT->Branch("strip_neg_n",&strip_neg_n,"strip_neg_n/I");
    outT->Branch("strip_n",&strip_n,"strip_n/I");	
	outT->Branch("pixel_pos_xi",pixel_pos_xi,"pixel_pos_xi[pixel_pos_n]/F");
	outT->Branch("pixel_pos_xi_eff",pixel_pos_xi_eff,"pixel_pos_xi_eff[pixel_pos_n]/F");
	outT->Branch("pixel_neg_xi",pixel_neg_xi,"pixel_neg_xi[pixel_neg_n]/F");
	outT->Branch("pixel_neg_xi_eff",pixel_neg_xi_eff,"pixel_neg_xi_eff[pixel_neg_n]/F");
	outT->Branch("strip_pos_xi",strip_pos_xi,"strip_pos_xi[strip_pos_n]/F");
	outT->Branch("strip_neg_xi",strip_neg_xi,"strip_neg_xi[strip_neg_n]/F");
    outT->Branch("mass_pixel_n",&mass_pixel_n,"mass_pixel_n/I");
	outT->Branch("mass_pixel_all",mass_pixel_all,"mass_pixel_all[mass_pixel_n]/F");
	
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
	  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++) foutVars[fvars[i]]=0.;
	  for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++) ioutVars[ivars[i]]=0;
	  
	  int runNumber = (ev.isData) ? ev.run : 297050; // era2017B

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
      
	  if(!boutVars["hasSLT"]) continue;

      //lepton selection
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);     
      SelectionTool::QualityFlags muId(SelectionTool::LOOSE); //TIGHT
      leptons = selector.selLeptons(leptons,muId,SelectionTool::MVA90,minLeptonPt,2.1);
      if(leptons.size()==0) continue;
	  ioutVars["nl"] = leptons.size();
      //if(leptons[0].id()!=13) continue;
	  
      //select jets
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);      
      std::vector<Jet> allJets = selector.getGoodJets(ev,30.,2.4,leptons,{});
      bool passJets(allJets.size()>=minJetMultiplicity);

      //cout << " passJets = " << passJets << endl;
      if(!passJets) continue;
	  foutVars["HT"] = 0;
	  ioutVars["nj"] = allJets.size();
	  std::vector<Jet>      bJets,lightJets;
      for(size_t ij=0; ij<allJets.size(); ij++) {
		  foutVars["HT"]+=allJets[ij].pt();
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
	  if(!passBJets) continue;
	  
      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
	  
	  // compute reconstructed tops:
	  neutrinoPzComputer.SetMET(met);
	  neutrinoPzComputer.SetLepton(leptons[0].p4());
	  float nupz=neutrinoPzComputer.Calculate();
	  TLorentzVector neutrinoP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(met.Pt()*met.Pt()+nupz*nupz));
	  
	  foutVars["chi2min"]=999; float invm_TOP = 1. / m_TOP, invm_W = 1. / m_W; int ttbar_idx[4] = {-1,-1,-1,-1};
	  for(int iblep=0;iblep<ioutVars["nj"];iblep++){
		  if(ev.j_btag[allJets[iblep].getJetIndex()]==0) continue; // not a b jet
		  float _chi2 = TMath::Power( invm_TOP*(allJets[iblep].p4() + leptons[0].p4() + neutrinoP4).M() - 1, 2 );
	  for(int ibhad=0;ibhad<ioutVars["nj"];ibhad++){
		  if(ibhad==iblep) continue;
		  if(ev.j_btag[allJets[ibhad].getJetIndex()]==0) continue; // not a b jet
	  for(int ij1=0;ij1<ioutVars["nj"];ij1++){
		  if(ij1==ibhad || ij1==iblep) continue;
	  for(int ij2=0;ij2<ioutVars["nj"];ij2++){
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

      //event weight
      float evWgt(1.0);
      
      //data specific: check event rates after selection
      if(ev.isData){
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
          ht.fill("ratevsrun",runBin,lumi,"inc");
        }else{
          //cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
      }

      //MC specific: compute event weight and get truth info
      if (!ev.isData) {

        float normWgt(normH? normH->GetBinContent(1) : 1.0);        
        TString period = lumi.assignRunPeriod();
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        EffCorrection_t selSF(1.0,0.0);// = lepEffH.getOfflineCorrection(leptons[0], period);
        EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});

        evWgt  = normWgt*puWgt*selSF.first*l1prefireProb.first;
        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);  
		
		// proton reconstruction information:
		// missing
		
		// truth information (if available)
		TLorentzVector top1, top2;
		for(int i=0;i<ev.ngtop;i++){
			if(ev.gtop_id[i]==6)  top1.SetPtEtaPhiM(ev.gtop_pt[i],ev.gtop_eta[i],ev.gtop_phi[i],ev.gtop_m[i]);
			if(ev.gtop_id[i]==-6) top2.SetPtEtaPhiM(ev.gtop_pt[i],ev.gtop_eta[i],ev.gtop_phi[i],ev.gtop_m[i]);
			if(ev.gtop_id[i]==2212){ // stable proton
				if(ev.gtop_pz[i]>(1-XIMAX)*EBEAM) {foutVars["xip_truth"] = 1 - ev.gtop_pz[i]/EBEAM; ioutVars["nxip"]++;}
				else if(ev.gtop_pz[i]<(XIMAX-1)*EBEAM) {foutVars["xin_truth"] = 1 + ev.gtop_pz[i]/EBEAM; ioutVars["nxin"]++;}
			}
		}
		if(top1.Pt() && top2.Pt()){ // check if found pair of top quarks
			float mtt = (top1+top2).M();
			float Ytt = (top1+top2).Rapidity();
			foutVars["gen_mtt"] = mtt;
			foutVars["gen_Ytt"] = Ytt;
		}
				
      }
	  else{
		foutVars["gen_mtt"] = -1; foutVars["gen_Ytt"] = -999;
	  }
      

      //ht.fill("nvtx",       ev.nvtx,        evWgt, "inc");
      //ht.fill("njets",      allJets.size(), evWgt, "inc");


      //lepton-b systems
      for(size_t ij=0; ij<allJets.size(); ij++) 
        {
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);
          if(!passBtag) continue;

          float mlb( (leptons[0]+allJets[ij]).M() );
          //std::vector<TString> tags={"inc",leptons[0].charge()>0 ? "plus" : "minus"};
          //ht.fill("mlb",mlb,evWgt,tags);
        }

      //fill data with roman pot information
	  pixel_pos_n=0;pixel_neg_n=0;
	  strip_pos_n=0;strip_neg_n=0;
	  multi_n=0; multi_pos_xi=0; multi_neg_xi=0;
      for(int irp=0; irp<50; irp++) {
		pixel_pos_xi[irp]=0;pixel_neg_xi[irp]=0;
		pixel_pos_xi_eff[irp]=0;pixel_neg_xi_eff[irp]=0;
		strip_pos_xi[irp]=0;strip_neg_xi[irp]=0;
	  }
	  
	//  if (ev.isData) {
	//	  const edm::EventID ev_id( ev.run, ev.lumi, ev.event );
	//	  const ctpps::conditions_t lhc_cond = lhc_conds.get( ev_id );
	//	  ioutVars["beamXangle"] = std::round(lhc_cond.crossing_angle);
	//  }
	  

	  
      for (int ift=0; ift<ev.nfwdtrk; ift++) { // nppstrk nfwdtrk

        const unsigned short pot_raw_id = ev.fwdtrk_pot[ift]; // ppstrk_pot fwdtrk_pot
		float xi = ev.fwdtrk_xi[ift];
		if(xi<XIMIN) continue; // cut on the event selection
//		float x  = ev.ppstrk_x[ift];
		if(ev.fwdtrk_vx[ift]==0) continue; // FIXME used for debug!!!

        //single pot reconstruction
        if(ev.fwdtrk_method[ift]==0){

		if(pot_raw_id==23){
			pixel_neg_xi[pixel_neg_n]=xi;
			pixel_neg_xi_eff[pixel_neg_n] = (ev.isData) ? 1 : PPS_eff.getEff(xi,pot_raw_id,runNumber);
			pixel_neg_n++;}
		else if(pot_raw_id==123){
			pixel_pos_xi[pixel_pos_n]=xi;
			pixel_pos_xi_eff[pixel_pos_n] = (ev.isData) ? 1 : PPS_eff.getEff(xi,pot_raw_id,runNumber);
			pixel_pos_n++;}
		else if(pot_raw_id==3){strip_neg_xi[strip_neg_n]=xi; strip_neg_n++;}
		else if(pot_raw_id==103){strip_pos_xi[strip_pos_n]=xi; strip_pos_n++;}
		else continue;
		
		}

        //multi reconstruction
        if(ev.fwdtrk_method[ift]==1){
		multi_n++;
		if (pot_raw_id>100) multi_pos_xi=xi;
		else multi_neg_xi=xi;
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
	  }	  

	  // Save output tree
	  if(skimtree) outT->Fill();
 
    } // END EVENT LOOP
  
  //close input file
  f->Close();
  
  //save histos to file
  cout << endl << "Writes " << fOut->GetName() << endl;  
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
