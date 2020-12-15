#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/FwdProtons.h"
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

void RunFwdProtons(const TString in_fname,
                      TString outname,
                      TH1F *normH, 
                      TH1F *genPU, 
                      TString era,
                      bool debug) 
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  float EBEAM(6500);
  float XIMIN(0.01);
  float XIMAX(0.2);
  bool SKIMME(true);
    
  //preselection cuts to apply
  float minLeptonPt(30);
  float minJetPt(15);
  size_t minJetMultiplicity(4);
  size_t minBJetMultiplicity(2);

  //const char* CMSSW_BASE = getenv("CMSSW_BASE");
  MiniEvent_t ev;  

  //CORRECTIONS: LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  std::map<Int_t,Float_t> lumiPerRun=lumi.lumiPerRun();
  
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

  // Fixel efficiency wrapper
  PPSEfficiencies PPSeff;
  if(isData) PPSeff.init();
  else PPSeff.init("2017C1");


  //PREPARE OUTPUT (BOOK SOME HISTOGRAMS)
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  
  TTree *outT=new TTree("tree","tree");
	  
  // Event info
  outT->Branch("run",&ev.run,"run/i");
  outT->Branch("event",&ev.event,"event/l");
  outT->Branch("lumi",&ev.lumi,"lumi/i");
  outT->Branch("nvtx",&ev.nvtx,"nvtx/I");
  //outT->Branch("nfill",&ev.nfill,"nfill/I");
  outT->Branch("beamXangle",&ev.beamXangle,"beamXangle/F");
  outT->Branch("betastar",&ev.betaStar,"betastar/F");
  
  //Variables (bools, ints, floats):
  TString bvars[]={"hasMTrigger","hasETrigger","hasEMTrigger","hasMETrigger","hasMMTrigger","hasEETrigger","hasSLT","hasDLT"};
  std::map<TString,bool> boutVars;
  for(size_t i=0; i<sizeof(bvars)/sizeof(TString); i++) boutVars[bvars[i]]=false;
  
  TString ivars[]={"nl","nj","nbj","nlj"};
  std::map<TString,Int_t> ioutVars;
  for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++) ioutVars[ivars[i]]=0;
  
  TString fvars[]={"HT"};
  std::map<TString,Float_t> foutVars;
  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++) foutVars[fvars[i]]=0.;

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
	
  //Protons
  int pixel_pos_n(0), pixel_neg_n(0), pixel_n(0);
  int strip_pos_n(0), strip_neg_n(0), strip_n(0);
  int proton_pos_n(0), proton_neg_n(0), proton_n(0);
  int multi_n(0);
  float pixel_pos_xi[ev.MAXPROTONS], pixel_neg_xi[ev.MAXPROTONS];
  float pixel_pos_eff[ev.MAXPROTONS], pixel_neg_eff[ev.MAXPROTONS];
  float pixel_pos_x[ev.MAXPROTONS], pixel_neg_x[ev.MAXPROTONS];
  float pixel_pos_y[ev.MAXPROTONS], pixel_neg_y[ev.MAXPROTONS];
  float strip_pos_xi[ev.MAXPROTONS], strip_neg_xi[ev.MAXPROTONS]; 
  float strip_pos_eff[ev.MAXPROTONS], strip_neg_eff[ev.MAXPROTONS]; 
  float strip_pos_x[ev.MAXPROTONS], strip_neg_x[ev.MAXPROTONS]; 
  float strip_pos_y[ev.MAXPROTONS], strip_neg_y[ev.MAXPROTONS]; 
  float proton_pos_xi[ev.MAXPROTONS], proton_neg_xi[ev.MAXPROTONS]; 
  float proton_pos_t[ev.MAXPROTONS], proton_neg_t[ev.MAXPROTONS]; 
  float multi_pos_xi, multi_neg_xi;
  
  outT->Branch("multi_n",&multi_n,"multi_n/I");
  outT->Branch("multi_pos_xi",&multi_pos_xi,"multi_pos_xi/F");
  outT->Branch("multi_neg_xi",&multi_neg_xi,"multi_neg_xi/F");
  outT->Branch("proton_pos_n",&proton_pos_n,"proton_pos_n/I");
  outT->Branch("proton_neg_n",&proton_neg_n,"proton_neg_n/I");
  outT->Branch("proton_n",&proton_n,"proton_n/i");
  outT->Branch("pixel_pos_n",&pixel_pos_n,"pixel_pos_n/I");
  outT->Branch("pixel_neg_n",&pixel_neg_n,"pixel_neg_n/I");
  outT->Branch("pixel_n",&pixel_n,"pixel_n/i");
  outT->Branch("strip_pos_n",&strip_pos_n,"strip_pos_n/I");
  outT->Branch("strip_neg_n",&strip_neg_n,"strip_neg_n/I");
  outT->Branch("strip_n",&strip_n,"strip_n/I");	
  outT->Branch("proton_pos_xi",proton_pos_xi,"proton_pos_xi[proton_pos_n]/F");
  outT->Branch("proton_neg_xi",proton_neg_xi,"proton_neg_xi[proton_neg_n]/F");
  outT->Branch("proton_pos_t",proton_pos_t,"proton_pos_t[proton_pos_n]/F");
  outT->Branch("proton_neg_t",proton_neg_t,"proton_neg_t[proton_neg_n]/F");
  outT->Branch("pixel_pos_xi",pixel_pos_xi,"pixel_pos_xi[pixel_pos_n]/F");
  outT->Branch("pixel_neg_xi",pixel_neg_xi,"pixel_neg_xi[pixel_neg_n]/F");
  outT->Branch("pixel_pos_eff",pixel_pos_eff,"pixel_pos_eff[pixel_pos_n]/F");
  outT->Branch("pixel_neg_eff",pixel_neg_eff,"pixel_neg_eff[pixel_neg_n]/F");
  outT->Branch("strip_pos_xi",strip_pos_xi,"strip_pos_xi[strip_pos_n]/F");
  outT->Branch("strip_neg_xi",strip_neg_xi,"strip_neg_xi[strip_neg_n]/F");
  outT->Branch("strip_pos_eff",strip_pos_eff,"strip_pos_eff[strip_pos_n]/F");
  outT->Branch("strip_neg_eff",strip_neg_eff,"strip_neg_eff[strip_neg_n]/F");
  outT->Branch("pixel_pos_x",pixel_pos_x,"pixel_pos_x[pixel_pos_n]/F");
  outT->Branch("pixel_pos_y",pixel_pos_y,"pixel_pos_y[pixel_pos_n]/F");
  outT->Branch("pixel_neg_x",pixel_neg_x,"pixel_neg_x[pixel_neg_n]/F");
  outT->Branch("pixel_neg_y",pixel_neg_y,"pixel_neg_y[pixel_neg_n]/F");
  outT->Branch("strip_pos_x",strip_pos_x,"strip_pos_x[strip_pos_n]/F");
  outT->Branch("strip_pos_y",strip_pos_y,"strip_pos_y[strip_pos_n]/F");
  outT->Branch("strip_neg_x",strip_neg_x,"strip_neg_x[strip_neg_n]/F");
  outT->Branch("strip_neg_y",strip_neg_y,"strip_neg_y[strip_neg_n]/F");
	  
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
  ht.addHist("norm",     new TH1F("h_norm",    ";Category; Events",4,0,4) );
  ht.getPlots()["norm"]->GetXaxis()->SetBinLabel(1,"SumUnWeighted");
  ht.getPlots()["norm"]->GetXaxis()->SetBinLabel(2,"Sumweighted");
  ht.getPlots()["norm"]->GetXaxis()->SetBinLabel(3,"accepted");
  ht.getPlots()["norm"]->GetXaxis()->SetBinLabel(4,"nentries");
  ht.getPlots()["norm"]->SetBinContent(1,counter->GetBinContent(1));
  ht.getPlots()["norm"]->SetBinContent(2,counter->GetBinContent(2));
  ht.getPlots()["norm"]->SetBinContent(3,counter->GetBinContent(3));
  ht.getPlots()["norm"]->SetBinContent(nentries,counter->GetBinContent(3));
  
  std::cout << "initialization done" << std::endl;

  //EVENT SELECTION WRAPPER (GETS LISTS OF PHYSICS OBJECTS FROM THE INPUT)
  SelectionTool selector(in_fname, false, triggerList);
  
  //EVENT LOOP
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      //if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)(iev+1)/(float)nentries); fflush(stdout); }
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
		
		
	  // Reset vars
      // ------------------------------------------------------------------------------------------------------- //
	  if(debug) cout << "Event #"<<iev<<" Reset vars" << endl;
	  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++) foutVars[fvars[i]]=0.;
	  for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++) ioutVars[ivars[i]]=0;
	  proton_pos_n=0;proton_neg_n=0;
	  pixel_pos_n=0;pixel_neg_n=0;
	  strip_pos_n=0;strip_neg_n=0;
	  multi_n=0; multi_pos_xi=0; multi_neg_xi=0;
      for(int irp=0; irp<ev.MAXPROTONS; irp++) {
		proton_pos_xi[irp]=0;proton_neg_xi[irp]=0;
		pixel_pos_xi[irp]=0;pixel_neg_xi[irp]=0;
		pixel_pos_x[irp]=0;pixel_neg_y[irp]=0;
		strip_pos_xi[irp]=0;strip_neg_xi[irp]=0;
		strip_pos_x[irp]=0;strip_neg_y[irp]=0;
	  }
      // ------------------------------------------------------------------------------------------------------- //

      //trigger
      boutVars["hasMTrigger"] = (selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits) ); 
	  boutVars["hasETrigger"] = (selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v",                                  ev.triggerBits) || 
                                selector.hasTriggerBit("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v",                     ev.triggerBits) ||
							    selector.hasTriggerBit("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v", ev.triggerBits));
      boutVars["hasMMTrigger"]= (selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                  ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",          ev.triggerBits));
	  boutVars["hasEETrigger"]= (selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev.triggerBits) );
	  boutVars["hasMETrigger"]= (selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits));
	  boutVars["hasEMTrigger"]= (selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits));
      
	  boutVars["hasSLT"] = selector.passSingleLeptonTrigger(ev); // boutVars["hasMTrigger"] || boutVars["hasETrigger"];
	  boutVars["hasDLT"] = (boutVars["hasMMTrigger"]||boutVars["hasEETrigger"]||boutVars["hasEMTrigger"]||boutVars["hasMETrigger"]);
      
	  if(debug) cout << "runNumber #"<<ev.run<<", event = " << ev.event<<", pass trigger (hasSLT) = "  << boutVars["hasSLT"] <<endl;
	  if(!boutVars["hasSLT"] && SKIMME) continue;

      //lepton selection
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);     
      SelectionTool::QualityFlags muId(SelectionTool::TIGHT); //TIGHT
      leptons = selector.selLeptons(leptons,muId,SelectionTool::MVA90,minLeptonPt,2.1);
      bool passLeptons(leptons.size()==1);
	  if(debug) cout << "passLeptons = "<<passLeptons<< "(nl="<<leptons.size()<<")"<<endl;
      if(!passLeptons && SKIMME) continue;
	  ioutVars["nl"] = leptons.size();
	  
	  
	  // ------ update latter with UL recommendations, id needed ------ 
      //btvSF.addBTagDecisions(ev);
      //if(!isData) btvSF.updateBTagDecisions(ev); 
	  
      //select jets 
      std::vector<Jet> allJets = selector.getGoodJets(ev,minJetPt,2.4,leptons,{});
	  foutVars["nj"] = allJets.size();
      bool passJets(foutVars["nj"]>=minJetMultiplicity);

	  if(debug) cout << "passJets = "<<passJets<< "(njets="<<allJets.size()<<")"<<endl;
      if(!passJets && SKIMME) continue;
	  std::vector<Jet>      bJets,lightJets;
      for(size_t ij=0; ij<allJets.size(); ij++) {
		  foutVars["HT"]+=allJets[ij].pt();
		  if(ev.j_btag[allJets[ij].getJetIndex()]) {
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
      }
	  else{ //MC specific: compute event weight and get truth info

        //float normWgt(normH? normH->GetBinContent(1) : 1.0);        
        TString period = lumi.assignRunPeriod();
        //float puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        //EffCorrection_t selSF(1.0,0.0); // = lepEffH.getOfflineCorrection(leptons[0], period);
		//EffCorrection_t trigSF(1.0,0.0);// = lepEffH.getTriggerCorrection(leptons,{},{},period);
        //EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});

        //evWgt  = normWgt*puWgt*trigSF.first*selSF.first*l1prefireProb.first;
        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);  
				
		// truth protons (if available)
	    if(debug) cout << "ngtop="<<ev.ngtop<<endl;
		for(int i=0;i<ev.ngtop;i++){
			if(ev.gtop_id[i]==2212){ // stable proton with xi < XIMAX (small momentum loss)
				if(ev.gtop_pz[i]>(1-XIMAX)*EBEAM) {proton_pos_xi[proton_pos_n++] = 1 - ev.gtop_pz[i]/EBEAM;}
				else if(ev.gtop_pz[i]<(XIMAX-1)*EBEAM) {proton_neg_xi[proton_neg_n] = 1 + ev.gtop_pz[i]/EBEAM; proton_neg_n++;}
			}
		}
		proton_n = proton_pos_n + proton_neg_n;
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
		if (pot_raw_id<100){ multi_pos_xi=xi;}
		else { multi_neg_xi=xi;}
		}
				
      } // loop over RP tracks
	  pixel_n=pixel_pos_n+pixel_neg_n;
	  strip_n=strip_pos_n+strip_neg_n;
	  
	  // add X,Y information:
	  int _pixel_neg_n=0, _pixel_pos_n=0, _strip_neg_n=0, _strip_pos_n=0;
	  for (int ift=0; ift<ev.nppstrk; ift++) {
		const unsigned short pot_raw_id = ev.ppstrk_pot[ift]; // ppstrk_pot fwdtrk_pot
	    if(pot_raw_id==123 && _pixel_neg_n<pixel_neg_n){
			pixel_neg_x[_pixel_neg_n]=ev.ppstrk_x[ift];
			pixel_neg_y[_pixel_neg_n]=ev.ppstrk_y[ift];
			pixel_neg_eff[_pixel_neg_n]=PPSeff.GetPixelEff(ev.ppstrk_x[ift],ev.ppstrk_y[ift],0);
			_pixel_neg_n++;}
		else if(pot_raw_id==23 && _pixel_pos_n<pixel_pos_n){
			pixel_pos_x[_pixel_pos_n]=ev.ppstrk_x[ift];
			pixel_pos_y[_pixel_pos_n]=ev.ppstrk_y[ift];
			pixel_pos_eff[_pixel_pos_n]=PPSeff.GetPixelEff(ev.ppstrk_x[ift],ev.ppstrk_y[ift],1);
			_pixel_pos_n++;}
		else if(pot_raw_id==103 && _strip_neg_n<strip_neg_n){
			strip_neg_x[_strip_neg_n]=ev.ppstrk_x[ift];
			strip_neg_y[_strip_neg_n]=ev.ppstrk_y[ift];
			strip_neg_eff[_strip_neg_n]=PPSeff.GetStripEff(ev.ppstrk_x[ift],ev.ppstrk_y[ift],0);
			_strip_neg_n++;}
		else if(pot_raw_id==3 && _strip_pos_n<strip_pos_n){
			strip_pos_x[_strip_pos_n]=ev.ppstrk_x[ift];
			strip_pos_y[_strip_pos_n]=ev.ppstrk_y[ift];
			strip_pos_eff[_strip_pos_n]=PPSeff.GetStripEff(ev.ppstrk_x[ift],ev.ppstrk_y[ift],1);
			_strip_pos_n++;}
		else continue;
	  }
	  
	  

	  //if(multi_n!=2 && SKIMME) continue;
	  //if(foutVars["chi2min"]>XXX && SKIMME) continue;

	  // Save output tree
	  outT->Fill();
 
    } // END EVENT LOOP
  
  //close input file
  f->Close();
  
  //save histos to file
  cout << endl << "Writes " << fOut->GetName() << " with " << outT->GetEntries() <<  " events." << endl;  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }  
  outT->Write();
  
  fOut->Close();
}

// PPSEfficiencies class:

float PPSEfficiencies::GetPixelEff(float x, float y, int pos){
	if(isData) return 0;
	if(pos>0) {
		bin = h_pixels45->FindBin(x,y);
		return h_pixels45->GetBinContent(bin);
	}
	bin = h_pixels56->FindBin(x,y);
	return h_pixels56->GetBinContent(bin);
}

float PPSEfficiencies::GetStripEff(float x, float y, int pos){
	if(isData) return 0;
	if(pos>0) {
		bin = h_strips45->FindBin(x,y);
		return h_strips45->GetBinContent(bin);
	}
	bin = h_strips56->FindBin(x,y);
	return h_strips56->GetBinContent(bin);
}

void PPSEfficiencies::init(TString era){
			
	f_strips = new TFile(path+"/Strips/StripsTracking/PreliminaryEfficiencies_October92019_1D2DMultiTrack.root");
	f_pixels = new TFile(path+"/Pixel/RPixTracking/pixelEfficiencies_radiation.root");

	h_strips45  = (TH2D *)f_strips->Get("Strips/2017/"+era+"/h45_"+era+"_all_2D");
	h_strips56  = (TH2D *)f_strips->Get("Strips/2017/"+era+"/h56_"+era+"_all_2D");
	h_pixels45  = (TH2D *)f_pixels->Get("Pixel/2017/"+era+"/h45_220_"+era+"_all_2D");
	h_pixels56  = (TH2D *)f_pixels->Get("Pixel/2017/"+era+"/h56_220_"+era+"_all_2D");
	if(!h_strips45 || !h_strips56) cout << "ERROR: didn't loaded strip efficiency histograms! " << endl;
	if(!h_pixels45 || !h_pixels56) cout << "ERROR: didn't loaded pixel efficiency histograms! " << endl;
	
}

void PPSEfficiencies::init(){
	isData=true;
}


