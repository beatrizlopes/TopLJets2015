#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveDiTau.h"
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

void RunExclusiveDiTau(const TString in_fname,
                      TString outname,
					  TH1F *genPU,
                      TString era,
                      bool debug,
					  std::string systVar) 
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  float EBEAM(6500);
  float XIMIN(0.03);
  float XIMAX(0.2);
  bool SKIMME(true);
  bool isData = in_fname.Contains("Data13TeV");
    
  //const char* CMSSW_BASE = getenv("CMSSW_BASE");
  MiniEvent_t ev;  

  //preselection cuts to apply
  float minLeptonPt(30);
  float minJetPt(25);  
  
  //CORRECTIONS: LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  
  //CORRECTIONS: LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper lepEffH(isData,era);

  //CORRECTIONS: L1-prefire 
  L1PrefireEfficiencyWrapper l1PrefireWR(in_fname.Contains("Data13TeV"),era);
  
  
  //CORRECTIONS: B-TAG CALIBRATION
  BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
  
  //JEC/JER
  JECTools jec(era);
	
  //READ TREE FROM FILE
  TFile *f = TFile::Open(in_fname);
  if(f==NULL || f->IsZombie()) {
    cout << "Corrupted or missing file " << in_fname << endl;
    return;
  }
  TH1 *counter=(TH1 *)f->Get("analysis/counter");
  if(!counter) {cout << "Corrupted or missing counter: \"analysis/counter\" " << endl;return;}
  TH1 *RPcount=(TH1 *)f->Get("analysis/RPcount");
  if(!RPcount) {cout << "Corrupted or missing RPcount: \"analysis/RPcount\" " << endl;return;}
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  if(!triggerList) {cout << "Corrupted or missing triggerList: \"analysis/triggerList\" " << endl;return;}
  TTree *t = (TTree*)f->Get("analysis/tree");
  if(!t) {cout << "Corrupted or missing tree: \"analysis/tree\" " << endl;return;}
  attachToMiniEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = min(1000,nentries); //restrict number of entries for testing
  std::cout << "...producing " << outname << " from " << nentries << " events" << std::endl;  

  //PREPARE OUTPUT (BOOK SOME HISTOGRAMS)
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  if(isData) RPcount->Write();
  
  // Proton variables:
  const int MAXPR = 10;
  Int_t nRecoProtCand;
  Float_t pr_xi[MAXPR], pr_t[MAXPR];//, pr_ThX[MAXPR], pr_ThY[MAXPR]; 
  Short_t pr_rpid[MAXPR], pr_ismultirp[MAXPR]; //,  pr_arm[MAXPR];
  
  // BOOK PROTON TREE (DATA ONLY)
  TTree *outPT=new TTree("protons","protons");
  outPT->Branch("Run",&ev.run,"Run/i");
  outPT->Branch("EventNum",&ev.event,"EventNum/l");
  outPT->Branch("LumiSection",&ev.lumi,"LumiSection/i");
  outPT->Branch("nvtx",&ev.nvtx,"nvtx/I");
  outPT->Branch("rho",&ev.rho,"rho/F");
  outPT->Branch("nchPV",&ev.nchPV,"nchPV/I");
  outPT->Branch("CrossingAngle",&ev.beamXangle,"CrossingAngle/F");
  
  outPT->Branch("nRecoProtCand",&nRecoProtCand,"nRecoProtCand/I");
  outPT->Branch("ProtCand_xi",pr_xi,"ProtCand_xi[nRecoProtCand]/F");
  outPT->Branch("ProtCand_t",pr_t,"ProtCand_t[nRecoProtCand]/F");
  //outPT->Branch("ProtCand_ThX",pr_ThX,"ProtCand_ThX[nRecoProtCand]/F");
  //outPT->Branch("ProtCand_ThY",pr_ThY,"ProtCand_ThY[nRecoProtCand]/F");
  outPT->Branch("ProtCand_rpid",pr_rpid,"ProtCand_rpid[nRecoProtCand]/S");
  //outPT->Branch("ProtCand_arm",pr_arm,"ProtCand_arm[nRecoProtCand]/S");
  outPT->Branch("ProtCand_ismultirp",pr_ismultirp,"ProtCand_ismultirp[nRecoProtCand]/S");
	

  //Variables (bools, ints, floats):
  TString bvars[]={"hasMTrigger","hasETrigger","hasEMTrigger","hasMETrigger","hasMMTrigger",
                   "hasEETrigger","hasSLT","hasDLT","isOS","isOF"};
  std::map<TString,bool> boutVars;
  for(size_t i=0; i<sizeof(bvars)/sizeof(TString); i++) boutVars[bvars[i]]=false;
  
  TString ivars[]={"nBjets","nJets","Lep0Id","Lep1Id"};
  std::map<TString,Int_t> ioutVars;
  for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++) ioutVars[ivars[i]]=0;
  
  TString fvars[]={
				   //generator info (MC only)
				   "gen_mll","gen_Yll","p1_xi_truth","p2_xi_truth",
				   
				   //lepton recontsruction
				   "Lep0Pt","Lep1Pt","Lep0Eta","Lep1Eta","Lep0Phi","Lep1Phi",
				   "InvMass","Yll","pTll","Acopl",
				   
				   //weights:
				   "selSF_wgt", "trigSF_wgt", "prefire_prob", "pu_wgt","ptag_wgt","ptag_wgt_err"
					
				  };
  std::map<TString,Float_t> foutVars;
  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++) foutVars[fvars[i]]=0.;
  
  // Lepton variables:
  const int MAXLEP = 10;
  Int_t nLepCand;
  Float_t LepCand_pt[MAXLEP], LepCand_eta[MAXLEP], LepCand_phi[MAXLEP], LepCand_e[MAXLEP]; 
  Short_t LepCand_charge[MAXLEP], LepCand_id[MAXLEP]; 
    
  // BOOK OUTPUT TREE
  TTree *outT=new TTree("tree","tree");
  outT->Branch("Run",&ev.run,"Run/i");
  outT->Branch("EventNum",&ev.event,"EventNum/l");
  outT->Branch("LumiSection",&ev.lumi,"LumiSection/i");
  outT->Branch("nvtx",&ev.nvtx,"nvtx/I");
  outT->Branch("rho",&ev.rho,"rho/F");
  outT->Branch("nchPV",&ev.nchPV,"nchPV/I");
  outT->Branch("zPV2",&ev.zPV2,"zPV2/F");
  outT->Branch("sumPVChPt",&ev.sumPVChPt,"sumPVChPt/F");
  outT->Branch("sumPVChPz",&ev.sumPVChPz,"sumPVChPz/F");
  outT->Branch("sumPVChHt",&ev.sumPVChHt,"sumPVChHt/F");
  outT->Branch("nchPV_v",ev.nchPV_v,"nchPV[8]/I");
  outT->Branch("sumPVChPt_v",ev.sumPVChPt_v,"sumPVChPt_v[8]/F");
  outT->Branch("sumPVChPz_v",ev.sumPVChPz_v,"sumPVChPz_v[8]/F");
  outT->Branch("sumPVChHt_v",ev.sumPVChHt_v,"sumPVChHt_v[8]/F");
  outT->Branch("CrossingAngle",&ev.beamXangle,"CrossingAngle/F");

  outT->Branch("nRecoProtCand",&nRecoProtCand,"nRecoProtCand/I");
  outT->Branch("ProtCand_xi",pr_xi,"ProtCand_xi[nRecoProtCand]/F");
  outT->Branch("ProtCand_t",pr_t,"ProtCand_t[nRecoProtCand]/F");
  //outT->Branch("ProtCand_ThX",pr_ThX,"ProtCand_ThX[nRecoProtCand]/F");
  //outT->Branch("ProtCand_ThY",pr_ThY,"ProtCand_ThY[nRecoProtCand]/F");
  outT->Branch("ProtCand_rpid",pr_rpid,"ProtCand_rpid[nRecoProtCand]/S");
  //outT->Branch("ProtCand_arm",pr_arm,"ProtCand_arm[nRecoProtCand]/S");
  outT->Branch("ProtCand_ismultirp",pr_ismultirp,"ProtCand_ismultirp[nRecoProtCand]/S");

  outT->Branch("nLepCand",&nLepCand,"nLepCand/I");
  outT->Branch("LepCand_pt",LepCand_pt,"LepCand_pt[nLepCand]/F");
  outT->Branch("LepCand_eta",LepCand_eta,"LepCand_eta[nLepCand]/F");
  outT->Branch("LepCand_phi",LepCand_phi,"LepCand_phi[nLepCand]/F");
  outT->Branch("LepCand_e",LepCand_e,"LepCand_e[nLepCand]/F");
  outT->Branch("LepCand_charge",LepCand_charge,"LepCand_charge[nLepCand]/S");
  outT->Branch("LepCand_id",LepCand_id,"LepCand_id[nLepCand]/S");
		
	
  float evt_weight=0;
  outT->Branch("weight",&evt_weight,"weight/F");
  
  ADDVAR(&ev.met_pt,"met_pt","/F",outT);
  ADDVAR(&ev.met_phi,"met_phi","/F",outT);
  ADDVAR(&ev.met_sig,"met_sig","/F",outT);

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
	
  outT->SetDirectory(fOut);
  ADDVAR(&(ioutVars["nJets"]),"nJets","/I",outPT);
  ADDVAR(&(ioutVars["nBjets"]),"nBjets","/I",outPT);
  outPT->Branch("nLepCand",&nLepCand,"nLepCand/I");
  

  //BOOK HISTOGRAMS  
  HistTool ht;
  ht.setNsyst(0);
  // normalization and event count
  ht.addHist("evt_count",    new TH1F("evt_count",   ";Selection Stage;Events",10,0,10));
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(1,"Total");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(2,"Sumweighted");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(3,"accepted");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(4,"nothing (data)");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(5,"nothing");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(6,"trigger");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(7,"#geq2 lep");
  ht.getPlots()["evt_count"]->SetBinContent(1,counter->GetBinContent(1));
  ht.getPlots()["evt_count"]->SetBinContent(2,counter->GetBinContent(2));
  ht.getPlots()["evt_count"]->SetBinContent(3,counter->GetBinContent(3));
	  
  std::cout << "initialization done" << std::endl;

  //EVENT SELECTION WRAPPER (GETS LISTS OF PHYSICS OBJECTS FROM THE INPUT)
  SelectionTool selector(in_fname, false, triggerList);

  // JEC/JER settings
  int sys = 0;
  if(systVar.find("jetUp")!=string::npos) sys = 1;
  if(systVar.find("jetDn")!=string::npos) sys = -1;
  if(sys==1){cout << "Running JEC/JER up variation"<<endl;}
  else if(sys==-1){cout << "Running JEC/JER down variation"<<endl;}
  else{cout << "Running nominal jet callibration"<<endl;}
	
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  LOOP OVER EVENTS  /////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
		
	  //Reset vars:
      for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++) foutVars[fvars[i]]=0.;
	  foutVars["pu_wgt"] = foutVars["ptag_wgt"] = foutVars["ptag_wgt_err"] = 1;
	  
      ////////////////////
      // EVENT WEIGHTS //
      ///////////////////
  	  evt_weight = (ev.g_nw>0 ? ev.g_w[0] : 1.0);
	  std::vector<double>plotwgts(1,1.);
	  
	  ht.fill("evt_count", 3, plotwgts); // count all events before any selection
		
      //assign randomly a run period
      TString period = lumi.assignRunPeriod();	  

      //////////////////
      // CORRECTIONS //
      /////////////////
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);
      jec.smearJetEnergies(ev);


      //trigger
      boutVars["hasMTrigger"] = (selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits) ); 
	  boutVars["hasETrigger"] = (selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v",                                  ev.triggerBits));
      boutVars["hasMMTrigger"]= (selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                  ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",          ev.triggerBits));
	  boutVars["hasEETrigger"]= (selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev.triggerBits) );
	  boutVars["hasMETrigger"]= (selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits));
	  boutVars["hasEMTrigger"]= (selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits) ||
								 selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits));
      
	  boutVars["hasSLT"] = (boutVars["hasMTrigger"] || boutVars["hasETrigger"]); 
	  boutVars["hasDLT"] = (boutVars["hasMMTrigger"]||boutVars["hasEETrigger"]||boutVars["hasEMTrigger"]||boutVars["hasMETrigger"]);
      
	  if(debug) cout << "runNumber #"<<ev.run<<", event = " << ev.event<<", pass trigger (hasSLT) = "  << boutVars["hasSLT"] <<endl;

      //lepton selection
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);  
      std::vector<Particle> selectedLeptons;	  
      SelectionTool::QualityFlags muId(SelectionTool::TIGHT);
      leptons = selector.selLeptons(leptons,muId,SelectionTool::MVA90,minLeptonPt,2.4);

      // selection of leptons
	  nLepCand = 0;
      for( size_t i_lept=0;i_lept<leptons.size();i_lept++) {	  
	    if (leptons[i_lept].pt()<30.) continue;
		if (leptons[i_lept].id()==11 && fabs(leptons[i_lept].eta())>2.1) continue;
		// if (leptons[i_lept].reliso()>0.10) continue;   //usually tighter
		selectedLeptons.push_back(leptons[i_lept]);
		LepCand_pt[nLepCand] = leptons[i_lept].pt();
		LepCand_eta[nLepCand] = leptons[i_lept].eta();
		LepCand_phi[nLepCand] = leptons[i_lept].phi();
		LepCand_e[nLepCand] = leptons[i_lept].e();
		LepCand_charge[nLepCand] = leptons[i_lept].charge();
		LepCand_id[nLepCand] = leptons[i_lept].id();
		nLepCand++;
	  }

      // selection of protons
      nRecoProtCand = 0;
	  for (int ift=0; ift<ev.nfwdtrk; ift++) {
		  float xi = ev.fwdtrk_xi[ift];
		  if(xi<XIMIN || xi>XIMAX) continue;
		  const unsigned short pot_raw_id = ev.fwdtrk_pot[ift];
		  
		  pr_xi[nRecoProtCand] = xi;
		  pr_t[nRecoProtCand] = ev.fwdtrk_t[ift];
		  //pr_ThX[nRecoProtCand] = ev.fwdtrk_thetax[ift];
		  //pr_ThY[nRecoProtCand] = ev.fwdtrk_thetay[ift];
		  
		  pr_rpid[nRecoProtCand] = pot_raw_id;
		  //pr_arm[nRecoProtCand] = pot_raw_id / 100;
		  pr_ismultirp[nRecoProtCand] = ev.fwdtrk_method[ift];
		  nRecoProtCand++;
      }
	  
	  // ---- EVENT SELECTION --------------------------------------------------------------
	  //if ( ev.isData && ((foutVars["p1_xi"] ==0 ) && (foutVars["p2_xi"] == 0)) && SKIMME )        continue; // ONLY events with >0 protons
	  ht.fill("evt_count", 4, plotwgts); // count events after selection of protons
	  
	  if(!boutVars["hasSLT"] && !boutVars["hasDLT"] && SKIMME)   continue; // events with electrons (id=11) or muons (id=13)
	  ht.fill("evt_count", 5, plotwgts); // count events after channel selection
		
	  if (nLepCand<2){
		  outPT->Fill();
	  }
	  
	  if (nLepCand<2 && SKIMME) continue; // more than 2 selected leptons
      ht.fill("evt_count", 6, plotwgts); // count events after selection on number of leptons (SHOULD BE SAME)
	  
	  
	  // Proceed with event processing
	  
	  // Lepton variables:
	  boutVars["isOS"] = (selectedLeptons[0].charge()==-selectedLeptons[1].charge());
	  boutVars["isOF"] = (selectedLeptons[0].id()+selectedLeptons[1].id());
      foutVars["InvMass"] =  (selectedLeptons[0].p4()+selectedLeptons[1].p4()).M();
      foutVars["Yll"] =  (selectedLeptons[0].p4()+selectedLeptons[1].p4()).Rapidity();
      foutVars["pTll"] = (selectedLeptons[0].p4()+selectedLeptons[1].p4()).Pt();
      foutVars["Acopl"] = 1 - selectedLeptons[0].p4().DeltaPhi(selectedLeptons[1].p4())/TMath::Pi();
	  
	  foutVars["Lep0Pt"] = selectedLeptons[0].pt();
	  foutVars["Lep0Eta"] = selectedLeptons[0].eta();
	  foutVars["Lep0Phi"] = selectedLeptons[0].phi();
	  ioutVars["Lep0Id"] = selectedLeptons[0].id();
	  
	  foutVars["Lep1Pt"] = selectedLeptons[1].pt();
	  foutVars["Lep1Eta"] = selectedLeptons[1].eta();
	  foutVars["Lep1Phi"] = selectedLeptons[1].phi();	  
	  ioutVars["Lep1Id"] = selectedLeptons[1].id();	  
					
      //select jets 
      std::vector<Jet> allJets = selector.getGoodJets(ev,minJetPt,2.4,leptons,{});

	  foutVars["HT"] = 0;
	  foutVars["nJets"] = allJets.size();
      for(size_t ij=0; ij<allJets.size(); ij++) {
		  foutVars["HT"]+=allJets[ij].pt();
		  if(ev.j_btag[allJets[ij].getJetIndex()]) {
			  ioutVars["nBjets"]++;
		  }
      }
	  	
       //MC specific: compute event weight and get truth info
      if(!ev.isData){

        TString period = lumi.assignRunPeriod();
        EffCorrection_t selSF1 = lepEffH.getOfflineCorrection(selectedLeptons[0], period);
        EffCorrection_t selSF2 = lepEffH.getOfflineCorrection(selectedLeptons[1], period);
		EffCorrection_t trigSF = lepEffH.getTriggerCorrection(selectedLeptons,{},{},period);
        EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});
		
		foutVars["selSF_wgt"] = selSF1.first*selSF2.first;
		foutVars["trigSF_wgt"] = trigSF.first;
		foutVars["prefire_prob"] = l1prefireProb.first;
		
        evt_weight  *= trigSF.first*selSF1.first*selSF2.first*l1prefireProb.first;
				
		// truth information (if available)
		TLorentzVector lep1, lep2;
	    if(debug) cout << "ngtop="<<ev.ngtop<<endl;
		for(int i=0;i<ev.ngtop;i++){
			if(ev.gtop_id[i]==11 || ev.gtop_id[i]==13)  {
				lep1.SetPtEtaPhiM(ev.gtop_pt[i],ev.gtop_eta[i],ev.gtop_phi[i],ev.gtop_m[i]);
			}
			if(ev.gtop_id[i]==-11 || ev.gtop_id[i]==-13)  {
				lep2.SetPtEtaPhiM(ev.gtop_pt[i],ev.gtop_eta[i],ev.gtop_phi[i],ev.gtop_m[i]);
			}
			if(ev.gtop_id[i]==2212){ // stable proton with xi < XIMAX (small momentum loss)
				if(ev.gtop_pz[i]>(1-XIMAX)*EBEAM) {foutVars["p1_xi_truth"] = 1 - ev.gtop_pz[i]/EBEAM;}
				else if(ev.gtop_pz[i]<(XIMAX-1)*EBEAM) {foutVars["p2_xi_truth"] = 1 + ev.gtop_pz[i]/EBEAM; }
			}
		}
		if(lep1.Pt() && lep2.Pt()){ // check if found pair of leptons
			float mll = (lep1+lep2).M();
			float Yll = (lep1+lep2).Rapidity();
			foutVars["gen_mll"] = mll;
			foutVars["gen_Yll"] = Yll;		
		} // end if dilep
	  } // end else isData (MC specific)
	  else{ foutVars["selSF_wgt"] = foutVars["trigSF_wgt"] = foutVars["prefire_prob"] = 1;}
	  
      

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
  if(isData) outPT->Write();
  
  fOut->Close();
}
