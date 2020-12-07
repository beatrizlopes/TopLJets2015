#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/LowMu2020.h"
#include "TopLJets2015/TopAnalysis/interface/NeutrinoEllipseCalculator.h"
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

//
void RunLowMu2020(const TString in_fname,
                      TString outname,
                      TH1F *genPU, 
                      TString era,
                      bool debug) 
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  //const char* CMSSW_BASE = getenv("CMSSW_BASE");
  MiniEvent_t ev;

  bool isData = in_fname.Contains("Data13TeV") || in_fname.Contains("SingleMuon") || in_fname.Contains("HighEGJet");
  //bool isPythia8 = in_fname.Contains("pythia8");
  
  //object selection cuts
  float minLeadJetPt(145);
  float minLeadLeptonPt(17);
  float minJetPt(25);
  float minLeptonPt(8);
  
   
  //CORRECTIONS: LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper lepEffH(isData,era);


  //CORRECTIONS: B-TAG CALIBRATION
  BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
  
  //JEC/JER
  JECTools jec(era);
 
  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;
	
  //READ TREE FROM FILE
  TFile *f = TFile::Open(in_fname);
  if(f==NULL || f->IsZombie()) {
    cout << "Corrupted or missing file " << in_fname << endl;
    return;
  }
  TH1 *counter=(TH1 *)f->Get("analysis/counter");
  if(!counter) {cout << "Corrupted or missing counter: \"analysis/counter\" " << endl;return;}
  TH2F *RPcount=(TH2F *)f->Get("analysis/RPcount");
  if(!RPcount) {cout << "Corrupted or missing RPcount: \"analysis/RPcount\" " << endl;return;}
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/tree");
  attachToMiniEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  
  if (debug) nentries = min(1000,nentries); //restrict number of entries for testing
  t->GetEntry(0);
  cout << "...producing " << outname << " from " << nentries << " events" << endl;  
    
  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  
  // BOOK PROTON TREE (DATA ONLY)
  TTree *outPT=new TTree("protons","protons");
  outPT->Branch("run",&ev.run,"run/i");
  outPT->Branch("event",&ev.event,"event/l");
  outPT->Branch("lumi",&ev.lumi,"lumi/i");
  outPT->Branch("nvtx",&ev.nvtx,"nvtx/I");
  outPT->Branch("rho",&ev.rho,"rho/F");
  outPT->Branch("nchPV",&ev.nchPV,"nchPV/I");
  outPT->Branch("beamXangle",&ev.beamXangle,"beamXangle/F");
  float m_protonVars_p1_xi=0, m_protonVars_p2_xi=0;
  outPT->Branch("p1_xi",&m_protonVars_p1_xi);
  outPT->Branch("p2_xi",&m_protonVars_p2_xi);

  // BOOK OUTPUT TREE
  TTree *outT=new TTree("tree","tree");
  outT->Branch("run",&ev.run,"run/i");
  outT->Branch("event",&ev.event,"event/l");
  outT->Branch("lumi",&ev.lumi,"lumi/i");
  outT->Branch("nvtx",&ev.nvtx,"nvtx/I");
  outT->Branch("rho",&ev.rho,"rho/F");
  outT->Branch("nchPV",&ev.nchPV,"nchPV/I");
  outT->Branch("beamXangle",&ev.beamXangle,"beamXangle/F");
	
  // Parton shower weights
  if(!isData){
	outT->Branch("g_npsw",    &ev.g_npsw,   "g_npsw/I");
	outT->Branch("g_psw",      ev.g_psw,    "g_psw[g_npsw]/F");
  }

  ADDVAR(&ev.met_pt,"met_pt","/F",outT);
  ADDVAR(&ev.met_phi,"met_phi","/F",outT);
  ADDVAR(&ev.met_sig,"met_sig","/F",outT);
	
  //Variables (bools, ints, floats):
  TString bvars[]={
					"HLT_HIMu15","HLT_HIEle15","HLT_HIPFJetFwd140","HLT_HIPFJet140",
					
					"passLepSel", "passJetSel", "passTopSel"
  };

				  
  TString ivars[]={
					"nbj", "nlj", "np"
  };
 
  TString fvars[]={ 
					"boson_m","boson_dphi","boson_pt","boson_Y",
					"ttm","ttdphi","ttpt","ttY",
					"j2m","j2dphi","j2pt","j2Y",
					"j3m","j3thrust","j3pt","j3Y",
					"mpp","Ypp","xip","xin",
					
					"HT","mWT",
					
					"weight","gen_wgt","pu_wgt","ptag_wgt","ptag_wgt_err"
  };
  
  int jet_n(0);
  float jet_pt[ev.MAXJET];
  float jet_deepcsv[ev.MAXJET];
  
  outT->Branch("jet_n",&jet_n,"jet_n/I");
  outT->Branch("jet_pt",jet_pt,"jet_pt[jet_n]/F");
  outT->Branch("jet_deepcsv",jet_deepcsv,"jet_deepcsv[jet_n]/F");
	
  int lep_n(0);
  int lep_id[ev.MAXLEP];
  int lep_q[ev.MAXLEP];
  float lep_pt[ev.MAXLEP];

  outT->Branch("lep_n",&lep_n,"lep_n/I");
  outT->Branch("lep_pt",lep_pt,"lep_pt[lep_n]/F");
  outT->Branch("lep_id",lep_id,"lep_id[lep_n]/I");
  outT->Branch("lep_q",lep_q,"lep_q[lep_n]/I");
  
  // fill custom variables
  std::map<TString,bool> boutVars;
  for(size_t i=0; i<sizeof(bvars)/sizeof(TString); i++){
    boutVars[bvars[i]]=false;
    ADDVAR(&(boutVars[bvars[i]]),bvars[i],"/O",outT);
  }

  std::map<TString,Int_t> ioutVars;
  for(size_t i=0; i<sizeof(ivars)/sizeof(TString); i++){
    ioutVars[ivars[i]]=0;
	ADDVAR(&(ioutVars[ivars[i]]),ivars[i],"/I",outT);
  }
	
  std::map<TString,Float_t> foutVars;
  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++){
    foutVars[fvars[i]]=0.;
	ADDVAR(&(foutVars[fvars[i]]),fvars[i],"/F",outT);
  }	

  //PPS extra info
  int pixel_pos_n(0), pixel_neg_n(0);
  float pixel_pos_xi[ev.MAXPROTONS], pixel_neg_xi[ev.MAXPROTONS];
  float strip_pos_xi, strip_neg_xi;

  outT->Branch("pixel_pos_n",&pixel_pos_n,"pixel_pos_n/I");
  outT->Branch("pixel_neg_n",&pixel_neg_n,"pixel_neg_n/I");
  outT->Branch("pixel_pos_xi",pixel_pos_xi,"pixel_pos_xi[pixel_pos_n]/F");
  outT->Branch("pixel_neg_xi",pixel_neg_xi,"pixel_neg_xi[pixel_neg_n]/F");
  outT->Branch("strip_pos_xi",&strip_pos_xi,"strip_pos_xi/F");
  outT->Branch("strip_neg_xi",&strip_neg_xi,"strip_neg_xi/F");
	
  outT->SetDirectory(fOut);
  
  //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  
  // normalization and event count
  ht.addHist("evt_count",    new TH1F("evt_count",   ";Selection Stage;Events",10,0,10));
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(1,"Total");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(2,"Sumweighted");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(3,"preselection");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(4,"trigger");
  ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(5,"lep or jet");
  ht.getPlots()["evt_count"]->SetBinContent(1,counter->GetBinContent(1));
  ht.getPlots()["evt_count"]->SetBinContent(2,counter->GetBinContent(2));
  ht.getPlots()["evt_count"]->SetBinContent(3,counter->GetBinContent(3));

  // proton count
  ht.addHist("pn_count",    new TH1F("pn_count",   ";;Events",10,0,10));
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(1,"00 (no protons)");
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(2,"10 (RP0 only)");
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(3,"01 (RP1 only)");
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(4,"11 (both RP)");
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(5,"20 (RP0 saturated, RP1 no hits)");
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(6,"02 (RP1 saturated, RP0 no hits)");
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(7,"21 (RP0 saturated, RP1 hit)");
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(8,"12 (RP1 saturated, RP0 hit)");
  ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(9,"22 (both saturated)");

  int nbin = 0;
  nbin=RPcount->FindBin(0,0); ht.getPlots()["pn_count"]->SetBinContent(1,RPcount->GetBinContent(nbin));
  nbin=RPcount->FindBin(1,0); ht.getPlots()["pn_count"]->SetBinContent(2,RPcount->GetBinContent(nbin));
  nbin=RPcount->FindBin(0,1); ht.getPlots()["pn_count"]->SetBinContent(3,RPcount->GetBinContent(nbin));
  nbin=RPcount->FindBin(1,1); ht.getPlots()["pn_count"]->SetBinContent(4,RPcount->GetBinContent(nbin));
  nbin=RPcount->FindBin(2,0); ht.getPlots()["pn_count"]->SetBinContent(5,RPcount->GetBinContent(nbin));
  nbin=RPcount->FindBin(0,2); ht.getPlots()["pn_count"]->SetBinContent(6,RPcount->GetBinContent(nbin));
  nbin=RPcount->FindBin(2,1); ht.getPlots()["pn_count"]->SetBinContent(7,RPcount->GetBinContent(nbin));
  nbin=RPcount->FindBin(1,1); ht.getPlots()["pn_count"]->SetBinContent(8,RPcount->GetBinContent(nbin));
  nbin=RPcount->FindBin(2,2); ht.getPlots()["pn_count"]->SetBinContent(9,RPcount->GetBinContent(nbin));

  std::cout << "init done" << std::endl;
  

  //EVENT SELECTION WRAPPER (GETS LISTS OF PHYSICS OBJECTS FROM THE INPUT)
  SelectionTool selector(in_fname, false, triggerList, SelectionTool::AnalysisType::LOWMU);

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }

      ////////////////////
      // EVENT WEIGHTS //
      ///////////////////
      float wgt(1.0);
      std::vector<double>plotwgts(1,wgt);
	  
      //trigger
      boutVars["HLT_HIMu15"]=(selector.hasTriggerBit("HLT_HIMu15_v",  ev.addTriggerBits) );   
      boutVars["HLT_HIEle15"]=(selector.hasTriggerBit("HLT_HIEle15_WPLoose_Gsf_v", ev.addTriggerBits) );   
      boutVars["HLT_HIPFJet140"]=(selector.hasTriggerBit("HLT_HIPFJet140_v", ev.addTriggerBits) );   
      boutVars["HLT_HIPFJetFwd140"]=(selector.hasTriggerBit("HLT_HIPFJetFwd140_v", ev.addTriggerBits) );   
	  
	  bool passTrigger = selector.passSingleLeptonTrigger(ev) || selector.passJetTrigger(ev);
      if(!passTrigger) continue;
      ht.fill("evt_count", 3, plotwgts);  // count events after trigger cut
	  
      //////////////////////
      // OBJECT SELECTION //
      //////////////////////
      //lepton selection
      std::vector<Particle> bare_leptons = selector.flaggedLeptons(ev, minLeptonPt, 2.5);   
      bare_leptons = selector.selLeptons(bare_leptons,SelectionTool::TIGHT,SelectionTool::MVA90,minLeptonPt,2.5);
      std::vector<Particle> leptons;
      for( size_t i_lept=0;i_lept<bare_leptons.size();i_lept++) {
		if (bare_leptons[i_lept].id()==11 && fabs(bare_leptons[i_lept].eta())>2.1) continue;
        //        if (leptons[i_lept].reliso()>0.10) continue;   //usually tighter
        leptons.push_back(bare_leptons[i_lept]);
      }	  
	  
      //if(leptons.size()==0) continue;
      //if(leptons[0].id()!=13) continue;

      //select jets
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);   
      jec.smearJetEnergies(ev);	  
      std::vector<Jet> jets = selector.getGoodJets(ev,minJetPt,4.7,leptons,{});
      
	  boutVars["passLepSel"] = selector.passSingleLeptonTrigger(ev) 
								&& (leptons.size()>0 && leptons[0].pt()>minLeadLeptonPt);
	  boutVars["passJetSel"] = selector.passJetTrigger(ev) 
							    && (jets.size()>0 && jets[0].pt()>minLeadJetPt);
	  
      if(!boutVars["passJetSel"] && !boutVars["passLepSel"]) continue;
      ht.fill("evt_count", 4, plotwgts);  // count events after object selection cut

      std::vector<Jet>      bJets,lightJets;
      for(size_t ij=0; ij<jets.size(); ij++) {
          if(jets[ij].flavor()==5) bJets.push_back(jets[ij]);
          else                     lightJets.push_back(jets[ij]);
      }
	  	  
      // selection of protons
      double p1_xi =0.; // proton in positive pot
      double p2_xi =0.; // proton in negative pot
	  pixel_pos_n = pixel_neg_n = 0;
	  strip_pos_xi = strip_neg_xi = 0.;
      for (int ift=0; ift<ev.nfwdtrk; ift++) {
        const unsigned short pot_raw_id = ev.fwdtrk_pot[ift];
        if(ev.fwdtrk_method[ift]==1){  // selecting only MultiRP protons
          if (pot_raw_id<100){ // positive z  (pot_raw_id=3)
            p1_xi = ev.fwdtrk_xi[ift];
          }
          else {   // negative z   (pot_raw_id=103)
            p2_xi = ev.fwdtrk_xi[ift];
          }
        }
		else{
          if (pot_raw_id==3) strip_pos_xi = ev.fwdtrk_xi[ift];
          else if (pot_raw_id==103) strip_neg_xi = ev.fwdtrk_xi[ift];
          else if (pot_raw_id==23) pixel_pos_xi[pixel_pos_n++] = ev.fwdtrk_xi[ift];
          else if (pot_raw_id==123) pixel_neg_xi[pixel_neg_n++] = ev.fwdtrk_xi[ift];
  	    }
      }

	  // Store proton pool at preselection
	  if(ev.isData){
		m_protonVars_p1_xi = p1_xi;
		m_protonVars_p2_xi = p2_xi;
		outPT->Fill();
  	  }	
	  
      //event weight
      float evWgt(1.0);


      //////////////////////
      // RECONSTRUCTION   //
      //////////////////////

	  foutVars["weight"] = foutVars["gen_wgt"] = 1;
      
      //MC specific: compute event weight
      if (!ev.isData) {

        //TString period = lumi.assignRunPeriod();
        //double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
		
		//EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons,{},{},period);
        //EffCorrection_t selSF = lepEffH.getOfflineCorrection(leptons[0], period);

        //evWgt  = trigSF.first*selSF.first;
        foutVars["gen_wgt"] = (ev.g_nw>0 ? ev.g_w[0] : 1.0);        
        evWgt *= foutVars["gen_wgt"];
      }
      
      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
	  	  
      //determine neutrino kinematics (in leptonic events)
	  float nupz = 0;
	  if(boutVars["passLepSel"]){
		  neutrinoPzComputer.SetMET(met);
	      neutrinoPzComputer.SetLepton(leptons[0].p4());
		  nupz = neutrinoPzComputer.Calculate();
	  }
      TLorentzVector neutrino(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));

	  // set boson kinematics
	  TLorentzVector boson(0,0,0,0);
	  if(leptons.size()==1) boson = (leptons[0] + neutrino);
	  if(leptons.size()==2) boson = (leptons[0] + leptons[1]);
      
      // compute scalar ht and jets
      float scalarht(0.);
      jet_n   = int(jets.size());
      for(size_t ij=0; ij<jets.size(); ij++) {
        jet_pt[ij] = jets[ij].Pt();
        scalarht += jets[ij].pt();
      }
      if(boutVars["passLepSel"]) scalarht += leptons[0].pt();
      scalarht += neutrino.Pt();

	  // check HF selection
  	  TLorentzVector top1(0,0,0,0), top2(0,0,0,0);
	  boutVars["passTopSel"] =  boutVars["passLepSel"] && bJets.size()>1 && lightJets.size()>1;
      if(boutVars["passTopSel"]){
		  top1 += leptons[0] + neutrino;
		  if(bJets[0].DeltaR(leptons[0])<bJets[1].DeltaR(leptons[0])) {
			  top1+=bJets[0];
			  top2=bJets[1]+lightJets[0]+lightJets[1];
		  }
		  else{
			  top1+=bJets[1];
			  top2=bJets[0]+lightJets[0]+lightJets[1];
		  }
	  }
	  
	  // leptons
      lep_n   = int(leptons.size());
	  for(size_t ij=0; ij<leptons.size(); ij++) {
        lep_pt[ij] = leptons[ij].Pt();
        lep_id[ij] = leptons[ij].id();
        lep_q[ij] = leptons[ij].charge();
	  }
	  
      //////////////////////
      // FILL KINEMATICS  //
      //////////////////////
	  ioutVars["nbj"] = bJets.size();
	  ioutVars["nlj"] = lightJets.size();
	  ioutVars["np"] = (p1_xi>0) + (p2_xi>0);
	  
	  foutVars["xip"] = p1_xi;
	  foutVars["xin"] = p2_xi;
	  if( (p1_xi>0) && (p2_xi>0)){
		foutVars["mpp"] = 13000.*sqrt(p1_xi*p2_xi);
		foutVars["Ypp"] = 0.5*TMath::Log(p1_xi/p2_xi);
	  }
	  else{ foutVars["mpp"] = 0; foutVars["Ypp"] = 999;}
	  foutVars["boson_m"] = boson.M();
	  foutVars["boson_pt"] = boson.Pt();
	  foutVars["boson_Y"] = boson.Rapidity();
	  foutVars["mWT"] = (lep_n==1) ? sqrt(2*leptons[0].Pt()*met.Pt()*(1-cos(leptons[0].DeltaPhi(met)))) : 0;

	  if(lep_n==1) foutVars["boson_dphi"] = leptons[0].DeltaPhi(neutrino);
	  else if(lep_n==2) foutVars["boson_dphi"] = leptons[0].DeltaPhi(leptons[1]);
	  else foutVars["boson_dphi"] = 999;
	  
	  TLorentzVector ttbar = top1 + top2;
	  foutVars["ttm"] = ttbar.M();
	  foutVars["ttpt"] = ttbar.Pt();
	  foutVars["ttY"] = ttbar.Rapidity();
	  foutVars["ttdphi"] = boutVars["passTopSel"] ? top1.DeltaPhi(top2) : 999;

      if(jet_n==2){
		TLorentzVector v_jets = (jets[0]+jets[1]);
	    foutVars["j2m"] = v_jets.M();
	    foutVars["j2pt"] = v_jets.Pt();
	    foutVars["j2Y"] = v_jets.Rapidity();
	    foutVars["j2dphi"] = jets[0].DeltaPhi(jets[1]);
	  }
	  else if(jet_n==3){
		TLorentzVector v_jets = (jets[0]+jets[1]+jets[2]);
	    foutVars["j3m"] = v_jets.M();
	    foutVars["j3pt"] = v_jets.Pt();
	    foutVars["j3Y"] = v_jets.Rapidity();
	    foutVars["j3thrust"] = 999;
	  }
	  else{
	    foutVars["j2m"] = foutVars["j3m"] = 0;
	    foutVars["j2pt"] = foutVars["j3pt"] = 0;
	    foutVars["j2Y"] = foutVars["j3Y"] = 999;
	    foutVars["j2dphi"] = foutVars["j3thrust"] = 999;
	  }								
      foutVars["HT"] = scalarht;
      foutVars["weight"] = evWgt;

	  // Save output tree
	  outT->Fill();
 
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
  outT->Write();
  if(isData) outPT->Write();
  fOut->Close();
}
