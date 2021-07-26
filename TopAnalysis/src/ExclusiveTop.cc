// ROOT includes
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TApplication.h"

// cpp includes
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>

// package includes
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveTop.h"
#include "TopLJets2015/TopAnalysis/interface/NeutrinoEllipseCalculator.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"
#include "TopLJets2015/TopAnalysis/interface/PPSEff.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"

//dilepton top kinematic reconstruction
#include "TopLJets2015/TopAnalysis/interface/KinematicReconstruction.h"

using namespace std;
#define ADDVAR(x,name,t,tree) tree->Branch(name,x,TString(name)+TString(t))

double m_TOP = 173.1;
double m_W   =  80.379;
double m_NU  =  0.;

// ---- SWITCHES ----------------------------------------------------------------
//#define DEBUG_ON                // activate in order to activate additional output for debugging and to limit entries
#define HISTOGRAMS_ON           // comment to avoid creating histograms in root file
#define SAVEERRORS_ON           // comment to avoid saving the quantities error in root file

// Returns the index of the lesser element of sum_squared vector
int dump_index(vector<double> sum_squared){
    int index=0;
    double a=0;
    double b=0;
    a=sum_squared[0];
    for(size_t i=1; i<sum_squared.size() ; i++){
        b=sum_squared[i];
        if(b<=a){
            a=b;
            index=i;
        }
    }
    //cout << "index " << index << endl;
    return index;
}

bool do_kin_reco(std::vector<Particle>& leptons, std::vector<Jet>& jets, std::vector<Jet>& bjets, TLorentzVector& met, double& k\
inReco_ttbar_mass, double& kinReco_ttbar_rapidity) {

    bool hasKinRecoSol = false;
    std::vector<TLorentzVector> allLeptons={};
    std::vector<TLorentzVector> allJets={};
    std::vector<TLorentzVector> allBJets={};

    for (size_t i=0; i < leptons.size(); i++){
      allLeptons.push_back(leptons[i]);
    }

    //Get jets indices                                                                                                                         
    std::vector<int> jetIndices;
    jetIndices.clear();
    for (size_t i=0; i < jets.size(); i++){
      jetIndices.push_back(i);
      allJets.push_back(jets[i]);
    }

    std::vector<int> bjetIndices;
    bjetIndices.clear();
    for (size_t i=0; i < bjets.size(); i++){
      bjetIndices.push_back(i);
      allBJets.push_back(bjets[i]);
    }

    if(jets.size()<2 || bjets.size()<1) {
      hasKinRecoSol = false;
     return hasKinRecoSol;
    }

    //if(debug) std::cout << "now doing ttbar kinematic reconstruction" << std::endl;                                                          
    const KinematicReconstruction* kinematicReconstruction(0);
    //For now random-number-based smearing                                                                                                     
    kinematicReconstruction = new KinematicReconstruction(1, true, false);
    int idx_l1 = leptons[0].charge()>0 ? 1 : 0 ;
    int idx_l2 = leptons[1].charge()>0 ? 1 : 0 ;
    KinematicReconstructionSolutions kinematicReconstructionSolutions  =
    kinematicReconstruction->solutions({idx_l1}, {idx_l2}, jetIndices, bjetIndices, allLeptons, allJets, met);
    const bool hasSolution = kinematicReconstructionSolutions.numberOfSolutions();
    if (hasSolution) {
      hasKinRecoSol = true;

      TLorentzVector recotop = kinematicReconstructionSolutions.solution().top();
      TLorentzVector recoantitop = kinematicReconstructionSolutions.solution().antiTop();
      TLorentzVector recottbar = recotop + recoantitop;

      //kinReco_ttbar_pt = recottbar.Pt();                                                                                                     
      kinReco_ttbar_mass = recottbar.M();
      kinReco_ttbar_rapidity = recotop.Rapidity();

    }
    return hasKinRecoSol;
}

// ---------------------------------------------------------------------------------------------------------------------
//      MAIN
// ---------------------------------------------------------------------------------------------------------------------
void RunExclusiveTop(TString filename,
                     TString outname,
                     Int_t channelSelection,
                     Int_t chargeSelection,
                     TH1F *normH,
                     TH1F *genPU,
                     TString era,
                     Bool_t debug,
					 std::string systVar,
					 int seed
					 )
{
    /////////////////////
    // INITIALIZATION //
    ///////////////////
    const char* CMSSW_BASE = getenv("CMSSW_BASE");
    // CTPPS reconstruction part
    //ctpps::XiReconstructor proton_reco;
    //proton_reco.feedDispersions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/dispersions.txt", CMSSW_BASE));
    
    //ctpps::AlignmentsFactory ctpps_aligns;
    //ctpps_aligns.feedAlignments(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/alignments_30jan2017.txt", CMSSW_BASE));
    
    //ctpps::LHCConditionsFactory lhc_conds;
    //lhc_conds.feedConditions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/xangle_tillTS2.csv", CMSSW_BASE));
    //lhc_conds.feedConditions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/xangle_afterTS2.csv", CMSSW_BASE));
    
    bool isTTbar = filename.Contains("_TTJets");
	bool isData = filename.Contains("Data13TeV") || filename.Contains("SingleMuon") || filename.Contains("SingleElectron");
	//bool isPythia8 = filename.Contains("pythia8");
    
    //PREPARE OUTPUT
    TString baseName=gSystem->BaseName(outname);
    TString dirName=gSystem->DirName(outname);
    TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
    fOut->cd();
    
    //READ TREE FROM FILE
    MiniEvent_t ev;
    TFile *f = TFile::Open(filename);
	TH1 *counter=(TH1 *)f->Get("analysis/counter");
    if(!counter) {cout << "Corrupted or missing counter: \"analysis/counter\" " << endl;return;}
	TH2F *RPcount=(TH2F *)f->Get("analysis/RPcount");
    if(!RPcount) {cout << "Corrupted or missing RPcount: \"analysis/RPcount\" " << endl;return;}
    TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
    TTree *t = (TTree*)f->Get("analysis/tree");
    attachToMiniEventTree(t,ev);
    Int_t nentries(t->GetEntriesFast());
    
    //debug = 1; // manual turn on DEBUG mode
#ifdef DEBUG_ON
//    if (debug) {
        int nentries_cap = 5000; // restricted number of entries for testing (10k+ for BTD)
        if (nentries>nentries_cap) nentries = nentries_cap;
        std::cout << "--- debug mode activated" << std::endl;
//    }
#endif
    
    t->GetEntry(0);
    
    std::cout << "--- producing " << outname << " from " << nentries << " events" << std::endl;
    
    //auxiliary to solve neutrino pZ using MET
    //MEzCalculator neutrinoPzComputer;
    
    //LUMINOSITY+PILEUP
    LumiTools lumi(era,genPU);
    
    //L1-prefire 
    L1PrefireEfficiencyWrapper l1PrefireWR(isData,era);

    //LEPTON EFFICIENCIES
    EfficiencyScaleFactorsWrapper lepEffH(isData,era);
    
    //B-TAG CALIBRATION
    //BTagSFUtil btvSF(era,"DeepCSV",BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
    BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",seed);
	
    // Proton correction class
    PPSEff *PPS_reco = new PPSEff(Form("%s/src/TopLJets2015/TopAnalysis/data/era2017/reco_charactersitics_version1.root", CMSSW_BASE));

    // BOOK PROTON TREE (DATA ONLY)
    fOut->cd();
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
	//if(!isData){
	//  outT->Branch("g_npsw",    &ev.g_npsw,   "g_npsw/I");
	//  outT->Branch("g_psw",      ev.g_psw,    "g_psw[g_npsw]/F");
	//}
	
	/// ISR/FSR Systematic variations, store as a vectors to account comp. split
	// total, [g2gg, g2qq, q2ga, x2xg] X [muR, cNS]. Q=uds, X=cb, NS=non-singular
	// see more info here: https://indico.cern.ch/event/746817/contributions/3101385/attachments/1702410/2742087/psweights_mseidel.pdf
	const int NPSRad_weights = 9;
	float isr_Up[NPSRad_weights], fsr_Up[NPSRad_weights], isr_Down[NPSRad_weights], fsr_Down[NPSRad_weights];
	outT->Branch("isr_Up",isr_Up,Form("isr_Up[%d]/F",NPSRad_weights));
	outT->Branch("fsr_Up",fsr_Up,Form("fsr_Up[%d]/F",NPSRad_weights));
	outT->Branch("isr_Down",isr_Down,Form("isr_Down[%d]/F",NPSRad_weights));
	outT->Branch("fsr_Down",fsr_Down,Form("fsr_Down[%d]/F",NPSRad_weights));
	map<string,int> PSmap;
    // from https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/PSweightsPythia/PythiaPSweightsSettings_cfi.py
	PSmap["nominal"] = 0;
	PSmap["isrDefHi"] = 6;
	PSmap["fsrDefHi"] = 7;
	PSmap["isrDefLo"] = 8;
	PSmap["fsrDefLo"] = 9;
	PSmap["fsr_G2GG_muR_dn"] = 14;
	PSmap["fsr_G2GG_muR_up"] = 15;
	PSmap["fsr_G2QQ_muR_dn"] = 16;
	PSmap["fsr_G2QQ_muR_up"] = 17;
	PSmap["fsr_Q2QG_muR_dn"] = 18;
	PSmap["fsr_Q2QG_muR_up"] = 19;
	PSmap["fsr_X2XG_muR_dn"] = 20;
	PSmap["fsr_X2XG_muR_up"] = 21;
	PSmap["fsr_G2GG_cNS_dn"] = 22;
	PSmap["fsr_G2GG_cNS_up"] = 23;
	PSmap["fsr_G2QQ_cNS_dn"] = 24;
	PSmap["fsr_G2QQ_cNS_up"] = 25;
	PSmap["fsr_Q2QG_cNS_dn"] = 26;
	PSmap["fsr_Q2QG_cNS_up"] = 27;
	PSmap["fsr_X2XG_cNS_dn"] = 28;
	PSmap["fsr_X2XG_cNS_up"] = 29;
	PSmap["isr_G2GG_muR_dn"] = 30;
	PSmap["isr_G2GG_muR_up"] = 31;
	PSmap["isr_G2QQ_muR_dn"] = 32;
	PSmap["isr_G2QQ_muR_up"] = 33;
	PSmap["isr_Q2QG_muR_dn"] = 34;
	PSmap["isr_Q2QG_muR_up"] = 35;
	PSmap["isr_X2XG_muR_dn"] = 36;
	PSmap["isr_X2XG_muR_up"] = 37;
	PSmap["isr_G2GG_cNS_dn"] = 38;
	PSmap["isr_G2GG_cNS_up"] = 39;
	PSmap["isr_G2QQ_cNS_dn"] = 40;
	PSmap["isr_G2QQ_cNS_up"] = 41;
	PSmap["isr_Q2QG_cNS_dn"] = 42;
	PSmap["isr_Q2QG_cNS_up"] = 43;
	PSmap["isr_X2XG_cNS_dn"] = 44;
	PSmap["isr_X2XG_cNS_up"] = 45;
	
	float met_pt, met_phi;
	ADDVAR(&met_pt,"met_pt","/F",outT);
    ADDVAR(&met_phi,"met_phi","/F",outT);
    ADDVAR(&ev.met_sig,"met_sig","/F",outT);	
	
    TString fvars[]={
        "ttbar_m", "ttbar_y",
        // quantities for all objects in event
        "nBjets", "nLightJets", "nJets", "cat",
        "l1_pt", "l1_eta", "l1_phi", "l1_m", "l1_E", "l1_y", "lepton1_isolation",
        "l2_pt", "l2_eta", "l2_phi", "l2_m", "l2_E", "l2_y", "lepton2_isolation",
        "p1_xi", "p2_xi", "ppsSF_wgt", "ppsSF_wgt_err",
		"p1_x","p2_x","p1_y","p2_y",
		"p1_220_x","p2_220_x","p1_220_y","p2_220_y",
        "mpp","ypp",
        "yvis","mll","min_dy","extra_rapidity","extra_rapidity_sum",
		"weight", "gen_wgt", "toppt_wgt", "EL_SF_wgt","MU_SF_wgt", "EL_trigSF_wgt","MU_trigSF_wgt", "L1Prefire_wgt",
		"EL_SF_wgt_err","MU_SF_wgt_err", "EL_trigSF_wgt_err","MU_trigSF_wgt_err", "pu_wgt", "ptag_wgt", "ptag_wgt_err","L1Prefire_wgt_err",
		"ren_err","fac_err",
		"pdf_as","pdf_hs",
        "bJet0_pt","bJet0_eta", "bJet0_phi", "bJet0_m", "bJet0_E","bJet0_y",
        "bJet1_pt","bJet1_eta", "bJet1_phi", "bJet1_m", "bJet1_E","bJet1_y",
        "lightJet0_pt", "lightJet0_eta", "lightJet0_phi", "lightJet0_m", "lightJet0_E",
        "lightJet1_pt", "lightJet1_eta", "lightJet1_phi", "lightJet1_m", "lightJet1_E"
        };
    
    std::map<TString,Float_t> outVars;
    for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++){
        outVars[fvars[i]]=0.;
        ADDVAR(&(outVars[fvars[i]]),fvars[i],"/F",outT);
    }
    ADDVAR(&(outVars["nJets"]),"nJets","/F",outPT);
    ADDVAR(&(outVars["nBjets"]),"nBjets","/F",outPT);

#ifdef HISTOGRAMS_ON
    //BOOK HISTOGRAMS
    HistTool ht;
    ht.setNsyst(0);
    ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
    ht.addHist("ch_tag",       new TH1F("ch_tag",      ";Channel Tag;Events",10,0,10));
    ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",55,-0.5,49.5));
    ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",15,-0.5,14.5));
    ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",10,-0.5,9.5));

    ht.addHist("mttbar_rec",   new TH1F("mttbar_rec",  ";M_{ttbar,rec} [GeV];Events",50,300,1200));
     
    ht.addHist("yttbar_rec", new TH1F("yttbar_rec",";Y_{ttbar,rec} ;Events",75,-2.5,2.5));
    	
    // normalization and event count
    ht.addHist("evt_count",    new TH1F("evt_count",   ";Selection Stage;Events",10,0,10));
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(1,"Total");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(2,"Sumweighted");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(3,"preselection");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(4,"=2 p (data)");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(5,"trigger");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(6,"event cleaning");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(7,"=2 lep");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(8,"#geq2 jets");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(9,"#geq2 bjets");
    ht.getPlots()["evt_count"]->SetBinContent(1,counter->GetBinContent(1));
    ht.getPlots()["evt_count"]->SetBinContent(2,counter->GetBinContent(2));
    ht.getPlots()["evt_count"]->SetBinContent(3,counter->GetBinContent(3));	

    // proton count
    ht.addHist("pn_count",    new TH1F("pn_count",   ";;Events",5,0,5));
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(1,"0 hits");
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(2,"RP0 hit");
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(3,"RP1 hit");
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(4,"both RP");
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(5,"pres + nB#geq1");
	int nbin = 0;
	nbin=RPcount->FindBin(0.5,0.5); ht.getPlots()["pn_count"]->SetBinContent(1,RPcount->GetBinContent(nbin));
	nbin=RPcount->FindBin(1.5,0.5); ht.getPlots()["pn_count"]->SetBinContent(2,RPcount->GetBinContent(nbin));
	nbin=RPcount->FindBin(0.5,1.5); ht.getPlots()["pn_count"]->SetBinContent(3,RPcount->GetBinContent(nbin));
	nbin=RPcount->FindBin(1.5,1.5); ht.getPlots()["pn_count"]->SetBinContent(4,RPcount->GetBinContent(nbin));
	ht.getPlots()["pn_count"]->SetBinContent(5,counter->GetBinContent(4));
#endif
    
    std::cout << "--- init done" << std::endl;
    
    //EVENT SELECTION WRAPPER
    SelectionTool selector(filename, false, triggerList, SelectionTool::AnalysisType::TOP);
    SelectionTool selector_nominal(filename, false, triggerList, SelectionTool::AnalysisType::TOP);
	// selector.minJetPt = 25;
    
	// JEC/JER settings
	int sys = 0;
	if(systVar.find("jerUp")!=string::npos) sys = 1;
	if(systVar.find("jerDn")!=string::npos) sys = -1;
	if(systVar.find("jecUp")!=string::npos) sys = 2;
	if(systVar.find("jecDn")!=string::npos) sys = -2;
    if(systVar.find("jecAbsUp")!=string::npos) sys = 3;
	if(systVar.find("jecAbsDn")!=string::npos) sys = -3;	
    if(systVar.find("jecRelUp")!=string::npos) sys = 4;
	if(systVar.find("jecRelDn")!=string::npos) sys = -4;	
    if(systVar.find("jecPileupUp")!=string::npos) sys = 5;
	if(systVar.find("jecPileupDn")!=string::npos) sys = -5;	
    if(systVar.find("jecFlavUp")!=string::npos) sys = 6;
	if(systVar.find("jecFlavDn")!=string::npos) sys = -6;	
    if(systVar.find("jecHighPtUp")!=string::npos) sys = 7;
	if(systVar.find("jecHighPtDn")!=string::npos) sys = -7;	
    if(systVar.find("jecTimeUp")!=string::npos) sys = 8;
	if(systVar.find("jecTimeDn")!=string::npos) sys = -8;
	if(sys>0){cout << "INFO: Running JEC/JER up variation"<<endl;}
	else if(sys<0){cout << "INFO: Running JEC/JER down variation"<<endl;}
	else{cout << "INFO: Running nominal jet callibration"<<endl;}

	// btagSF settings
	string optionlf = "central", optionhf = "central";
	if(systVar.find("btaghfUp")!=string::npos) {
		cout << "INFO: Running btaghf up variation, seed = "<< seed << endl;
		optionhf = "up";}
	else if(systVar.find("btaghfDn")!=string::npos) {
		cout << "INFO: Running btaghf down variation, seed = "<< seed << endl;
		optionhf = "down";}
        else if(systVar.find("btaglfUp")!=string::npos) {
                cout << "INFO: Running btaglf up variation, seed = "<< seed << endl;
                optionlf = "up";}
        else if(systVar.find("btaglfDn")!=string::npos) {
                cout << "INFO: Running btaglf down variation, seed = "<< seed << endl;
                optionlf = "down";}
	else{cout << "INFO: Running nominal btag SF, seed = "<< seed << endl;}
    
	// Proton reconstuction systematics
	float pps45_err = 0, pps56_err = 0;
	if(systVar.find("pps45Up")!=string::npos) {
		cout << "INFO: Running PPS xi reco. up variation"<<endl;
		pps45_err = 1;}
	else if(systVar.find("pps45Dn")!=string::npos) {
		cout << "INFO: Running PPS xi reco. down variation"<<endl;
		pps45_err = -1;}
	else if(systVar.find("pps56Up")!=string::npos) {
		cout << "INFO: Running PPS xi reco. up variation"<<endl;
		pps56_err = 1;}
	else if(systVar.find("pps56Dn")!=string::npos) {
		cout << "INFO: Running PPS xi reco. down variation"<<endl;
		pps56_err = -1;}
	else{cout << "INFO: Running nominal PPS reco"<<endl;}
	
	// MET uncertanty https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/interface/MET.h#L152-L167
	int met_err = 0;
	if(systVar.find("UnclusteredEnUp")!=string::npos) {
		cout << "INFO: Running MET UnclusteredEn up variation"<<endl;
		met_err = 10;}
	else if(systVar.find("UnclusteredEnDn")!=string::npos) {
		cout << "INFO: Running MET UnclusteredEn down variation"<<endl;
		met_err = 11;}
	else{cout << "INFO: Running nominal MET "<<endl;}
	
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  LOOP OVER EVENTS  /////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    for (Int_t iev=0;iev<nentries;iev++) {
        t->GetEntry(iev);
        if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
        
        //      int fill_number = run_to_fill.getFillNumber(ev.run);
        //      proton_reco.setAlignmentConstants(pots_align.getAlignmentConstants(fill_number));
        
        //assign randomly a run period
        TString period = lumi.assignRunPeriod();
        
        ////////////////////
        // EVENT WEIGHTS //
        ///////////////////
        float wgt(1.0);
        std::vector<double>plotwgts(1,wgt);
        
		////////////////////
		// EVENT cleaning //
		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2018_2017_data_and_MC_UL
		// list of met_filterBits defined in:
		// cat ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/python/miniAnalyzer_cfi.py
		////////////////////
		//bool passMETfilters = selector.passMETFilters(ev);
		bool passMETfilters = ( ((ev.met_filterBits >> 0) & 0x1) && //Flag_goodVertices
		                ((ev.met_filterBits >> 1) & 0x1) && //Flag_globalSuperTightHalo2016Filter
			            ((ev.met_filterBits >> 2) & 0x1) && //Flag_HBHENoiseFilter
                        ((ev.met_filterBits >> 3) & 0x1) && //Flag_HBHENoiseIsoFilter
                        ((ev.met_filterBits >> 4) & 0x1) && //Flag_EcalDeadCellTriggerPrimitiveFilter
                        ((ev.met_filterBits >> 5) & 0x1) && //Flag_eeBadScFilter
                        ((ev.met_filterBits >> 6) & 0x1) && //Flag_ecalBadCalibFilter          
                        ((ev.met_filterBits >> 7) & 0x1) && //Flag_BadPFMuonFilter          
                        ((ev.met_filterBits >> 8) & 0x1) ); //Flag_BadPFMuonDzFilter				
        //////////////////
        // CORRECTIONS //
        /////////////////
        btvSF.addBTagDecisions(ev);
        if(!isData) btvSF.updateBTagDecisions(ev, optionhf, optionlf); 
        //btvSF.updateBTagDecisions(ev); // check it !!!
        
        //////////////////////////
        // RECO LEVEL SELECTION //
        //////////////////////////
        TString chTag = selector.flagFinalState(ev, {}, {}, sys); // writes the name in chTag, last argument is JEC/JER systematics
		selector_nominal.flagFinalState(ev, {}, {}, 0); // selector tool with nominal event selection
        // ch
		Int_t ch_tag = 0;
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 3, plotwgts); // count all events before any selection
        if      (chTag=="EM")    ch_tag = 1;
        else if (chTag=="MM")    ch_tag = 2;
        else if (chTag=="EE")    ch_tag = 3;
        else if (chTag=="E")     ch_tag = 4;
        else if (chTag=="M")     ch_tag = 5;
        else                     ch_tag = 9;
        ht.fill("ch_tag", ch_tag, plotwgts);
#endif
        std::vector<Particle> &leptons     = selector.getSelLeptons();
		//std::vector<Particle> allPhotons=selector.getSelPhotons();
        std::vector<Jet>      &allJets        = selector.getJets();
        std::vector<Jet>      &nominalJets        = selector_nominal.getJets();
        std::vector<Jet>      jets,bJets,lightJets;
        std::vector<Particle> selectedLeptons;
        TLorentzVector extra_system(0,0,0,0);
        float extra_rapidity(0.), extra_rapidity_sum(0.);

        double p1_xi =0.; // proton in positive pot
        double p2_xi =0.; // proton in negative pot
        double p1_x = 0, p1_y =0.; // proton near track position in positive pot
	    double p2_x = 0, p2_y =0.; // proton near track position in positive pot
        double p1_220_x = 0, p1_220_y =0.; // proton near track position in positive pot
	    double p2_220_x = 0, p2_220_y =0.; // proton near track position in positive pot
        double mpp=-99., ypp =-99.;  //reconstructed pp quantities
        double mttbar=-99., yttbar =-99.;  //reconstructed ttbar quantities
 
        // high level variables for BDT
        double yvis = 0;
        double mll = 0;
        double min_dy = 0;


        //  selection of lightJets and bJets
        for(size_t ij=0; ij<allJets.size(); ij++) {
			
            double pt=allJets[ij].Pt();
            if(pt<30.) continue;

			int idx=allJets[ij].getJetIndex();		
			int jid=ev.j_id[idx];
			bool passLoosePu((jid>>2)&0x1);    
			if(!passLoosePu) {continue;}
			
			jets.push_back(allJets[ij]);
            if(allJets[ij].flavor()==5) bJets.push_back(allJets[ij]);
            else { 
                lightJets.push_back(allJets[ij]);
                extra_system+=allJets[ij];
                extra_rapidity_sum+=fabs(allJets[ij].Rapidity());
            }
        }
		
        if(lightJets.size()>0) extra_rapidity=extra_system.Rapidity();
        else extra_rapidity=0;
		
        // met selection:
		met_pt=ev.met_pt;
		met_phi=ev.met_phi;
        TLorentzVector met(0,0,0,0);
		
		// propogate JER to MET:
		TLorentzVector tmp_met(0,0,0,0);
		tmp_met.SetPtEtaPhiM(met_pt,0.,met_phi,0.);
		float newx=tmp_met.Px(), newy=tmp_met.Py();
		for(size_t ij=0; ij<allJets.size(); ij++) {
			int idx=allJets[ij].getJetIndex();
			int jid=ev.j_id[idx];
			bool passLoosePu((jid>>2)&0x1);    
			if(!passLoosePu) {continue;}
			newx+=allJets[ij].Px()*(1-ev.j_rawsf[idx]);
			newy+=allJets[ij].Py()*(1-ev.j_rawsf[idx]);
		}
		tmp_met.SetPxPyPzE(newx,newy,0,sqrt(newx*newx+newy*newy));
		met_pt=tmp_met.Pt();
		met_phi=tmp_met.Phi();
			
		if(sys){ // if apply JEC/JER variation, propogate the difference to MET
            tmp_met.SetPtEtaPhiM(met_pt,0.,met_phi,0.);
			float dx = 0, dy =0;
			for(size_t ij=0; ij<nominalJets.size(); ij++) {
				int idx=nominalJets[ij].getJetIndex();		
				int jid=ev.j_id[idx];
				bool passLoosePu((jid>>2)&0x1);    
				if(!passLoosePu) {continue;}
				dx+=nominalJets[ij].Px();
				dy+=nominalJets[ij].Py();
			}
			for(size_t ij=0; ij<allJets.size(); ij++) {
				int idx=allJets[ij].getJetIndex();		
				int jid=ev.j_id[idx];
				bool passLoosePu((jid>>2)&0x1);    
				if(!passLoosePu) {continue;}
				dx-=allJets[ij].Px();
				dy-=allJets[ij].Py();
			}
			float newx=tmp_met.Px()+dx, newy=tmp_met.Py()+dy;
			tmp_met.SetPxPyPzE(newx,newy,0,sqrt(newx*newx+newy*newy));
			met_pt=tmp_met.Pt();
			met_phi=tmp_met.Phi();
		}
		if(met_err){
			met_pt=ev.met_ptShifted[met_err];
			met_phi=ev.met_phiShifted[met_err];
		}
		met.SetPtEtaPhiM(met_pt,0.,met_phi,0.);	
        
        // selection of leptons
        for( size_t i_lept=0;i_lept<leptons.size();i_lept++) {
            if (leptons[i_lept].pt()<20.) continue;
			if (leptons[i_lept].id()==11 && fabs(leptons[i_lept].eta())>2.4) continue;
            //        if (leptons[i_lept].reliso()>0.10) continue;   //usually tighter
            selectedLeptons.push_back(leptons[i_lept]);
        }

        // selection of protons and reco unc.
		// Note fwdtrk_xiError is deprecated, using reco from 
		// https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsRecoCharacteristics
        for (int ift=0; ift<ev.nfwdtrk; ift++) {
            const unsigned short pot_raw_id = ev.fwdtrk_pot[ift];
            if(ev.fwdtrk_method[ift]==1){  // selecting only MultiRP protons
              if (pot_raw_id<100){ // positive z  (pot_raw_id=3)
                    p1_xi = ev.fwdtrk_xi[ift] + pps45_err*PPS_reco->getRecoErr(ev.fwdtrk_xi[ift],0,ev.run);
		    p1_x = ev.fwdtrk_NearX[ift];
		    p1_y = ev.fwdtrk_NearY[ift];
		    p1_220_x = ev.fwdtrk_FarX[ift];
		    p1_220_y = ev.fwdtrk_FarY[ift];

              }
              else {   // negative z   (pot_raw_id=103)
                p2_xi = ev.fwdtrk_xi[ift] + pps56_err*PPS_reco->getRecoErr(ev.fwdtrk_xi[ift],1,ev.run);
		p2_x = ev.fwdtrk_NearX[ift];
		p2_y = ev.fwdtrk_NearY[ift];
		p2_220_x = ev.fwdtrk_FarX[ift];
		p2_220_y = ev.fwdtrk_FarY[ift];
              }
            }
        }
	
        outVars["nJets"]=jets.size();
        outVars["nBjets"]=bJets.size();

        //compute mpp and ypp
        if(p1_xi>0 && p2_xi>0) {
            mpp=13000*sqrt(p1_xi*p2_xi);
            ypp=0.5*TMath::Log(p1_xi/p2_xi);
        }
		
        // Store proton pool at preselection
		if(ev.isData){
			m_protonVars_p1_xi = p1_xi;
			m_protonVars_p2_xi = p2_xi;
			outPT->Fill();
		}

    // ---- EVENT SELECTION --------------------------------------------------------------

	if(chTag!="EE" && chTag!="MM" && chTag!="EM")   continue; // events with 2 LEPTONS (electrons (id=11) or muons (id=13))
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 4, plotwgts); // count events after channel selection
#endif	
	//if(ev.isData && !passMETfilters)   continue; // event cleaning
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 5, plotwgts); // count events after channel selection
#endif
    if(selectedLeptons.size()!=2) continue; // ONLY events with =2 selected leptons
    if(selectedLeptons[0].Pt()<30 || fabs(leptons[0].Eta())>2.1) continue; //pt is Trigger safe
    if(selectedLeptons[0].charge()*selectedLeptons[1].charge()>0) continue; //exclude same sign
    if((selectedLeptons[0]+selectedLeptons[1]).M()<20) continue; // ONLY events with M(ll)>20
    if(fabs((selectedLeptons[0]+selectedLeptons[1]).M()-91)<15 && chTag!="EM") continue; // Zpeak cut

#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 6, plotwgts); // count events after selection on number of leptons and Z peak cut
#endif
        if ( jets.size()  < 2 )        continue; // ONLY events with at least 4 jets
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 7, plotwgts); // count events after selection on number of jets
#endif
        if ( bJets.size() < 2 )        continue; // ONLY events with at least 2 BJets
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 8, plotwgts); // count events after selection on number of Bjets
#endif

#ifdef HISTOGRAMS_ON
        ht.fill("puwgtctr",0,plotwgts);
#endif  
        outVars["weight"] = outVars["gen_wgt"] = outVars["toppt_wgt"] = 1;
		outVars["EL_SF_wgt"] = outVars["MU_SF_wgt"] = 1;
		outVars["EL_trigSF_wgt"] = outVars["MU_trigSF_wgt"] = 1;

        outVars["pu_wgt"] = outVars["ptag_wgt"] =  outVars["L1Prefire_wgt"] = outVars["ppsSF_wgt"] = 1;
        
		outVars["EL_SF_wgt_err"] = outVars["MU_SF_wgt_err"] = 0;
		outVars["EL_trigSF_wgt_err"] = outVars["MU_trigSF_wgt_err"] = 0;

		outVars["L1Prefire_wgt_err"] =  0;
		outVars["ppsSF_wgt_err"] = outVars["ptag_wgt_err"] = 0;
		
        outVars["ren_err"] = outVars["fac_err"] = 0;
	    
		outVars["pdf_as"] = 0;
		outVars["pdf_hs"] = 0;

		for(int ips=0;ips<NPSRad_weights;ips++){
			isr_Up[ips] = 0; fsr_Up[ips] = 0;
			isr_Down[ips] = 0; fsr_Down[ips] = 0;
		}		
		
        if (!ev.isData) {
            wgt  = (normH? normH->GetBinContent(1) : 1.0);          // norm weight
            double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);     // pu weight
            std::vector<double>puPlotWgts(1,puWgt);
            
#ifdef HISTOGRAMS_ON
            ht.fill("puwgtctr",1,puPlotWgts);
#endif
            
            leptons=selectedLeptons;

            // lepton trigger*selection weights (update the code later)
            EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons,{},{},"");
            EffCorrection_t  sel1SF = lepEffH.getOfflineCorrection(leptons[0], period);
            EffCorrection_t  sel2SF = lepEffH.getOfflineCorrection(leptons[1], period);
            
			if(leptons[0].id()==11){
				outVars["EL_trigSF_wgt"] *= trigSF.first;
				outVars["EL_trigSF_wgt_err"] *= trigSF.second;
				if(!outVars["EL_trigSF_wgt"]) cout << "WARNING: EL_trigSF_wgt = 0, check your selection! " << endl;
				outVars["EL_SF_wgt"] *= sel1SF.first;
				outVars["EL_SF_wgt_err"] *= sel1SF.second;
			}
			if(leptons[0].id()==13){
				outVars["MU_trigSF_wgt"] *= trigSF.first;
				outVars["MU_trigSF_wgt_err"] *= trigSF.second;
				if(!outVars["MU_trigSF_wgt"]) cout << "WARNING: MU_trigSF_wgt = 0, check your selection! " << endl;
				outVars["MU_SF_wgt"] *= sel1SF.first;
				outVars["MU_SF_wgt_err"] *= sel1SF.second;
			}
            if(leptons[1].id()==11){
                outVars["EL_trigSF_wgt"] *= trigSF.first;
                outVars["EL_trigSF_wgt_err"] *= trigSF.second;
                if(!outVars["EL_trigSF_wgt"]) cout << "WARNING: EL_trigSF_wgt = 0, check your selection! " << endl;
                outVars["EL_SF_wgt"] *= sel2SF.first;
                outVars["EL_SF_wgt_err"] *= sel2SF.second;
            }
            if(leptons[1].id()==13){
                outVars["MU_trigSF_wgt"] *= trigSF.first;
                outVars["MU_trigSF_wgt_err"] *= trigSF.second;
                if(!outVars["MU_trigSF_wgt"]) cout << "WARNING: MU_trigSF_wgt = 0, check your selection! " << endl;
                outVars["MU_SF_wgt"] *= sel2SF.first;
                outVars["MU_SF_wgt_err"] *= sel2SF.second;
            }
			
            wgt *= outVars["EL_SF_wgt"]*outVars["MU_SF_wgt"];
            wgt *= outVars["EL_trigSF_wgt"]*outVars["MU_trigSF_wgt"];

			//L1 pre-fire
			EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});
			outVars["L1Prefire_wgt"] = l1prefireProb.first;
			outVars["L1Prefire_wgt_err"] = l1prefireProb.second;

			wgt *= outVars["L1Prefire_wgt"];
            
            //top pt weighting
            outVars["toppt_wgt"] = 1.0;
            if(isTTbar) {
                for (int igen=0; igen<ev.ngtop; igen++) {
                    if(abs(ev.gtop_id[igen])!=6) continue;
                    outVars["toppt_wgt"] *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
                }
            }
			// The recommendation is not to apply this weight but use it as a systematic uncertainty.
			//wgt *= outVars["toppt_wgt"]; 
			
			// generator weights
			outVars["gen_wgt"] = (ev.g_nw>0 ? ev.g_w[0] : 1.0);
			
            wgt *= outVars["gen_wgt"];                              // generator level weights
					
            plotwgts[0]=wgt;                                        //update weight for plotter
			outVars["weight"] = wgt;
			
			// Systematic uncertainties (convert to 1 +/- form):
			if(ev.g_nw>5 && ev.g_w[1]){ // scale variations
				outVars["fac_err"] = (ev.g_w[2]-ev.g_w[3])/(ev.g_w[2]+ev.g_w[3]);
				outVars["ren_err"] =  (ev.g_w[4]-ev.g_w[5])/(ev.g_w[4]+ev.g_w[5]);
				// protect against broken weights
				if(abs(outVars["fac_err"])>1) outVars["fac_err"]=0;
				if(abs(outVars["ren_err"])>1) outVars["ren_err"]=0;
			}
			
			// pdf variations 
			float sig_mean = 0; float Nmem=0;
			for (int i=7;i<107;i++){
				if (ev.g_w[i]==0) continue;
				sig_mean+=ev.g_w[i];
				Nmem++;
			}
			if(Nmem){
			  sig_mean /=Nmem; // eq 22 in https://arxiv.org/pdf/1510.03865.pdf
			
			  float sig_var = 0;
			  for (int i=7;i<107;i++){
				if (ev.g_w[i]==0) continue;
				sig_var+=(ev.g_w[i]-sig_mean)*(ev.g_w[i]-sig_mean);
			  }			
			  outVars["pdf_hs"] = sqrt(sig_var/(Nmem-1))/sig_mean; // eq 21 in https://arxiv.org/pdf/1510.03865.pdf
			  
			  // for pdf_as use eq 27 in https://arxiv.org/pdf/1510.03865.pdf
			  // r=0.75 used to adjust to 68% for d(as)=0.002 as suggested in eq 29 
			  outVars["pdf_as"] = 0.75*0.5*(ev.g_w[108]-ev.g_w[107])/sig_mean; 
			}

			// Py8 PS variations 
			int nW = ev.g_npsw;
			float PSnom_weight = (nW>PSmap["nominal"]) ? ev.g_psw[PSmap["nominal"]] : 0;
			if(nW==46 && PSnom_weight){ // 46 weights of ISR and FSR, take only the first ones
				isr_Up[0] = (ev.g_psw[PSmap["isrDefLo"]])/PSnom_weight-1;
				fsr_Up[0] = (ev.g_psw[PSmap["fsrDefLo"]])/PSnom_weight-1;
				isr_Down[0] = 1-(ev.g_psw[PSmap["isrDefHi"]])/PSnom_weight;
				fsr_Down[0] = 1-(ev.g_psw[PSmap["fsrDefHi"]])/PSnom_weight;
					
				//G2GG muR
				isr_Up[1] = (ev.g_psw[PSmap["isr_G2GG_muR_up"]])/PSnom_weight-1;
				fsr_Up[1] = (ev.g_psw[PSmap["fsr_G2GG_muR_up"]])/PSnom_weight-1;
				isr_Down[1] = 1-(ev.g_psw[PSmap["isr_G2GG_muR_dn"]])/PSnom_weight;
				fsr_Down[1] = 1-(ev.g_psw[PSmap["fsr_G2GG_muR_dn"]])/PSnom_weight;

				//G2QQ muR
				isr_Up[2] = (ev.g_psw[PSmap["isr_G2QQ_muR_up"]])/PSnom_weight-1;
				fsr_Up[2] = (ev.g_psw[PSmap["fsr_G2QQ_muR_up"]])/PSnom_weight-1;
				isr_Down[2] = 1-(ev.g_psw[PSmap["isr_G2QQ_muR_dn"]])/PSnom_weight;
				fsr_Down[2] = 1-(ev.g_psw[PSmap["fsr_G2QQ_muR_dn"]])/PSnom_weight;

				//Q2QG muR
				isr_Up[3] = (ev.g_psw[PSmap["isr_Q2QG_muR_up"]])/PSnom_weight-1;
				fsr_Up[3] = (ev.g_psw[PSmap["fsr_Q2QG_muR_up"]])/PSnom_weight-1;
				isr_Down[3] = 1-(ev.g_psw[PSmap["isr_Q2QG_muR_dn"]])/PSnom_weight;
				fsr_Down[3] = 1-(ev.g_psw[PSmap["fsr_Q2QG_muR_dn"]])/PSnom_weight;

				//X2XG muR
				isr_Up[4] = (ev.g_psw[PSmap["isr_X2XG_muR_up"]])/PSnom_weight-1;
				fsr_Up[4] = (ev.g_psw[PSmap["fsr_X2XG_muR_up"]])/PSnom_weight-1;
				isr_Down[4] = 1-(ev.g_psw[PSmap["isr_X2XG_muR_dn"]])/PSnom_weight;
				fsr_Down[4] = 1-(ev.g_psw[PSmap["fsr_X2XG_muR_dn"]])/PSnom_weight;

				//G2GG cNS
				isr_Up[5] = (ev.g_psw[PSmap["isr_G2GG_cNS_up"]])/PSnom_weight-1;
				fsr_Up[5] = (ev.g_psw[PSmap["fsr_G2GG_cNS_up"]])/PSnom_weight-1;
				isr_Down[5] = 1-(ev.g_psw[PSmap["isr_G2GG_cNS_dn"]])/PSnom_weight;
				fsr_Down[5] = 1-(ev.g_psw[PSmap["fsr_G2GG_cNS_dn"]])/PSnom_weight;

				//G2QQ cNS
				isr_Up[6] = (ev.g_psw[PSmap["isr_G2QQ_cNS_up"]])/PSnom_weight-1;
				fsr_Up[6] = (ev.g_psw[PSmap["fsr_G2QQ_cNS_up"]])/PSnom_weight-1;
				isr_Down[6] = 1-(ev.g_psw[PSmap["isr_G2QQ_cNS_dn"]])/PSnom_weight;
				fsr_Down[6] = 1-(ev.g_psw[PSmap["fsr_G2QQ_cNS_dn"]])/PSnom_weight;

				//Q2QG cNS
				isr_Up[7] = (ev.g_psw[PSmap["isr_Q2QG_cNS_up"]])/PSnom_weight-1;
				fsr_Up[7] = (ev.g_psw[PSmap["fsr_Q2QG_cNS_up"]])/PSnom_weight-1;
				isr_Down[7] = 1-(ev.g_psw[PSmap["isr_Q2QG_cNS_dn"]])/PSnom_weight;
				fsr_Down[7] = 1-(ev.g_psw[PSmap["fsr_Q2QG_cNS_dn"]])/PSnom_weight;

				//X2XG cNS
				isr_Up[8] = (ev.g_psw[PSmap["isr_X2XG_cNS_up"]])/PSnom_weight-1;
				fsr_Up[8] = (ev.g_psw[PSmap["fsr_X2XG_cNS_up"]])/PSnom_weight-1;
				isr_Down[8] = 1-(ev.g_psw[PSmap["isr_X2XG_cNS_dn"]])/PSnom_weight;
				fsr_Down[8] = 1-(ev.g_psw[PSmap["fsr_X2XG_cNS_dn"]])/PSnom_weight;			
			}
			if(nW==24 && PSnom_weight){ // in signal no ISR weights

				fsr_Up[0] = (ev.g_psw[PSmap["fsrDefLo"]-4])/PSnom_weight-1;
				fsr_Down[0] = 1-(ev.g_psw[PSmap["fsrDefHi"]-3])/PSnom_weight;

				//G2GG muR
				fsr_Up[1] = (ev.g_psw[PSmap["fsr_G2GG_muR_up"]-6])/PSnom_weight-1;
				fsr_Down[1] = 1-(ev.g_psw[PSmap["fsr_G2GG_muR_dn"]-6])/PSnom_weight;

				//G2QQ muR
				fsr_Up[2] = (ev.g_psw[PSmap["fsr_G2QQ_muR_up"]-6])/PSnom_weight-1;
				fsr_Down[2] = 1-(ev.g_psw[PSmap["fsr_G2QQ_muR_dn"]-6])/PSnom_weight;

				//Q2QG muR
				fsr_Up[3] = (ev.g_psw[PSmap["fsr_Q2QG_muR_up"]-6])/PSnom_weight-1;
				fsr_Down[3] = 1-(ev.g_psw[PSmap["fsr_Q2QG_muR_dn"]-6])/PSnom_weight;

				//X2XG muR
				fsr_Up[4] = (ev.g_psw[PSmap["fsr_X2XG_muR_up"]-6])/PSnom_weight-1;
				fsr_Down[4] = 1-(ev.g_psw[PSmap["fsr_X2XG_muR_dn"]-6])/PSnom_weight;

				//G2GG cNS
				fsr_Up[5] = (ev.g_psw[PSmap["fsr_G2GG_cNS_up"]-6])/PSnom_weight-1;
				fsr_Down[5] = 1-(ev.g_psw[PSmap["fsr_G2GG_cNS_dn"]-6])/PSnom_weight;

				//G2QQ cNS
				fsr_Up[6] = (ev.g_psw[PSmap["fsr_G2QQ_cNS_up"]-6])/PSnom_weight-1;
				fsr_Down[6] = 1-(ev.g_psw[PSmap["fsr_G2QQ_cNS_dn"]-6])/PSnom_weight;

				//Q2QG cNS
				fsr_Up[7] = (ev.g_psw[PSmap["fsr_Q2QG_cNS_up"]-6])/PSnom_weight-1;
				fsr_Down[7] = 1-(ev.g_psw[PSmap["fsr_Q2QG_cNS_dn"]-6])/PSnom_weight;

				//X2XG cNS
				fsr_Up[8] = (ev.g_psw[PSmap["fsr_X2XG_cNS_up"]-6])/PSnom_weight-1;
				fsr_Down[8] = 1-(ev.g_psw[PSmap["fsr_X2XG_cNS_dn"]-6])/PSnom_weight;	
			}
		
        } // end is MC
        
        //if (ev.isData) {
        //    const edm::EventID ev_id( ev.run, ev.lumi, ev.event );
            // LHC information retrieval from LUT
        //    const ctpps::conditions_t lhc_cond = lhc_conds.get( ev_id );
        //    const double xangle = lhc_cond.crossing_angle;
        //    for (int ift=0; ift<ev.nfwdtrk; ift++) {
                // only look at strips!
        //        const unsigned short pot_raw_id = 100*ev.fwdtrk_arm[ift]+/*10*ev.fwdtrk_station[ift]+*/ev.fwdtrk_pot[ift];
        //        const ctpps::alignment_t align = ctpps_aligns.get( ev_id, pot_raw_id );
        //        double xi, xi_error;
        //        proton_reco.reconstruct(xangle, pot_raw_id, ev.fwdtrk_x[ift]/10.+align.x_align, xi, xi_error);
        //    }
        //}
      
        // ----- START RECONSTRUCTION OF TTBAR -------------------------------------------------------
        if(bJets.size()>=2)        {    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //  in theory this selection was made already
            do_kin_reco(leptons, jets, bJets, met, mttbar, yttbar);
        
            //high level variables for BDT
            yvis = (bJets[0]+bJets[1]+leptons[0]+leptons[1]).Rapidity();
            mll = (leptons[0]+leptons[1]).M();
            min_dy= min(fabs((bJets[0]+leptons[0]).Rapidity())-fabs((bJets[1]+leptons[1]).Rapidity()),  fabs((bJets[0]+leptons[1]).Rapidity())-fab\
s((bJets[1]+leptons[0]).Rapidity()));

#ifdef HISTOGRAMS_ON
            //control histograms
            ht.fill("nvtx",   ev.nvtx,      plotwgts);
            ht.fill("nbjets", bJets.size(), plotwgts);
            ht.fill("njets",  jets.size(),  plotwgts);

            ht.fill("mttbar_rec", mttbar,   plotwgts);
            ht.fill("yttbar_rec", yttbar,   plotwgts);

#endif

            //   --- write the variables which will be written in tree ---
            outVars["ttbar_y"]=yttbar;
            outVars["ttbar_m"]=mttbar;
            outVars["yvis"]=yvis;
            outVars["mll"]=mll;
            outVars["min_dy"]=min_dy;
            outVars["extra_rapidity"]=extra_rapidity;
            outVars["extra_rapidity_sum"]=extra_rapidity_sum;
            
            // quantities for all objects in event
            outVars["nLightJets"]=lightJets.size();
            outVars["cat"]=float(ch_tag);

            outVars["l1_pt"]=leptons[0].Pt();
            outVars["l1_eta"]=leptons[0].Rapidity();
            outVars["l1_phi"]=leptons[0].Phi();
            outVars["l1_m"]=leptons[0].M();
            outVars["l1_E"]=leptons[0].E();
            outVars["lepton1_isolation"]=ev.l_relIso[0];
            outVars["l2_pt"]=leptons[1].Pt();
            outVars["l2_eta"]=leptons[1].Rapidity();
            outVars["l2_phi"]=leptons[1].Phi();
            outVars["l2_m"]=leptons[1].M();
            outVars["l2_E"]=leptons[1].E();
            outVars["lepton2_isolation"]=ev.l_relIso[1];

            outVars["p1_xi"]=p1_xi;
            outVars["p2_xi"]=p2_xi;

            outVars["mpp"]=mpp;
            outVars["ypp"]=ypp;
            
            outVars["p1_x"]=p1_x;
            outVars["p2_x"]=p2_x;
            outVars["p1_y"]=p1_y;
            outVars["p2_y"]=p2_y;
            outVars["p1_220_x"]=p1_220_x;
            outVars["p2_220_x"]=p2_220_x;
            outVars["p1_220_y"]=p1_220_y;
            outVars["p2_220_y"]=p2_220_y;
			
            if (lightJets.size()==0) {
            outVars["lightJet0_pt"]=  0. ;
            outVars["lightJet0_eta"]= 0. ;
            outVars["lightJet0_phi"]= 0. ;
            outVars["lightJet0_m"]=   0. ;
            outVars["lightJet0_E"]=   0. ;
            outVars["lightJet1_pt"]=  0. ;
            outVars["lightJet1_eta"]= 0. ;
            outVars["lightJet1_phi"]= 0. ;
            outVars["lightJet1_m"]=   0. ;
            outVars["lightJet1_E"]=   0. ;
            }
            if (lightJets.size()==1) {
            outVars["lightJet0_pt"]=lightJets[1].Pt();
            outVars["lightJet0_eta"]=lightJets[1].Rapidity();
            outVars["lightJet0_phi"]=lightJets[1].Phi();
            outVars["lightJet0_m"]=lightJets[1].M();
            outVars["lightJet0_E"]=lightJets[1].E();
            outVars["lightJet1_pt"]=  0. ;
            outVars["lightJet1_eta"]= 0. ;
            outVars["lightJet1_phi"]= 0. ;
            outVars["lightJet1_m"]=   0. ;
            outVars["lightJet1_E"]=   0. ;
            }
            else if (lightJets.size()>=2) {
                outVars["lightJet0_pt"]=  lightJets[0].Pt();
                outVars["lightJet0_eta"]= lightJets[0].Rapidity();
                outVars["lightJet0_phi"]= lightJets[0].Phi();
                outVars["lightJet0_m"]=   lightJets[0].M();
                outVars["lightJet0_E"]=   lightJets[0].E();
                outVars["lightJet1_pt"]=  lightJets[1].Pt();
                outVars["lightJet1_eta"]= lightJets[1].Rapidity();
                outVars["lightJet1_phi"]= lightJets[1].Phi();
                outVars["lightJet1_m"]=   lightJets[1].M();
                outVars["lightJet1_E"]=   lightJets[1].E();
            }

            outVars["bJet0_pt"]=bJets[0].Pt();
            outVars["bJet0_eta"]=bJets[0].Rapidity();
            outVars["bJet0_phi"]=bJets[0].Phi();
            outVars["bJet0_m"]=bJets[0].M();
            outVars["bJet0_E"]=bJets[0].E();
            outVars["bJet1_pt"]=bJets[1].Pt();
            outVars["bJet1_eta"]=bJets[1].Rapidity();
            outVars["bJet1_phi"]=bJets[1].Phi();
            outVars["bJet1_m"]=bJets[1].M();
            outVars["bJet1_E"]=bJets[1].E();
            if (bJets.size()==3) {
                outVars["bJet2_pt"]=  bJets[2].Pt();
                outVars["bJet2_eta"]= bJets[2].Rapidity();
                outVars["bJet2_phi"]= bJets[2].Phi();
                outVars["bJet2_m"]=   bJets[2].M();
                outVars["bJet2_E"]=   bJets[2].E();
                outVars["bJet3_pt"]=  0. ;
                outVars["bJet3_eta"]= 0. ;
                outVars["bJet3_phi"]= 0. ;
                outVars["bJet3_m"]=   0. ;
                outVars["bJet3_E"]=   0. ;
            } else if (bJets.size()>=4) {
                outVars["bJet2_pt"]=  bJets[2].Pt();
                outVars["bJet2_eta"]= bJets[2].Rapidity();
                outVars["bJet2_phi"]= bJets[2].Phi();
                outVars["bJet2_m"]=   bJets[2].M();
                outVars["bJet2_E"]=   bJets[2].E();
                outVars["bJet3_pt"]=  bJets[3].Pt();
                outVars["bJet3_eta"]= bJets[3].Rapidity();
                outVars["bJet3_phi"]= bJets[3].Phi();
                outVars["bJet3_m"]=   bJets[3].M();
                outVars["bJet3_E"]=   bJets[3].E();
            } else { // bJets.size()==2
                outVars["bJet2_pt"]=  0. ;
                outVars["bJet2_eta"]= 0. ;
                outVars["bJet2_phi"]= 0. ;
                outVars["bJet2_m"]=   0. ;
                outVars["bJet2_E"]=   0. ;
                outVars["bJet3_pt"]=  0. ;
                outVars["bJet3_eta"]= 0. ;
                outVars["bJet3_phi"]= 0. ;
                outVars["bJet3_m"]=   0. ;
                outVars["bJet3_E"]=   0. ;
            }
          
            //---------------------------------
            // FILL TREE
            outT->Fill();

        } // end of if(bJets.size()>=2 && lightJets.size()>=2)
    } // end of loop over events
        std::cout << std::endl;
		std::cout << "saved " << outT->GetEntries() << " events " << std::endl;
        //close input file
        f->Close();

        //save histos to file
        fOut->cd();

    #ifdef HISTOGRAMS_ON
        for (auto& it : ht.getPlots())  {
            it.second->SetDirectory(fOut); it.second->Write();
        }
        for (auto& it : ht.get2dPlots())  {
            it.second->SetDirectory(fOut); it.second->Write();
        }
    #endif
        outT->Write();
        if(isData) outPT->Write();
        fOut->Close();
}  // end of RunExclusiveTop()

// --- THAT'S ALL FOLKS ---

