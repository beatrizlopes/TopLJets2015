#include <iostream>
#include <algorithm>

/*
Code was developed for exclusive top analysis by Enrico Robutti and Michael Pitt
Code mix pileup protons with MC simulation and it was used for 2017 data only (special strip treatment)

Usage:
mergePileupProtons /eos/cms/store/group/phys_top/TTbarCentralExclProd/ntuples/mc/excl_ttbar_semilep_QED_xa150_era2017_preTS2.root /eos/cms/store/group/phys_top/TTbarCentralExclProd/ntuples/data/
NOTE: for the events in /eos/cms/store/group/phys_top/TTbarCentralExclProd/ntuples/mc/ 
  are already after the proton mixing was applied (be aware of wrong normalization if running twice)
*/

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TMath.h>
#include "TSystem.h"

#include "TopLJets2015/TopAnalysis/interface/protonTrackRatios.h"
#include "TopLJets2015/TopAnalysis/interface/PPSEff.h"

#include <time.h>
using namespace std;

int main(int argc, char* argv[])
{
	
  // CMSSW settings  
  const char* CMSSW_BASE = getenv("CMSSW_BASE");
  TString data_path = Form("%s/src/TopLJets2015/TopAnalysis/data/era2017/", CMSSW_BASE);

  // Parameters
  int rndSeed = 1234567890; 
  int nEventsToMix = -1;
  int nMCEventsToSkip = 0;
  
  TRandom3 *rand_gen = new TRandom3(rndSeed);
  
  // Check arguments
  if (argc < 3 || argc > 5) {
    cout << "Two to four arguments required!" << endl
         << "Usage: " << argv[0] << " <inputMCFileName> <inputDataPath\
> [<nMCEventsToMix>=all] [<nMCEventsToSkip>=0]" << endl;
    cout << "Exampe: " << argv[0] << " /eos/cms/store/group/phys_top/TTbarCentralExclProd/ntuples/mc/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root \
/eos/cms/store/group/phys_top/TTbarCentralExclProd/ntuples/data/" << endl;
    return 1;
  }
  
  // Check input files
  if(gSystem->AccessPathName(argv[1])){
	  cout << "ERROR! Missing input file: " << argv[1] << endl;
	  return 0;
  }
  if(gSystem->AccessPathName(argv[2])){
	  cout << "ERROR! wrong path to data files: " << argv[2] << endl;
	  return 0;
  }  
  
  string inMCFileName = argv[1];
  string inPUFileName = argv[2];
  string outFileName = inMCFileName.substr(inMCFileName.find_last_of('/') + 1, inMCFileName.find_last_of('.') - inMCFileName.find_last_of('/') - 1) + "_enriched.root";
  bool isSignal = TString(inMCFileName.c_str()).Contains("excl_ttbar") || TString(inMCFileName.c_str()).Contains("ExclusiveTTbar");
  
  if (argc >= 4) {
    TString arg3 = TString(argv[3]);
    if (arg3.IsDec())
      nEventsToMix = arg3.Atoi();
    else
      cout << "Argument 3 must be an integer number. All input events will be mixed." << endl;
  }
  if (argc >= 5) { 
    TString arg4 = TString(argv[4]);
    if (arg4.IsDec())
      nMCEventsToSkip = arg4.Atoi();
    else
      cout << "Argument 4 must be an integer number. No pileup events will be s\
kipped." << endl;
  }
   
   // Proton efficiency class
   PPSEff *MultiRP_eff = new PPSEff(Form("%s/pixelEfficiencies_multiRP.root",data_path.Data()));
   //PPSEff *Strip_eff = new PPSEff(Form("%s/PreliminaryEfficiencies_July132020_1D2DMultiTrack.root",data_path.Data()));
   //PPSEff *Strip_eff = new PPSEff(Form("%s/Alt2017ERestrictedRuns_PreliminaryEfficiencies_March122021_1D2DMultiTrack.root",data_path.Data()));
   PPSEff *Strip_eff = new PPSEff(Form("%s/PreliminaryEfficiencies_March302021_1D2DMultiTrack.root",data_path.Data()));
   
   // ---------------------------------------------------------------------------------------------------------------------------------- //
   // ---------------------------------------------------------------------------------------------------------------------------------- //
   // ---------------------------------------------------------------------------------------------------------------------------------- //
   // NEW code, create pools from the data files, calculate relative cross-sections and nvtx distributions (era,xangle)
   // Read data, prerare PU trees:
   float era_lumi[] = {2.367138592,8.680676647,4.142765586,1.447798525,13.247272608};
   TString era[] = {"2017B","2017C","2017D","2017E","2017F"}; int n_era = (sizeof(era)/sizeof(TString));
   int xangle[] = {120, 130, 140, 150}; int n_xa = (sizeof(xangle)/sizeof(int));
   
   const int n_PUregionsMAX = n_era * n_xa; // used in signal normalization
   
   // Extra normalization factor used for signal since we have dedicated signal samples per era/xangle
   float signal_fraction_regions[n_PUregionsMAX], extra_signal_normalization=0;
   if(isSignal){
	   
	  //Get extra normalization
      for(int i_era=0;i_era<n_era;i_era++){
	    TString name_el = Form("%s/SingleElectron_%s.root",inPUFileName.c_str(),era[i_era].Data());
	    TString name_mu = Form("%s/SingleMuon_%s.root",inPUFileName.c_str(),era[i_era].Data());
	      for(int i_xa=0;i_xa<n_xa;i_xa++){ // create proton pools
		    // calculate fraction of Xangle at preselection fpr ALL REGIONS:
		    TChain * _ch2 = new TChain("tree"); _ch2->Add(name_el); _ch2->Add(name_mu);
		    signal_fraction_regions[i_era*n_xa+i_xa] = _ch2->GetEntries(Form("beamXangle==%d",xangle[i_xa]));
	      }
	    // normalize properly per selected crossing-angle (sometimes data contains unselected values lke 141,142...)
	    float total_event_per_era = 0; for(int ii=0;ii<n_xa;ii++) total_event_per_era+=signal_fraction_regions[i_era*n_xa+ii];
	    for(int ii=0;ii<n_xa;ii++) signal_fraction_regions[i_era*n_xa+ii] *= (era_lumi[i_era]/29.885651958)/total_event_per_era;
      }

	  // reduce the number of regions in case if the sample is a simulated signal
	  
	  // set crossing angle
	  int i_xa = 0;
	  if(TString(inMCFileName.c_str()).Contains("xa120")) {xangle[0] = 120; i_xa=0;}
	  if(TString(inMCFileName.c_str()).Contains("xa130")) {xangle[0] = 130; i_xa=1;}
	  if(TString(inMCFileName.c_str()).Contains("xa140")) {xangle[0] = 140; i_xa=2;}
	  if(TString(inMCFileName.c_str()).Contains("xa150")) {xangle[0] = 150; i_xa=3;}  

	  // set list of periods
	  if(TString(inMCFileName.c_str()).Contains("postTS2")){
		   n_era = 2;
		   era_lumi[0] = era_lumi[3]; era[0] = "2017E";
		   era_lumi[1] = era_lumi[4]; era[1] = "2017F";
		   for(int i_era=3;i_era<5;i_era++)
		     extra_signal_normalization += signal_fraction_regions[i_era*n_xa+i_xa];
	  }
	  else {
		  n_era = 3; // for preTS2 first 3 entries are the same
		  for(int i_era=0;i_era<3;i_era++)
		    extra_signal_normalization += signal_fraction_regions[i_era*n_xa+i_xa];
	  }
	  
	  // set to 1 the crossing angle counter
	  n_xa=1; 
	  cout << "Signal extra normalization = " << extra_signal_normalization << endl;
   }

   const int n_PUregions = n_era * n_xa;
   // puleup histograms (used for reweighting)
   TH1F * pu_weights[n_PUregions];
   TFile * _puweights = new TFile("pu_weights.root","recreate"); // dummy file for protons pools to be accosiated with.

   cout << "Reads data from " <<    inPUFileName << endl;
   TTree * PUpr[n_PUregions];  // puleup proton trees (for each era and xangle)
   
   // calculate fraction of the total luminosity
   float total_lumi=0;
   for(int i_era=0;i_era<n_era;i_era++) total_lumi+=era_lumi[i_era];
   
   float norm_weight[n_PUregions], norm_weight_err[n_PUregions], fraction_regions[n_PUregions]; int counter_regions[n_PUregions];
   float norm_weight_0p[n_PUregions], norm_weight_0p_err[n_PUregions]; // for signal events
   float norm_weight_1pRP0[n_PUregions], norm_weight_1pRP0_err[n_PUregions]; // for signal events
   float norm_weight_1pRP1[n_PUregions], norm_weight_1pRP1_err[n_PUregions]; // for signal events
   for(int i_era=0;i_era<n_era;i_era++){
	   printf ("\rProcessing %s ", era[i_era].Data());
	   TString name_el = Form("%s/SingleElectron_%s.root",inPUFileName.c_str(),era[i_era].Data());
	   TString name_mu = Form("%s/SingleMuon_%s.root",inPUFileName.c_str(),era[i_era].Data());
	   
	   // Read 2 proton hits and total number of events
	   int n_p2_e = int(((TH1F*)TFile::Open(name_el)->Get("evt_count"))->GetBinContent(4));
	   int n_p2_m = int(((TH1F*)TFile::Open(name_mu)->Get("evt_count"))->GetBinContent(4));
	   int n_e = int(((TH1F*)TFile::Open(name_el)->Get("evt_count"))->GetBinContent(3));
	   int n_m = int(((TH1F*)TFile::Open(name_mu)->Get("evt_count"))->GetBinContent(3));
	   int n_p2 = (n_p2_e+n_p2_m), n = (n_e+n_m);
	   
	   // Get systematics for n_proton==2 event fraction from sub-selection of nBjet>0 events:
	   int n_e_sys = int(((TH1F*)TFile::Open(name_el)->Get("pn_count"))->GetBinContent(5));
	   int n_m_sys = int(((TH1F*)TFile::Open(name_mu)->Get("pn_count"))->GetBinContent(5));
	   TChain * _ch_protons = new TChain("protons");
	   _ch_protons->Add(name_el);
	   _ch_protons->Add(name_mu);
	   int n_sys = (n_e_sys+n_m_sys), n_p2_sys = _ch_protons->GetEntries("nBjets>0");
	   
	   // Get 1 proton hits (depend on the arm)
	   int n_p1_eRP0 = int(((TH1F*)TFile::Open(name_el)->Get("pn_count"))->GetBinContent(2));
	   int n_p1_mRP0 = int(((TH1F*)TFile::Open(name_mu)->Get("pn_count"))->GetBinContent(2));
	   int n_p1_eRP1 = int(((TH1F*)TFile::Open(name_el)->Get("pn_count"))->GetBinContent(3));
	   int n_p1_mRP1 = int(((TH1F*)TFile::Open(name_mu)->Get("pn_count"))->GetBinContent(3));
	   int n_p1_RP0 = (n_p1_eRP0+n_p1_mRP0), n_p1_RP1 = (n_p1_eRP1+n_p1_mRP1);
		
	   // events with exactly 0 pu tracks
	   int n_p0 = n - n_p1_RP0 - n_p1_RP1 - n_p2; 

	   for(int i_xa=0;i_xa<n_xa;i_xa++){ // create proton pools
		   TChain * _ch = new TChain("protons");
		   _ch->Add(name_el);
		   _ch->Add(name_mu);
  		   _puweights->cd();
		   PUpr[i_era*n_xa+i_xa] = (_ch->CopyTree(Form("beamXangle==%d",xangle[i_xa])))->CloneTree();
		   pu_weights[i_era*n_xa+i_xa] = new TH1F(Form("data_puw_%d",i_era*n_xa+i_xa),";nvtx;w",100,0,100);
		   PUpr[i_era*n_xa+i_xa]->Draw(Form("nvtx>>data_puw_%d",i_era*n_xa+i_xa));
		   counter_regions[i_era*n_xa+i_xa] = PUpr[i_era*n_xa+i_xa]->GetEntries();
		   pu_weights[i_era*n_xa+i_xa]->Scale(1/float(counter_regions[i_era*n_xa+i_xa]));
		   // calculate fraction of Xangle at preselection:
		   TChain * _ch2 = new TChain("tree"); _ch2->Add(name_el); _ch2->Add(name_mu);
		   fraction_regions[i_era*n_xa+i_xa] = _ch2->GetEntries(Form("beamXangle==%d",xangle[i_xa]));
	       norm_weight[i_era*n_xa+i_xa] = n_p2/float(n); // probability of 2 tracks
	       norm_weight_err[i_era*n_xa+i_xa] = (n_p2_sys/float(n_sys)) / norm_weight[i_era*n_xa+i_xa];
		   		   
		   // probabilities for 1 track in signal events
	       norm_weight_1pRP0[i_era*n_xa+i_xa] = n_p1_RP0/float(n); // probability of 0 tracks in RP0
	       norm_weight_1pRP0_err[i_era*n_xa+i_xa] = 0.95; // 5% flat

	       norm_weight_1pRP1[i_era*n_xa+i_xa] = n_p1_RP1/float(n); // probability of 0 tracks in RP1
	       norm_weight_1pRP1_err[i_era*n_xa+i_xa] = 0.95; // 5% flat

		   // probabilities for 0 track in signal events
	       norm_weight_0p[i_era*n_xa+i_xa] = (n_p0)/float(n); // probability of 0 tracks in both arms
	       norm_weight_0p_err[i_era*n_xa+i_xa] = 0.95; // 5% flat
		   
	   }
	    // normalize properly per selected crossing-angle (sometimes data contains unselected values like 100,110,...)
	   float total_event_per_era = 0; for(int ii=0;ii<n_xa;ii++) total_event_per_era+=fraction_regions[i_era*n_xa+ii];
	   for(int ii=0;ii<n_xa;ii++) fraction_regions[i_era*n_xa+ii] *= (era_lumi[i_era]/total_lumi)/total_event_per_era;
   }
   
   // List obtained fractions:
   cout << "INFO List obtained fractions for all "<<n_PUregions<<" SRs, (for 120,130,140,150), including statistics:"<<endl;
   for(int i_era=0;i_era<n_era;i_era++){
	   cout << "era " << era[i_era]<<": ";
	   for(int ii=0;ii<n_xa-1;ii++) cout << fraction_regions[i_era*n_xa+ii]<<" ("<<counter_regions[i_era*n_xa+ii]<<"),";
	   cout << fraction_regions[i_era*n_xa+n_xa-1] <<" ("<<counter_regions[i_era*n_xa+n_xa-1]<<")"<< endl;
   }   
   
   // variables to be used to read from PU pools
    unsigned int poll_run[n_PUregions]; float poll_p1_xi[n_PUregions],poll_p2_xi[n_PUregions];
	_puweights->cd();
	for(int i=0;i<n_PUregions;i++){
		PUpr[i]->SetBranchAddress("run",&poll_run[i]);
		PUpr[i]->SetBranchAddress("p1_xi",&poll_p1_xi[i]);
		PUpr[i]->SetBranchAddress("p2_xi",&poll_p2_xi[i]);
	}
    cout << "done "<< endl;
   //
   // ---------------------------------------------------------------------------------------------------------------------------------- //
   // ---------------------------------------------------------------------------------------------------------------------------------- //
   // ---------------------------------------------------------------------------------------------------------------------------------- //
   
   //Save relative fractions of periods:
   //_puweights->cd();
   //cout << "write " <<   _puweights->GetName()  << endl;
   //for(int i=0;i<n_PUregions;i++) pu_weights[i]->Write();
   //_puweights->Write();
   //_puweights->Close();
   
   
   
  // read MC file    
  TFile *oldfile = new TFile(inMCFileName.c_str());
  TTree * chMCEvents = (TTree *)oldfile->Get("tree");
  cout << "Store MC nvtx distribution " << endl;
  TH1F * mc_pu = new TH1F("mc_pu",";nvtx;w",100,0,100);
  chMCEvents->Draw("nvtx>>mc_pu","","norm");
  TH1F * evt_count = (TH1F*)oldfile->Get("evt_count"); // event counter with SumWeights
  
  //list of branches to update:
  unsigned int run; int nvtx;
  float beamXangle, p1_xi, p2_xi, weight, ppsSF_wgt, ppsSF_wgt_err, pu_wgt, ptag_wgt, ptag_wgt_err; 
  float p1_x = 0, p1_y = 0, p2_x = 0, p2_y = 0;  
  //float p1_220_x = 0, p1_220_y = 0, p2_220_x = 0, p2_220_y = 0;  
  chMCEvents->SetBranchAddress("run",&run);
  chMCEvents->SetBranchAddress("beamXangle",&beamXangle);
  chMCEvents->SetBranchAddress("p1_xi",&p1_xi);
  chMCEvents->SetBranchAddress("p2_xi",&p2_xi);
  if(isSignal){
	  chMCEvents->SetBranchAddress("p1_x",&p1_x);
	  chMCEvents->SetBranchAddress("p1_y",&p1_y);
	  chMCEvents->SetBranchAddress("p2_x",&p2_x);
	  chMCEvents->SetBranchAddress("p2_y",&p2_y);
	  //chMCEvents->SetBranchAddress("p1_220_x",&p1_220_x);
	  //chMCEvents->SetBranchAddress("p1_220_y",&p1_220_y);
	  //chMCEvents->SetBranchAddress("p2_220_x",&p2_220_x);
	  //chMCEvents->SetBranchAddress("p2_220_y",&p2_220_y);
  }
  chMCEvents->SetBranchAddress("weight",&weight);
  chMCEvents->SetBranchAddress("pu_wgt",&pu_wgt);
  chMCEvents->SetBranchAddress("ppsSF_wgt",&ppsSF_wgt);
  chMCEvents->SetBranchAddress("ppsSF_wgt_err",&ppsSF_wgt_err);
  chMCEvents->SetBranchAddress("ptag_wgt",&ptag_wgt);
  chMCEvents->SetBranchAddress("ptag_wgt_err",&ptag_wgt_err);
  chMCEvents->SetBranchAddress("nvtx",&nvtx);
  chMCEvents->SetBranchStatus("*",1); // activate all branches to copy
  
  int nMCEntries = chMCEvents->GetEntries();
  cout << nMCEntries << " events read from file(s) " << inMCFileName << endl;
  if (nEventsToMix == -1)
    nEventsToMix = nMCEntries;

  // Get pileup proton ratios
  protonTrackRatios ptr;
  int nLines = ptr.readFromFile(Form("%s/protonRatios_2017.dat",data_path.Data()));
  if (nLines <= 0) {
    cout << "No true-zero-track ratio read from file! No mixing performed." << endl;
    return 4;
  }
  
  // Open output mixed file
  TFile* fOutMixed = new TFile(outFileName.c_str(), "RECREATE");
  if(!fOutMixed->IsOpen()) {
    cout<< "Error opening output file mixedPUProtons. Abort." << endl;
    return 5;
  }
 
  // Create new output with updated PPS values:
  TTree* tMCMixed = chMCEvents->CloneTree(0);
  
  int signal_protons = 0; // additional variables used for signal to indecate if all/part/none of the protons are from pileup
  tMCMixed->Branch("signal_protons",&signal_protons);
   
  int times,timed;
  times=time(NULL);  
  cout << "Loop over MC entries and add pileup protons" << endl;
  // Loop over MC entries and mix pileup protons
  int iMCEntry;
  for (iMCEntry = 0; (iMCEntry < nEventsToMix) && (iMCEntry + nMCEventsToSkip < nMCEntries); iMCEntry++) {
    if(iMCEntry%1000==0) printf ("\r [%3.0f%%] done", 100.*(float)iMCEntry/(float)nEventsToMix);
    
	// ------------- old approach --------------------- //
	// sample region index using PDF and get protons from the pool
	//int i_reg=0; float cdf_=0, rndm = rand_gen->Rndm();
	//while (cdf_<rndm) {cdf_+= fraction_regions[i_reg];i_reg++;}
	//i_reg--; // region index should start from 0

	// ------------- new approach --------------------- //
	// samples region index using flat PDF, and asign extra weight to event
    int i_reg = rand_gen->Rndm() * n_PUregions;
	weight *= (fraction_regions[i_reg]*float(n_PUregions));

	// -------------------------------------------------- //
	
	int i_event = rand_gen->Rndm()*counter_regions[i_reg];  
	PUpr[i_reg]->GetEntry(i_event);
	
	// Add additional selection of proton pool here:
	// while(condition) {
	//	i_event = rand_gen->Rndm()*counter_regions[i_reg];
	//	PUpr[i_reg]->GetEntry(i_event);
	//}
	
	// asign the protons to the MC event
    chMCEvents->GetEntry(iMCEntry + nMCEventsToSkip);
	
	// fix run number and crossing-angle from proton pool
	run = poll_run[i_reg];
	beamXangle = xangle[i_reg % n_xa];
	
	// proton efficiency implementation (xi of protons that fail reco. will be set to zero)
	ppsSF_wgt = 1.; ppsSF_wgt_err=0;
	if(p1_xi>0){
		float SF = Strip_eff->getEff(p1_x,p1_y,0,run) * MultiRP_eff->getEff(p1_x,p1_y,0,run);
		if(rand_gen->Rndm()>SF) p1_xi=0;
		else ppsSF_wgt *= SF;
		ppsSF_wgt_err += MultiRP_eff->getRelEffErrSq(p1_x,p1_y,0,run);
	}
	if(p2_xi>0){
		float SF = Strip_eff->getEff(p2_x,p2_y,1,run) * MultiRP_eff->getEff(p2_x,p2_y,1,run);
		if(rand_gen->Rndm()>SF) p2_xi=0;
		else ppsSF_wgt *= SF;
		ppsSF_wgt_err += MultiRP_eff->getRelEffErrSq(p2_x,p2_y,1,run);
	}
	ppsSF_wgt_err = sqrt(ppsSF_wgt_err);
	
	// mixing with the pileup protons (note the different treatment for the signal)
	if(isSignal && p1_xi>0 && p2_xi>0){ // two accepted signal protons
		ptag_wgt = norm_weight_0p[i_reg];
		ptag_wgt *= ptr.trueZeroTracksRatio(run, beamXangle, 0) * ptr.trueZeroTracksRatio(run, beamXangle, 1);
		ptag_wgt_err = norm_weight_0p_err[i_reg];
		weight *= ptag_wgt;
		//ppsSF_wgt = MultiRP_eff->getEff(p1_x,p1_y,0,run) *  MultiRP_eff->getEff(p2_x,p2_y,1,run);
		//ppsSF_wgt_err = MultiRP_eff->getRelEffErrSq(p1_x,p1_y,0,run) + 
		//                MultiRP_eff->getRelEffErrSq(p2_x,p2_y,1,run);
		//ppsSF_wgt_err = sqrt(ppsSF_wgt_err);
		// strip efficiency
		//ppsSF_wgt *= Strip_eff->getEff(p1_x,p1_y,0,run) * Strip_eff->getEff(p2_x,p2_y,1,run);
		signal_protons = 11;
		//weight *= ppsSF_wgt;
	}
	else if(isSignal && p1_xi==0 && p2_xi>0){ // one signal proton in arm1
		p1_xi = poll_p1_xi[i_reg];
		ptag_wgt = norm_weight_1pRP0[i_reg];
		ptag_wgt *= ptr.trueZeroTracksRatio(run, beamXangle, 1);
		ptag_wgt_err = norm_weight_1pRP0_err[i_reg];
		weight *= ptag_wgt;
		//ppsSF_wgt = MultiRP_eff->getEff(p2_x,p2_y,1,run);
		//ppsSF_wgt_err = MultiRP_eff->getRelEffErrSq(p2_x,p2_y,1,run);
		// strip efficiency
		//ppsSF_wgt *= Strip_eff->getEff(p2_x,p2_y,1,run);
		signal_protons = 1;
		//weight *= ppsSF_wgt;
	}
	else if(isSignal && p1_xi>0 && p2_xi==0){ // one signal proton in arm0
		p2_xi = poll_p2_xi[i_reg];
		ptag_wgt = norm_weight_1pRP1[i_reg];
		ptag_wgt *= ptr.trueZeroTracksRatio(run, beamXangle, 0);
		ptag_wgt_err = norm_weight_1pRP1_err[i_reg];
		weight *= ptag_wgt;
		//ppsSF_wgt = MultiRP_eff->getEff(p1_x,p1_y,0,run);
		//ppsSF_wgt_err = MultiRP_eff->getRelEffErrSq(p1_x,p1_y,0,run);
		// strip efficiency
		//ppsSF_wgt *= Strip_eff->getEff(p1_x,p1_y,0,run);
		signal_protons = 10;
		//weight *= ppsSF_wgt;
	}
	else{ // no signal protons in the acceptance region
		p1_xi = poll_p1_xi[i_reg];
		p2_xi = poll_p2_xi[i_reg];
		ptag_wgt = norm_weight[i_reg];
		ptag_wgt_err = norm_weight_err[i_reg];
		weight *= ptag_wgt;
		signal_protons = 0;
	}

	// Fix ptag weight from w_sys/w_nom to 1 +/- err
	ptag_wgt_err = TMath::Abs(1 - ptag_wgt_err);

    // puleup reweighting
	float w_mc = mc_pu->GetBinContent(nvtx+1);
	pu_wgt = (w_mc) ? pu_weights[i_reg]->GetBinContent(nvtx+1)/w_mc : 0;
	weight *= pu_wgt;
	
	// Add extra weight to signal since we have simulation for each era/xangle
	if(isSignal) {
		weight *= extra_signal_normalization;
	}
	   
	// add extra weight in case of signal events:
    //pup_2tkRatio = ptr.twoTracksRatio(run, beamXangle);
    //for (int iArm = 0; iArm < 2; iArm++) {
      //pup_0tkRatio[iArm] = ptr.zeroTracksRatio(run, beamXangle, iArm);
      //pup_true0tkRatio[iArm] = ptr.trueZeroTracksRatio(run, beamXangle, iArm);
    //}
    tMCMixed->Fill();
  }

  cout << "\nDone. " << iMCEntry << " events processed." << endl;
      timed=time(NULL);
    times=timed-times;
    cout << "time from start to end = " << times << endl;
	
  // Write mixed tree to file and close everything
  fOutMixed->cd();
  cout << "Writes " << fOutMixed->GetName() << endl;
  tMCMixed->Write();
  if(nMCEventsToSkip==0) evt_count->Write();
  fOutMixed->Close();
  
  return iMCEntry;

}
