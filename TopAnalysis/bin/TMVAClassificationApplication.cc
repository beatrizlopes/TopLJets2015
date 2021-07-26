/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained classifiers
/// within an analysis module
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Exectuable: TMVAClassificationApplication
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker
/// Developed further by Matteo Pisano

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;
using namespace std;

void TMVAClassificationApplication( char * _fname, TString myMethodList = "")
{
   
   TString fname(_fname);
   
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Cut optimisation
   //
   // Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost


   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST correspond in name and type to those given in the weight file(s) used
   Float_t var1, var2;
   Float_t var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, var15;
   reader->AddVariable( "mpp", &var1 );
   reader->AddVariable( "kinreco_mtt", &var2 );
   reader->AddVariable( "ypp",                &var3 );
   reader->AddVariable( "kinreco_ytt",                &var4 );
   reader->AddVariable( "yvis", &var5 );
   reader->AddVariable( "deltarll", &var6 );
   reader->AddVariable( "mll",                &var7 );
   reader->AddVariable( "fabs(b2phi-b1phi)", &var8 );
   reader->AddVariable( "ysum", &var9 );
   reader->AddVariable( "min_dy", &var10 );
   reader->AddVariable( "nljets", &var11 );
   reader->AddVariable( "metpt", &var12 );
   reader->AddVariable( "E2", &var13 );
   reader->AddVariable( "extra_ysum", &var14 );
   reader->AddVariable( "extrasystem_y", &var15 );

   // Spectator variables declared in the training have to be added to the reader, too
   //Float_t spec1,spec2;
   //reader->AddSpectator( "spec1 := var1*2",   &spec1 );
   //reader->AddSpectator( "spec2 := var1*3",   &spec2 );

   //Float_t Category_cat1, Category_cat2, Category_cat3;
   //if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
   //   reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
   //   reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
   //   reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   //}

   // Book the MVA methods
   //const char* CMSSW_BASE = getenv("CMSSW_BASE");
   //TString dir    = Form("%s/src/TopLJets2015/TopAnalysis/data/weights/ExclusiveTop/", CMSSW_BASE);
   //TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
		  std::cout << it->second << " , " << TString((it->first).c_str()) << std::endl;
         TString methodName = TString((it->first).c_str()) + TString(" method");
         TString weightfile = "/eos/user/b/bribeiro/exclusiveTTBarDilepton/BDTTrainingResults/TMVAClassification_BDT_new.weights.xml";
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 10000;

   TH1F *histBdt(0);

   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -1.0, 0.5 );

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input(0);
   //TString fname = "../NewRoberto/output_xi_errors_data.root";
   //TString fname = "output_xi_errors_sgn.root";
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   else {
      TFile::SetCacheFileDir(".");
      input = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD"); // if not: download from ROOT server
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TFile input_file(fname);
   TTree* theTree = (TTree*)input_file.Get("tree");
   TH1F * evt_count = (TH1F*)input_file.Get("evt_count");
   
   Float_t userVar_nljets, userVar_ytt, userVar_mtt, userVar_yvis, userVar_mll,
      userVar_mindy, userVar_extray, userVar_extraysum,
      userVar_b1phi, userVar_b2phi, userVar_met,
      userVar_b1y, userVar_b2y, userVar_b1E, userVar_b2E,
      userVar_l1y, userVar_l2y, userVar_l1E, userVar_l2E,
      userVar_l1eta, userVar_l2eta,
      userVar_l1phi, userVar_l2phi, userVar_xi0, userVar_xi1;


   std::vector<Float_t> vecVar(9); // vector for EvaluateMVA tests
   // Create output file with the tree:
   TFile* fOutMVA = new TFile(fname.ReplaceAll(".root","_MVA.root"), "RECREATE");
   if(!fOutMVA->IsOpen()) {
    std::cout<< "Error opening output file mixedPUProtons. Abort." << std::endl;
    return;
   }   
   // Create new output with MVA outputs:
   TTree* theTreeNew = theTree->CloneTree(0);
   float wBDT;
   theTreeNew->Branch("wBDT",&wBDT);
   
   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (int ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

      userVar_nljets=theTree->GetLeaf("nLightJets")->GetValue(0);
      userVar_ytt=theTree->GetLeaf("ttbar_y")->GetValue(0);
      userVar_mtt=theTree->GetLeaf("ttbar_m")->GetValue(0);
      userVar_yvis=theTree->GetLeaf("yvis")->GetValue(0);
      userVar_mll=theTree->GetLeaf("mll")->GetValue(0);
      userVar_mindy=theTree->GetLeaf("min_dy")->GetValue(0);
      userVar_extray=theTree->GetLeaf("extra_rapidity")->GetValue(0);
      userVar_extraysum=theTree->GetLeaf("extra_rapidity_sum")->GetValue(0);
      userVar_b1phi=theTree->GetLeaf("bJet0_phi")->GetValue(0);
      userVar_b2phi=theTree->GetLeaf("bJet1_phi")->GetValue(0);
      userVar_met=theTree->GetLeaf("met_pt")->GetValue(0);
      userVar_b1y=theTree->GetLeaf("bJet0_y")->GetValue(0);
      userVar_b2y=theTree->GetLeaf("bJet1_y")->GetValue(0);
      userVar_b1E=theTree->GetLeaf("bJet0_E")->GetValue(0);
      userVar_b2E=theTree->GetLeaf("bJet1_E")->GetValue(0);
 
      userVar_l1y=theTree->GetLeaf("l1_y")->GetValue(0);
      userVar_l2y=theTree->GetLeaf("l2_y")->GetValue(0);
      userVar_l1E=theTree->GetLeaf("l1_E")->GetValue(0);
      userVar_l2E=theTree->GetLeaf("l2_E")->GetValue(0);
      userVar_l1eta=theTree->GetLeaf("l1_eta")->GetValue(0);
      userVar_l2eta=theTree->GetLeaf("l2_eta")->GetValue(0);
      userVar_l1phi=theTree->GetLeaf("l1_phi")->GetValue(0);
      userVar_l2phi=theTree->GetLeaf("l2_phi")->GetValue(0);  

      userVar_xi0=theTree->GetLeaf("p1_xi")->GetValue(0);
      userVar_xi1=theTree->GetLeaf("p2_xi")->GetValue(0);

      //mpp
      var1 = 13000*sqrt(userVar_xi0*userVar_xi1);
      //mtt
      var2 = userVar_mtt;
      //ypp
      var3 = 0.5*TMath::Log(userVar_xi0/userVar_xi1);
      //ytt
      var4 = userVar_ytt;
      //yvis
      var5 = userVar_yvis;
      //drll
      var6 = sqrt( pow(userVar_l1eta-userVar_l2eta,2) + pow(userVar_l1phi-userVar_l2phi,2) );
      //mll
      var7 = userVar_mll;
      //deltaphi(bb)
      var8 = fabs(userVar_b1phi-userVar_b2phi);
      //ysum 
      var9 = fabs(userVar_b1y)+fabs(userVar_b2y)+fabs(userVar_l1y)+fabs(userVar_l2y);
      //min_dy
      var10 = userVar_mindy;
      //nljets
      var11 = userVar_nljets;
      //metpt
      var12 = userVar_met;
      //E2 
      var13 = pow(userVar_b1E+userVar_b2E+userVar_l1E+userVar_l2E+userVar_met,2);
      //extra ysum
      var14 = userVar_extraysum;
      //extray
      var15 = userVar_extray;  
      // Return the MVA outputs and fill into histograms

      if (Use["BDT"          ]){ wBDT = reader->EvaluateMVA("BDT method");  histBdt->Fill(wBDT);}
      
      if(wBDT>999) theTreeNew->Fill();
       
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
   
   // Write output file:
   fOutMVA->cd();
   theTreeNew->Write();
   if(evt_count) evt_count->Write();
   fOutMVA->Write();
   fOutMVA->Close();

   // Get efficiency for cuts classifier
  
   // Write histograms

   TFile *target  = new TFile( "TMVApp.root","RECREATE" );
   if (Use["BDT"          ])   histBdt    ->Write();
   
   target->Close();

   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}

int main( int argc, char** argv )
{
   if (argc < 2){
	   std::cout << "missing file!" << std::endl
	   << "Usage: " << argv[0] << " <inputMCFileName> [methods]" << std::endl;
    std::cout << "Exampe: " << argv[0] << " output_xi_errors_sgn.root BDT" << std::endl;
    return 1;
  }
  
   TString methodList;
   for (int i=2; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAClassificationApplication(argv[1],methodList);
   return 0;
}
