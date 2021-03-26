#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "TApplication.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TSystem.h"

using namespace std;

/*
The code was writen by Matteo Pisano for exclusive top analysis with 2017 CMS+PPS data
*/


//Funzione che calcola l'indice corrispondente all'entrata minima del vettore di input.
int dump_index(vector<double> sum_squared, int entries){
  int index=0;
  double a=0;
  double b=0;
  a=sum_squared[0];
  for(int i=1; i<entries; i++){
    b=sum_squared[i];
    if(b<=a){
      a=b;
      index=i;
    }
  }
  //cout << "index " << index << endl;
  return index;
}
double m_inv_t;

int main(int argc, char* argv[]){
  //TApplication app("app",NULL,NULL);
	
  int nMCEventsToProcess = -1;
  int nMCEventsToSkip = 0;

  // Check arguments
  if (argc < 2 || argc > 4) {
    cout << "ERROR: missing arguments" << endl
         << "Usage: " << argv[0] << " inputMCFileName [<nMCEventsToProcess>=all] [<nMCEventsToSkip>=0]" << endl;
    cout << "Exampe: " << argv[0] << " /eos/cms/store/group/phys_top/TTbarCentralExclProd/ntuples/mc/excl_ttbar_semilep_QED_xa120_era2017_preTS2.root" << endl;
    return 1;
  }
  

  // Check input files
  if(gSystem->AccessPathName(argv[1])){
	  cout << "ERROR: missing input file " << argv[1] << endl;
	  return 0;
  }

  TString infile(argv[1]);
  bool isData = infile.Contains("Data") || infile.Contains("SingleElectron") || infile.Contains("SingleMuon");
  
  //Preparing output file
  TString outfile = TString(infile.Tokenize("/")->Last()->GetName()).ReplaceAll(".root","_chi2.root");
  if (argc > 2) nMCEventsToProcess = TString(argv[2]).Atoi(); 
  if (argc > 3) {
	  nMCEventsToSkip = TString(argv[3]).Atoi(); 
	  outfile = outfile.ReplaceAll(".root",Form("_%d_%d.root",nMCEventsToSkip,nMCEventsToProcess));
  }  

  
  TFile output(outfile,"RECREATE","");
  TTree tree_out("tree","");
  double tot_dist=0;
  double entry1,entry2,entry3,entry4,entry5,entry6,entry7,entry8,entry9,entry10,entry11;
  double entry12,entry13,entry14,entry15,entry16,entry17,entry18,entry19,entry20,entry21;
  double entry22,entry23,entry24,entry25,entry26,entry27,entry28,entry29,entry30,entry31;
  double entry32,entry33,entry34,entry35,entry36,entry37,entry38,entry39,entry40,entry41;
  double entry42,entry43,entry44,entry45,entry46,entry47,entry48,entry49;
  tree_out.Branch("tbar_m",&entry1,"tbar_m/D");
  tree_out.Branch("t_pt",&entry2,"t_pt/D");
  tree_out.Branch("ttbar_pt",&entry3,"ttbar_pt/D");
  tree_out.Branch("l_pt",&entry4,"l_pt/D");
  tree_out.Branch("met_pt",&entry5,"met_pt/D");
  tree_out.Branch("chisquareAnalyticMattFit",&entry6,"chisquareAnalyticMattFit/D");
  tree_out.Branch("lightJet3_pt",&entry7,"lightJet3_pt/D");
  tree_out.Branch("lightJet4_pt",&entry8,"lightJet4_pt/D");
  tree_out.Branch("nJets",&entry9,"nJets/D");
  tree_out.Branch("t_eta",&entry10,"t_eta/D");
  tree_out.Branch("tbar_eta",&entry11,"tbar_eta/D");
  tree_out.Branch("t_phi",&entry12,"t_phi/D");
  tree_out.Branch("tbar_phi",&entry13,"tbar_phi/D");
  tree_out.Branch("bJet0_phi",&entry14,"bJet0_phi/D");
  tree_out.Branch("bJet1_phi",&entry15,"bJet1_phi/D");
  tree_out.Branch("bJet0_eta",&entry16,"bJet0_eta/D");
  tree_out.Branch("bJet1_eta",&entry17,"bJet1_eta/D");
  tree_out.Branch("bJet2_phi",&entry41,"bJet2_phi/D");
  tree_out.Branch("bJet3_phi",&entry42,"bJet3_phi/D");
  tree_out.Branch("bJet2_eta",&entry43,"bJet2_eta/D");
  tree_out.Branch("bJet3_eta",&entry44,"bJet3_eta/D");
  tree_out.Branch("lightJet0_E",&entry18,"lightJet0_E/D");
  tree_out.Branch("lightJet1_E",&entry19,"lightJet1_E/D");
  tree_out.Branch("lightJet2_E",&entry20,"lightJet2_E/D");
  tree_out.Branch("lightJet3_E",&entry21,"lightJet3_E/D");
  tree_out.Branch("l_E",&entry22,"l_E/D");
  tree_out.Branch("ttbar_E",&entry23,"ttbar_E/D");
  tree_out.Branch("bJet0_E",&entry24,"bJet0_E/D");
  tree_out.Branch("bJet1_E",&entry25,"bJet1_E/D");
  tree_out.Branch("bJet2_E",&entry26,"bJet2_E/D");
  tree_out.Branch("bJet3_E",&entry27,"bJet3_E/D");
  tree_out.Branch("bJet0_m",&entry28,"bJet0_m/D");
  tree_out.Branch("bJet1_m",&entry29,"bJet1_m/D");
  tree_out.Branch("bJet2_m",&entry34,"bJet2_m/D");
  tree_out.Branch("bJet3_m",&entry35,"bJet3_m/D");
  tree_out.Branch("lightJet0_m",&entry30,"lightJet0_m/D");
  tree_out.Branch("lightJet1_m",&entry31,"lightJet1_m/D");
  tree_out.Branch("lightJet2_m",&entry32,"lightJet2_m/D");
  tree_out.Branch("lightJet3_m",&entry33,"lightJet3_m/D");
  tree_out.Branch("ttbar_m",&entry36,"ttbar_m/D");
  tree_out.Branch("p_nu",&entry37,"p_nu/D");
  tree_out.Branch("p_l",&entry38,"p_l/D");
  tree_out.Branch("nLightJets",&entry39,"nLightJets/D");
  tree_out.Branch("lepton_isolation",&entry40,"lepton_isolation/D");
  tree_out.Branch("mean_ljets_delta_R",&entry45,"mean_ljets_delta_R/D");
  tree_out.Branch("min_ljets_delta_R",&entry46,"min_ljets_delta_R/D");
  tree_out.Branch("weight",&entry47,"weight/D");
  tree_out.Branch("p1_xi",&entry48,"p1_xi/D");
  tree_out.Branch("p2_xi",&entry49,"p2_xi/D");
  
  //MC weights variaitons to pass to a new output:
  unsigned int run; float beamXangle, pu_wgt, toppt_wgt, ptag_wgt_err, L1Prefire_wgt_err, ppsSF_wgt_err, trigSF_wgt_err, selSF_wgt_err; 
  float ren_err, fac_err;
  float pdf_as, pdf_hs;
  const int NPSRad_weights = 9;
  float isr_Up[NPSRad_weights], fsr_Up[NPSRad_weights], isr_Down[NPSRad_weights], fsr_Down[NPSRad_weights];

  
  int signal_protons =0, cat=0;
  tree_out.Branch("run",&run,"run/i");
  tree_out.Branch("beamXangle",&beamXangle);
  tree_out.Branch("pu_wgt",&pu_wgt);
  tree_out.Branch("toppt_wgt",&toppt_wgt);
  tree_out.Branch("ptag_wgt_err",&  ptag_wgt_err);
  tree_out.Branch("ppsSF_wgt_err",&  ppsSF_wgt_err);
  tree_out.Branch("trigSF_wgt_err",&trigSF_wgt_err);
  tree_out.Branch("selSF_wgt_err",&selSF_wgt_err);
  tree_out.Branch("L1Prefire_wgt_err",&L1Prefire_wgt_err);

  tree_out.Branch("pdf_as_err",&pdf_as);
  tree_out.Branch("pdf_hs_err",&pdf_hs);
  for(int i=0;i<NPSRad_weights;i++){
    tree_out.Branch(Form("isr%d_Up",i),&isr_Up[i]);
    tree_out.Branch(Form("fsr%d_Up",i),&fsr_Up[i]);
    tree_out.Branch(Form("isr%d_Down",i),&isr_Down[i]);
    tree_out.Branch(Form("fsr%d_Down",i),&fsr_Down[i]);
  }

  tree_out.Branch("ren_err",&ren_err);
  tree_out.Branch("fac_err",&fac_err);
  tree_out.Branch("signal_protons",&signal_protons);
  tree_out.Branch("cat",&cat);

  
  //Input files  
  //TFile file("../NewExcTopOUTPUT/bass.root");  
  //TFile file("./data/data.root");
  //TFile file("./mc/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
  TFile file(argv[1]);
  TTree* tree = (TTree*) file.Get("tree");
  TH1F * evt_count = (TH1F*)file.Get("evt_count");
  
  //Istogrammi che rappresentano i risultati
  
  TH1D histo("histo","histo", 300, 80.42 ,80.43);
  TH1D histo1("histo1","histo1", 300, 80.42 ,80.43);
  TH1D histo2("histo2","histo2", 300, 173. ,173.2);
  TH1D histo3("histo3","histo3", 300, 173. ,173.2);
  TH1D histo5("histo5","histo5", 300, 0 ,1000);
  TH1D histo4("histo4","histo4", 300, -1.5 ,1.5);
  TH1D histo6("histo6","histo6", 300, -0.1 ,0.1);
  TH1D histo7("histo7","histo7", 300, -0.1 ,0.1);
  TH1D histo8("histo8","histo8", 3000, -0.1 ,0.1);
  TH1D histo9("histo9","histo9", 300, -300 ,300);
  TH1D histo10("histo10","histo10", 300, -300 ,300);
  TH1D histo11("histo11","histo11", 300, 0 ,300);
  TH1D histo12("histo12","histo12", 2000, -6 ,6);
  //double hsump[5266];
  int occhio=0;
  TGraph profile;
  //Ciclo sugli eventi contenuti nel file di input
  
  int nevents = tree->GetEntries();
  std::cout << nevents << " events read from file(s) " << argv[1] << std::endl;
  if (nMCEventsToProcess == -1)
    nMCEventsToProcess = nevents;

  for(int h=0; (h<nMCEventsToProcess) && ((h + nMCEventsToSkip) < nevents); h++){ //68000
    if(h%1000==0) printf ("\r [%3.0f%%] done", 100.*(float)h/(float)nMCEventsToProcess);
	
    TMatrixD chi2(1,1);
    
    tree->GetEvent(h + nMCEventsToSkip);
    //int o=tree->GetEvent(h + nMCEventsToSkip);

    //Input variables
    
    double phi_b1=tree->GetLeaf("bJet0_phi")->GetValue(0);             
    //double phi_b2=tree->GetLeaf("bJet1_phi")->GetValue(0);
    
    //double phi_l=tree->GetLeaf("l_phi")->GetValue(0);                                
    //double phi_nu=tree->GetLeaf("nu_phi")->GetValue(0);                    
    //double phi_q1=tree->GetLeaf("lightJet0_phi")->GetValue(0);                  
    //double phi_q2=tree->GetLeaf("lightJet1_phi")->GetValue(0);                        
    //double eta_b1=tree->GetLeaf("bJet0_eta")->GetValue(0);
    //double eta_b2=tree->GetLeaf("bJet1_eta")->GetValue(0);
    //double eta_l=tree->GetLeaf("l_eta")->GetValue(0);
    //double eta_nu=tree->GetLeaf("nu_eta")->GetValue(0);
    //double eta_q1=tree->GetLeaf("lightJet0_eta")->GetValue(0);
    //double eta_q2=tree->GetLeaf("lightJet1_eta")->GetValue(0);
    double pt_b1=tree->GetLeaf("bJet0_pt")->GetValue(0);                  //adronico
    //double pt_b2=tree->GetLeaf("bJet1_pt")->GetValue(0);                  //leptonico
    //double pt_l=tree->GetLeaf("l_pt")->GetValue(0);                  
    double pt_nu=tree->GetLeaf("nu_pt")->GetValue(0);                              
    //double pt_q1=tree->GetLeaf("lightJet0_pt")->GetValue(0);                      
    //double pt_q2=tree->GetLeaf("lightJet1_pt")->GetValue(0);
    //double m_l=tree->GetLeaf("l_m")->GetValue(0);;   
    double m_b1=tree->GetLeaf("bJet0_m")->GetValue(0);
    double m_b2=tree->GetLeaf("bJet1_m")->GetValue(0);
    double m_q1= tree->GetLeaf("lightJet0_m")->GetValue(0);
    double m_q2= tree->GetLeaf("lightJet1_m")->GetValue(0);
    double pt_b3=tree->GetLeaf("bJet2_pt")->GetValue(0);
    double phi_b3=tree->GetLeaf("bJet2_phi")->GetValue(0);
    double eta_b3=tree->GetLeaf("bJet2_eta")->GetValue(0); 
    double m_b3=tree->GetLeaf("bJet2_m")->GetValue(0);
    double pt_b4=tree->GetLeaf("bJet3_pt")->GetValue(0);
    double phi_b4=tree->GetLeaf("bJet3_phi")->GetValue(0);
    double eta_b4=tree->GetLeaf("bJet3_eta")->GetValue(0);
    double m_b4=tree->GetLeaf("bJet3_m")->GetValue(0);
    double pt_q3=tree->GetLeaf("lightJet2_pt")->GetValue(0);
    double phi_q3=tree->GetLeaf("lightJet2_phi")->GetValue(0);
    double eta_q3=tree->GetLeaf("lightJet2_eta")->GetValue(0);
    double m_q3=tree->GetLeaf("lightJet2_m")->GetValue(0);
    double pt_q4=tree->GetLeaf("lightJet3_pt")->GetValue(0);
    double phi_q4=tree->GetLeaf("lightJet3_phi")->GetValue(0);
    double eta_q4=tree->GetLeaf("lightJet3_eta")->GetValue(0);
    double m_q4=tree->GetLeaf("lightJet3_m")->GetValue(0);
    //double pt_t_gen=tree->GetLeaf("gen_t_pt")->GetValue(0);
    //double pt_tbar_gen=tree->GetLeaf("gen_tbar_pt")->GetValue(0);
    double m_ttbar_gen=tree->GetLeaf("gen_ttbar_m")->GetValue(0);
    double xi_1=tree->GetLeaf("p1_xi")->GetValue(0);
    double xi_2=tree->GetLeaf("p2_xi")->GetValue(0); 
    double s=13000*13000;  
    int is_h= tree->GetLeaf("t_isHadronic")->GetValue(0);
    //double pt_b_gen=tree->GetLeaf("gen_b_pt")->GetValue(0);
    //double eta_b_gen=tree->GetLeaf("gen_b_eta")->GetValue(0);
    //double phi_b_gen=tree->GetLeaf("gen_b_phi")->GetValue(0);
    //double pt_bbar_gen=tree->GetLeaf("gen_bbar_pt")->GetValue(0);
    //double eta_bbar_gen=tree->GetLeaf("gen_bbar_eta")->GetValue(0);
    //double phi_bbar_gen=tree->GetLeaf("gen_bbar_phi")->GetValue(0);
    TF1 theta("theta", "2*atan(exp(-x))",0, 1e10);
    double e_b1_x=tree->GetLeaf("e_bJet_had_px")->GetValue(0);
    double e_b1_y=tree->GetLeaf("e_bJet_had_py")->GetValue(0);
    double e_b1_z=tree->GetLeaf("e_bJet_had_pz")->GetValue(0);    
    double e_b2_x=tree->GetLeaf("e_bJet_lep_px")->GetValue(0);
    double e_b2_y=tree->GetLeaf("e_bJet_lep_py")->GetValue(0);
    double e_b2_z=tree->GetLeaf("e_bJet_lep_pz")->GetValue(0);   
    double e_l_x=tree->GetLeaf("e_l_px")->GetValue(0);
    double e_l_y=tree->GetLeaf("e_l_py")->GetValue(0);
    double e_l_z=tree->GetLeaf("e_l_pz")->GetValue(0);
    double e_nu_x=6.*pt_nu;
    double e_nu_y=6.*pt_nu;
    double e_nu_z=6.*pt_nu;     
    double e_q1_x=tree->GetLeaf("e_lightJet0_px")->GetValue(0);
    double e_q1_y=tree->GetLeaf("e_lightJet0_py")->GetValue(0);
    double e_q1_z=tree->GetLeaf("e_lightJet0_pz")->GetValue(0);   
    double e_q2_x=tree->GetLeaf("e_lightJet1_px")->GetValue(0);
    double e_q2_y=tree->GetLeaf("e_lightJet1_py")->GetValue(0);
    double e_q2_z=tree->GetLeaf("e_lightJet1_pz")->GetValue(0);

    
    ////Traduzione Variabili di input in componenti x,y,z (nomi definiti nel pdf di riferimento)
  
    double bar_p_x_b1= tree->GetLeaf("bJet_had_px")->GetValue(0);   //b1 è adronico
    double bar_p_y_b1= tree->GetLeaf("bJet_had_py")->GetValue(0);
    double bar_p_z_b1= tree->GetLeaf("bJet_had_pz")->GetValue(0);
    double bar_p_x_b2= tree->GetLeaf("bJet_lep_px")->GetValue(0);   //b2 è leptonico
    double bar_p_y_b2= tree->GetLeaf("bJet_lep_py")->GetValue(0);
    double bar_p_z_b2= tree->GetLeaf("bJet_lep_pz")->GetValue(0);
    double bar_p_x_b3= pt_b3*cos(phi_b3);   //b1 è adronico
    double bar_p_y_b3= pt_b3*sin(phi_b3);
    double bar_p_z_b3= pt_b3/tan(theta.Eval(eta_b3));
    double bar_p_x_b4= pt_b4*cos(phi_b4);   //b2 è leptonico
    double bar_p_y_b4= pt_b4*sin(phi_b4);
    double bar_p_z_b4= pt_b4/tan(theta.Eval(eta_b4));
    double bar_p_x_l= tree->GetLeaf("l_px")->GetValue(0);
    double bar_p_y_l= tree->GetLeaf("l_py")->GetValue(0);
    double bar_p_z_l= tree->GetLeaf("l_pz")->GetValue(0);
    double bar_p_x_nu= tree->GetLeaf("nu_px")->GetValue(0);
    double bar_p_y_nu= tree->GetLeaf("nu_py")->GetValue(0);
    double bar_p_z_nu= tree->GetLeaf("nu_pz")->GetValue(0);;
    double bar_p_x_q1= tree->GetLeaf("lightJet0_px")->GetValue(0);
    double bar_p_y_q1= tree->GetLeaf("lightJet0_py")->GetValue(0);
    double bar_p_z_q1= tree->GetLeaf("lightJet0_pz")->GetValue(0);
    double bar_p_x_q2= tree->GetLeaf("lightJet1_px")->GetValue(0);
    double bar_p_y_q2= tree->GetLeaf("lightJet1_py")->GetValue(0);
    double bar_p_z_q2= tree->GetLeaf("lightJet0_pz")->GetValue(0);
    double bar_p_x_q3= pt_q3*cos(phi_q3);
    double bar_p_y_q3= pt_q3*sin(phi_q3);;
    double bar_p_z_q3= pt_q3/tan(theta.Eval(eta_q3));
    double bar_p_x_q4= pt_q4*cos(phi_q4);
    double bar_p_y_q4= pt_q4*cos(phi_q4);
    double bar_p_z_q4= pt_q4/tan(theta.Eval(eta_q4));
    double bar_xi_1= xi_1;
    double bar_xi_2= xi_2;

    int max=50; //Setta il numero di iterazioni primo dello stop (set number of iterations)

    if(xi_1==0 || xi_2==0) max=0;
	
	//Read weights and variables to pass to the new output:
	cat = int(tree->GetLeaf("cat")->GetValue(0));
	run = tree->GetLeaf("run")->GetValue(0);
	beamXangle = tree->GetLeaf("beamXangle")->GetValue(0);
	pu_wgt = tree->GetLeaf("pu_wgt")->GetValue(0);
	toppt_wgt = tree->GetLeaf("toppt_wgt")->GetValue(0);
	ptag_wgt_err = tree->GetLeaf("ptag_wgt_err")->GetValue(0);
	ppsSF_wgt_err = tree->GetLeaf("ppsSF_wgt_err")->GetValue(0);
	trigSF_wgt_err = tree->GetLeaf("trigSF_wgt_err")->GetValue(0);
	selSF_wgt_err =tree->GetLeaf("selSF_wgt_err")->GetValue(0); 
	L1Prefire_wgt_err =tree->GetLeaf("L1Prefire_wgt_err")->GetValue(0); 
	
	pdf_as =tree->GetLeaf("pdf_as")->GetValue(0); 
	pdf_hs =tree->GetLeaf("pdf_hs")->GetValue(0); 

	// add guard from nans, and large values (broken weights)
	if(pdf_as!=pdf_as) pdf_as=0;
	if(pdf_hs!=pdf_hs) pdf_hs=0;
	
	for(int i=0;i<NPSRad_weights;i++){
		isr_Up[i]=tree->GetLeaf("isr_Up")->GetValue(i);
		fsr_Up[i]=tree->GetLeaf("fsr_Up")->GetValue(i);
		isr_Down[i]=tree->GetLeaf("isr_Down")->GetValue(i);
		fsr_Down[i]=tree->GetLeaf("fsr_Down")->GetValue(i);
	}	

	ren_err =tree->GetLeaf("ren_err")->GetValue(0); 
	fac_err =tree->GetLeaf("fac_err")->GetValue(0); 
	signal_protons = (!isData) ? tree->GetLeaf("signal_protons")->GetValue(0) : 0; 
  
    ////Definisco il vettore di traslazione
    
    TMatrixD f(27,1);
    TMatrixD f_iniziale(27,1);
    f.Zero();
    TMatrixD C(20,20); //Matrice di covarianza
    TMatrixD C_sum(20,20);
    TMatrixD InvC(20,20);
    TMatrixD C_star_i(27,27);
    TMatrixD P_star_i(27,27);
    TMatrixD Q_star_i(27,1);
    C_star_i.Zero();
    C_sum.Zero();
    
    for(int iteration=0; iteration<max; iteration++){

      ///////////////////////////////////////////////////////////////////////NEW
      //double dump=1.;
      //int step=0;
      //int escape=0;
      ///////////////////////////////////////////////////////////////////////END/NEW

    

    if(iteration==0){
    C[0][0]=pow(abs(e_b1_x),2);
    C[1][1]=pow(abs(e_b1_y),2);
    C[2][2]=pow(abs(e_b1_z),2);
    C[3][3]=pow(abs(e_b2_x),2);
    C[4][4]=pow(abs(e_b2_y),2);
    C[5][5]=pow(abs(e_b2_z),2);
    C[6][6]=pow(abs(e_l_x),2);
    C[7][7]=pow(abs(e_l_y),2);
    C[8][8]=pow(abs(e_l_z),2);
    C[9][9]=pow(abs(e_nu_x),2);
    C[10][10]=pow(abs(e_nu_y),2);
    C[11][11]=pow(abs(e_nu_z),2);
    C[12][12]=pow(abs(e_q1_x),2);
    C[13][13]=pow(abs(e_q1_y),2);
    C[14][14]=pow(abs(e_q1_z),2);
    C[15][15]=pow(abs(e_q2_x),2);
    C[16][16]=pow(abs(e_q2_y),2);
    C[17][17]=pow(abs(e_q2_z),2);
    //C[18][18]=pow(0.0164*bar_xi_1+0.00129,2);
    //C[19][19]=pow(0.0152*bar_xi_2+0.00130,2);
    // preTS2
	//C[18][18]=pow(8585.*pow(bar_xi_1,5)-2896.*pow(bar_xi_1,4)+374.*pow(bar_xi_1,3)-23.07*pow(bar_xi_1,2)+0.745*bar_xi_1-0.0067,2);
	//C[19][19]=pow(111.*pow(bar_xi_2,4)-31.41*pow(bar_xi_2,3)+2.823*pow(bar_xi_2,2)-0.014*bar_xi_2+0.00146,2);
	// postTS2
	C[18][18]=pow(3424.*pow(bar_xi_1,5)-1273.*pow(bar_xi_1,4)+182.*pow(bar_xi_1,3)-12.51*pow(bar_xi_1,2)+0.478*bar_xi_1-0.0043,2);
	C[19][19]=pow(10.7*pow(bar_xi_2,4)+1.637*pow(bar_xi_2,3)-1.075*pow(bar_xi_2,2)+0.0176*bar_xi_2-0.00176,2);
    
    }
    
    InvC=C.Invert();
    

    double m_t=173.1; //Massa del top
    double m_tbar=173.1; //Massa del top
    double m_W=80.426; //Massa del W
    ////

    //// Definisco le matrici trattate nel pdf di riferimento

    TMatrixD C_star(27,27);
    TMatrixD K_star(27,27);
    TMatrixD P_star(27,27);
    TMatrixD Q_star(27,1);

    C_star.Zero();  //Riempio di zeri tutte le matrici.
    K_star.Zero();
    P_star.Zero();
    Q_star.Zero();

    ////Definisco il vettore contenente le quantità da fittare;
    TMatrixD x(27,1);
    TMatrixD fittedsolutions(27,1);
    x.Zero();
    fittedsolutions.Zero();

 
    //// Riempio gli oggetti con gli input

    //Riempizione C_star
    for(int row=0; row<20; row++){
      for(int col=0; col<20; col++){
	C_star[row][col]=InvC[row][col];
      }
    }

    if(iteration==0){
      C_star_i=C_star;
    }

    
    //if(h==37) C_star_i.Print(); 

    //Riempizione K_star
    K_star[0][20]=1.;
    K_star[1][21]=1.;
    K_star[3][20]=1.;
    K_star[4][21]=1.;
    K_star[6][20]=1.;
    K_star[7][21]=1.;
    K_star[9][20]=1.;
    K_star[10][21]=1.;
    K_star[12][20]=1.;
    K_star[13][21]=1.;
    K_star[15][20]=1.;
    K_star[16][21]=1.;

    K_star[20][0]=1.;
    K_star[21][1]=1.;
    K_star[20][3]=1.;
    K_star[21][4]=1.;
    K_star[20][6]=1.;
    K_star[21][7]=1.;
    K_star[20][9]=1.;
    K_star[21][10]=1.;
    K_star[20][12]=1.;
    K_star[21][13]=1.;
    K_star[20][15]=1.;
    K_star[21][16]=1.;

    //Riempizione Q_star

    //Per calcolare le entrate di Q_star devo prima calcolare
    //i valori della massa invariante dei sistemi adronici
    //e leptonici usando i valori di p non fittati.

    double mod_bar_p_b1=sqrt(pow(bar_p_x_b1,2)+pow(bar_p_y_b1,2)+pow(bar_p_z_b1,2)); //Moduli trivettori p
    double mod_bar_p_b2=sqrt(pow(bar_p_x_b2,2)+pow(bar_p_y_b2,2)+pow(bar_p_z_b2,2));
    double mod_bar_p_b3=sqrt(pow(bar_p_x_b3,2)+pow(bar_p_y_b3,2)+pow(bar_p_z_b3,2));
    double mod_bar_p_b4=sqrt(pow(bar_p_x_b4,2)+pow(bar_p_y_b4,2)+pow(bar_p_z_b4,2));
    double mod_bar_p_l=sqrt(pow(bar_p_x_l,2)+pow(bar_p_y_l,2)+pow(bar_p_z_l,2));
    double mod_bar_p_nu=sqrt(pow(bar_p_x_nu,2)+pow(bar_p_y_nu,2)+pow(bar_p_z_nu,2));
    double mod_bar_p_q1=sqrt(pow(bar_p_x_q1,2)+pow(bar_p_y_q1,2)+pow(bar_p_z_q1,2));
    double mod_bar_p_q2=sqrt(pow(bar_p_x_q2,2)+pow(bar_p_y_q2,2)+pow(bar_p_z_q2,2));
    double mod_bar_p_q3=sqrt(pow(bar_p_x_q3,2)+pow(bar_p_y_q3,2)+pow(bar_p_z_q3,2));
    double mod_bar_p_q4=sqrt(pow(bar_p_x_q4,2)+pow(bar_p_y_q4,2)+pow(bar_p_z_q4,2));
  
    double bar_m_lep=sqrt(2*mod_bar_p_b2*mod_bar_p_l+2*mod_bar_p_b2*mod_bar_p_nu+2*mod_bar_p_l*mod_bar_p_nu-2*bar_p_x_b2*bar_p_x_l-2*bar_p_x_b2*bar_p_x_nu-2*bar_p_x_l*bar_p_x_nu-2*bar_p_y_b2*bar_p_y_l-2*bar_p_y_b2*bar_p_y_nu-2*bar_p_y_l*bar_p_y_nu-2*bar_p_z_b2*bar_p_z_l-2*bar_p_z_b2*bar_p_z_nu-2*bar_p_z_l*bar_p_z_nu);
  
    double bar_m_ad=sqrt(2*mod_bar_p_b1*mod_bar_p_q1+2*mod_bar_p_b1*mod_bar_p_q2+2*mod_bar_p_q1*mod_bar_p_q2-2*bar_p_x_b1*bar_p_x_q1-2*bar_p_x_b1*bar_p_x_q2-2*bar_p_x_q1*bar_p_x_q2-2*bar_p_y_b1*bar_p_y_q1-2*bar_p_y_b1*bar_p_y_q2-2*bar_p_y_q1*bar_p_y_q2-2*bar_p_z_b1*bar_p_z_q1-2*bar_p_z_b1*bar_p_z_q2-2*bar_p_z_q1*bar_p_z_q2);

    double bar_m_lepnu = sqrt(2*mod_bar_p_l*mod_bar_p_nu-2*bar_p_x_l*bar_p_x_nu-2*bar_p_y_l*bar_p_y_nu-2*bar_p_z_l*bar_p_z_nu);
    double bar_m_q1q2 = sqrt(2*mod_bar_p_q1*mod_bar_p_q2-2*bar_p_x_q1*bar_p_x_q2-2*bar_p_y_q1*bar_p_y_q2-2*bar_p_z_q1*bar_p_z_q2);

    TLorentzVector test_b1,test_b2,test_l,test_nu,test_q1,test_q2,test_b3,test_b4,test_q3,test_q4;
    test_b1.SetPxPyPzE(bar_p_x_b1,bar_p_y_b1,bar_p_z_b1,mod_bar_p_b1);
    test_b2.SetPxPyPzE(bar_p_x_b2,bar_p_y_b2,bar_p_z_b2,mod_bar_p_b2);
    test_b3.SetPxPyPzE(bar_p_x_b3,bar_p_y_b3,bar_p_z_b3,mod_bar_p_b3);
    test_b4.SetPxPyPzE(bar_p_x_b4,bar_p_y_b4,bar_p_z_b4,mod_bar_p_b4);
    test_l.SetPxPyPzE(bar_p_x_l,bar_p_y_l,bar_p_z_l,mod_bar_p_l);
    test_nu.SetPxPyPzE(bar_p_x_nu,bar_p_y_nu,bar_p_z_nu,mod_bar_p_nu);
    test_q1.SetPxPyPzE(bar_p_x_q1,bar_p_y_q1,bar_p_z_q1,mod_bar_p_q1);
    test_q2.SetPxPyPzE(bar_p_x_q2,bar_p_y_q2,bar_p_z_q2,mod_bar_p_q2);
    test_q3.SetPxPyPzE(bar_p_x_q3,bar_p_y_q3,bar_p_z_q3,mod_bar_p_q3);
    test_q4.SetPxPyPzE(bar_p_x_q4,bar_p_y_q4,bar_p_z_q4,mod_bar_p_q4);
    TLorentzVector ttbar= test_b1+test_b2+test_l+test_nu+test_q1+test_q2;
    TLorentzVector lep=test_b2+test_l+test_nu;
    TLorentzVector ad=test_b1+test_q1+test_q2;
    TLorentzVector lepnu=test_l+test_nu;
    TLorentzVector q1q2= test_q1+test_q2;
    bar_m_lep=lep.M();
    bar_m_ad=ad.M();
    bar_m_lepnu=lepnu.M();
    bar_m_q1q2=q1q2.M();

    
    //if(iteration==999) cout << bar_m_ad << endl;
    Q_star[22][0]=(pow(bar_m_lep,2)-m_t*m_t)/(m_t*m_t);
    Q_star[23][0]=(pow(bar_m_ad,2)-m_tbar*m_tbar)/(m_tbar*m_tbar);
    Q_star[24][0]=(pow(bar_m_lepnu,2)-m_W*m_W)/(m_W*m_W);
    Q_star[25][0]=(pow(bar_m_q1q2,2)-m_W*m_W)/(m_W*m_W);

    
    Q_star[26][0]=ttbar.M()*ttbar.M()/s-bar_xi_1*bar_xi_2;
    //Q_star[25][0]=0; 

    if(iteration==0) Q_star_i=Q_star;
    

    //Riempizione f
    
    f[0][0]=bar_p_x_b1;
    f[1][0]=bar_p_y_b1;
    f[2][0]=bar_p_z_b1;
    f[3][0]=bar_p_x_b2;
    f[4][0]=bar_p_y_b2;
    f[5][0]=bar_p_z_b2;
    f[6][0]=bar_p_x_l;
    f[7][0]=bar_p_y_l;
    f[8][0]=bar_p_z_l;
    f[9][0]=bar_p_x_nu;
    f[10][0]=bar_p_y_nu;
    f[11][0]=bar_p_z_nu;
    f[12][0]=bar_p_x_q1;
    f[13][0]=bar_p_y_q1;
    f[14][0]=bar_p_z_q1;
    f[15][0]=bar_p_x_q2;
    f[16][0]=bar_p_y_q2;
    f[17][0]=bar_p_z_q2;
    f[18][0]=bar_xi_1;
    f[19][0]=bar_xi_2;

    if(iteration==0){
    f_iniziale[0][0]=bar_p_x_b1;
    f_iniziale[1][0]=bar_p_y_b1;
    f_iniziale[2][0]=bar_p_z_b1;
    f_iniziale[3][0]=bar_p_x_b2;
    f_iniziale[4][0]=bar_p_y_b2;
    f_iniziale[5][0]=bar_p_z_b2;
    f_iniziale[6][0]=bar_p_x_l;
    f_iniziale[7][0]=bar_p_y_l;
    f_iniziale[8][0]=bar_p_z_l;
    f_iniziale[9][0]=bar_p_x_nu;
    f_iniziale[10][0]=bar_p_y_nu;
    f_iniziale[11][0]=bar_p_z_nu;
    f_iniziale[12][0]=bar_p_x_q1;
    f_iniziale[13][0]=bar_p_y_q1;
    f_iniziale[14][0]=bar_p_z_q1;
    f_iniziale[15][0]=bar_p_x_q2;
    f_iniziale[16][0]=bar_p_y_q2;
    f_iniziale[17][0]=bar_p_z_q2;
    f_iniziale[18][0]=bar_xi_1;
    f_iniziale[19][0]=bar_xi_2;
    }
    

    //Riempizione P_star;

    P_star[3][22]=2*(mod_bar_p_l*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_x_b2/mod_bar_p_b2-bar_p_x_l-bar_p_x_nu)/(m_t*m_t);
    P_star[4][22]=2*(mod_bar_p_l*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_y_b2/mod_bar_p_b2-bar_p_y_l-bar_p_y_nu)/(m_t*m_t);
    P_star[5][22]=2*(mod_bar_p_l*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_z_b2/mod_bar_p_b2-bar_p_z_l-bar_p_z_nu)/(m_t*m_t);
    P_star[6][22]=2*(mod_bar_p_b2*bar_p_x_l/mod_bar_p_l+mod_bar_p_nu*bar_p_x_l/mod_bar_p_l-bar_p_x_b2-bar_p_x_nu)/(m_t*m_t);
    P_star[7][22]=2*(mod_bar_p_b2*bar_p_y_l/mod_bar_p_l+mod_bar_p_nu*bar_p_y_l/mod_bar_p_l-bar_p_y_b2-bar_p_y_nu)/(m_t*m_t);
    P_star[8][22]=2*(mod_bar_p_b2*bar_p_z_l/mod_bar_p_l+mod_bar_p_nu*bar_p_z_l/mod_bar_p_l-bar_p_z_b2-bar_p_z_nu)/(m_t*m_t);
    P_star[9][22]=2*(mod_bar_p_b2*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu-bar_p_x_b2-bar_p_x_l)/(m_t*m_t);
    P_star[10][22]=2*(mod_bar_p_b2*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu-bar_p_y_b2-bar_p_y_l)/(m_t*m_t);
    P_star[11][22]=2*(mod_bar_p_b2*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu-bar_p_z_b2-bar_p_z_l)/(m_t*m_t);
    P_star[0][23]=2*(mod_bar_p_q1*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_x_b1/mod_bar_p_b1-bar_p_x_q1-bar_p_x_q2)/(m_tbar*m_tbar);
    P_star[1][23]=2*(mod_bar_p_q1*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_y_b1/mod_bar_p_b1-bar_p_y_q1-bar_p_y_q2)/(m_tbar*m_tbar);
    P_star[2][23]=2*(mod_bar_p_q1*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_z_b1/mod_bar_p_b1-bar_p_z_q1-bar_p_z_q2)/(m_tbar*m_tbar);
    P_star[12][23]=2*(mod_bar_p_b1*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_b1-bar_p_x_q2)/(m_tbar*m_tbar);
    P_star[13][23]=2*(mod_bar_p_b1*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_b1-bar_p_y_q2)/(m_tbar*m_tbar);
    P_star[14][23]=2*(mod_bar_p_b1*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_b1-bar_p_z_q2)/(m_tbar*m_tbar);
    P_star[15][23]=2*(mod_bar_p_b1*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_x_q2/mod_bar_p_q2-bar_p_x_b1-bar_p_x_q1)/(m_tbar*m_tbar);
    P_star[16][23]=2*(mod_bar_p_b1*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_y_q2/mod_bar_p_q2-bar_p_y_b1-bar_p_y_q1)/(m_tbar*m_tbar);
    P_star[17][23]=2*(mod_bar_p_b1*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_z_q2/mod_bar_p_q2-bar_p_z_b1-bar_p_z_q1)/(m_tbar*m_tbar);
    P_star[6][24]=2*(mod_bar_p_nu*bar_p_x_l/mod_bar_p_l-bar_p_x_nu)/(m_W*m_W);
    P_star[7][24]=2*(mod_bar_p_nu*bar_p_y_l/mod_bar_p_l-bar_p_y_nu)/(m_W*m_W);
    P_star[8][24]=2*(mod_bar_p_nu*bar_p_z_l/mod_bar_p_l-bar_p_z_nu)/(m_W*m_W);
    P_star[9][24]=2*(mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu-bar_p_x_l)/(m_W*m_W);
    P_star[10][24]=2*(mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu-bar_p_y_l)/(m_W*m_W);
    P_star[11][24]=2*(mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu-bar_p_z_l)/(m_W*m_W);
    P_star[12][25]=2*(mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_q2)/(m_W*m_W);
    P_star[13][25]=2*(mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_q2)/(m_W*m_W);
    P_star[14][25]=2*(mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_q2)/(m_W*m_W);
    P_star[15][25]=2*(mod_bar_p_q1*bar_p_x_q2/mod_bar_p_q2-bar_p_x_q1)/(m_W*m_W);
    P_star[16][25]=2*(mod_bar_p_q1*bar_p_y_q2/mod_bar_p_q2-bar_p_y_q1)/(m_W*m_W);
    P_star[17][25]=2*(mod_bar_p_q1*bar_p_z_q2/mod_bar_p_q2-bar_p_z_q1)/(m_W*m_W);
    P_star[0][26]=2*(mod_bar_p_b2*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_l*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_nu*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_q1*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_x_b1/mod_bar_p_b1-bar_p_x_b2-bar_p_x_l-bar_p_x_nu-bar_p_x_q1-bar_p_x_q2)/s;
    P_star[1][26]=2*(mod_bar_p_b2*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_l*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_nu*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_q1*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_y_b1/mod_bar_p_b1-bar_p_y_b2-bar_p_y_l-bar_p_y_nu-bar_p_y_q1-bar_p_y_q2)/s;
    P_star[2][26]=2*(mod_bar_p_b2*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_l*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_nu*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_q1*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_z_b1/mod_bar_p_b1-bar_p_z_b2-bar_p_z_l-bar_p_z_nu-bar_p_z_q1-bar_p_z_q2)/s;
    P_star[3][26]=2*(mod_bar_p_b1*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_l*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_q1*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_q2*bar_p_x_b2/mod_bar_p_b2-bar_p_x_b1-bar_p_x_l-bar_p_x_nu-bar_p_x_q1-bar_p_x_q2)/s;
    P_star[4][26]=2*(mod_bar_p_b1*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_l*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_q1*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_q2*bar_p_y_b2/mod_bar_p_b2-bar_p_y_b1-bar_p_y_l-bar_p_y_nu-bar_p_y_q1-bar_p_y_q2)/s;
    P_star[5][26]=2*(mod_bar_p_b1*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_l*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_q1*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_q2*bar_p_z_b2/mod_bar_p_b2-bar_p_z_b1-bar_p_z_l-bar_p_z_nu-bar_p_z_q1-bar_p_z_q2)/s;
    P_star[6][26]=2*(mod_bar_p_b1*bar_p_x_l/mod_bar_p_l+mod_bar_p_b2*bar_p_x_l/mod_bar_p_l+mod_bar_p_nu*bar_p_x_l/mod_bar_p_l+mod_bar_p_q1*bar_p_x_l/mod_bar_p_l+mod_bar_p_q2*bar_p_x_l/mod_bar_p_l-bar_p_x_b1-bar_p_x_b2-bar_p_x_nu-bar_p_x_q1-bar_p_x_q2)/s;
    P_star[7][26]=2*(mod_bar_p_b1*bar_p_y_l/mod_bar_p_l+mod_bar_p_b2*bar_p_y_l/mod_bar_p_l+mod_bar_p_nu*bar_p_y_l/mod_bar_p_l+mod_bar_p_q1*bar_p_y_l/mod_bar_p_l+mod_bar_p_q2*bar_p_y_l/mod_bar_p_l-bar_p_y_b1-bar_p_y_b2-bar_p_y_nu-bar_p_y_q1-bar_p_y_q2)/s;
    P_star[8][26]=2*(mod_bar_p_b1*bar_p_z_l/mod_bar_p_l+mod_bar_p_b2*bar_p_z_l/mod_bar_p_l+mod_bar_p_nu*bar_p_z_l/mod_bar_p_l+mod_bar_p_q1*bar_p_z_l/mod_bar_p_l+mod_bar_p_q2*bar_p_z_l/mod_bar_p_l-bar_p_z_b1-bar_p_z_b2-bar_p_z_nu-bar_p_z_q1-bar_p_z_q2)/s;
    P_star[9][26]=2*(mod_bar_p_b1*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_b2*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_q1*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_q2*bar_p_x_l/mod_bar_p_nu-bar_p_x_b1-bar_p_x_b2-bar_p_x_l-bar_p_x_q1-bar_p_x_q2)/s;
    P_star[10][26]=2*(mod_bar_p_b1*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_b2*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_q1*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_q2*bar_p_y_l/mod_bar_p_nu-bar_p_y_b1-bar_p_y_b2-bar_p_y_l-bar_p_y_q1-bar_p_y_q2)/s;
    P_star[11][26]=2*(mod_bar_p_b1*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_b2*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_q1*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_q2*bar_p_z_l/mod_bar_p_nu-bar_p_z_b1-bar_p_z_b2-bar_p_z_l-bar_p_z_q1-bar_p_z_q2)/s;
    P_star[12][26]=2*(mod_bar_p_b1*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_b2*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_l*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_nu*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_b1-bar_p_x_b2-bar_p_x_l-bar_p_x_nu-bar_p_x_q2)/s;
    P_star[13][26]=2*(mod_bar_p_b1*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_b2*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_l*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_nu*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_b1-bar_p_y_b2-bar_p_y_l-bar_p_y_nu-bar_p_y_q2)/s;
    P_star[14][26]=2*(mod_bar_p_b1*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_b2*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_l*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_nu*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_b1-bar_p_z_b2-bar_p_z_l-bar_p_z_nu-bar_p_z_q2)/s;
    P_star[15][26]=2*(mod_bar_p_b1*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_b2*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_l*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_nu*bar_p_x_q1/mod_bar_p_q2+mod_bar_p_q2*bar_p_x_q2/mod_bar_p_q2-bar_p_x_b1-bar_p_x_b2-bar_p_x_l-bar_p_x_nu-bar_p_x_q1)/s;
    P_star[16][26]=2*(mod_bar_p_b1*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_b2*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_l*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_nu*bar_p_y_q1/mod_bar_p_q2+mod_bar_p_q2*bar_p_y_q2/mod_bar_p_q2-bar_p_y_b1-bar_p_y_b2-bar_p_y_l-bar_p_y_nu-bar_p_y_q1)/s;
    P_star[17][26]=2*(mod_bar_p_b1*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_b2*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_l*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_nu*bar_p_z_q1/mod_bar_p_q2+mod_bar_p_q2*bar_p_z_q2/mod_bar_p_q2-bar_p_z_b1-bar_p_z_b2-bar_p_z_l-bar_p_z_nu-bar_p_z_q1)/s;
    P_star[18][26]=-bar_xi_2;
    P_star[19][26]=-bar_xi_1;
    
    
    
    


    P_star[22][3]=2*(mod_bar_p_l*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_x_b2/mod_bar_p_b2-bar_p_x_l-bar_p_x_nu)/(m_t*m_t);
    P_star[22][4]=2*(mod_bar_p_l*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_y_b2/mod_bar_p_b2-bar_p_y_l-bar_p_y_nu)/(m_t*m_t);
    P_star[22][5]=2*(mod_bar_p_l*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_z_b2/mod_bar_p_b2-bar_p_z_l-bar_p_z_nu)/(m_t*m_t);
    P_star[22][6]=2*(mod_bar_p_b2*bar_p_x_l/mod_bar_p_l+mod_bar_p_nu*bar_p_x_l/mod_bar_p_l-bar_p_x_b2-bar_p_x_nu)/(m_t*m_t);
    P_star[22][7]=2*(mod_bar_p_b2*bar_p_y_l/mod_bar_p_l+mod_bar_p_nu*bar_p_y_l/mod_bar_p_l-bar_p_y_b2-bar_p_y_nu)/(m_t*m_t);
    P_star[22][8]=2*(mod_bar_p_b2*bar_p_z_l/mod_bar_p_l+mod_bar_p_nu*bar_p_z_l/mod_bar_p_l-bar_p_z_b2-bar_p_z_nu)/(m_t*m_t);
    P_star[22][9]=2*(mod_bar_p_b2*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu-bar_p_x_b2-bar_p_x_l)/(m_t*m_t);
    P_star[22][10]=2*(mod_bar_p_b2*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu-bar_p_y_b2-bar_p_y_l)/(m_t*m_t);
    P_star[22][11]=2*(mod_bar_p_b2*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu-bar_p_z_b2-bar_p_z_l)/(m_t*m_t);
    P_star[23][0]=2*(mod_bar_p_q1*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_x_b1/mod_bar_p_b1-bar_p_x_q1-bar_p_x_q2)/(m_tbar*m_tbar);
    P_star[23][1]=2*(mod_bar_p_q1*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_y_b1/mod_bar_p_b1-bar_p_y_q1-bar_p_y_q2)/(m_tbar*m_tbar);
    P_star[23][2]=2*(mod_bar_p_q1*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_z_b1/mod_bar_p_b1-bar_p_z_q1-bar_p_z_q2)/(m_tbar*m_tbar);
    P_star[23][12]=2*(mod_bar_p_b1*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_b1-bar_p_x_q2)/(m_tbar*m_tbar);
    P_star[23][13]=2*(mod_bar_p_b1*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_b1-bar_p_y_q2)/(m_tbar*m_tbar);
    P_star[23][14]=2*(mod_bar_p_b1*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_b1-bar_p_z_q2)/(m_tbar*m_tbar);
    P_star[23][15]=2*(mod_bar_p_b1*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_x_q2/mod_bar_p_q2-bar_p_x_b1-bar_p_x_q1)/(m_tbar*m_tbar);
    P_star[23][16]=2*(mod_bar_p_b1*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_y_q2/mod_bar_p_q2-bar_p_y_b1-bar_p_y_q1)/(m_tbar*m_tbar);
    P_star[23][17]=2*(mod_bar_p_b1*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_z_q2/mod_bar_p_q2-bar_p_z_b1-bar_p_z_q1)/(m_tbar*m_tbar);
    P_star[24][6]=2*(mod_bar_p_nu*bar_p_x_l/mod_bar_p_l-bar_p_x_nu)/(m_W*m_W);
    P_star[24][7]=2*(mod_bar_p_nu*bar_p_y_l/mod_bar_p_l-bar_p_y_nu)/(m_W*m_W);
    P_star[24][8]=2*(mod_bar_p_nu*bar_p_z_l/mod_bar_p_l-bar_p_z_nu)/(m_W*m_W);
    P_star[24][9]=2*(mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu-bar_p_x_l)/(m_W*m_W);
    P_star[24][10]=2*(mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu-bar_p_y_l)/(m_W*m_W);
    P_star[24][11]=2*(mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu-bar_p_z_l)/(m_W*m_W);
    P_star[25][12]=2*(mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_q2)/(m_W*m_W);
    P_star[25][13]=2*(mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_q2)/(m_W*m_W);
    P_star[25][14]=2*(mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_q2)/(m_W*m_W);
    P_star[25][15]=2*(mod_bar_p_q1*bar_p_x_q2/mod_bar_p_q2-bar_p_x_q1)/(m_W*m_W);
    P_star[25][16]=2*(mod_bar_p_q1*bar_p_y_q2/mod_bar_p_q2-bar_p_y_q1)/(m_W*m_W);
    P_star[25][17]=2*(mod_bar_p_q1*bar_p_z_q2/mod_bar_p_q2-bar_p_z_q1)/(m_W*m_W);
    P_star[26][0]=2*(mod_bar_p_b2*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_l*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_nu*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_q1*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_x_b1/mod_bar_p_b1-bar_p_x_b2-bar_p_x_l-bar_p_x_nu-bar_p_x_q1-bar_p_x_q2)/s;
    P_star[26][1]=2*(mod_bar_p_b2*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_l*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_nu*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_q1*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_y_b1/mod_bar_p_b1-bar_p_y_b2-bar_p_y_l-bar_p_y_nu-bar_p_y_q1-bar_p_y_q2)/s;
    P_star[26][2]=2*(mod_bar_p_b2*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_l*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_nu*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_q1*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_z_b1/mod_bar_p_b1-bar_p_z_b2-bar_p_z_l-bar_p_z_nu-bar_p_z_q1-bar_p_z_q2)/s;
    P_star[26][3]=2*(mod_bar_p_b1*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_l*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_q1*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_q2*bar_p_x_b2/mod_bar_p_b2-bar_p_x_b1-bar_p_x_l-bar_p_x_nu-bar_p_x_q1-bar_p_x_q2)/s;
    P_star[26][4]=2*(mod_bar_p_b1*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_l*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_q1*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_q2*bar_p_y_b2/mod_bar_p_b2-bar_p_y_b1-bar_p_y_l-bar_p_y_nu-bar_p_y_q1-bar_p_y_q2)/s;
    P_star[26][5]=2*(mod_bar_p_b1*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_l*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_q1*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_q2*bar_p_z_b2/mod_bar_p_b2-bar_p_z_b1-bar_p_z_l-bar_p_z_nu-bar_p_z_q1-bar_p_z_q2)/s;
    P_star[26][6]=2*(mod_bar_p_b1*bar_p_x_l/mod_bar_p_l+mod_bar_p_b2*bar_p_x_l/mod_bar_p_l+mod_bar_p_nu*bar_p_x_l/mod_bar_p_l+mod_bar_p_q1*bar_p_x_l/mod_bar_p_l+mod_bar_p_q2*bar_p_x_l/mod_bar_p_l-bar_p_x_b1-bar_p_x_b2-bar_p_x_nu-bar_p_x_q1-bar_p_x_q2)/s;
    P_star[26][7]=2*(mod_bar_p_b1*bar_p_y_l/mod_bar_p_l+mod_bar_p_b2*bar_p_y_l/mod_bar_p_l+mod_bar_p_nu*bar_p_y_l/mod_bar_p_l+mod_bar_p_q1*bar_p_y_l/mod_bar_p_l+mod_bar_p_q2*bar_p_y_l/mod_bar_p_l-bar_p_y_b1-bar_p_y_b2-bar_p_y_nu-bar_p_y_q1-bar_p_y_q2)/s;
    P_star[26][8]=2*(mod_bar_p_b1*bar_p_z_l/mod_bar_p_l+mod_bar_p_b2*bar_p_z_l/mod_bar_p_l+mod_bar_p_nu*bar_p_z_l/mod_bar_p_l+mod_bar_p_q1*bar_p_z_l/mod_bar_p_l+mod_bar_p_q2*bar_p_z_l/mod_bar_p_l-bar_p_z_b1-bar_p_z_b2-bar_p_z_nu-bar_p_z_q1-bar_p_z_q2)/s;
    P_star[26][9]=2*(mod_bar_p_b1*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_b2*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_q1*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_q2*bar_p_x_l/mod_bar_p_nu-bar_p_x_b1-bar_p_x_b2-bar_p_x_l-bar_p_x_q1-bar_p_x_q2)/s;
    P_star[26][10]=2*(mod_bar_p_b1*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_b2*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_q1*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_q2*bar_p_y_l/mod_bar_p_nu-bar_p_y_b1-bar_p_y_b2-bar_p_y_l-bar_p_y_q1-bar_p_y_q2)/s;
    P_star[26][11]=2*(mod_bar_p_b1*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_b2*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_q1*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_q2*bar_p_z_l/mod_bar_p_nu-bar_p_z_b1-bar_p_z_b2-bar_p_z_l-bar_p_z_q1-bar_p_z_q2)/s;
    P_star[26][12]=2*(mod_bar_p_b1*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_b2*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_l*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_nu*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_b1-bar_p_x_b2-bar_p_x_l-bar_p_x_nu-bar_p_x_q2)/s;
    P_star[26][13]=2*(mod_bar_p_b1*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_b2*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_l*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_nu*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_b1-bar_p_y_b2-bar_p_y_l-bar_p_y_nu-bar_p_y_q2)/s;
    P_star[26][14]=2*(mod_bar_p_b1*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_b2*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_l*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_nu*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_b1-bar_p_z_b2-bar_p_z_l-bar_p_z_nu-bar_p_z_q2)/s;
    P_star[26][15]=2*(mod_bar_p_b1*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_b2*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_l*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_nu*bar_p_x_q1/mod_bar_p_q2+mod_bar_p_q2*bar_p_x_q2/mod_bar_p_q2-bar_p_x_b1-bar_p_x_b2-bar_p_x_l-bar_p_x_nu-bar_p_x_q1)/s;
    P_star[26][16]=2*(mod_bar_p_b1*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_b2*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_l*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_nu*bar_p_y_q1/mod_bar_p_q2+mod_bar_p_q2*bar_p_y_q2/mod_bar_p_q2-bar_p_y_b1-bar_p_y_b2-bar_p_y_l-bar_p_y_nu-bar_p_y_q1)/s;
    P_star[26][17]=2*(mod_bar_p_b1*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_b2*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_l*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_nu*bar_p_z_q1/mod_bar_p_q2+mod_bar_p_q2*bar_p_z_q2/mod_bar_p_q2-bar_p_z_b1-bar_p_z_b2-bar_p_z_l-bar_p_z_nu-bar_p_z_q1)/s;
    P_star[26][18]=-bar_xi_2;
    P_star[26][19]=-bar_xi_1;

    if(iteration==0) P_star_i=P_star;
    
    //// Determino il valore di x
    
    x=-1.*((C_star+K_star+P_star).Invert())*0.5*(2.*K_star*f+Q_star);
    
    fittedsolutions=x+f;

    TMatrixD xT(1,27);
    xT.Transpose(x);

    
    chi2=xT*C_star*x;
    
    TGraph gr;
    double chi_min=chi2[0][0];
    double err=0;
    double err1=0;

    //Aggiorno gli errori sulle quantità post fit approssimando il chi2 con una parabola
    
    for(int j=0; j<20; j++){
      for(int i=0; i<3; i++){
	TMatrixD draw(27,1);
	draw.Zero();
	draw[j][0]= fittedsolutions[j][0]+((double) i)*fittedsolutions[j][0]/50.;
	x[j][0]=draw[j][0]-f[j][0];
	xT.Transpose(x);
	chi2=xT*C_star*x;
	gr.SetPoint(i,draw[j][0],chi2[0][0]);
	if(i==2){
	  double x1,x2,x3,y1,y2,y3,a,b,c,sol1,sol2;
	  gr.GetPoint(0,x1,y1);
	  gr.GetPoint(1,x2,y2);
	  gr.GetPoint(2,x3,y3);
	  
	  a=(y3-y1-(y2-y1)*(x3-x1)/(x2-x1))/(x3*x3-x1*x1-(x2*x2-x1*x1)*(x3-x1)/(x2-x1));
	  b=(y2-y1-a*(x2*x2-x1*x1))/(x2-x1);
	  c=y1-a*x1*x1-b*x1;
	  sol1=(-b+sqrt(b*b-4.*a*(c-chi_min-1)))/(2.*a);
	  sol2=(-b-sqrt(b*b-4.*a*(c-chi_min-1)))/(2.*a);
	  err=abs(sol2-sol1)*0.5;
	  if(err!=err) err=0;
	}
      }
    
      
      C[j][j]=abs(err*err);
      if (C[j][j]<2.e-16)  C[j][j]=2.5e-16;
      x[j][0]=fittedsolutions[j][0]-f[j][0];
      err1+=abs(err/fittedsolutions[j][0]);
      //double k=(double) iteration +0.;
      if(j==1 && h==0) {TMatrixD x1(27,1);TMatrixD x1T(1,27);TMatrixD x2T(1,27);TMatrixD x3T(1,27); x1=fittedsolutions-f_iniziale; x1T.Transpose(x1); chi2=x1T*(C_star_i+0.*P_star_i)*x1+0.*x2T.Transpose(fittedsolutions)*K_star*fittedsolutions+0.*x3T.Transpose(Q_star_i)*x1;}
      
      
    }
    err1=0;
 
    chi2[0][0]=chi_min;
    
    //cout << h << endl;
   
    
     
    bar_p_x_b1= fittedsolutions[0][0];   //b1 è adronico
    bar_p_y_b1= fittedsolutions[1][0];
    bar_p_z_b1= fittedsolutions[2][0];
    bar_p_x_b2= fittedsolutions[3][0];   //b2 è leptonico
    bar_p_y_b2= fittedsolutions[4][0];
    bar_p_z_b2= fittedsolutions[5][0];
    bar_p_x_l= fittedsolutions[6][0];
    bar_p_y_l= fittedsolutions[7][0];
    bar_p_z_l= fittedsolutions[8][0];
    bar_p_x_nu= fittedsolutions[9][0];
    bar_p_y_nu= fittedsolutions[10][0];
    bar_p_z_nu= fittedsolutions[11][0];
    bar_p_x_q1= fittedsolutions[12][0];
    bar_p_y_q1= fittedsolutions[13][0];
    bar_p_z_q1= fittedsolutions[14][0];
    bar_p_x_q2= fittedsolutions[15][0];
    bar_p_y_q2= fittedsolutions[16][0];
    bar_p_z_q2= fittedsolutions[17][0];
    bar_xi_1= fittedsolutions[18][0];
    bar_xi_2= fittedsolutions[19][0];
    

    ////////////////////////////////////////////////////##################################
    test_b1.SetPxPyPzE(bar_p_x_b1,bar_p_y_b1,bar_p_z_b1,mod_bar_p_b1);
    test_b2.SetPxPyPzE(bar_p_x_b2,bar_p_y_b2,bar_p_z_b2,mod_bar_p_b2);
    test_l.SetPxPyPzE(bar_p_x_l,bar_p_y_l,bar_p_z_l,mod_bar_p_l);
    test_nu.SetPxPyPzE(bar_p_x_nu,bar_p_y_nu,bar_p_z_nu,mod_bar_p_nu);
    test_q1.SetPxPyPzE(bar_p_x_q1,bar_p_y_q1,bar_p_z_q1,mod_bar_p_q1);
    test_q2.SetPxPyPzE(bar_p_x_q2,bar_p_y_q2,bar_p_z_q2,mod_bar_p_q2);
    ttbar= test_b1+test_b2+test_l+test_nu+test_q1+test_q2;
    lep=test_b2+test_l+test_nu;
    ad=test_b1+test_q1+test_q2;
    lepnu=test_l+test_nu;
    q1q2= test_q1+test_q2;
    TLorentzVector top= test_b2+test_l+test_nu;
    TLorentzVector tbar=test_b1+test_q1+test_q2;
    TLorentzVector ttsystem=top+tbar;
    /////////////////////////////////////////////////////##################################

    double p_x_t;
    double p_x_tbar;
    double p_y_t;
    double p_y_tbar;
    double p_z_t;
    double p_z_tbar;
    double E_t;
    double E_tbar;
    double p_t;
    //double p_tbar;
    //double m_inv_tbar;  
    //double m_lepnu;     
    //double m_q1q2;     
    double m_inv_ttbar;
    //double pt_ttbar; 
    //double pt_tbar;   

     //Costruisco quantità per il calcolo della massa inv
    p_x_tbar=fittedsolutions[0][0]+fittedsolutions[12][0]+fittedsolutions[15][0];
    p_y_t=fittedsolutions[4][0]+fittedsolutions[7][0]+fittedsolutions[10][0];
    p_y_tbar=fittedsolutions[1][0]+fittedsolutions[13][0]+fittedsolutions[16][0];
    p_z_t=fittedsolutions[5][0]+fittedsolutions[8][0]+fittedsolutions[11][0];
    p_z_tbar=fittedsolutions[2][0]+fittedsolutions[14][0]+fittedsolutions[17][0];

    E_t=sqrt(fittedsolutions[3][0]*fittedsolutions[3][0]+fittedsolutions[4][0]*fittedsolutions[4][0]+fittedsolutions[5][0]*fittedsolutions[5][0])+sqrt(fittedsolutions[6][0]*fittedsolutions[6][0]+fittedsolutions[7][0]*fittedsolutions[7][0]+fittedsolutions[8][0]*fittedsolutions[8][0])+sqrt(fittedsolutions[9][0]*fittedsolutions[9][0]+fittedsolutions[10][0]*fittedsolutions[10][0]+fittedsolutions[11][0]*fittedsolutions[11][0]);

    
    E_tbar=sqrt(fittedsolutions[0][0]*fittedsolutions[0][0]+fittedsolutions[1][0]*fittedsolutions[1][0]+fittedsolutions[2][0]*fittedsolutions[2][0])+sqrt(fittedsolutions[12][0]*fittedsolutions[12][0]+fittedsolutions[13][0]*fittedsolutions[13][0]+fittedsolutions[14][0]*fittedsolutions[14][0])+sqrt(fittedsolutions[15][0]*fittedsolutions[15][0]+fittedsolutions[16][0]*fittedsolutions[16][0]+fittedsolutions[17][0]*fittedsolutions[17][0]);

     p_x_t=fittedsolutions[3][0]+fittedsolutions[6][0]+fittedsolutions[9][0];       

    //double p_l=sqrt(fittedsolutions[6][0]*fittedsolutions[6][0]+fittedsolutions[7][0]*fittedsolutions[7][0]+fittedsolutions[8][0]*fittedsolutions[8][0]);
    //double p_nu=sqrt(fittedsolutions[9][0]*fittedsolutions[9][0]+fittedsolutions[10][0]*fittedsolutions[10][0]+fittedsolutions[11][0]*fittedsolutions[11][0]);
    //double p_q1=sqrt(fittedsolutions[12][0]*fittedsolutions[12][0]+fittedsolutions[13][0]*fittedsolutions[13][0]+fittedsolutions[14][0]*fittedsolutions[14][0]);
    //double p_q2=sqrt(fittedsolutions[15][0]*fittedsolutions[15][0]+fittedsolutions[16][0]*fittedsolutions[16][0]+fittedsolutions[17][0]*fittedsolutions[17][0]);

    double p_x_ttbar=p_x_t+p_x_tbar;
    double p_y_ttbar=p_y_t+p_y_tbar;
    double p_z_ttbar=p_z_t+p_z_tbar;
    p_t=sqrt(p_x_t*p_x_t+p_y_t*p_y_t+p_z_t*p_z_t);
    //p_tbar=sqrt(p_x_tbar*p_x_tbar+p_y_tbar*p_y_tbar+p_z_tbar*p_z_tbar);
    double p_ttbar=sqrt(p_x_ttbar*p_x_ttbar+p_y_ttbar*p_y_ttbar+p_z_ttbar*p_z_ttbar);

    
    //double E_ttbar_1 = sqrt(p_t*p_t+173.1*173.1)+sqrt(p_tbar*p_tbar+173.1*173.1);
    double E_ttbar=E_t+E_tbar;
    m_inv_ttbar = sqrt(-p_ttbar*p_ttbar+E_ttbar*E_ttbar);
    //cout << E_ttbar_1-E_ttbar << endl;
    //double pt_t=sqrt(p_x_t*p_x_t+p_y_t*p_y_t);
    //pt_tbar=sqrt(p_x_tbar*p_x_tbar+p_y_tbar*p_y_tbar);
    //pt_ttbar=sqrt(p_x_ttbar*p_x_ttbar+p_y_ttbar*p_y_ttbar);
    

    //m_lepnu = sqrt(2*p_l*p_nu-2*fittedsolutions[6][0]*fittedsolutions[9][0]-2*fittedsolutions[7][0]*fittedsolutions[10][0]-2*fittedsolutions[8][0]*fittedsolutions[11][0]);
    //m_q1q2 = sqrt(2*p_q1*p_q2-2*fittedsolutions[12][0]*fittedsolutions[15][0]-2*fittedsolutions[13][0]*fittedsolutions[16][0]-2*fittedsolutions[14][0]*fittedsolutions[17][0]);

    m_inv_t = sqrt(-p_t*p_t+E_t*E_t);
    //m_inv_tbar = sqrt(-p_tbar*p_tbar+E_tbar*E_tbar);
    //double y=0.5*log((E_ttbar+p_z_ttbar)/(E_ttbar-p_z_ttbar));
    
    //double p_y_t_gen=pt_t_gen*sin(tree->GetLeaf("gen_t_phi")->GetValue(0));
    //double p_y_tbar_gen=pt_tbar_gen*sin(tree->GetLeaf("gen_tbar_phi")->GetValue(0));
    double p_x_b2_gen=tree->GetLeaf("gen_b_pt")->GetValue(0)*cos(tree->GetLeaf("gen_b_phi")->GetValue(0));
    double p_x_b2_bar_gen=tree->GetLeaf("gen_bbar_pt")->GetValue(0)*cos(tree->GetLeaf("gen_bbar_phi")->GetValue(0));
    double eta_ttbar_gen=tree->GetLeaf("gen_ttbar_eta")->GetValue(0);
    
    
    if(iteration==49 && is_h==0 ) histo9.Fill((bar_p_x_b1-p_x_b2_bar_gen));
    if(iteration==49 && is_h==1 ) histo9.Fill((bar_p_x_b1-p_x_b2_gen));			      
    if(iteration==49 && is_h==0) histo10.Fill((pt_b1*cos(phi_b1)-p_x_b2_bar_gen));
    if(iteration==49 && is_h==1) histo10.Fill((pt_b1*cos(phi_b1)-p_x_b2_gen));

    //if(iteration==49 && is_h==1 && intent==0 && no_save!=10) histo9.Fill(bar_p_x_b2-p_x_b2_bar_gen);
    //if(iteration==49 && is_h==0 && intent==0 && no_save!=10) histo9.Fill(bar_p_x_b2-p_x_b2_gen);			      
    //if(iteration==49 && is_h==1 && intent==0 && no_save!=10) histo10.Fill(pt_b2*cos(phi_b2)-p_x_b2_bar_gen);
    //if(iteration==49 && is_h==0 && intent==0 && no_save!=10) histo10.Fill(pt_b2*cos(phi_b2)-p_x_b2_gen);
    
    
    
    if(iteration==49) {x=fittedsolutions-f_iniziale; xT.Transpose(x); chi2=xT*C_star_i*x;}
    if(iteration==49) histo.Fill(lepnu.M());
    if(iteration==49) histo1.Fill(q1q2.M());
    if(iteration==49) histo2.Fill(top.M());
    if(iteration==49) histo3.Fill(tbar.M());
    if(iteration==49) histo4.Fill((ttbar.M()-m_ttbar_gen)/m_ttbar_gen);
    if(iteration==49) histo5.Fill(chi2[0][0]);
    if(iteration==49) histo6.Fill(p_x_t+p_x_tbar);
    if(iteration==49) histo7.Fill(p_y_t+p_y_tbar);
    if(iteration==49) histo8.Fill((m_inv_ttbar*m_inv_ttbar/s-bar_xi_1*bar_xi_2)/(bar_xi_1*bar_xi_2));
    //if(iteration==49 && is_h==0 && intent==0) histo10.Fill(pt_tbar-pt_tbar_gen);
    //if(iteration==49 && is_h==1 && intent==0) histo10.Fill(pt_tbar-pt_t_gen);
    //if(iteration==49 && is_h==0 && intent==0) histo9.Fill((sqrt(p_x_t*p_x_t+p_y_t*p_y_t)-pt_t_gen));
    //if(iteration==49 && is_h==1 && intent==0) histo9.Fill((sqrt(p_x_t*p_x_t+p_y_t*p_y_t)-pt_tbar_gen));
    if(iteration==49) histo11.Fill(sqrt(bar_p_x_nu*bar_p_x_nu+bar_p_y_nu*bar_p_y_nu));
    if(iteration==49) histo12.Fill((ttsystem.Eta()-eta_ttbar_gen));
    if(iteration==49 && chi2[0][0]<500.) profile.SetPoint(h, chi2[0][0], (ttbar.M()-m_ttbar_gen)/m_ttbar_gen );
    if(iteration==49){
      entry1=tbar.M();
      entry2=top.Pt();
      entry3=ttbar.Pt();  
      entry4=test_l.Pt();    
      entry5=test_nu.Pt();
      entry6=chi2[0][0];
      entry7=test_q3.Pt();
      entry8=test_q4.Pt();
      entry9=tree->GetLeaf("nJets")->GetValue(0);
      entry10=top.Eta();
      entry11=tbar.Eta();
      entry12=top.Phi();
      entry13=tbar.Phi();
      entry14=test_b1.Phi();
      entry15=test_b2.Phi();
      entry16=test_b1.Eta();
      entry17=test_b2.Eta();
      entry18=test_q1.E();
      entry19=test_q2.E();
      entry20=test_q3.E();
      entry21=test_q4.E();
      entry22=test_l.E();
      entry23=ttbar.E();
      entry24=test_b1.E();
      entry25=test_b2.E();
      entry26=test_b3.E();
      entry27=test_b4.E();
      entry41=test_b3.Phi();
      entry42=test_b4.Phi();
      entry43=test_b3.Eta();
      entry44=test_b4.Eta();
      entry28=m_b1;
      entry29=m_b2;
      entry30=m_q1;
      entry31=m_q2;
      entry32=m_q3;
      entry33=m_q4;
      entry34=m_b3;
      entry35=m_b4;
      entry36=ttbar.M();
      entry37=sqrt(bar_p_x_nu*bar_p_x_nu+bar_p_y_nu*bar_p_y_nu);
      entry38=sqrt(bar_p_x_l*bar_p_x_l+bar_p_y_l*bar_p_y_l);
      entry39=tree->GetLeaf("nLightJets")->GetValue(0);
      entry40=tree->GetLeaf("lepton_isolation")->GetValue(0);
      vector<double> ang_dist;
      ang_dist.push_back(test_q1.DeltaR(test_q2));
      if(test_q3.Pt()!=0 && test_q4.Pt()!=0) {tot_dist=(test_q1.DeltaR(test_q2)+test_q1.DeltaR(test_q3)+test_q1.DeltaR(test_q4)+test_q2.DeltaR(test_q3)+test_q2.DeltaR(test_q4)+test_q3.DeltaR(test_q4))/6.; ang_dist.push_back(test_q1.DeltaR(test_q3)); ang_dist.push_back(test_q1.DeltaR(test_q4)); ang_dist.push_back(test_q2.DeltaR(test_q3)); ang_dist.push_back(test_q2.DeltaR(test_q4));ang_dist.push_back(test_q3.DeltaR(test_q4));}
      if(test_q3.Pt()==0 && test_q4.Pt()!=0)  {tot_dist=(test_q1.DeltaR(test_q2)+test_q1.DeltaR(test_q4)+test_q2.DeltaR(test_q4))/3.; ang_dist.push_back(test_q1.DeltaR(test_q4)); ang_dist.push_back(test_q2.DeltaR(test_q4));}
      if(test_q3.Pt()!=0 && test_q4.Pt()==0)  {tot_dist=(test_q1.DeltaR(test_q2)+test_q1.DeltaR(test_q3)+test_q2.DeltaR(test_q3))/3.; ang_dist.push_back(test_q1.DeltaR(test_q3)); ang_dist.push_back(test_q2.DeltaR(test_q3));}
      if(test_q3.Pt()==0 && test_q4.Pt()==0)  tot_dist=test_q1.DeltaR(test_q2);
      entry45=tot_dist;
      entry46=ang_dist[dump_index(ang_dist,ang_dist.size())];
      entry47=tree->GetLeaf("weight")->GetValue(0);
      entry48=tree->GetLeaf("p1_xi")->GetValue(0);
      entry49=tree->GetLeaf("p2_xi")->GetValue(0);
      tree_out.Fill();
    }
    
    }
  }

  cout << "occhio vale: " << occhio << endl; 
  TCanvas c0;
  histo.Draw("E");
  TCanvas c1;
  histo1.Draw("E");
  TCanvas c2;
  histo2.Draw("E");
  TCanvas c3;
  histo3.Draw("E");
  TCanvas c4;
  histo4.Draw("E");
  TCanvas c5;
  histo5.Draw("E");
  TCanvas c6;
  histo6.Draw("E");
  TCanvas c7;
  histo7.Draw("E");
  TCanvas c8;
  histo8.Draw("E");
  TCanvas c9;
  histo9.Draw("E");
  TCanvas c10;
  histo10.Draw("E");
  TCanvas c11;
  histo11.Draw("E");
  TCanvas c12;
  histo12.Draw("E");
  TCanvas c13;
  profile.Draw("AP");
  TFile ciao("histo_data_final.root","RECREATE");
  ciao.cd();
  histo.Write();
  histo1.Write();
  histo2.Write();
  histo3.Write();
  histo4.Write();
  histo5.Write();
  histo6.Write();
  histo7.Write();
  histo8.Write();
  histo9.Write();
  histo10.Write();
  histo11.Write();
  histo12.Write();
  ciao.Close();
  output.cd();
  if(evt_count) evt_count->Write();
  cout << "Writes " << output.GetName() << endl;
  output.Write();
  output.Close();
  //app.Run("true");

  //////
  return 0;
}
