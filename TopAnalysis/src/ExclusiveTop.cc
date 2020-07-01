#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
//#include "TLorentzVector.h"  // !!!!!!!!!!!! IT SHOULD GO WITH ", NOT <
#include <TGraphAsymmErrors.h>
#include "TMatrixD.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TApplication.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveTop.h"
#include "TopLJets2015/TopAnalysis/interface/NeutrinoEllipseCalculator.h"
//#include "TopLJets2015/CTPPSAnalysisTools/interface/LHCConditionsFactory.h"
//#include "TopLJets2015/CTPPSAnalysisTools/interface/AlignmentsFactory.h"
//#include "TopLJets2015/CTPPSAnalysisTools/interface/XiReconstructor.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

/*
Modifications form the original code:
0. copy neutrino calculator from https://github.com/betchart/analytic-nu 
1. comment CTPPSAnalysisTools dependences, and comment strip reconstruction
2. change met_XX to float (not vector)
3. BTagSFUtil btvSF() different from the branch
4. attachToMiniEventTree() remove systematics
5. ADDVAR(&(outVars[fvars[i]]),fvars[i],"/F",outT); // with /F instead of F
6. rename input tree to "analysis/tree"
*/

using namespace std;
#define ADDVAR(x,name,t,tree) tree->Branch(name,x,TString(name)+TString(t))

double m_TOP = 173.1;
double m_W   =  80.379;
double m_NU  =  0.;

// ---- SWITCHES ----------------------------------------------------------------
#define HISTOGRAMS_ON           // comment to avoid creating histograms in root file
#define SAVEERRORS_ON           // comment to avoid saving the quantities error in root file
//#define KINFIT_ON             // comment to avoid calculating KinFit
//#define MattFIT_ON            // comment to avoid calculating MattFit
//#define analyticMattFit_ON    // comment to avoid calculating analyticMattFit


#ifdef MattFIT_ON
    // b1 is the hadronic b, b2 is the leptonic one
    double phi_b1, phi_b2, phi_l, phi_nu, phi_q1, phi_q2;
    double eta_b1, eta_b2, eta_l, eta_nu, eta_q1, eta_q2;
    double pt_b1, pt_b2, pt_l, pt_nu, pt_q1, pt_q2;
    double m_l     = 0.;
    double m_b     = 0.;
    double m_q1    = 0.;
    double m_q2    = 0.;
    double m_t     = m_TOP;
    double m_tbar  = m_TOP;
    double m_W_ad  = m_W;
    double m_W_lep = m_W;
    double lambda_Value = 1000.;                 // 0.01 = initial values of lambdas
    double f_value = 0.; // value of f (chistar of MattFit). Global because it is computed inside chistar and then read outside

    double e_b1, e_b2, e_l, e_nu, e_q1, e_q2;   //error on pt
    double e_t, e_tbar, e_W_ad, e_W_lep;        //error on masses
    double e_phi_b1, e_phi_b2, e_phi_l, e_phi_nu, e_phi_q1, e_phi_q2;
    double e_eta_b1, e_eta_b2, e_eta_l, e_eta_nu, e_eta_q1, e_eta_q2;

    TF1 H1("H1","[0]+[1]+[2]+[3]+[4]+[5]",-1.e9,1.e9); //somma impulsi asse x
    TF1 H2("H2","[0]+[1]+[2]+[3]+[4]+[5]",-1.e9,1.e9); //somma impulsi asse y
    TF1 theta("theta", "2*atan(exp(-x))",0, 1e10);
#endif

TF1 theta("theta", "2*atan(exp(-x))",0, 1e10);
// Returns the index of the lesser element of sum_squared vector
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

// PARTICLE_err returns the uncertainty on the particle: p (if err_type==0), eta (if err_type==1), phi (if err_type==2).
// The uncertainties are extracted from tables
double bjet_err(double p,double eta,double phi, int err_type){
 eta=abs(eta);

 double c_p,r_p,n_p,c_e,r_e,n_e,c_ph,r_ph,n_ph;

 if(eta>=0 && eta<0.0870)
 {c_p=0.08790000;r_p=0.90500000; n_p=2.93000000;c_e=1.60700000;r_e=0.00000000; n_e=0.00611000;c_ph=1.67550000;r_ph=0.00000000; n_ph=0.00849000;}
 if(eta>=0.0870 && eta<0.1740)
 {c_p=0.08460000;r_p=0.96600000; n_p=2.14000000;c_e=1.62100000;r_e=0.00000000; n_e=0.00554;c_ph=1.70500000;r_ph=0.00000000; n_ph=0.00791000;}
 if(eta>=0.1740 && eta<0.2610)
 {c_p=0.07890000;r_p=0.99900000; n_p=1.74000000;c_e=1.61160000;r_e=0.00000000; n_e=0.00618;c_ph=1.70810000;r_ph=0.00000000; n_ph=0.0078;}
 if(eta>=0.2610 && eta<0.3480)
 {c_p=0.08150000;r_p=0.96100000; n_p=2.08000000;c_e=1.62740000;r_e=0.00000000; n_e=0.00577;c_ph=1.70130000;r_ph=0.00000000; n_ph=0.00766;}
 if(eta>=0.3480 && eta<0.4350)
 {c_p=0.07630000;r_p=1.02400000; n_p=1.21000000;c_e=1.63090000;r_e=0.00000000; n_e=0.00587;c_ph=1.69770000;r_ph=0.00000000; n_ph=0.00784;}
 if(eta>=0.4350 && eta<0.5220)
 {c_p=0.07600000;r_p=1.00900000; n_p=1.65000000;c_e=1.63880000;r_e=0.00000000; n_e=0.00648;c_ph=1.71110000;r_ph=0.00000000; n_ph=0.00787;}
 if(eta>=0.5220 && eta<0.6090)
 {c_p=0.06880000;r_p=1.05500000; n_p=0.92000000;c_e=1.64200000;r_e=0.00000000; n_e=0.00607;c_ph=1.70990000;r_ph=0.00000000; n_ph=0.00762;}
 if(eta>=0.6090 && eta<0.6960)
 {c_p=0.07200000;r_p=1.05170000; n_p=0.00000000;c_e=1.65260000;r_e=0.00000000; n_e=0.00656;c_ph=1.72850000;r_ph=0.00000000; n_ph=0.00750;}
 if(eta>=0.6960 && eta<0.7830)
 {c_p=0.07540000;r_p=1.06800000; n_p=0.00000000;c_e=1.65360000;r_e=0.00000000; n_e=0.00687;c_ph=1.71680000;r_ph=0.00000000; n_ph=0.00779;}
 if(eta>=0.7830 && eta<0.8700)
 {c_p=0.07510000;r_p=1.07970000; n_p=0.00000000;c_e=1.65300000;r_e=0.00000000; n_e=0.00764;c_ph=1.74680000;r_ph=0.00000000; n_ph=0.00744;}
 if(eta>=0.8700 && eta<0.9570)
 {c_p=0.07980000;r_p=1.07840000; n_p=0.00000000;c_e=1.66150000;r_e=0.00000000; n_e=0.00782;c_ph=1.75840000;r_ph=0.00000000; n_ph=0.00758;}
 if(eta>=0.9570 && eta<1.0440)
 {c_p=0.07430000;r_p=1.11190000; n_p=0.00000000;c_e=1.67650000;r_e=0.00000000; n_e=0.00785;c_ph=1.77180000;r_ph=0.00000000; n_ph=0.00798;}
 if(eta>=1.0440 && eta<1.1310)
 {c_p=0.07570000;r_p=1.14090000; n_p=0.00000000;c_e=1.68270000;r_e=0.00000000; n_e=0.00816;c_ph=1.79710000;r_ph=0.00000000; n_ph=0.00812;}
 if(eta>=1.1310 && eta<1.2180)
 {c_p=0.07550000;r_p=1.15240000; n_p=0.00000000;c_e=1.73900000;r_e=0.00000000; n_e=0.00810;c_ph=1.81990000;r_ph=0.00000000; n_ph=0.00846;}
 if(eta>=1.2180 && eta<1.3050)
 {c_p=0.08250000;r_p=1.15410000; n_p=0.00000000;c_e=1.74910000;r_e=0.00000000; n_e=0.01052;c_ph=1.86830000;r_ph=0.00000000; n_ph=0.00910;}
 if(eta>=1.3050 && eta<1.3920)
 {c_p=0.08770000;r_p=1.18100000; n_p=0.00000000;c_e=1.74100000;r_e=0.00000000; n_e=0.01346;c_ph=1.91170000;r_ph=0.00000000; n_ph=0.01070;}
 if(eta>=1.3920 && eta<1.4790)
 {c_p=0.10310000;r_p=1.14400000; n_p=0.00000000;c_e=1.77200000;r_e=0.00000000; n_e=0.01026;c_ph=1.99920000;r_ph=0.00000000; n_ph=0.01069;}
 if(eta>=1.4790 && eta<1.5660)
 {c_p=0.10260000;r_p=1.15100000; n_p=0.00000000;c_e=1.78600000;r_e=0.00000000; n_e=0.01088;c_ph=1.99100000;r_ph=0.00000000; n_ph=0.00985;}
 if(eta>=1.5660 && eta<1.6530)
 {c_p=0.09090000;r_p=1.16300000; n_p=0.00000000;c_e=1.82400000;r_e=0.00000000; n_e=0.01008;c_ph=1.97000000;r_ph=0.00000000; n_ph=0.00874;}
 if(eta>=1.6530 && eta<1.7400)
 {c_p=0.08850000;r_p=1.15600000; n_p=0.00000000;c_e=1.84900000;r_e=0.00000000; n_e=0.00872;c_ph=1.57000000;r_ph=0.00000000; n_ph=0.00761;}
 if(eta>=1.7400 && eta<1.8300)
 {c_p=0.08930000;r_p=1.12800000; n_p=0.00000000;c_e=1.91200000;r_e=0.00000000; n_e=0.00922;c_ph=1.93300000;r_ph=0.00000000; n_ph=0.00690;}
 if(eta>=1.8300 && eta<1.9300)
 {c_p=0.07960000;r_p=1.13400000; n_p=0.00000000;c_e=1.88400000;r_e=0.00000000; n_e=0.01102;c_ph=1.89500000;r_ph=0.00000000; n_ph=0.00665;}
 if(eta>=1.9300 && eta<2.0430)
 {c_p=0.06530000;r_p=1.13630000; n_p=0.00000000;c_e=1.89500000;r_e=0.00000000; n_e=0.01168;c_ph=1.85700000;r_ph=0.00000000; n_ph=0.00646;}
 if(eta>=2.0430 && eta<2.1720)
 {c_p=0.06240000;r_p=1.13800000; n_p=0.00000000;c_e=1.90200000;r_e=0.00000000; n_e=0.01183;c_ph=1.87300000;r_ph=0.00000000; n_ph=0.00647;}
 if(eta>=2.1720 && eta<2.3220)
 {c_p=0.06530000;r_p=1.12400000; n_p=0.00000000;c_e=1.91100000;r_e=0.00000000; n_e=0.01127;c_ph=1.87300000;r_ph=0.00000000; n_ph=0.00540;}
 if(eta>=2.3220)
 {c_p=0.05040000;r_p=1.20200000; n_p=0.50000000;c_e=2.28000000;r_e=0.00000000; n_e=0.01053;c_ph=1.98800000;r_ph=0.00000000; n_ph=0.00887;}

 double theta=2*atan(exp(-eta));

 if(err_type==0){return sqrt(c_p*c_p*p*p*sin(theta)*sin(theta)+r_p*r_p*p*sin(theta)+n_p*n_p);}
 else if (err_type==1){return sqrt(c_e*c_e/(p*p)+r_e*r_e/p+n_e*n_e);}
 else if (err_type==2){return sqrt(c_ph*c_ph/(p*p)+r_ph*r_ph/p+n_ph*n_ph);}
 else { return 0; } // THIS WILL CAUSE THE FIT TO CRASH --> err_type must be 0,1,2
}

double ljet_err(double p,double eta,double phi, int err_type){
 eta=abs(eta);

 double c_p,r_p,n_p,c_e,r_e,n_e,c_ph,r_ph,n_ph;

 if(eta>=0 && eta<0.0870)
 {c_p=0.07390000;r_p=0.90400000; n_p=2.39000000;c_e=1.48700000;r_e=0.00000000; n_e=0.00986000;c_ph=1.54610000;r_ph=0.00000000; n_ph=0.01155000;}
 if(eta>=0.0870 && eta<0.1740)
 {c_p=0.07390000;r_p=0.92500000; n_p=1.90000000;c_e=1.45050000;r_e=0.00000000; n_e=0.01024;c_ph=1.48680000;r_ph=0.00000000; n_ph=0.01224000;}
 if(eta>=0.1740 && eta<0.2610)
 {c_p=0.06720000;r_p=0.94300000; n_p=2.02000000;c_e=1.44220000;r_e=0.00000000; n_e=0.01047;c_ph=1.50530000;r_ph=0.00000000; n_ph=0.01159;}
 if(eta>=0.2610 && eta<0.3480)
 {c_p=0.06570000;r_p=0.97100000; n_p=1.73000000;c_e=1.48270000;r_e=0.00000000; n_e=0.01009;c_ph=1.49650000;r_ph=0.00000000; n_ph=0.01204;}
 if(eta>=0.3480 && eta<0.4350)
 {c_p=0.06690000;r_p=0.94300000; n_p=2.09000000;c_e=1.47810000;r_e=0.00000000; n_e=0.01026;c_ph=1.51630000;r_ph=0.00000000; n_ph=0.01184;}
 if(eta>=0.4350 && eta<0.5220)
 {c_p=0.06690000;r_p=0.95200000; n_p=1.89000000;c_e=1.47910000;r_e=0.00000000; n_e=0.01065;c_ph=1.50480000;r_ph=0.00000000; n_ph=0.01221;}
 if(eta>=0.5220 && eta<0.6090)
 {c_p=0.06980000;r_p=0.93500000; n_p=2.02000000;c_e=1.46160000;r_e=0.00000000; n_e=0.01103;c_ph=1.50120000;r_ph=0.00000000; n_ph=0.01193;}
 if(eta>=0.6090 && eta<0.6960)
 {c_p=0.06340000;r_p=0.98100000; n_p=1.80000000;c_e=1.47420000;r_e=0.00000000; n_e=0.01095;c_ph=1.52640000;r_ph=0.00000000; n_ph=0.01159;}
 if(eta>=0.6960 && eta<0.7830)
 {c_p=0.05730000;r_p=1.03000000; n_p=1.11000000;c_e=1.46330000;r_e=0.00000000; n_e=0.01123;c_ph=1.51400000;r_ph=0.00000000; n_ph=0.01193;}
 if(eta>=0.7830 && eta<0.8700)
 {c_p=0.06720000;r_p=0.9920000; n_p=1.73000000;c_e=1.48330000;r_e=0.00000000; n_e=0.01149;c_ph=1.54820000;r_ph=0.00000000; n_ph=0.01122;}
 if(eta>=0.8700 && eta<0.9570)
 {c_p=0.06450000;r_p=1.06100000; n_p=0.9300000;c_e=1.5212000;r_e=0.00000000; n_e=0.01132;c_ph=1.56030000;r_ph=0.00000000; n_ph=0.01175;}
 if(eta>=0.9570 && eta<1.0440)
 {c_p=0.06170000;r_p=1.07700000; n_p=1.44000000;c_e=1.48860000;r_e=0.00000000; n_e=0.01218;c_ph=1.57330000;r_ph=0.00000000; n_ph=0.01135;}
 if(eta>=1.0440 && eta<1.1310)
 {c_p=0.06420000;r_p=1.11000000; n_p=0.72000000;c_e=1.48070000;r_e=0.00000000; n_e=0.01239;c_ph=1.57390000;r_ph=0.00000000; n_ph=0.01221;}
 if(eta>=1.1310 && eta<1.2180)
 {c_p=0.05980000;r_p=1.15800000; n_p=0.54000000;c_e=1.52920000;r_e=0.00000000; n_e=0.01204;c_ph=1.61650000;r_ph=0.00000000; n_ph=0.01161;}
 if(eta>=1.2180 && eta<1.3050)
 {c_p=0.05320000;r_p=1.21410000; n_p=0.00000000;c_e=1.52750000;r_e=0.00000000; n_e=0.01453;c_ph=1.64260000;r_ph=0.00000000; n_ph=0.01244;}
 if(eta>=1.3050 && eta<1.3920)
 {c_p=0.05830000;r_p=1.24180000; n_p=0.00000000;c_e=1.47370000;r_e=0.00000000; n_e=0.01841;c_ph=1.68600000;r_ph=0.00000000; n_ph=0.01441;}
 if(eta>=1.3920 && eta<1.4790)
 {c_p=0.06730000;r_p=1.25480000; n_p=0.00000000;c_e=1.52920000;r_e=0.00000000; n_e=0.01432;c_ph=1.71230000;r_ph=0.00000000; n_ph=0.01507;}
 if(eta>=1.4790 && eta<1.5660)
 {c_p=0.05920000;r_p=1.25950000; n_p=0.00000000;c_e=1.57180000;r_e=0.00000000; n_e=0.01358;c_ph=1.73310000;r_ph=0.00000000; n_ph=0.01307;}
 if(eta>=1.5660 && eta<1.6530)
 {c_p=0.03990000;r_p=1.23260000; n_p=0.00000000;c_e=1.56400000;r_e=0.00000000; n_e=0.01267;c_ph=1.73490000;r_ph=0.00000000; n_ph=0.01179;}
 if(eta>=1.6530 && eta<1.7400)
 {c_p=0.02820000;r_p=1.21750000; n_p=0.00000000;c_e=1.61360000;r_e=0.00000000; n_e=0.01261;c_ph=1.71360000;r_ph=0.00000000; n_ph=0.01118;}
 if(eta>=1.7400 && eta<1.8300)
 {c_p=0.03860000;r_p=1.17960000; n_p=0.00000000;c_e=1.67590000;r_e=0.00000000; n_e=0.01174;c_ph=1.68610000;r_ph=0.00000000; n_ph=0.01100;}
 if(eta>=1.8300 && eta<1.9300)
 {c_p=0.03810000;r_p=1.13200000; n_p=0.00000000;c_e=1.63300000;r_e=0.00000000; n_e=0.01365;c_ph=1.67390000;r_ph=0.00000000; n_ph=0.01044;}
 if(eta>=1.9300 && eta<2.0430)
 {c_p=0.0000000;r_p=1.12480000; n_p=0.00000000;c_e=1.65090000;r_e=0.00000000; n_e=0.01340;c_ph=1.680100000;r_ph=0.00000000; n_ph=0.01056;}
 if(eta>=2.0430 && eta<2.1720)
 {c_p=0.00000000;r_p=1.11740000; n_p=0.0000000;c_e=1.63300000;r_e=0.00000000; n_e=0.01481;c_ph=1.65640000;r_ph=0.00000000; n_ph=0.01084;}
 if(eta>=2.1720 && eta<2.3220)
 {c_p=0.00000000;r_p=1.08890000; n_p=1.51040000;c_e=1.63500000;r_e=0.00000000; n_e=0.01430;c_ph=1.66830000;r_ph=0.00000000; n_ph=0.01113;}
 if(eta>=2.3220)
 {c_p=0.00000000;r_p=1.06820000; n_p=2.65000000;c_e=2.04200000;r_e=0.00000000; n_e=0.01414;c_ph=1.75290000;r_ph=0.00000000; n_ph=0.01311;}

 double theta=2*atan(exp(-eta));


 if(err_type==0){return sqrt(c_p*c_p*p*p*sin(theta)*sin(theta)+r_p*r_p*p*sin(theta)+n_p*n_p);}
 else if(err_type==1){return sqrt(c_e*c_e/(p*p)+r_e*r_e/p+n_e*n_e);}
 else if(err_type==2){return sqrt(c_ph*c_ph/(p*p)+r_ph*r_ph/p+n_ph*n_ph);}
 else { return 0; } // THIS WILL CAUSE THE FIT TO CRASH --> err_type must be 0,1,2
}

double mu_err(double p,double eta,double phi, int err_type){
 eta=abs(eta);

 double c_p,r_p,n_p,c_e,r_e,n_e,c_ph,r_ph,n_ph;

 if(eta>=0 && eta<0.10)
 {c_p=0.00530000;r_p=0.0012329; n_p=2.0001434;c_e=0.00305000;r_e=0.00000000; n_e=0.00043240;c_ph=0.00304000;r_ph=0.00018100; n_ph=0.0000719;}
 if(eta>=0.10 && eta<0.20)
 {c_p=0.00553000;r_p=0.00122959; n_p=0.00013670;c_e=0.00350000;r_e=0.00000000; n_e=0.00038330;c_ph=0.0029900;r_ph=0.00015000; n_ph=0.00007190;}
 if(eta>=0.20 && eta<0.30)
 {c_p=0.00609000;r_p=0.00126893; n_p=0.00013220;c_e=0.00269000;r_e=0.00000000; n_e=0.00033830;c_ph=0.003040000;r_ph=0.00023400; n_ph=0.00006820;}
 if(eta>=0.30 && eta<0.40)
 {c_p=0.00687000;r_p=0.00132409; n_p=0.00012760;c_e=0.00273000;r_e=0.00000000; n_e=0.00030640;c_ph=0.00313000;r_ph=0.00024200; n_ph=0.00006700;}
 if(eta>=0.40 && eta<0.50)
 {c_p=0.00699000;r_p=0.00136460; n_p=0.00013320;c_e=0.00281000;r_e=0.00000000; n_e=0.00028470;c_ph=0.003310000;r_ph=0.00022300; n_ph=0.00006760;}
 if(eta>=0.50 && eta<0.60)
 {c_p=0.00742000;r_p=0.00138092; n_p=0.00012850;c_e=0.00242000;r_e=0.00000000; n_e=0.00028330;c_ph=0.003400000;r_ph=0.00018000; n_ph=0.00006940;}
 if(eta>=0.60 && eta<0.70)
 {c_p=0.00788000;r_p=0.00139850; n_p=0.00012410;c_e=0.00276000;r_e=0.00000000; n_e=0.00028300;c_ph=0.00344000;r_ph=0.000268000; n_ph=0.00006640;}
 if(eta>=0.70 && eta<0.80)
 {c_p=0.00832000;r_p=0.00143528; n_p=0.00012380;c_e=0.00326000;r_e=0.00000000; n_e=0.00027360;c_ph=0.00350000;r_ph=0.00022000; n_ph=0.00007160;}
 if(eta>=0.80 && eta<0.90)
 {c_p=0.00930000;r_p=0.00152113; n_p=0.00012440;c_e=0.00325000;r_e=0.00000000; n_e=0.00025360;c_ph=0.00328000;r_ph=0.00044500; n_ph=0.00006070;}
 if(eta>=0.90 && eta<1.00)
 {c_p=0.01098000;r_p=0.00180097; n_p=0.00014770;c_e=0.00326000;r_e=0.00000000; n_e=0.00024060;c_ph=0.00395000;r_ph=0.00030500; n_ph=0.00007610;}
 if(eta>=1.0000 && eta<1.100)
 {c_p=0.01188000;r_p=0.00184714; n_p=0.00014360;c_e=0.00387000;r_e=0.00000000; n_e=0.00022460;c_ph=0.00411000;r_ph=0.00027000; n_ph=0.00008110;}
 if(eta>=1.10 && eta<1.20)
 {c_p=0.01375000;r_p=0.00184140; n_p=0.00012330;c_e=0.00390000;r_e=0.00000000; n_e=0.00021530;c_ph=0.00434000;r_ph=0.00042400; n_ph=0.00006740;}
 if(eta>=1.20 && eta<1.30)
 {c_p=0.01432000;r_p=0.00197793; n_p=0.00013660;c_e=0.00390000;r_e=0.00000000; n_e=0.00020530;c_ph=0.00378000;r_ph=0.00066300; n_ph=0.00005580;}
 if(eta>=1.30 && eta<1.40)
 {c_p=0.01465000;r_p=0.00208521; n_p=0.00014840;c_e=0.00389000;r_e=0.00000000; n_e=0.0002131;c_ph=0.00397000;r_ph=0.00059300; n_ph=0.00007000;}
 if(eta>=1.40 && eta<1.50)
 {c_p=0.01419000;r_p=0.00207764; n_p=0.00015210;c_e=0.00403000;r_e=0.00000000; n_e=0.00022100;c_ph=0.00411000;r_ph=0.00074000; n_ph=0.00004600;}
 if(eta>=1.50 && eta<1.60)
 {c_p=0.01334000;r_p=0.00207255; n_p=0.00016100;c_e=0.00400000;r_e=0.00000000; n_e=0.00021970;c_ph=0.00376000;r_ph=0.00082600; n_ph=0.00004900;}
 if(eta>=1.60 && eta<1.70)
 {c_p=0.01304000;r_p=0.00222486; n_p=0.00018980;c_e=0.00404000;r_e=0.00000000; n_e=0.0002242;c_ph=0.00379000;r_ph=0.00083800; n_ph=0.00006660;}
 if(eta>=1.70 && eta<1.80)
 {c_p=0.01277000;r_p=0.00263811; n_p=0.00027250;c_e=0.00396000;r_e=0.00000000; n_e=0.00025370;c_ph=0.00429000;r_ph=0.00083000; n_ph=0.00008600;}
 if(eta>=1.80 && eta<1.90)
 {c_p=0.01489000;r_p=0.00313961; n_p=0.00033100;c_e=0.00417000;r_e=0.00000000; n_e=0.00027520;c_ph=0.00419000;r_ph=0.00084000; n_ph=0.00011500;}
 if(eta>=1.90 && eta<2.00)
 {c_p=0.01473000;r_p=0.00366521; n_p=0.00045600;c_e=0.00437000;r_e=0.00000000; n_e=0.00029680;c_ph=0.003920000;r_ph=0.00104000; n_ph=0.00011900;}
 if(eta>=2.00 && eta<2.1000)
 {c_p=0.01667000;r_p=0.00440120; n_p=0.00058100;c_e=0.00429000;r_e=0.00000000; n_e=0.00031940;c_ph=0.00405000;r_ph=0.00110000; n_ph=0.00014900;}
 if(eta>=2.100 && eta<2.200)
 {c_p=0.01641000;r_p=0.00495142; n_p=0.00074700;c_e=0.00479000;r_e=0.00000000; n_e=0.00037380;c_ph=0.0053900;r_ph=0.00087000; n_ph=0.00019500;}
 if(eta>=2.200 && eta<2.30)
 {c_p=0.01710000;r_p=0.00568498; n_p=0.00094500;c_e=0.00526000;r_e=0.00000000; n_e=0.00042750;c_ph=0.00470000;r_ph=0.00129000; n_ph=0.00018900;}
 if(eta>=2.30)
 {c_p=0.01590000;r_p=0.00643950; n_p=0.00130400;c_e=0.00449000;r_e=0.00000000; n_e=0.00052200;c_ph=0.00370000;r_ph=0.00176000; n_ph=0.00019800;}

 double theta=2*atan(exp(-eta));

 if(err_type==0){return p*p*sin(theta)*sin(theta)*sqrt(c_p*c_p/(p*p*sin(theta)*sin(theta))+r_p*r_p/(p*sin(theta))+n_p*n_p);}
 else if(err_type==1){return sqrt(c_e*c_e/(p*p)+r_e*r_e/p+n_e*n_e);}
 else if(err_type==2){return sqrt(c_ph*c_ph/(p*p)+r_ph*r_ph/p+n_ph*n_ph);}
 else { return 0; } // THIS WILL CAUSE THE FIT TO CRASH --> err_type must be 0,1,2
}

double e_err(double p,double eta,double phi, int err_type){
 eta=abs(eta);

 double c_p,r_p,n_p,c_e,r_e,n_e,c_ph,r_ph,n_ph;

 if(eta>=0 && eta<0.1740)
 {c_p=0.00635000;r_p=0.07110000; n_p=0.26800000;c_e=0.00460000;r_e=0.00000000; n_e=0.00045350;c_ph=0.00311000;r_ph=0.00120600; n_ph=0.00009290;}
 if(eta>=0.1740 && eta<0.2610)
 {c_p=0.00330000;r_p=0.08990000; n_p=0.12600000;c_e=0.00381000;r_e=0.00000000; n_e=0.00038600;c_ph=0.00362000;r_ph=0.00120300; n_ph=0.00008700;}
 if(eta>=0.2610 && eta<0.3480)
 {c_p=0.00566000;r_p=0.07850000; n_p=0.22700000;c_e=0.00327000;r_e=0.00000000; n_e=0.00035180;c_ph=0.00514000;r_ph=0.00091000; n_ph=0.00012590;}
 if(eta>=0.3480 && eta<0.4350)
 {c_p=0.00519000;r_p=0.08990000; n_p=0.19800000;c_e=0.00232000;r_e=0.00000000; n_e=0.00033220;c_ph=0.00368000;r_ph=0.00127600; n_ph=0.00010300;}
 if(eta>=0.4350 && eta<0.5220)
 {c_p=0.00448000;r_p=0.08830000; n_p=0.27800000;c_e=0.00267000;r_e=0.00000000; n_e=0.00031120;c_ph=0.00355000;r_ph=0.00128800; n_ph=0.00009500;}
 if(eta>=0.5220 && eta<0.6090)
 {c_p=0.00360000;r_p=0.09150000; n_p=0.18800000;c_e=0.002660;r_e=0.00000000; n_e=0.00031790;c_ph=0.00000000;r_ph=0.00159900; n_ph=0.00000000;}
 if(eta>=0.6090 && eta<0.6960)
 {c_p=0.00370000;r_p=0.09020000; n_p=0.23600000;c_e=0.00295000;r_e=0.00000000; n_e=0.00032990;c_ph=0.00418000;r_ph=0.00130200; n_ph=0.00010700;}
 if(eta>=0.6960 && eta<0.7830)
 {c_p=0.00200000;r_p=0.09860000; n_p=0.24800000;c_e=0.00312000;r_e=0.00000000; n_e=0.00033050;c_ph=0.00456000;r_ph=0.00116000; n_ph=0.00014610;}
 if(eta>=0.7830 && eta<0.8700)
 {c_p=0.00320000;r_p=0.10460000; n_p=0.28600000;c_e=0.00392000;r_e=0.00000000; n_e=0.00030260;c_ph=0.00280000;r_ph=0.00156000; n_ph=0.00012400;}
 if(eta>=0.8700 && eta<0.9570)
 {c_p=0.00000000;r_p=0.10680000; n_p=0.30600000;c_e=0.00405000;r_e=0.00000000; n_e=0.00028930;c_ph=0.00220000;r_ph=0.00174000; n_ph=0.00012900;}
 if(eta>=0.9570 && eta<1.0440)
 {c_p=0.00000000;r_p=0.10910000; n_p=0.48700000;c_e=0.00369000;r_e=0.00000000; n_e=0.00028420;c_ph=0.00400000;r_ph=0.00173000; n_ph=0.00014800;}
 if(eta>=1.0440 && eta<1.1310)
 {c_p=0.00000000;r_p=0.13000000; n_p=0.58700000;c_e=0.00384000;r_e=0.00000000; n_e=0.00029050;c_ph=0.00000000;r_ph=0.00204600; n_ph=0.00015000;}
 if(eta>=1.1310 && eta<1.2180)
 {c_p=0.00000000;r_p=0.14650000; n_p=0.81900000;c_e=0.00474000;r_e=0.00000000; n_e=0.00027940;c_ph=0.00000000;r_ph=0.00209400; n_ph=0.00018100;}
 if(eta>=1.2180 && eta<1.3050)
 {c_p=0.00000000;r_p=0.13800000; n_p=0.90300000;c_e=0.00444000;r_e=0.00000000; n_e=0.00026780;c_ph=0.00400000;r_ph=0.00188000; n_ph=0.00024000;}
 if(eta>=1.3050 && eta<1.3920)
 {c_p=0.00000000;r_p=0.14290000; n_p=0.98800000;c_e=0.00408000;r_e=0.00000000; n_e=0.00029510;c_ph=0.00490000;r_ph=0.00209000; n_ph=0.00024400;}
 if(eta>=1.3920 && eta<1.4790)
 {c_p=0.00000000;r_p=0.20160000; n_p=0.89000000;c_e=0.00400000;r_e=0.00000000; n_e=0.00029060;c_ph=0.00400000;r_ph=0.00249000; n_ph=0.00023100;}
 if(eta>=1.4790 && eta<1.6530)
 {c_p=0.02770000;r_p=0.24400000; n_p=0.66000000;c_e=0.00451000;r_e=0.00000000; n_e=0.00027320;c_ph=0.00360000;r_ph=0.00262000; n_ph=0.00026200;}
 if(eta>=1.6530 && eta<1.7400)
 {c_p=0.01440000;r_p=0.15800000; n_p=1.07000000;c_e=0.00497000;r_e=0.00000000; n_e=0.00027700;c_ph=0.00540000;r_ph=0.00265000; n_ph=0.00034100;}
 if(eta>=1.7400 && eta<1.8300)
 {c_p=0.00620000;r_p=0.19300000; n_p=0.59000000;c_e=0.00371000;r_e=0.00000000; n_e=0.00030640;c_ph=0.00600000;r_ph=0.00279000; n_ph=0.00035800;}
 if(eta>=1.8300 && eta<1.9300)
 {c_p=0.01240000;r_p=0.15400000; n_p=0.66000000;c_e=0.00444000;r_e=0.00000000; n_e=0.00030570;c_ph=0.00950000;r_ph=0.00267000; n_ph=0.00033500;}
 if(eta>=1.9300 && eta<2.0430)
 {c_p=0.00960000;r_p=0.17100000; n_p=0.37000000;c_e=0.00497000;r_e=0.00000000; n_e=0.00033390;c_ph=0.00410000;r_ph=0.00362000; n_ph=0.00022500;}
 if(eta>=2.0430 && eta<2.1720)
 {c_p=0.01890000;r_p=0.01000000; n_p=0.73500000;c_e=0.00466000;r_e=0.00000000; n_e=0.00036940;c_ph=0.01090000;r_ph=0.00287000; n_ph=0.00032600;}
 if(eta>=2.1720 && eta<2.3220)
 {c_p=0.01130000;r_p=0.14700000; n_p=0.00000000;c_e=0.00472000;r_e=0.00000000; n_e=0.00045540;c_ph=0.00400000;r_ph=0.00446000; n_ph=0.00000000;}
 if(eta>=2.3220)
 {c_p=0.01460000;r_p=0.12350000; n_p=0.00000000;c_e=0.00360000;r_e=0.00000000; n_e=0.00059400;c_ph=0.01550000;r_ph=0.00370000; n_ph=0.00035000;}

 double theta=2*atan(exp(-eta));


 if(err_type==0){return sqrt(c_p*c_p*(p*p*sin(theta)*sin(theta))+r_p*r_p*(p*sin(theta))+n_p*n_p);}
 else if(err_type==1){return sqrt(c_e*c_e/(p*p)+r_e*r_e/p+n_e*n_e);}
 else if(err_type==2){return sqrt(c_ph*c_ph/(p*p)+r_ph*r_ph/p+n_ph*n_ph);}
 else { return 0; } // THIS WILL CAUSE THE FIT TO CRASH --> err_type must be 0,1,2

}


//
void RunExclusiveTop(TString filename,
                     TString outname,
                     Int_t channelSelection,
                     Int_t chargeSelection,
                     TH1F *normH,
                     TH1F *genPU,
                     TString era,
                     Bool_t debug)
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

    //  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
    bool isTTbar = 1;

    //PREPARE OUTPUT
    TString baseName=gSystem->BaseName(outname);
    TString dirName=gSystem->DirName(outname);
    TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
    fOut->cd();

    //READ TREE FROM FILE
    MiniEvent_t ev;
    TFile *f = TFile::Open(filename);
    TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
    TTree *t = (TTree*)f->Get("analysis/tree");
    attachToMiniEventTree(t,ev);
    Int_t nentries(t->GetEntriesFast());

    //debug = 1; // manual turn on DEBUG mode
    if (debug) nentries = 5000; //restrict number of entries for testing
    t->GetEntry(0);

    std::cout << "--- producing " << outname << " from " << nentries << " events" << std::endl;

    //auxiliary to solve neutrino pZ using MET
    MEzCalculator neutrinoPzComputer;

    //LUMINOSITY+PILEUP
    LumiTools lumi(era,genPU);

    //LEPTON EFFICIENCIES
    //EfficiencyScaleFactorsWrapper lepEffH(filename.Contains("Data13TeV"),era);

    //B-TAG CALIBRATION
    //BTagSFUtil btvSF(era,"DeepCSV",BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
    BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);

    //JEC/JER
    JECTools jec(era);

    // BOOK OUTPUT TREE
    TTree *outT=new TTree("tree","tree");
    outT->Branch("run",&ev.run,"run/i");
    outT->Branch("event",&ev.event,"event/l");
    outT->Branch("lumi",&ev.lumi,"lumi/i");
    outT->Branch("nvtx",&ev.nvtx,"nvtx/I");

    ADDVAR(&ev.met_pt,"met_pt","/F",outT);
    ADDVAR(&ev.met_phi,"met_phi","/F",outT);
    ADDVAR(&ev.met_sig,"met_sig","/F",outT);
    TString fvars[]={
        "t_pt","t_eta", "t_phi", "t_m", "t_charge", "t_isHadronic",
        "tbar_pt","tbar_eta", "tbar_phi", "tbar_m",
        "ttbar_pt","ttbar_eta", "ttbar_phi", "ttbar_m", "ttbar_E",
        "gen_ttbar_pt","gen_ttbar_eta", "gen_ttbar_phi", "gen_ttbar_m", "gen_ttbar_E",
        "gen_t_pt","gen_t_eta", "gen_t_phi", "gen_t_m",
        "gen_tbar_pt","gen_tbar_eta", "gen_tbar_phi", "gen_tbar_m",
        "l_pt", "l_eta", "l_phi", "l_m", "l_E",
        "nu_pt", "nu_eta", "nu_phi",
        "p1_xi", "p2_xi",
#ifdef SAVEERRORS_ON
        "e_l_px", "e_l_py", "e_l_pz",
        "e_met_px", "e_met_py", "e_met_pxpy",
        "e_bJet0_px", "e_bJet0_py", "e_bJet0_pz",
        "e_bJet1_px", "e_bJet1_py", "e_bJet1_pz",
        "e_lightJet0_px", "e_lightJet0_py", "e_lightJet0_pz",
        "e_lightJet1_px", "e_lightJet1_py", "e_lightJet1_pz",
#endif
#ifdef KINFIT_ON
        "nu_pt_uncorrected", "nu_eta_uncorrected", "nu_phi_uncorrected",
        "nu_pt_min", "nu_eta_min", "nu_phi_min", "nu_pt_max", "nu_eta_max", "nu_phi_max",
#endif
        "chisquareKinFit", "chisquareMattFit", "chisquareAnalyticMattFit",
        "bJet0_pt","bJet0_eta", "bJet0_phi", "bJet0_m", "bJet0_E",
        "bJet1_pt","bJet1_eta", "bJet1_phi", "bJet1_m", "bJet1_E",
        "bJet2_pt","bJet2_eta", "bJet2_phi", "bJet2_m", "bJet2_E",
        "bJet3_pt","bJet3_eta", "bJet3_phi", "bJet3_m", "bJet3_E",

        "lightJet0_pt", "lightJet0_eta", "lightJet0_phi", "lightJet0_m", "lightJet0_E",
        "lightJet1_pt", "lightJet1_eta", "lightJet1_phi", "lightJet1_m", "lightJet1_E",
        "lightJet2_pt", "lightJet2_eta", "lightJet2_phi", "lightJet2_m", "lightJet2_E",
        "lightJet3_pt", "lightJet3_eta", "lightJet3_phi", "lightJet3_m", "lightJet3_E",

        "bJet0_px","bJet0_py", "bJet0_pz",
        "bJet1_px","bJet1_py", "bJet1_pz",
        "lightJet0_px", "lightJet0_py", "lightJet0_pz",
        "lightJet1_px", "lightJet1_py", "lightJet1_pz",

        "gen_b_pt","gen_b_eta", "gen_b_phi", "gen_b_m",
        "gen_bbar_pt","gen_bbar_eta", "gen_bbar_phi", "gen_bbar_m",
        "nBjets", "nLightJets", "nJets", "ht", "lepton_isolation"
        //           "lightJets_pt", "lightJets_eta", "bJets_pt", "bJets_eta"
        // added later as vectors:   lightJets_pt and bJets_pt lightJets_eta and bJets_eta
    };

    std::map<TString,Float_t> outVars;
    for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++){
        outVars[fvars[i]]=0.;
        ADDVAR(&(outVars[fvars[i]]),fvars[i],"/F",outT);
    }
    vector<double> vLightJet_pt, vLightJet_eta;
    vector<double> vBJet_pt, vBJet_eta;

    outT->Branch("bJets_pt", &vBJet_pt);
    outT->Branch("bJets_eta", &vBJet_eta);
    outT->Branch("lightJets_pt", &vLightJet_pt);
    outT->Branch("lightJets_eta", &vLightJet_eta);

#ifdef HISTOGRAMS_ON
    //BOOK HISTOGRAMS
    HistTool ht;
    ht.setNsyst(0);
    ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
    ht.addHist("ch_tag",       new TH1F("ch_tag",      ";Channel Tag;Events",10,0,10));
    ht.addHist("evt_count",    new TH1F("evt_count",   ";Selection Stage;Events",10,0,10));
    ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",55,-0.5,49.5));
    ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",15,-0.5,14.5));
    ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",10,-0.5,9.5));
    ht.addHist("ht",           new TH1F("ht",          ";H_{T} [GeV];Events",500,0,2500));

    ht.addHist("mttbar_rec",   new TH1F("mttbar_rec",  ";M_{ttbar,rec} [GeV];Events",50,300,1200));
    ht.addHist("mttbar_gen",   new TH1F("mttbar_gen",  ";M_{ttbar,gen} [GeV];Events",50,300,1200));
    ht.addHist("mttbar_res",   new TH1F("mttbar_res",  ";M_{ttbar,rec}-M_{ttbar,gen} [GeV];Events",100,-500,500));

    ht.addHist("ptttbar_rec", new TH1F("ptttbar_rec",";Pt_{ttbar,rec} [GeV];Events",50,0,100));
    ht.addHist("ptttbar_gen", new TH1F("ptttbar_gen",";Pt_{ttbar,gen} [GeV];Events",50,0,100));
    ht.addHist("ptttbar_res", new TH1F("ptttbar_res",";Pt_{ttbar,rec}-Pt_{ttbar,gen} [GeV];Events",50,-150,150));

    ht.addHist("yttbar_rec", new TH1F("yttbar_rec",";Y_{ttbar,rec} ;Events",75,-2.5,2.5));
    ht.addHist("yttbar_gen", new TH1F("yttbar_gen",";Y_{ttbar,gen} ;Events",75,-2.5,2.5));
    ht.addHist("yttbar_res", new TH1F("yttbar_res",";Y_{ttbar,rec}-Y_{ttbar,gen} ;Events",75,-1.5,1.5));

    ht.addHist("mt_res"   ,   new TH1F("mt_res",  ";M_{t,rec}-M_{t,gen} [GeV];Events",50,-700,700));
    ht.addHist("mtbar_res",   new TH1F("mtbar_res",  ";M_{tbar,rec}-M_{tbar,gen} [GeV];Events",50,-700,700));
    ht.addHist("ptt_res",     new TH1F("ptt_res",";Pt_{t,rec}-Pt_{t,gen} [GeV];Events",50,-150,150));
    ht.addHist("pttbar_res",  new TH1F("pttbar_res",";Pt_{tbar,rec}-Pt_{ttbar,gen} [GeV];Events",50,-150,150));
    ht.addHist("yt_res",    new TH1F("yt_res",";Y_{t,rec}-Y_{t,gen} ;Events",75,-2.5,2.5));
    ht.addHist("ytbar_res", new TH1F("ytbar_res",";Y_{tbar,rec}-Y_{tbar,gen} ;Events",75,-2.5,2.5));

    ht.addHist("mtop_res_hadronic",   new TH1F("mtop_res_hadronic",  ";M_{top,rec}-M_{top,gen} [GeV];Events",100,-500,500));
    ht.addHist("pttop_res_hadronic",  new TH1F("pttop_res_hadronic",";Pt_{top,rec}-Pt_{top,gen} [GeV];Events",50,-150,150));
    ht.addHist("ytop_res_hadronic",   new TH1F("ytop_res_hadronic",";Y_{top,rec}-Y_{top,gen} ;Events",75,-2.5,2.5));
    ht.addHist("mtop_res_leptonic",   new TH1F("mtop_res_leptonic",  ";M_{top,rec}-M_{top,gen} [GeV];Events",100,-500,500));
    ht.addHist("pttop_res_leptonic",  new TH1F("pttop_res_leptonic",";Pt_{top,rec}-Pt_{top,gen} [GeV];Events",50,-150,150));
    ht.addHist("ytop_res_leptonic",   new TH1F("ytop_res_leptonic",";Y_{top,rec}-Y_{top,gen} ;Events",75,-2.5,2.5));
#endif

    std::cout << "--- init done" << std::endl;


    //EVENT SELECTION WRAPPER
    SelectionTool selector(filename, false, triggerList);

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

        //////////////////
        // CORRECTIONS //
        /////////////////
        btvSF.addBTagDecisions(ev);
        btvSF.updateBTagDecisions(ev);
        jec.smearJetEnergies(ev);

        //////////////////////////
        // RECO LEVEL SELECTION //
        //////////////////////////
        TString chTag = selector.flagFinalState(ev); // writes the name in chTag
        // ch
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 0, plotwgts); // count all events before any selection
        Int_t ch_tag = 0;
        if      (chTag=="EM")    ch_tag = 1;
        else if (chTag=="MM")    ch_tag = 2;
        else if (chTag=="EE")    ch_tag = 3;
        else if (chTag=="E")     ch_tag = 4;
        else if (chTag=="M")     ch_tag = 5;
        else                     ch_tag = 9;
        ht.fill("ch_tag", ch_tag, plotwgts);
#endif
        if(chTag!="E" && chTag!="M" )   continue; // events with electrons (id=11) or muons (id=13)
//        if(chTag!="M")                continue; // ONLY events with muons (id=13), not electrons (id=11)
//        if(chTag!="E")                continue; // ONLY events with electrons (id=11), not muons (id=13)
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 1, plotwgts); // count events after channel selection
#endif
        std::vector<Particle> &leptons     = selector.getSelLeptons();
        std::vector<Jet>      &jets        = selector.getJets();
        std::vector<Jet>      bJets,lightJets;
        std::vector<Particle> selectedLeptons;

        //  selection of lightJets and bJets
        for(size_t ij=0; ij<jets.size(); ij++) {
            if(jets[ij].flavor()==5) bJets.push_back(jets[ij]);
            else                     lightJets.push_back(jets[ij]);
        }

        // selection of leptons
        for( size_t i_lept=0;i_lept<leptons.size();i_lept++) {
            if (leptons[i_lept].pt()<25.) continue;
            //        if (leptons[i_lept].reliso()>0.10) continue;   //usually tighter
            selectedLeptons.push_back(leptons[i_lept]);
        }
        // ---- EVENT SELECTION --------------------------------------------------------------
        if (selectedLeptons.size()!=1) continue; // ONLY events with 1 selected lepton
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 2, plotwgts); // count events after selection on number of leptons (SHOULD BE SAME)
#endif
        if ( jets.size()  < 4 )        continue; // ONLY events with at least 4 jets
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 3, plotwgts); // count events after selection on number of jets
#endif
        if ( bJets.size() < 2 )        continue; // ONLY events with at least 2 BJets
#ifdef HISTOGRAMS_ON
       ht.fill("evt_count", 4, plotwgts); // count events after selection on number of Bjets
       ht.fill("puwgtctr",0,plotwgts);
#endif

        if (!ev.isData) {
            wgt  = (normH? normH->GetBinContent(1) : 1.0);          // norm weight
            double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);     // pu weight
            std::vector<double>puPlotWgts(1,puWgt);

#ifdef HISTOGRAMS_ON
            ht.fill("puwgtctr",1,puPlotWgts);
#endif

            // lepton trigger*selection weights
            EffCorrection_t trigSF(1.0,0.0);//lepEffH.getTriggerCorrection(leptons,{},{},period);
            EffCorrection_t  selSF(1.0,0.0); //lepEffH.getOfflineCorrection(leptons[0], period);
            wgt *= puWgt*trigSF.first*selSF.first;

            //top pt weighting
            double topptsf = 1.0;
            if(isTTbar) {
                for (int igen=0; igen<ev.ngtop; igen++) {
                    if(abs(ev.gtop_id[igen])!=6) continue;
                    topptsf *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
                }
            }
            wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);                   // generator level weights
            plotwgts[0]=wgt;                                        //update weight for plotter
        }
		


// ----- START RECONSTRUCTION OF TTBAR -------------------------------------------------------
        bool isThadronic = 0;
        if(bJets.size()>=2 && lightJets.size()>=2)        {    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            //determine the neutrino kinematics
            TLorentzVector met(0,0,0,0);
            met.SetPtEtaPhiM(ev.met_pt,0.,ev.met_phi,0.);
            neutrinoPzComputer.SetMET(met);
            neutrinoPzComputer.SetLepton(leptons[0].p4());
            float nupz=neutrinoPzComputer.Calculate();
            TLorentzVector neutrinoP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));

            // compute scalar ht
            float scalarht(0.);
            for(size_t ij=0; ij<jets.size(); ij++) {
                    scalarht += jets[ij].pt();
            }
            scalarht += leptons[0].pt();
            scalarht += neutrinoP4.Pt();

            // create candidates for the t and tbar by coupling the jets and leptons
            // T1 is hadronic, T2 is leptonic
//            TLorentzVector firstT1Candidate(bJets[0].p4()+lightJets[0].p4()+lightJets[1].p4());
//            TLorentzVector firstT2Candidate(bJets[1].p4()+leptons[0].p4()+neutrinoP4 );
//            TLorentzVector secondT1Candidate(bJets[1].p4()+lightJets[0].p4()+lightJets[1].p4());
//            TLorentzVector secondT2Candidate(bJets[0].p4()+leptons[0].p4()+neutrinoP4 );
//            double mass_dev_1st = TMath::Power( firstT1Candidate.M()-m_TOP ,2) + TMath::Power( firstT2Candidate.M()-m_TOP ,2);
//            double mass_dev_2nd = TMath::Power( secondT1Candidate.M()-m_TOP ,2) + TMath::Power( secondT2Candidate.M()-m_TOP ,2);
            TLorentzVector t_rec;
            TLorentzVector tbar_rec;
            TLorentzVector bJet_Leptonic; // this is the bJet of the leptonic top which will be used in the kinematic fit
            TLorentzVector bJet_Hadronic;
            // the combination whose mass is closer to the top one is chosen
/*            if ( mass_dev_1st < mass_dev_2nd ) {
                if (leptons[0].charge()>0.) { //the T2 is the one decaying in leptons: check the lepton charge to know if it's t or tbar
                    tbar_rec = firstT1Candidate;
                    t_rec = firstT2Candidate;
                    bJet_Leptonic = bJets[1].p4();
                    bJet_Hadronic = bJets[0].p4();
                    isThadronic = 0; //in this case the top has decayed into leptons
                } else {
                    t_rec = firstT1Candidate;
                    tbar_rec = firstT2Candidate;
                    bJet_Leptonic = bJets[1].p4();
                    bJet_Hadronic = bJets[0].p4();
                    isThadronic = 1;
                }
            } else {
                if (leptons[0].charge()>0.) {
                    tbar_rec = secondT1Candidate;
                    t_rec = secondT2Candidate;
                    bJet_Leptonic = bJets[0].p4();
                    bJet_Hadronic = bJets[1].p4();
                    isThadronic = 0;
                } else {
                    t_rec = secondT1Candidate;
                    tbar_rec = secondT2Candidate;
                    bJet_Leptonic = bJets[0].p4();
                    bJet_Hadronic = bJets[1].p4();
                    isThadronic = 1;
                }
            }
*/

            vector <double> mistake;        // mass difference between top candidate and true mass
            vector <int> sequence;          // indexes of:   bJetHad, lightJet1, lightJet2, bJetLep
            vector <vector <int>> marker;   // list of all sequences
            TLorentzVector t_rec_had;
            TLorentzVector t_rec_lep;

            // the combination whose mass is closer to the top one is chosen

            // candidate bJets[0]+lightJets[0]+lightJets[1]+bJets[1]
            t_rec_had = bJets[0].p4()+lightJets[0].p4()+lightJets[1].p4();   // b1+q1+q2
            t_rec_lep = bJets[1].p4()+leptons[0].p4()+neutrinoP4;            // b2+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(0);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(1);
            marker.push_back(sequence);
            sequence.clear();

            // candidate bJets[0]+lightJets[0]+lightJets[1]+bJets[2]
            t_rec_had = bJets[0].p4()+lightJets[0].p4()+lightJets[1].p4();   // b1+q1+q2
            t_rec_lep = bJets[2].p4()+leptons[0].p4()+neutrinoP4;            // b3+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(0);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(2);
            marker.push_back(sequence);
            sequence.clear();

            // candidate bJets[0]+lightJets[0]+lightJets[1]+bJets[3]
            t_rec_had = bJets[0].p4()+lightJets[0].p4()+lightJets[1].p4();   // b1+q1+q2
            t_rec_lep = bJets[3].p4()+leptons[0].p4()+neutrinoP4;            // b4+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(0);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(3);
            marker.push_back(sequence);
            sequence.clear();

            // candidate bJets[1]+lightJets[0]+lightJets[1]+bJets[2]
            t_rec_had = bJets[1].p4()+lightJets[0].p4()+lightJets[1].p4();   // b2+q1+q2
            t_rec_lep = bJets[2].p4()+leptons[0].p4()+neutrinoP4;            // b3+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(1);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(2);
            marker.push_back(sequence);
            sequence.clear();

            // candidate bJets[0]+lightJets[0]+lightJets[1]+bJets[1]
            t_rec_had = bJets[2].p4()+lightJets[0].p4()+lightJets[1].p4();   // b3+q1+q2
            t_rec_lep = bJets[3].p4()+leptons[0].p4()+neutrinoP4;            // b4+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(2);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(3);
            marker.push_back(sequence);
            sequence.clear();

// ------- symmetric
            // candidate bJets[1]+lightJets[0]+lightJets[1]+bJets[0]
            t_rec_had = bJets[1].p4()+lightJets[0].p4()+lightJets[1].p4();   // b2+q1+q2
            t_rec_lep = bJets[0].p4()+leptons[0].p4()+neutrinoP4;            // b1+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(1);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(0);
            marker.push_back(sequence);
            sequence.clear();

            // candidate bJets[2]+lightJets[0]+lightJets[1]+bJets[0]
            t_rec_had = bJets[2].p4()+lightJets[0].p4()+lightJets[1].p4();   // b3+q1+q2
            t_rec_lep = bJets[0].p4()+leptons[0].p4()+neutrinoP4;            // b1+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(2);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(0);
            marker.push_back(sequence);
            sequence.clear();

            // candidate bJets[3]+lightJets[0]+lightJets[1]+bJets[0]
            t_rec_had = bJets[3].p4()+lightJets[0].p4()+lightJets[1].p4();   // b4+q1+q2
            t_rec_lep = bJets[0].p4()+leptons[0].p4()+neutrinoP4;            // b1+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(3);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(0);
            marker.push_back(sequence);
            sequence.clear();

            // candidate bJets[2]+lightJets[0]+lightJets[1]+bJets[1]
            t_rec_had = bJets[2].p4()+lightJets[0].p4()+lightJets[1].p4();   // b3+q1+q2
            t_rec_lep = bJets[1].p4()+leptons[0].p4()+neutrinoP4;            // b2+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(2);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(1);
            marker.push_back(sequence);
            sequence.clear();

            // candidate bJets[3]+lightJets[0]+lightJets[1]+bJets[2]
            t_rec_had = bJets[3].p4()+lightJets[0].p4()+lightJets[1].p4();   // b4+q1+q2
            t_rec_lep = bJets[2].p4()+leptons[0].p4()+neutrinoP4;            // b3+l+nu
            mistake.push_back(abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP));
            sequence.push_back(3);
            sequence.push_back(0);
            sequence.push_back(1);
            sequence.push_back(2);
            marker.push_back(sequence);
            sequence.clear();

            int correct_index=dump_index(mistake,mistake.size());   // index of best candidate
            vector <int> correct_sequence=marker[correct_index];    // sequence of best candidate

            if (leptons[0].charge()>0.) {
                t_rec    = bJets[ correct_sequence[3] ].p4()+leptons[0].p4()+neutrinoP4;
                tbar_rec = bJets[ correct_sequence[0] ].p4()+lightJets[ correct_sequence[1] ].p4()+lightJets[ correct_sequence[2] ].p4();
                bJet_Leptonic = bJets[ correct_sequence[3] ].p4();
                bJet_Hadronic = bJets[ correct_sequence[0] ].p4();
                isThadronic = 0;
            } else {
                t_rec    = bJets[ correct_sequence[0] ].p4()+leptons[0].p4()+neutrinoP4;
                tbar_rec = bJets[ correct_sequence[3] ].p4()+lightJets[ correct_sequence[1] ].p4()+lightJets[ correct_sequence[2] ].p4();
                bJet_Leptonic = bJets[ correct_sequence[0] ].p4();
                bJet_Hadronic = bJets[ correct_sequence[3] ].p4();
                isThadronic = 1;
            }

            TF1 err_theta_fracdeltaeta("err_theta","abs(-2*(exp(-x)/(1+pow(exp(-x),2))))",-1e10,1e10);

            // compute uncertainties
            //---   b1 = bJet_Hadronic (e_bJet0_p)   , b2=bJet_Leptonic
            double pt_b1     =   bJet_Hadronic.Pt();
            double eta_b1    =   bJet_Hadronic.Rapidity();
            double phi_b1    =   bJet_Hadronic.Phi();

            double pt_b2     =   bJet_Leptonic.Pt();
            double eta_b2    =   bJet_Leptonic.Rapidity();
            double phi_b2    =   bJet_Leptonic.Phi();

            double pt_q1     =   lightJets[ correct_sequence[1] ].p4().Pt();
            double eta_q1    =   lightJets[ correct_sequence[1] ].p4().Rapidity();
            double phi_q1    =   lightJets[ correct_sequence[1] ].p4().Phi();

            double pt_q2     =   lightJets[ correct_sequence[2] ].p4().Pt();
            double eta_q2    =   lightJets[ correct_sequence[2] ].p4().Rapidity();
            double phi_q2    =   lightJets[ correct_sequence[2] ].p4().Phi();

            double pt_l     =   leptons[0].p4().Pt();
            double eta_l    =   leptons[0].p4().Rapidity();
            double phi_l    =   leptons[0].p4().Phi();
            double m_l      =   leptons[0].p4().M();

            double pt_nu     =   neutrinoP4.Pt();
            double eta_nu    =   neutrinoP4.Rapidity();
            double phi_nu    =   neutrinoP4.Phi();
            double m_nu      =   neutrinoP4.M();

            // xyz components
            double bar_p_x_b1   = pt_b1*cos(phi_b1); // bJet_Hadronic.Px();  //b1 è adronico
            double bar_p_y_b1   = pt_b1*sin(phi_b1);
            double bar_p_z_b1   = pt_b1/tan(theta.Eval(eta_b1));
            double bar_p_x_b2   = pt_b2*cos(phi_b2); // bJet_Leptonic.Px();  //b2 è leptonico
            double bar_p_y_b2   = pt_b2*sin(phi_b2);
            double bar_p_z_b2   = pt_b2/tan(theta.Eval(eta_b2));
            double bar_p_x_l    = pt_l*cos(phi_l);   //   leptons[0].p4().Px();
            double bar_p_y_l    = pt_l*sin(phi_l);
            double bar_p_z_l    = pt_l/tan(theta.Eval(eta_l));
            double bar_p_x_nu   = pt_nu*cos(phi_nu);    // neutrinoP4.Px();
            double bar_p_y_nu   = pt_nu*sin(phi_nu);
            double bar_p_z_nu   = pt_nu/tan(theta.Eval(eta_nu));
            double bar_p_x_q1   = pt_q1*cos(phi_q1);    // lightJets[ correct_sequence[1] ].p4().Px();
            double bar_p_y_q1   = pt_q1*sin(phi_q1);;
            double bar_p_z_q1   = pt_q1/tan(theta.Eval(eta_q1));
            double bar_p_x_q2   = pt_q2*cos(phi_q2);    // lightJets[ correct_sequence[2] ].p4().Px();
            double bar_p_y_q2   = pt_q2*cos(phi_q2);
            double bar_p_z_q2   = pt_q2/tan(theta.Eval(eta_q2));

            // error xyz components
            double e_b1_eta=bjet_err(abs(pt_b1/sin(theta.Eval(eta_b1))),eta_b1,phi_b1,1);
            double e_b1_p=sqrt(pow(bjet_err(abs(pt_b1/sin(theta.Eval(eta_b1))),eta_b1,phi_b1,0)/(sin(theta.Eval(eta_b1))),2)+pow(pt_b1/pow(sin(theta.Eval(eta_b1)),2)*cos(theta.Eval(eta_b1))*err_theta_fracdeltaeta.Eval(eta_b1)*e_b1_eta,2));
            double e_b1_phi=bjet_err(abs(pt_b1/sin(theta.Eval(eta_b1))),eta_b1,phi_b1,2);
            double e_b1_theta=e_b1_eta*err_theta_fracdeltaeta.Eval(eta_b1);

            double e_b2_eta=bjet_err(abs(pt_b2/sin(theta.Eval(eta_b2))),eta_b2,phi_b2,1);
            double e_b2_p=sqrt(pow(bjet_err(abs(pt_b2/sin(theta.Eval(eta_b2))),eta_b2,phi_b2,0)/(sin(theta.Eval(eta_b2))),2)+pow(pt_b2/pow(sin(theta.Eval(eta_b2)),2)*cos(theta.Eval(eta_b2))*err_theta_fracdeltaeta.Eval(eta_b2)*e_b2_eta,2));
            double e_b2_phi=bjet_err(abs(pt_b2/sin(theta.Eval(eta_b2))),eta_b1,phi_b2,2);
            double e_b2_theta=e_b2_eta*err_theta_fracdeltaeta.Eval(eta_b2);

            double e_q1_eta=ljet_err(abs(pt_q1/sin(theta.Eval(eta_q1))),eta_q1,phi_q1,1);
            double e_q1_p=sqrt(pow(ljet_err(abs(pt_q1/sin(theta.Eval(eta_q1))),eta_q1,phi_q1,0)/(sin(theta.Eval(eta_q1))),2)+pow(pt_q1/pow(sin(theta.Eval(eta_q1)),2)*cos(theta.Eval(eta_q1))*err_theta_fracdeltaeta.Eval(eta_q1)*e_q1_eta,2));
            double e_q1_phi=ljet_err(abs(pt_q1/sin(theta.Eval(eta_q1))),eta_q1,phi_q1,2);
            double e_q1_theta=e_q1_eta*err_theta_fracdeltaeta.Eval(eta_q1);

            double e_q2_eta=ljet_err(abs(pt_q2/sin(theta.Eval(eta_q2))),eta_q2,phi_q2,1);
            double e_q2_p=sqrt(pow(ljet_err(abs(pt_q2/sin(theta.Eval(eta_q2))),eta_q2,phi_q2,0)/(sin(theta.Eval(eta_q2))),2)+pow(pt_q2/pow(sin(theta.Eval(eta_q2)),2)*cos(theta.Eval(eta_q2))*err_theta_fracdeltaeta.Eval(eta_q2)*e_q2_eta,2));
            double e_q2_phi=ljet_err(abs(pt_q2/sin(theta.Eval(eta_q2))),eta_q2,phi_q2,2);
            double e_q2_theta=e_q2_eta*err_theta_fracdeltaeta.Eval(eta_q2);

            double e_l_eta,e_l_p,e_l_phi,e_l_theta;
            if(m_l>=0.1){
              e_l_eta=mu_err(abs(pt_l/sin(theta.Eval(eta_l))),eta_l,phi_l,1);
              e_l_p=sqrt(pow(mu_err(abs(pt_l/sin(theta.Eval(eta_l))),eta_l,phi_l,0)/(sin(theta.Eval(eta_l))),2)+pow(pt_l/pow(sin(theta.Eval(eta_l)),2)*cos(theta.Eval(eta_l))*err_theta_fracdeltaeta.Eval(eta_l)*e_l_eta,2));
              e_l_phi=mu_err(abs(pt_l/sin(theta.Eval(eta_l))),eta_l,phi_l,2);
              e_l_theta=e_l_eta*err_theta_fracdeltaeta.Eval(eta_l);
            } else {
              e_l_eta=e_err(abs(pt_l/sin(theta.Eval(eta_l))),eta_l,phi_l,1);
              e_l_p=sqrt(pow(e_err(abs(pt_l/sin(theta.Eval(eta_l))),eta_l,phi_l,0)/(sin(theta.Eval(eta_l))),2)+pow(pt_l/pow(sin(theta.Eval(eta_l)),2)*cos(theta.Eval(eta_l))*err_theta_fracdeltaeta.Eval(eta_l)*e_l_eta,2));
              e_l_phi=e_err(abs(pt_l/sin(theta.Eval(eta_l))),eta_l,phi_l,2);
              e_l_theta=e_l_eta*err_theta_fracdeltaeta.Eval(eta_l);
            }

            // compute x,y,z components of uncertainties
            double e_b1_x=sqrt(pow(e_b1_p*sin(theta.Eval(eta_b1))*cos(phi_b1),2)+pow(pt_b1/sin(theta.Eval(eta_b1))*cos(theta.Eval(eta_b1))*cos(phi_b1)*e_b1_theta,2)+pow(pt_b1/sin(theta.Eval(eta_b1))*sin(theta.Eval(eta_b1))*sin(phi_b1)*e_b1_phi,2));
              double e_b1_y=sqrt(pow(e_b1_p*sin(theta.Eval(eta_b1))*sin(phi_b1),2)+pow(pt_b1/sin(theta.Eval(eta_b1))*cos(theta.Eval(eta_b1))*sin(phi_b1)*e_b1_theta,2)+pow(pt_b1/sin(theta.Eval(eta_b1))*sin(theta.Eval(eta_b1))*cos(phi_b1)*e_b1_phi,2));
              double e_b1_z=sqrt(pow(e_b1_p*cos(theta.Eval(eta_b1)),2)+pow(pt_b1/sin(theta.Eval(eta_b1))*sin(theta.Eval(eta_b1))*e_b1_theta,2));

            double e_b2_x=sqrt(pow(e_b2_p*sin(theta.Eval(eta_b2))*cos(phi_b2),2)+pow(pt_b2/sin(theta.Eval(eta_b2))*cos(theta.Eval(eta_b2))*cos(phi_b2)*e_b2_theta,2)+pow(pt_b2/sin(theta.Eval(eta_b2))*sin(theta.Eval(eta_b2))*sin(phi_b2)*e_b2_phi,2));
            double e_b2_y=sqrt(pow(e_b2_p*sin(theta.Eval(eta_b2))*sin(phi_b2),2)+pow(pt_b2/sin(theta.Eval(eta_b2))*cos(theta.Eval(eta_b2))*sin(phi_b2)*e_b2_theta,2)+pow(pt_b2/sin(theta.Eval(eta_b2))*sin(theta.Eval(eta_b2))*cos(phi_b2)*e_b2_phi,2));
            double e_b2_z=sqrt(pow(e_b2_p*cos(theta.Eval(eta_b2)),2)+pow(pt_b2/sin(theta.Eval(phi_b2))*sin(theta.Eval(eta_b2))*e_b2_theta,2));

            double e_l_x=sqrt(pow(e_l_p*sin(theta.Eval(eta_l))*cos(phi_l),2)+pow(pt_l/sin(theta.Eval(eta_l))*cos(theta.Eval(eta_l))*cos(phi_l)*e_l_theta,2)+pow(pt_l/sin(theta.Eval(eta_l))*sin(theta.Eval(eta_l))*sin(phi_l)*e_l_phi,2));
            double e_l_y=sqrt(pow(e_l_p*sin(theta.Eval(eta_l))*sin(phi_l),2)+pow(pt_b1/sin(theta.Eval(eta_l))*cos(theta.Eval(eta_l))*sin(phi_l)*e_l_theta,2)+pow(pt_l/sin(theta.Eval(eta_l))*sin(theta.Eval(eta_l))*cos(phi_l)*e_l_phi,2));;
            double e_l_z=sqrt(pow(e_l_p*cos(theta.Eval(eta_l)),2)+pow(pt_l/sin(theta.Eval(eta_l))*sin(theta.Eval(eta_l))*e_l_theta,2));

            double e_nu_x=6.*neutrinoP4.Pt(); // 6 times because six is big enough. We set a very big error on neutrino
            double e_nu_y=6.*neutrinoP4.Pt();
            double e_nu_z=6.*neutrinoP4.Pt();

            double e_q1_x=sqrt(pow(e_q1_p*sin(theta.Eval(eta_q1))*cos(phi_q1),2)+pow(pt_q1/sin(theta.Eval(eta_q1))*cos(theta.Eval(eta_q1))*cos(phi_q1)*e_q1_theta,2)+pow(pt_q1/sin(theta.Eval(eta_q1))*sin(theta.Eval(eta_q1))*sin(phi_b1)*e_q1_phi,2));
            double e_q1_y=sqrt(pow(e_q1_p*sin(theta.Eval(eta_q1))*sin(phi_q1),2)+pow(pt_q1/sin(theta.Eval(eta_q1))*cos(theta.Eval(eta_q1))*sin(phi_q1)*e_q1_theta,2)+pow(pt_q1/sin(theta.Eval(eta_q1))*sin(theta.Eval(eta_q1))*cos(phi_q1)*e_q1_phi,2));
            double e_q1_z=sqrt(pow(e_q1_p*cos(theta.Eval(eta_q1)),2)+pow(pt_q1/sin(theta.Eval(eta_q1))*sin(theta.Eval(eta_q1))*e_q1_theta,2));

            double e_q2_x=sqrt(pow(e_q2_p*sin(theta.Eval(eta_q2))*cos(phi_q2),2)+pow(pt_q1/sin(theta.Eval(eta_q2))*cos(theta.Eval(eta_q2))*cos(phi_q2)*e_q2_theta,2)+pow(pt_q2/sin(theta.Eval(eta_q2))*sin(theta.Eval(eta_q2))*sin(phi_b2)*e_q2_phi,2));
            double e_q2_y=sqrt(pow(e_q2_p*sin(theta.Eval(eta_q2))*sin(phi_q2),2)+pow(pt_q2/sin(theta.Eval(eta_q2))*cos(theta.Eval(eta_q2))*sin(phi_q2)*e_q2_theta,2)+pow(pt_q2/sin(theta.Eval(eta_q2))*sin(theta.Eval(eta_q2))*cos(phi_q2)*e_q2_phi,2));
            double e_q2_z=sqrt(pow(e_q2_p*cos(theta.Eval(eta_q2)),2)+pow(pt_q2/sin(theta.Eval(eta_q2))*sin(theta.Eval(eta_q2))*e_q2_theta,2));

            // ----- START KINEMATIC FIT -----------------------------------------------------------------
            TLorentzVector lepton = leptons[0].p4();  // <-- MOVE EARLIER IN CODE !!!!!!!!!!!!!!!
            double chisquareKinFit  = -1.;
            double chisquareMattFit = -1.;
            double chisquareAnalyticMattFit = -1.;

#ifdef KINFIT_ON
            TLorentzVector t_rec_uncorr      = t_rec;
            TLorentzVector tbar_rec_uncorr   = tbar_rec;
            TLorentzVector neutrinoP4_uncorr = neutrinoP4;

            TLorentzVector *neutrinoP4_new = new TLorentzVector(); // neutrino after KinFit correction
            TLorentzVector *neutrinoP4_min = new TLorentzVector(); // min value of neutrino after KinFit correction
            TLorentzVector *neutrinoP4_max = new TLorentzVector(); // max value of neutrino after KinFit correction

            // apply kinematic fit on leptonic top and re-define it with the corrected MET, all other quantities are left unchanged
            NeutrinoEllipseCalculator n( bJet_Leptonic , lepton , m_TOP , m_W, m_NU ); // mass of neutrino is always set to 0.
            TMatrixD ma=n.getNeutrinoEllipse();
            TMatrixD V_0(3,3); // kinematic values of MET
            V_0[0][0]=0.;            V_0[0][1]=0.;            V_0[0][2]= neutrinoP4.Px(); //met x  ( p_x )
            V_0[1][0]=0.;            V_0[1][1]=0.;            V_0[1][2]= neutrinoP4.Py(); //met y  ( p_y )
            V_0[2][0]=0.;            V_0[2][1]=0.;            V_0[2][2]= neutrinoP4.Pz(); //met z  ( p_z )
            TMatrixD sigma(3,3);
            sigma[0][0]=0.3*abs(neutrinoP4.Px())+1.e-3;   // err met x
            sigma[0][1]=0.;                               // correlations
            sigma[0][2]=0.;
            sigma[1][0]=0.;
            sigma[1][1]=0.3*abs(neutrinoP4.Py())+1.e-3;   // err met y
            sigma[1][2]=0.;
            sigma[2][0]=0.;
            sigma[2][1]=0.;
            sigma[2][2]=1.e9;                             // err met z

            TMatrixD sigmaINVERTED(3,3);
            sigmaINVERTED = sigma.Invert();

            n.KinFit(neutrinoP4_new, neutrinoP4_max, neutrinoP4_min, sigmaINVERTED, V_0 , n.GetH(), chisquareKinFit);
            neutrinoP4 = *neutrinoP4_new;
#endif
#ifdef MattFIT_ON
            const int nPar = 24;
            TMinuit minuit(nPar);
            minuit.SetPrintLevel(-1); // quiet (also suppresse all warnings) https://root.cern.ch/root/html534/src/TMinuit.cxx.html#aN1PoB
            minuit.SetFCN(chistar);
            minuit.SetErrorDef(1.);

            //Temporary redefinition of errors making them huge
            e_b1=1.e9;
            e_b2=1.e9;
            e_l=1.e9;
            e_nu=1.e9;
            e_q1=1.e9;
            e_q2=1.e9;
            e_t=1.e9;
            e_tbar=1.e9;
            e_W_ad=1.e9;
            e_W_lep=1.e9;
            e_phi_b1=1.e9;
            e_phi_b2=1.e9;
            e_phi_l=1.e9;
            e_phi_nu=1.e9;
            e_phi_q1=1.e9;
            e_phi_q2=1.e9;
            e_eta_b1=1.e9;
            e_eta_b2=1.e9;
            e_eta_l=1.e9;
            e_eta_nu=1.e9;
            e_eta_q1=1.e9;
            e_eta_q2=1.e9;

            // !!!!!!!!!!!!!  THIS CAN BE SIMPLIFIED, these variables are not necessary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // ASSIGN VALUES TO ALL QUANTITIES
            // par[0]  = bJet_Hadronic.Pt();
            pt_b1  = bJet_Hadronic.Pt();
            pt_b2  = bJet_Leptonic.Pt();
            pt_l   = lepton.Pt();
            pt_nu  = neutrinoP4.Pt();
            pt_q1  = lightJets[0].Pt();
            pt_q2  = lightJets[1].Pt();
            phi_b1 = bJet_Hadronic.Phi();
            phi_b2 = bJet_Leptonic.Phi();
            phi_l  = lepton.Phi();
            phi_nu = neutrinoP4.Phi();
            phi_q1 = lightJets[0].Phi();
            phi_q2 = lightJets[1].Phi();
            eta_b1 = bJet_Hadronic.Rapidity();
            eta_b2 = bJet_Leptonic.Rapidity();
            eta_l  = lepton.Rapidity();
            eta_nu = neutrinoP4.Rapidity();
            eta_q1 = lightJets[0].Rapidity();
            eta_q2 = lightJets[1].Rapidity();

            // map of    par[]     vector
            // {"jb1_pt","jb2_pt","l_pt","nu_pt","j_q1_pt","j_q2_pt","t_m","t_bar_m","W_ad_m","W_lep_m","phi_b1","phi_b2","phi_l","phi_nu","phi_j_q1","phi_j_q2","eta_b1","eta_b2","eta_l","eta_nu","eta_j_q1","eta_j_q2","lam_1","lam_2"};
            //    0        1         2      3       4         5         6     7         8         9        10       11       12       13       14         15         16       17      18      19        20         21       22      23
            double par[nPar] = { pt_b1,pt_b2,pt_l,pt_nu, pt_q1, pt_q2,m_t,m_tbar,m_W_ad,m_W_lep,phi_b1,phi_b2,phi_l,phi_nu,phi_q1,phi_q2,eta_b1,eta_b2,eta_l,eta_nu,eta_q1,eta_q2,lambda_Value,lambda_Value};
            double step[nPar]={0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1}; // computational parameter whatever
            double min[nPar]={pt_b1-e_b1,pt_b2-e_b2,pt_l-e_l,pt_nu-e_nu,pt_q1-e_q1,pt_q2-e_q2,m_t-e_t,m_tbar-e_tbar,m_W_ad-e_W_ad,m_W_lep-e_W_lep,phi_b1-e_phi_b1,phi_b2-e_phi_b2,phi_l-e_phi_l,phi_nu-e_phi_nu,phi_q1-e_phi_q1,phi_q2-e_phi_q2,eta_b1-e_eta_b1,eta_b2-e_eta_b2,eta_l-e_eta_l,eta_nu-e_eta_nu,eta_q1-e_eta_q1,eta_q2-e_eta_q2,0,0}; // 0 in both min and max errors (with no point) means no limit
            double max[nPar]={pt_b1+e_b1,pt_b2+e_b2,pt_l+e_l,pt_nu+e_nu,pt_q1+e_q1,pt_q2+e_q2,m_t+e_t,m_tbar+e_tbar,m_W_ad+e_W_ad,m_W_lep+e_W_lep,phi_b1+e_phi_b1,phi_b2+e_phi_b2,phi_l+e_phi_l,phi_nu+e_phi_nu,phi_q1+e_phi_q1,phi_q2+e_phi_q2,eta_b1+e_eta_b1,eta_b2+e_eta_b2,eta_l+e_eta_l,eta_nu+e_eta_nu,eta_q1+e_eta_q1,eta_q2+e_eta_q2,0,0};// 0 in both min and max errors (with no point) means no limit
            string cpar[nPar]={"jb1","jb2","l","nu","j_q1","j_q2","t","t_bar","W_ad","W_lep","phi_b1","phi_b2","phi_l","phi_nu","phi_j_q1","phi_j_q2","eta_b1","eta_b2","eta_l","eta_nu","eta_j_q1","eta_j_q2","lam_1","lam_2"};

            // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ASSIGN VALUES TO ALL ERRORS, now are only indicative !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            e_b1=2*0.05*pt_b1;
            e_b2=2*0.05*pt_b2;
            e_l=0.01*pt_l;
            e_nu=0.3*pt_nu;
            e_q1=2*0.05*pt_q1;
            e_q2=2*0.05*pt_q2;
            e_t=0.48;                   //0.76;
            e_tbar=0.48;                 //0.76;
            e_W_ad=2.082;
            e_W_lep=2.082;
            e_phi_b1=0.087;
            e_phi_b2=0.087;
            e_phi_l=0.087;
            e_phi_nu=0.087;
            e_phi_q1=0.087;
            e_phi_q2=0.087;
            e_eta_b1=0.001;
            e_eta_b2=0.001;
            e_eta_l=0.001;
            e_eta_nu=0.001;
            e_eta_q1=0.001;
            e_eta_q2=0.001;

            for(int i=0; i<nPar; i++){
                minuit.DefineParameter(i,cpar[i].c_str(),par[i],step[i],min[i],max[i]);
            }
            minuit.Migrad();

            // map of    par[]     vector
            // {"jb1_pt","jb2_pt","l_pt","nu_pt","j_q1_pt","j_q2_pt","t_m","t_bar_m","W_ad_m","W_lep_m","phi_b1","phi_b2","phi_l","phi_nu","phi_j_q1","phi_j_q2","eta_b1","eta_b2","eta_l","eta_nu","eta_j_q1","eta_j_q2","lam_1","lam_2"};
            //    0        1         2      3       4         5         6     7         8         9        10       11       12       13       14         15         16       17      18      19        20         21       22      23

            // ------- REWRITE CORRECTED QUANTITIES OF ALL OBJECTS ----------------------------
            double corrected_pt, corrected_eta, corrected_phi;
            double erroron_corrected_pt, erroron_corrected_eta, erroron_corrected_phi;
            double mass; // variables in which to temporary store

            // rewrite bJet_Hadronic with corrected quantities
            minuit.GetParameter(0,corrected_pt,erroron_corrected_pt);
            minuit.GetParameter(16,corrected_eta,erroron_corrected_eta);
            minuit.GetParameter(10,corrected_phi,erroron_corrected_phi);
            mass=bJet_Hadronic.M();
            bJet_Hadronic.SetPtEtaPhiM( corrected_pt, corrected_eta, corrected_phi, mass);

            // rewrite bJet_Leptonic with corrected quantities
            minuit.GetParameter(1,corrected_pt,erroron_corrected_pt);
            minuit.GetParameter(17,corrected_eta,erroron_corrected_eta);
            minuit.GetParameter(11,corrected_phi,erroron_corrected_phi);
            mass=bJet_Leptonic.M();
            bJet_Leptonic.SetPtEtaPhiM( corrected_pt, corrected_eta, corrected_phi, mass);

            // rewrite lepton with corrected quantities
            minuit.GetParameter(2,corrected_pt,erroron_corrected_pt);
            minuit.GetParameter(18,corrected_eta,erroron_corrected_eta);
            minuit.GetParameter(12,corrected_phi,erroron_corrected_phi);
            mass=lepton.M();
            lepton.SetPtEtaPhiM( corrected_pt, corrected_eta, corrected_phi, mass);

            // rewrite neutrino with corrected quantities
            minuit.GetParameter(3,corrected_pt,erroron_corrected_pt);
            minuit.GetParameter(19,corrected_eta,erroron_corrected_eta);
            minuit.GetParameter(13,corrected_phi,erroron_corrected_phi);
            mass=neutrinoP4.M();
            neutrinoP4.SetPtEtaPhiM( corrected_pt, corrected_eta, corrected_phi, mass);

            // rewrite lightJets[0] with corrected quantities
            minuit.GetParameter(4,corrected_pt,erroron_corrected_pt);
            minuit.GetParameter(20,corrected_eta,erroron_corrected_eta);
            minuit.GetParameter(14,corrected_phi,erroron_corrected_phi);
            mass=lightJets[0].M();
            lightJets[0].SetPtEtaPhiM( corrected_pt, corrected_eta, corrected_phi, mass);

            // rewrite lightJets[1] with corrected quantities
            minuit.GetParameter(5,corrected_pt,erroron_corrected_pt);
            minuit.GetParameter(21,corrected_eta,erroron_corrected_eta);
            minuit.GetParameter(14,corrected_phi,erroron_corrected_phi);
            mass=lightJets[1].M();
            lightJets[1].SetPtEtaPhiM( corrected_pt, corrected_eta, corrected_phi, mass);


            chisquareMattFit = f_value; // f_value is assigned in last iteration of chistar inside of Migrad ( f_value is a global variable :( )
#endif
#ifdef analyticMattFit_ON
            TMatrixD covarianceMatrix(18,18);
            covarianceMatrix.Zero();
            covarianceMatrix[0][0]=pow(e_b1_x,2);   //  pow(abs(2*0.05* bJet_Hadronic.Px() ),2); // error on bJet_Hadronic px
            covarianceMatrix[1][1]=pow(e_b1_y,2);   //  pow(abs(2*0.05* bJet_Hadronic.Py() ),2); // error on bJet_Hadronic py
            covarianceMatrix[2][2]=pow(e_b1_z,2);   //  pow(abs(2*0.05* bJet_Hadronic.Pz() ),2); // error on bJet_Hadronic pz
            covarianceMatrix[3][3]=pow(e_b2_x,2);   //  pow(abs(2*0.05* bJet_Leptonic.Px() ),2); // error on bJet_Leptonic px
            covarianceMatrix[4][4]=pow(e_b2_y,2);   //  pow(abs(2*0.05* bJet_Leptonic.Px() ),2); // error on bJet_Leptonic py
            covarianceMatrix[5][5]=pow(e_b2_z,2);   //  pow(abs(2*0.05* bJet_Leptonic.Py() ),2); // error on bJet_Leptonic pz
            covarianceMatrix[6][6]=pow(e_l_x,2);    // error on Lepton px
            covarianceMatrix[7][7]=pow(e_l_y,2);    // error on Lepton py
            covarianceMatrix[8][8]=pow(e_l_z,2);    // error on Lepton pz
            covarianceMatrix[9][9]=pow(e_nu_x,2);   //pow(abs(0.3*    neutrinoP4.Px()    ),2); // error on Neutrino px
            covarianceMatrix[10][10]=pow(e_nu_y,2); //pow(abs(0.3*  neutrinoP4.Py()    ),2); // error on Neutrino py
            covarianceMatrix[11][11]=pow(e_nu_z,2); // error on Neutrino pz -> set to be huge
            covarianceMatrix[9][10]=0.;             //ev.e_met_pxpy[0];
            covarianceMatrix[10][9]=0.;             //ev.e_met_pxpy[0];
      //    covarianceMatrix[xxx][xxx]=pow(abs(1.e9*bar_p_z_nu),2);  \\ covariance pxpy neutrino

            covarianceMatrix[12][12]=pow(e_q1_x,2);//  pow(abs(2*0.05* lightJets[0].Px() ),2); // error on lightJet0 px
            covarianceMatrix[13][13]=pow(e_q1_y,2);//  pow(abs(2*0.05* lightJets[0].Py() ),2); // error on lightJet0 py
            covarianceMatrix[14][14]=pow(e_q1_z,2);//  pow(abs(2*0.05* lightJets[0].Pz() ),2); // error on lightJet0 pz
            covarianceMatrix[15][15]=pow(e_q2_x,2);//  pow(abs(2*0.05* lightJets[1].Px() ),2); // error on lightJet1 px
            covarianceMatrix[16][16]=pow(e_q2_y,2);//  pow(abs(2*0.05* lightJets[1].Py() ),2); // error on lightJet1 py
            covarianceMatrix[17][17]=pow(e_q2_z,2);//  pow(abs(2*0.05* lightJets[1].Pz() ),2); // error on lightJet1 pz

/*            covarianceMatrix[0][9]=-1.*ev.e_j_px[0]*ev.e_met_px[0];
            covarianceMatrix[9][0]=-1.*ev.e_j_px[0]*ev.e_met_px[0];
            covarianceMatrix[1][10]=-1.*ev.e_j_py[0]*ev.e_met_py[0];
            covarianceMatrix[10][1]=-1.*ev.e_j_py[0]*ev.e_met_py[0];
            covarianceMatrix[3][9]=-1.*ev.e_j_px[1]*ev.e_met_px[0];
            covarianceMatrix[9][3]=-1.*ev.e_j_px[1]*ev.e_met_px[0];
            covarianceMatrix[4][10]=-1.*ev.e_j_py[1]*ev.e_met_py[0];
            covarianceMatrix[10][4]=-1.*ev.e_j_py[1]*ev.e_met_py[0];
            covarianceMatrix[6][9]=-1.*abs(0.01*lepton.Px())*ev.e_met_px[0];
            covarianceMatrix[9][6]=-1.*abs(0.01*lepton.Px())*ev.e_met_px[0];
            covarianceMatrix[7][10]=-1.*abs(0.01*lepton.Py())*ev.e_met_py[0];
            covarianceMatrix[10][7]=-1.*abs(0.01*lepton.Py())*ev.e_met_py[0];
            covarianceMatrix[12][9]=-1.*ev.e_j_px[2]*ev.e_met_px[0];
            covarianceMatrix[9][12]=-1.*ev.e_j_px[2]*ev.e_met_px[0];
            covarianceMatrix[13][10]=-1.*ev.e_j_py[2]*ev.e_met_py[0];
            covarianceMatrix[10][13]=-1.*ev.e_j_py[2]*ev.e_met_py[0];
            covarianceMatrix[15][9]=-1.*ev.e_j_px[3]*ev.e_met_px[0];
            covarianceMatrix[9][15]=-1.*ev.e_j_px[3]*ev.e_met_px[0];
            covarianceMatrix[16][10]=-1.*ev.e_j_py[3]*ev.e_met_py[0];
            covarianceMatrix[10][16]=-1.*ev.e_j_py[3]*ev.e_met_py[0];
*/
            // !!! LAUNCH KINEMATIC FIT AND RE-WRITE CORRECTED VARIABLES!
            analyticMattFit( &bJet_Hadronic, &bJet_Leptonic, &lepton, &neutrinoP4, &lightJets[0], &lightJets[1], covarianceMatrix, &chisquareAnalyticMattFit );
#endif

            //reconstruction of t and tbar with corrected met
            if (isThadronic==1) {
                t_rec    = bJet_Hadronic + lightJets[0].p4() + lightJets[1].p4();
                tbar_rec = bJet_Leptonic + lepton + neutrinoP4;
            } else {
                t_rec    = bJet_Leptonic + lepton + neutrinoP4;
                tbar_rec = bJet_Hadronic + lightJets[0].p4() + lightJets[1].p4();
            }

            //reconstructed ttbar system
            TLorentzVector ttbarSystem_rec = t_rec + tbar_rec;

            //generated t and tbar
            TLorentzVector t_gen(0,0,0,0);
            TLorentzVector tbar_gen(0,0,0,0);
            TLorentzVector ttbarSystem_gen(0,0,0,0);
            TLorentzVector b_gen(0,0,0,0);
            TLorentzVector bbar_gen(0,0,0,0);
            Bool_t firstB = 1;
            Bool_t firstBbar = 1;
            for (int igen=0; igen<ev.ngtop; igen++) {
                if( ev.gtop_id[igen]==6 ) {
                    t_gen.SetPtEtaPhiM(ev.gtop_pt[igen],ev.gtop_eta[igen],ev.gtop_phi[igen],ev.gtop_m[igen]);
                } else if( ev.gtop_id[igen]==-6 ) {
                    tbar_gen.SetPtEtaPhiM(ev.gtop_pt[igen],ev.gtop_eta[igen],ev.gtop_phi[igen],ev.gtop_m[igen]);
                } else if( ev.gtop_id[igen]==5 && firstB==1 ) {
                    b_gen.SetPtEtaPhiM(ev.gtop_pt[igen],ev.gtop_eta[igen],ev.gtop_phi[igen],ev.gtop_m[igen]);
                    firstB = 0;
                } else if( ev.gtop_id[igen]==-5 && firstBbar==1 ) {
                    bbar_gen.SetPtEtaPhiM(ev.gtop_pt[igen],ev.gtop_eta[igen],ev.gtop_phi[igen],ev.gtop_m[igen]);
                    firstBbar = 0;
                }
            }
            ttbarSystem_gen = t_gen+tbar_gen;

#ifdef HISTOGRAMS_ON
            //control histograms
            ht.fill("nvtx",   ev.nvtx,      plotwgts);
            ht.fill("nbjets", bJets.size(), plotwgts);
            ht.fill("njets",  jets.size(),  plotwgts);

            // create separate histograms for resolution of tops if reconstructed from jets or leptons
            if (isThadronic) {
                ht.fill("mtop_res_hadronic", t_rec.M()-t_gen.M(),   plotwgts);
                ht.fill("mtop_res_leptonic", tbar_rec.M()-tbar_gen.M(),   plotwgts);

                ht.fill("pttop_res_hadronic", t_rec.Pt()-t_gen.Pt(),   plotwgts);
                ht.fill("pttop_res_leptonic", tbar_rec.Pt()-tbar_gen.Pt(),   plotwgts);

                ht.fill("ytop_res_hadronic", t_rec.Rapidity()-t_gen.Rapidity(),   plotwgts);
                ht.fill("ytop_res_leptonic", tbar_rec.Rapidity()-tbar_gen.Rapidity(),   plotwgts);
            } else {
                ht.fill("mtop_res_leptonic", t_rec.M()-t_gen.M(),   plotwgts);
                ht.fill("mtop_res_hadronic", tbar_rec.M()-tbar_gen.M(),   plotwgts);

                ht.fill("pttop_res_leptonic", t_rec.Pt()-t_gen.Pt(),   plotwgts);
                ht.fill("pttop_res_hadronic", tbar_rec.Pt()-tbar_gen.Pt(),   plotwgts);

                ht.fill("ytop_res_leptonic", t_rec.Rapidity()-t_gen.Rapidity(),   plotwgts);
                ht.fill("ytop_res_hadronic", tbar_rec.Rapidity()-tbar_gen.Rapidity(),   plotwgts);
            }

            ht.fill("mttbar_gen", ttbarSystem_gen.M(),   plotwgts);
            ht.fill("mttbar_rec", ttbarSystem_rec.M(),   plotwgts);
            ht.fill("mttbar_res", ttbarSystem_rec.M() - ttbarSystem_gen.M(), plotwgts);

            ht.fill("yttbar_rec", ttbarSystem_rec.Rapidity(),   plotwgts);
            ht.fill("yttbar_gen", ttbarSystem_gen.Rapidity(),   plotwgts);
            ht.fill("yttbar_res", ttbarSystem_rec.Rapidity()-ttbarSystem_gen.Rapidity(),   plotwgts);

            ht.fill("ptttbar_rec", ttbarSystem_rec.Pt(),   plotwgts);
            ht.fill("ptttbar_gen", ttbarSystem_gen.Pt(),   plotwgts);
            ht.fill("ptttbar_res", ttbarSystem_rec.Pt()-ttbarSystem_gen.Pt(),   plotwgts);

            ht.fill("mt_res", t_rec.M() - t_gen.M(), plotwgts);
            ht.fill("mtbar_res", tbar_rec.M() - tbar_gen.M(), plotwgts);
            ht.fill("yt_res", t_rec.Rapidity()-t_gen.Rapidity(),   plotwgts);
            ht.fill("ytbar_res", tbar_rec.Rapidity()-tbar_gen.Rapidity(),   plotwgts);
            ht.fill("ptt_res", t_rec.Pt()-tbar_gen.Pt(),   plotwgts);
            ht.fill("pttbar_res", tbar_rec.Pt()-tbar_gen.Pt(),   plotwgts);

            ht.fill("ht", scalarht , plotwgts);
#endif

            //   --- write the variables which will be written in tree ---
            outVars["t_pt"]=t_rec.Pt();
            outVars["t_eta"]=t_rec.Rapidity();
            outVars["t_phi"]=t_rec.Phi();
            outVars["t_m"]=t_rec.M();
            outVars["t_charge"]=1;
            outVars["t_isHadronic"]= isThadronic;

            outVars["tbar_pt"]=tbar_rec.Pt();
            outVars["tbar_eta"]=tbar_rec.Rapidity();
            outVars["tbar_phi"]=tbar_rec.Phi();
            outVars["tbar_m"]=tbar_rec.M();

            outVars["ttbar_pt"]=ttbarSystem_rec.Pt();
            outVars["ttbar_eta"]=ttbarSystem_rec.Rapidity();
            outVars["ttbar_phi"]=ttbarSystem_rec.Phi();
            outVars["ttbar_m"]=ttbarSystem_rec.M();
            outVars["ttbar_E"]=ttbarSystem_rec.E();

            outVars["gen_t_pt"]=t_gen.Pt();
            outVars["gen_t_eta"]=t_gen.Rapidity();
            outVars["gen_t_phi"]=t_gen.Phi();
            outVars["gen_t_m"]=t_gen.M();

            outVars["gen_tbar_pt"]=tbar_gen.Pt();
            outVars["gen_tbar_eta"]=tbar_gen.Rapidity();
            outVars["gen_tbar_phi"]=tbar_gen.Phi();
            outVars["gen_tbar_m"]=tbar_gen.M();

            outVars["gen_ttbar_pt"]=ttbarSystem_gen.Pt();
            outVars["gen_ttbar_eta"]=ttbarSystem_gen.Rapidity();
            outVars["gen_ttbar_phi"]=ttbarSystem_gen.Phi();
            outVars["gen_ttbar_m"]=ttbarSystem_gen.M();
            outVars["gen_ttbar_E"]=ttbarSystem_gen.E();

            outVars["l_pt"]=leptons[0].Pt();
            outVars["l_eta"]=leptons[0].Rapidity();
            outVars["l_phi"]=leptons[0].Phi();
            outVars["l_m"]=leptons[0].M();
            outVars["l_E"]=leptons[0].E();
            outVars["lepton_isolation"]=ev.l_relIso[0];

            outVars["nu_pt"]=neutrinoP4.Pt();
            outVars["nu_eta"]=neutrinoP4.Rapidity();
            outVars["nu_phi"]=neutrinoP4.Phi();
			
			for (int ift=0; ift<ev.nfwdtrk; ift++) {
				const unsigned short pot_raw_id = ev.fwdtrk_pot[ift];
				if(ev.fwdtrk_method[ift]==1){
				if (pot_raw_id<100){ outVars["p1_xi"]=ev.fwdtrk_xi[ift];}
				else { outVars["p2_xi"] = ev.fwdtrk_xi[ift];}
				}
			}
            //outVars["p1_xi"]=ev.fwdtrk_xi[0];
            //outVars["p2_xi"]=ev.fwdtrk_xi[1];

#ifdef KINFIT_ON
            outVars["nu_pt_uncorrected"] =neutrinoP4_uncorr.Pt();
            outVars["nu_eta_uncorrected"]=neutrinoP4_uncorr.Rapidity();
            outVars["nu_phi_uncorrected"]=neutrinoP4_uncorr.Phi();
            outVars["nu_pt_min"] =neutrinoP4_min->Pt();
            outVars["nu_eta_min"]=neutrinoP4_min->Rapidity();
            outVars["nu_phi_min"]=neutrinoP4_min->Phi();
            outVars["nu_pt_max"] =neutrinoP4_max->Pt();
            outVars["nu_eta_max"]=neutrinoP4_max->Rapidity();
            outVars["nu_phi_max"]=neutrinoP4_max->Phi();
#endif
            outVars["chisquareKinFit"]=chisquareKinFit;
            outVars["chisquareMattFit"]=chisquareMattFit;
            outVars["chisquareAnalyticMattFit"]=chisquareAnalyticMattFit;

            outVars["lightJet0_pt"]=lightJets[0].Pt();
            outVars["lightJet0_eta"]=lightJets[0].Rapidity();
            outVars["lightJet0_phi"]=lightJets[0].Phi();
            outVars["lightJet0_m"]=lightJets[0].M();
            outVars["lightJet0_E"]=lightJets[0].E();
            outVars["lightJet1_pt"]=lightJets[1].Pt();
            outVars["lightJet1_eta"]=lightJets[1].Rapidity();
            outVars["lightJet1_phi"]=lightJets[1].Phi();
            outVars["lightJet1_m"]=lightJets[1].M();
            outVars["lightJet1_E"]=lightJets[1].E();
            if (lightJets.size()==3) {
                outVars["lightJet2_pt"]=  lightJets[2].Pt();
                outVars["lightJet2_eta"]= lightJets[2].Rapidity();
                outVars["lightJet2_phi"]= lightJets[2].Phi();
                outVars["lightJet2_m"]=   lightJets[2].M();
                outVars["lightJet2_E"]=   lightJets[2].E();
                outVars["lightJet3_pt"]=  0. ;
                outVars["lightJet3_eta"]= 0. ;
                outVars["lightJet3_phi"]= 0. ;
                outVars["lightJet3_m"]=   0. ;
                outVars["lightJet3_E"]=   0. ;
            } else if (lightJets.size()>=4) {
                outVars["lightJet2_pt"]=  lightJets[2].Pt();
                outVars["lightJet2_eta"]= lightJets[2].Rapidity();
                outVars["lightJet2_phi"]= lightJets[2].Phi();
                outVars["lightJet2_m"]=   lightJets[2].M();
                outVars["lightJet2_E"]=   lightJets[2].E();
                outVars["lightJet3_pt"]=  lightJets[3].Pt();
                outVars["lightJet3_eta"]= lightJets[3].Rapidity();
                outVars["lightJet3_phi"]= lightJets[3].Phi();
                outVars["lightJet3_m"]=   lightJets[3].M();
                outVars["lightJet3_E"]=   lightJets[3].E();
            } else { // lightJets.size()==2
                outVars["lightJet2_pt"]=  0. ;
                outVars["lightJet2_eta"]= 0. ;
                outVars["lightJet2_phi"]= 0. ;
                outVars["lightJet2_m"]=   0. ;
                outVars["lightJet2_E"]=   0. ;
                outVars["lightJet3_pt"]=  0. ;
                outVars["lightJet3_eta"]= 0. ;
                outVars["lightJet3_phi"]= 0. ;
                outVars["lightJet3_m"]=   0. ;
                outVars["lightJet3_E"]=   0. ;
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

            outVars["bJet0_px"]=bar_p_x_b1; //bJets[0].Px();
            outVars["bJet0_py"]=bar_p_y_b1;
            outVars["bJet0_pz"]=bar_p_z_b1;
            outVars["bJet1_px"]=bar_p_x_b2;
            outVars["bJet1_py"]=bar_p_y_b2;
            outVars["bJet1_pz"]=bar_p_z_b2;
            outVars["lightJet0_px"]=bar_p_x_q1;
            outVars["lightJet0_py"]=bar_p_y_q1;
            outVars["lightJet0_pz"]=bar_p_z_q1;
            outVars["lightJet1_px"]=bar_p_x_q2;
            outVars["lightJet1_py"]=bar_p_y_q2;
            outVars["lightJet1_pz"]=bar_p_z_q2;

#ifdef SAVEERRORS_ON
            outVars["e_l_px"]=0.01*lepton.Px() ;
            outVars["e_l_py"]=0.01*lepton.Py();
            outVars["e_l_pz"]=0.01*lepton.Pz();
            outVars["e_met_px"]=ev.e_met_px;
            outVars["e_met_py"]=ev.e_met_py;
            outVars["e_met_pxpy"]=ev.e_met_pxpy;
            outVars["e_bJet0_px"]=ev.e_j_px[0];
            outVars["e_bJet0_py"]=ev.e_j_py[0];
            outVars["e_bJet0_pz"]=ev.e_j_pz[0];
            outVars["e_bJet1_px"]=ev.e_j_px[1];
            outVars["e_bJet1_py"]=ev.e_j_py[1];
            outVars["e_bJet1_pz"]=ev.e_j_pz[1];
            outVars["e_lightJet0_px"]=ev.e_j_px[2];
            outVars["e_lightJet0_py"]=ev.e_j_py[2];
            outVars["e_lightJet0_pz"]=ev.e_j_pz[2];
            outVars["e_lightJet1_px"]=ev.e_j_px[3];
            outVars["e_lightJet1_py"]=ev.e_j_py[3];
            outVars["e_lightJet1_pz"]=ev.e_j_pz[3];
/*          "e_l_px", "e_l_py", "e_l_pz",
            "e_met_px", "e_met_py", "e_met_pxpy",
            "e_bJet0_px", "e_bJet0_py", "e_bJet0_pz",
            "e_bJet1_px", "e_bJet1_py", "e_bJet1_pz",
            "e_lightJet0_px", "e_lightJet0_py", "e_lightJet0_pz",
            "e_lightJet1_px", "e_lightJet1_py", "e_lightJet1_pz", */
#endif
            outVars["gen_b_pt"]=b_gen.Pt();
            outVars["gen_b_eta"]=b_gen.Rapidity();
            outVars["gen_b_phi"]=b_gen.Phi();
            outVars["gen_b_m"]=b_gen.M();

            outVars["gen_bbar_pt"]=bbar_gen.Pt();
            outVars["gen_bbar_eta"]=bbar_gen.Rapidity();
            outVars["gen_bbar_phi"]=bbar_gen.Phi();
            outVars["gen_bbar_m"]=bbar_gen.M();

            outVars["ht"]=scalarht;
            outVars["nJets"]=jets.size();
            outVars["nBjets"]=bJets.size();
            outVars["nLightJets"]=lightJets.size();

            for(size_t ij=0; ij<jets.size() || ij<20; ij++)  {
                if(jets[ij].flavor()==5) {
                    vBJet_pt.push_back( jets[ij].pt() );
                    vBJet_eta.push_back( jets[ij].eta() );
                } else  {
                    vLightJet_pt.push_back( jets[ij].pt() );
                    vLightJet_eta.push_back( jets[ij].eta() );
                }
            }
            /*           outT->Branch("bJets_pt", &vBJet_pt);
             outT->Branch("bJets_eta", &vBJet_eta);
             outT->Branch("lightJets_pt", &vLightJet_pt);
             outT->Branch("lightJets_eta", &vLightJet_eta);
             */

            // FILL TREE
            outT->Fill();

            // CLEAR MEMORY FROM VECTORS
            vLightJet_pt.clear();
            vLightJet_eta.clear();
            vBJet_pt.clear();
            vBJet_eta.clear();
#ifdef KINFIT_ON
            // CLEAR MEMORY FROM POINTERS
            neutrinoP4_new->Delete();
            neutrinoP4_max->Delete();
            neutrinoP4_min->Delete();
#endif
            /*           vector<double>().swap(vLightJet_pt);
             vector<double>().swap(vLightJet_eta);
             vector<double>().swap(vBJet_pt);
             vector<double>().swap(vBJet_eta);
             */
        } // end of if(bJets.size()>=2 && lightJets.size()>=2)
    } // end of loop over events
    cout << endl;
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
    fOut->Close();
}  // end of RunExclusiveTop()



#ifdef MattFIT_ON
void chistar(int& npar, double *deriv, double& f, double *par, int flag){  // function to be minimized by Minuit
    double parx[6];
    double pary[6];
    double parz[6];

    parx[0]=par[0]*cos(par[10]); //b1
    parx[1]=par[1]*cos(par[11]); //b2
    parx[2]=par[2]*cos(par[12]); //l
    parx[3]=par[3]*cos(par[13]); //nu
    parx[4]=par[4]*cos(par[14]); //q1
    parx[5]=par[5]*cos(par[15]); //q2
    pary[0]=par[0]*sin(par[10]);
    pary[1]=par[1]*sin(par[11]);
    pary[2]=par[2]*sin(par[12]);
    pary[3]=par[3]*sin(par[13]);
    pary[4]=par[4]*sin(par[14]);
    pary[5]=par[5]*sin(par[15]);
    parz[0]=par[0]/tan(theta.Eval(par[16]));
    parz[1]=par[1]/tan(theta.Eval(par[17]));
    parz[2]=par[2]/tan(theta.Eval(par[18]));
    parz[3]=par[3]/tan(theta.Eval(par[19]));
    parz[4]=par[4]/tan(theta.Eval(par[20]));
    parz[5]=par[5]/tan(theta.Eval(par[21]));

    //Calcolo quantità adroniche (figlie top adronico)
    double p_tot_ad=sqrt(pow(parx[0]+parx[4]+parx[5],2)+pow(pary[0]+pary[4]+pary[5],2)+pow(parz[0]+parz[4]+parz[5],2));
    double E_tot_ad= sqrt(parx[0]*parx[0]+pary[0]*pary[0]+parz[0]*parz[0]+m_b*m_b)+sqrt(parx[4]*parx[4]+pary[4]*pary[4]+parz[4]*parz[4]+m_q1*m_q1)+sqrt(parx[5]*parx[5]+pary[5]*pary[5]+parz[5]*parz[5]+m_q2*m_q2);
    double m_inv_ad=sqrt(E_tot_ad*E_tot_ad-p_tot_ad*p_tot_ad);

    //Calcolo quantità leptoniche (figlie top leptonico)
    double p_tot_lep=sqrt(pow(parx[1]+parx[2]+parx[3],2)+pow(pary[1]+pary[2]+pary[3],2)+pow(parz[1]+parz[2]+parz[3],2));
    double E_tot_lep= sqrt(parx[1]*parx[1]+pary[1]*pary[1]+parz[1]*parz[1]+m_b*m_b)+sqrt(parx[2]*parx[2]+pary[2]*pary[2]+parz[2]*parz[2]+m_l*m_l)+sqrt(parx[3]*parx[3]+pary[3]*pary[3]+parz[3]*parz[3]+0.*0.);
    double m_inv_lep=sqrt(E_tot_lep*E_tot_lep-p_tot_lep*p_tot_lep);

    //Calcolo W adronico
    double p_W_ad = sqrt(pow(parx[4]+parx[5],2)+pow(pary[4]+pary[5],2)+pow(parz[4]+parz[5],2));
    double E_W_ad=sqrt(parx[4]*parx[4]+pary[4]*pary[4]+parz[4]*parz[4]+m_q1*m_q1)+sqrt(parx[5]*parx[5]+pary[5]*pary[5]+parz[5]*parz[5]+m_q2*m_q2);
    double m_inv_W_ad=sqrt(E_W_ad*E_W_ad-p_W_ad*p_W_ad);

    //Calcolo W leptonico
    double p_W_lep = sqrt(pow(parx[2]+parx[3],2)+pow(pary[2]+pary[3],2)+pow(parz[2]+parz[3],2));
    double E_W_lep= sqrt(parx[2]*parx[2]+pary[2]*pary[2]+parz[2]*parz[2]+m_l*m_l)+sqrt(parx[3]*parx[3]+pary[3]*pary[3]+parz[3]*parz[3]+0.*0.);
    double m_inv_W_lep=sqrt(E_W_lep*E_W_lep-p_W_lep*p_W_lep);

    f=pow(m_inv_lep-par[6],2)/pow(1.3,2)+pow(m_inv_ad-par[7],2)/pow(1.3,2)+pow(m_inv_W_lep-par[9],2)/pow(2.082,2)+pow(m_inv_W_ad-par[8],2)/pow(2.082,2); //f=chistar parte su top e W
    f+=pow(pt_b1-par[0],2)/pow(e_b1,2);
    f+=pow(pt_b2-par[1],2)/pow(e_b2,2);
    f+=pow(pt_l-par[2],2)/pow(e_l,2);
    f+=pow(pt_nu-par[3],2)/pow(e_nu,2);
    f+=pow(pt_q1-par[4],2)/pow(e_q1,2);
    f+=pow(pt_q2-par[5],2)/pow(e_q2,2);
    f+=pow(m_t-par[6],2)/pow(e_t,2);
    f+=pow(m_tbar-par[7],2)/pow(e_tbar,2);
    f+=pow(m_W_ad-par[8],2)/pow(e_W_ad,2);
    f+=pow(m_W_lep-par[9],2)/pow(e_W_lep,2);
    f+=pow(phi_b1-par[10],2)/pow(e_phi_b1,2);
    f+=pow(phi_b2-par[11],2)/pow(e_phi_b2,2);
    f+=pow(phi_l-par[12],2)/pow(e_phi_l,2);
    f+=pow(phi_nu-par[13],2)/pow(e_phi_nu,2);
    f+=pow(phi_q1-par[14],2)/pow(e_phi_q1,2);
    f+=pow(phi_q2-par[15],2)/pow(e_phi_q2,2);
    f+=pow(eta_b1-par[16],2)/pow(e_eta_b1,2);
    f+=pow(eta_b2-par[17],2)/pow(e_eta_b2,2);
    f+=pow(eta_l-par[18],2)/pow(e_eta_l,2);
    f+=pow(eta_nu-par[19],2)/pow(e_eta_nu,2);
    f+=pow(eta_q1-par[20],2)/pow(e_eta_q1,2);
    f+=pow(eta_q2-par[21],2)/pow(e_eta_q2,2);

    f_value = f;

    //Parte sui moltiplicatori di lagrange
    H1.SetParameter(0,parx[0]);
    H1.SetParameter(1,parx[1]);
    H1.SetParameter(2,parx[2]);
    H1.SetParameter(3,parx[3]);
    H1.SetParameter(4,parx[4]);
    H1.SetParameter(5,parx[5]);
    H2.SetParameter(0,pary[0]);
    H2.SetParameter(1,pary[1]);
    H2.SetParameter(2,pary[2]);
    H2.SetParameter(3,pary[3]);
    H2.SetParameter(4,pary[4]);
    H2.SetParameter(5,pary[5]);

    double d1=H1.Eval(0);
    double d2=H2.Eval(0);

    f+=abs(2*par[22]*(d1)); //Aggiungo il moltiplicatore di lagrange x
    f+=abs(2*par[23]*(d2)); //Aggiungo il moltiplicatore di lagrange y

//    cout << "chisquare dentro chistar = "<< f << endl;
//    f_value = f;
}
#endif

#ifdef analyticMattFit_ON
void analyticMattFit (TLorentzVector* bJet_Had, TLorentzVector* bJet_Lep, TLorentzVector* lep, TLorentzVector* nu, TLorentzVector* lightJet0, TLorentzVector* lightJet1, TMatrixD covarianceMatrix, double* chiSquare) {

    double phi_b1= bJet_Had->Phi(); // bJet1 is hadronic
    double phi_b2= bJet_Lep->Phi();
    double phi_l=  lep->Phi();
    double phi_nu= nu->Phi();
    double phi_q1= lightJet0->Phi();
    double phi_q2= lightJet1->Phi();
    double eta_b1= bJet_Had->Rapidity();
    double eta_b2= bJet_Lep->Rapidity();
    double eta_l=  lep->Rapidity();
    double eta_nu= nu->Rapidity();
    double eta_q1= lightJet0->Rapidity();
    double eta_q2= lightJet1->Rapidity();
    double pt_b1=  bJet_Had->Pt();
    double pt_b2=  bJet_Lep->Pt();
    double pt_l=   lep->Pt();
    double pt_nu=  nu->Pt();
    double pt_q1=  lightJet0->Pt();
    double pt_q2=  lightJet1->Pt();

    //    double pt_t_gen=tree->GetLeaf("gen_t_ptF")->GetValue(0);
    //    double pt_tbar_gen=tree->GetLeaf("gen_tbar_ptF")->GetValue(0);
    //    double m_ttbar_gen=tree->GetLeaf("gen_ttbar_mF")->GetValue(0);

    //   int is_h= tree->GetLeaf("t_isHadronicF")->GetValue(0);

    TF1 theta("theta", "2*atan(exp(-x))",0, 1e10);

    ////Variabili in input (nomi definiti nel pdf di riferimento
    double bar_p_x_b1= pt_b1*cos(phi_b1);   //b1 è adronico
    double bar_p_y_b1= pt_b1*sin(phi_b1);
    double bar_p_z_b1= pt_b1/tan(theta.Eval(eta_b1));
    double bar_p_x_b2= pt_b2*cos(phi_b2);   //b2 è leptonico
    double bar_p_y_b2= pt_b2*sin(phi_b2);
    double bar_p_z_b2= pt_b2/tan(theta.Eval(eta_b2));
    double bar_p_x_l= pt_l*cos(phi_l);
    double bar_p_y_l= pt_l*sin(phi_l);
    double bar_p_z_l= pt_l/tan(theta.Eval(eta_l));
    double bar_p_x_nu= pt_nu*cos(phi_nu);
    double bar_p_y_nu= pt_nu*sin(phi_nu);
    double bar_p_z_nu= pt_nu/tan(theta.Eval(eta_nu));
    double bar_p_x_q1= pt_q1*cos(phi_q1);
    double bar_p_y_q1= pt_q1*sin(phi_q1);;
    double bar_p_z_q1= pt_q1/tan(theta.Eval(eta_q1));
    double bar_p_x_q2= pt_q2*cos(phi_q2);
    double bar_p_y_q2= pt_q2*cos(phi_q2);
    double bar_p_z_q2= pt_q2/tan(theta.Eval(eta_q2));

    double constraint_violation_1=0.;
    double constraint_violation_2=0.;
    vector<double> constraint_violation_1_vector;
    vector<double> constraint_violation_2_vector;
    vector<double> sum_squared;
    vector<double> dumping;

    int maxIter=15;//15;    // Setta il numero di iterazioni prima dello stop
    int maxStep=5000;       // Setta il numero di step in cui si suddivide lo spostamento verso il punto di convergenza
    double m_t = m_TOP; //Massa del top

    ////Definisco il vettore di traslazione
    TMatrixD f(24,1);
    f.Zero();

    for(int iteration=0; iteration<maxIter; iteration++){
        double dump=1.;
        int step=0;
        int escape=0;

        TMatrixD C(18,18); //Matrice di covarianza
        C.Zero();
        C = covarianceMatrix;

        //// Definisco le matrici trattate nel pdf di riferimento
        TMatrixD C_star(24,24);
        TMatrixD K_star(24,24);
        TMatrixD P_star(24,24);
        TMatrixD Q_star(24,1);
        C_star.Zero();  //Riempio di zeri tutte le matrici.
        K_star.Zero();
        P_star.Zero();
        Q_star.Zero();

        //Insert in matrix only


        ////Definisco il vettore contenente le quantità da fittare;
        TMatrixD x(24,1);
        TMatrixD fittedsolutions(24,1);
        x.Zero();
        fittedsolutions.Zero();

        //// Riempio gli oggetti con gli input
        TMatrixD InvC(18,18);
        // cout << "ciao" << endl;
        InvC=C.Invert();
        // cout << "ciao" << endl;

        //Riempizione C_star
        for(int row=0; row<18; row++){
            for(int col=0; col<18; col++){
                C_star[row][col]=InvC[row][col];
            }
        }

        //Riempizione K_star
        K_star[0][18]=1.;
        K_star[1][19]=1.;
        K_star[3][18]=1.;
        K_star[4][19]=1.;
        K_star[6][18]=1.;
        K_star[7][19]=1.;
        K_star[9][18]=1.;
        K_star[10][19]=1.;
        K_star[12][18]=1.;
        K_star[13][19]=1.;
        K_star[15][18]=1.;
        K_star[16][19]=1.;
        K_star[18][0]=1.;
        K_star[19][1]=1.;
        K_star[18][3]=1.;
        K_star[19][4]=1.;
        K_star[18][6]=1.;
        K_star[19][7]=1.;
        K_star[18][9]=1.;
        K_star[19][10]=1.;
        K_star[18][12]=1.;
        K_star[19][13]=1.;
        K_star[18][15]=1.;
        K_star[19][16]=1.;

        //Riempizione Q_star
        //Per calcolare le entrate di Q_star devo prima calcolare
        //i valori della massa invariante dei sistemi adronici
        //e leptonici usando i valori di p non fittati.
        double mod_bar_p_b1=sqrt(pow(bar_p_x_b1,2)+pow(bar_p_y_b1,2)+pow(bar_p_z_b1,2)); //Moduli trivettori p
        double mod_bar_p_b2=sqrt(pow(bar_p_x_b2,2)+pow(bar_p_y_b2,2)+pow(bar_p_z_b2,2));
        double mod_bar_p_l=sqrt(pow(bar_p_x_l,2)+pow(bar_p_y_l,2)+pow(bar_p_z_l,2));
        double mod_bar_p_nu=sqrt(pow(bar_p_x_nu,2)+pow(bar_p_y_nu,2)+pow(bar_p_z_nu,2));
        double mod_bar_p_q1=sqrt(pow(bar_p_x_q1,2)+pow(bar_p_y_q1,2)+pow(bar_p_z_q1,2));
        double mod_bar_p_q2=sqrt(pow(bar_p_x_q2,2)+pow(bar_p_y_q2,2)+pow(bar_p_z_q2,2));

        double bar_m_lep = sqrt(2*mod_bar_p_b2*mod_bar_p_l+2*mod_bar_p_b2*mod_bar_p_nu+2*mod_bar_p_l*mod_bar_p_nu-2*bar_p_x_b2*bar_p_x_l-2*bar_p_x_b2*bar_p_x_nu-2*bar_p_x_l*bar_p_x_nu-2*bar_p_y_b2*bar_p_y_l-2*bar_p_y_b2*bar_p_y_nu-2*bar_p_y_l*bar_p_y_nu-2*bar_p_z_b2*bar_p_z_l-2*bar_p_z_b2*bar_p_z_nu-2*bar_p_z_l*bar_p_z_nu);
        double bar_m_ad = sqrt(2*mod_bar_p_b1*mod_bar_p_q1+2*mod_bar_p_b1*mod_bar_p_q2+2*mod_bar_p_q1*mod_bar_p_q2-2*bar_p_x_b1*bar_p_x_q1-2*bar_p_x_b1*bar_p_x_q2-2*bar_p_x_q1*bar_p_x_q2-2*bar_p_y_b1*bar_p_y_q1-2*bar_p_y_b1*bar_p_y_q2-2*bar_p_y_q1*bar_p_y_q2-2*bar_p_z_b1*bar_p_z_q1-2*bar_p_z_b1*bar_p_z_q2-2*bar_p_z_q1*bar_p_z_q2);
        double bar_m_lepnu = sqrt(2*mod_bar_p_l*mod_bar_p_nu-2*bar_p_x_l*bar_p_x_nu-2*bar_p_y_l*bar_p_y_nu-2*bar_p_z_l*bar_p_z_nu);
        double bar_m_q1q2 = sqrt(2*mod_bar_p_q1*mod_bar_p_q2-2*bar_p_x_q1*bar_p_x_q2-2*bar_p_y_q1*bar_p_y_q2-2*bar_p_z_q1*bar_p_z_q2);

        Q_star[20][0]=pow(bar_m_lep,2)-m_t*m_t;
        Q_star[21][0]=pow(bar_m_ad,2)-m_t*m_t;
        Q_star[22][0]=pow(bar_m_lepnu,2)-m_W*m_W;
        Q_star[23][0]=pow(bar_m_q1q2,2)-m_W*m_W;

        //Riempizione f
        if(iteration==0){
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
        }

        //Riempizione P_star
        P_star[3][20]=2*(mod_bar_p_l*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_x_b2/mod_bar_p_b2-bar_p_x_l-bar_p_x_nu);
        P_star[4][20]=2*(mod_bar_p_l*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_y_b2/mod_bar_p_b2-bar_p_y_l-bar_p_y_nu);
        P_star[5][20]=2*(mod_bar_p_l*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_z_b2/mod_bar_p_b2-bar_p_z_l-bar_p_z_nu);
        P_star[6][20]=2*(mod_bar_p_b2*bar_p_x_l/mod_bar_p_l+mod_bar_p_nu*bar_p_x_l/mod_bar_p_l-bar_p_x_b2-bar_p_x_nu);
        P_star[7][20]=2*(mod_bar_p_b2*bar_p_y_l/mod_bar_p_l+mod_bar_p_nu*bar_p_y_l/mod_bar_p_l-bar_p_y_b2-bar_p_y_nu);
        P_star[8][20]=2*(mod_bar_p_b2*bar_p_z_l/mod_bar_p_l+mod_bar_p_nu*bar_p_z_l/mod_bar_p_l-bar_p_z_b2-bar_p_z_nu);
        P_star[9][20]=2*(mod_bar_p_b2*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu-bar_p_x_b2-bar_p_x_l);
        P_star[10][20]=2*(mod_bar_p_b2*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu-bar_p_y_b2-bar_p_y_l);
        P_star[11][20]=2*(mod_bar_p_b2*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu-bar_p_z_b2-bar_p_z_l);
        P_star[0][21]=2*(mod_bar_p_q1*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_x_b1/mod_bar_p_b1-bar_p_x_q1-bar_p_x_q2);
        P_star[1][21]=2*(mod_bar_p_q1*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_y_b1/mod_bar_p_b1-bar_p_y_q1-bar_p_y_q2);
        P_star[2][21]=2*(mod_bar_p_q1*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_z_b1/mod_bar_p_b1-bar_p_z_q1-bar_p_z_q2);
        P_star[12][21]=2*(mod_bar_p_b1*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_b1-bar_p_x_q2);
        P_star[13][21]=2*(mod_bar_p_b1*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_b1-bar_p_y_q2);
        P_star[14][21]=2*(mod_bar_p_b1*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_b1-bar_p_z_q2);
        P_star[15][21]=2*(mod_bar_p_b1*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_x_q2/mod_bar_p_q2-bar_p_x_b1-bar_p_x_q1);
        P_star[16][21]=2*(mod_bar_p_b1*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_y_q2/mod_bar_p_q2-bar_p_y_b1-bar_p_y_q1);
        P_star[17][21]=2*(mod_bar_p_b1*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_z_q2/mod_bar_p_q2-bar_p_z_b1-bar_p_z_q1);
        P_star[6][22]=2*(mod_bar_p_nu*bar_p_x_l/mod_bar_p_l-bar_p_x_nu);
        P_star[7][22]=2*(mod_bar_p_nu*bar_p_y_l/mod_bar_p_l-bar_p_y_nu);
        P_star[8][22]=2*(mod_bar_p_nu*bar_p_z_l/mod_bar_p_l-bar_p_z_nu);
        P_star[9][22]=2*(mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu-bar_p_x_l);
        P_star[10][22]=2*(mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu-bar_p_y_l);
        P_star[11][22]=2*(mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu-bar_p_z_l);
        P_star[12][23]=2*(mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_q2);
        P_star[13][23]=2*(mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_q2);
        P_star[14][23]=2*(mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_q2);
        P_star[15][23]=2*(mod_bar_p_q1*bar_p_x_q2/mod_bar_p_q2-bar_p_x_q1);
        P_star[16][23]=2*(mod_bar_p_q1*bar_p_y_q2/mod_bar_p_q2-bar_p_y_q1);
        P_star[17][23]=2*(mod_bar_p_q1*bar_p_z_q2/mod_bar_p_q2-bar_p_z_q1);

        P_star[20][3]=2*(mod_bar_p_l*bar_p_x_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_x_b2/mod_bar_p_b2-bar_p_x_l-bar_p_x_nu); // symmetric matrix
        P_star[20][4]=2*(mod_bar_p_l*bar_p_y_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_y_b2/mod_bar_p_b2-bar_p_y_l-bar_p_y_nu);
        P_star[20][5]=2*(mod_bar_p_l*bar_p_z_b2/mod_bar_p_b2+mod_bar_p_nu*bar_p_z_b2/mod_bar_p_b2-bar_p_z_l-bar_p_z_nu);
        P_star[20][6]=2*(mod_bar_p_b2*bar_p_x_l/mod_bar_p_l+mod_bar_p_nu*bar_p_x_l/mod_bar_p_l-bar_p_x_b2-bar_p_x_nu);
        P_star[20][7]=2*(mod_bar_p_b2*bar_p_y_l/mod_bar_p_l+mod_bar_p_nu*bar_p_y_l/mod_bar_p_l-bar_p_y_b2-bar_p_y_nu);
        P_star[20][8]=2*(mod_bar_p_b2*bar_p_z_l/mod_bar_p_l+mod_bar_p_nu*bar_p_z_l/mod_bar_p_l-bar_p_z_b2-bar_p_z_nu);
        P_star[20][9]=2*(mod_bar_p_b2*bar_p_x_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu-bar_p_x_b2-bar_p_x_l);
        P_star[20][10]=2*(mod_bar_p_b2*bar_p_y_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu-bar_p_y_b2-bar_p_y_l);
        P_star[20][11]=2*(mod_bar_p_b2*bar_p_z_nu/mod_bar_p_nu+mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu-bar_p_z_b2-bar_p_z_l);
        P_star[21][0]=2*(mod_bar_p_q1*bar_p_x_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_x_b1/mod_bar_p_b1-bar_p_x_q1-bar_p_x_q2);
        P_star[21][1]=2*(mod_bar_p_q1*bar_p_y_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_y_b1/mod_bar_p_b1-bar_p_y_q1-bar_p_y_q2);
        P_star[21][2]=2*(mod_bar_p_q1*bar_p_z_b1/mod_bar_p_b1+mod_bar_p_q2*bar_p_z_b1/mod_bar_p_b1-bar_p_z_q1-bar_p_z_q2);
        P_star[21][12]=2*(mod_bar_p_b1*bar_p_x_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_b1-bar_p_x_q2);
        P_star[21][13]=2*(mod_bar_p_b1*bar_p_y_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_b1-bar_p_y_q2);
        P_star[21][14]=2*(mod_bar_p_b1*bar_p_z_q1/mod_bar_p_q1+mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_b1-bar_p_z_q2);
        P_star[21][15]=2*(mod_bar_p_b1*bar_p_x_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_x_q2/mod_bar_p_q2-bar_p_x_b1-bar_p_x_q1);
        P_star[21][16]=2*(mod_bar_p_b1*bar_p_y_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_y_q2/mod_bar_p_q2-bar_p_y_b1-bar_p_y_q1);
        P_star[21][17]=2*(mod_bar_p_b1*bar_p_z_q2/mod_bar_p_q2+mod_bar_p_q1*bar_p_z_q2/mod_bar_p_q2-bar_p_z_b1-bar_p_z_q1);
        P_star[22][6]=2*(mod_bar_p_nu*bar_p_x_l/mod_bar_p_l-bar_p_x_nu);
        P_star[22][7]=2*(mod_bar_p_nu*bar_p_y_l/mod_bar_p_l-bar_p_y_nu);
        P_star[22][8]=2*(mod_bar_p_nu*bar_p_z_l/mod_bar_p_l-bar_p_z_nu);
        P_star[22][9]=2*(mod_bar_p_l*bar_p_x_nu/mod_bar_p_nu-bar_p_x_l);
        P_star[22][10]=2*(mod_bar_p_l*bar_p_y_nu/mod_bar_p_nu-bar_p_y_l);
        P_star[22][11]=2*(mod_bar_p_l*bar_p_z_nu/mod_bar_p_nu-bar_p_z_l);
        P_star[23][12]=2*(mod_bar_p_q2*bar_p_x_q1/mod_bar_p_q1-bar_p_x_q2);
        P_star[23][13]=2*(mod_bar_p_q2*bar_p_y_q1/mod_bar_p_q1-bar_p_y_q2);
        P_star[23][14]=2*(mod_bar_p_q2*bar_p_z_q1/mod_bar_p_q1-bar_p_z_q2);
        P_star[23][15]=2*(mod_bar_p_q1*bar_p_x_q2/mod_bar_p_q2-bar_p_x_q1);
        P_star[23][16]=2*(mod_bar_p_q1*bar_p_y_q2/mod_bar_p_q2-bar_p_y_q1);
        P_star[23][17]=2*(mod_bar_p_q1*bar_p_z_q2/mod_bar_p_q2-bar_p_z_q1);

        double p_x_t;
        double p_x_tbar;
        double p_y_t;
        double p_y_tbar;
        double p_z_t;
        double p_z_tbar;
        double E_t;
        double E_tbar;
        double p_t;
        double p_tbar;
        double m_inv_t;

        //// Determino il valore di x --> correction on quantities
        x=-1.*((C_star+K_star+P_star).Invert())*0.5*(2.*K_star*f+Q_star);

        ///////Verifica dei vincoli cinematici (step method).
        while(1){
            fittedsolutions=dump*x+f;

            p_x_t=fittedsolutions[0][0]+fittedsolutions[6][0]+fittedsolutions[9][0];        //Costruisco quantità per il calcolo della massa inv
            p_x_tbar=fittedsolutions[3][0]+fittedsolutions[12][0]+fittedsolutions[15][0];
            p_y_t=fittedsolutions[1][0]+fittedsolutions[7][0]+fittedsolutions[10][0];
            p_y_tbar=fittedsolutions[4][0]+fittedsolutions[13][0]+fittedsolutions[16][0];
            p_z_t=fittedsolutions[2][0]+fittedsolutions[8][0]+fittedsolutions[11][0];
            p_z_tbar=fittedsolutions[5][0]+fittedsolutions[14][0]+fittedsolutions[17][0];

            E_t=sqrt(fittedsolutions[0][0]*fittedsolutions[0][0]+fittedsolutions[1][0]*fittedsolutions[1][0]+fittedsolutions[2][0]*fittedsolutions[2][0])+sqrt(fittedsolutions[6][0]*fittedsolutions[6][0]+fittedsolutions[7][0]*fittedsolutions[7][0]+fittedsolutions[8][0]*fittedsolutions[8][0])+sqrt(fittedsolutions[9][0]*fittedsolutions[9][0]+fittedsolutions[10][0]*fittedsolutions[10][0]+fittedsolutions[11][0]*fittedsolutions[11][0]);

            E_tbar=sqrt(fittedsolutions[3][0]*fittedsolutions[3][0]+fittedsolutions[4][0]*fittedsolutions[4][0]+fittedsolutions[5][0]*fittedsolutions[5][0])+sqrt(fittedsolutions[12][0]*fittedsolutions[12][0]+fittedsolutions[13][0]*fittedsolutions[13][0]+fittedsolutions[14][0]*fittedsolutions[14][0])+sqrt(fittedsolutions[15][0]*fittedsolutions[15][0]+fittedsolutions[16][0]*fittedsolutions[16][0]+fittedsolutions[17][0]*fittedsolutions[17][0]);

            p_t=sqrt(p_x_t*p_x_t+p_y_t*p_y_t+p_z_t*p_z_t);
            p_tbar=sqrt(p_x_tbar*p_x_tbar+p_y_tbar*p_y_tbar+p_z_tbar*p_z_tbar);

            m_inv_t = sqrt(-p_t*p_t+E_t*E_t);
            double m_inv_tbar = sqrt(-p_tbar*p_tbar+E_tbar*E_tbar);

            constraint_violation_1= abs(m_inv_t-m_TOP);      //Calcolo della violazione del vincolo sulla massa invariante del sistema adronico.
            constraint_violation_2= abs(m_inv_tbar-m_TOP);   //Calcolo della violazione del vincolo sulla massa invariante del sistema leptonico.

            constraint_violation_1_vector.push_back(constraint_violation_1);
            constraint_violation_2_vector.push_back(constraint_violation_2);
            sum_squared.push_back(sqrt(constraint_violation_1*constraint_violation_1+0*constraint_violation_2*constraint_violation_2)); //Somma in quadratura dell'errore sui due vincoli.
            dumping.push_back(dump);

            dump=2*(1.-double (step/(double(maxStep)) ));                //Selezione dello step.

            if(escape==1){
                escape=0;
                constraint_violation_1_vector.clear();
                constraint_violation_2_vector.clear();
                sum_squared.clear();
                dumping.clear();
                break;
            }

            if(step==maxStep){
                //cerco il minimo della somma in quadratura dei vincoli e uso minimo valore per dumping
                dump=dumping[dump_index(sum_squared,maxStep+1)];
                escape=1;
            }

            step++;
        }

        //fittedsolutions=0.5*x+f; // 0.5 is the dump on the correction in order to make the correction converge

        TMatrixD xT(1,24);
        xT.Transpose(x);

        TMatrixD chi2(1,1);
        chi2=xT*C_star*x;

        bar_p_x_b1= fittedsolutions[0][0];   //b1 è adronico
        bar_p_y_b1= fittedsolutions[1][0];
        bar_p_z_b1= fittedsolutions[2][0];
        bar_p_x_b2= fittedsolutions[3][0];   //b2 è leptonico
        bar_p_y_b2= fittedsolutions[4][0];
        bar_p_z_b2= fittedsolutions[5][0];
        bar_p_x_l=  fittedsolutions[6][0];
        bar_p_y_l=  fittedsolutions[7][0];
        bar_p_z_l=  fittedsolutions[8][0];
        bar_p_x_nu= fittedsolutions[9][0];
        bar_p_y_nu= fittedsolutions[10][0];
        bar_p_z_nu= fittedsolutions[11][0];
        bar_p_x_q1= fittedsolutions[12][0];
        bar_p_y_q1= fittedsolutions[13][0];
        bar_p_z_q1= fittedsolutions[14][0];
        bar_p_x_q2= fittedsolutions[15][0];
        bar_p_y_q2= fittedsolutions[16][0];
        bar_p_z_q2= fittedsolutions[17][0];

        *chiSquare = chi2[0][0];

        /* ===============================================================================
        // =========================       DO NOT REWRITE     ============================
        // =============================================================================== */
        if( iteration==(maxIter-1) ) { // rewrite the particles at the last iteration
            double m = 0;
            m=bJet_Had->M();
            bJet_Had->SetPxPyPzE(fittedsolutions[0][0],fittedsolutions[1][0],fittedsolutions[2][0], sqrt(fittedsolutions[0][0]*fittedsolutions[0][0] +fittedsolutions[1][0]*fittedsolutions[1][0]+ fittedsolutions[2][0]*fittedsolutions[2][0] + m*m));

            m=bJet_Lep->M();
            bJet_Lep->SetPxPyPzE(fittedsolutions[3][0],fittedsolutions[4][0],fittedsolutions[5][0], sqrt(fittedsolutions[3][0]*fittedsolutions[3][0] +fittedsolutions[4][0]*fittedsolutions[4][0]+ fittedsolutions[5][0]*fittedsolutions[5][0] + m*m));

            m=lep->M();
            lep->SetPxPyPzE(fittedsolutions[6][0],fittedsolutions[7][0],fittedsolutions[8][0], sqrt(fittedsolutions[6][0]*fittedsolutions[6][0] +fittedsolutions[7][0]*fittedsolutions[7][0]+ fittedsolutions[8][0]*fittedsolutions[8][0] + m*m));

            m=0;
            nu->SetPxPyPzE(fittedsolutions[9][0],fittedsolutions[10][0],fittedsolutions[11][0], sqrt(fittedsolutions[9][0]*fittedsolutions[9][0] +fittedsolutions[10][0]*fittedsolutions[10][0]+ fittedsolutions[11][0]*fittedsolutions[11][0] + m*m));

            m=lightJet0->M();
            lightJet0->SetPxPyPzE(fittedsolutions[12][0],fittedsolutions[13][0],fittedsolutions[14][0], sqrt(fittedsolutions[12][0]*fittedsolutions[12][0] +fittedsolutions[13][0]*fittedsolutions[13][0]+ fittedsolutions[14][0]*fittedsolutions[14][0] + m*m));

            m=lightJet1->M();
            lightJet1->SetPxPyPzE(fittedsolutions[15][0],fittedsolutions[16][0],fittedsolutions[17][0], sqrt(fittedsolutions[15][0]*fittedsolutions[15][0] +fittedsolutions[16][0]*fittedsolutions[16][0]+ fittedsolutions[17][0]*fittedsolutions[17][0] + m*m));
        } // */
    }
}
#endif


// --- THAT'S ALL FOLKS ---



