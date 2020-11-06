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
					 std::string systVar)
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
	bool isData = filename.Contains("Data13TeV");
	bool isPythia8 = filename.Contains("pythia8");
    
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
    MEzCalculator neutrinoPzComputer;
    
    //LUMINOSITY+PILEUP
    LumiTools lumi(era,genPU);
    
    //LEPTON EFFICIENCIES
    EfficiencyScaleFactorsWrapper lepEffH(isData,era);
    
    //B-TAG CALIBRATION
    //BTagSFUtil btvSF(era,"DeepCSV",BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
    BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
	
    //JEC/JER
    JECTools jec(era);
    
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
	
    TString fvars[]={
        "t_pt","t_eta", "t_phi", "t_m", "t_charge", "t_isHadronic",
        "tbar_pt","tbar_eta", "tbar_phi", "tbar_m",
        "ttbar_pt","ttbar_eta", "ttbar_phi", "ttbar_m", "ttbar_E",
        
        // quantities of objects used to build the ttbar
        "l_px", "l_py", "l_pz",
        "nu_px", "nu_py", "nu_pz", "nu_complex",
        "bJet_had_px","bJet_had_py", "bJet_had_pz",
        "bJet_lep_px","bJet_lep_py", "bJet_lep_pz",
        "lightJet0_px", "lightJet0_py", "lightJet0_pz",
        "lightJet1_px", "lightJet1_py", "lightJet1_pz",
        
#ifdef SAVEERRORS_ON
        "e_l_px", "e_l_py", "e_l_pz",
        "e_nu_px", "e_nu_py", "e_nu_pz", "e_nu_pxpy",
        "e_bJet_had_px", "e_bJet_had_py", "e_bJet_had_pz",
        "e_bJet_lep_px", "e_bJet_lep_py", "e_bJet_lep_pz",
        "e_lightJet0_px", "e_lightJet0_py", "e_lightJet0_pz",
        "e_lightJet1_px", "e_lightJet1_py", "e_lightJet1_pz",
#endif
        
#ifdef DEBUG_ON
        // test variables
        "test1", "test2", "test3", "test4", "test5",
#endif
        
        // quantities for all objects in event
        "nBjets", "nLightJets", "nJets", "ht", "cat",

        "l_pt", "l_eta", "l_phi", "l_m", "l_E", "lepton_isolation",
        "nu_pt", "nu_eta", "nu_phi",
        "p1_xi", "p2_xi",
		"weight", "gen_wgt", "toppt_wgt", "selSF_wgt", "trigSF_wgt",
		"selSF_wgt_err", "trigSF_wgt_err", "pu_wgt", "ptag_wgt", "ptag_wgt_err",
		"ren_err","fac_err","scale_err",
		"isr_err","fsr_err",

        "bJet0_pt","bJet0_eta", "bJet0_phi", "bJet0_m", "bJet0_E",
        "bJet1_pt","bJet1_eta", "bJet1_phi", "bJet1_m", "bJet1_E",
        "bJet2_pt","bJet2_eta", "bJet2_phi", "bJet2_m", "bJet2_E",
        "bJet3_pt","bJet3_eta", "bJet3_phi", "bJet3_m", "bJet3_E",
        
        "lightJet0_pt", "lightJet0_eta", "lightJet0_phi", "lightJet0_m", "lightJet0_E",
        "lightJet1_pt", "lightJet1_eta", "lightJet1_phi", "lightJet1_m", "lightJet1_E",
        "lightJet2_pt", "lightJet2_eta", "lightJet2_phi", "lightJet2_m", "lightJet2_E",
        "lightJet3_pt", "lightJet3_eta", "lightJet3_phi", "lightJet3_m", "lightJet3_E",
        
        // quantities of generated objects (MC truth)
        "gen_ttbar_pt","gen_ttbar_eta", "gen_ttbar_phi", "gen_ttbar_m", "gen_ttbar_E",
        "gen_t_pt","gen_t_eta", "gen_t_phi", "gen_t_m",
        "gen_tbar_pt","gen_tbar_eta", "gen_tbar_phi", "gen_tbar_m",
        "gen_b_pt","gen_b_eta", "gen_b_phi", "gen_b_m",
        "gen_bbar_pt","gen_bbar_eta", "gen_bbar_phi", "gen_bbar_m"
        };
    std::map<TString,Float_t> outVars;
    for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++){
        outVars[fvars[i]]=0.;
        ADDVAR(&(outVars[fvars[i]]),fvars[i],"/F",outT);
    }
    ADDVAR(&(outVars["nJets"]),"nJets","/I",outPT);
    ADDVAR(&(outVars["nBjets"]),"nBjets","/I",outPT);

#ifdef HISTOGRAMS_ON
    //BOOK HISTOGRAMS
    HistTool ht;
    ht.setNsyst(0);
    ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
    ht.addHist("ch_tag",       new TH1F("ch_tag",      ";Channel Tag;Events",10,0,10));
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
	
    // normalization and event count
    ht.addHist("evt_count",    new TH1F("evt_count",   ";Selection Stage;Events",10,0,10));
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(1,"Total");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(2,"Sumweighted");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(3,"preselection");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(4,"=2 p (data)");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(5,"trigger");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(6,"=1 lep");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(7,"#geq4 jets");
    ht.getPlots()["evt_count"]->GetXaxis()->SetBinLabel(8,"#geq2 bjets");
    ht.getPlots()["evt_count"]->SetBinContent(1,counter->GetBinContent(1));
    ht.getPlots()["evt_count"]->SetBinContent(2,counter->GetBinContent(2));
    ht.getPlots()["evt_count"]->SetBinContent(3,counter->GetBinContent(3));	

    // proton count
    ht.addHist("pn_count",    new TH1F("pn_count",   ";;Events",5,0,5));
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(1,"0 hits");
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(2,"RP0 hit");
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(3,"RP1 hit");
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(4,"both RP");
	ht.getPlots()["pn_count"]->GetXaxis()->SetBinLabel(5,"both RP nB#geq1");
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
	// selector.minJetPt = 25;
    
	// JEC/JER settings
	int sys = 0;
	if(systVar.find("jerUp")!=string::npos) sys = 1;
	if(systVar.find("jerDn")!=string::npos) sys = -1;
	if(systVar.find("jecUp")!=string::npos) sys = 2;
	if(systVar.find("jecDn")!=string::npos) sys = -2;
	if(sys>0){cout << "INFO: Running JEC/JER up variation"<<endl;}
	else if(sys<0){cout << "INFO: Running JEC/JER down variation"<<endl;}
	else{cout << "INFO: Running nominal jet callibration"<<endl;}

	// btagSF settings
	string option = "central";
	if(systVar.find("btagUp")!=string::npos) {
		cout << "INFO: Running btagSF up variation"<<endl;
		option = "up";}
	else if(systVar.find("btagDn")!=string::npos) {
		cout << "INFO: Running btagSF down variation"<<endl;
		option = "down";}
	else{cout << "INFO: Running nominal btag SF"<<endl;}

		
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
        if(!isData) btvSF.updateBTagDecisions(ev, option, option); 
        //btvSF.updateBTagDecisions(ev); // check it !!!
        jec.smearJetEnergies(ev);
        
        //////////////////////////
        // RECO LEVEL SELECTION //
        //////////////////////////
        TString chTag = selector.flagFinalState(ev, {}, {}, sys); // writes the name in chTag, last argument is JEC/JER systematics
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
        std::vector<Jet>      &jets        = selector.getJets();
        std::vector<Jet>      bJets,lightJets;
        std::vector<Particle> selectedLeptons;
        double p1_xi =0.; // proton in positive pot
        double p2_xi =0.; // proton in negative pot
        
        //  selection of lightJets and bJets
        for(size_t ij=0; ij<jets.size(); ij++) {
            if(jets[ij].flavor()==5) bJets.push_back(jets[ij]);
            else                     lightJets.push_back(jets[ij]);
        }
        
        // selection of leptons
        for( size_t i_lept=0;i_lept<leptons.size();i_lept++) {
            if (leptons[i_lept].pt()<30.) continue;
			if (leptons[i_lept].id()==11 && fabs(leptons[i_lept].eta())>2.1) continue;
            //        if (leptons[i_lept].reliso()>0.10) continue;   //usually tighter
            selectedLeptons.push_back(leptons[i_lept]);
        }
		
        // selection of protons
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
    }
	
        outVars["nJets"]=jets.size();
        outVars["nBjets"]=bJets.size();
 
		// Store proton pool at preselection
		if(ev.isData){
			m_protonVars_p1_xi = p1_xi;
			m_protonVars_p2_xi = p2_xi;
			outPT->Fill();
		}

    // ---- EVENT SELECTION --------------------------------------------------------------

	if(chTag!="E" && chTag!="M" )   continue; // events with electrons (id=11) or muons (id=13)
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 4, plotwgts); // count events after channel selection
#endif
	if (selectedLeptons.size()!=1) continue; // ONLY events with 1 selected lepton
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 5, plotwgts); // count events after selection on number of leptons (SHOULD BE SAME)
#endif
        if ( jets.size()  < 4 )        continue; // ONLY events with at least 4 jets
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 6, plotwgts); // count events after selection on number of jets
#endif
        if ( bJets.size() < 2 )        continue; // ONLY events with at least 2 BJets
#ifdef HISTOGRAMS_ON
        ht.fill("evt_count", 7, plotwgts); // count events after selection on number of Bjets
#endif

#ifdef HISTOGRAMS_ON
        ht.fill("puwgtctr",0,plotwgts);
#endif  
        outVars["weight"] = outVars["gen_wgt"] = outVars["toppt_wgt"] = outVars["selSF_wgt"] = outVars["trigSF_wgt"] = 1;
        outVars["pu_wgt"] = outVars["ptag_wgt"] = outVars["ptag_wgt_err"] = 1;
		
        outVars["ren_err"] = outVars["scale_err"] = outVars["fac_err"] = 0;
		
		
        if (!ev.isData) {
            wgt  = (normH? normH->GetBinContent(1) : 1.0);          // norm weight
            double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);     // pu weight
            std::vector<double>puPlotWgts(1,puWgt);
            
#ifdef HISTOGRAMS_ON
            ht.fill("puwgtctr",1,puPlotWgts);
#endif
            
            // lepton trigger*selection weights (update the code later)
            EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons,{},{},period);
            EffCorrection_t  selSF = lepEffH.getOfflineCorrection(leptons[0], period);
			outVars["trigSF_wgt"] = trigSF.first;
			outVars["trigSF_wgt_err"] = trigSF.second;
			if(!outVars["trigSF_wgt"]) cout << "WARNING: trigSF = 0, check your selection! " << endl;
			outVars["selSF_wgt"] = selSF.first;
			outVars["selSF_wgt_err"] = selSF.second;
			if(!outVars["trigSF_wgt"]) cout << "WARNING: trigSF = 0, check your selection! " << endl;
            wgt *= outVars["trigSF_wgt"]*outVars["selSF_wgt"];
            
            //top pt weighting
            outVars["toppt_wgt"] = 1.0;
            if(isTTbar) {
                for (int igen=0; igen<ev.ngtop; igen++) {
                    if(abs(ev.gtop_id[igen])!=6) continue;
                    outVars["toppt_wgt"] *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
                }
            }
			wgt *= outVars["toppt_wgt"];
			
			// generator weights
			outVars["gen_wgt"] = (ev.g_nw>0 ? ev.g_w[0] : 1.0);
			
            wgt *= outVars["gen_wgt"];                              // generator level weights
			
            plotwgts[0]=wgt;                                        //update weight for plotter
			outVars["weight"] = wgt;
			
			// Systematic uncertainties:
			if(ev.g_nw>5){ // scale variations
				outVars["ren_err"] = (ev.g_w[2]-ev.g_w[3])/ev.g_w[1];
				outVars["fac_err"] = (ev.g_w[4]-ev.g_w[5])/ev.g_w[1];
				outVars["scale_err"] = (ev.g_w[6]-ev.g_w[7])/ev.g_w[1];
			}
			outVars["isr_err"] = outVars["fsr_err"] = 0;
			if(isPythia8){ // Py8 PS variations 
//https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/PSweightsPythia/PythiaPSweightsSettings_cfi.py
				int nW = ev.g_npsw;
				if(nW==46){ // 46 weights of ISR and FSR, take only the first ones
					outVars["isr_err"] = (ev.g_psw[8]-ev.g_psw[6])/ev.g_psw[0];
					outVars["fsr_err"] = (ev.g_psw[9]-ev.g_psw[7])/ev.g_psw[0];
				}
				if(nW==24){ // in signal no ISR weights
					outVars["fsr_err"] = (ev.g_psw[5]-ev.g_psw[6])/ev.g_psw[0];
				}
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
        if(bJets.size()>=2 && lightJets.size()>=2)        {    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //  in theory this selection was made already
            
            //determine the neutrino kinematics
            TLorentzVector met(0,0,0,0);
            met.SetPtEtaPhiM(ev.met_pt,0.,ev.met_phi,0.);
            neutrinoPzComputer.SetMET(met);
            neutrinoPzComputer.SetLepton(leptons[0].p4());
            float nupz=neutrinoPzComputer.Calculate();
            TLorentzVector neutrino(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
            
            // compute scalar ht
            float scalarht(0.);
            for(size_t ij=0; ij<jets.size(); ij++) {
                scalarht += jets[ij].pt();
            }
            scalarht += leptons[0].pt();
            scalarht += neutrino.Pt();

            TLorentzVector t_rec_had(0,0,0,0);
            TLorentzVector t_rec_lep(0,0,0,0);
            TLorentzVector bJet_had(0,0,0,0);
            TLorentzVector bJet_lep(0,0,0,0);
            TLorentzVector lightJet0(0,0,0,0);
            TLorentzVector lightJet1(0,0,0,0);

            TLorentzVector t_rec(0,0,0,0);
            TLorentzVector tbar_rec(0,0,0,0);

            // the combination of object whose mass is closer to the top one is chosen
            vector <int> combination;          // indexes of:   bJetHad, lightJet1, lightJet2, bJetLep
            vector <vector <int> > marker;  // list of all combinations
            vector <double> mistake;        // mass differences between top candidate and true mass for each combination

            // create list of all combinations
            // for the time being only the first 2 LightJets are used
            if (bJets.size()>=2) {
                marker.push_back( {0,0,1,1} );
                marker.push_back( {1,0,1,0} );
            }
            if (bJets.size()>=3) {
                marker.push_back( {0,0,1,2} );
                marker.push_back( {2,0,1,0} );
                marker.push_back( {1,0,1,2} );
                marker.push_back( {2,0,1,1} );
            }
            if (bJets.size()==4) {
                marker.push_back( {0,0,1,3} );
                marker.push_back( {3,0,1,0} );
                marker.push_back( {1,0,1,3} );
                marker.push_back( {3,0,1,1} );
                marker.push_back( {2,0,1,3} );
                marker.push_back( {3,0,1,2} );
            }

            // build t_rec_had and t_rec_lep for every combination,
            // compute difference between top candidate and true mass for every combination
            // and store it in the mistake vector
            for (size_t i_comb=0; i_comb<marker.size(); i_comb++) {
                t_rec_had = bJets[ marker[i_comb][0] ].p4() +lightJets[ marker[i_comb][1] ].p4()+ lightJets[  marker[i_comb][2] ].p4();   // b+q+q
                t_rec_lep = bJets[ marker[i_comb][3] ].p4() +leptons[0].p4() +neutrino;            // b+l+nu
                mistake.push_back( abs(t_rec_had.M()-m_TOP)+abs(t_rec_lep.M()-m_TOP) );
            }

            int correct_index=dump_index(mistake);   // index of best candidate
            combination = marker[correct_index];     // combination of best candidate
			
            // build the best top candidates
            bJet_had  = bJets[ combination[0] ].p4();
            bJet_lep  = bJets[ combination[3] ].p4();
            lightJet0 = lightJets[ combination[1] ].p4();
            lightJet1 = lightJets[  combination[2] ].p4();
			
            t_rec_had = bJet_had + lightJet0 + lightJet1;   // b+q+q
            t_rec_lep = bJet_lep + leptons[0].p4() + neutrino;            // b+l+nu

            bool isThadronic = 0;
            if (leptons[0].charge()>0.) {
                isThadronic = 0; // the top quark decayed in leptons, the tbar decayed in jets
                t_rec    = t_rec_lep;
                tbar_rec = t_rec_had;
            } else {
                isThadronic = 1; // the top quark decayed in jets, the tbar decayed in leptons
                t_rec    = t_rec_had;
                tbar_rec = t_rec_lep;
            }
            
            //reconstructed ttbar system
            TLorentzVector ttbarSystem_rec = t_rec_had + t_rec_lep;

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


#ifdef SAVEERRORS_ON
            // MATTEO'S CODE ----------- compute UNCERTAINTIES ---------------------------------------------
            TF1 theta("theta", "2*atan(exp(-x))",0, 1e10);
            TF1 err_theta_fracdeltaeta("err_theta","abs(-2*(exp(-x)/(1+pow(exp(-x),2))))",-1e10,1e10);

            // translate my quantities in Matteo's quantities
            double pt_b1     =   bJet_had.Pt();
            double eta_b1    =   bJet_had.Eta();    // .Rapidity();
            double phi_b1    =   bJet_had.Phi();

            double pt_b2     =   bJet_lep.Pt();
            double eta_b2    =   bJet_lep.Eta();    // .Rapidity();
            double phi_b2    =   bJet_lep.Phi();

            double pt_q1     =   lightJet0.Pt();
            double eta_q1    =   lightJet0.Eta();    // .Rapidity();
            double phi_q1    =   lightJet0.Phi();

            double pt_q2     =   lightJet1.Pt();
            double eta_q2    =   lightJet1.Eta();    // .Rapidity();
            double phi_q2    =   lightJet1.Phi();

            double pt_l     =   leptons[0].p4().Pt();
            double eta_l    =   leptons[0].p4().Eta();    // .Rapidity();
            double phi_l    =   leptons[0].p4().Phi();
            double m_l      =   leptons[0].p4().M();

//            double pt_nu     =   neutrino.Pt();
//            double eta_nu    =   neutrino.Eta();    // .Rapidity();
//            double phi_nu    =   neutrino.Phi();
//            double m_nu      =   neutrino.M();

            // compute uncertainties on eta,p,phi,theta quantities
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

            // propagate uncertainties to x,y,z quantities
             double e_b1_x=sqrt(pow(e_b1_p*sin(theta.Eval(eta_b1))*cos(phi_b1),2)+pow(pt_b1/sin(theta.Eval(eta_b1))*cos(theta.Eval(eta_b1))*cos(phi_b1)*e_b1_theta,2)+pow(pt_b1/sin(theta.Eval(eta_b1))*sin(theta.Eval(eta_b1))*sin(phi_b1)*e_b1_phi,2));
//             if ( e_b1_x>4.*bJet_Hadronic.Px()) e_b1_x=4.*bJet_Hadronic.Px(); // if the error diverges set it to a maximum
               double e_b1_y=sqrt(pow(e_b1_p*sin(theta.Eval(eta_b1))*sin(phi_b1),2)+pow(pt_b1/sin(theta.Eval(eta_b1))*cos(theta.Eval(eta_b1))*sin(phi_b1)*e_b1_theta,2)+pow(pt_b1/sin(theta.Eval(eta_b1))*sin(theta.Eval(eta_b1))*cos(phi_b1)*e_b1_phi,2));
//             if ( e_b1_y>4.*bJet_Hadronic.Py()) e_b1_y=4.*bJet_Hadronic.Py(); // if the error diverges set it to a maximum
               double e_b1_z=sqrt(pow(e_b1_p*cos(theta.Eval(eta_b1)),2)+pow(pt_b1/sin(theta.Eval(eta_b1))*sin(theta.Eval(eta_b1))*e_b1_theta,2));
//             if ( e_b1_z>4.*bJet_Hadronic.Pz()) e_b1_z=4.*bJet_Hadronic.Pz(); // if the error diverges set it to a maximum

             double e_b2_x=sqrt(pow(e_b2_p*sin(theta.Eval(eta_b2))*cos(phi_b2),2)+pow(pt_b2/sin(theta.Eval(eta_b2))*cos(theta.Eval(eta_b2))*cos(phi_b2)*e_b2_theta,2)+pow(pt_b2/sin(theta.Eval(eta_b2))*sin(theta.Eval(eta_b2))*sin(phi_b2)*e_b2_phi,2));
//             if ( e_b2_x>4.*bJet_Leptonic.Px()) e_b2_x=4.*bJet_Leptonic.Px(); // if the error diverges set it to a maximum
             double e_b2_y=sqrt(pow(e_b2_p*sin(theta.Eval(eta_b2))*sin(phi_b2),2)+pow(pt_b2/sin(theta.Eval(eta_b2))*cos(theta.Eval(eta_b2))*sin(phi_b2)*e_b2_theta,2)+pow(pt_b2/sin(theta.Eval(eta_b2))*sin(theta.Eval(eta_b2))*cos(phi_b2)*e_b2_phi,2));
//             if ( e_b2_y>4.*bJet_Leptonic.Py()) e_b2_y=4.*bJet_Leptonic.Py(); // if the error diverges set it to a maximum
             double e_b2_z=sqrt(pow(e_b2_p*cos(theta.Eval(eta_b2)),2)+pow(pt_b2/sin(theta.Eval(phi_b2))*sin(theta.Eval(eta_b2))*e_b2_theta,2));
//             if ( e_b2_z>4.*bJet_Leptonic.Pz()) e_b2_z=4.*bJet_Leptonic.Pz(); // if the error diverges set it to a maximum

             double e_l_x=sqrt(pow(e_l_p*sin(theta.Eval(eta_l))*cos(phi_l),2)+pow(pt_l/sin(theta.Eval(eta_l))*cos(theta.Eval(eta_l))*cos(phi_l)*e_l_theta,2)+pow(pt_l/sin(theta.Eval(eta_l))*sin(theta.Eval(eta_l))*sin(phi_l)*e_l_phi,2));
//             if ( e_l_x>4.*leptons[0].p4().Px() ) e_l_x=4.*leptons[0].p4().Px(); // if the error diverges set it to a maximum
             double e_l_y= sqrt(pow(e_l_p*sin(theta.Eval(eta_l))*sin(phi_l),2)+pow(pt_b1/sin(theta.Eval(eta_l))*cos(theta.Eval(eta_l))*sin(phi_l)*e_l_theta,2)+pow(pt_l/sin(theta.Eval(eta_l))*sin(theta.Eval(eta_l))*cos(phi_l)*e_l_phi,2));
//             if ( e_l_y>4.*leptons[0].p4().Py() ) e_l_y=4.*leptons[0].p4().Py(); // if the error diverges set it to a maximum
             double e_l_z=sqrt(pow(e_l_p*cos(theta.Eval(eta_l)),2)+pow(pt_l/sin(theta.Eval(eta_l))*sin(theta.Eval(eta_l))*e_l_theta,2));
//             if ( e_l_z>4.*leptons[0].p4().Pz() ) e_l_z=4.*leptons[0].p4().Pz(); // if the error diverges set it to a maximum

             double e_nu_x=6.*neutrino.Pt(); // 6 times because six is big enough. We set a very big error on neutrino
             double e_nu_y=6.*neutrino.Pt();
             double e_nu_z=6.*neutrino.Pt();

             double e_q1_x=sqrt(pow(e_q1_p*sin(theta.Eval(eta_q1))*cos(phi_q1),2)+pow(pt_q1/sin(theta.Eval(eta_q1))*cos(theta.Eval(eta_q1))*cos(phi_q1)*e_q1_theta,2)+pow(pt_q1/sin(theta.Eval(eta_q1))*sin(theta.Eval(eta_q1))*sin(phi_b1)*e_q1_phi,2));
             double e_q1_y=sqrt(pow(e_q1_p*sin(theta.Eval(eta_q1))*sin(phi_q1),2)+pow(pt_q1/sin(theta.Eval(eta_q1))*cos(theta.Eval(eta_q1))*sin(phi_q1)*e_q1_theta,2)+pow(pt_q1/sin(theta.Eval(eta_q1))*sin(theta.Eval(eta_q1))*cos(phi_q1)*e_q1_phi,2));
             double e_q1_z=sqrt(pow(e_q1_p*cos(theta.Eval(eta_q1)),2)+pow(pt_q1/sin(theta.Eval(eta_q1))*sin(theta.Eval(eta_q1))*e_q1_theta,2));

             double e_q2_x=sqrt(pow(e_q2_p*sin(theta.Eval(eta_q2))*cos(phi_q2),2)+pow(pt_q1/sin(theta.Eval(eta_q2))*cos(theta.Eval(eta_q2))*cos(phi_q2)*e_q2_theta,2)+pow(pt_q2/sin(theta.Eval(eta_q2))*sin(theta.Eval(eta_q2))*sin(phi_b2)*e_q2_phi,2));
             double e_q2_y=sqrt(pow(e_q2_p*sin(theta.Eval(eta_q2))*sin(phi_q2),2)+pow(pt_q2/sin(theta.Eval(eta_q2))*cos(theta.Eval(eta_q2))*sin(phi_q2)*e_q2_theta,2)+pow(pt_q2/sin(theta.Eval(eta_q2))*sin(theta.Eval(eta_q2))*cos(phi_q2)*e_q2_phi,2));
             double e_q2_z=sqrt(pow(e_q2_p*cos(theta.Eval(eta_q2)),2)+pow(pt_q2/sin(theta.Eval(eta_q2))*sin(theta.Eval(eta_q2))*e_q2_theta,2));
            // END OF MATTEO'S CODE -------------------------------------------------------------------------------------------

            // translate Matteo's quantities back (maybe fill some TLorentzVector??)

            TVector3 errXYZ_bJet_had(e_b1_x,e_b1_y,e_b1_z);
            TVector3 errXYZ_bJet_lep(e_b2_x,e_b2_y,e_b2_z);
            TVector3 errXYZ_lightJet0(e_q1_x,e_q1_y,e_q1_z);
            TVector3 errXYZ_lightJet1(e_q2_x,e_q2_y,e_q2_z);
            TVector3 errXYZ_lepton(e_l_x,e_l_y,e_l_z);
            TVector3 errXYZ_neutrino(e_nu_x,e_nu_y,e_nu_z);
/*
            TVector3 errPOLAR_bJet_had(0,0,0);      errPOLAR_bJet_had.SetPtEtaPhi();
            TVector3 errPOLAR_bJet_lep(0,0,0);
            TVector3 errPOLAR_lightJet0(0,0,0);
            TVector3 errPOLAR_lightJet1(0,0,0);
            TVector3 errPOLAR_lepton(0,0,0);
            TVector3 errPOLAR_neutrino(0,0,0);
*/
#endif

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
            outVars["t_pt"]     = t_rec.Pt();
            outVars["t_eta"]    = t_rec.Rapidity();
            outVars["t_phi"]    = t_rec.Phi();
            outVars["t_m"]      = t_rec.M();
            outVars["t_charge"] =1;
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

            /* quantities of objects used to build the ttbar */
            outVars["l_px"]=leptons[0].p4().Px();
            outVars["l_py"]=leptons[0].p4().Py();
            outVars["l_pz"]=leptons[0].p4().Pz();
            outVars["nu_px"]=neutrino.Px();
            outVars["nu_py"]=neutrino.Py();
            outVars["nu_pz"]=neutrino.Pz();
            outVars["nu_complex"] = neutrinoPzComputer.IsComplex() ? 1 : 0;
            outVars["bJet_had_px"]=bJet_had.Px();
            outVars["bJet_had_py"]=bJet_had.Py();
            outVars["bJet_had_pz"]=bJet_had.Pz();
            outVars["bJet_lep_px"]=bJet_lep.Px();
            outVars["bJet_lep_py"]=bJet_lep.Py();
            outVars["bJet_lep_pz"]=bJet_lep.Pz();
            outVars["lightJet0_px"]=lightJet0.Px();
            outVars["lightJet0_py"]=lightJet0.Py();
            outVars["lightJet0_pz"]=lightJet0.Pz();
            outVars["lightJet1_px"]=lightJet1.Px();
            outVars["lightJet1_py"]=lightJet1.Py();
            outVars["lightJet1_pz"]=lightJet1.Pz();

#ifdef SAVEERRORS_ON
            outVars["e_l_px"]       = errXYZ_lepton.Px();
            outVars["e_l_py"]       = errXYZ_lepton.Py();
            outVars["e_l_pz"]       = errXYZ_lepton.Pz();
            outVars["e_nu_px"]      = errXYZ_neutrino.Px();
            outVars["e_nu_py"]      = errXYZ_neutrino.Py();
            outVars["e_nu_pz"]      = errXYZ_neutrino.Pz();
            outVars["e_nu_pxpy"]    = 0.;
            outVars["e_bJet_had_px"]= errXYZ_bJet_had.Px();
            outVars["e_bJet_had_py"]= errXYZ_bJet_had.Py();
            outVars["e_bJet_had_pz"]= errXYZ_bJet_had.Pz();
            outVars["e_bJet_lep_px"]= errXYZ_bJet_lep.Px();
            outVars["e_bJet_lep_py"]= errXYZ_bJet_lep.Py();
            outVars["e_bJet_lep_pz"]= errXYZ_bJet_lep.Pz();
            outVars["e_lightJet0_px"]=errXYZ_lightJet0.Px();
            outVars["e_lightJet0_py"]=errXYZ_lightJet0.Py();
            outVars["e_lightJet0_pz"]=errXYZ_lightJet0.Pz();
            outVars["e_lightJet1_px"]=errXYZ_lightJet1.Px();
            outVars["e_lightJet1_py"]=errXYZ_lightJet1.Py();
            outVars["e_lightJet1_pz"]=errXYZ_lightJet1.Pz();
#endif
            
#ifdef DEBUG_ON
            // test variables
            outVars["test1"]       = (e_l_p-errXYZ_lepton.Mag())/e_l_p;
            outVars["test2"]       = (e_l_eta-errXYZ_lepton.Eta())/e_l_eta;
            outVars["test3"]       = (e_l_phi-errXYZ_lepton.Phi())/e_l_phi;
            outVars["test4"]       = (e_l_theta-errXYZ_lepton.Theta())/e_l_theta;
            outVars["test5"]       = -1.;
#endif
            
            // quantities for all objects in event
            outVars["nLightJets"]=lightJets.size();
            outVars["ht"]=scalarht;
            outVars["cat"]=float(ch_tag);

            outVars["l_pt"]=leptons[0].Pt();
            outVars["l_eta"]=leptons[0].Rapidity();
            outVars["l_phi"]=leptons[0].Phi();
            outVars["l_m"]=leptons[0].M();
            outVars["l_E"]=leptons[0].E();
            outVars["lepton_isolation"]=ev.l_relIso[0];

            outVars["nu_pt"]=neutrino.Pt();
            outVars["nu_eta"]=neutrino.Rapidity();
            outVars["nu_phi"]=neutrino.Phi();

            outVars["p1_xi"]=p1_xi;
            outVars["p2_xi"]=p2_xi;
			
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
            // quantities of generated objects (MC truth)
            outVars["gen_ttbar_pt"]=ttbarSystem_gen.Pt();
            outVars["gen_ttbar_eta"]=ttbarSystem_gen.Rapidity();
            outVars["gen_ttbar_phi"]=ttbarSystem_gen.Phi();
            outVars["gen_ttbar_m"]=ttbarSystem_gen.M();
            outVars["gen_ttbar_E"]=ttbarSystem_gen.E();

            outVars["gen_t_pt"]=t_gen.Pt();
            outVars["gen_t_eta"]=t_gen.Rapidity();
            outVars["gen_t_phi"]=t_gen.Phi();
            outVars["gen_t_m"]=t_gen.M();
            outVars["gen_tbar_pt"]=tbar_gen.Pt();
            outVars["gen_tbar_eta"]=tbar_gen.Rapidity();
            outVars["gen_tbar_phi"]=tbar_gen.Phi();
            outVars["gen_tbar_m"]=tbar_gen.M();

            outVars["gen_b_pt"]=b_gen.Pt();
            outVars["gen_b_eta"]=b_gen.Rapidity();
            outVars["gen_b_phi"]=b_gen.Phi();
            outVars["gen_b_m"]=b_gen.M();
            outVars["gen_bbar_pt"]=bbar_gen.Pt();
            outVars["gen_bbar_eta"]=bbar_gen.Rapidity();
            outVars["gen_bbar_phi"]=bbar_gen.Phi();
            outVars["gen_bbar_m"]=bbar_gen.M();

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

