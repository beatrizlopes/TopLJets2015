#include <vector>

//#include <TLorentzVector.h>

#include "TopLJets2015/TopAnalysis/interface/KinematicReconstruction_MeanSol.h"
//#include "gammatt/include/classes.h"
//#include "gammatt/include/analysisUtils.h"



KinematicReconstruction_MeanSol::KinematicReconstruction_MeanSol(const double& topm, const double& Wmmass, const double& Wpmass):
sum_weight_(0),
max_sum_weight_(0),
mass_top_(topm),
mass_Wm_(Wmmass),
mass_Wp_(Wpmass)
{}


KinematicReconstruction_MeanSol::~KinematicReconstruction_MeanSol()
{}



void KinematicReconstruction_MeanSol::clear()
{
    v_top_.clear();
    v_topbar_.clear();
    //MOD
    v_wm_.clear();
    v_wp_.clear();
    v_n_.clear();
    v_nbar_.clear();
    v_weight_.clear();
    sum_weight_=0;
    max_sum_weight_=0;
}



void KinematicReconstruction_MeanSol::add(const TLorentzVector& top, const TLorentzVector& topbar, const TLorentzVector& n, const TLorentzVector& nbar,const double& weight,const double& mbl_weight)
{
    v_top_.push_back(top);
    v_topbar_.push_back(topbar);
    v_n_.push_back(n);
    v_nbar_.push_back(nbar);
    v_weight_.push_back(weight);
    
    sum_weight_ = sum_weight_ + weight;
    max_sum_weight_ = max_sum_weight_ + mbl_weight;
}



void KinematicReconstruction_MeanSol::add(const TLorentzVector& top, const TLorentzVector& topbar, const TLorentzVector& wm, const TLorentzVector& wp, const TLorentzVector& n,const TLorentzVector& nbar,const double& weight)
{
    v_top_.push_back(top);
    v_topbar_.push_back(topbar);
    v_wm_.push_back(wm);
    v_wp_.push_back(wp);
    v_n_.push_back(n);
    v_nbar_.push_back(nbar);
    v_weight_.push_back(weight);
    
    sum_weight_ = sum_weight_ + weight;
    max_sum_weight_ = max_sum_weight_ +weight;
}



void KinematicReconstruction_MeanSol::getMeanVect(TLorentzVector& lv, const std::vector<TLorentzVector>& vlv, const double& mass)const
{
    double px_sum=0;
    double py_sum=0;
    double pz_sum=0;
    double px=0;
    double py=0;
    double pz=0;
    
    for(int i=0;i<((int)vlv.size());++i)
    {
          px_sum=px_sum+v_weight_.at(i)*vlv.at(i).Px();
          py_sum=py_sum+v_weight_.at(i)*vlv.at(i).Py();
          pz_sum=pz_sum+v_weight_.at(i)*vlv.at(i).Pz();
    }
    
    px=px_sum/sum_weight_;
    py=py_sum/sum_weight_;
    pz=pz_sum/sum_weight_;
    
    lv.SetXYZM(px,py,pz,mass);
}



void KinematicReconstruction_MeanSol::getMeanSol(TLorentzVector& top, TLorentzVector& topbar, TLorentzVector& wm, TLorentzVector& wp, TLorentzVector& n, TLorentzVector& nbar)const
{
    this->getMeanVect(top,v_top_,mass_top_);
    this->getMeanVect(topbar,v_topbar_,mass_top_);
    this->getMeanVect(wm,v_wm_,mass_Wm_);
    this->getMeanVect(wp,v_wp_,mass_Wp_);
    this->getMeanVect(n,v_n_,0);
    this->getMeanVect(nbar,v_nbar_,0);
}



double KinematicReconstruction_MeanSol::getSumWeight()const
{
    return sum_weight_; // for 1 weight
//     return max_sum_weight_; // for 2 weights
    
}



int KinematicReconstruction_MeanSol::getNsol()const
{
    return (int)v_top_.size();
}



void KinematicReconstruction_MeanSol::meanSolution(TLorentzVector& top, TLorentzVector& antiTop, TLorentzVector& wm, TLorentzVector& wp, TLorentzVector& neutrino, TLorentzVector& antiNeutrino)const
{
    TLorentzVector topTemp = top;
    TLorentzVector antiTopTemp = antiTop;
    TLorentzVector wmTemp = wm;
    TLorentzVector wpTemp = wp;
    TLorentzVector neutrinoTemp = neutrino;
    TLorentzVector antiNeutrinoTemp = antiNeutrino;
    
    this->getMeanVect(topTemp, v_top_, mass_top_);
    this->getMeanVect(antiTopTemp, v_topbar_, mass_top_);
    this->getMeanVect(wmTemp,v_wm_,mass_Wm_);
    this->getMeanVect(wpTemp,v_wp_,mass_Wp_);
    this->getMeanVect(neutrinoTemp, v_n_, 0);
    this->getMeanVect(antiNeutrinoTemp, v_nbar_, 0);
    
    top = topTemp;
    antiTop = antiTopTemp;
    wm = wmTemp;
    wp = wpTemp;
    neutrino = neutrinoTemp;
    antiNeutrino = antiNeutrinoTemp;
}




