#ifndef _KinematicReconstruction_MeanSol_h
#define _KinematicReconstruction_MeanSol_h

#include <vector>

#include <TLorentzVector.h>


class KinematicReconstruction_MeanSol{
    
public:
    
    KinematicReconstruction_MeanSol(const double& topm, const double& Wmmass, const double& Wpmass);
    ~KinematicReconstruction_MeanSol();  
    void add(const TLorentzVector& top, const TLorentzVector& topbar, const TLorentzVector& n, const TLorentzVector& nbar, const double& weight, const double& mbl_weight);
    void add(const TLorentzVector& top, const TLorentzVector& topbar, const TLorentzVector& wm, const TLorentzVector& wp, const TLorentzVector& n, const TLorentzVector& nbar, const double& weight);
    
    void getMeanVect(TLorentzVector& lv, const std::vector<TLorentzVector>& vlv, const double& mass)const;
    void getMeanSol(TLorentzVector& top, TLorentzVector& topbar, TLorentzVector& wm, TLorentzVector& wp, TLorentzVector& n, TLorentzVector& nbar)const;
    double getSumWeight()const;
    int getNsol()const;
    void clear();
    
    
    // FIXME: right now this is only using getMeanSol(), since the whole class is based on TLorentzVector
    /// Get Lorentz vectors of (anti-)top and (anti-)neutrino for mean solution
    void meanSolution(TLorentzVector& top, TLorentzVector& antiTop, TLorentzVector& wm, TLorentzVector& wp, TLorentzVector& neutrino, TLorentzVector& antiNeutrino)const;
    
    
    
private:
    
    std::vector<TLorentzVector> v_top_;
    std::vector<TLorentzVector> v_topbar_;
    std::vector<TLorentzVector> v_wm_;
    std::vector<TLorentzVector> v_wp_;
    std::vector<TLorentzVector> v_n_;
    std::vector<TLorentzVector> v_nbar_;
    
    std::vector<double> v_weight_;
    double sum_weight_;  
    double max_sum_weight_;
    
    const double mass_top_;
    const double mass_Wm_;
    const double mass_Wp_;
};



#endif


