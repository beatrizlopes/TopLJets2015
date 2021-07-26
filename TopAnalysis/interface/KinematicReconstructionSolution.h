#ifndef _KinematicReconstructionSolution_h
#define _KinematicReconstructionSolution_h

#include <map>
#include <vector>
#include <TLorentzVector.h>


class KinematicReconstructionSolution{
    
public:
    
    /// Enumeration for all defined weights of a kinematic reconstruction solution
    enum WeightType{defaultForMethod, neutrinoEnergy, averagedSumSmearings_mlb, undefinedWeight};
    
    
    
    /// Default constructor
    KinematicReconstructionSolution();
    
    /// Constructor for setting values of the solution
    KinematicReconstructionSolution(const std::vector<TLorentzVector>* const allLeptons,
                                    const std::vector<TLorentzVector>* const allJets,
                                    const int leptonIndex, const int antiLeptonIndex,
                                    const int bjetIndex, const int antiBjetIndex,
                                    const TLorentzVector& top, const TLorentzVector& antiTop,
                                    const TLorentzVector& wMinus, const TLorentzVector& wPlus,
                                    const TLorentzVector& neutrino, const TLorentzVector& antiNeutrino,
                                    const double& reconstructedTopMass,
                                    const int numberOfBtags,
                                    const std::map<WeightType, double>& m_weight,
                                    const bool isNoSmearSol_ = false);
    
    /// Copy constructor
    KinematicReconstructionSolution(const KinematicReconstructionSolution& kinematicReconstructionSolution);
    
    /// Destructor
    ~KinematicReconstructionSolution(){}
    
    
    
    /// Assignment operator (needed to allow std::vector<KinematicReconstructionSolution> with const data members)
    const KinematicReconstructionSolution& operator=(const KinematicReconstructionSolution& rhs){return rhs;}
    
    /// Print all stored quantities
    void print()const;
    
    /// Check if this is really filled with a solution, or whether it is a dummy
    bool dummy()const{return !allLeptons_;}
    
    
    
    const TLorentzVector& lepton()const{return allLeptons_->at(leptonIndex_);}
    const TLorentzVector& antiLepton()const{return allLeptons_->at(antiLeptonIndex_);}
    const TLorentzVector& bjet()const{return allJets_->at(bjetIndex_);}
    const TLorentzVector& antiBjet()const{return allJets_->at(antiBjetIndex_);}
    
    const TLorentzVector& top()const{return top_;}
    const TLorentzVector& antiTop()const{return antiTop_;}
    const TLorentzVector& wMinus()const{return wbar_;}
    const TLorentzVector& wPlus()const{return w_;}
    const TLorentzVector& neutrino()const{return neutrino_;}
    const TLorentzVector& antiNeutrino()const{return antiNeutrino_;}
    
    TLorentzVector ttbar()const{return top_ + antiTop_;}
    
    // FIXME: do these two make sense? these are just the sum of reco objects, not kinReco Ws
    //TLorentzVector wMinus()const{return this->lepton() + antiNeutrino_;}
    //TLorentzVector wPlus()const{return this->antiLepton() + neutrino_;}
    
    int leptonIndex()const{return leptonIndex_;}
    int antiLeptonIndex()const{return antiLeptonIndex_;}
    int bjetIndex()const{return bjetIndex_;}
    int antiBjetIndex()const{return antiBjetIndex_;}
    const double& reconstructedTopMass()const{return reconstructedTopMass_;}
    int numberOfBtags()const{return numberOfBtags_;}
    const double& weight(const WeightType weightType =defaultForMethod)const{return m_weight_.at(weightType);}
    
    const std::map<WeightType, double>& weightMap()const{return m_weight_;}
    
    bool isNoSmearSol()const{return isNoSmearSol_;}
    
    
    
private:
    
    const std::vector<TLorentzVector>* const allLeptons_;
    const std::vector<TLorentzVector>* const allJets_;
    
    // Indices of input quantities
    const int leptonIndex_;
    const int antiLeptonIndex_;
    const int bjetIndex_;
    const int antiBjetIndex_;
    
    // Reconstructed quantities
    const TLorentzVector top_;
    const TLorentzVector antiTop_;
    //MOD
    const TLorentzVector w_;
    const TLorentzVector wbar_;
    const TLorentzVector neutrino_;
    const TLorentzVector antiNeutrino_;
    
    const double reconstructedTopMass_;
    
    // Number of b-tagged jets for solution
    const int numberOfBtags_;
    
    // Weights
    const std::map<WeightType, double> m_weight_;
    
    // Is there solution without smearing
    const bool isNoSmearSol_;
};





class KinematicReconstructionSolutions{
    
public:
    
    /// Empty constructor
    KinematicReconstructionSolutions();
    
    /// Destructor
    ~KinematicReconstructionSolutions(){}
    
    
    
    /// Add a vector of solutions
    void addSolutions(const std::vector<KinematicReconstructionSolution>& solutions);
    
    ///  Add a solution
    void addSolution(const KinematicReconstructionSolution& solution);
    
    
    /// Number of all solutions
    size_t numberOfSolutions()const{return v_solution_.size();}
    
    /// Number of 2 b-tag solutions
    size_t numberOfSolutionsTwoBtags()const{return v_solutionTwoBtags_.size();}
    
    /// Number of 1 b-tag solutions
    size_t numberOfSolutionsOneBtag()const{return v_solutionOneBtag_.size();}
    
    /// Number of 0 b-tag solutions
    size_t numberOfSolutionsNoBtags()const{return v_solutionNoBtags_.size();}
    
    
    
    /// Access from all solutions the one selected with solutionNumber, ranked by decreasing specific weight
    const KinematicReconstructionSolution& solution(const KinematicReconstructionSolution::WeightType weightType =KinematicReconstructionSolution::defaultForMethod,
                                                    const size_t solutionNumber =0)const;
    
    /// Access from solutions with 2 b-tags the one selected with solutionNumber, ranked by decreasing specific weight
    const KinematicReconstructionSolution& solutionTwoBtags(const KinematicReconstructionSolution::WeightType weightType =KinematicReconstructionSolution::defaultForMethod,
                                                            const size_t solutionNumber =0)const;
    
    /// Access from solutions with 1 b-tag the one selected with solutionNumber, ranked by decreasing specific weight
    const KinematicReconstructionSolution& solutionOneBtag(const KinematicReconstructionSolution::WeightType weightType =KinematicReconstructionSolution::defaultForMethod,
                                                           const size_t solutionNumber =0)const;
    
    /// Access from solutions with 0 b-tags the one selected with solutionNumber, ranked by decreasing specific weight
    const KinematicReconstructionSolution& solutionNoBtags(const KinematicReconstructionSolution::WeightType weightType =KinematicReconstructionSolution::defaultForMethod,
                                                           const size_t solutionNumber =0)const;
    
    
    
private:
    
    /// Insert solution index for specific weight in vector for all solutions, ranked by weight
    void insertIndex(const size_t solutionIndex,
                     const double weight, std::vector<size_t>& v_index)const;
    
    /// Insert solution index for specific weight in vector for b-tag categorised solutions, ranked by weight
    void insertIndexByCategory(const std::vector<size_t>& v_solutionIndex,
                               const double weight, std::vector<size_t>& v_solutionIndexByCategory)const;
    
    
    
    /// Vector containing all solutions
    std::vector<KinematicReconstructionSolution> v_solution_;
    
    /// Vector containing indices of the solutions with 2 b-tags stored in v_solution_
    std::vector<size_t> v_solutionTwoBtags_;
    
    /// Vector containing indices of the solutions with 1 b-tag stored in v_solution_
    std::vector<size_t> v_solutionOneBtag_;
    
    /// Vector containing indices of the solutions with 0 b-tags stored in v_solution_
    std::vector<size_t> v_solutionNoBtags_;
    
    
    
    /// Map associating specific weight type to vector containing indices of all solutions, ordered for this weight
    std::map<KinematicReconstructionSolution::WeightType, std::vector<size_t>> m_weightIndex_;
    
    /// Map associating specific weight type to vector containing indices of 2 b-tag solutions, ordered for this weight
    std::map<KinematicReconstructionSolution::WeightType, std::vector<size_t>> m_weightIndexTwoBtags_;
    
    /// Map associating specific weight type to vector containing indices of 1 b-tag solutions, ordered for this weight
    std::map<KinematicReconstructionSolution::WeightType, std::vector<size_t>> m_weightIndexOneBtag_;
    
    /// Map associating specific weight type to vector containing indices of 0 b-tag solutions, ordered for this weight
    std::map<KinematicReconstructionSolution::WeightType, std::vector<size_t>> m_weightIndexNoBtags_;
};





#endif


