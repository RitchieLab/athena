/* 
 * File:   NNSolution.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 4:38 PM
 */

#ifndef _NNSOLUTION_H
#define	_NNSOLUTION_H

#include "Solution.h"
#include "TerminalSymbCreator.h"
#include "SolutionCalculator.h"

class NNSolution: public Solution{
    
public:
   
    /// Constructor
    NNSolution(){set_name("Neural Network");}
 
    virtual Solution* clone();
    
    /// returns the genotypes present in the solution
    virtual vector<int> get_genotypes(bool dummy_encoded=true);
    
    /// returns covariates present in the solution
    virtual vector<int> get_covariates();
    
    /// outputs a more human-readable version of the network
    virtual void output_clean(std::ostream& os, data_manage::Dataholder& data,
      bool map_used,  bool ott_dummy);
    
    /// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
    virtual void adjust_score_out(Dataset* train_set, Dataset* test_set);
    
    /// Adjusts output of the scores when needed for training set only
    virtual void adjust_score_out(Dataset* train_set);
    
    /// Adjusts score passed and returns value
    virtual float adjust_score_out(float score, int nIndsTested, float sstotal);
    
private:

    /// adjusts indexes back to original values if dummy-encoded
    int adjust_dummy_encoding(int genotype);
    
    float alter_score(float mse, int total_inds, float sstotal);
    
    int calc_inds(Dataset* set, float& sstotal);
    
};


#endif	/* _NNSOLUTION_H */

