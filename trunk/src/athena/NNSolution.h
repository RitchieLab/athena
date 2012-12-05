/*
Copyright Marylyn Ritchie 2011

This file is part of ATHENA.

ATHENA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ATHENA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ATHENA.  If not, see <http://www.gnu.org/licenses/>.
*/
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
      bool map_used,  bool ott_dummy, bool continmap_used);
    
    /// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
    virtual void adjust_score_out(Dataset* train_set, Dataset* test_set);
    
    /// Adjusts output of the scores when needed for training set only
    virtual void adjust_score_out(Dataset* train_set);
    
    /// Adjusts score passed and returns value
    virtual float adjust_score_out(float score, int nIndsTested, float sstotal);
    
    void set_gram_depth(int dep){gram_depth = dep;}
    int get_gram_depth(){return gram_depth;}
    
    void set_nn_depth(int dep){nn_depth = dep;}
    int get_nn_depth(){return nn_depth;}
    
private:

    /// adjusts indexes back to original values if dummy-encoded
    int adjust_dummy_encoding(int genotype);
    
    float alter_score(float mse, int total_inds, float sstotal);
    
    int calc_inds(Dataset* set, float& sstotal);
    
    int nn_depth, gram_depth;
    
};


#endif	/* _NNSOLUTION_H */

