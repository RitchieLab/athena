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
 * File:   NNSolutionCreator.h
 * Author: dudeksm
 *
 * Created on December 1, 2008, 4:18 PM
 */

#ifndef _NNSOLUTIONCREATOR_H
#define	_NNSOLUTIONCREATOR_H

#include "SolutionCreator.h"
#include "TerminalSymbCreator.h"
#include "SolutionCalculator.h"
#include "NNSolution.h"
#include "NNLog.h"
#include <set>

class NNSolutionCreator: public SolutionCreator{
    
public:
   
    /// Constructor
    NNSolutionCreator();
    
    /// Alternative constructor
    NNSolutionCreator(vector<string>& symbols);
    
    void initialize();
    
    /// creates solution from vector of strings
    virtual void establish_solution(vector<string>& symbols, Dataset* set);
 
     /// creates solution from vector of strings
    virtual void establish_solution(vector<string>& symbols);
 
    /// returns fitness score through evaluation of solution
    virtual float evaluate(Dataset* set);
    
    /// optimize solution by running back propagation
    int optimizeSolution(std::vector<std::string>& symbols, Dataset* set);
    
    /// returns optimized score
    float getOptimizedScore(){return optimized_score;}
       
     virtual Solution* create_new_solution(){
         NNSolution* sol = new NNSolution;
        return sol;}
    
    /// frees memory associated with constants
    void free_solution(){
      for(unsigned int i=0; i<constants.size(); i++){
        delete constants[i];
      }
      constants.clear();
    }
    
    virtual void restrict(vector<string>& vars){}
    
    /// outputs each individual with score of the network 
    virtual float evaluate_with_output(Dataset* set, ostream& os);
    
    /// writes a dot compatible text file representing the network
    virtual void graphical_output(ostream& os, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, bool continmap_used);
 
    vector<int> getGeneIndexes();
    vector<int> getCovarIndexes();
    
    /// Returns worst score
    inline float get_worst(){return calculator->get_worst();}
 
    inline unsigned int get_num_genes(){return getGeneIndexes().size();}
    
    inline unsigned int get_num_covars(){return getCovarIndexes().size();}

    inline int getNumIndsEvaluated(){return nIndsEvaluated;}

    /// Returns symbol that corresponds to start of optimization symbols  
    string getStartOptSymbol(){return startopt;}
    
    std::set<string> getOptIncluded(){return optsymbols;}

    char getLeftOptBound(){return left_opt_bound;}
    char getRightOptBound(){return right_opt_bound;}
    std::set<string> getOptArgSymbols(){return optargsymbols;}

    void detailed_logging();
    unsigned int get_detailed_log();

protected:
    
    void compress_operator(vector<TerminalSymbol*> & postfix_stack,
      vector<TerminalSymbol*>& new_stack);
    
    /// evaluates single individual and returns value for that individual
    virtual float evaluate_ind(Individual* ind);
    
    /// resets number of genotypes and covariates when necessary
    void set_variables(Dataset* set);
    
    /// returns true when ind has complete data for solution
    bool use_ind(Individual* ind, Dataset* set);
    
    TerminalSymbCreator term_holder;
    
    float optimized_score;
    bool terminals_set;
    unsigned int nn_terminal_size, nn_depth;
    vector<TerminalSymbol *> postfix_stack;
    int nIndsEvaluated;
    string startopt;
    std::set<string> optsymbols, optargsymbols;
    
    char left_opt_bound, right_opt_bound;
    
    // vectors for holding list of covariates and genotypes in current solution
    // these can be checked against an ind and the individual can be skipped 
    // when needed
    map<TerminalSymbol*, int> covars, genos;
    
    // vector holds pointers to constants that can be destroyed after all evaluations
    vector<TerminalSymbol* > constants;
};


#endif	/* _NNSOLUTIONCREATOR_H */

