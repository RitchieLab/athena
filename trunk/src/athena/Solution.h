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
 * File:   Solution.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 3:35 PM
 */

#ifndef _SOLUTION_H
#define	_SOLUTION_H

#include <Dataset.h>
#include <string>
#include <vector>
#include <iostream>
#include <Dataholder.h>

using namespace data_manage;
///
/// Solution is the abstract base class for solutions in the HEMANN
/// system.  Solutions are the structure created by the algorithms.
/// An example is NeuralNetwork.  Future possible solutions would be
/// support vector machines, etc.
///
class Solution{
    
public:
 
    /// default constructor
    Solution(){solution_name = "";}
    
    /// Named constructor
    Solution(std::string name){set_name(name);}
    
    /// Destructor
    virtual ~Solution(){}
    
    /// Clones the current solution
    virtual Solution* clone();
    
    /// return number of symbols
    unsigned int get_num_symbols(){return symbols.size();}
    
    /// operator [] for accessing symbols
    inline std::string& operator[](unsigned int index){return symbols[index];}
   
    /// sets solution namem
    void set_name(std::string name){solution_name = name;}
    
    /// returns solution name (type)
    std::string get_name(){return solution_name;}
    
    /// set symbols
    void set_symbols(std::vector<std::string>& sym){ symbols=sym;}
    
    /// returns symbols
    std::vector<std::string>& get_symbols(){return symbols;}
    
    /// returns fitness
    float fitness(){return sol_fitness;}
    
    /// set fitness
    void fitness(float fit){sol_fitness = fit;}
    
    /// returns test score
    float testval(){return test_score;}
    
    /// sets test score
    void testval(float val){test_score = val;}
    
    /// adjusts output so that 
    virtual void adjust_dummy_encoding(){};
    
    /// outputs as a text file
    virtual void output_solution(std::ostream& os);
    
    /// outputs as a graphviz compatible file
    virtual void output_graph(std::ostream& os);

    /// outputs a more human-readable version of the network
    virtual void output_clean(std::ostream& os, data_manage::Dataholder& data, 
      bool map_used, bool ott_dummy, bool continmap_used);

    void copy(Solution* other);
    
    /// returns the genotypes present in the solution
    virtual std::vector<int> get_genotypes(bool dummy_encoded = true){
        std::vector<int> blank;
        return blank;}
    
    /// returns covariates present in the solution
    virtual std::vector<int> get_covariates(){
        std::vector<int> blank;
        return blank;}
    
    /// Adjusts output of the scores when needed (e.g. meansquared to rsquared)
    virtual void adjust_score_out(Dataset* train_set, Dataset* test_set){}
    
    /// Adjusts output of the scores when needed
    virtual void adjust_score_out(Dataset* train_set){}
    
    /// Adjusts score passed and returns value
    virtual float adjust_score_out(float score, int nIndsTested, float constant){return score;}  

protected:
    
    std::vector<std::string> symbols;
    float sol_fitness, test_score;
      
private:
    
    std::string solution_name;
};


#endif	/* _SOLUTION_H */

