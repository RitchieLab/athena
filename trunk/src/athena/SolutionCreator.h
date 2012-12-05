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
 * File:   SolutionCreator.h
 * Author: dudeksm
 *
 * Created on December 1, 2008, 4:10 PM
 */

#ifndef _SOLUTIONCREATOR_H
#define	_SOLUTIONCREATOR_H

#include "Solution.h"
#include "CalculatorFactory.h"
#include "AlgorithmLog.h"
#include "Structs.h"
#include <Dataset.h>
#include <set>

#include <iostream>

using namespace data_manage;
///
/// Solution is the abstract base class for solutions in the HEMANN
/// system.  Solutions are the structure created by the algorithms.
/// An example is NeuralNetwork.  Future possible solutions would be
/// support vector machines, etc.
///
class SolutionCreator{
    
public:

    virtual ~SolutionCreator(){}
    
    /// creates solution from a vector of strings
    virtual void establish_solution(std::vector<std::string>& symbols, Dataset* set)=0;

    /// creates solution from a vector of strings
    virtual void establish_solution(std::vector<std::string>& symbols)=0;
    
    /// returns fitness score through evaluation of solution
    virtual float evaluate(Dataset* set)=0;
    
    /// adds a solution calculator for determining fitness
    virtual void set_calculator(std::string calc_type){
        if(calculator != NULL)
          delete calculator;
        calculator = CalculatorFactory::create_calculator(calc_type);
    }
    
    virtual void free_solution()=0;
    
    /// returns solution
    Solution* get_solution(){return sol;}
    
    /// returns fitness
    float fitness(){return sol_fitness;}
    
    /// set fitness
    void fitness(float fit){sol_fitness = fit;}
    
    /// optimize solution
    virtual int optimizeSolution(std::vector<std::string>& symbols, Dataset* set)=0;
    
    /// returns vector of optimized values for use in adapting the original
    vector<float> getOptimizedValues(){return opt_values;}
    
    /// returns optimized score
    virtual float getOptimizedScore()=0;
    
    /// returns vector of structs containing optimized values stored as strings
    vector<symbVector> getOptimizedSymbols(){return opt_val_symbols;}
    
    /// returns a blank solution of appropriate type
    virtual Solution* create_new_solution()=0;
    
    virtual void set_calculator_constant(float constant){calculator->set_constant(constant);}
    virtual float get_calculator_constant(){return calculator->get_constant();}
    
    virtual bool max_best(){return calculator->max_best();}
    
    virtual float get_worst(){return calculator->get_worst();}
    
    virtual void restrict(vector<string>& restrictions)=0;
    
    virtual float evaluate_with_output(Dataset* set, ostream& os)=0;
    
    /// writes a graphical or file that can be converted to a graphic representation of the solution
    virtual void graphical_output(ostream& os, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, bool continmap_used)=0;
    
    virtual std::string graphicExt(){return graphic_extension;}
    
    virtual unsigned int get_num_genes()=0;
    
    virtual unsigned int get_num_covars()=0;
    
    virtual vector<int> getGeneIndexes()=0;
    virtual vector<int> getCovarIndexes()=0;
    
    virtual string getStartOptSymbol()=0;
    virtual std::set<string> getOptIncluded()=0;

    virtual char getLeftOptBound()=0;
    virtual char getRightOptBound()=0;
    virtual std::set<string> getOptArgSymbols()=0;
    
    virtual int getNumIndsEvaluated()=0;
    
    virtual void detailed_logging()=0;
    virtual unsigned int get_detailed_log()=0;
    
    std::string calculator_name(){return calculator->get_name();}
    
protected:
    Solution* sol;
    SolutionCalculator * calculator;
    std::string graphic_extension;
    float sol_fitness;
    vector<float> opt_values;
    vector<symbVector> opt_val_symbols;
};



#endif	/* _SOLUTIONCREATOR_H */

