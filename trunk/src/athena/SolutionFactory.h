/* 
 * File:   SolutionFactory.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 4:24 PM
 */

#ifndef _SOLUTIONFACTORY_H
#define	_SOLUTIONFACTORY_H

using namespace std;

#include <map>
#include "SolutionCreator.h"
#include "HemannExcept.h"

///
/// Creates and returns solution based on name passed
///
class SolutionFactory{
    
public:
    
    static SolutionCreator* create_solution(string solution_name);
    static SolutionCreator* create_solution(string solution_name, vector<string>& vars);
    
private:
    
    enum SolutionType{
        /// For no match
        MissingSolutionType,
        /// For neural network
        NNSolutionType,
        /// For neural network requiring all variables to be fit
        NNSolutionAllType,
        /// For neural network restricting the networks to only one occurrence of each variable
        NNSolutionOnceType,
        /// For Symbolic Regression 
        SymRegressSolutionType
    };
    
    static map<string, SolutionType> SolutionMap;
    
    /// Sets map for use in Solution creation
    static void setSolutionMap();
    
};

#endif	/* _SOLUTIONFACTORY_H */

