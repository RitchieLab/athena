/* 
 * File:   AlgorithmFactory.h
 * Author: dudeksm
 *
 * Created on November 10, 2008, 2:09 PM
 */

#ifndef _ALGORITHMFACTORY_H
#define	_ALGORITHMFACTORY_H

#include "Algorithm.h"

class AlgorithmFactory{
    
public:
    /// Creates and returns the Algorithm
    static Algorithm* create_algorithm(std::string alg_name);
    
    /// Enumeration for algorithm types
    enum AlgorithmType{
        /// Enum for missing
        NoAlgorithm,
        /// Enum for GENN
        GENNAlgorithm,
        /// Enum for Symbolic Regression
        GESymbRegAlgorithm
    };
    
private:
    
    /// Sets map for Algorithm creation
    static void set_algorithm_map();
    
    static std::map<std::string, AlgorithmType> AlgorithmMap;
    
};



#endif	/* _ALGORITHMFACTORY_H */


