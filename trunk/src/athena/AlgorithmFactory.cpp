#include "AlgorithmFactory.h"
#include "GESymbReg.h"

std::map<string, AlgorithmFactory::AlgorithmType> AlgorithmFactory::AlgorithmMap;


///
/// Creates an algorithm based on the name.  Caller
/// is responsible for freeing the memory.
/// @param alg_name
/// @return new Algorithm
/// @throws AthenaExcept when no algorithm with that name
///
Algorithm* AlgorithmFactory::create_algorithm(string alg_name){
    
    if(AlgorithmMap.empty()){
        set_algorithm_map();
    }
    
    Algorithm* newAlgorithm = NULL;
    switch(AlgorithmMap[alg_name]){
        case NoAlgorithm:
            throw AthenaExcept(alg_name + " is not a valid Algorithm name.");
            break;
        case GENNAlgorithm:
            newAlgorithm = new GENNAlg;
            break;
        case GESymbRegAlgorithm:
            newAlgorithm = new GESymbReg;
            break;
    }
    
    return newAlgorithm;
}


///
/// Establishes map for creating algorithms
/// 
void AlgorithmFactory::set_algorithm_map(){
    AlgorithmMap["GENN"] = GENNAlgorithm;
    AlgorithmMap["GESYMBREG"] = GESymbRegAlgorithm;
}

