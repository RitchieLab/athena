#include "AlgorithmFactory.h"
#include "GESymbReg.h"
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

