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
Algorithm* AlgorithmFactory::createAlgorithm(string algName){
		
		if(AlgorithmMap.empty()){
				setAlgorithmMap();
		}
		
		Algorithm* newAlgorithm = NULL;
		switch(AlgorithmMap[algName]){
				case NoAlgorithm:
						throw AthenaExcept(algName + " is not a valid Algorithm name.");
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
void AlgorithmFactory::setAlgorithmMap(){
		AlgorithmMap["GENN"] = GENNAlgorithm;
		AlgorithmMap["GESYMBREG"] = GESymbRegAlgorithm;
}

