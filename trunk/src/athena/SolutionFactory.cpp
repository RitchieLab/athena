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
#include "SolutionFactory.h"
#include "NNSolutionCreator.h"
#include "NNSolutionCreatorIncludeAll.h"
#include "NNSolutionCreatorIncludeOnce.h"
#include "SymRegressSolutionCreator.h"

#include <iostream>

map<string, SolutionFactory::SolutionType> SolutionFactory::solutionMap;

///
/// Function that creates a solution based on the name 
/// @param solutionName name of solution
/// @return pointer to new Solution object
/// @throws AthenaExcept if not a valid solution name
///
SolutionCreator* SolutionFactory::createSolution(string solutionName){
		
	if(solutionMap.empty()){
		setSolutionMap();
	}

	SolutionCreator* newSolution;
	switch(solutionMap[solutionName]){
			case MissingSolutionType:
					throw AthenaExcept("No solution matching " + solutionName);
					break;
			case SymRegressSolutionType:
					newSolution = new SymRegressSolutionCreator;
					break;
			case NNSolutionType:
					newSolution = new NNSolutionCreator;
					break;
			case NNSolutionAllType:
					newSolution = new NNSolutionCreatorIncludeAll;
					break;
			case NNSolutionOnceType:
					newSolution = new NNSolutionCreatorIncludeOnce;
					break;
			default:
					throw AthenaExcept("No solution matching " + solutionName); 
	}
	
	return newSolution;
}

///
/// Functions that creates a solution based on the name and passes the string vector
/// for restricting the solutions
/// @param solutionName name of solution
/// @param vars vector that contains strings restricting solution
/// @return pointer to new Solution object
/// @throws AthenaExcept if not a valid solution name
///
SolutionCreator* SolutionFactory::createSolution(string solutionName, vector<string>& vars){
	SolutionCreator* newCreator = createSolution(solutionName);
	newCreator->restrict(vars);
	return newCreator;
}


///
/// Establishes the map for use in creating solutions
/// @return 
///
void SolutionFactory::setSolutionMap(){
	solutionMap["NN"]=NNSolutionType;
	solutionMap["NNALL"]=NNSolutionAllType;
	solutionMap["NNONCE"] = NNSolutionOnceType;
	solutionMap["SYMBREG"]=SymRegressSolutionType;
}

