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

map<string, SolutionFactory::SolutionType> SolutionFactory::SolutionMap;

///
/// Function that creates a solution based on the name 
/// @param solution_name name of solution
/// @return pointer to new Solution object
/// @throws AthenaExcept if not a valid solution name
///
SolutionCreator* SolutionFactory::create_solution(string solution_name){
    
  if(SolutionMap.empty()){
    setSolutionMap();
  }

  SolutionCreator* new_solution;
  switch(SolutionMap[solution_name]){
      case MissingSolutionType:
          throw AthenaExcept("No solution matching " + solution_name);
          break;
      case SymRegressSolutionType:
          new_solution = new SymRegressSolutionCreator;
          break;
      case NNSolutionType:
          new_solution = new NNSolutionCreator;
          break;
      case NNSolutionAllType:
          new_solution = new NNSolutionCreatorIncludeAll;
          break;
      case NNSolutionOnceType:
          new_solution = new NNSolutionCreatorIncludeOnce;
          break;
      default:
          throw AthenaExcept("No solution matching " + solution_name); 
  }
  
  return new_solution;
}

///
/// Functions that creates a solution based on the name and passes the string vector
/// for restricting the solutions
/// @param solution_name name of solution
/// @param vars vector that contains strings restricting solution
/// @return pointer to new Solution object
/// @throws AthenaExcept if not a valid solution name
///
SolutionCreator* SolutionFactory::create_solution(string solution_name, vector<string>& vars){
  SolutionCreator* new_creator = create_solution(solution_name);
  new_creator->restrict(vars);
  return new_creator;
}


///
/// Establishes the map for use in creating solutions
/// @return 
///
void SolutionFactory::setSolutionMap(){
  SolutionMap["NN"]=NNSolutionType;
  SolutionMap["NNALL"]=NNSolutionAllType;
  SolutionMap["NNONCE"] = NNSolutionOnceType;
  SolutionMap["SYMBREG"]=SymRegressSolutionType;
}

