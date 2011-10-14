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
/// @throws HEMannExcept if not a valid solution name
///
SolutionCreator* SolutionFactory::create_solution(string solution_name){
    
  if(SolutionMap.empty()){
    setSolutionMap();
  }

  SolutionCreator* new_solution;
  switch(SolutionMap[solution_name]){
      case MissingSolutionType:
          throw HemannExcept("No solution matching " + solution_name);
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
          throw HemannExcept("No solution matching " + solution_name); 
  }
  
  return new_solution;
}

///
/// Functions that creates a solution based on the name and passes the string vector
/// for restricting the solutions
/// @param solution_name name of solution
/// @param vars vector that contains strings restricting solution
/// @return pointer to new Solution object
/// @throws HEMannExcept if not a valid solution name
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

