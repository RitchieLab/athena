#include "NNSolutionCreatorIncludeOnce.h"
#include <set>

#include <iostream>
using namespace std;

///
/// Establishes solution.  If not all variables are in set, doesn't bother
/// evaluating model at all.  Sets flag so that worst fitness will be returned.
///
void NNSolutionCreatorIncludeOnce::establish_solution(vector<string>& symbols, Dataset* set){
  
  NNSolutionCreator::establish_solution(symbols, set);
    
  if(!vars_set){
    required_symbols(req_symbols, set);
    vars_set = true;
  }

  // check genos
  all_vars_included = true;

  // must have same number in required and genos in this solution
  if((req_genos.size() != genos.size()) || req_covars.size() != covars.size()){
    all_vars_included = false;
    return;
  }
  // need to make sure that no geno or covar is duplicated
  map<TerminalSymbol*, int>::iterator iter;
  
  for(iter = genos.begin(); iter != genos.end(); iter++){
    if(iter->second > 1){
      all_vars_included = false;
      return;
    }
  }
  
  for(iter = covars.begin(); iter != covars.end(); iter++){
    if(iter->second > 1){
      all_vars_included = false;
      return;
    }
  }

}

///
/// returns fitness score through evaluation of solution
///
float NNSolutionCreatorIncludeOnce::evaluate(Dataset* set){

  if(all_vars_included)
    return NNSolutionCreator::evaluate(set);
  else
    return calculator->get_worst();
}


///
/// Takes list of symbols and Dataset to establish list of terminals that must 
/// be part of the neural network for it to receive a score
///
void NNSolutionCreatorIncludeOnce::required_symbols(vector<string>& symbols, Dataset* set){
  
  if(!terminals_set)
    NNSolutionCreator::set_variables(set);
    
  // use symbols to get pointers that can be used to check for inclusion of all variables
  req_covars.clear();
  req_genos.clear();
  vector<string>::iterator iter;
  TerminalSymbol* var;
  for(iter=symbols.begin(); iter != symbols.end(); iter++){
    var = term_holder.get_term(*iter);
    if(var->get_term_type() == TerminalSymbol::Genotype)
      req_genos.push_back(var);
    else if(var->get_term_type() == TerminalSymbol::Covariate)
      req_covars.push_back(var);
  }
}
