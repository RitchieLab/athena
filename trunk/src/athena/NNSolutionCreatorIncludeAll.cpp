//NNSolutionCreatorIncludeAll.cpp

#include <iostream>
using namespace std;

#include "NNSolutionCreatorIncludeAll.h"

///
/// Establishes solution.  If not all variables are in set, doesn't bother
/// evaluating model at all.  Sets flag so that worst fitness will be returned.
///
void NNSolutionCreatorIncludeAll::establish_solution(vector<string>& symbols, Dataset* set){
  
  NNSolutionCreator::establish_solution(symbols, set);
    
  if(!vars_set){
    required_symbols(req_symbols, set);
    vars_set = true;
  }

  // check genos
  all_vars_included = true;
  vector<TerminalSymbol*>::iterator iter;
  for(iter=req_genos.begin(); iter != req_genos.end(); *iter++){
    if(genos.find(*iter) == genos.end()){
      all_vars_included = false;
      return;
    }
  }
  
  // check covars
  for(iter=req_covars.begin(); iter != req_covars.end(); *iter++){
    if(covars.find(*iter) == covars.end()){
      all_vars_included = false;
      return;
    }
  }
}

///
/// returns fitness score through evaluation of solution
///
float NNSolutionCreatorIncludeAll::evaluate(Dataset* set){

  if(all_vars_included)
    return NNSolutionCreator::evaluate(set);
  else
    return calculator->get_worst();
}


///
/// Takes list of symbols and Dataset to establish list of terminals that must 
/// be part of the neural network for it to receive 
///
void NNSolutionCreatorIncludeAll::required_symbols(vector<string>& symbols, Dataset* set){
  
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
