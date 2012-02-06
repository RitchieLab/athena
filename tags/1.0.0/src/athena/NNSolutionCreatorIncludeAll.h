// NNSolutionCreatorIncludeAll.h

#ifndef _NNSOLUTIONCREATORINCLUDEALL_H
#define	_NNSOLUTIONCREATORINCLUDEALL_H

#include "NNSolutionCreator.h"

///
/// Inherits from the NNSolutionCreator class.  
/// Adds functionality allowing which assigns worst fitness to any network
/// that doesn't contain all the variables in the grammar.
///

class NNSolutionCreatorIncludeAll: public NNSolutionCreator{

  public:
  
    NNSolutionCreatorIncludeAll(){vars_set = false;}
  
    /// creates solution from vector of strings
    void establish_solution(vector<string>& symbols, Dataset* set);
    
    /// returns fitness score through evaluation of solution
    float evaluate(Dataset* set);
    
    /// passes vector of symbols that must be part of the neural network
    void required_symbols(vector<string>& symbols, Dataset* set);
    
    /// gets list of required symbols
    void restrict(vector<string>& vars){req_symbols = vars;}
    
    
  private:
    
    vector<TerminalSymbol*> req_genos, req_covars;
    vector<string> req_symbols;
    bool all_vars_included, vars_set;
    
};

#endif
