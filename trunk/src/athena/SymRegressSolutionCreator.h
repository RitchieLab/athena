// SymRegressSolutionCreator.h

#ifndef _SYMBREGRESSSOLUTIONCREATOR_H
#define	_SYMBREGRESSSOLUTIONCREATOR_H

#include "NNSolutionCreator.h"

///
/// Inherits from the NNSolutionCreator class.  
/// Adds functionality allowing which assigns worst fitness to any network
/// that doesn't contain all the variables in the grammar.
///

class SymRegressSolutionCreator: public NNSolutionCreator{

  public:
    /// creates solution from vector of strings
    virtual void establish_solution(vector<string>& symbols, Dataset* set);
 
     /// creates solution from vector of strings
    virtual void establish_solution(vector<string>& symbols);

    /// writes a dot compatible text file representing the network
    virtual void graphical_output(ostream& os, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy){}
    
  protected:
  
     /// evaluates single individual and returns value for that individual
    float evaluate_ind(Individual* ind);
    
};

#endif
