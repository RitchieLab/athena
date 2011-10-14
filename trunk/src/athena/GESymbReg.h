/* 
 * File:   GESymbReg.h
 * Author: dudeksm
 *
 * Created on August 10, 2010, 2:03 PM
 */

#ifndef _GESYMBREG_H
#define	_GESYMBREG_H

#include "GENNAlg.h"

///
/// Inherits from the GENNAlg class.
/// Contains code for running GE to produce symbolic regression equations
/// 
 
class GESymbReg: public GENNAlg{
  public:
  
    /// Set the parameters for the algorithm
    void set_params(AlgorithmParams& alg_params, int numExchanges, int numGenos, int numContin);  
    
    /// Runs a step of the algorithm
    virtual void step();
  
  private:
    
    /// Sets GA for run
    void set_ga_params();
  
};


#endif

