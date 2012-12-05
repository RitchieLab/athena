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

