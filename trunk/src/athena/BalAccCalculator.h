/* 
 * File:   BalAccCalculator.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 2:46 PM
 */
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

#ifndef _BALACCCALCULATOR_H
#define	_BALACCCALCULATOR_H

#include "SolutionCalculator.h"

#include <iostream>

///
/// Calculates balanced accuracy
///
class BalAccCalculator: public SolutionCalculator{
   
public:
    
    BalAccCalculator();
    
    /// resets calculator for new analysis set
    void reset();
    
    /// adds evaluation score to results
    void add_ind_score(float score, float status);
    
    /// returns balanced accuracy
    float get_score(){
      return .5 * (float(caseright) / (caseright + casewrong) + 
      float(controlright) /(controlright + controlwrong));
    }
   bool max_best(){return true;}
   
   /// returns worst score
   float get_worst(){return 0.0;}
private:
    float balanced_acc;
    unsigned int caseright, casewrong, controlright, controlwrong;
};


#endif	/* _BALACCCALCULATOR_H */

