/* 
 * File:   RSquared.h
 * Author: dudeksm
 *
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

#ifndef _RSQUAREDCALCULATOR_H
#define	_RSQUAREDCALCULATOR_H

#include "SolutionCalculator.h"


#include <iostream>

///
/// Calculates mean squared error
///
class RSquaredCalculator: public SolutionCalculator{
    
    public:
    
    RSquaredCalculator();
    
    /// resets calculator for new analysis set
    void reset();
    
    /// adds evaluation score to results
    void add_ind_score(float score, float status);
    
    /// returns r-squared
    float get_score(){
        return 1-(squared_error_total / sstotal);
    }
    
    /// returns false so that the best is the smallest
    bool max_best(){return true;}
    
    /// returns worst score
    float get_worst(){return -100;}
    
    /// Sets the SSTotal which is the square of the differences of the mean status and each ind status
    void set_constant(float constant){sstotal = constant;}
    
private:
    
    float squared_error_total, sstotal;
    unsigned int total_inds_tested;
    
    
};




#endif	/* _MEANSQUAREDERRCALCULATOR_H */

