/* 
 * File:   MeanSquaredErrCalculator.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 3:45 PM
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

#ifndef _MEANSQUAREDERRCALCULATOR_H
#define	_MEANSQUAREDERRCALCULATOR_H

#include "SolutionCalculator.h"


#include <iostream>
#include <vector>

///
/// Calculates mean squared error
///
class MeanSquaredErrCalculator: public SolutionCalculator{
    
    public:
    
    MeanSquaredErrCalculator();
    
    /// resets calculator for new analysis set
    void reset();
    
    /// adds evaluation score to results
    void add_ind_score(float score, float status);
    
    /// returns balanced accuracy
    float get_score(){
        if(total_inds_tested == 0)
          return get_worst();
        else
          return squared_error_total / total_inds_tested;
    }
    
    /// returns false so that the best is the smallest
    bool max_best(){return false;}
    
    /// returns worst score
    float get_worst(){return 1000000;}
    /// sets the sstotal for use in calculating rsquared
    void set_constant(float constant){sstotal = constant;}
    /// returns the sstotal for use in calculating rsquared
    float get_constant();

    /// returns rsquared for display to user
    float outputScore(float sc){return 1-(sc/sstotal);}

    
private:
    
    float squared_error_total, sstotal;
    std::vector<float> stat_total;
    unsigned int total_inds_tested;
    
    
};




#endif	/* _MEANSQUAREDERRCALCULATOR_H */

