/* 
 * File:   SolutionCalculator.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 2:38 PM
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

#ifndef _SOLUTIONCALCULATOR_H
#define	_SOLUTIONCALCULATOR_H

#include <string>

///
/// Base class for calculation of solution scores.
/// For example, calculates Balanced accuracy for binary outcomes and
/// mean squared error for continuous outcomes.
///

class SolutionCalculator{
    
public:
    
    SolutionCalculator(){name = "Calculator";}
    
    virtual ~SolutionCalculator(){}

    virtual void add_ind_score(float score, float status)=0;
    
    virtual float get_score()=0;
    
    virtual void reset()=0;
    
    virtual bool max_best()=0;
    
    virtual float get_worst()=0;
    
    virtual float get_constant(){return 0.0;}
    
    /// Used when a sub class needs a constant value for calculations as in RSquared
    virtual void set_constant(float constant){}
    
    std::string get_name(){return name;}
    
protected:
    std::string name;
    
};



#endif	/* _SOLUTIONCALCULATOR_H */

