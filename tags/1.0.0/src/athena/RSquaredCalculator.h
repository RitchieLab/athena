/* 
 * File:   RSquared.h
 * Author: dudeksm
 *
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

