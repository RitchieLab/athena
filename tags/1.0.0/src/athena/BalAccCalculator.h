/* 
 * File:   BalAccCalculator.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 2:46 PM
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

