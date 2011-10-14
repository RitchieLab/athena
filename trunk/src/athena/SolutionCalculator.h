/* 
 * File:   SolutionCalculator.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 2:38 PM
 */

#ifndef _SOLUTIONCALCULATOR_H
#define	_SOLUTIONCALCULATOR_H

///
/// Base class for calculation of solution scores.
/// For example, calculates Balanced accuracy for binary outcomes and
/// mean squared error for continuous outcomes.
///

class SolutionCalculator{
    
public:
    
    virtual ~SolutionCalculator(){}

    virtual void add_ind_score(float score, float status)=0;
    
    virtual float get_score()=0;
    
    virtual void reset()=0;
    
    virtual bool max_best()=0;
    
    virtual float get_worst()=0;
    
    virtual float get_constant(){return 0.0;}
    
    /// Used when a sub class needs a constant value for calculations as in RSquared
    virtual void set_constant(float constant){}
    
private:
    
};



#endif	/* _SOLUTIONCALCULATOR_H */

