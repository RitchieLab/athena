/* 
 * File:   CalculatorFactory.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 4:00 PM
 */

#ifndef _CALCULATORFACTORY_H
#define	_CALCULATORFACTORY_H

#include "SolutionCalculator.h"
#include "HemannExcept.h"
#include <map>
#include <string>

///
/// Returns calculator depending on string passed
///
class CalculatorFactory{
    
public:
        
    static SolutionCalculator* create_calculator(std::string calc_name);
    
private:
    
    static void set_calc_map();
    
    enum CalcType{
        NoCalcType,
        BalanceCalcType,
        MeanSquaredErrType,
        RSquaredType
    };
    
    static std::map<std::string, CalcType> CalcMap;

};



#endif	/* _CALCULATORFACTORY_H */

