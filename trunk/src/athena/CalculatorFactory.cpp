#include "CalculatorFactory.h"
#include "CalculatorList.h"

#include <iostream>
using namespace std;

std::map<std::string, CalculatorFactory::CalcType> CalculatorFactory::CalcMap;


///
/// Creates and returns calculator
/// @param calc_type
/// @return calculator
/// @throws HemannExcept when no matching calculator to create
///
SolutionCalculator* CalculatorFactory::create_calculator(string calc_name){
   if(CalcMap.empty()){
      set_calc_map();
   }
   
   SolutionCalculator* newSolution;
   
   switch(CalcMap[calc_name]){
       case NoCalcType:
           throw HemannExcept(calc_name + " is not a valid calculation");
           break;
       case MeanSquaredErrType:
           newSolution = new MeanSquaredErrCalculator;
           break;
       case BalanceCalcType:
           newSolution = new BalAccCalculator;
           break;
       case RSquaredType:
           newSolution = new RSquaredCalculator;
           break;
       default:
           throw HemannExcept(calc_name + " is not a valid calculation");
           break;
   }
   
   return newSolution;
}


///
/// Sets the Calculator map
///
void CalculatorFactory::set_calc_map(){
    CalcMap["BALANCEDACC"] = BalanceCalcType;
    CalcMap["RSQUARED"] = MeanSquaredErrType;
    
}

