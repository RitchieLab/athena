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
#include "CalculatorFactory.h"
#include "CalculatorList.h"

#include <iostream>
using namespace std;

std::map<std::string, CalculatorFactory::CalcType> CalculatorFactory::CalcMap;


///
/// Creates and returns calculator
/// @param calc_type
/// @return calculator
/// @throws AthenaExcept when no matching calculator to create
///
SolutionCalculator* CalculatorFactory::create_calculator(string calc_name){
   if(CalcMap.empty()){
      set_calc_map();
   }
   
   SolutionCalculator* newSolution;
   
   switch(CalcMap[calc_name]){
       case NoCalcType:
           throw AthenaExcept(calc_name + " is not a valid calculation");
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
           throw AthenaExcept(calc_name + " is not a valid calculation");
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

