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
#include "BalAccCalculator.h"


BalAccCalculator::BalAccCalculator(){
    reset();
    name = "Balanced Accuracy";
}


void BalAccCalculator::reset(){
    caseright=0;
    casewrong=0;
    controlright=0;
    controlwrong=0;
}

#include<iostream>
using namespace std;

///
/// Adds score to running total within object
/// @param score
///
void BalAccCalculator::add_ind_score(float score, float stat){
   
    unsigned int result = score > 0.5?1:0;
    unsigned int status = (unsigned int)stat;
    
    if(result != status){
        if(status)
            casewrong++;
        else
            controlwrong++;
    }
    else if(status)
        caseright++;
    else
        controlright++;
}

