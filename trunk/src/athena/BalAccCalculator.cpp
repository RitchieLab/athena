#include "BalAccCalculator.h"


BalAccCalculator::BalAccCalculator(){
    reset();
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

