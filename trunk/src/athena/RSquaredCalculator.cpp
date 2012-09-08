#include "RSquaredCalculator.h"

RSquaredCalculator::RSquaredCalculator(){
    reset();
    sstotal = 1;
    name = "R-Squared";
}

void RSquaredCalculator::reset(){
   total_inds_tested = 0;
   squared_error_total = 0.0;
}

///
/// Adds score to running total within object
/// @param score
///
void RSquaredCalculator::add_ind_score(float score, float stat){
    float difference = score - stat;
    squared_error_total += difference * difference;
    total_inds_tested++;
}
