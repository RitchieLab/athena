//MeanSquaredErrCalculatorTest.cpp

#include "../MeanSquaredErrCalculator.h"
#include <gtest/gtest.h>

#include <iostream>
using namespace std;

// check that it is correctly indicating
// if max is best
TEST(MeanSquaredErrCalculatorTest, MaxBest){
  MeanSquaredErrCalculator meanerr;
  EXPECT_FALSE(meanerr.max_best()) << "Max is worst for mean squared error";
}

// Check worst score
TEST(MeanSquaredErrCalculatorTest, GetWorst){
  MeanSquaredErrCalculator meanerr;
  EXPECT_EQ(1000000,meanerr.get_worst()) << "Worst should be zero for balanced accuracy calculator";
}

// Add individuals and calculate the balanced accuracy
// Also check sstotal calculation
TEST(MeanSquaredErrCalculatorTest, GetScore){
  MeanSquaredErrCalculator meanerr;
  
  srand(7);
  
  for(int i=0; i<50; i++){
    meanerr.add_ind_score(rand()/float(RAND_MAX),0);
  }
  
  for(int i=0; i<50; i++){
    meanerr.add_ind_score(rand()/float(RAND_MAX),0);
  }
  
  for(int i=0; i<40; i++){
    meanerr.add_ind_score(rand()/float(RAND_MAX),1);
  }
  
  for(int i=0; i<30; i++){
    meanerr.add_ind_score(rand()/float(RAND_MAX),1);
  }
  
  // check calculation
  EXPECT_GT(0.00001, fabs(meanerr.get_score() - 0.332651)) << 
    "Mean squared error score is incorrect";
  EXPECT_GT(0.0001, fabs(meanerr.get_constant() - 41.1765)) <<
    "sstotal value incorrect within mean squared error calculation";

  // test reset functionality
  meanerr.reset();
  EXPECT_EQ(1000000, meanerr.get_score()) << "Reset mean squared error should be 1000000 for worst possible";
  EXPECT_EQ(0, meanerr.get_constant()) << "Constant (sstotal) should be 0";

}


