//BalAccCalculatorTest.cpp

#include "../BalAccCalculator.h"
#include <gtest/gtest.h>

#include <iostream>
using namespace std;

// check that it is correctly indicating
// if max is best
TEST(BalAccCalculatorTest, MaxBest){
  BalAccCalculator balacc;
  EXPECT_TRUE(balacc.max_best()) << "Max is best for balanced accuracy";
}

// Check worst score
TEST(BalAccCalculatorTest, GetWorst){
  BalAccCalculator balacc;
  EXPECT_EQ(0.0,balacc.get_worst()) << "Worst should be zero for balanced accuracy calculator";
}

// Add individuals and calculate the balanced accuracy
TEST(BalAccCalculatorTest, GetScore){
  BalAccCalculator balacc;
  
  for(int i=0; i<50; i++){
    balacc.add_ind_score(0,0);
  }
  for(int i=0; i<50; i++){
    balacc.add_ind_score(1,0);
  }
  for(int i=0; i<40; i++){
    balacc.add_ind_score(1,1);
  }
  for(int i=0; i<30; i++){
    balacc.add_ind_score(0,1);
  }
  // check calculation
  EXPECT_GT(0.00001, fabs(balacc.get_score() - 0.535714)) << 
    "Balanced accuracy score is incorrect";
}

