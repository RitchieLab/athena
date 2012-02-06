//DatasetTest.cpp

#include "../Dataset.h"
#include "../Stringmanip.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace data_manage;
using namespace std;

class DatasetTest:public testing::Test{
  protected:
  virtual void SetUp(){
  
    string ids[10] = {"tmp883", "tmp8291", "tmp1", "ind77", "ind701", "ind1010",
      "ind7577", "tmp5", "tmp99", "ind200"};
    float stats[10] = {0.2, 0.3, 0.22, 0.89, 0.48, 0.57, 0.62, 0.84, 0.05, 0.5};  
      
    // fill holder with individuals
    for(int i=0; i<10; i++){
      Individual *ind = new Individual;
      ind->set_status(stats[i]);
      ind->set_id(Stringmanip::itos(i+1));
      set1.add_ind(ind);
      
      Individual *indID = new Individual;
      indID->set_status(stats[i]);
      indID->set_id(ids[i]);
      indID->add_covariate(0.5/(i+1));
      indID->add_genotype(i%2);
      indID->add_genotype(i%3);
      set2.add_ind(indID);
    }
  }
  
  virtual void TearDown(){
    for(unsigned int i=0; i<set1.num_inds(); i++){
      delete set1[i];
    }
    for(unsigned int i=0; i<set2.num_inds(); i++){
      delete set2[i];
    }
  }
  
  Dataset set1, set2;
};

// Constructor test 
TEST_F(DatasetTest, Constructor){
  Dataset set3(-9999, 3);
  
  EXPECT_EQ(3, set3.get_missing_genotype());
  EXPECT_GT(0.000001, fabs(set3.get_missing_covalue()-(-9999)));
}

// Check num covariates
TEST_F(DatasetTest, NumCovariates){
  EXPECT_EQ(1, int(set2.num_covariates()));
  EXPECT_EQ(0, int(set1.num_covariates()));
}

// Check number of genotypes
TEST_F(DatasetTest, NumGenotypes){
  EXPECT_EQ(2, int(set2.num_genos()));
  EXPECT_EQ(0, int(set1.num_genos()));
}


// Test the ss total calculation (used in r-squared calculations)
TEST_F(DatasetTest, SSTotalCalc){
  set2.calc_sstotal();
  EXPECT_GT(0.000001, fabs(set2.get_sstotal()-0.68741));
}

