//CVIntervalTest.cpp

#include "../CVInterval.h"
#include "../Stringmanip.h"
#include <gtest/gtest.h>
#include <iostream>
using namespace std;

using namespace data_manage;

class CVIntervalTest:public testing::Test{
  protected:
  virtual void SetUp(){
  
    string ids[10] = {"tmp883", "tmp8291", "tmp1", "ind77", "ind701", "ind1010",
      "ind7577", "tmp5", "tmp99", "ind200"};
    float stats[10] = {0.2, 0.3, 0.22, 0.89, 0.48, 0.57, 0.62, 0.84, 0.05, 0.5};  
      
    // fill sets with individuals for tests
    for(int i=0; i<10; i++){
      Individual *ind = new Individual;
      ind->set_status(stats[i]);
      ind->set_id(ids[i]);
      set1.add_ind(ind);
    }
    
    for(int i=0; i<5; i++){
      Individual *indID = new Individual;
      indID->set_status(stats[i]);
      indID->set_id(Stringmanip::itos(i+1));
      indID->add_covariate(0.5/(i+1));
      indID->add_genotype(i%2);
      indID->add_genotype(i%3);
      set2.add_ind(indID);
    }
    
    for(int i=0; i<5; i++){
      Individual *indID = new Individual;
      indID->set_status(stats[i]);
      indID->set_id(Stringmanip::itos(i+1+5));
      indID->add_covariate(0.5/(i+1));
      indID->add_genotype(i%2);
      indID->add_genotype(i%3);
      set3.add_ind(indID);      
    }

  }
  
  virtual void TearDown(){
    for(unsigned int i=0; i<set1.num_inds(); i++){
      delete set1[i];
    }
    for(unsigned int i=0; i<set2.num_inds(); i++){
      delete set2[i];
    }
    for(unsigned int i=0; i<set3.num_inds(); i++){
      delete set3[i];
    }
  }
  
  Dataset set1, set2, set3;
};


// Test adding sets and retrieving
TEST_F(CVIntervalTest, AddSet){
  CVInterval cvint;
  cvint.add_set(set1);
  cvint.add_set(set2);
  cvint.add_set(set3);
  
  EXPECT_EQ(3, int(cvint.num_sets()));
  Dataset& set = cvint.get_training();
  EXPECT_EQ(10, int(set.num_inds()));
  set = cvint.get_testing();
  EXPECT_EQ(5, int(set.num_inds()));
  set = cvint.get_validation();
  EXPECT_EQ(5, int(set.num_inds()));
  
}


// Test for setting up blank sets
TEST_F(CVIntervalTest, NumSets){
  CVInterval cvint;
  cvint.num_sets(3);
  EXPECT_EQ(3, int(cvint.num_sets()));
  Dataset& set = cvint.get_training();
  EXPECT_EQ(0, int(set.num_inds()));
  set = cvint.get_testing();
  EXPECT_EQ(0, int(set.num_inds()));
  set = cvint.get_validation();
  EXPECT_EQ(0, int(set.num_inds()));  
  
}

// Get sets by index
TEST_F(CVIntervalTest, IndexSets){

  CVInterval cvint;
  cvint.add_set(set1);
  cvint.add_set(set2);
  cvint.add_set(set3);
  
  EXPECT_EQ(3, int(cvint.num_sets()));
  Dataset& set = cvint.get_set(0);
  EXPECT_EQ(10, int(set.num_inds()));
  set = cvint.get_set(1);
  EXPECT_EQ(5, int(set.num_inds()));
  set = cvint.get_set(2);
  EXPECT_EQ(5, int(set.num_inds()));  

}


