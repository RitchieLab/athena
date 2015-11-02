//ScalContinuousTest.cpp

#include "../ScaleContinuous.h"
#include "../Stringmanip.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace data_manage;
using namespace std;

class ScalContinuousTest:public testing::Test{
  protected:
  virtual void SetUp(){
  
    string ids[10] = {"tmp883", "tmp8291", "tmp1", "ind77", "ind701", "ind1010",
      "ind7577", "tmp5", "tmp99", "ind200"};
    float stats[10] = {0.2, 0.3, 0.22, 0.89, 0.48, 0.57, 0.62, 0.84, 0.05, 0.5};  
    float covars[10] = {0.99, -1.3, 2.8, 4.5, -2, 0.83, 1.9, -1.1, 0.05, -4};
    
    // fill holder with individuals
    for(int i=0; i<10; i++){
      Individual *ind = new Individual;
      ind->set_status(stats[i]);
      ind->set_id(Stringmanip::itos(i+1));
      holder1.add_ind(ind);
      
      Individual *indID = new Individual;
      indID->set_status(stats[i]);
      indID->set_id(ids[i]);
      indID->add_covariate(0.5/(i+1));
      // add 12 snps to match the sample map file
      for(int snp=0; snp<12;snp++){
        indID->add_genotype((snp+i)%3);
      }
      indID->add_covariate(covars[i]);
      holder2.add_ind(indID);
    }
  }
  
  Dataholder holder1, holder2;
};



// test for correct normalization of status values
TEST_F(ScalContinuousTest, AdjustStatus){

  ScaleContinuous scale;
  scale.adjust_status(&holder2);

  EXPECT_GT(0.00001, fabs(0.224719 - holder2.get_ind(0)->get_status())) <<
    "ind (index 0) should have status 0.224719";
  EXPECT_GT(0.00001, fabs(0.539326 - holder2.get_ind(4)->get_status())) <<
    "ind (index 0) should have status 0.539326";  
  EXPECT_GT(0.00001, fabs(0.561798 - holder2.get_ind(9)->get_status())) <<
    "ind (index 0) should have status 0.561798";

  EXPECT_GT(0.00001, fabs(0.2 - scale.get_original_status(holder2.get_ind(0)->get_status()))) <<
    "ind (index 0) should have a calculated original status of 0.2";
  EXPECT_GT(0.00001, fabs(0.48 - scale.get_original_status(holder2.get_ind(4)->get_status()))) <<
    "ind (index 0) should have a calculated original status of 0.48"; 
  EXPECT_GT(0.00001, fabs(0.5 - scale.get_original_status(holder2.get_ind(9)->get_status()))) <<
    "ind (index 0) should have a calculated original status of 0.5"; 
}

// check covariate conversion
TEST_F(ScalContinuousTest, AdjustCovar){
  ScaleContinuous scale;

  scale.adjust_contin(&holder2, 1);
  
  EXPECT_GT(0.00001, fabs(0.587059 - holder2.get_ind(0)->get_covariate(1))) <<
    "ind (index 0) should have status 0.587059";
  EXPECT_GT(0.00001, fabs(0.235294 - holder2.get_ind(4)->get_covariate(1))) <<
    "ind (index 0) should have status 0.235294";  
  EXPECT_GT(0.00001, fabs(0 - holder2.get_ind(9)->get_covariate(1))) <<
    "ind (index 0) should have status 0";
  
}

