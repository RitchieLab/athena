//CrossValidatorTest.cpp

#include "../CrossValidator.h"
#include "../Dataholder.h"
#include "../MDRFileHandler.h"
#include <gtest/gtest.h>

using namespace data_manage;
using namespace std;

class CrossValidatorTest:public testing::Test{
  protected:
  virtual void SetUp(){
    // load dataholder
    MDRFileHandler mdrf;
    mdrf.parse_file("sample_data/small_sample.txt", &holder, -1, -9, false);
    srand(7);
  }
  Dataholder holder;
};

// check 2-way split
TEST_F(CrossValidatorTest, SplitData){
  
  CrossValidator cv;
  
  unsigned int num_cv=5;
  
  // split data and then check
  CVSet set = cv.split_data(num_cv, &holder);
  
  EXPECT_EQ(num_cv, set.num_intervals());
  
  // check one set to see size of datasets in training and testing
  for(unsigned int i=0; i<num_cv; i++){
    CVInterval inter = set.get_interval(i);
    EXPECT_EQ(80, int(inter.get_training().num_inds())) << "Training set from interval (index " << i << ") should have 80 inds in it";
    EXPECT_EQ(20, int(inter.get_testing().num_inds())) << "Testing set from interval (index " << i << ") should have 20 inds in it";
  }
  
}

