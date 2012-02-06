//IndividualTest.cpp

#include "../Individual.h"
#include <gtest/gtest.h>
#include <vector>
using namespace data_manage;
using namespace std;

class IndividualTest:public testing::Test{
  protected:
  virtual void SetUp(){
    ind2.set_status(0.5);    
    ind2.add_genotype(2);
    ind2.add_genotype(1);
    ind2.add_genotype(0);
    ind3 = ind2;
    ind3.add_covariate(0.2);
    ind3.add_covariate(0.1);
    
  }
  Individual ind1, ind2, ind3;
};


TEST_F(IndividualTest, AddGenotype){
  ind1.add_genotype(2);
  ind1.add_genotype(1);
  EXPECT_EQ(2, ind1.get_genotype(0));
  EXPECT_EQ(1, ind1.get_genotype(1));
}

TEST_F(IndividualTest, SetAllGenotypes){
  vector<char> genos;
  genos.push_back(1);
  genos.push_back(2);
  genos.push_back(0);
  
  ind1.set_all_genotypes(genos);
  
  EXPECT_EQ(3, int(ind1.num_genotypes()));
  EXPECT_EQ(2, ind1.get_genotype(1));
  
  ind2.set_all_genotypes(genos);
  EXPECT_EQ(3, int(ind2.num_genotypes()));
  EXPECT_EQ(2, ind2.get_genotype(1));
}


TEST_F(IndividualTest, DefaultConstructor){
  EXPECT_EQ(0, int(ind1.num_genotypes()));
}


TEST_F(IndividualTest, CopyTest){
  Individual copyind = ind3;
  ASSERT_EQ(ind3.num_covariates(), copyind.num_covariates());
  for(unsigned int i=0; i<ind3.num_covariates(); i++){
    EXPECT_EQ(ind3.get_covariate(i), copyind.get_covariate(i));
  }
  ASSERT_EQ(ind3.num_genotypes(), copyind.num_genotypes());
  for(unsigned int i=0; i<ind3.num_genotypes(); i++){
    EXPECT_EQ(ind3.get_genotype(i), copyind.get_genotype(i));
  }
  
}


TEST_F(IndividualTest, CovariateTest){
  EXPECT_EQ(2, int(ind3.num_covariates()));
  EXPECT_GT(0.0000001, fabs(ind3.get_covariate(0)-0.2));
}
