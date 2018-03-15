//DataholderTest.cpp

#include "../Dataholder.h"
#include "../Stringmanip.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace data_manage;
using namespace std;

class DataholderTest:public testing::Test{
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
      holder1.add_ind(ind);
      
      Individual *indID = new Individual;
      indID->set_status(stats[i]);
      indID->set_id(ids[i]);
      indID->add_covariate(0.5/(i+1));
      indID->add_genotype(i%2);
      indID->add_genotype(i%3);
      holder2.add_ind(indID);
    }
  }
  
  Dataholder holder1, holder2;
};


// Check overloaded [] operator
TEST_F(DataholderTest, OverloadedOp){
  string id("tmp883");
  EXPECT_EQ(id, holder2[0]->get_id());
  id = "ind200";
  EXPECT_EQ(id, holder2[9]->get_id());
  id = "tmp5";
  EXPECT_EQ(id, holder2[7]->get_id());
}

// Test covar name
TEST_F(DataholderTest, CovarName){
  holder2.add_covar_name("covar1");
  EXPECT_EQ("covar1", holder2.get_covar_name(0));
}


// Test number of inds
TEST_F(DataholderTest, NumInds){
  EXPECT_EQ(10, int(holder1.num_inds()));
  EXPECT_EQ(10, int(holder2.num_inds()));
}

TEST_F(DataholderTest, GenoName){
  holder2.add_geno_name("GENO1");
  holder2.add_geno_name("GENO2");
  EXPECT_EQ("GENO1", holder2.get_geno_name(0));
  EXPECT_EQ("GENO2", holder2.get_geno_name(1));
}

// check that default snp creation is correct
TEST_F(DataholderTest, DefaultSNPS){
  holder2.add_default_snps();
  EXPECT_EQ("1", holder2.get_geno_name(0));
  EXPECT_EQ("2", holder2.get_geno_name(1));
  
  // check snp indexes returned
  EXPECT_EQ(1, int(holder2.get_geno_index("2")));
}

// check default covariates
TEST_F(DataholderTest, DefaultCovars){
  holder2.add_default_covars();
  EXPECT_EQ("1", holder2.get_covar_name(0));
}

// check getting ind by ID
TEST_F(DataholderTest, IndividualByID){
  EXPECT_EQ(0, holder2.get_ind_by_id("tmp883")->get_genotype(0));
  EXPECT_EQ(1, holder2.get_ind_by_id("tmp8291")->get_genotype(0));
}


