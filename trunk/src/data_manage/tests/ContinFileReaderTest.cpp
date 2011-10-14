//ContinFileReaderTest.cpp

#include <gtest/gtest.h>
#include <vector>
#include "../ContinFileReader.h"
#include "../Stringmanip.h"
#include "../Dataholder.h"
using namespace data_manage;
using namespace std;

class ContinFileReaderTest:public testing::Test{
  protected:
  virtual void SetUp(){
  
    string ids[10] = {"tmp883", "tmp8291", "tmp1", "ind77", "ind701", "ind1010",
      "ind7577", "tmp5", "tmp99", "ind200"};
  
    // fill holder with individuals
    for(int i=0; i<10; i++){
      Individual ind;
      ind.set_status(i*0.1);
      ind.set_id(Stringmanip::itos(i+1));
      holder.add_ind(ind);
      
      Individual indID;
      indID.set_status(i*0.1);
      indID.set_id(ids[i]);
      holderID.add_ind(indID);
    }
  }
  
  Dataholder holder, holderID;
};


// test basic reader
TEST_F(ContinFileReaderTest, ReadContinFile){

  ContinFileReader cfilereader;
  
  cfilereader.read_contin_file("sample_data/contin_sample.txt", &holder, -9999, 
    false);

  // check on some of the covariates that have been read in
  EXPECT_GT(0.0000001, fabs(holder.get_ind(0)->get_covariate(0)-0.3));
  EXPECT_GT(0.0000001, fabs(holder.get_ind(1)->get_covariate(1)-(-14.5)));
  EXPECT_GT(0.0000001, fabs(holder.get_ind(9)->get_covariate(1)-9));
  EXPECT_GT(0.0000001, fabs(holder.get_ind(7)->get_covariate(0)-0.2));
  
  // check using ID strings
  EXPECT_GT(0.0000001, fabs(holder.get_ind_by_id("1")->get_covariate(0)-0.3));
  EXPECT_GT(0.0000001, fabs(holder.get_ind_by_id("5")->get_covariate(1)-1.2));
  EXPECT_GT(0.0000001, fabs(holder.get_ind_by_id("10")->get_covariate(1)-9));

}


// read file with id
TEST_F(ContinFileReaderTest, ReadContinFileID){
  
  ContinFileReader cfilereader;
  cfilereader.read_contin_file("sample_data/contin_sample_id.txt", &holderID, -9999,
    true);
  
  // check covariates by ID
  EXPECT_GT(0.0000001, fabs(holderID.get_ind_by_id("tmp883")->get_covariate(0)-0.3));
  EXPECT_GT(0.0000001, fabs(holderID.get_ind_by_id("tmp8291")->get_covariate(1)-(-14.5)));
  EXPECT_GT(0.0000001, fabs(holderID.get_ind_by_id("ind200")->get_covariate(1)-9));
  EXPECT_GT(0.0000001, fabs(holderID.get_ind_by_id("tmp5")->get_covariate(0)-0.2));  

  // check by position
  EXPECT_GT(0.0000001, fabs(holderID.get_ind(0)->get_covariate(0)-0.3));
  EXPECT_GT(0.0000001, fabs(holderID.get_ind(1)->get_covariate(1)-(-14.5)));
  EXPECT_GT(0.0000001, fabs(holderID.get_ind(9)->get_covariate(1)-9));
  EXPECT_GT(0.0000001, fabs(holderID.get_ind(7)->get_covariate(0)-0.2));  

}
