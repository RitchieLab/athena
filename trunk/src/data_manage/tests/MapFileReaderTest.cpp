//MapFileReaderTest.cpp

#include "../MapFileReader.h"
#include "../Stringmanip.h"
#include <gtest/gtest.h>
#include <iostream>

using namespace data_manage;
using namespace std;

class MapFileReaderTest:public testing::Test{
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
      // add 12 snps to match the sample map file
      for(int snp=0; snp<12;snp++){
        indID->add_genotype((snp+i)%3);
      }
      holder2.add_ind(indID);
    }
  }

  
  Dataholder holder1, holder2;
};


// read map file
TEST_F(MapFileReaderTest, MapFileRead){
  MapFileReader mapread;
  mapread.parse_map_file("sample_data/sample12.map", &holder2);
  
  EXPECT_EQ(2, int(holder2.get_geno_index("rs3"))) << "get_geno_index for rs3 is incorrect";
  EXPECT_EQ(10, int(holder2.get_geno_index("rs11"))) << "get_geno_index for rs11 is incorrect";
  EXPECT_EQ("rs9", holder2.get_geno_name(8)) << "get_geno_name for index 8 is incorrect";
  EXPECT_EQ("rs1", holder2.get_geno_name(0)) << "get_geno_name for index 0 is incorrect";
  EXPECT_EQ("rs12", holder2.get_geno_name(11)) << "get_geno_name for index 11 is incorrect";
  
}
