//FileHandlerTest.cpp

// Tests all filehandler implementations.  First will only test
// the MDRFileHandler but additional formats will be supported in 
// the future.

#include "../MDRFileHandler.h"

#include <gtest/gtest.h>
using namespace data_manage;

template <class T> FileHandler* CreateFileHandler();

template <> FileHandler* CreateFileHandler<MDRFileHandler>(){
  return new MDRFileHandler;
}

// Define test fixture class template
template <class T> class FileHandlerTest: public testing::Test{

  protected:
    // constructor calls the factory function to create a 
    // FileHandler implemented by T
    FileHandlerTest():handler_(CreateFileHandler<T>()){}
    
    virtual ~FileHandlerTest(){ delete handler_;}
    
    FileHandler* const handler_;
};

#if GTEST_HAS_TYPED_TEST

using testing::Types;

typedef Types<MDRFileHandler> Implementations;

TYPED_TEST_CASE(FileHandlerTest, Implementations);

// test the standard read on file with all data present
TYPED_TEST(FileHandlerTest, TextNoMissing){
  Dataholder holder;
  
  // basic parsing
  this->handler_->parse_file("sample_data/small_sample.txt", &holder,
    -1, -9999, false);
  
  // check what is held in holder
  EXPECT_EQ(100, int(holder.num_inds())) << "Incorrect number of individuals";
  // genotypes will be empty because no default snps set yet
  EXPECT_EQ(0, int(holder.num_genotypes())) << "Should be 0 genotypes reported";
  EXPECT_FALSE(holder.any_missing_genos()) << "Should be no missing data";
  EXPECT_EQ(63, int(holder.num_genos())) << "Should have read 63 SNPs from file";
  EXPECT_EQ(2, int(holder.get_max_locus_value())) << "Max locus should be 2";
  EXPECT_EQ(3, holder.get_missing_genotype()) << "Missing genotype should now be 3";
  // check some genotypes from file
  EXPECT_EQ(0, int(holder.get_genotype(0, 0))) << "First genotype should be 0";
  // check 48th genotype
  EXPECT_EQ(2, int(holder.get_genotype(0, 47))) << "48th genotype (first ind) should be 2";
  EXPECT_EQ(2, int(holder.get_genotype(0, 62))) << "Last genotype (first ind) should be 2";
  EXPECT_EQ(1, int(holder.get_genotype(99, 0))) << "First genotype (last ind) should be 1";
  EXPECT_EQ(1, int(holder.get_genotype(99, 47))) << "48th genotype (last ind) should be 1";
  EXPECT_EQ(1, int(holder.get_genotype(99, 62))) << "Last genotype (last ind) should be 1";  
  // check the 54th individual
  EXPECT_EQ(1, int(holder.get_genotype(53, 0))) << "First genotype (54th ind) should be 1";
  EXPECT_EQ(2, int(holder.get_genotype(53, 47))) << "48th genotype (54th ind) should be 2";
  EXPECT_EQ(2, int(holder.get_genotype(53, 62))) << "Last genotype (54th ind) should be 2";
  
}

// test the standard read on file with all data present
TYPED_TEST(FileHandlerTest, TextMissing){
  Dataholder holder;
  
  // basic parsing
  this->handler_->parse_file("sample_data/small_sample_miss.txt", &holder,
    -1, -9999, false);
  
  // check what is held in holder
  EXPECT_EQ(100, int(holder.num_inds())) << "Incorrect number of individuals";
  // genotypes will be empty because no default snps set yet
  EXPECT_EQ(0, int(holder.num_genotypes())) << "Should be 0 genotypes reported";
  EXPECT_TRUE(holder.any_missing_genos()) << "Should be no missing data";
  EXPECT_EQ(63, int(holder.num_genos())) << "Should have read 63 SNPs from file";
  EXPECT_EQ(2, int(holder.get_max_locus_value())) << "Max locus should be 2";
  EXPECT_EQ(3, holder.get_missing_genotype()) << "Missing genotype should now be 3";
  
  // check some genotypes from file
  EXPECT_EQ(0, int(holder.get_genotype(0, 0))) << "First genotype should be 0";
  // check 9th genotype
  EXPECT_EQ(3, int(holder.get_genotype(0, 48))) << "9th genotype (first ind) should be 3";
  EXPECT_EQ(2, int(holder.get_genotype(0, 62))) << "Last genotype (first ind) should be 2";
  // check last individual
  EXPECT_EQ(1, int(holder.get_genotype(99, 0))) << "First genotype (last ind) should be 1";
  EXPECT_EQ(1, int(holder.get_genotype(99, 47))) << "48th genotype (last ind) should be 1";
  EXPECT_EQ(1, int(holder.get_genotype(99, 62))) << "Last genotype (last ind) should be 1";  
  // check the 54th individual
  EXPECT_EQ(1, int(holder.get_genotype(53, 0))) << "First genotype (54th ind) should be 1";
  EXPECT_EQ(2, int(holder.get_genotype(53, 47))) << "48th genotype (54th ind) should be 2";
  EXPECT_EQ(2, int(holder.get_genotype(53, 62))) << "Last genotype (54th ind) should be 2";  
}


#endif // for #if GTEST_HAS_TYPED_TEST

