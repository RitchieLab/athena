//StringmanipTest.cpp

#include "../Stringmanip.h"
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

using namespace data_manage;
using namespace std;

// test conversion of the numbers to strings
TEST(StringmanipTest, itos){

  int sevenint = 7;
  EXPECT_EQ("7", Stringmanip::itos(sevenint)) << "Integer of 7 in itos";
  sevenint = -7;
  EXPECT_EQ("-7", Stringmanip::itos(sevenint)) << "Integer of -7 in itos";
  long sevenlong = 7;
  EXPECT_EQ("7", Stringmanip::itos(sevenlong)) << "Long of 7 in itos";
  unsigned int sevenuint = 7;
  EXPECT_EQ("7", Stringmanip::itos(sevenuint)) << "Unsigned int of 7 in itos";
  double sevendouble = 7;
  EXPECT_EQ("7", Stringmanip::itos(sevendouble)) << "Double of 7 in itos";
  sevendouble = 7.7;
  EXPECT_EQ("7.7", Stringmanip::itos(sevendouble)) << "Double of 7.7 in itos";
  float sevenfloat = 7.083;
  EXPECT_EQ("7.083", Stringmanip::itos(sevenfloat)) << "Float of 7.083 in itos";
}

TEST(StringmanipTest, CheckTrueFalse){
  EXPECT_EQ(true, Stringmanip::check_true_false("TRUE")) << "Checking TRUE in check_true_false";
  EXPECT_EQ(false, Stringmanip::check_true_false("FALSE")) << "Checking FALSE in check_true_false";
  EXPECT_EQ(true, Stringmanip::check_true_false("ON")) <<
    "Checking ON for true value in check_true_false";
  EXPECT_EQ(false, Stringmanip::check_true_false("AKIELAO")) <<
    "Checking nonsense string in check_true_false";
}

TEST(StringmanipTest, StringToDouble){
  double num = 1.2839;
  EXPECT_EQ(num, Stringmanip::stodouble("1.2839")) << "stodouble called for '1.2839'";
  num = -1.35;
  EXPECT_EQ(num, Stringmanip::stodouble("-1.35")) << "stodouble called for '-1.35'";
}

TEST(StringmanipTest, StringToUint){
  unsigned int num = 2938;
  EXPECT_EQ(num, Stringmanip::stouint("2938")) << "stouint called for 2938";
}

TEST(StringmanipTest, ToUpper){
  string teststr = "all lower";
  EXPECT_EQ("ALL LOWER", Stringmanip::to_upper(teststr)) << "All lower case converted in to_upper";
  teststr = "mIxED caSE";
  EXPECT_EQ("MIXED CASE", Stringmanip::to_upper(teststr)) << "Mixed case converted in to_upper";
}

TEST(StringmanipTest, IsNumber){
  string isnum = "-938.2458";
  EXPECT_TRUE(Stringmanip::is_number(isnum)) << "Checking is_number " << isnum << endl;
  string notnum = "44a.8392";
  EXPECT_FALSE(Stringmanip::is_number(notnum)) << "Checking is_number " << notnum << endl;
}

TEST(StringmanipTest, SplitString){
  string tosplit = "This should split on spaces";
  vector<string> tokens;
  Stringmanip::split_string(tosplit, " ", " ", tokens);
  EXPECT_EQ(5, int(tokens.size())) << "Should be 5 tokens on split string for test string";
  EXPECT_EQ("should", tokens[1]);
  EXPECT_EQ("spaces", tokens[4]);
}
