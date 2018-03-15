//CVSetTest.cpp

#include "../CVSet.h"
#include "../Stringmanip.h"
#include <gtest/gtest.h>
#include <iostream>
using namespace std;
using namespace data_manage;

class CVSetTest:public testing::Test{
  protected:
  virtual void SetUp(){
    // create interval for use in test
    cvint1.num_sets(3);
    cvint2.num_sets(2);
    cvint3.num_sets(1);
    cvset.add_interval(cvint1);
    cvset.add_interval(cvint2);
    cvset.add_interval(cvint3);
  }
  CVSet cvset;
  CVInterval cvint1, cvint2, cvint3;

};


TEST_F(CVSetTest, NumIntervals){
  EXPECT_EQ(3, int(cvset.num_intervals())) << "Should have 3 intervals in set";
}

TEST_F(CVSetTest, AddInteval){
  CVSet set;
  CVInterval newint(3);
  EXPECT_EQ(0, int(set.num_intervals())) << "Set should be empty";
  set.add_interval(newint);
  EXPECT_EQ(1, int(set.num_intervals())) << "Set should have 1 interval";
  set.add_interval(cvint1);
  EXPECT_EQ(2, int(set.num_intervals())) << "Set should have 2 intervals";
}


TEST_F(CVSetTest, GetInterval){
  CVInterval& inter = cvset.get_interval(1);
  EXPECT_EQ(2, int(inter.num_sets())) << "Interval is not correct one added";
  inter = cvset.get_interval(0);
  EXPECT_EQ(3, int(inter.num_sets())) << "Interval is not correct one added";
}


