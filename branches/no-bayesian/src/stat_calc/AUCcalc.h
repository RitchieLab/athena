#ifndef AUC_CALC_H
#define AUC_CALC_H

#include <vector>

namespace stat
{


struct TestResult{
	float score;
	int status;
};

/// sort TestResult in descending order
bool TestResultSorter(TestResult const& ltr, TestResult const& rtr);


class AUCCalc
{
public:
  AUCCalc();
  ~AUCCalc();

  /// returns the AUC for the values passed
  static float calculateAUC(std::vector<TestResult>& results);
  
  /// calculate area of trapezoid
  static float trapArea(int x1, int x2, int y1, int y2);
};

}

#endif /*AUC_CALC_H*/