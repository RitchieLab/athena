#include "AUCcalc.h"
#include<algorithm>
#include<iostream>
using namespace std;

namespace stat
{

bool TestResultSorter(TestResult const& ltr, TestResult const& rtr){
	return ltr.score > rtr.score;
}


///
/// Returns the area under the curve for a ROC curve.
/// @param results TestResult vector
/// @return area under the curve
///
float AUCCalc::calculateAUC(std::vector<TestResult>& results){
 
 	int fp=0;
 	int tp=0;
 	float fpPrev = 0.0;
 	float tpPrev = 0.0;
 	float area = 0.0;
 	float fPrev = 1e99;
 	// sort results in descending order
 	sort(results.begin(), results.end(), TestResultSorter);
 	std::vector<TestResult>::iterator resultIter=results.begin();
 	std::vector<TestResult>::iterator endIter=results.end();
 	
 	while(resultIter != endIter){
 		if(fPrev != resultIter->score){
 			area += trapArea(fp,fpPrev,tp,tpPrev);
 			fPrev = resultIter->score;
 			fpPrev = fp;
 			tpPrev = tp;
 		}
 		if(resultIter->status==1){
 			tp++;
 		}
 		else{
 			fp++;
 		}
 		++resultIter;
 	}
 	area += trapArea(fp,fpPrev, tp,tpPrev);

 	// A / (PXN) 
 	// P is number of positive examples (equals tp)
 	// N is numer of negative examples (equals fp)
 	// scale from P x N onto the unit square
 	area = area / (tp*fp); 
 	return area;
 }


///
/// Returns the trapezoidal area
/// @param x1
/// @param x2
/// @param y1
/// @param y2
/// @return area under the curve
///
float AUCCalc::trapArea(int x1, int x2, int y1, int y2){
	int base = x1-x2;
	float avgHeight = float(y1+y2)/2;
	return base * avgHeight;
}


}
