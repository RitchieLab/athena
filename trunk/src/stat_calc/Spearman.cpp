#include "Spearman.h"
#include<algorithm>
#include<cmath>
using namespace std;

namespace stat
{

///
/// Returns the coefficient of the ranks of the PairedValues
/// Note that the vectors passed will be sorted 
/// @param pairs PairedValues
/// @return coefficient
///
float Spearman::calculate(vector<PairedValues>& pairs){
 
	// rank the two vectors by value and assign rank values
	rank(pairs,0);
	rank(pairs,1);
	
	// calculate coefficient
	// calculation adapted from ftp://ftp.biostat.wisc.edu/pub/chappell/800/hw/spearman.pdf
	float xiSum=0.0, yiSum=0.0, xi, yi;
	
	for(vector<PairedValues>::iterator iter=pairs.begin(); iter != pairs.end(); iter++){
		xi = iter->ranks[0];
		yi = iter->ranks[1];	
		
		xiSum += xi;
		yiSum += yi;
	}
	
	unsigned int size = pairs.size();
	float xMean=xiSum/size;
	float yMean=yiSum/size;
	float numerator=0.0;
	float denominator=0.0;
	float xDiff,yDiff,xDiffDenom=0.0,yDiffDenom=0.0;
	
	for(vector<PairedValues>::iterator iter=pairs.begin(); iter != pairs.end(); iter++){
		xi = iter->ranks[0];
		yi = iter->ranks[1];	
		xDiff = xi-xMean;
		yDiff = yi-yMean;
		numerator += xDiff * yDiff;
		xDiffDenom += xDiff * xDiff;
		yDiffDenom += yDiff * yDiff;
	}
	denominator = xDiffDenom * yDiffDenom;
	return numerator / sqrt(denominator);
}



void Spearman::rank(vector<PairedValues>& pairs, int index){

	sort(pairs.begin(), pairs.end(), RankValues(index));
	
	float rank=1.0;
	vector<PairedValues>::iterator startIter;
	int rankNum;
	
	for(vector<PairedValues>::iterator iter=pairs.begin(); iter != pairs.end(); iter++){
		float rankTotal = rank;
		rankNum=1;
		startIter = iter;		
		while(iter != pairs.end() && iter->values[index] == (iter+1)->values[index]){
			rank++;
			rankTotal += rank;
			iter++;
			rankNum++;
		}
		
		float rankValue = rankTotal/rankNum;
		
		for(;startIter != iter+1; startIter++){
			startIter->ranks[index] = rankValue;
		}
		rank++;
		
	}
	
}	
	
}