#include "MDR.h"

///
/// Calculate balanced accuracy using MDR
/// @param dSet Dataset
/// @param refSet Dataset for calculating conditional probability tables
/// @param variables
///
float MDR::calcBalAccuracy(Dataset* dSet, Dataset* refSet, vector<IndividualTerm*> &variables){

// cout << "IN MDR: ";
// for(size_t i=0; i<variables.size(); i++){
// cout << variables[i]->getLabel() << " ";
// }
// cout << endl;
		
	vector<vector<int> > totals, conditTotals;
	totalVariables(dSet, refSet, variables, totals, conditTotals);
	
	int fp=0, fn=0, tp=0, tn=0;
	int nLevels = totals.size();
	size_t nVariableCombinations = totals[0].size();
	
	for(size_t i=0; i<nVariableCombinations; i++){
		if(conditTotals[0][i] > conditTotals[1][i]){
			tn += totals[0][i];
			fn += totals[1][i];
		}
		else{
			tp += totals[1][i];
			fp += totals[0][i];
		}
	}
	
	float balacc = 0.5 * (float(tp)/(tp+fn) + float(tn)/(tn+fp));
// cout << "MDR::balacc=" << balacc << endl;
	return balacc;
}


///
/// Calculates totals for each combination of parent values 
/// @param dSet Dataset 
/// @param variables Variables for the phenotype
/// @param totals Out parameter for totals for the set for each combination
/// @param conditTotals Out parameter for conditional table totals
///
void MDR::totalVariables(Dataset* dSet, Dataset* refSet, vector<IndividualTerm*> &variables,
	vector<vector<int> >& totals, vector<vector<int> >& conditTotals){
	
	// set number of levels for each parent
	vector<int> nLevels(variables.size(), 0);
	int nl = 1;
	for(size_t i=0; i<variables.size(); i++){
		nLevels[i] = variables[i]->getNumLevels(dSet);
		nl *= nLevels[i];
	}
	
	vector<int> cumulativeLevels(variables.size(), 1);
	for(unsigned int i=1; i<cumulativeLevels.size(); i++){
		cumulativeLevels[i] = cumulativeLevels[i-1] * nLevels[i-1];
	}		
	
	vector<int> inner(nl, 0);
	totals.assign(dSet->getNumStatusLevels(), inner);
	
	deque<float> args;
	unsigned int nVariables = variables.size();
	Individual* ind;
	
	for(unsigned int i=0; i < dSet->numInds(); i++){
		ind = (*dSet)[i];
		IndividualTerm::setInd(ind);

		int value = 0;	
		for(unsigned int j=0; j < nVariables; j++){
			value += variables[j]->evaluate(args) * cumulativeLevels[j];
		}
		totals[ind->getStatus()][value]++;
	}
	
	if(refSet == dSet){
		conditTotals = totals;
	}
	else{
		conditTotals.assign(dSet->getNumStatusLevels(), inner);
		// sets are different so conditional table will be based on reference set
		for(unsigned int i=0; i < refSet->numInds(); i++){
			ind = (*refSet)[i];
			IndividualTerm::setInd(ind);

			int value = 0;	
			for(unsigned int j=0; j < nVariables; j++){
				value += variables[j]->evaluate(args) * cumulativeLevels[j];
			}
			conditTotals[ind->getStatus()][value]++;
		}
	}

}
