//BestModelSelector.cpp
#include "BestModelSelector.h"
#include<cmath>

using namespace stat;


BestModelSelector::BestModelSelector(){
	correlationThreshold=0.8;
	crossValThreshold=2;
}


///
/// Counts occurrences of variables in each model
/// @param vars
/// @param varTotals
/// @return set containing all variables found in model
///
set<int> BestModelSelector::countVariables(vector<int>& vars,  map<int, int>& varTotals){
		set<int> multVars;
		
		vector<int>::iterator varIter;
		map<int, int>::iterator mapIter;
		
		for(varIter = vars.begin(); varIter != vars.end(); ++varIter){
			mapIter = varTotals.find(*varIter);
			if(mapIter == varTotals.end()){
				varTotals[*varIter]=0;
			}
			// increment total if variable hasn't already appeared in this model
			if(multVars.find(*varIter)==multVars.end()){
				varTotals[*varIter] += 1;
			}
			multVars.insert(*varIter);
		}	
		return multVars;
}


///
/// Increment variable cross-validation scores where needed due to correlated
/// signals
///
void BestModelSelector::addSignalCV(set<int>& varMult,
	map<int, int>& varTotals, map<int, vector<int> >& correlationSignals){

	vector<int>::iterator varIter, linkedIter; 
	map<int, int>::iterator mapIter;
	map<int, vector<int> >::iterator signalIter;

	// if they appear in the set, they have already been counted
	set<int> counted=varMult;

	for(set<int>::iterator varIter=varMult.begin(); varIter != varMult.end(); varIter++){
		signalIter = correlationSignals.find(*varIter);
		// if no correlated signal continue to next one
		if(signalIter != correlationSignals.end()){
			for(linkedIter = signalIter->second.begin(); linkedIter != signalIter->second.end();
				++linkedIter){
				// only increment the linked variables when they haven't appeared in the model
				// otherwise they have already been counted and then add it to set as counted	
				if(counted.find(*linkedIter)==counted.end()){
					varTotals[*linkedIter] += 1;
					counted.insert(*linkedIter);
				}
			}			
		}
	}


}


void BestModelSelector::selectBestVariables(std::vector<Solution*> & solutions, Dataholder * data){

	includeGenos.clear();
	includeContins.clear();
	
	vector<int> genos, contins;
	
	// CV totals
	map<int, int> genoCVTotals, continCVTotals;
	
	vector<int>::iterator varIter, linkedIter;
	map<int, int>::iterator mapIter;
	set<int> uniqueContins, uniqueGenos;
	
	// each set contains the variables for one model (removes duplicates)
	vector<set<int> > modelContinMult, modelGenoMult;
	
	for(vector<Solution*>::iterator iter=solutions.begin(); iter != solutions.end(); ++iter){
		// don't worry about dummy encoding as will just use the variables as they exist
		// they will be converted at end for display
		genos = (*iter)->getGenotypes(false);
		contins = (*iter)->getCovariates();
		
		for_each (genos.begin(), genos.end(), data_manage::Utility::decrement);
		for_each (contins.begin(), contins.end(), data_manage::Utility::decrement);
		
		modelGenoMult.push_back(countVariables(genos, genoCVTotals));
		modelContinMult.push_back(countVariables(contins, continCVTotals));
		uniqueContins.insert(modelContinMult.back().begin(), modelContinMult.back().end());
		uniqueGenos.insert(modelGenoMult.back().begin(), modelGenoMult.back().end());
	}
	genos.clear();
	copy(uniqueGenos.begin(), uniqueGenos.end(), back_inserter(genos));
	contins.clear();
	copy(uniqueContins.begin(), uniqueContins.end(), back_inserter(contins));
	
	// when continuous variables are present calculate the correlation coefficients
	// and then group those surpassing threshold.  count the CV for these signals
	// in the list of models
	genoCorrelations.clear();
	continCorrelations.clear();
	map<int, vector<int> > continCorrelationSignals = calculateCorrelation(contins, data, true);
	map<int, vector<int> > genoCorrelationSignals = calculateCorrelation(genos, data, false);
	map<int, vector<int> >::iterator signalIter;
	int setCount=0;
	
	for(vector<Solution*>::iterator iter=solutions.begin(); iter != solutions.end(); ++iter){	
		set<int>& multContin = modelContinMult[setCount];
// 		genos = (*iter)->getGenotypes();
// 		contins = (*iter)->getCovariates();
		
		addSignalCV(multContin, continCVTotals, continCorrelationSignals);
		set<int>& multGeno = modelGenoMult[setCount];
		addSignalCV(multGeno, genoCVTotals, genoCorrelationSignals);
		setCount++;
	}

	
	// final step is to create vectors containing all genotypes and continuous variables
	// that exceed the cross validation threshold 
	map<int, int>::iterator totalIter;
	for(totalIter = genoCVTotals.begin(); totalIter != genoCVTotals.end(); totalIter++){
		if(totalIter->second >= crossValThreshold){
			includeGenos.push_back(totalIter->first);
		}
	}
	for(totalIter = continCVTotals.begin(); totalIter != continCVTotals.end(); totalIter++){
		if(totalIter->second >= crossValThreshold){
			includeContins.push_back(totalIter->first);
		}
	}
	
}


///
/// Calculate coefficients
///
map<int, vector<int> > BestModelSelector::calculateCorrelation(vector<int>& vars, 
	Dataholder* data, bool continUsed){
	
	map<int, vector<int> > signals;
	
	CorrelationScore corr;
	for(vector<int>::iterator iter=vars.begin(); iter != vars.end(); ++iter){
		vector<PairedValues> pairs;
		corr.var1 = *iter;
		for(vector<int>::iterator innerIter=iter+1; innerIter != vars.end(); ++innerIter){
			if(continUsed)
				fillContinPairs(pairs, *data, *iter, *innerIter);
			else
				fillGenoPairs(pairs, *data, *iter, *innerIter);
			corr.var2 = *innerIter;
			corr.score = Spearman::calculate(pairs);
			if(abs(corr.score) >= correlationThreshold){
				signals[*iter].push_back(*innerIter);
				signals[*innerIter].push_back(*iter);
			}
			if(continUsed){
				corr.name1 = data->getCovarName(corr.var1);
				corr.name2 = data->getCovarName(corr.var2);
				continCorrelations.push_back(corr);
			}
			else{
				corr.name1 = data->getGenoName(corr.var1);
				corr.name2 = data->getGenoName(corr.var2);
				genoCorrelations.push_back(corr);
			}
		}
	}
	return signals;
}


void BestModelSelector::fillGenoPairs(vector<PairedValues>& pairs, Dataholder& data, int var1,
	int var2){
	pairs.clear();
	unsigned int nInds = data.numInds();
	float missingVal = data.getMissingGenotype();
	PairedValues pair;
	for(unsigned int i=0; i < nInds; i++){
		Individual * ind = data[i];
		if(ind->getGenotype(var1) != missingVal && ind->getGenotype(var2) != missingVal){
			pair.setPair(ind->getGenotype(var1), ind->getGenotype(var2));
			pairs.push_back(pair);
		}
	}	
}


void BestModelSelector::fillContinPairs(vector<PairedValues>& pairs, Dataholder& data, int var1,
	int var2){

	pairs.clear();
	unsigned int nInds = data.numInds();
	float missingVal = data.getMissingCoValue();
	PairedValues pair;
	for(unsigned int i=0; i < nInds; i++){
		Individual * ind = data[i];
		if(ind->getCovariate(var1) != missingVal && ind->getCovariate(var2) != missingVal){
			pair.setPair(ind->getCovariate(var1), ind->getCovariate(var2));
			pairs.push_back(pair);
		}
	}
	
}


///
/// Output operator
///
ostream& operator<< (ostream& os, BestModelSelector & selector){
	os << "Correlation threshold=" << selector.getCorrThreshold() << endl;
	std::vector<BestModelSelector::CorrelationScore>::iterator iter;
	for(iter = selector.getContinCorrelations().begin(); iter != selector.getContinCorrelations().end(); ++iter){
		os << iter->name1 << " " << iter->name2 << ":\t" << iter->score << endl;
	}
	for(iter = selector.getGenoCorrelations().begin(); iter != selector.getGenoCorrelations().end(); ++iter){
		os << iter->name1 << " " << iter->name2 << ":\t" << iter->score << endl;
	}
	return os;
}

