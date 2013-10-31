#include "GAParetoRankSelector.h"
#include <ga/gaconfig.h>

#include <ga/GAGenome.h>
#include <ga/GASelector.h>
#include <ga/garandom.h>
#include <ga/gaerror.h>
#include "GE1DArrayGenome.h"

void GAParetoRankSelector::update(){

	n = pop->size();
	
	inds.clear();
	complexitySet.clear();
	
	if(which == GASelectionScheme::RAW){
		pop->sort(gaFalse, GAPopulation::RAW);
		sortBasis = GAPopulation::RAW;
	}
	else{
		pop->sort(gaFalse, GAPopulation::SCALED);
		sortBasis = GAPopulation::SCALED;
	}
	
	int totalValidInds=0;
	
	nComplexities=0;
	
	// iterate through population
	// establish the initial deques containing 
	// the different complexity sizes
	// skip any solutions that are invalid
	if(pop->order() == GAPopulation::HIGH_IS_BEST){
		for(int i=n-1; i>=0; --i){
			if(!addInd(i))
				break;
			totalValidInds++;
		}
	}
	else{
		for(int i=0; i<n; ++i){
			if(!addInd(i))
				break;
			totalValidInds++;
		}
	}
	
	// determine number that will be selected from these (at most 1/2)
	// any beyond this number will be random networks 
	if(totalValidInds < n/2){
		nSelect = totalValidInds;
	}
	else{
		nSelect = n/2;
	}

	// start filling up the list of individuals to be selected
	selectedInds.clear();
	
	// first inds are the best for each complexity (from lowest to most complex)
	// this is the pareto front for the population
	while(selectedInds.size() < nSelect){
		for(std::set<int>::iterator setIter=complexitySet.begin(); setIter != complexitySet.end();){
			selectedInds.push_back(inds[*setIter].front());
			inds[*setIter].pop_front();
			if(inds[*setIter].empty()){
				std::set<int>::iterator deleteIter = setIter;
				++setIter;
			// remove from complexities when no more  individuals of that complexity
				complexitySet.erase(deleteIter);
			}
			else{
				++setIter;
			}
			if(selectedInds.size() == nSelect)
				break;
		}
	}

	*currInd=0;
	*currSelectedInd=0;
}


///
/// Adds index to correct position in map if network is valid
/// @param i index of individual in population
/// @return false if network invalid
///
bool GAParetoRankSelector::addInd(int i){

	GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(pop->individual(i, sortBasis));
	if(!genome.isValid()){
		return false;
	}	
	int complexity = genome.getComplexity();
	std::map<int, std::deque<int> >::iterator iter = inds.find(complexity);
	
	std::pair<std::map<int,deque<int> >::iterator, bool> ret;
	
	if(iter == inds.end()){
		deque<int> empty;
		ret = inds.insert(std::pair<int,deque<int> >(complexity,empty));
		iter=ret.first;		
		nComplexities++;
		complexitySet.insert(complexity);
	}
	iter->second.push_back(i);
	return true;
}


// return individuals from population until all selected
// then return sensibly initialized new individuals
GAGenome &
GAParetoRankSelector::select() const {

	if(*currSelectedInd >= selectedInds.size()){
		// return new genome
 		(*currInd)++;
// 		pop->individual(*currInd, (which == SCALED ? 
// 			  GAPopulation::SCALED : GAPopulation::RAW)).initialize();
		GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(pop->individual(*currInd, (which == SCALED ? 
			  GAPopulation::SCALED : GAPopulation::RAW)));
		genome.initialize();
		genome.establish();		  
		return pop->individual(*currInd, (which == SCALED ? 
			  GAPopulation::SCALED : GAPopulation::RAW));
	}
	else{
		int individual = selectedInds[*currSelectedInd];
  	(*currSelectedInd)++;
		return pop->individual(individual, (which == SCALED ? 
			  GAPopulation::SCALED : GAPopulation::RAW));
	}
}

