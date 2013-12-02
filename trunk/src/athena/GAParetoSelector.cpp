#include "GAParetoSelector.h"
#include <ga/gaconfig.h>

#include <ga/GAGenome.h>
#include <ga/GASelector.h>
#include <ga/garandom.h>
#include <ga/gaerror.h>
#include "GE1DArrayGenome.h"

// #include "GEObjective.h"

void GAParetoSelector::update(){
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

	for(int i=0; i<n; ++i){
		if(!addInd(i))
			break;
		totalValidInds++;
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
	
	map<int,paretoValue> paretos;

	// first inds are the best for each complexity (from lowest to most complex)
	// this is the pareto front for the population
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
			paretoValue pareto(1, *setIter, *setIter);
			paretos[*setIter]=pareto;
			++setIter;
		}
	}
	
	// after the counter reaches the target
	// additional individuals are selected
	// based on the percent chance of selection
	set<paretoValue> paretoProbs;	
	// every complexity size gets a percentage based on 1-(this size/total selected)
	// then need to adjust for percentages
	while(selectedInds.size() < nSelect){		
		float paretoTotal=0.0;
		
		for(std::set<int>::iterator setIter=complexitySet.begin(); setIter != complexitySet.end();
			++setIter){
			paretos[*setIter].probability = 1 - (paretos[*setIter].selected/(float)selectedInds.size());
			paretoTotal += paretos[*setIter].probability;
		}
		
		float currTotal=0.0;
		for(std::set<int>::iterator setIter=complexitySet.begin(); setIter != complexitySet.end();
			++setIter){
			paretos[*setIter].probability = paretos[*setIter].probability/paretoTotal + currTotal;
			currTotal = paretos[*setIter].probability;
			paretoProbs.insert(paretos[*setIter]);
		}
	
		// select a random value (use the GA library random number generator)
		paretoValue randProb;
		randProb.probability = GARandomFloat();
		set<paretoValue>::iterator selectedPareto = paretoProbs.upper_bound(randProb);
		// take next best individual from the selected complexity
		selectedInds.push_back(inds[selectedPareto->complexity].front());
		inds[selectedPareto->complexity].pop_front();
		paretos[selectedPareto->complexity].selected += 1;
		if(inds[selectedPareto->complexity].empty()){
			complexitySet.erase(selectedPareto->complexity);
			paretoProbs.clear();
		}
	}	
	*currInd=1;
	*currSelectedInd=0;
}


///
/// Adds index to correct position in map if network is valid
/// @param i index of individual in population
/// @return false if network invalid
///
bool GAParetoSelector::addInd(int i){

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
GAParetoSelector::select() const {
	
	if(*currSelectedInd >= selectedInds.size()){
		// return new genome
 		int individual = (*currInd)++;
		pop->individual(individual, (which == SCALED ? 
			  GAPopulation::SCALED : GAPopulation::RAW)).initialize();
		GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(pop->individual(individual, (which == SCALED ? 
			  GAPopulation::SCALED : GAPopulation::RAW)));
		genome.establish();  
		return pop->individual(individual, (which == SCALED ? 
			  GAPopulation::SCALED : GAPopulation::RAW));
	}
	else{
		int individual = selectedInds[*currSelectedInd];
  	(*currSelectedInd)++;
		return pop->individual(individual, (which == SCALED ? 
			  GAPopulation::SCALED : GAPopulation::RAW));
	}
}

