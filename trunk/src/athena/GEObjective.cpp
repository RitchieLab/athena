/*
Copyright Marylyn Ritchie 2011

This file is part of ATHENA.

ATHENA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ATHENA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ATHENA.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "GEObjective.h"
#include <math.h>
#include "Structs.h"

AthenaGrammarSI* GEObjective::mapper = NULL;
data_manage::Dataset* GEObjective::set = NULL;
data_manage::Dataset* GEObjective::referenceSet = NULL;
SolutionCreator* GEObjective::solCreator = NULL;
unsigned int GEObjective::maxGenSize = 250;
bool GEObjective::additionalLogging = false;

int GEObjective::rank=0;

///
/// Objective function
/// @param g GAGenome which is being evaluated
/// @return fitness
///
float GEObjective::GEObjectiveFunc(GAGenome& g){

	 GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);

	//Assign genotype to mapper
	mapper->setGenotype(genome); 
	 Phenotype const *phenotype=mapper->getPhenotype();
	 
	 float fitness;
	 genome.setValid(phenotype->getValid());
	 if(phenotype->getValid()){
	 
			 unsigned int phenoSize=(*phenotype).size();
			 vector<string> symbols(phenoSize, "");     
			 
			for(unsigned int i=0; i<phenoSize; ++i){
					symbols[i] = *((*phenotype)[i]);
// cout << symbols[i] << " ";
			}
// cout << endl;

			try{
				solCreator->establishSolution(symbols, set);
				}catch(AthenaExcept& ae){
					fitness = solCreator->getWorst();
	  			 genome.clearScores();
// cout << "failed: fitness=" << fitness << endl;
	  			 return fitness;
				}

			// alter genome to match any variable changes
			if(solCreator->anyChangedVariables()){
				// change them in the mapper
				mapper->changeVariables(genome, solCreator->getChangedVariables());		
			}

			fitness = solCreator->evaluate(set);
// cout << "success: fitness=" << fitness << endl;
			if(additionalLogging){
				solCreator->detailedLogging();
				genome.setDepth(solCreator->getDetailedLog());
				genome.setGramDepth(mapper->buildDerivationTree());
			}
			
			solCreator->freeSolution();

			genome.setEffectiveSize(mapper->getGenotype()->getEffectiveSize());   
			genome.setNumCovars(solCreator->getNumCovars());
			genome.setNumGenes(solCreator->getNumGenes());
			genome.addGenos(solCreator->getGeneIndexes());
			genome.addCovars(solCreator->getCovarIndexes());
			genome.setNumIndsEvaluated(solCreator->getNumIndsEvaluated());
			genome.setNumNodes(solCreator->getNumNodes());
			genome.setComplexity(solCreator->getComplexity());
			
			// when set 
			if(genome.getNumIndsEvaluated() != int(set->numInds())){
				genome.setSSTotal(solCreator->getCalculatorConstant());
			}
			else{
				genome.setSSTotal(set->getSSTotal());
			}
/////////////   DEBUG 
// if(solCreator->anyChangedVariables()){			
// mapper->setGenotype(genome); 
// Phenotype const *phenotype=mapper->getPhenotype();			
// for(unsigned int i=0; i<phenoSize; ++i){
// symbols[i] = *((*phenotype)[i]);
// cout << symbols[i] << " ";
// }
// cout << endl; cout << "-------------------------------------------------" << endl;
// }
///////////   END DEBUG 
	 }
	 else{
				// set fitness to worst score initially
		 fitness = solCreator->getWorst();
		 genome.clearScores();
	 }

	return fitness;
}


///
/// Performs basic tasks (except evaluation) on genome
///  Necessary for conducting certain types of crossover
/// @param g GAGenome which is being evaluated
/// @return none
///
void GEObjective::GEObjectiveInit(GAGenome& g){

	 GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);

	//Assign genotype to mapper
	mapper->setGenotype(genome); 
	 Phenotype const *phenotype=mapper->getPhenotype();
	 
	 float fitness;
	 genome.setValid(phenotype->getValid());
	
	
	 if(phenotype->getValid()){
			genome.setEffectiveSize(mapper->getGenotype()->getEffectiveSize());
	 }
	 else{
				// set fitness to worst score initially
		 fitness = solCreator->getWorst();
		 genome.clearScores();
	 }

}


///
/// Outputs symbols from genome to output
///
void GEObjective::outputSymbols(GAGenome& g, ostream& os){
	GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);

	//Assign genotype to mapper
	mapper->setGenotype(genome);
	 Phenotype const *phenotype=mapper->getPhenotype();  
			 unsigned int phenoSize=(*phenotype).size();
			 vector<string> symbols(phenoSize, "");     

			for(unsigned int i=0; i<phenoSize; ++i){
					symbols[i] = *((*phenotype)[i]);        
						os << symbols[i] << " ";
			}   
		os << endl;
}


///
/// Objective function which pass output stream to evaluator so that individual evaluations
/// can be displayed.
/// @param g GAGenome which is being evaluated
/// @param os
/// @return fitness
///
float GEObjective::GEObjectiveFuncOut(GAGenome& g, ostream& os){

	 GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);
	 
	//Assign genotype to mapper
	mapper->setGenotype(genome);

	 Phenotype const *phenotype=mapper->getPhenotype();
	 
	 float fitness;

	 if(phenotype->getValid()){
 
			 unsigned int phenoSize=(*phenotype).size();
			 vector<string> symbols(phenoSize, "");     

			for(unsigned int i=0; i<phenoSize; ++i){
					symbols[i] = *((*phenotype)[i]);
			}

			solCreator->establishSolution(symbols, set);
			fitness = solCreator->evaluateWithOutput(set, os);
			solCreator->freeSolution();
			
			genome.setEffectiveSize(mapper->getGenotype()->getEffectiveSize());   

	 }
	 else{
				// set fitness to worst score initially
		 fitness = solCreator->getWorst();
	 }
	 return fitness;
}

///
/// Calculates fitness on model supplied
/// 
///
void GEObjective::calcFitness(Solution* sol){

	solCreator->establishSolution(sol->getSymbols(), set);
	double fitness = solCreator->evaluate(set); 
	sol->fitness(fitness);

	solCreator->freeSolution();
}


///
/// Calculates fitness on model supplied and writes each individual
/// to the output stream
/// 
///
void GEObjective::calcFitnessOut(Solution* sol, ostream& os){

	solCreator->establishSolution(sol->getSymbols(), set);
	double fitness = solCreator->evaluateWithOutput(set, os);
	sol->fitness(fitness);

	solCreator->freeSolution();
}

///
/// Return values for final output
/// @param g GAGenome to analyze
///
vector<std::string> GEObjective::calcAdditionalFinalOutput(Solution* sol){ 
	try{
		solCreator->establishSolution(sol->getSymbols(), set);
	}catch(AthenaExcept& ae){
		// fails so return empty strings
		vector<string> tmp;
		return tmp;
	}
	solCreator->evaluateForOutput(set);
	solCreator->freeSolution();
		
	return solCreator->getAdditionalFinalOutput();	
}



///
/// Return values for final output
/// @param g GAGenome to analyze
///
vector<std::string> GEObjective::getAdditionalFinalOutput(GAGenome& g){

	 GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);
	 
	//Assign genotype to mapper
	mapper->setGenotype(genome);

	 Phenotype const *phenotype=mapper->getPhenotype();
	 
	 solCreator->clearAdditionalOutput();
// 	vector<std::string> outputValues;

	 if(phenotype->getValid()){
// cout << "valid" << endl;
			 unsigned int phenoSize=(*phenotype).size();
			 vector<string> symbols(phenoSize, "");     

			for(unsigned int i=0; i<phenoSize; ++i){
					symbols[i] = *((*phenotype)[i]);
// cout << symbols[i] << " ";
			}
		  try{
// cout << endl;
				solCreator->establishSolution(symbols, set);
			}catch(AthenaExcept& ae){
				// return empty string as this isn't a valid network
				vector<std::string> tmp;
				return tmp;
			}
			solCreator->evaluateForOutput(set, referenceSet);
			solCreator->freeSolution();
			
			genome.setEffectiveSize(mapper->getGenotype()->getEffectiveSize());   

	 }
// cout << "returning " << solCreator->getAdditionalFinalOutput()[0] << " " << solCreator->getAdditionalFinalOutput()[1] << endl;
	 return solCreator->getAdditionalFinalOutput();	
}


///
/// sets the Dataset for objective function to work with
///
void GEObjective::setDataset(data_manage::Dataset* ds){
	set = ds; 
	if(!set->isCaseControl() && solCreator->getCalculator()->requiresCaseControl()){
		throw AthenaExcept(solCreator->getCalculator()->getName() + " requires a case-control dataset");
	}
	solCreator->setCalculatorConstant(ds);
}

///
/// sets the Dataset for objective function to work with
///
void GEObjective::setRefDataset(data_manage::Dataset* ds){
	referenceSet = ds; 
// 	if(!set->isCaseControl() && solCreator->getCalculator()->requiresCaseControl()){
// 		throw AthenaExcept(solCreator->getCalculator()->getName() + " requires a case-control dataset");
// 	}
// 	solCreator->setCalculatorConstant(ds);
}

///
/// Optimizes current model using process provided by SolutionCreator
/// @param g GAGenome to optimize
/// @returns optimized score
///
void GEObjective::optimizeSolution(GAGenome& g){

	float oldScore, optScore;

	GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);
	//Assign genotype to mapper
	vector<AthenaGrammarSI::codonBlocks> blocks = mapper->setGenotypeOpt(genome, 
	  solCreator->getOptIncluded(), solCreator->getStartOptSymbol(),
	  solCreator->singleOpt());
	Phenotype const *phenotype=mapper->getPhenotype();	
// cout << endl;
// for(unsigned int i=0; i<blocks.size(); i++){
// cout << "i=" << i << " block.start=" << blocks[i].start << " block.end=" << blocks[i].end << endl;
// }
// unsigned int phenoSize=(*phenotype).size();
// vector<string> symbols(phenoSize, "");     			 
// for(unsigned int i=0; i<phenoSize; ++i){
// symbols[i] = *((*phenotype)[i]);
// cout << symbols[i] << " ";
// }
// cout << endl;

	// run optimization on the model if it is a valid solution
	if(phenotype->getValid()){
		unsigned int phenoSize=(*phenotype).size();
		vector<string> symbols(phenoSize, "");     

		oldScore = genome.score();
		// network is not optimizable
		if(oldScore == solCreator->getWorst()){
// cout << "worst score " << endl;
			return;
		}

		for(unsigned int i=0; i<phenoSize; ++i){
		 symbols[i] = *((*phenotype)[i]);
		}

		int numEpochsTrained = solCreator->optimizeSolution(symbols, set);
	 
		// after optimization, get the new constant list
		vector<symbVector> newWeights = solCreator->getOptimizedSymbols();
// cout << "number of newWeights=" << newWeights.size() << endl;
// cout << "genome.score=" << genome.score() << endl;

		optScore = solCreator->getOptimizedScore();
// cout << "optScore=" << optScore << endl;

		// skip the optimization if it isn't improving 
		if(optScore < genome.score()){
		// go through original list and copy to new genome 
		// use vector to create new genome list
			vector<int> newCodons;
			unsigned int currBlock = 0;    
		// replace existing operator chunks with new values and set that to be the genome
		// copy entire genome including possibly unused portion at end after effective size
			for(int i=0; i<genome.size();){
				if(currBlock < blocks.size() && i==blocks[currBlock].start){
				// set i to be equal to end of original block so that copying will continue
				// from the correct position
					i=blocks[currBlock].end;
					vector<int> tempCodons = mapper->translateOptValue(newWeights[currBlock]);
				// in mapper have function that returns codon list for a specified value
					newCodons.insert(newCodons.end(), tempCodons.begin(), tempCodons.end());
					currBlock++;

				}
				else{
					newCodons.push_back(genome.gene(i));
					i++;
				}
			}

		// resize if new codon size is able to fit within the max genome
			if(newCodons.size() <= maxGenSize){
				genome.resize(newCodons.size());
				int i=0;
				for(vector<int>::iterator iter=newCodons.begin(); iter != newCodons.end(); ++iter){
					genome.gene(i++, *iter);
				}
			}

			genome.setNumEpochsTrained(numEpochsTrained);

		// as last step evaluate new genome and set score in it to be new value
			genome.score(GEObjective::GEObjectiveFunc(genome));
// cout << "orig=" << oldScore << " k2 score=" << genome.score() << endl;
		}
	}  	
// exit(1);
}

///
/// Add constants to solution creator.
///
void GEObjective::addConstants(std::vector<std::string> constants){
	solCreator->addConstants(constants);
}



