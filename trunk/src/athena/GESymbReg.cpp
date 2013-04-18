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
/* 
 * File:   GESymbReg.cpp
 * Author: dudeksm
 *
 * Created on August 10, 2010, 2:03 PM
 */
#include "GESymbReg.h"
#include "GEObjective.h"
#include "Terminals.h"
#include "GENNGrammarAdjuster.h"
#include <ga/ga.h>
#include <ctime>
 
 ///
/// Sets parameters within the algorithm for the 
/// analysis.
/// @param algParam AlgorithmParams
/// @param numExchanges total number of times the best model will be passed to other algorithms
/// when more than one algorithm
/// @param numGenos number of Genotypes in set
/// @param numContin number of Continuous variables in set
/// 
void GESymbReg::setParams(AlgorithmParams& algParam, int numExchanges, int numGenos, int numContin){
			
		map<string, string>::iterator mapIter;
		
		for(mapIter = algParam.params.begin(); mapIter != algParam.params.end(); 
			mapIter++){       
				switch(paramMap[mapIter->first]){
						case noMatchParam:
								throw AthenaExcept("No match for parameter " + mapIter->first +
												"in Algorithm GE Symbolic Regression");
								break;
						case minSizeParam:
								minSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case maxSizeParam:
								maxSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case tailRatioParam:
								tailRatio = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case growRateParam:
								growRate = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case maxDepthParam:
								maxDepth = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case tailSizeParam:
								tailSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case sensibleInitParam:
								sensibleInit = Stringmanip::check_true_false(mapIter->second);
								break;
						case popSizeParam:
								popSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case probCrossParam:
								probCross = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case probMutParam:
								probMut = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case gramFileParam:
								grammarFile = mapIter->second;
								break;
						case stepSizeParam:
								stepSize = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case calcType:
								calculatorName = Stringmanip::to_upper(mapIter->second);
								break;        
						case useEffectiveXO:
								effectiveXO = Stringmanip::check_true_false(mapIter->second);
								break;
						case useAllSnps:
								useAllVars = Stringmanip::check_true_false(mapIter->second);
								break;
						case useAllCovariates:
								useAllCovars = Stringmanip::check_true_false(mapIter->second);
								break;
						case requireAll:
								requireAllVars = Stringmanip::check_true_false(mapIter->second);
								break;
						case requireAllOnce:
								requireAllVarsOnce = Stringmanip::check_true_false(mapIter->second);
								break;
						case bioInitFract:
								initBioFract = Stringmanip::stringToNumber<double>(mapIter->second);     
								break;
						case restrictVarGens:
								ngensVarRestrict = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case bioModelSelection:
								if(bioModelSelectionMap.find(Stringmanip::to_upper(mapIter->second)) == 
										bioModelSelectionMap.end())
									throw AthenaExcept("No match for bio model selection type " + mapIter->second);
								else
									biofilterSelectorType = bioModelSelectionMap[Stringmanip::to_upper(mapIter->second)];
								break;
						case gaSelection:
								if(gaSelectorMap.find(Stringmanip::to_upper(mapIter->second)) ==
									gaSelectorMap.end())
									throw AthenaExcept("No match for GA selection type " + mapIter->second);
								else
									gaSelector = gaSelectorMap[Stringmanip::to_upper(mapIter->second)];
									break;
						case doubleTournF:
								doubleTourneyF = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case doubleTournD:
								doubleTourneyD = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
						case doubleTournFitFirst:
								fitFirst = Stringmanip::check_true_false(mapIter->second);
								break;
						case blockCrossGens:
								ngensBlockCross = Stringmanip::stringToNumber<unsigned int>(mapIter->second);
								break;
						case resetVarsAtMigration:
								resetRestrictedAtMigration = Stringmanip::check_true_false(mapIter->second);
								break;
#ifdef ATHENA_BLOAT_CONTROL
						case prunePlantFract:
								pruneAndPlantFract = Stringmanip::stringToNumber<double>(mapIter->second);
								break;
#endif
						case bpstart:
								throw AthenaExcept("No match for parameter " + mapIter->first +
												" in Algorithm GE Symbolic Regression");
								break;
						case bpfreq:    
								throw AthenaExcept("No match for parameter " + mapIter->first +
												" in Algorithm GE Symbolic Regression");
								break;
						default:
								throw AthenaExcept("No match for parameter " + mapIter->first +
												" in Algorithm GE Symbolic Regression");               
				}
		}
		
		if(stepSize > numGenerations){
				stepSize = numGenerations;
		}
		
		numGenerations = stepSize * numExchanges;
		numGenotypes = numGenos;
		numContinuous = numContin;

		setGAParams();   
}



///
/// Sets parameters for use with GAlib
/// @throws AthenaExcept on error
///
void GESymbReg::setGAParams(){   
		GARandomSeed(randSeed);
		srand(randSeed);
	
		// first free ga memory if already run
		freeMemory();
		
		//Initialize the GE mapper
	 //Set maximum number of wrapping events per mapping
	 mapper.setMaxWraps(wrapEvents);
		adjuster.readGrammarFile(grammarFile);
		if(useAllVars){
			adjuster.includeAllVars(numGenotypes, numContinuous);
		}
		else{
			adjuster.expandVariables();
			if(dummyEncoded)
				adjuster.doubleGenotypeGrammar();
		 }
		adjuster.setMapper(mapper);

		setMapperPrefs(mapper);
	  
	  if(!requireAllVars && !requireAllVarsOnce)
			GEObjective::setSolutionType("SYMBREG", calculatorName);
		else if(requireAllVars){
			vector<string> varStrings = adjuster.getVariables();
			GEObjective::setSolutionType("SYMBREGALL", calculatorName, varStrings);
		}
		else if(requireAllVarsOnce){
			vector<string> varStrings = adjuster.getVariables();
			GEObjective::setSolutionType("SYMBREGONCE", calculatorName, varStrings);
		}

	 if(calculatorName.compare("RSQUARED")==0){
		 pop.setConvertScores(true);
	 }

	 // add mapper to objective function
	 GEObjective::setMapper(&mapper);
	 InitGEgenome::setMinSize(minSize);
	 InitGEgenome::setMaxSize(maxSize);
	 InitGEgenome::setMapper(&mapper);
}



///
/// Runs an indicated interval for the algorithm. This interval is set
/// in the stepSize variable.  
///
int GESymbReg::step(){

		GE1DArrayGenome::myRank = myRank;

		for(unsigned int i=0; i < stepSize; i++){
				if(!ga->done()){
				 if(logTypeSelected!=LogNone){
						geLog->addGeneration();
				 }         
				 ga->step();
				 restrictStepsDone++;
 
					if(ngensVarRestrict && restrictStepsDone == ngensVarRestrict){
					 // have to convert the networks back to original mapping
					 GEObjective::setMapper(&mapper);
					 GE1DArrayGenome::setMapper(&mapper);
					 convertNetworks(restrictMapper, mapper);
				 }
				 
				 // check for need to change crossovers
				 if(ngensBlockCross && restrictStepsDone == ngensBlockCross){
						resetCrossover();
				 }
 
				 fillLog();
				}
		}

		#ifdef PARALLEL
			// if restricted variables has been used and are still in effect
			// need to convert all networks back to original grammar and exchange
			// then construct new grammar restricted to only variables in the 
			// new population and continue to process
			if(ngensVarRestrict && restrictStepsDone < ngensVarRestrict){
			// have to convert the networks back to original mapping
				 GEObjective::setMapper(&mapper);
				 GE1DArrayGenome::setMapper(&mapper);
				 convertNetworks(restrictMapper, mapper);
			}
			SendAndReceiveStruct(totalNodes, myRank);

			if(ngensVarRestrict && restrictStepsDone < ngensVarRestrict){
				// after transfer construct new grammar and set for use
				setRestrictedGrammar(resetRestrictedAtMigration);
			}

		#endif
		
		// only need to fill population at this point not at end of each generation
		fillPopulation();
		
		return 0;
}



