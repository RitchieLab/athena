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
			
		GENNAlg::setParams(algParam, numExchanges, numGenos, numContin);	
			
		map<string, string>::iterator mapIter;
		
		for(mapIter = algParam.params.begin(); mapIter != algParam.params.end(); 
			mapIter++){       
				switch(paramMap[mapIter->first]){
						case noMatchParam:
								throw AthenaExcept("No match for parameter " + mapIter->first +
												"in Algorithm GE Symbolic Regression");
								break;
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
	 expandVariables();
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

	 if(calculatorName.compare("RSQUARED")==0 || calculatorName.compare("MEANABSOLUTE")==0){
//    if(calculatorName.compare("RSQUARED")==0){
		 pop.setConvertScores(true);
	 }

	 // add mapper to objective function
	 GEObjective::setMapper(&mapper);
	 setInitParams();
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
			popMigrator.sendAndReceiveStruct(totalNodes, myRank, ga);

			if(ngensVarRestrict && restrictStepsDone < ngensVarRestrict){
				// after transfer construct new grammar and set for use
				setRestrictedGrammar(resetRestrictedAtMigration);
			}

		#endif
		
		// only need to fill population at this point not at end of each generation
		fillPopulation();
		
		return 0;
}



