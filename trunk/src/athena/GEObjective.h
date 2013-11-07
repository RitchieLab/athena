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
 * File:   GEObjective.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 3:50 PM
 */

#ifndef _GEOBJECTIVE_H
#define	_GEOBJECTIVE_H

#include <GE/ge.h>
#include <ga/ga.h>
#include <Dataset.h>
#include "SolutionFactory.h"
#include "GE1DArrayGenome.h"
#include "Population.h"
#include "AthenaGrammarSI.h"

///
/// Provides objective function for use with GALIb and GE library
///
class GEObjective{
		
public:
		
		/// Used to assign * fitness values in GE
		static float GEObjectiveFunc(GAGenome& g);
		
		/// Used to output individual evaluations
		static float GEObjectiveFuncOut(GAGenome& g, ostream& os);
		
		/// Outputs symbols from genome to output
		static void outputSymbols(GAGenome& g, ostream& os);
		
		/// sets the mapper to use
		static void setMapper(AthenaGrammarSI* m){
			mapper = m;
			mapper->setLeftOptBound(solCreator->getLeftOptBound());
			mapper->setRightOptBound(solCreator->getRightOptBound());
			mapper->setArgSymbols(solCreator->getOptArgSymbols());
		}
		
		/// sets the Dataset for objective function to work with
		static void setDataset(data_manage::Dataset* ds);
		
		/// sets Solution type for objective function
		static void setSolutionType(std::string solutionName, std::string calculatorName){
				solCreator = SolutionFactory::createSolution(solutionName);
				solCreator->setCalculator(calculatorName);
		}
		
		/// Alternative method for creating Solution
		static void setSolutionType(std::string solutionName, std::string calculatorName,
			vector<string>& vars){
			solCreator = SolutionFactory::createSolution(solutionName, vars);
			solCreator->setCalculator(calculatorName);
		}
		
		static std::string calculatorName(){return solCreator->calculatorName();}
		
		static Solution* getBlankSolution(){return solCreator->createNewSolution();}
		
		static bool maxBest(){return solCreator->maxBest();}
		
		static std::string getGraphicalExt(){return solCreator->graphicExt();}
		
		static void outputModel(ostream& os, Solution* sol, data_manage::Dataholder* data,
			bool mapUsed, bool ottDummy, bool continMapUsed){
			solCreator->establishSolution(sol->getSymbols());
			solCreator->graphicalOutput(os, data, mapUsed, ottDummy, continMapUsed);
			solCreator->freeSolution();
		}
		
		static void outputEquation(ostream& os, Solution* sol, data_manage::Dataholder* data,
			bool mapUsed, bool ottDummy, bool continMapUsed){
			solCreator->establishSolutionEquation(sol->getSymbols());
			solCreator->equationOutput(os, data, mapUsed, ottDummy, continMapUsed);
			solCreator->freeSolution();
		}
		
		static void GEObjectiveInit(GAGenome& g);

		/// Return additional column names for output
		static vector<std::string>  getAdditionalOutputNames(){return solCreator->getAdditionalOutputNames();}
		
		/// Return values for final output
		static vector<std::string> getAdditionalFinalOutput(GAGenome& g);
		
		
		/// Optimizes current model using process provided by SolutionCreator
		static void optimizeSolution(GAGenome& g);
		
		/// Sets maximum genome size
		static void setMaxGenomeSize(unsigned int maxSize){maxGenSize = maxSize;}
	 
		// returns worst score
		static float getWorstScore(){return solCreator->getWorst();}
		
		static void addLogging(bool val){additionalLogging=val;}

		static void setrank(int r){rank=r;}
		static int rank;
private:

		static AthenaGrammarSI* mapper;
		static data_manage::Dataset* set;
		static SolutionCreator* solCreator;
		static unsigned int maxGenSize;
		static bool additionalLogging;  
};

#endif	/* _GEOBJECTIVE_H */

