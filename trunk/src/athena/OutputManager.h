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
 * File:   OutputManager.h
 * Author: dudeksm
 *
 * Created on December 29, 2008, 3:37 PM
 */

#ifndef _OUTPUTMANAGER_H
#define	_OUTPUTMANAGER_H

#include "Population.h"
#include "AthenaExcept.h"
#include <vector>
#include <fstream>
#include "Algorithm.h"


///
/// Ouputs best models from the analysis as a summary
/// and as individual files showing the results.
///

class OutputManager{
		
public:
		
		/// sets basename for output
		void setBasename(std::string base){basename = base;}
		
		/// creates new files for writing
		void setFiles(bool mapFileUsed, std::string fitnessName, 
			std::vector<std::string> additionalHeaders);
		
		/// outputs summary of best models
		void outputSummary(Population& pop, int currPop, data_manage::Dataholder& data, 
			std::vector<std::string> extraColumns,
			bool mapFileUsed=false, bool dummyEncoded=true, bool continMapUsed=false,
			std::string fitnessName=" ");
		
		void outputBest(Solution* bestmodel, data_manage::Dataholder& data,
			std::vector<std::string>& extraColumns,
			bool mapFileUsed=false, bool dummyEncoded=true, bool continMapUsed=false,
			std::string fitnessName=" ");
		
		/// outputs a file for each best model
		void outputBestModels(Population& pop, int nmodels, int currPop, std::string scaleInfo,
			data_manage::Dataholder& data, bool mapUsed, bool ottDummy, bool continMapUsed);
		
		/// returns a stream for writing
		std::ostream& getStream(std::string filename);
		
		/// closes the provided stream
		void closeStream(){if(logStream.is_open()){ logStream.close();}}
		
		void outputInds(std::istream &is, std::string base);
		
		/// output graphical representation as defined in algorithm
		void outputGraphic(Algorithm* alg,  Population& pop, int currPop, std::string basename, int numModels,
			data_manage::Dataholder& data, bool mapUsed, bool ottDummy, bool continMapUsed);
		
private:
		
		std::string basename;
		std::ofstream logStream;
		std::vector<std::string> addHeaders;
	
		
		struct indOutScores{
			std::string output;
			double diff;
		};
		
		struct scoreComp {
		  bool operator() (const indOutScores& lhs, const indOutScores& rhs) const
		  {return lhs.diff>rhs.diff;}
		};
		
};

#endif	/* _OUTPUTMANAGER_H */

