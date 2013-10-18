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
// NNModelLog.h

#ifndef _NNMODELLOG_H
#define	_NNMODELLOG_H

#include "NNSolution.h"
#include <fstream>

///
/// Outputs each model in the population and lists the details kept in
/// in the solutions about the model.  The output resembles the .best files
/// but with more information.  Some of the information will be grammar depth,
/// network depth, fitness and the terminal representation of the network
///

class NNModelLog{

public:
		NNModelLog(){detailed=false;}
		
		~NNModelLog();

		void openLog(std::string filename, std::string fitnessName, bool varsOnly);
		
		void closeLog();
		
		void writeSolution(NNSolution & solution, int generation, int rank);
		
		void setDetailed(bool include){detailed=include;}
		
		bool getDetailed(){return detailed;}
		
		void addGeneration(int generation);
		
		void writeVariables(NNSolution& solution, bool ottDummyEncoded);
		
private:
		
		bool detailed;
		
		void header(ostream& os, std::string fitnessName);
		
		std::ofstream logStream;

};

#endif
