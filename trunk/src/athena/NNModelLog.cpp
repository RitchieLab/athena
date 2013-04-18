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
//NNModelLog.cpp

#include "NNModelLog.h"
#include <iomanip>

///
/// Destructor
///
NNModelLog::~NNModelLog(){
		if(logStream.is_open()){
				closeLog();
		}
}

void NNModelLog::header(ostream& os, std::string fitnessName){
	os << "GEN\tRANK\t" << fitnessName  << 
		" FITNESS\tGRAM_DEPTH\tNN_DEPTH\tNUM_G\tNUM_C";
	if(detailed)
		os << "\tMODEL";
	os << "\n";
}

void NNModelLog::openLog(std::string filename, std::string fitnessName){
	logStream.open(filename.c_str(), ios::out);
	if(!logStream.is_open()){
		throw AthenaExcept(filename + " unable to open for writing results");
	}
	header(logStream, fitnessName);
}

///
/// Closes output stream
///
void NNModelLog::closeLog(){
		logStream.close();
}
		
///
/// Add solution to log file
///
void NNModelLog::writeSolution(NNSolution & solution, int generation, int rank){
		logStream << generation << "\t";
		logStream << rank << "\t";
		logStream << solution.fitness() << "\t";
		logStream << solution.getGramDepth() << "\t";
		logStream << solution.getNNDepth() << "\t";
		logStream << solution.getGenotypes().size() << "\t";
		logStream << solution.getCovariates().size() << "\t";
		if(detailed)
				solution.outputSolution(logStream);
		else
				logStream << "\n";
}
