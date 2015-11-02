/*
Copyright Marylyn Ritchie 2014

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
//BayesModelLog.cpp

#include "BayesModelLog.h"
#include <iomanip>

///
/// Destructor
///
BayesModelLog::~BayesModelLog(){
		if(logStream.is_open()){
				closeLog();
		}
}

void BayesModelLog::header(ostream& os, std::string fitnessName){
	os << "GEN\tRANK\t" << fitnessName  << 
		" FITNESS\tCOMPLEXITY\tNUM_G\tNUM_C";
	if(detailed)
		os << "\tMODEL";
	os << "\n";
}

void BayesModelLog::openLog(std::string filename, std::string fitnessName, bool varsOnly){
	logStream.open(filename.c_str(), ios::out);
	if(!logStream.is_open()){
		throw AthenaExcept(filename + " unable to open for writing results");
	}
	if(!varsOnly)
		header(logStream, fitnessName);
}

void BayesModelLog::addGeneration(int generation){
	logStream << "Gen: " << generation << "\n";
}

///
/// Closes output stream
///
void BayesModelLog::closeLog(){
		logStream.close();
}

///
/// Writes variables to log 
/// @param solution to write
/// @param whether solution has genotypes converted to Ott representation 
/// @param data for converting to variable names
///
void BayesModelLog::writeVariables(BayesSolution& solution, bool ottDummyConverted){
	vector<int> genotypes=solution.getGenotypes(ottDummyConverted);
	vector<int> covariates=solution.getCovariates();
	
	vector<int>::iterator iter;
	
	iter=genotypes.begin();
	if(iter != genotypes.end()){
		logStream << "G" << *iter;
		++iter;
	}
	for(; iter != genotypes.end(); iter++){
		logStream << " G" << *iter;
	}
	iter=covariates.begin();
	if(iter != covariates.end()){
		logStream << "C" << *iter;
		++iter;
	}	
	for(iter=covariates.begin(); iter != covariates.end(); iter++){
		logStream << " C"<< *iter;
	}
	if(!genotypes.empty() || !covariates.empty())
		logStream << "\n";	
}

///
/// Add solution to log file
///
void BayesModelLog::writeSolution(BayesSolution & solution, int generation, int rank){
		logStream << generation << "\t";
		logStream << rank << "\t";
		logStream << solution.fitness() << "\t";
		logStream << solution.getComplexity() << "\t";
		logStream << solution.getGenotypes().size() << "\t";
		logStream << solution.getCovariates().size() << "\t";
		if(detailed)
				solution.outputSolution(logStream);
		else
				logStream << "\n";
}
