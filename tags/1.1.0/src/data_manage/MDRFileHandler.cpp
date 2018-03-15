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
#include "MDRFileHandler.h"
#include <fstream>
#include <sstream>
#include "Stringmanip.h"
#include <math.h>
#include <float.h>

#include <iostream>
#include <cstdlib>


namespace data_manage
{

MDRFileHandler::MDRFileHandler()
{
	dummyID = 1;
}

MDRFileHandler::~MDRFileHandler()
{
}



///
///Parses the data file and fills this object
///@param datafile name of data file
///@param holder DataHolder to fill with
///data from the datafile
///@param missingValue value of missing data in current set
///@param statusMissingValue value that indicates missing status for an individual
///@param containds_id true when file has first column as ID value
///@return none
///@throws DataExcept on error
///
void MDRFileHandler::parseFile(string filename, Dataholder * holder,
				int missingValue, float statusMissingValue, bool containsID){

	int maxLocusValue = 0;
	ifstream dataStream(filename.c_str(), ios::in);

	if(!dataStream.is_open()){
		throw DataExcept("Error:  Unable to open " + filename + "\n");
	}

	string line, out, indID;
	string::size_type lastNumberPos;
	double indStatus;
	int geno;
	char missingGeno = 3;

	getline(dataStream, line);

	// format is staus followed by genotypes separated by white space
	bool anyMissing = false;
	// assume case/control until find an individual that isn't
	holder->setCasecontrol(true);

	do{
		if(line.find_first_of("0123456789") == string::npos || line.find("#")==0){
			getline(dataStream, line);
			continue;
		}
		// remove windows carriage return
		lastNumberPos = line.find_last_of("0123456789");
		line = line.substr(0,lastNumberPos+1);

		stringstream ss(line);
		Individual ind;
		if(containsID){
				ss >> indID;
		}
		else{
				indID = Stringmanip::numberToString(dummyID++);
		}
		ind.setID(indID);
		ss >> indStatus;

		// check to see if ind should be skipped because status is missing
		if(fabs(indStatus - statusMissingValue) > DBL_EPSILON){

			if(fabs(indStatus-0) > DBL_EPSILON && fabs(indStatus-1) > DBL_EPSILON){
				holder->setCasecontrol(false);
			}

			ind.setStatus(indStatus);
			while(ss >> geno){
				if(geno == missingValue){
					ind.addGenotype(missingGeno);
					anyMissing = true;
				}
				else if(geno > 2){
				throw DataExcept("All genotypes must be less than 3 or equal the missing value " +
						Stringmanip::numberToString(missingValue) + "\nIndividual " +
						Stringmanip::numberToString(dummyID++) +
						" has the genotype " + Stringmanip::numberToString(geno));
				}
				else{
					if(geno > maxLocusValue){
						maxLocusValue = geno;
					}
					ind.addGenotype(geno);
				}
			}
			holder->addInd(ind);
		}
		else{
			cout << "Skipping individual " << indID << " with status " << indStatus <<
	  		" (STATUSMISSINGVALUE is set to " << statusMissingValue << ")" << endl;
		}

		getline(dataStream, line);
	}while(!dataStream.eof());

	dataStream.close();

	holder->setMaxLocusValue(maxLocusValue);
	holder->anyMissingGenos(anyMissing);
	holder->setMissingGenotype(missingGeno);
}


///
///Parses the data file and fills this object
///@param datafile name of data file
///@param holder DataHolder to fill with
///data from the datafile
///@param missingValue value of missing data in current set
///@param containds_id true when file has first column as ID value
///@return none
///@throws DataExcept on error
///
void MDRFileHandler::parseFile(std::string trainFile, std::string testFile, Dataholder* holder,
	  int missingValue, float statusMissingValue, bool containsID){

	parseFile(trainFile, holder, missingValue, statusMissingValue, containsID);
	// now the number of inds in the set will specify the split point for CV generation
	holder->setTestSplit(holder->numInds());

	// parse the testing file and add to dataholder
	parseFile(testFile, holder, missingValue, statusMissingValue, containsID);

}



///
///Writes the data file to a text file using
///standard MDR format text file format.
///@param datafile name of data file
///@param dataholder DataHolder containing
///data to write to file
///@return none
///
void MDRFileHandler::writeFile(string filename,
	Dataholder* dataholder){

	std::ofstream outStream(filename.c_str(), ios::out);
	if(!outStream.is_open()){
		throw DataExcept("ERROR: Unable to open " + filename + " for writing\n");
	}

	unsigned int numInds = dataholder->numInds();
	unsigned int lastLocus = dataholder->numGenos();
	unsigned int currLoc;

	Individual* currIndividual;

	for(unsigned int currInd=0; currInd < numInds; currInd++){
		currIndividual = dataholder->getInd(currInd);
		outStream << Stringmanip::numberToString(currIndividual->getStatus())	<< " ";
		for(currLoc=0; currLoc < lastLocus; currLoc++){
			outStream << Stringmanip::numberToString(currIndividual->getGenotype(currLoc)) << " ";

		}
		outStream << Stringmanip::numberToString(currIndividual->getGenotype(currLoc)) << endl;
	}

	outStream.close();
}

}
