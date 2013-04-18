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
#include "ContinFileReader.h"
#include "Stringmanip.h"
#include <fstream>
#include <sstream>
#include <iostream>

namespace data_manage
{

ContinFileReader::ContinFileReader()
{
	dummyID = 1;
}

ContinFileReader::~ContinFileReader()
{
}

///
/// Reads continuous variables from file and stores in Dataholder
/// @param filename string
/// @param holder Dataholder
/// @param missingValue value that identifies absence of data
/// @param containsID true when first column is an ID that should match the
/// ID in the corresponding genotype data file
///
void ContinFileReader::readContinFile(string filename, Dataholder* holder,
				float missingValue, bool containsID){

		holder->setMissingCoValue(missingValue);

		std::ifstream cStream(filename.c_str(), ios::in);

		if(!cStream.is_open()){
			throw DataExcept("ERROR: Unable to open " + filename + "\n");
		}

		string line, indID;
		float value;
		unsigned int currInd=0;
		unsigned int numCovars=0;

		while(!cStream.eof()){
			getline(cStream, line);

			if(line.find_first_of("0123456789") == string::npos){
				continue;
			}

			stringstream ss(line);
			
			if(containsID){
					ss >> indID;
			}
			else{
					indID = Stringmanip::numberToString(dummyID++);
			}
		 Individual * ind;
		 try{
			 ind = holder->getIndByID(indID);
		 }
		 catch(DataExcept& de){
			 cout << "Skipping individual " << indID << " while reading " << filename << endl;  
			 continue;
		 }
			
			while(ss >> value){
				ind->addCovariate(value);
			}
			if(numCovars > 0 && numCovars != ind->numCovariates()){
				throw DataExcept("ERROR: in file " +  filename +" individual " + indID +" has " 
					+ Stringmanip::numberToString(ind->numCovariates()) + " values but previous line has " 
					+ Stringmanip::numberToString(numCovars));
			}
			else{
				 numCovars = ind->numCovariates();
			}
			currInd++;
		}

		cStream.close();
		
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
void ContinFileReader::readContinFile(std::string trainFile, std::string testFile, Dataholder* holder,
	  float missingValue, bool containsID){
	  
	// parse training file and add to dataholder  
	readContinFile(trainFile, holder, missingValue, containsID);
	
	// parse the testing file and add to dataholder
	readContinFile(testFile, holder, missingValue, containsID);
	
}

}
