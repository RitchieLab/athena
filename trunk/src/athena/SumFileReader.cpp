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

#include "SumFileReader.h"
#include "GEObjective.h"
#include <fstream>
#include <sstream>


using namespace std;

///
/// Read summary file designated and store ATHENA model representations
/// in population
///
///
void SumFileReader::readSumFile(std::string filename){

  ifstream sumStream(filename.c_str(), ios::in);
  if(!sumStream.is_open()){
    throw AthenaExcept("Failed in attempt to open file " + filename
       + "\n");
  }
  modelPop.clear();

  string line;
  // create new solutions using GEObjective and fill the population
  // search for start of internal representations
  while(sumStream.good()){
    getline(sumStream,line);
    size_t found = line.find("Internal");
    if(found != string::npos)
      break;
  }

  string modelnum, modelstring;
  Solution * solPtr;
  // grab if not empty and create model
  while(sumStream.good()){
    getline(sumStream,line);
    stringstream ss(line);
    ss >> modelnum;
    if(line.find_first_of("0123456789") == string::npos)
      continue;
    // construct solution
    solPtr = GEObjective::getBlankSolution();
    // convert string into vector
    vector<string> symbols;
    while(ss >> modelstring){
    	symbols.push_back(modelstring);
    }
    solPtr->setSymbols(symbols);
    modelPop.push_back(solPtr);
  }
  sumStream.close();
}
