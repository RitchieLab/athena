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
#include "OutputSet.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
using namespace std;

namespace data_manage
{

///
/// outputs designated set
///
void OutputSet::outputSet(string name, Dataset& set){

		ofstream outfile;
		outfile.open(name.c_str(), ios::out);    
		outfile << set;
		outfile.close();
	 
}

///
/// outputs set and adds cv number to name
///
void OutputSet::outputCV(string name, Dataset& set, int cv){

	stringstream ss;
	ss << name << "." << cv << ".txt";
	outputSet(ss.str(), set);

}

}
