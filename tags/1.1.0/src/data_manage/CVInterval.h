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
#ifndef CVINTERVAL_H_
#define CVINTERVAL_H_

#include "Dataset.h"

namespace data_manage
{

///
/// A crossvalidation interval with a dataset for a training, testing, and
/// an optional validation set.
///

class CVInterval
{
public:
	CVInterval();
	CVInterval(unsigned int num);
	~CVInterval();

	Dataset& get_set(unsigned int index){return sets[index];}

	unsigned int numSets(){return sets.size();}

	void numSets(unsigned int num);

	void addSet(Dataset& set){sets.push_back(set);}

	Dataset& getTraining(){return sets[0];}
	Dataset& getTesting(){return sets[1];}
	Dataset& getValidation(){return sets[2];}

private:

	std::vector<Dataset> sets;

};

}

#endif /*CVINTERVAL_H_*/
