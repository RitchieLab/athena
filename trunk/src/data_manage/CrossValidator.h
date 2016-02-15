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
#ifndef CROSSVALIDATOR_H_
#define CROSSVALIDATOR_H_

#include <cstdlib>

#include "CVSet.h"
#include "Dataholder.h"

namespace data_manage
{


class TRandom
{
	public:
	int operator()(int n)
	{
		return rand() % n;
	}

};

///
///  Performs crossvalidation on a dataholder to return datasets with the
///

class CrossValidator
{
public:
	CrossValidator();
	~CrossValidator();

	/// splits data into testing and training sets for the indicated number of intervals
	CVSet splitData(unsigned int numCrossVal, Dataholder* holder, bool multiCatData=false);

	/// save splits with individual IDs
	void saveSplits(std::string filename);

	/// load splits with individual IDs
	CVSet loadSplits(std::string filename, Dataholder* holder);

	/// splits data by multiple category sets
	CVSet splitDataMultiCat(unsigned int numCrossVal, Dataholder* holder);

private:

	void distributeInds(unsigned int numSplits, std::vector<Individual*>& inds,
			std::vector<std::vector<Individual*> >& splits);

	void shuffleInds(std::vector<Individual*> & inds);

	void statusBin(Dataholder* holder, std::vector<Individual*>& affected,
			std::vector<Individual*>& unaffected);

	void statusBin(Dataholder* holder, vector<vector<Individual*> >& categories);

	void statusBinAffOnly(Dataholder* holder, std::vector<Individual*>& affected,
			std::vector<Individual*>& unaffected);

	CVSet splitByNum(Dataholder* holder);

	CVSet createSet(unsigned int num_cv, Dataholder* holder);

	vector<vector<Individual*> > splits;

};

}

#endif /*CROSSVALIDATOR_H_*/
