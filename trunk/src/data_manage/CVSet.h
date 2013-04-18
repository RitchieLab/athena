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
#ifndef CVSET_H_
#define CVSET_H_

#include "CVInterval.h"

namespace data_manage
{

class CVSet
{
public:
	CVSet();
	~CVSet();

	/// Get number of intervals in set
	unsigned int numIntervals(){return intervals.size();}

	/// Get specific interval from set
	CVInterval&  getInterval(unsigned int index){return intervals[index];}

	/// Add interval to set
	void addInterval(CVInterval& interval){intervals.push_back(interval);}


private:
	std::vector<CVInterval> intervals;

};

}

#endif /*CVSET_H_*/
