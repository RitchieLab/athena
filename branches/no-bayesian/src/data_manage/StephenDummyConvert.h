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
#ifndef _STEPHENDUMMYCONVERT_H_
#define _STEPHENDUMMYCONVERT_H_

#include "DummyConvert.h"

namespace data_manage
{


///
/// Converts genotypes in individuals to the ott dummy encoding.
/// Each single genotype is replaced by a pair of genotypes
///
class StephenDummyConvert: public DummyConvert{
public:
	StephenDummyConvert();
	~StephenDummyConvert();

	/// converts all genotypes to the ott-dummy encoding
	void convertGenotypes(Dataholder* holder);

private:

};

}

#endif /*_STEPHENDUMMYCONVERT_H_*/
