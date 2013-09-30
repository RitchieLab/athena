//DummyConvert.h

/*
Copyright Marylyn Ritchie 2013

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
#ifndef _DUMMYCONVERT_H_
#define _DUMMYCONVERT_H_

#include "Dataholder.h"

namespace data_manage
{


///
/// Converts genotypes in individuals to the ott dummy encoding.
/// Each single genotype is replaced by a pair of genotypes
///
class DummyConvert
{
public:

	DummyConvert(){multVars=false;}
	virtual ~DummyConvert(){}

	/// converts all genotypes to the dummy encoding
	virtual void convertGenotypes(Dataholder* holder){};

	/// are genotypes expanded to 2 variables in converted set?
	bool getMultipleVars(){return multVars;}
	
	void setMultipleVars(bool m){multVars=m;}

private:
	bool multVars;

};

}

#endif /*_DUMMYCONVERT_H_*/
