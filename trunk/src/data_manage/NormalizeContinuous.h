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
#ifndef NORMALIZECONTINUOUS_H_
#define NORMALIZECONTINUOUS_H_

#include "ScaleData.h"


namespace data_manage
{

class NormalizeContinuous : public data_manage::ScaleData
{
public:
	NormalizeContinuous();
	~NormalizeContinuous();

	void adjustContin(Dataholder* holder, unsigned int varIndex);
	
	void adjustContin(Dataholder* holder);
	
	/// Adjust continuous status values
	void adjustStatus(Dataholder* holder);

	/// Returns original value for a scaled status value
	float getOriginalStatus(float status){return (status*stdDev)+meanVal;}
	
private:
	
	float stdDev, meanVal;
	
};

}

#endif /*NORMALIZECONTINUOUS_H_*/
