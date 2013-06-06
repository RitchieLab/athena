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
#ifndef SCALECONTINUOUS_H_
#define SCALECONTINUOUS_H_

#include "ScaleData.h"

namespace data_manage
{

class ScaleContinuous : public data_manage::ScaleData
{
public:
	ScaleContinuous();
	~ScaleContinuous();

	void adjustContin(Dataholder* holder, unsigned int varIndex);
	
	void adjustContin(Dataholder* holder);
	
	/// Scales status by dividing all 
	void adjustStatus(Dataholder* holder);

	/// Returns original value for a scaled status value
	float getOriginalStatus(float status){return status*(statMax-statMin)+statMin;}
	
	/// Writes scaling information to ostream
	string outputScaleInfo();
	
private:
	
	float statMin, statMax;
	std::map<std::string, float> maxGroupValues;
	

};

}

#endif /*SCALECONTINUOUS_H_*/
