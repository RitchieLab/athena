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
#ifndef SCALEDATA_H_
#define SCALEDATA_H_

#include "Dataholder.h"
#include <string>

namespace data_manage
{

///
/// Abstract base class that performs scaling of data
///
class ScaleData
{
public:
	ScaleData();
	virtual ~ScaleData();

	virtual void adjustContin(Dataholder* holder){return;}

	virtual void adjustContin(Dataholder* holder, unsigned int varIndex){return;}
	
	virtual void adjustStatus(Dataholder* holder){return;}
	
	virtual void setContinGroups(Dataholder* holder){return;}
	
	virtual float getOriginalStatus(float status){return status;}

	void sigmoidStatus(Dataholder* holder);
	
	virtual std::string outputScaleInfo(){return "";}

};

}

#endif /*DATAADJUST_H_*/
