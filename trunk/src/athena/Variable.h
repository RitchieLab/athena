/*
Copyright Marylyn Ritchie 2015

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

#ifndef _VARIABLE_H
#define	_VARIABLE_H

#include <Dataholder.h>

///
/// Provides objective function for use with GALIb and GE library
///
class Variable{

public:
	Variable(int idx);
	virtual ~Variable(){}
	virtual float getValue(data_manage::Individual* ind)=0;
	int getIndex(){return index;}
	virtual bool isGeno(){return gType;}
	virtual int getNumLevels(){return nLevels;}
	void setNumLevels(int nl){nLevels=nl;}
	int getMissingVal(){return nLevels-1;}
	virtual string getName(data_manage::Dataholder* holder)=0;

protected:
	int index, nLevels;
	bool gType;
};

class GenoVariable :public Variable{
	public:
		GenoVariable(int index);
		virtual float getValue(data_manage::Individual* ind);
		virtual string getName(data_manage::Dataholder* holder);
// 		virtual float getMissingVal();
};

class ConVariable :public Variable{
	public:
		ConVariable(int index);
		virtual float getValue(data_manage::Individual* ind);
		virtual string getName(data_manage::Dataholder* holder);
// 		virtual float getMissingVal();
};

#endif