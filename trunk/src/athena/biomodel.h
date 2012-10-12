//
// C++ Implementation: biomodel.h
//
// Description: Encapsulates any functionality duplicated in Snp/Snp or Gene/Gene models
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
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
#ifndef __BIOFILTER_BIOMODEL_H
#define __BIOFILTER_BIOMODEL_H

namespace Biofilter {

class BioModel {
public:
	BioModel();
	BioModel(float implIdx);
	virtual ~BioModel();

	virtual float ImplicationIndex() const;
protected:
	float implicationIndex;
};

inline
BioModel::BioModel() { }

inline
BioModel::BioModel(float implIndex) : implicationIndex(implIndex) { }

inline
BioModel::~BioModel() { }
	
inline
float BioModel::ImplicationIndex() const {
	return implicationIndex;
}

}

#endif

