/*
Copyright Marylyn Ritchie 2014

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
/* 
 * File:   MDR.h
 * Author: dudeksm
 *
 * Created on August 22, 2014, 3:18 PM
 */
 
#ifndef _MDR_H
#define _MDR_H

#include <Dataset.h>
#include "Terminals.h"

class MDR{

public:

	static float calcBalAccuracy(Dataset* dSet, Dataset* refSet, vector<IndividualTerm*> &variables);


private:
 	MDR();
 
	static void totalVariables(Dataset* dSet, Dataset* refSet, vector<IndividualTerm*> &variables,
		vector<vector<int> >& totals, vector<vector<int> >& conditTotals);
 
};
 
#endif
 
 
 