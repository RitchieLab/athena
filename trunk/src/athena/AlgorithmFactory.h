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
/* 
 * File:   AlgorithmFactory.h
 * Author: dudeksm
 *
 * Created on November 10, 2008, 2:09 PM
 */

#ifndef _ALGORITHMFACTORY_H
#define	_ALGORITHMFACTORY_H

#include "Algorithm.h"

class AlgorithmFactory{
		
public:
		/// Creates and returns the Algorithm
		static Algorithm* createAlgorithm(std::string algName);
		
		/// Enumeration for algorithm types
		enum AlgorithmType{
				/// Enum for missing
				NoAlgorithm,
				/// Enum for GENN
				GENNAlgorithm,
				/// Enum for Symbolic Regression
				GESymbRegAlgorithm,
				/// Enum for Bayesian
				GEBayesAlgorithm,
				/// Enum for Bayesian Discriminant
				GEDiscrimBayesAlgorithm
		};
		
private:
		
		/// Sets map for Algorithm creation
		static void setAlgorithmMap();
		
		static std::map<std::string, AlgorithmType> AlgorithmMap;
		
};



#endif	/* _ALGORITHMFACTORY_H */


