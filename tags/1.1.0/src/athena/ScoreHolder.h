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
 * File:   ScoreHolder.h
 * Author: dudeksm
 *
 * Created on September 3, 2014, 10:18 AM
 */

#include <map> 
#include "Terminals.h"
  
class ScoreHolder{
 
public:
	ScoreHolder();
	
	void clear();
	
	struct ScoreNode{
	
		ScoreNode(){
			sc = ScoreHolder::notFound();
		}
		
		std::map<IndividualTerm*, ScoreNode> scores;
		float sc;
	};
	
	ScoreNode* getScore(vector<IndividualTerm*>& terms);
	
	inline static float notFound(){return noValue;}

private:

	void initialize();

	int count, maximumCount;
	static float noValue;
	std::map<IndividualTerm*, ScoreNode> holder;
	ScoreNode none;
	
 
 };