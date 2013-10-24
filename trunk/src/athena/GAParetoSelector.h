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
#ifndef _GAPARETOSELECTOR_H
#define	_GAPARETOSELECTOR_H

#include <ga/GASelector.h>
#include <ga/GAPopulation.h>

#include <deque>
#include <map>
#include <set>

class GAGenome;
class GAPopulation;

class paretoValue{
	public: 
	
	paretoValue(){
		probability=0.0;
		selected=0;
		complexity=0;
	}
	
	paretoValue(int s, float p, int c){
		probability=p;
		selected=s;
		complexity=c;
	}

	bool operator< (const paretoValue& other) const {
		return probability < other.probability;
  }

	float probability;
	int complexity, selected;
};


class GAParetoSelector : public GASelectionScheme {
public:
  GADefineIdentity("GAParetoSelector", GAID::GAParetoSelection);

  GAParetoSelector(int w=GASelectionScheme::RAW) : 
  GASelectionScheme(w), n(0), nComplexities(0), cInd(0), nSelect(0), cSelectedInd(0) 
    { currInd = &cInd;
      currSelectedInd=&cSelectedInd;}
    
  GAParetoSelector(const GAParetoSelector& orig) 
    { copy(orig); }
    
  GAParetoSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
    
  virtual ~GAParetoSelector() {}
  
  virtual GASelectionScheme* clone() const
    { return new GAParetoSelector; }
    
  virtual void copy(const GASelectionScheme& orig) {
    GASelectionScheme::copy(orig);
//  GAParetoSelector& sel = 
//       DYN_CAST(GAParetoSelector&,orig);
    const GAParetoSelector& sel = dynamic_cast<const GAParetoSelector&>(orig);
		n=sel.n;
		nSelect=sel.n;
		nComplexities=sel.nComplexities;
		selectedInds=sel.selectedInds;
		inds=sel.inds;
		complexitySet=sel.complexitySet;
  }
  
  virtual GAGenome& select() const;
  virtual void update();

protected:

	bool addInd(int i);

	// key is complexity and value is deque of indexes into population
	std::map<int, std::deque<int> > inds;
	// list of inds that are selected
	std::deque<int> selectedInds;

	// set contains list of current valid complexities for selection
	std::set<int> complexitySet;

  int n, nComplexities, cInd;
  int * currInd;
  unsigned int *currSelectedInd;
  unsigned int nSelect, cSelectedInd;
  GAPopulation::SortBasis sortBasis;
};

#endif