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
#ifndef _GAPARETORANKSELECTOR_H
#define	_GAPARETORANKSELECTOR_H

#include <ga/GASelector.h>
#include <ga/GAPopulation.h>

#include <deque>
#include <map>
#include <set>

class GAGenome;
class GAPopulation;


class GAParetoRankSelector : public GASelectionScheme {
public:
  GADefineIdentity("GAParetoSelector", GAID::GAParetoSelection);

  GAParetoRankSelector(int w=GASelectionScheme::RAW) : 
  GASelectionScheme(w), n(0), nComplexities(0), cInd(0), nSelect(0), cSelectedInd(0) 
    { currInd = &cInd;
      currSelectedInd=&cSelectedInd;}
    
  GAParetoRankSelector(const GAParetoRankSelector& orig) 
    { copy(orig); }
    
  GAParetoRankSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
    
  virtual ~GAParetoRankSelector() {}
  
  virtual GASelectionScheme* clone() const
    { return new GAParetoRankSelector; }
    
  virtual void copy(const GASelectionScheme& orig) {
    GASelectionScheme::copy(orig);
    const GAParetoRankSelector& sel = dynamic_cast<const GAParetoRankSelector&>(orig);
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