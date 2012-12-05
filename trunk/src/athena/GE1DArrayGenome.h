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
// GE1DArrayGenome.h -*- C++ -*-
#ifndef _GE1DARRAYGENOME_H
#define _GE1DARRAYGENOME_H

#include<ga/ga.h>
#include<iostream>
#include<vector>
#include "AthenaGrammarSI.h"

class GE1DArrayGenome:public GA1DArrayGenome<int>{
public:

  GE1DArrayGenome(unsigned int len);
  GE1DArrayGenome(const GE1DArrayGenome& source);
  virtual ~GE1DArrayGenome();

  virtual GE1DArrayGenome* clone(GAGenome::CloneMethod method) const;
  virtual void copy(const GAGenome& source);
  virtual int equal(const GAGenome& source) const;

  void copy(const GE1DArrayGenome& orig,
    int r, int x, unsigned int& l);

  void copy(GE1DArrayGenome& orig,
    int r, int x, unsigned int l);

  GE1DArrayGenome operator=(const GAGenome& source);

  int getEffectiveSize() const;
  void setEffectiveSize(const int newEffSize);
  
  unsigned int getNumGenes()const;
  void setNumGenes(const unsigned int numGenes);
  
  unsigned int getNumCovars()const;
  void setNumCovars(const unsigned int numCovars);
  
  unsigned int getDepth() const;
  void setDepth(const unsigned int depth);
  
  unsigned int getGramDepth() const;
  void setGramDepth(const unsigned int depth);
  
  static int effCrossover(const GAGenome& p1,
			  const GAGenome& p2,
			  GAGenome* c1, 
			  GAGenome* c2);
		  
	static int blockCrossover(const GAGenome& p1,
  		  const GAGenome& p2,
			  GAGenome* c1, 
			  GAGenome* c2);

#ifdef ATHENA_BLOAT_CONTROL
  static int prune_and_plant(GAGenome* start, GAGenome* planted);
#endif

  static int codonMutator(GAGenome& g, float pmut);

  float getTestValue();
  void setTestValue(float val);

  int output(std::ostream & os) const;
  
  void add_geno(int g);
  void add_genos(std::vector<int> g);
  std::vector<int>& get_genos();
  void add_covar(int c);
  void add_covars(std::vector<int> c);
  std::vector<int>& get_covars();
  
  inline bool isValid() const {return validnn;}
  inline void setValid(bool val){validnn= val;}
  
  inline int getNumEpochsTrained(){return numEpochsTrained;}
  inline void setNumEpochsTrained(int e){numEpochsTrained = e;}
  
  inline void setNumIndsEvaluated(int n){numIndsEvaluated = n;}
  inline int getNumIndsEvaluated(){return numIndsEvaluated;}
  
  inline void setSSTotal(float s){ssTotal = s;}
  inline float getSSTotal(){return ssTotal;}
  
  void clearScores();
  
  static void setMapper(AthenaGrammarSI* m){mapper = m;}
  
  static int myrank;

private:

  void helpCopy(const GE1DArrayGenome& source);
  int helpCompare(const GE1DArrayGenome& source) const;

  unsigned int effSize, numGenes, numCovars, netDepth, gramDepth;
  float testVal;
  int genome_id;
  bool validnn;
  int numEpochsTrained;
  int numIndsEvaluated;
  float ssTotal;
  std::vector<int> genos, covars;
  
  static AthenaGrammarSI* mapper;
  
};

#endif

