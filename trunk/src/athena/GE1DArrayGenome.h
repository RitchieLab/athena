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

	typedef void (*Establishinator)(GAGenome &);

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
	
	unsigned int getNumGenes()const{return numGenes;} 
	void setNumGenes(const unsigned int numGenes);
	
	unsigned int getNumCovars()const{return numCovars;} 
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
	static int pruneAndPlant(GAGenome* start, GAGenome* planted);
#endif

	static int codonMutator(GAGenome& g, float pmut);

	float getTestValue();
	void setTestValue(float val);

	int output(std::ostream & os) const;
	
	void addGeno(int g);
	void addGenos(std::vector<int> g);
	std::vector<int>& getGenos();
	void addCovar(int c);
	void addCovars(std::vector<int> c);
	std::vector<int>& getCovars();
	
	inline bool isValid() const {return validnn;}
	inline void setValid(bool val){validnn= val;}
	
	inline int getNumEpochsTrained(){return numEpochsTrained;}
	inline void setNumEpochsTrained(int e){numEpochsTrained = e;}
	
	inline void setNumIndsEvaluated(int n){numIndsEvaluated = n;}
	inline int getNumIndsEvaluated(){return numIndsEvaluated;}
	
	inline void setSSTotal(float s){ssTotal = s;}
	inline float getSSTotal(){return ssTotal;}
	unsigned int getNumNodes()const;
	void setNumNodes(const unsigned int nNodes);
	
	void setComplexity(const int c);
	int getComplexity(){return complexity;}
	
	Establishinator establishinator() const {return estab;}
  Establishinator establishinator(Establishinator f) { return(estab=f); }
  void establish();
	
	void clearScores();
	
	float getIndivScore(unsigned int indIndex){return indivScores[indIndex];}
	void setIndivScore(unsigned int indIndex, float val){indivScores[indIndex]=val;}
	void resizeIndivScores(unsigned int size){indivScores.resize(size);}
	std::vector<float> & getIndivScores(){return indivScores;}
	void clearIndivScores(){indivScores.clear();}
	
	static void setMapper(AthenaGrammarSI* m){mapper = m;}
	
	static int myRank;

private:

	void helpCopy(const GE1DArrayGenome& source);
	int helpCompare(const GE1DArrayGenome& source) const;
		
	bool validnn;	
	float ssTotal, testVal;
	unsigned int effSize, numGenes, numCovars, netDepth, gramDepth;
	int genomeID, numEpochsTrained, numIndsEvaluated, numNodes, complexity;
	
	std::vector<int> genos, covars;
	std::vector<float> indivScores;
	Establishinator estab;		// establishes newly initialized genomes that don't need full evaulation function
	static AthenaGrammarSI* mapper;
	
};

#endif

