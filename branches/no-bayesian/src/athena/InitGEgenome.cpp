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
#include "InitGEgenome.h"

// GEGrammarSI* InitGEgenome::mapper = NULL;
AthenaGrammarSI* InitGEgenome::mapper = NULL;
unsigned int InitGEgenome::minSize = 50;
unsigned int InitGEgenome::maxSize = 200;

int InitGEgenome::rank = 0;

/// 
/// Random initializer.  Called from within GA library.
/// Adapted from an example in the GE library.
/// @param g GAGenome to initialize
///
void InitGEgenome::initFuncRandom(GAGenome &g){
	GA1DArrayGenome<unsigned int>& genome = static_cast<GA1DArrayGenome<unsigned int> &>(g);
	int n=GARandomInt(minSize,maxSize);
	genome.resize(n);
	for(int i=0; i<n; i++){
		genome.gene(i, GARandomInt(0,2560000));
	}
}


///
/// Sensible initializer for the GA library used as 
/// part of GE.  Adapted from an example in the GE library.
/// @param g GAGenome to initialize
///
void InitGEgenome::initFuncSI(GAGenome& g){

	GA1DArrayGenome<unsigned int> &genome=
		static_cast<GA1DArrayGenome<unsigned int> &>(g);
	
	if(!mapper->init()){
		cerr << "Error using sensible initialisation.\n";
		cerr << "Execution aborted.\n";
		exit(0);
	}
	
	unsigned int genomeSize = mapper->getGenotype()->size();

	genome.resize(genomeSize);
	
	if(genomeSize > (unsigned int)genome.size()){
		throw AthenaExcept("Sensible initialization is producing models larger than the maximum genome size.  Try reducing the MAXDEPTH or increasing the MAXSIZE");
	}
	
	int i=0;
 
 	// Now copy genotype onto genome
	Genotype::const_iterator genIt=(mapper->getGenotype())->begin();
	while(genIt!=(mapper->getGenotype())->end()){
	  genome.gene(i, *genIt);
	  genIt++;
	  i++;
	}
}
