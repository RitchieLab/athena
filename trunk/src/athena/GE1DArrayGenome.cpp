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
// GA1DArrayGenome.cpp -*- C++ -*-
#ifndef _GE1DARRAYGENOME_CPP
#define _GE1DARRAYGENOME_CPP

#include"GE1DArrayGenome.h"
#include"InitGEgenome.h"

int GE1DArrayGenome::myrank=-1;

using namespace std;

AthenaGrammarSI* GE1DArrayGenome::mapper = NULL;


GE1DArrayGenome::GE1DArrayGenome(unsigned int len)
  : GA1DArrayGenome<int>(len), effSize(0), numGenes(0), numCovars(0),  testVal(0), validnn(false), 
    numEpochsTrained(0), numIndsEvaluated(0), ssTotal(0.0), netDepth(0), gramDepth(0)
{
}

GE1DArrayGenome::GE1DArrayGenome(const GE1DArrayGenome& source)
  : GA1DArrayGenome<int>(50)
{
  this->copy(source);
}

GE1DArrayGenome::~GE1DArrayGenome()
{
}

GE1DArrayGenome* GE1DArrayGenome::clone(GAGenome::CloneMethod) const
{
  return new GE1DArrayGenome(*this);
}

void GE1DArrayGenome::copy(const GAGenome& source)
{
  if (&source != this)
  {
    GA1DArrayGenome<int>::copy(source);
    const GE1DArrayGenome& ge1DArrayGenome = 
      static_cast<const GE1DArrayGenome&>(source);
    this->helpCopy(ge1DArrayGenome);
  }
}

void GE1DArrayGenome::copy(const GE1DArrayGenome& orig,
  int r, int x, unsigned int& l){
  
  const GA1DArrayGenome<int>& g = DYN_CAST(const GA1DArrayGenome<int>&, orig);
  
  GA1DArrayGenome<int>::copy(g,r,x,l);
  
  
}


void GE1DArrayGenome::copy(GE1DArrayGenome& orig,
  int r, int x, unsigned int l){
  
  const GA1DArrayGenome<int>& g = DYN_CAST(const GA1DArrayGenome<int>&, orig);
  
  GA1DArrayGenome<int>::copy(g,r,x,l);
}


int GE1DArrayGenome::equal(const GAGenome& source) const
{
  if (GA1DArrayGenome<int>::equal(source))
  {
    return 1;
  }

  const GE1DArrayGenome& ge1DArrayGenome = static_cast<const GE1DArrayGenome&>(source);
  return this->helpCompare(ge1DArrayGenome);
}

GE1DArrayGenome GE1DArrayGenome::operator=(const GAGenome& source)
{
    this->copy(source);
    return *this;
}

void GE1DArrayGenome::helpCopy(const GE1DArrayGenome& source)
{
  effSize = source.effSize;
  testVal = source.testVal;
  validnn = source.validnn;
  numGenes = source.numGenes;
  numCovars = source.numCovars;
  netDepth = source.netDepth;
  gramDepth = source.gramDepth;
  genos = source.genos;
  covars = source.covars;
  numEpochsTrained = source.numEpochsTrained;
  numIndsEvaluated = source.numIndsEvaluated;
  ssTotal = source.ssTotal;
}

int GE1DArrayGenome::helpCompare(const GE1DArrayGenome& source) const
{
  if (&source == this)
  {
    return 1;
  }

  if (effSize == source.effSize)
  {
    return 1;
  }

  return 0;
}

int GE1DArrayGenome::getEffectiveSize() const{
	return effSize;
}
void GE1DArrayGenome::setEffectiveSize(const int newEffSize){
	effSize=newEffSize;
        _altscore=newEffSize;
}

float GE1DArrayGenome::getTestValue(){
  return testVal;
}

void GE1DArrayGenome::setTestValue(float val){
  testVal = val;
}

void GE1DArrayGenome::add_geno(int g)
{genos.push_back(g);}


void GE1DArrayGenome::add_genos(vector<int> g){
  genos=g;
}

vector<int>& GE1DArrayGenome::get_genos()
{return genos;}
  
void GE1DArrayGenome::add_covar(int c)
{covars.push_back(c);}

void GE1DArrayGenome::add_covars(vector<int> c){
  covars = c;
}

vector<int>& GE1DArrayGenome::get_covars()
{return covars;}


unsigned int GE1DArrayGenome::getNumGenes()const{
  return numGenes;
}

void GE1DArrayGenome::setNumGenes(const unsigned int nGenes){
  numGenes = nGenes;
}

unsigned int GE1DArrayGenome::getNumCovars()const{
  return numCovars;
}

void GE1DArrayGenome::setNumCovars(const unsigned int nCovars){
  numCovars = nCovars;
}

unsigned int GE1DArrayGenome::getDepth()const{
  return netDepth;
}

void GE1DArrayGenome::setDepth(const unsigned int depth){
  netDepth = depth;
}

unsigned int GE1DArrayGenome::getGramDepth()const{
  return gramDepth;
}

void GE1DArrayGenome::setGramDepth(const unsigned int depth){
  gramDepth = depth;
}


// Allow crossover only on effective length of individuals
int GE1DArrayGenome::effCrossover(const GAGenome &p1, 
			       const GAGenome &p2, 
			       GAGenome *c1, 
			       GAGenome *c2){
			       
	const GE1DArrayGenome &mom = DYN_CAST(const GE1DArrayGenome &, p1);
	const GE1DArrayGenome &dad = DYN_CAST(const GE1DArrayGenome &, p2);		       
	int nc=0;
	unsigned int momsite, momlen;
  unsigned int dadsite, dadlen;

	if(c1 && c2){
    GE1DArrayGenome &sis=DYN_CAST(GE1DArrayGenome &, *c1);
    GE1DArrayGenome &bro=DYN_CAST(GE1DArrayGenome &, *c2);	       
    
    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE &&
       bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() || 
	 sis.length() != bro.length() ||
	 sis.length() != mom.length()){
	GAErr(GA_LOC, mom.className(), "effective cross", gaErrSameLengthReqd);
	return nc;
      }
      momsite = dadsite = GARandomInt(0, mom.length()-1);
      momlen = dadlen = mom.length() - momsite;
    }
    else if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE ||
	    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      GAErr(GA_LOC, mom.className(), "effective cross", gaErrSameBehavReqd);
      return nc;
    }
    else{
    	if(mom.getEffectiveSize()>mom.length() || mom.getEffectiveSize() == 0){
		  // Wrapping was used, choose random location
    		momsite = GARandomInt(0, mom.length()-1);
    	}
    	else{
		    // Choose point from effective length
  		  momsite = GARandomInt(0, mom.getEffectiveSize()-1);
	    }
  	  if(dad.getEffectiveSize()>dad.length() || dad.getEffectiveSize() == 0){
	  	  // Wrapping was used, choose random location
		    dadsite = GARandomInt(0,dad.length()-1);
  	  }
	    else{
		    // Choose point from effective length
		    dadsite = GARandomInt(0, dad.getEffectiveSize()-1);
  	  }
      momlen = mom.length() - momsite;
      dadlen = dad.length() - dadsite;
      sis.resize(momsite+dadlen);
      bro.resize(dadsite+momlen);
    }
  
    sis.copy(mom, 0, 0, momsite);
    sis.copy(dad, momsite, dadsite, dadlen);
    bro.copy(dad, 0, 0, dadsite);
    bro.copy(mom, dadsite, momsite, momlen);
    
    nc=2;
  }
  else if(c1 || c2){
    GE1DArrayGenome &sis = (c1 ? 
			       DYN_CAST(GE1DArrayGenome &, *c1) : 
			       DYN_CAST(GE1DArrayGenome &, *c2));

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() || sis.length() != mom.length()){
	GAErr(GA_LOC, mom.className(), "effective cross", gaErrSameLengthReqd);
	return nc;
      }
      momsite = dadsite = GARandomInt(0, mom.length()-1);
      momlen = dadlen = mom.length() - momsite;
    }
    else{
    	if(mom.getEffectiveSize()>mom.length()){
		  // Wrapping was used, choose random location
    		momsite = GARandomInt(0, mom.length()-1);
    	}
    	else{
		    // Choose point from effective length
  		  momsite = GARandomInt(0, mom.getEffectiveSize()-1);
	    }
  	  if(dad.getEffectiveSize()>dad.length()){
	  	  // Wrapping was used, choose random location
		    dadsite = GARandomInt(0,dad.length()-1);
  	  }
	    else{
		    // Choose point from effective length
		    dadsite = GARandomInt(0, dad.getEffectiveSize()-1);
  	  } 
      sis.resize(momsite+dadlen);
    }
    
    if(GARandomBit()){
      sis.copy(mom, 0, 0, momsite);
      sis.copy(dad, momsite, dadsite, dadlen);
    }
    else{
      sis.copy(dad, 0, 0, dadsite);
      sis.copy(mom, dadsite, momsite, momlen);
    }

    nc = 1;
  }
  
  return nc;
}


#ifdef ATHENA_BLOAT_CONTROL
///
/// Performs prune and plant on genome passed
///
int GE1DArrayGenome::prune_and_plant(GAGenome* startGenome,
  GAGenome* planted){
 
  GAGenome * cloned = startGenome->clone();
  const GE1DArrayGenome &orig=DYN_CAST(GE1DArrayGenome &, *cloned);
  int inc=0;
  
  GE1DArrayGenome &prune=DYN_CAST(GE1DArrayGenome &, *startGenome);
  GE1DArrayGenome &plant=DYN_CAST(GE1DArrayGenome &, *planted);

  mapper->setGenotype(orig);
  Phenotype const *phenotype=mapper->getPhenotype();  

  // if not valid, make pruned a copy and use sensible init to create a planted one
  if(!phenotype->getValid()){
    unsigned int origl = orig.length();
    prune.copy(orig,0,0,origl);
    InitGEgenome::initFuncSI(plant);
    delete cloned;
    return inc=1;
  }
 
  // when valid look for an <expr> to use as basis for prune and plant
  // establish codons
  mapper->establishCodons(orig);
  
  string ruleString = "<expr>";

  unsigned int origsite = mapper->getMatchingCodon(ruleString);

  int origBlockLen = mapper->determineBlockLength(origsite);

  // create block to replace original <expr> in pruned one
  // set max depth to 5 and then run initialization to create small subtree
  unsigned int origMaxDepth = mapper->getMaxDepth();
  mapper->setMaxDepth(5);
  // use plant to hold the new subtree temporarily
  InitGEgenome::initFuncSI(plant);

  // replace original <expr> with one that uses only <v> or Concat<num>
  // as created in the plant holder
  prune.resize(orig.length()-origBlockLen+plant.length());
  prune.copy(orig,0,0,origsite);
  prune.copy(plant,int(origsite),0,(unsigned int)plant.length());
  unsigned int endpoint = orig.length()-(origsite+origBlockLen);
  prune.copy(orig,int(origsite+plant.length()),int(origsite+origBlockLen),endpoint);
 
  //resize the planted one with block length 
  //and copy over the codons to it
  plant.resize(origBlockLen);
  endpoint = (unsigned int)origBlockLen;
  plant.copy(orig, 0, origsite,endpoint);
 
  mapper->setMaxDepth(origMaxDepth);
  delete cloned;
  return inc=1;
}
#endif



///
/// Crosses over at same nonterminal point to and then swap the codons
/// for each genome at that block.  
///
int GE1DArrayGenome::blockCrossover(const GAGenome& p1,
  		  const GAGenome& p2,
			  GAGenome* c1, 
			  GAGenome* c2){
  // when one or the other is invalid (perform effCrossover instead)
  const GE1DArrayGenome &mom = DYN_CAST(const GE1DArrayGenome &, p1);
	const GE1DArrayGenome &dad = DYN_CAST(const GE1DArrayGenome &, p2);		       

  // when neither is valid use standard effective crossover
  if(!mom.isValid() && !dad.isValid()){
    return effCrossover(p1,p2,c1,c2);
  }
  
	int nc=0;
	int momsite;
  int dadsite;
	
	if(c1 && c2){
    GE1DArrayGenome &sis=DYN_CAST(GE1DArrayGenome &, *c1);
    GE1DArrayGenome &bro=DYN_CAST(GE1DArrayGenome &, *c2);

    // if either chromosome is of fixed size return error message
    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE ||
	    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      GAErr(GA_LOC, mom.className(), "block-preserve cross", gaErrSameBehavReqd);
      return nc;
    }
    else{
      
      // check for case where one is valid and other is not
      if(!mom.isValid()){
        unsigned int dadl = dad.length();
        dadsite = GARandomInt(0, dad.getEffectiveSize()-1);
        unsigned int msite = GARandomInt(0, mom.length()-1);        
        unsigned int dadlen = dad.length() - dadsite;
        sis.resize(msite + dadlen);
        sis.copy(mom, 0, 0, msite); // create new genome from dad and mom
        sis.copy(dad, msite, dadsite, dadlen);
        bro.copy(dad,0,0,dadl); // copy dad over as it is functional network
        return nc = 2;
      }
      else if(!dad.isValid()){
        unsigned int moml = mom.length();
        momsite = GARandomInt(0, mom.getEffectiveSize()-1);
        unsigned int dsite = GARandomInt(0, dad.length()-1);
        unsigned int momlen = mom.length() - momsite;
        sis.copy(mom,0,0,moml); // copy mom over
        bro.resize(dsite + momlen);
        bro.copy(dad,0,0,dsite); // create new genome
        bro.copy(mom, dsite, momsite, momlen);
        return nc = 2;
      }
    
      // select random site on the mother chromosome
      momsite = GARandomInt(0, mom.getEffectiveSize()-1);
      
      // determine the non-terminal at that location
      mapper->establishCodons(mom);
      int momBlockLen = mapper->determineBlockLength(momsite);

      // if not complete block copy original over (or can do effCrossover with these?)
      if(momBlockLen < 0){
        unsigned int moml = mom.length();
        unsigned int dadl = dad.length();
        sis.copy(mom,0,0,moml);
        bro.copy(dad,0,0,dadl);
      }
      else{
        string ruleString = mapper->getRuleString(momsite);

        // have valid block from mom so find corresponding block in dad
        mapper->establishCodons(dad);
        dadsite = mapper->getMatchingCodon(ruleString);

        int dadBlockLen = -1;
        if(dadsite >= 0)
          dadBlockLen = mapper->determineBlockLength(dadsite);
        
        if(dadBlockLen < 0){
          unsigned int moml = mom.length();
          unsigned int dadl = dad.length();
          sis.copy(mom,0,0,moml);
          bro.copy(dad,0,0,dadl);          
        }
        else{
          // everything ok so copy to the new genomes, swapping the blocks in the parents
          sis.resize(mom.length()-momBlockLen+dadBlockLen);
          unsigned int msite = momsite;
          sis.copy(mom,0,0, msite);
          unsigned int dBLockLen = dadBlockLen;
          unsigned int mBlockLen = momBlockLen;
          sis.copy(dad, momsite, dadsite, dBLockLen);
          unsigned int endpoint = mom.length()-(msite+mBlockLen);
          sis.copy(mom,momsite+dadBlockLen, momsite+momBlockLen, endpoint);
          
          // copy other genome
          bro.resize(dad.length()-dadBlockLen+momBlockLen);
          unsigned int dsite = dadsite;
          bro.copy(dad,0,0,dsite);
          bro.copy(mom,dadsite,momsite,mBlockLen);
          endpoint = dad.length()-(dadsite+dadBlockLen);
          bro.copy(dad, dadsite+momBlockLen, dadsite+dadBlockLen, endpoint);          
        }
        
      }
      nc = 2;
      }     
    }
    else if(c1 || c2){
      GE1DArrayGenome &sis = (c1 ? 
			  DYN_CAST(GE1DArrayGenome&, *c1) : 
				DYN_CAST(GE1DArrayGenome&, *c2));
				
      // when fixed size required return an error
      if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE){
	        GAErr(GA_LOC, mom.className(), "one-point cross", gaErrSameLengthReqd);
	        return nc;
      }		  
		  else{
		    // determine which will be base genome      
        bool change = GARandomBit();
        
        
        // when one parent is not valid 
        // copy the other parent as new child so that blocks are maintained
        if(!dad.isValid()){
            sis.resize(mom.length());
            unsigned int moml = mom.length();
            sis.copy(mom,0,0,moml);
            return nc=1;
        }
        if(!mom.isValid()){
            sis.resize(dad.length());
            unsigned int dadl = dad.length();
            sis.copy(dad,0,0,dadl);
            return nc=1;
        }
        
        const GE1DArrayGenome* startParent = (change?&mom:&dad);
        const GE1DArrayGenome* otherParent = (change?&dad:&mom);
        
        int startsite = GARandomInt(0, startParent->getEffectiveSize()-1);
        // determine the non-terminal at that location
        mapper->establishCodons(*startParent);
        int startBlockLen = mapper->determineBlockLength(startsite);
        // if not complete block copy original over (or can do effCrossover with these?)
        if(startBlockLen < 0){
          unsigned int startl = startParent->length();
          sis.copy(*startParent,0,0,startl);
        }
        else{
          string ruleString = mapper->getRuleString(startsite);
          // have valid block from mom so find corresponding block in dad
          mapper->establishCodons(*otherParent);
          int othersite = mapper->getMatchingCodon(ruleString);
          // when no match found copy and return starting parent
          if(othersite < 0){
              sis.resize(startParent->length());
              unsigned int startl = startParent->length();
              sis.copy(*startParent, 0, 0, startl);
              return nc=1;
          }
 
          int otherBlockLen = mapper->determineBlockLength(othersite);
        
          if(otherBlockLen < 0){
            unsigned int startl = startParent->length();
            sis.copy(*startParent,0,0,startl);
          }
          else{
            // everything ok so copy to the new genomes, swapping the blocks in the parents
            sis.resize(startParent->length()-startBlockLen+otherBlockLen);
            unsigned int ssite = startsite;
            sis.copy(*startParent,0,0, ssite);
            unsigned int oBLockLen = otherBlockLen;
            unsigned int sBlockLen = startBlockLen;
            sis.copy(*otherParent, startsite, othersite, oBLockLen);
            unsigned int endpoint = startParent->length()-(ssite+sBlockLen);
            sis.copy(*startParent,startsite+otherBlockLen, startsite+startBlockLen, endpoint);
          }
        }
              
		  }
    
      nc = 1;
    }
  return nc;
}



// Calculate number of mutations based on length and do those
int GE1DArrayGenome::codonMutator(GAGenome & g, float pmut){
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;
  
  // if pmut is zero return 0 for no mutations done
  if(pmut <= 0.0)
    return 0;
  
  float nmut = genome.length() * pmut;
  register int i,n;
  
  int length = genome.length()-1;
  
  // when total is less than one have to check each codon
  if(nmut < 1.0){
  	for(i=length; i>=0; i--){
  		if(GAFlipCoin(pmut))
	  		genome.gene(i, GARandomInt(0,INT_MAX));
  	}
  }
  else{
  	for(n=0; n<nmut; n++)
  		genome.gene(GARandomInt(0, length), GARandomInt(0,INT_MAX));
  }
}
 
int
GE1DArrayGenome::output(ostream & os) const
{
    int len = length();
    os << "OUTPUT Genome ";
    for(int i=0; i<len; i++){
      os << (unsigned int)gene(i) << endl;
    }
  return os.fail() ? 1 : 0;
}


void GE1DArrayGenome::clearScores(){
    std::vector<int> temp;
    setNumCovars(0);
    setNumGenes(0);
    add_genos(temp);
    add_covars(temp);
    setNumIndsEvaluated(0);
    setDepth(0);
    setGramDepth(0);
}

#endif

