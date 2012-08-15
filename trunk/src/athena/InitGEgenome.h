/* 
 * File:   InitGEgenome.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 11:32 AM
 */

#ifndef _INITGEGENOME_H
#define	_INITGEGENOME_H

#include <GE/ge.h>
#include <ga/ga.h>
#include "AthenaGrammarSI.h"
#include "AthenaExcept.h"

///
/// Contains functions for initializing genomes in the GA libaray when
/// used with the GE library.
///

class InitGEgenome{
    
public:
    
    /// conducts random intialization of genomes
    static void initFuncRandom(GAGenome &g);
    
    /// conducts sensible initialization of genomes
    static void initFuncSI(GAGenome &g);
    
    /// set minimum size in random initialization
    static void setMinSize(unsigned int msize){min_size = msize;}
    
    /// set maximum size in random initialization
    static void setMaxSize(unsigned int msize){max_size = msize;}
    
    /// sets the mapper to use
    static void setMapper(AthenaGrammarSI* m){mapper = m;}
    
    /// sets rank
    static void setrank(int r){rank=r;}
    
private:
    
    static AthenaGrammarSI* mapper;
    static unsigned int min_size, max_size;
    
    static int rank;
    
};


#endif	/* _INITGEGENOME_H */

