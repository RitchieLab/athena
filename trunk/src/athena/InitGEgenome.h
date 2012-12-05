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

