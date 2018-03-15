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
 * File:   GenomeTransfer.h
 * Author: dudeksm
 *
 * Created on Fri May  9 11:45:50 CDT 2014 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_CXX_MPI

#ifndef _GENOMETRANSFER_H
#define	_GENOMETRANSFER_H

#include "Config.h"
#include "Structs.h"
#include <ga/ga.h>
#include <mpi.h>

struct genomeMPI{
	float genomeParams[8];
	int codons[MAX_GENOME_SIZE];
};
		
typedef struct genomeMPI structMPI;


class GenomeTransfer{
		
public:

		GenomeTransfer();

		/// update population with migration
		void updateWithMigration(structMPI* mpiGenomes, int totalNodes, int myRank, GASimpleGA* ga);
		
		/// exhange best variables between populations
		void exchangeBestVariables(int totalNodes, int myRank, vector<int>& genotypes,
		  vector<int>& contins);
		
		/// pass structures with genomes  
		void sendAndReceiveStruct(int totalNodes, int myRank, GASimpleGA* ga);
		
		/// check for any nodes completed
		int nodesCompleted(int complete);
		
		/// set ranks of process
		void setRank(int rank);
		
private:
		int genomeInfo, genomeArray;
		
};

#endif	/* _GENOMETRANSFER_H */

#endif /* HAVE_CXX_MPI */
