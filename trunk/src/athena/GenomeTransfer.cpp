#include "GenomeTransfer.h"
#ifdef HAVE_CXX_MPI
#include "Algorithm.h"
#include "InitGEgenome.h"
#include "GEObjective.h"

GenomeTransfer::GenomeTransfer(){
	genomeInfo = 112;
	genomeArray=113;
}


void GenomeTransfer::setRank(int rank){
				InitGEgenome::setrank(rank);
				GEObjective::setrank(rank);
}


///
/// Returns sum of all the fitness goal checks on all nodes.  A 0 means no node has reached
/// the fitness.  A value > 0 means that many populations found a good enough fit.
/// @param rank
/// @param complete
///
int GenomeTransfer::nodesCompleted(int complete){
	int sum_completed=0;
	MPI_Allreduce(&complete, &sum_completed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	return sum_completed;
}

///
/// Alternate MPI messaging taking advantage of MPI_Allgather functionality using struct
/// @param totalNodes
/// @param myRank
///
void GenomeTransfer::sendAndReceiveStruct(int totalNodes, int myRank, GASimpleGA* ga){

	structMPI * send = new structMPI;
	
	GE1DArrayGenome& genome = (GE1DArrayGenome&)ga->statistics().bestIndividual();
	// package genome Info
	send->genomeParams[0] = genome.length();
	send->genomeParams[1] = genome.score();
	send->genomeParams[2] = genome.getEffectiveSize();
	send->genomeParams[3] = genome.getNumGenes();
	send->genomeParams[4] = genome.getNumCovars();
	send->genomeParams[5] = genome.getNumIndsEvaluated();
	send->genomeParams[6] = genome.getSSTotal();
	send->genomeParams[7] = genome.getNumNodes();
	
	// package codons
	for(int i=0; i<genome.length(); i++){
		send->codons[i] = genome.gene(i);
	}
	
	// prepare receiving array
	structMPI * recv = new structMPI[totalNodes];
	
	MPI_Allgather (send, sizeof(*send), MPI_BYTE, recv, sizeof(*send), MPI_BYTE, MPI_COMM_WORLD);
	
	updateWithMigration(recv, totalNodes, myRank, ga);
	
	delete send;
	delete [] recv;
	
}

///
/// Broadcast best genotypes and continuous variables to all other nodes
///
void GenomeTransfer::exchangeBestVariables(int totalNodes, int myRank, vector<int>& genotypes,
  vector<int>& contins){

  int nVars = 1000;
  int * variables = new int[nVars];

  vector<int>::iterator iter;
  if(myRank==0){
    int currVar=0;
    variables[currVar++]=int(genotypes.size());
    for(iter=genotypes.begin(); iter!=genotypes.end(); iter++){
      variables[currVar++]=*iter;
    }
    variables[currVar++]=int(contins.size());
    for(iter=contins.begin(); iter!=contins.end(); iter++){
      variables[currVar++]=*iter;
    }    
  }
  MPI_Bcast(variables, nVars, MPI_INT, 0, MPI_COMM_WORLD);
  
  genotypes.clear();
  contins.clear();
  
  int currVar=0;
  int n = variables[currVar++];
  for(int i=0; i<n; i++){
    genotypes.push_back(variables[currVar++]);
  }
  n=variables[currVar++];
  for(int i=0; i<n; i++){
    contins.push_back(variables[currVar++]);
  }
}



///
/// Incorporates migration into population
/// @param genomes 
/// @param totalNodes
/// @param myRank
///
void GenomeTransfer::updateWithMigration(structMPI* mpiGenomes, int totalNodes, int myRank,
	GASimpleGA* ga){
	GAPopulation pop(ga->population());
	
	for(int node=0; node < totalNodes; node++){
		if(myRank==node){
			continue;
		}
		
		GAGenome *tmpInd = ga->population().individual(0).clone();
		GE1DArrayGenome& genome = (GE1DArrayGenome&)*tmpInd;
		int len = mpiGenomes[node].genomeParams[0];
		genome.length(len);
		genome.setEffectiveSize(mpiGenomes[node].genomeParams[2]);
		genome.setNumGenes(mpiGenomes[node].genomeParams[3]);
		genome.setNumCovars(mpiGenomes[node].genomeParams[4]);
		genome.setNumIndsEvaluated(mpiGenomes[node].genomeParams[5]);
		genome.setSSTotal(mpiGenomes[node].genomeParams[6]);
		genome.setNumNodes(mpiGenomes[node].genomeParams[7]);
		for(int i=0; i<len; i++){
			genome.gene(i, mpiGenomes[node].codons[i]);
		}
		genome.score(mpiGenomes[node].genomeParams[1]);

		pop.add(genome);
		
		delete tmpInd;
	}

	// remove worst individuals from population
	for(int i=0; i < totalNodes-1; i++)
		pop.destroy();
		
	ga->population(pop);

}

#endif /* HAVE_CXX_MPI */
