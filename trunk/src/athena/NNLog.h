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
//NNLog.h

#ifndef _NNLOG_H
#define	_NNLOG_H

#include "AlgorithmLog.h"
#include <vector>

struct indModel{
	std::vector<int> snps;
	float fitness;
};

bool compareIndModels(indModel first, indModel second);

///
/// Stores data for display on algorithm processing for investigating results.
///
class NNLog: public AlgorithmLog{

	public:
	
		NNLog(int numSnps):AlgorithmLog(numSnps) {
			genNumber=0;
			genTotals newGen(numSnps, genNumber);
			gens.push_back(newGen);
			genIndex = 0;
			totalSnps = numSnps;
			maxBest = true;
		}
	
		/// Outputs data stored in this log
		void outputLog(std::ostream& os);

		/// Sets max as best (true) or false for min score is best
		void setMaxBest(bool val){maxBest = val;}

		/// Adds a generation worth of data information
		inline void addGeneration(){
				gens.clear();
				genNumber++;
				genTotals newGen(totalSnps,genNumber);
				gens.push_back(newGen);
				genIndex=0;
		}


		inline unsigned int getCurrentGen(){return genNumber;}

		inline void completeGen(){
			gens[genIndex].avgScores();
		}

		inline void addNetwork(){gens[genIndex].numValidNN++;}

		inline void addFitness(float fitness, std::vector<int> snps){
			gens[genIndex].totalFitness += fitness;
			if(fitness > gens[genIndex].maxFitness){
				gens[genIndex].maxFitness = fitness;
				if(maxBest)
					gens[genIndex].bestModelSnps = snps;
			}
			if(fitness < gens[genIndex].minFitness){
				gens[genIndex].minFitness = fitness;
				if(!maxBest)
					gens[genIndex].bestModelSnps = snps;
			}
			addSnps(snps);
			gens[genIndex].allFitness.push_back(fitness);
		}
		
		inline void addFitness(float fitness){
			gens[genIndex].totalFitness += fitness;
			if(fitness > gens[genIndex].maxFitness)
				gens[genIndex].maxFitness = fitness;
			if(fitness < gens[genIndex].minFitness)
				gens[genIndex].minFitness = fitness;
			gens[genIndex].allFitness.push_back(fitness);
		}

		inline void addEpochs(int nepochs){
			if(nepochs < gens[genIndex].minEpochs)
				gens[genIndex].minEpochs = nepochs;
			if(nepochs > gens[genIndex].maxEpochs)
				gens[genIndex].maxEpochs = nepochs;
			gens[genIndex].totalEpochs += nepochs;
		}
		
		inline void addNumGenos(int nGenos){
			gens[genIndex].numGenotypes += nGenos;
		}
		
		inline void addNumCovars(int nCovars){
			gens[genIndex].numCovariates += nCovars;
		}
		
		inline void addNNSize(int size){
			gens[genIndex].totalSize += size;
		}
		
		inline void addNNDepth(int depth){
			gens[genIndex].totalDepth += depth;
		}
		
		/// add the snps present in the network
		inline void addSnps(std::vector<int> snps){
			std::vector<int>::iterator iter;
			for(iter = snps.begin(); iter != snps.end(); iter++){
				gens[genIndex].snpTotals[*iter]++;
			}
			gens[genIndex].allModels.push_back(snps);
		}
		
		inline void setBestMod(std::vector<int> snps){
			gens[genIndex].bestModelSnps = snps;
		}

		/// output snp sizes for each
		void outputSNPSizes(std::ostream& os, unsigned int totalPopSize);
		
		/// output snp sizes headers
		void outputSNPHeaders(std::ostream& os);
		
		/// output fitnesses for all models in each generation
		void outputFitness(std::ostream& os, unsigned int totalPopSize);
		
		void outputFitnessHeaders(std::ostream& os);
		
		void outputMainHeaders(std::ostream& os);

		#ifdef PARALLEL
			void sendLog(); // for slaves
			void receiveLogs(int nprocs); // for master
			void sendDetailedLog(); // for slaves
			void receiveDetailedLogs(int nprocs); // for master
			void sendReceiveLogs(int nprocs, int myrank);
		#endif

	private:
		struct genTotals{   
			genTotals(int numSnps, int gen_num){
				generationNumber=gen_num;
				numValidNN = 0;
				numGenotypes = 0;
				numCovariates = 0;
				totalFitness = 0.0;
				totalSize = 0;
				totalDepth=0;
				maxFitness = -100000000;
				minFitness = 100000000;
				avgGenos = 0.0;
				avgCovars = 0.0;
				avgFitness = 0.0;
				avgSize = 0.0;
				snpTotals.assign(numSnps,0);
				minEpochs = 100000000;
				maxEpochs = 0;
				avgEpochs = 0.0;
				totalEpochs = 0;
			};
			
			void avgScores(){
				avgGenos = numGenotypes / float(numValidNN);
				avgCovars = numCovariates / float(numValidNN);
				avgFitness = totalFitness / float(numValidNN);
				avgSize = totalSize / float(numValidNN);
				avgEpochs = totalEpochs / float(numValidNN);
				avgDepth = totalDepth / float(numValidNN);
			}
			
			int numValidNN, numGenotypes, numCovariates, totalSize, totalEpochs, minEpochs, 
				maxEpochs, totalDepth, generationNumber;
			float avgGenos, avgCovars, avgFitness, avgSize, totalFitness, maxFitness, minFitness,
				avgEpochs, avgDepth;
			std::vector<int> snpTotals, bestModelSnps;
			std::vector<std::vector<int> > allModels;
			std::vector<float> allFitness;
		};
			 
		#ifdef PARALLEL
			void fillMasterLog(float* recData, int basicSize, float* snps, int snpSize,
				int nprocs);
			void fillSendBuffer(float* sendData);
			void packageDetailed(int& totalFitSize, float*& fitSend, int& totalSnpSize, int*& snpSend);
			void mergeDetailed(int totalFitSize, float* fitRec, int totalSnpSize, int* snpRec, 
				std::vector<std::vector<indModel> >& indModels);
		#endif

		std::vector<genTotals> gens;
		unsigned int genIndex;
		int totalSnps, genNumber;
		bool maxBest;
};

#endif
