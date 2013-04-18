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
//NNLog.cpp

#include "NNLog.h"
#include <algorithm>
#include <iomanip>

using namespace std;

bool compareIndModels(indModel first, indModel second){
	if(first.fitness > second.fitness)
		return true;
	else
		return false;
}

///
/// Ouputs log information to stream 
/// os ostream
///
void NNLog::outputLog(ostream& os){

	os.precision(3);
	for(unsigned int gen=0; gen < gens.size(); gen++){
		os << setw(5) << left << gens[gen].generationNumber << setw(11) << gens[gen].numValidNN << setw(10) << gens[gen].avgSize 
			<< setw(9) << gens[gen].avgDepth << setw(8) << gens[gen].avgFitness 
			<< setw(9) << gens[gen].maxFitness << setw(9) << gens[gen].minFitness
			<< setw(9) << gens[gen].avgGenos << setw(7) << gens[gen].avgCovars 
			<< setw(11) << gens[gen].avgEpochs << setw(10) << gens[gen].maxEpochs <<  " ";
		
		if(!gens[gen].bestModelSnps.empty()){  
				os << gens[gen].bestModelSnps[0];
				for(unsigned int snp=1; snp < gens[gen].bestModelSnps.size(); snp++){
					os << "_" << gens[gen].bestModelSnps[snp];
				}
				os << " ";
		}
		
		for(int snp=0; snp < totalSnps; snp++){
			os << gens[gen].snpTotals[snp] << " ";
		}  
		os << endl;
	}
}

void NNLog::outputMainHeaders(ostream& os){
			os.precision(3);
 
	// output headers
	os << setw(5) << left << "Gen" << setw(11) << "ValidNN" << setw(10) << "AvgSize" << setw(9) 
		<<  "AvgDepth" << setw(8) << "AvgFit" 
		<< setw(9) << "MaxFit" << setw(9) << "MinFit" 
		<< setw(9) << "AvgGeno" << setw(7) << "AvgCov "
		<< setw(11) << "AvgEpochs" << setw(10) << "MaxEpoch"
		<< setw(15) << "Bestmod";
 
	for(int snp=0; snp < totalSnps; snp++){
		os << snp << " ";
	}
	os << endl;
}

void NNLog::outputSNPHeaders(std::ostream& os){
	// first row is generation number
	for(size_t i=0; i<gens.size(); i++){
		os << "g" << i << " ";
	}
	os << endl;    
}

///
/// output snp sizes for each
/// @param os ostream
///
void NNLog::outputSNPSizes(std::ostream& os, unsigned int totalPopSize){
	// in each row place the number of snps in each model at that ranking
	for(unsigned int currInd=0; currInd < totalPopSize; currInd++){
		for(size_t i=0; i<gens.size(); i++){
			if(currInd >= gens[i].allModels.size()){
				os << "NA ";
				continue;
			}
			os << gens[i].allModels[currInd].size() << " ";
		}
		os << endl;
	}
}

void NNLog::outputFitnessHeaders(std::ostream& os){
	// first row is generation number
	for(size_t i=0; i<gens.size(); i++){
		os << "g" << i << " ";
	}
	os << endl;  
}


///
/// output fitnesses for each
/// @param os ostream
///
void NNLog::outputFitness(std::ostream& os, unsigned int totalPopSize){
	// first row is generation number
	for(size_t i=0; i<gens.size(); i++){
		os << "g" << i << " ";
	}
	os << endl;
	
	// in each row place the number of snps in each model at that ranking
	for(unsigned int currInd=0; currInd < totalPopSize; currInd++){
		for(size_t i=0; i<gens.size(); i++){
			if(currInd >= gens[i].allFitness.size()){
				os << "NA ";
				continue;
			}
			os << gens[i].allFitness[currInd]<< " ";
		}
		os << endl;
	}
}


#ifdef PARALLEL

	///
	/// Send and receive all log information using MPI_Gather.  The log information
	/// is gather by the master.
	///
	void NNLog::sendReceiveLogs(int nprocs, int myRank){
		float* rbuf=NULL;
		int bestModelSize = 256;
		int sendBufferSize = 12 + bestModelSize;
		
		float* sendBuffer = new float[sendBufferSize];
		if(myRank == 0){
				rbuf = new float[nprocs*sendBufferSize];
		}
		
		fill_send_buffer(sendBuffer);
		
		MPI_Gather(sendBuffer, sendBufferSize, MPI_FLOAT, rbuf, sendBufferSize, 
				MPI_FLOAT, 0, MPI_COMM_WORLD);
				
		delete [] sendBuffer;
		
		// receive snp totals
		int sendSize = gens[0].snpTotals.size();
		int totalSize = gens.size() * sendSize;
		int * sendSnps = new int[totalSize];
		float* rbufSnps=NULL;
		
		if(myRank==0){
				rbufSnps = new float[totalSize * nprocs];        
		}
		
		for(unsigned int i=0; i < gens.size(); i++){
			int base = i * sendSize;
			for(unsigned int j=0; j < gens[i].snpTotals.size(); j++){
				sendSnps[base+j] = gens[i].snpTotals[j];
			}
		}
		
		MPI_Gather(sendSnps, totalSize, MPI_INT, rbufSnps, totalSize, MPI_INT,
				0, MPI_COMM_WORLD);
		
		delete [] sendSnps;
		
		if(myRank==0){
				fill_master_log(rbuf, sendBufferSize, rbufSnps, totalSize, nprocs);
				delete [] rbuf;
				delete [] rbufSnps;
		}
		
	}

	void NNLog::fillMasterLog(float* recData, int basicSize, float* snps, int snpSize,
		int nprocs){
		
		int i=0, currIndex;
		
		// add to the master log
		// can skip master info as that is already in the master node
		for(int proc=1; proc < nprocs; proc++){
				currIndex = proc * basicSize;
				gens[i].numValidNN += recData[currIndex++];
				gens[i].numGenotypes += recData[currIndex++];
				gens[i].numCovariates += recData[currIndex++];
				gens[i].totalFitness +=  recData[currIndex++];
				gens[i].totalSize +=  recData[currIndex++];    
				float maxFit = recData[currIndex++];
				float minFit = recData[currIndex++];

				float minEpochs = recData[currIndex++];
				float maxEpochs = recData[currIndex++];
				gens[i].totalEpochs += recData[currIndex++];
				gens[i].totalDepth += recData[currIndex++];
				int modSize = recData[currIndex++];
				vector<int> snps(modSize,0);
				for(int j=0; j < modSize; j++){
					snps[j] = recData[currIndex++];
				}
				
				if(maxFit > gens[i].maxFitness){
					gens[i].maxFitness = maxFit;
					gens[i].bestModelSnps = snps;
				}
				if(minFit < gens[i].minFitness){
					gens[i].minFitness = minFit;
				}
				
				if(minEpochs < gens[i].minEpochs){
					gens[i].minEpochs = minEpochs;
				}
				if(maxEpochs > gens[i].maxEpochs){
					gens[i].maxEpochs = maxEpochs;
				} 
		}

		gens[i].avgScores();        
		
		// add snp totals to logging information
		for(int proc=1; proc < nprocs; proc++){
				currIndex = proc * snpSize;
				for(int snp=0; snp<snpSize; snp++){
						gens[i].snpTotals[snp] += snps[currIndex++];
				}
		}
	}



	void NNLog::fillSendBuffer(float* sendData){ 
		int currIndex=0;
		for(unsigned int i=0; i < gens.size(); i++){
			sendData[currIndex++] = gens[i].numValidNN;
			sendData[currIndex++] = gens[i].numGenotypes;
			sendData[currIndex++] = gens[i].numCovariates;
			sendData[currIndex++] = gens[i].totalFitness;
			sendData[currIndex++] = gens[i].totalSize;
			sendData[currIndex++] = gens[i].maxFitness;
			sendData[currIndex++] = gens[i].minFitness; 
			
			sendData[currIndex++] = gens[i].minEpochs;
			sendData[currIndex++] = gens[i].maxEpochs;
			sendData[currIndex++] = gens[i].totalEpochs;
			sendData[currIndex++] = gens[i].totalDepth;
			
			sendData[currIndex++] = gens[i].bestModelSnps.size();
			for(unsigned int j=0; j<gens[i].bestModelSnps.size(); j++){
				sendData[currIndex++] = gens[i].bestModelSnps[j];
			}
			sendData[currIndex] = -1000;
		}
	}


	///
	/// Send log information back to master
	///
	void NNLog::sendLog(){

	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

		// need to calculate total size for first part
		int totalSize = 0;
		for(unsigned int i=0; i< gens.size(); i++){
			totalSize += 12; // 10 info + 1 for size of best model
			totalSize += gens[i].bestModelSnps.size();
		}   
		// send total size to master
		MPI_Send(&totalSize, 1, MPI_INT, 0, 144, MPI_COMM_WORLD);
		
		// create send array
		float * sendData = new float[totalSize];
		int currIndex = 0;
		for(unsigned int i=0; i < gens.size(); i++){
			sendData[currIndex++] = gens[i].numValidNN;
			sendData[currIndex++] = gens[i].numGenotypes;
			sendData[currIndex++] = gens[i].numCovariates;
			sendData[currIndex++] = gens[i].totalFitness;
			sendData[currIndex++] = gens[i].totalSize;
			sendData[currIndex++] = gens[i].maxFitness;
			sendData[currIndex++] = gens[i].minFitness; 
			
			sendData[currIndex++] = gens[i].minEpochs;
			sendData[currIndex++] = gens[i].maxEpochs;
			sendData[currIndex++] = gens[i].totalEpochs;
			sendData[currIndex++] = gens[i].totalDepth;
			
			sendData[currIndex++] = gens[i].bestModelSnps.size();
			for(unsigned int j=0; j<gens[i].bestModelSnps.size(); j++){
				sendData[currIndex++] = gens[i].bestModelSnps[j];
			}
		}
		
		MPI_Send(sendData, totalSize, MPI_FLOAT, 0, 145, MPI_COMM_WORLD);
		delete [] sendData;
		
		// Send snp index information 
		int sendSize = gens[0].snpTotals.size();
		totalSize = gens.size() * sendSize;
		int * sendSnps = new int[totalSize];
		
		for(unsigned int i=0; i < gens.size(); i++){
			int base = i * sendSize;
			for(unsigned int j=0; j < gens[i].snpTotals.size(); j++){
				sendSnps[base+j] = gens[i].snpTotals[j];
			}
		}
		
		MPI_Send(sendSnps, totalSize, MPI_INT, 0, 146, MPI_COMM_WORLD);
		delete [] sendSnps;
	}



 ///
	/// Receive log information from slaves and store in generations in master.
	/// After receiving all logs recompute the averages.
	/// @param nprocs Number of slave processes to receive
	///
	void NNLog::receiveLogs(int nprocs){ 
		
		int nProcsDone = 1; // one for master

		MPI_Status status;
		vector<int> totalSizes(nprocs, 0);
		int totalSize;

		while(nProcsDone < nprocs){
			MPI_Recv(&totalSize, 1, MPI_INT, MPI_ANY_SOURCE, 144, MPI_COMM_WORLD, &status);
			totalSizes[status.MPI_SOURCE]=totalSize;
			nProcsDone++;
		}
		
		// receive until all slaves have reported log information
		for(int proc=1; proc < nprocs; proc++){
		
			totalSize = totalSizes[proc];
			float * recData = new float[totalSize];
			
			MPI_Recv(recData, totalSize, MPI_FLOAT, proc, 145, MPI_COMM_WORLD, &status);
			
			int currIndex=0, modSize, maxEpochs, minEpochs;
			float maxFit, minFit;
			
			for(unsigned int i=0; i < gens.size(); i++){
				gens[i].numValidNN += recData[currIndex++];
				gens[i].numGenotypes += recData[currIndex++];
				gens[i].numCovariates += recData[currIndex++];
				gens[i].totalFitness +=  recData[currIndex++];
				gens[i].totalSize +=  recData[currIndex++];    
				maxFit = recData[currIndex++];
				minFit = recData[currIndex++];

				minEpochs = recData[currIndex++];
				maxEpochs = recData[currIndex++];
				gens[i].totalEpochs += recData[currIndex++];
				gens[i].totalDepth += recData[currIndex++];
				
				modSize = recData[currIndex++];
				vector<int> snps(modSize,0);
				for(int j=0; j < modSize; j++){
					snps[j] = recData[currIndex++];
				}
				
				if(maxFit > gens[i].maxFitness){
					gens[i].maxFitness = maxFit;
					gens[i].bestModelSnps = snps;
				}
				if(minFit < gens[i].minFitness){
					gens[i].minFitness = minFit;
				}
				
				if(minEpochs < gens[i].minEpochs){
					gens[i].minEpochs = minEpochs;
				}
				if(maxEpochs > gens[i].maxEpochs){
					gens[i].maxEpochs = maxEpochs;
				} 
			}
			delete [] recData;
		}

		for(unsigned int i=0; i < gens.size(); i++){
			gens[i].avgScores();
		}      

		 // Send snp index information 
		int sendSize = gens[0].snpTotals.size();
		totalSize = gens.size() * sendSize;
		int * recSnps = new int[totalSize];
		
		nProcsDone = 1;
		while(nProcsDone < nprocs){
		
			MPI_Recv(recSnps, totalSize, MPI_INT, MPI_ANY_SOURCE, 146, MPI_COMM_WORLD, &status);
		
			for(unsigned int i=0; i < gens.size(); i++){
				int base = i * sendSize;
				for(unsigned int j=0; j < gens[i].snpTotals.size(); j++){
					gens[i].snpTotals[j] += recSnps[base+j];
				}
			}
			nProcsDone++;
		}
		
		delete [] recSnps;
	}



///
/// Sends detailed information to master such as fitness of each model 
/// and the number of snps in each model.
///
void NNLog::sendDetailedLog(){
	int fitnessSizeSignal = 199;
	int fitnessSignal = 200;
	int snpSignal = 202;
	
	float *fitnessSend=NULL;
	int totalModels, totalSnps, *snpSend=NULL; 
	packageDetailed(totalModels, fitnessSend, totalSnps, snpSend);

	int * messageSizes = new int[2];
	messageSizes[0] = totalModels;
	messageSizes[1] = totalSnps;
	MPI_Send(messageSizes, 2, MPI_INT, 0, fitnessSizeSignal, MPI_COMM_WORLD);
	delete [] messageSizes;
	// send to master
	MPI_Send(fitnessSend, totalModels, MPI_FLOAT, 0, fitnessSignal, MPI_COMM_WORLD);
	delete [] fitnessSend;
	// send to master
	MPI_Send(snpSend, totalSnps, MPI_INT, 0, snpSignal, MPI_COMM_WORLD);  
	delete [] snpSend;
}


///
/// Receives detailed information from slaves 
/// @param nprocs Number of processors running
///
void NNLog::receiveDetailedLogs(int nprocs){

	int fitnessSizeSignal = 199;
	int fitnessSignal = 200;
	int snpSignal = 202;
	
	vector<vector<indModel> > indModels;
	
	// need to package master data first
	float *fitnessSend=NULL;
	int totalModels, totalSnps, *snpSend=NULL;  
	packageDetailed(totalModels, fitnessSend, totalSnps, snpSend);
	
	// need to clear master fitness and snp information
	for(size_t i=0; i<gens.size(); i++){
		gens[i].allFitness.clear();
		gens[i].allModels.clear();
	}

	mergeDetailed(totalModels, fitnessSend, totalSnps, snpSend, indModels);
	delete [] fitnessSend;
	delete [] snpSend;

	MPI_Status status;
	
	vector<int> totalFitSizes(nprocs, 0);

	// stores models from all nodes before merging and sorting
	vector<indModel> temp;
	vector<vector<indModel> > models(nprocs, temp);
	
	int * messageSizes = new int[2];
	
	for(int proc=1; proc < nprocs; proc++){  
		MPI_Recv(messageSizes, 2, MPI_INT, proc, fitnessSizeSignal, MPI_COMM_WORLD, &status);
		float * fitRec = new float[messageSizes[0]];
		MPI_Recv(fitRec, messageSizes[0], MPI_FLOAT, proc, fitnessSignal, MPI_COMM_WORLD, &status);
		int * snp_rec = new int[messageSizes[1]];
		MPI_Recv(snp_rec, messageSizes[1], MPI_FLOAT, proc, snpSignal, MPI_COMM_WORLD, &status);
		mergeDetailed(messageSizes[0], fitRec, messageSizes[1], snp_rec, indModels);
		delete [] fitRec;
		delete [] snp_rec;
	}
	
	// need to sort the models by fitness and then transfer into log 
	for(vector<vector<indModel> >::iterator modIter = indModels.begin(); modIter != indModels.end(); 
		modIter++){
		sort(modIter->begin(), modIter->end(), compareIndModels);
	}
	
	for(size_t currGen=0; currGen < gens.size(); currGen++){
		// add sorted information into log
		for(vector<indModel>::iterator modIter=indModels[currGen].begin(); modIter != indModels[currGen].end();
			modIter++){
				gens[currGen].allFitness.push_back(modIter->fitness);
				gens[currGen].allModels.push_back(modIter->snps);
		}
	}
	
}


///
/// Merges information from a single node into vector of models that will
/// be sorted and then inserted back into master log
/// @param totalFitSize total size of fitSend array
/// @param fitRec Array with fitness values
/// @param totalSnpSize total size of snpSend array
/// @param snp_rec Array with snps in models
///
void NNLog::mergeDetailed(int totalFitSize, float* fitRec, int totalSnpSize,
	int* snp_rec, vector<vector<indModel> >& indModels){
	vector<indModel> temp;
	vector<vector<indModel> > tempmodels;
 
	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
 
	float splitGen = 77777;
	// fill with fitness values
	for(int currFit=0; currFit < totalFitSize; currFit++){
		if(fitRec[currFit] == splitGen){
			tempmodels.push_back(temp);
			temp.clear();
		}
		else{
			indModel indtemp;
			indtemp.fitness = fitRec[currFit];
			temp.push_back(indtemp);
		} 
	}
	int currGen=0, currMod=0;
	vector<int> currSnps;

	// fill with snp information
	// it is an array of ints with each terminated on a -1
	// end of a generation is marked with a -2  
	for(int currSnp=0; currSnp < totalSnpSize; currSnp++){
		if(snp_rec[currSnp] == -2){
			currGen++;
			currMod = 0;
		}
		else if(snp_rec[currSnp] == -1){
			tempmodels[currGen][currMod].snps = currSnps;
			currSnps.clear();
			currMod++;
		}
		else{
			currSnps.push_back(snp_rec[currSnp]);
		}
		
	}
	
	temp.clear();
	// add the current models to the overall model
	for(size_t g=0; g<gens.size(); g++){
		if(indModels.size() <= g)
			indModels.push_back(temp);
		indModels[g].insert(indModels[g].end(), tempmodels[g].begin(), tempmodels[g].end()); 
	}
}


///
/// Packages data from a single node for transfer
/// @param totalFitSize Out parameter for total size of fitSend array
/// @param fitSend Array that will be allocated and filled
/// @param totalSnpSize Out parameter for total size of snpSend array
/// @param snpSend Array that will allocated and filled with snps in models
///
void NNLog::packageDetailed(int& totalFitSize, float*& fitSend, int& totalSnpSize,
	int*&snpSend){

int myRank;
MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  

	float splitGen = 77777;
	// determine overall size 
	totalFitSize=0;
	for(size_t i=0; i < gens.size(); i++){
		totalFitSize += int(gens[i].allFitness.size());
	}
	// add space for indicator to mark generation end
	totalFitSize += int(gens.size());
	fitSend = new float[totalFitSize];   
	int currFit=0;
	
	// fill float array to send
	for(size_t i=0; i<gens.size(); i++){  
		for(vector<float>::iterator fiter=gens[i].allFitness.begin(); fiter != gens[i].allFitness.end();
			fiter++){
			fitSend[currFit++] = *fiter;
		}
		fitSend[currFit++] = splitGen;
	}
	totalSnpSize=0;
	
	for(size_t i=0; i<gens.size(); i++){  
		for(size_t j=0; j<gens[i].allModels.size(); j++){
			totalSnpSize += int(gens[i].allModels[j].size());
			totalSnpSize++;
		}
		totalSnpSize++;
	}
	
	snpSend = new int[totalSnpSize];
	
	int currSnp =0;
	
	// prepare snp list for models and send
	// it is an array of ints with each terminated on a -1
	// end of a generation is marked with a -2  
	for(size_t i=0; i<gens.size(); i++){  
		for(size_t j=0; j<gens[i].allModels.size(); j++){
			for(vector<int>::iterator snpIter=gens[i].allModels[j].begin(); snpIter != gens[i].allModels[j].end();
				snpIter++){
				snpSend[currSnp++] = *snpIter;
			}
			snpSend[currSnp++] = -1;
		}
		snpSend[currSnp++] = -2;
	}    
}
	
#endif
