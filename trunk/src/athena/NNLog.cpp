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
void NNLog::output_log(ostream& os){
  
//   os.precision(3);
//  
//   // output headers
//   os << setw(5) << left << "Gen" << setw(11) << "ValidNN" << setw(10) << "AvgSize" << setw(8) 
//     <<  "AvgDepth" << setw(8) << "AvgFit" 
//     << setw(9) << "MaxFit" << setw(9) << "MinFit" 
//     << setw(9) << "AvgGeno" << setw(7) << "AvgCov "
//     << setw(11) << "AvgEpochs" << setw(10) << "MaxEpoch"
//     << setw(15) << "Bestmod";
//  
//   for(int snp=0; snp < total_snps; snp++){
//     os << snp << " ";
//   }
//   os << endl;

  os.precision(3);
  for(unsigned int gen=0; gen < gens.size(); gen++){
    os << setw(5) << left << gens[gen].generationNumber << setw(11) << gens[gen].numValidNN << setw(10) << gens[gen].avgSize 
      << setw(9) << gens[gen].avgDepth << setw(8) << gens[gen].avgFitness 
      << setw(9) << gens[gen].maxFitness << setw(9) << gens[gen].minFitness
      << setw(9) << gens[gen].avgGenos << setw(7) << gens[gen].avgCovars 
      << setw(11) << gens[gen].avgEpochs << setw(10) << gens[gen].maxEpochs <<  " ";
      
    os << gens[gen].best_model_snps[0];
    for(unsigned int snp=1; snp < gens[gen].best_model_snps.size(); snp++){
      os << "_" << gens[gen].best_model_snps[snp];
    }
    os << " ";
    
    for(int snp=0; snp < total_snps; snp++){
      os << gens[gen].snp_totals[snp] << " ";
    }  
    os << endl;
  }
}

void NNLog::output_main_headers(ostream& os){
      os.precision(3);
 
  // output headers
  os << setw(5) << left << "Gen" << setw(11) << "ValidNN" << setw(10) << "AvgSize" << setw(9) 
    <<  "AvgDepth" << setw(8) << "AvgFit" 
    << setw(9) << "MaxFit" << setw(9) << "MinFit" 
    << setw(9) << "AvgGeno" << setw(7) << "AvgCov "
    << setw(11) << "AvgEpochs" << setw(10) << "MaxEpoch"
    << setw(15) << "Bestmod";
 
  for(int snp=0; snp < total_snps; snp++){
    os << snp << " ";
  }
  os << endl;
}

void NNLog::output_snp_headers(std::ostream& os){
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
void NNLog::output_snp_sizes(std::ostream& os, unsigned int totalPopSize){
// first row is generation number
//   for(size_t i=0; i<gens.size(); i++){
//     os << "g" << i << " ";
//   }
//   os << endl;
  
  // in each row place the number of snps in each model at that ranking
  for(unsigned int curr_ind=0; curr_ind < totalPopSize; curr_ind++){
    for(size_t i=0; i<gens.size(); i++){
      if(curr_ind >= gens[i].allmodels.size()){
        os << "NA ";
        continue;
      }
      os << gens[i].allmodels[curr_ind].size() << " ";
    }
    os << endl;
  }
}

void NNLog::output_fitness_headers(std::ostream& os){
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
void NNLog::output_fitness(std::ostream& os, unsigned int totalPopSize){
  // first row is generation number
  for(size_t i=0; i<gens.size(); i++){
    os << "g" << i << " ";
  }
  os << endl;
  
  // in each row place the number of snps in each model at that ranking
  for(unsigned int curr_ind=0; curr_ind < totalPopSize; curr_ind++){
    for(size_t i=0; i<gens.size(); i++){
      if(curr_ind >= gens[i].allfitness.size()){
        os << "NA ";
        continue;
      }
      os << gens[i].allfitness[curr_ind]<< " ";
    }
    os << endl;
  }
}


#ifdef PARALLEL

  ///
  /// Send and receive all log information using MPI_Gather.  The log information
  /// is gather by the master.
  ///
  void NNLog::SendReceiveLogs(int nprocs, int myrank){
    float* rbuf=NULL;
    int best_model_size = 256;
    int send_buffer_size = 12 + best_model_size;
    
    float* send_buffer = new float[send_buffer_size];
    if(myrank == 0){
        rbuf = new float[nprocs*send_buffer_size];
    }
    
    fill_send_buffer(send_buffer);
    
    MPI_Gather(send_buffer, send_buffer_size, MPI_FLOAT, rbuf, send_buffer_size, 
        MPI_FLOAT, 0, MPI_COMM_WORLD);
        
    delete [] send_buffer;
    
    // receive snp totals
    int send_size = gens[0].snp_totals.size();
    int total_size = gens.size() * send_size;
    int * send_snps = new int[total_size];
    float* rbuf_snps=NULL;
    
    if(myrank==0){
        rbuf_snps = new float[total_size * nprocs];        
    }
    
    for(unsigned int i=0; i < gens.size(); i++){
      int base = i * send_size;
      for(unsigned int j=0; j < gens[i].snp_totals.size(); j++){
        send_snps[base+j] = gens[i].snp_totals[j];
      }
    }
    
    MPI_Gather(send_snps, total_size, MPI_INT, rbuf_snps, total_size, MPI_INT,
        0, MPI_COMM_WORLD);
    
    delete [] send_snps;
    
    if(myrank==0){
        fill_master_log(rbuf, send_buffer_size, rbuf_snps, total_size, nprocs);
        delete [] rbuf;
        delete [] rbuf_snps;
    }
    
  }

  void NNLog::fill_master_log(float* rec_data, int basic_size, float* snps, int snp_size,
    int nprocs){
    
    int i=0, curr_index;
    
    // add to the master log
    // can skip master info as that is already in the master node
    for(int proc=1; proc < nprocs; proc++){
        curr_index = proc * basic_size;
//         for(int curr_index=start; curr_index<finish; curr_index++){
        gens[i].numValidNN += rec_data[curr_index++];
        gens[i].numGenotypes += rec_data[curr_index++];
        gens[i].numCovariates += rec_data[curr_index++];
        gens[i].totalFitness +=  rec_data[curr_index++];
        gens[i].totalSize +=  rec_data[curr_index++];    
        float maxfit = rec_data[curr_index++];
        float minfit = rec_data[curr_index++];

        float minepochs = rec_data[curr_index++];
        float maxepochs = rec_data[curr_index++];
        gens[i].totalEpochs += rec_data[curr_index++];
        gens[i].totalDepth += rec_data[curr_index++];
        int modsize = rec_data[curr_index++];
        vector<int> snps(modsize,0);
        for(int j=0; j < modsize; j++){
          snps[j] = rec_data[curr_index++];
        }
        
        if(maxfit > gens[i].maxFitness){
          gens[i].maxFitness = maxfit;
          gens[i].best_model_snps = snps;
        }
        if(minfit < gens[i].minFitness){
          gens[i].minFitness = minfit;
        }
        
        if(minepochs < gens[i].minEpochs){
          gens[i].minEpochs = minepochs;
        }
        if(maxepochs > gens[i].maxEpochs){
          gens[i].maxEpochs = maxepochs;
        } 
//       }
//       delete [] rec_data;
    }

//     for(unsigned int i=0; i < gens.size(); i++){
      gens[i].avgScores();
//     }         
    
    // add snp totals to logging information
    for(int proc=1; proc < nprocs; proc++){
        curr_index = proc * snp_size;
        for(int snp=0; snp<snp_size; snp++){
            gens[i].snp_totals[snp] += snps[curr_index++];
        }
    }
    
  }

  void NNLog::fill_send_buffer(float* send_data){ 
    int curr_index=0;
    for(unsigned int i=0; i < gens.size(); i++){
      send_data[curr_index++] = gens[i].numValidNN;
      send_data[curr_index++] = gens[i].numGenotypes;
      send_data[curr_index++] = gens[i].numCovariates;
      send_data[curr_index++] = gens[i].totalFitness;
      send_data[curr_index++] = gens[i].totalSize;
      send_data[curr_index++] = gens[i].maxFitness;
      send_data[curr_index++] = gens[i].minFitness; 
      
      send_data[curr_index++] = gens[i].minEpochs;
      send_data[curr_index++] = gens[i].maxEpochs;
      send_data[curr_index++] = gens[i].totalEpochs;
      send_data[curr_index++] = gens[i].totalDepth;
      
      send_data[curr_index++] = gens[i].best_model_snps.size();
      for(unsigned int j=0; j<gens[i].best_model_snps.size(); j++){
        send_data[curr_index++] = gens[i].best_model_snps[j];
      }
      send_data[curr_index] = -1000;
    }
    
  }


  ///
  /// Send log information back to master
  ///
  void NNLog::sendLog(){

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // need to calculate total size for first part
    int total_size = 0;
    for(unsigned int i=0; i< gens.size(); i++){
      total_size += 12; // 10 info + 1 for size of best model
      total_size += gens[i].best_model_snps.size();
    }   
    // send total size to master
    MPI_Send(&total_size, 1, MPI_INT, 0, 144, MPI_COMM_WORLD);
    
    // create send array
    float * send_data = new float[total_size];
    int curr_index = 0;
    for(unsigned int i=0; i < gens.size(); i++){
      send_data[curr_index++] = gens[i].numValidNN;
      send_data[curr_index++] = gens[i].numGenotypes;
      send_data[curr_index++] = gens[i].numCovariates;
      send_data[curr_index++] = gens[i].totalFitness;
      send_data[curr_index++] = gens[i].totalSize;
      send_data[curr_index++] = gens[i].maxFitness;
      send_data[curr_index++] = gens[i].minFitness; 
      
      send_data[curr_index++] = gens[i].minEpochs;
      send_data[curr_index++] = gens[i].maxEpochs;
      send_data[curr_index++] = gens[i].totalEpochs;
      send_data[curr_index++] = gens[i].totalDepth;
      
      send_data[curr_index++] = gens[i].best_model_snps.size();
      for(unsigned int j=0; j<gens[i].best_model_snps.size(); j++){
        send_data[curr_index++] = gens[i].best_model_snps[j];
      }
    }
    
    MPI_Send(send_data, total_size, MPI_FLOAT, 0, 145, MPI_COMM_WORLD);
    delete [] send_data;
    
    // Send snp index information 
    int send_size = gens[0].snp_totals.size();
    total_size = gens.size() * send_size;
    int * send_snps = new int[total_size];
    
    for(unsigned int i=0; i < gens.size(); i++){
      int base = i * send_size;
      for(unsigned int j=0; j < gens[i].snp_totals.size(); j++){
        send_snps[base+j] = gens[i].snp_totals[j];
      }
    }
    
    MPI_Send(send_snps, total_size, MPI_INT, 0, 146, MPI_COMM_WORLD);
    delete [] send_snps;
    
  }

 ///
  /// Receive log information from slaves and store in generations in master.
  /// After receiving all logs recompute the averages.
  /// @param nprocs Number of slave processes to receive
  ///
  void NNLog::receiveLogs(int nprocs){ 
    
    int nprocs_done = 1; // one for master

    MPI_Status status;
    vector<int> total_sizes(nprocs, 0);
    int total_size;

    while(nprocs_done < nprocs){
      MPI_Recv(&total_size, 1, MPI_INT, MPI_ANY_SOURCE, 144, MPI_COMM_WORLD, &status);
      total_sizes[status.MPI_SOURCE]=total_size;
      nprocs_done++;
    }
    
    // receive until all slaves have reported log information
    for(int proc=1; proc < nprocs; proc++){
    
      total_size = total_sizes[proc];
      float * rec_data = new float[total_size];
      
      MPI_Recv(rec_data, total_size, MPI_FLOAT, proc, 145, MPI_COMM_WORLD, &status);
      
      int curr_index=0, modsize, maxepochs, minepochs;
      float maxfit, minfit;
      
      for(unsigned int i=0; i < gens.size(); i++){
        gens[i].numValidNN += rec_data[curr_index++];
        gens[i].numGenotypes += rec_data[curr_index++];
        gens[i].numCovariates += rec_data[curr_index++];
        gens[i].totalFitness +=  rec_data[curr_index++];
        gens[i].totalSize +=  rec_data[curr_index++];    
        maxfit = rec_data[curr_index++];
        minfit = rec_data[curr_index++];

        minepochs = rec_data[curr_index++];
        maxepochs = rec_data[curr_index++];
        gens[i].totalEpochs += rec_data[curr_index++];
        gens[i].totalDepth += rec_data[curr_index++];
        
        modsize = rec_data[curr_index++];
        vector<int> snps(modsize,0);
        for(int j=0; j < modsize; j++){
          snps[j] = rec_data[curr_index++];
        }
        
        if(maxfit > gens[i].maxFitness){
          gens[i].maxFitness = maxfit;
          gens[i].best_model_snps = snps;
        }
        if(minfit < gens[i].minFitness){
          gens[i].minFitness = minfit;
        }
        
        if(minepochs < gens[i].minEpochs){
          gens[i].minEpochs = minepochs;
        }
        if(maxepochs > gens[i].maxEpochs){
          gens[i].maxEpochs = maxepochs;
        } 
      }
      delete [] rec_data;
    }

    for(unsigned int i=0; i < gens.size(); i++){
      gens[i].avgScores();
    }      

     // Send snp index information 
    int send_size = gens[0].snp_totals.size();
    total_size = gens.size() * send_size;
    int * rec_snps = new int[total_size];
    
    nprocs_done = 1;
    while(nprocs_done < nprocs){
    
      MPI_Recv(rec_snps, total_size, MPI_INT, MPI_ANY_SOURCE, 146, MPI_COMM_WORLD, &status);
    
      for(unsigned int i=0; i < gens.size(); i++){
        int base = i * send_size;
        for(unsigned int j=0; j < gens[i].snp_totals.size(); j++){
          gens[i].snp_totals[j] += rec_snps[base+j];
        }
      }
      nprocs_done++;
    }
    
    delete [] rec_snps;
  }

///
/// Sends detailed information to master such as fitness of each model 
/// and the number of snps in each model.
///
void NNLog::sendDetailedLog(){
  int fitness_size_signal = 199;
  int fitness_signal = 200;
  int snp_signal = 202;
  
  float *fitness_send=NULL;
  int total_models, total_snps, *snp_send=NULL; 
  packageDetailed(total_models, fitness_send, total_snps, snp_send);

  int * message_sizes = new int[2];
  message_sizes[0] = total_models;
  message_sizes[1] = total_snps;
  MPI_Send(message_sizes, 2, MPI_INT, 0, fitness_size_signal, MPI_COMM_WORLD);
  delete [] message_sizes;
  // send to master
  MPI_Send(fitness_send, total_models, MPI_FLOAT, 0, fitness_signal, MPI_COMM_WORLD);
  delete [] fitness_send;
  // send to master
  MPI_Send(snp_send, total_snps, MPI_INT, 0, snp_signal, MPI_COMM_WORLD);  
  delete [] snp_send;
}


///
/// Receives detailed information from slaves 
/// @param nprocs Number of processors running
///
void NNLog::receiveDetailedLogs(int nprocs){

  int fitness_size_signal = 199;
  int fitness_signal = 200;
  int snp_signal = 202;
  
  vector<vector<indModel> > indmodels;
  
  // need to package master data first
  float *fitness_send=NULL;
  int total_models, total_snps, *snp_send=NULL;  
  packageDetailed(total_models, fitness_send, total_snps, snp_send);
  
  // need to clear master fitness and snp information
  for(size_t i=0; i<gens.size(); i++){
    gens[i].allfitness.clear();
    gens[i].allmodels.clear();
  }

  mergeDetailed(total_models, fitness_send, total_snps, snp_send, indmodels);
  delete [] fitness_send;
  delete [] snp_send;

  MPI_Status status;
  
  vector<int> total_fit_sizes(nprocs, 0);

  // stores models from all nodes before merging and sorting
  vector<indModel> temp;
  vector<vector<indModel> > models(nprocs, temp);
  
  int * message_sizes = new int[2];
  
  for(int proc=1; proc < nprocs; proc++){  
    MPI_Recv(message_sizes, 2, MPI_INT, proc, fitness_size_signal, MPI_COMM_WORLD, &status);
    float * fit_rec = new float[message_sizes[0]];
    MPI_Recv(fit_rec, message_sizes[0], MPI_FLOAT, proc, fitness_signal, MPI_COMM_WORLD, &status);
    int * snp_rec = new int[message_sizes[1]];
    MPI_Recv(snp_rec, message_sizes[1], MPI_FLOAT, proc, snp_signal, MPI_COMM_WORLD, &status);
    mergeDetailed(message_sizes[0], fit_rec, message_sizes[1], snp_rec, indmodels);
    delete [] fit_rec;
    delete [] snp_rec;
  }
  
  // need to sort the models by fitness and then transfer into log 
  for(vector<vector<indModel> >::iterator moditer = indmodels.begin(); moditer != indmodels.end(); 
    moditer++){
    sort(moditer->begin(), moditer->end(), compareIndModels);
  }
  
  for(size_t curr_gen=0; curr_gen < gens.size(); curr_gen++){
    // add sorted information into log
    for(vector<indModel>::iterator moditer=indmodels[curr_gen].begin(); moditer != indmodels[curr_gen].end();
      moditer++){
        gens[curr_gen].allfitness.push_back(moditer->fitness);
        gens[curr_gen].allmodels.push_back(moditer->snps);
    }
  }
  
}


///
/// Merges information from a single node into vector of models that will
/// be sorted and then inserted back into master log
/// @param total_fit_size total size of fit_send array
/// @param fit_rec Array with fitness values
/// @param total_snp_size total size of snp_send array
/// @param snp_rec Array with snps in models
///
void NNLog::mergeDetailed(int total_fit_size, float* fit_rec, int total_snp_size,
  int* snp_rec, vector<vector<indModel> >& indmodels){
  vector<indModel> temp;
  vector<vector<indModel> > tempmodels;
 
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
 
  float split_gen = 77777;
  // fill with fitness values
  for(int curr_fit=0; curr_fit < total_fit_size; curr_fit++){
    if(fit_rec[curr_fit] == split_gen){
      tempmodels.push_back(temp);
      temp.clear();
    }
    else{
      indModel indtemp;
      indtemp.fitness = fit_rec[curr_fit];
      temp.push_back(indtemp);
    } 
  }
  int curr_gen=0, curr_mod=0;
  vector<int> currsnps;

  // fill with snp information
  // it is an array of ints with each terminated on a -1
  // end of a generation is marked with a -2  
  for(int curr_snp=0; curr_snp < total_snp_size; curr_snp++){
    if(snp_rec[curr_snp] == -2){
      curr_gen++;
      curr_mod = 0;
    }
    else if(snp_rec[curr_snp] == -1){
      tempmodels[curr_gen][curr_mod].snps = currsnps;
      currsnps.clear();
      curr_mod++;
    }
    else{
      currsnps.push_back(snp_rec[curr_snp]);
    }
    
  }
  
  temp.clear();
  // add the current models to the overall model
  for(size_t g=0; g<gens.size(); g++){
    if(indmodels.size() <= g)
      indmodels.push_back(temp);
    indmodels[g].insert(indmodels[g].end(), tempmodels[g].begin(), tempmodels[g].end()); 
  }
}


///
/// Packages data from a single node for transfer
/// @param total_fit_size Out parameter for total size of fit_send array
/// @param fit_send Array that will be allocated and filled
/// @param total_snp_size Out parameter for total size of snp_send array
/// @param snp_send Array that will allocated and filled with snps in models
///
void NNLog::packageDetailed(int& total_fit_size, float*& fit_send, int& total_snp_size,
  int*&snp_send){

int myrank;
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  

  float split_gen = 77777;
  // determine overall size 
  total_fit_size=0;
  for(size_t i=0; i < gens.size(); i++){
    total_fit_size += int(gens[i].allfitness.size());
  }
  // add space for indicator to mark generation end
  total_fit_size += int(gens.size());
  fit_send = new float[total_fit_size];   
  int curr_fit=0;
  
  // fill float array to send
  for(size_t i=0; i<gens.size(); i++){  
    for(vector<float>::iterator fiter=gens[i].allfitness.begin(); fiter != gens[i].allfitness.end();
      fiter++){
      fit_send[curr_fit++] = *fiter;
    }
    fit_send[curr_fit++] = split_gen;
  }
  total_snp_size=0;
  
  for(size_t i=0; i<gens.size(); i++){  
    for(size_t j=0; j<gens[i].allmodels.size(); j++){
      total_snp_size += int(gens[i].allmodels[j].size());
      total_snp_size++;
    }
    total_snp_size++;
  }
  
  snp_send = new int[total_snp_size];
  
  int curr_snp =0;
  
  // prepare snp list for models and send
  // it is an array of ints with each terminated on a -1
  // end of a generation is marked with a -2  
  for(size_t i=0; i<gens.size(); i++){  
    for(size_t j=0; j<gens[i].allmodels.size(); j++){
      for(vector<int>::iterator snpiter=gens[i].allmodels[j].begin(); snpiter != gens[i].allmodels[j].end();
        snpiter++){
        snp_send[curr_snp++] = *snpiter;
      }
      snp_send[curr_snp++] = -1;
    }
    snp_send[curr_snp++] = -2;
  }    
}
 
  
 
  
#endif
