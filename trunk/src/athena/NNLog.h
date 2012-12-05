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
  
    NNLog(int num_snps):AlgorithmLog(num_snps) {
      gen_number=0;
      genTotals newgen(num_snps, gen_number);
      gens.push_back(newgen);
      gen_index = 0;
      total_snps = num_snps;
      maxbest = true;
    }
  
    /// Outputs data stored in this log
    void output_log(std::ostream& os);

    /// Sets max as best (true) or false for min score is best
    void set_max_best(bool val){maxbest = val;}

    /// Adds a generation worth of data information
    inline void add_generation(){
//       genTotals newgen(total_snps);
//       gens.push_back(newgen);
//       gen_index++;
        // change to re-use generation 
        // will output after every generation
        gens.clear();
        gen_number++;
        genTotals newgen(total_snps,gen_number);
        gens.push_back(newgen);
        gen_index=0;
    }


    inline unsigned int get_current_gen(){return gen_number;}

    inline void complete_gen(){
      gens[gen_index].avgScores();
    }

    inline void add_network(){gens[gen_index].numValidNN++;}

    inline void add_fitness(float fitness, std::vector<int> snps){
      gens[gen_index].totalFitness += fitness;
      if(fitness > gens[gen_index].maxFitness){
        gens[gen_index].maxFitness = fitness;
        if(maxbest)
          gens[gen_index].best_model_snps = snps;
      }
      if(fitness < gens[gen_index].minFitness){
        gens[gen_index].minFitness = fitness;
        if(!maxbest)
          gens[gen_index].best_model_snps = snps;
      }
      add_snps(snps);
      gens[gen_index].allfitness.push_back(fitness);
    }
    
    inline void add_fitness(float fitness){
      gens[gen_index].totalFitness += fitness;
      if(fitness > gens[gen_index].maxFitness)
        gens[gen_index].maxFitness = fitness;
      if(fitness < gens[gen_index].minFitness)
        gens[gen_index].minFitness = fitness;
      gens[gen_index].allfitness.push_back(fitness);
    }

    inline void add_epochs(int nepochs){
      if(nepochs < gens[gen_index].minEpochs)
        gens[gen_index].minEpochs = nepochs;
      if(nepochs > gens[gen_index].maxEpochs)
        gens[gen_index].maxEpochs = nepochs;
      gens[gen_index].totalEpochs += nepochs;
    }
    
    inline void add_num_genos(int nGenos){
      gens[gen_index].numGenotypes += nGenos;
    }
    
    inline void add_num_covars(int nCovars){
      gens[gen_index].numCovariates += nCovars;
    }
    
    inline void add_nn_size(int size){
      gens[gen_index].totalSize += size;
    }
    
    inline void add_nn_depth(int depth){
      gens[gen_index].totalDepth += depth;
    }
    
    /// add the snps present in the network
    inline void add_snps(std::vector<int> snps){
      std::vector<int>::iterator iter;
      for(iter = snps.begin(); iter != snps.end(); iter++){
        gens[gen_index].snp_totals[*iter]++;
      }
      gens[gen_index].allmodels.push_back(snps);
    }
    
    inline void set_best_mod(std::vector<int> snps){
      gens[gen_index].best_model_snps = snps;
    }

    /// output snp sizes for each
    void output_snp_sizes(std::ostream& os, unsigned int totalPopSize);
    
    /// output snp sizes headers
    void output_snp_headers(std::ostream& os);
    
    /// output fitnesses for all models in each generation
    void output_fitness(std::ostream& os, unsigned int totalPopSize);
    
    void output_fitness_headers(std::ostream& os);
    
    void output_main_headers(std::ostream& os);

    #ifdef PARALLEL
      void sendLog(); // for slaves
      void receiveLogs(int nprocs); // for master
      void sendDetailedLog(); // for slaves
      void receiveDetailedLogs(int nprocs); // for master
      void SendReceiveLogs(int nprocs, int myrank);
    #endif

  private:
    struct genTotals{   
      genTotals(int num_snps, int gen_num){
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
        snp_totals.assign(num_snps,0);
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
      std::vector<int> snp_totals, best_model_snps;
      std::vector<std::vector<int> > allmodels;
      std::vector<float> allfitness;
    };
       
    #ifdef PARALLEL
      void fill_master_log(float* rec_data, int basic_size, float* snps, int snp_size,
        int nprocs);
      void fill_send_buffer(float* send_data);
      void packageDetailed(int& total_fit_size, float*& fit_send, int& total_snp_size, int*& snp_send);
      void mergeDetailed(int total_fit_size, float* fit_rec, int total_snp_size, int* snp_rec, 
        std::vector<std::vector<indModel> >& indmodels);
    #endif

    std::vector<genTotals> gens;
    unsigned int gen_index;
    int total_snps, gen_number;
    bool maxbest;
};

#endif
