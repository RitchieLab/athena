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
 * File:   athena.cpp
 * Author: dudeksm
 *
 * Created on November 10, 2008, 2:10 PM
 */

//#include <stdlib.h>

#include "ConfigFileReader.h"
#include <Dataholder.h>
#include "MDRFileHandler.h"
#include "ContinFileReader.h"
#include "MapFileReader.h"
#include "ContinMapFileReader.h"
#include "CrossValidator.h"
#include "OttDummyConvert.h"
#include "AlgorithmFactory.h"
#include "OutputManager.h"
#include <iostream>
#include <sstream>
#include "OutputSet.h"
#include <ScaledDataFactory.h>
#include <StephenDummyConvert.h>
#include <time.h>
#include <random_func.h>

#ifdef PARALLEL
#include "TransferData.h"
#endif

void exit_app(AthenaExcept& he, int myrank);
void exit_app(DataExcept& de, int myrank);
void adjust_seed(Config& config, int orig_seed, int cv, int nproc, int myrank);
std::string time_diff(double dif);

int main(int argc, char** argv) {

  int nproc = 1; // only one processor when not running in parallel

 bool mapfile_used = false, continmap_used = false;
 int myrank = 0;
   
#ifdef PARALLEL
 

  // set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
#endif /* end PARALLEL code block */
   
    string version_date = "10/12/12";
    string exec_name = "ATHENA";
     time_t start,end;
     
    if(argc < 2){
        AthenaExcept he("\n\tATHENA\n\t" + version_date + 
            "\n\n\tUsage: ATHENA <config>\n\n");
        exit_app(he, myrank);
    }
    else{
#ifdef PARALLEL
  if(myrank==0){
#endif
        time (&start);
        cout << endl << "\t" << exec_name << ":\t" << version_date << endl << endl;
#ifdef PARALLEL
        }
#endif
    }
    
    string configfile = argv[1];
    ConfigFileReader configread;
    Config config;
    ScaleData* scaler = NULL;

    // read config file
    try{
      config = configread.read_config(configfile);
    }catch(AthenaExcept he){
        exit_app(he, myrank);
    }
    
    // fill dataholder with data
    data_manage::Dataholder data;

    try{
      // read in genotype data
        data_manage::MDRFileHandler mdr_reader;
        if(config.getTrainFile().size() ==0)
          mdr_reader.parse_file(config.getDataSetName(), &data, config.getMissingValue(),
               config.getStatusMissingValue(),config.getIDinData());
        else
          mdr_reader.parse_file(config.getTrainFile(), config.getTestFile(), 
            &data, config.getMissingValue(), config.getStatusMissingValue(), 
            config.getIDinData());
        
        // read in continuous data if any
        if(config.getContinTrainFile().size() > 0){
          data_manage::ContinFileReader contin_reader;
          contin_reader.read_contin_file(config.getContinTrainFile(), config.getContinTestFile(),
            &data, config.getContinMiss(), config.getIDinData());          
        }
        else if(config.getContinFileName().size() > 0){
            data_manage::ContinFileReader contin_reader;
            contin_reader.read_contin_file(config.getContinFileName(), &data,
                    config.getContinMiss(), config.getIDinData());
        }
        
        // if present read map file
        if(config.getMapName().size() > 0){
            mapfile_used = true;
            data_manage::MapFileReader map_reader;
            map_reader.parse_map_file(config.getMapName(), &data);
        }
        else{
            mapfile_used = false;
            data.add_default_snps();
        }
        
        if(config.getContinMapName().size() > 0){
            continmap_used = true;
            data_manage::ContinMapFileReader continmap_reader;
            continmap_reader.parse_map_file(config.getContinMapName(), &data);
        }
        else{
            continmap_used = false;
            data.add_default_covars();
        }
       
        // convert data to ott dummy representation if needed
        if(config.getOttEncoded()){
          data_manage::OttDummyConvert ott;
          ott.convert_genotypes(&data);
        }
        else if(config.getEncodeType() == Config::StephenDummy){
          data_manage::StephenDummyConvert st;
          st.convert_genotypes(&data);
        }

        // alter continuous variables and status value
        scaler = data_manage::ScaledDataFactory::create_scaler(config.getStatusAdjust());
        scaler->adjust_status(&data);
        for(unsigned int c=0; c < data.num_covariates(); c++){
          scaler->adjust_contin(&data, c);
        }

    }catch(AthenaExcept he){
        exit_app(he, myrank);
    }
 
 
 	// restart code needs to be placed here
 	// establish original splits and load in the saved states of the random 
 	// number generators
 	
    // set random seed  before splitting
    data_manage::rand_seed(config.getRandSeed());

    // construct crossvalidation sets to use in running algorithm
    CrossValidator cvmaker;
    CVSet cv_set;

	if(config.getSplitFile()==""){
	    cv_set = cvmaker.split_data(config.getNumCV(), &data);
#ifdef PARALLEL
		if(myrank==0)
#endif
	    cvmaker.save_splits(config.getOutputName() + ".cvsplit");
	}
	else{
	    try{
//       config = configread.read_config(configfile);
      cv_set = cvmaker.load_splits(config.getSplitFile(), &data);
    }catch(DataExcept de){
        exit_app(de, myrank);
    }
		cv_set = cvmaker.load_splits(config.getSplitFile(), &data);
// cvmaker.save_splits(config.getOutputName() + ".cvsplit2");
	}

    // run crossvalidations and store the populations
    int num_cv = cv_set.num_intervals();
    
// #ifdef PARALLEL
//   // when running in parallel need to change the seeds for slaves
//   // after splitting data so that all algorithms will not be running same seeds
//   if(myrank != 0){
//     int newseed = config.getRandSeed() + myrank * 25;
//     if(newseed > RAND_MAX)
//       newseed -= RAND_MAX;
//     config.setRandSeed(newseed);
//   }
// #endif /* end parallel seed adjustment for slaves */
   
//    	adjust_seed(Configuration& config, int cv);
   
    // create algorithm
    vector<AlgorithmParams> alg_params= config.getAlgorithmParams();
    Algorithm* alg = AlgorithmFactory::create_algorithm(alg_params[0].name);
#ifdef PARALLEL
    alg->setRank(myrank);
    alg->setTotalNodes(nproc);
#endif
    alg->setrand(config.getRandSeed());
    alg->set_dummy_encoding(config.getOttEncoded());
    alg->setLogType(config.getLogType());
    
    alg->set_params(alg_params[0], config.getNumExchanges(), 
      data.num_genos(), data.num_covariates());

    /// store results in a population vector
    vector<Population> pops;
    string cvfilename = "cv";

    OutputManager writer;
    writer.setBasename(config.getOutputName());
    int curr_cv=config.getStartCV()-1;
	if(curr_cv==0)
		writer.setFiles(mapfile_used, alg->get_fitness_name());
    int original_seed = config.getRandSeed();
    for(; curr_cv < num_cv; curr_cv++){
    	adjust_seed(config, original_seed,curr_cv, nproc, myrank);
    	alg->setrand(config.getRandSeed());
#ifdef PARALLEL
  if(myrank==0){
#endif
      cout << "Beginning Cross-validation " << curr_cv + 1 << "...";
      cout.flush();
#ifdef PARALLEL
  }
#endif

      alg->set_dataset(&(cv_set.get_interval(curr_cv).get_training()));

#ifdef PARALLEL
    if(myrank==0){  // only have master output cv when desired
#endif
      // output cv interval when requested
      if(config.getCVOutput()){
        OutputSet oset;
        oset.outputCV(cvfilename, cv_set.get_interval(curr_cv).get_training(), curr_cv+1);
      }
#ifdef PARALLEL
  } /* end check for master writing CV splits */
#endif
      alg->startLog(data.num_genos());

    if(config.getBioFilterFile().size() > 0){
      alg->getBioModels(config.getBioFilterFile(), config.getBioFileType(), &data);
    }
    else if(config.getBioGeneFile().size() > 0){
      alg->getBioModelsArchive(config.getBioGeneFile(), config.getBioArchiveFile(), &data);
    }
      // prepare logs -- these are written as job progresses
      alg->prepareLog(config.getOutputName(), curr_cv+1);
      try{
          alg->initialize();
    }catch(AthenaExcept& ae){
        exit_app(ae, myrank);
    }
      
      for(int step=0; step < config.getNumExchanges(); step++)
         alg->step();

		alg->CloseLog();
      // prepare and write logs -- algorithm will output as much information as requested
//       alg->writeLog(config.getOutputName(), curr_cv+1);

      if(num_cv > 1)
        alg->test_solution(&(cv_set.get_interval(curr_cv).get_testing()), nproc);
      
      // check population values
//       pops.push_back(alg->getPopulation());
      Population pop = alg->getPopulation();
	
      int curr_proc = 0;
#ifdef PARALLEL
  if(myrank==0){
    int lastproc = 1;
    if(config.outputAllNodesBest())
      lastproc = nproc;
    for(curr_proc=0; curr_proc < lastproc; curr_proc++){
#endif
      if(config.getIndOutput()){
        stringstream ss;
        ss << config.getOutputName() << "." << curr_cv+1 << "." << curr_proc+1 << ".ind_results.txt";
        ostream & os = writer.getStream(ss.str());
        os << "Ind ID\tPredicted\tOriginal\n";
        alg->output_ind_evals(&(cv_set.get_interval(curr_cv).get_training()), os, curr_proc);
        if(num_cv > 1)
          alg->output_ind_evals(&(cv_set.get_interval(curr_cv).get_testing()), os, curr_proc);
        writer.closeStream();
      }
#ifdef PARALLEL
  }
#endif 
    cout << " Completed" << endl;
    alg->finishLog(config.getOutputName(),curr_cv+1);
    int nmodels=1;
#ifdef PARALLEL
    if(myrank==0){

      if(config.outputAllNodesBest())
        nmodels = nproc;
#endif
	// update output when needed
    if(pop.getConvertScores()){
      if(num_cv > 1)
 //        for(int curr_cv=0; curr_cv < num_cv; curr_cv++){ 
          pop.convert_scores(&(cv_set.get_interval(curr_cv).get_training()), 
            &(cv_set.get_interval(curr_cv).get_testing()));
//         }
      else
        pop.convert_scores(&(cv_set.get_interval(0).get_training()));
    }
    writer.outputSummary(pop, curr_cv, data, mapfile_used, config.getOttEncoded(), continmap_used,
        alg->get_fitness_name());
        
    switch(config.getSummaryOnly()){
      case Config::False:
        writer.outputGraphic(alg, pop, curr_cv, config.getOutputName(), nmodels, data, 
         mapfile_used, config.getOttEncoded(), continmap_used);
      case Config::Best:
        writer.outputBestModels(pop, nmodels, curr_cv,scaler->output_scale_info(), data, 
          mapfile_used, config.getOttEncoded(), continmap_used); 
      default:
        ;
    }
    
#ifdef PARALLEL
} /* end of output */
#endif
		
    }
//     int nmodels = 1;
//     for(int curr_cv=1; curr_cv <= num_cv; curr_cv++){
//         alg->finishLog(config.getOutputName(),curr_cv);
//     }

// #ifdef PARALLEL
//     if(myrank==0){
// 
//       if(config.outputAllNodesBest())
//         nmodels = nproc;
// #endif

    // update output when needed
//     if(pops[0].getConvertScores()){
//       if(num_cv > 1)
//         for(int curr_cv=0; curr_cv < num_cv; curr_cv++){ 
//           pops[curr_cv].convert_scores(&(cv_set.get_interval(curr_cv).get_training()), 
//             &(cv_set.get_interval(curr_cv).get_testing()));
//         }
//       else
//         pops[0].convert_scores(&(cv_set.get_interval(0).get_training()));
//     }

    // Output results
//     writer.setBasename(config.getOutputName());
    
//     writer.outputSummary(pops, data, mapfile_used, config.getOttEncoded(), continmap_used,
//         alg->get_fitness_name());
    
//     switch(config.getSummaryOnly()){
//       case Config::False:
//         writer.outputGraphic(alg, pops, config.getOutputName(), nmodels, data, 
//          mapfile_used, config.getOttEncoded(), continmap_used);
//       case Config::Best:
//         writer.outputBestModels(pops, nmodels, scaler->output_scale_info(), data, 
//           mapfile_used, config.getOttEncoded(), continmap_used); 
//       default:
//         ;
//     }
//     cout << endl;
    
#ifdef PARALLEL
    }   /* ends master processing of output */
    MPI_Finalize();
 #endif
    delete scaler;

#ifdef PARALLEL
    if(myrank==0){
#endif
    time(&end);
    double dif = difftime (end,start);
    cout << "\n\tAnalysis took " << time_diff(dif) << endl << endl;
#ifdef PARALLEL
    }
#endif
    

    return (EXIT_SUCCESS);
}

///
/// Returns string formatted to indicate passage of time
/// @param dif time in seconds
///
std::string time_diff(double dif){
    double elapsed;
    string period;
    if(dif < 60){
        elapsed = dif;
        period = " seconds";
    }
    else if(dif < 3600){
        elapsed = dif/60;
        period = " minutes";
    }
    else{
        elapsed = dif/3600;
        period = " hours";
    }
    return Stringmanip::itos(elapsed) + period;
}

///
/// Outputs message in exception and exits program
/// @param he AthenaExcept
///
void exit_app(AthenaExcept& he, int myrank){
    if(myrank==0)
      cout << he.what() << endl << endl;;
#ifdef PARALLEL
  MPI_Finalize();
#endif
  exit(EXIT_FAILURE);    
} 


///
/// Outputs message in exception and exits program
/// @param he AthenaExcept
///
void exit_app(DataExcept& de, int myrank){
    if(myrank==0)
      cout << de.what() << endl << endl;;
#ifdef PARALLEL
  MPI_Finalize();
#endif
  exit(EXIT_FAILURE);    
} 

///
/// Adjust seed based on current CV and number of processors in run
/// @param config Configuration
/// @param cv Current cross-validation
/// @param nproc Total number of processors
/// @param myrank Rank of this process
///
void adjust_seed(Config& config, int orig_seed, int cv, int nproc, int myrank){
  int newseed = orig_seed + nproc * cv + myrank;
  config.setRandSeed(newseed);
}

