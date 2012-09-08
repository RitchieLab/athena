/* 
 * File:   GENNAlg.h
 * Author: dudeksm
 *
 * Created on November 10, 2008, 2:03 PM
 */

#ifndef _GENNALG_H
#define	_GENNALG_H

#include "Algorithm.h"
#include "InitGEgenome.h"
#include "AlgorithmFactory.h"
#include "NNLog.h"
#include "NNModelLog.h"
#include "GENNGrammarAdjuster.h"
#include "BioFilterModelCollection.h"

#ifdef PARALLEL
#define MAX_GENOME_SIZE 100000
#endif


class GENNAlg:public Algorithm{
    
public:

    GENNAlg();
    
    ~GENNAlg();
    
    /// Set the parameters for the algorithm
    virtual void set_params(AlgorithmParams& alg_params, int numExchanges, int numGenos, int numContin);
    
    /// Run algorithm 
    void run();
    
    /// Runs a step of the algorithm
    virtual void step();
    
    /// Sets random seed
    void setrand(unsigned int seed){randSeed = seed;}
    
    void set_dataset(Dataset* new_set);
    
    /// Sets testing values for best solutions
    void test_solution(Dataset* test_set, int nproc);
    
    /// Outputs individual evaluations to stream
    void output_ind_evals(Dataset* set, ostream& os, int model);
    
    /// Initializes the algorithm for each new dataset to test
    void initialize();
    
    /// Writes output to stream
    void writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
      bool map_used, bool ott_dummy, bool continmap_used);
    
    /// Returns extension for graphical representation
    std::string getGraphicalFileExt();
    
    /// Saves the log
    void saveLog();
    
    /// Starts log
    void startLog(int num_snps);
   
    /// Writes log
    void writeLog();
    
    /// Clears log
    void clearLogs();
    
    /// Close log
    void CloseLog();
    
    /// Finish log and pulls together model information
    void finishLog(std::string basename, int cv);
    
    /// Prepares log files
    void prepareLog(std::string basename, int cv);
    
    /// Returns covariates and snps in best network 
    vector<string> getBestVariables();
    
    /// Retrieves the models from BioFilter and stores the information in the algorithm
    void getBioModels(std::string filename, std::string biofiletype, data_manage::Dataholder* holder);
    
    /// Retrieves the models from the BioFilter archive files and stores information in algorithm
    void getBioModelsArchive(string genegeneFile, string archiveFile, data_manage::Dataholder* holder);
    
    #ifdef PARALLEL
      struct genome_mpi{
        float genomeParams[7];
        int codons[MAX_GENOME_SIZE];
      };
    
      typedef struct genome_mpi struct_mpi;
      
      virtual void setRank(int rank);
      void updateWithMigration(float* stats, int* codons, int totalNodes, int myRank, int max_length=0);
      void updateWithMigration(struct_mpi * genomes, int totalNodes, int myRank);
      void SendAndReceive(int totalNodes, int myRank);
      void SendAndReceiveStruct(int totalNodes, int myRank);
    #endif
    
    void tempoOutputName(string outname){oname = outname;
      templogos = new ofstream;
      templogos->open(oname.c_str(), ios::out);
    }
    
protected:
    
    /// Sets GA for run
    void set_ga_params();
    
    /// Sets default values for parameters
    void initialize_params();   
    
    void fill_population();
    
    void fillLog();
    
    void outputGenome(GAGenome& g);
    
    void convertNetworks(AthenaGrammarSI& currentMapper, AthenaGrammarSI& newMapper);
    
    void set_mapper_prefs(AthenaGrammarSI& hemannMapper);
    
    void ResetCrossover();
    
    void setRestrictedGrammar(bool clearVariables);
    
    void runBackPropagation();
    
    void setBioModels(BioFilterModelCollection& collection, data_manage::Dataholder* holder);
    
    enum GENNParams{
        noMatchParam,
        minSizeParam,
        maxSizeParam,
        tailRatioParam,
        growRateParam,
        maxDepthParam,
        tailSizeParam,
        sensibleInitParam,
        popSizeParam,
        numGenParam,
        probCrossParam,
        probMutParam,
        gramFileParam,
        stepSize,
        calcType,
        useEffectiveXO,
        useAllSnps,
        useAllCovariates,
        requireAll,
        requireAllOnce,
        bioInitFract,
        restrictVarGens,
        bioModelSelection,
        blockCrossGens,
        resetVarsAtMigration,
        bpfreq,
        bpstart,
        gaSelection,
        doubleTournF,
        doubleTournD,
        doubleTournFitFirst,
        prunePlantFract
    };
    
    
    enum BioSelectionType{
      NoMatchSelection,
      rouletteSelect,
      orderedSelect
    };
    
    enum GASelectionType{
      NoMatchSelector,
      RouletteWheelSelection,
      DoubleTournamentSelection
    };
  
    
    void free_memory();
    
    void output_alg_inds();
    
    bool parameters_set;
    
    void AnalyzePopulation();
    std::string oname;
    std::ofstream * templogos;
    
    std::map<std::string, GENNParams> param_map;
    std::map<std::string, BioSelectionType> BioModelSelectionMap;
    std::map<std::string, GASelectionType> GASelectorMap;
    
    // BioFilter parameters
    BioSelectionType biofilter_selector_type;
    
    // GA Selection parameters
    GASelectionType gaSelector;
    bool fitfirst;
    int doubletourneyF;
    float doubletourneyD;
    
    //Random initialization parameters
    unsigned int minSize, maxSize;
    
    //Sensible initialization parameters
    float tailRatio, growRate;
    unsigned int maxDepth, tailSize;
    bool sensibleInit;
    
    //Prune and plant parameters
    float pruneAndPlantFract;
    
    //General algorithm parameters
    bool effectiveXO, useAllVars, useAllCovars, requireAllVars, requireAllVarsOnce, maxbest,
      reset_restricted_at_migration;
    unsigned int wrapEvents, randSeed;
    std::string grammarFile, calculatorName, main_log_filename, fitness_log_filename, 
        snpname_log_filename;
    unsigned int pop_size, num_generations, step_size, ngens_var_restrict, restrict_steps_done,
      ngens_block_cross;
    double prob_cross, prob_mut, init_bio_fract;
    int num_genotypes, num_continuous, bp_first_gen, bp_freq_gen, bp_next_opt;
    
    NNLog* gelog;
    NNModelLog* modellog;
    
    // GE parameters
    AthenaGrammarSI mapper, restrictMapper;
    
    // Genetic algorithm
    GASimpleGA* ga;
    GENNGrammarAdjuster adjuster;
    
    #ifdef PARALLEL
      int genomeInfo, genomeArray;
    #endif
    
};


#endif	/* _GENNALG_H */

