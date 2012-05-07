
#include "GEObjective.h"
#include "GENNAlg.h"
#include "GEObjective.h"
#include "Terminals.h"
#include "GENNGrammarAdjuster.h"
#include <ga/ga.h>
#include <ctime>
#include <set>

///
/// Constructor
///
GENNAlg::GENNAlg(){
    initialize_params();
}

///
/// Destructor
///
GENNAlg::~GENNAlg(){
    if(gelog != NULL){
        delete gelog;
        gelog=NULL;
    }
    free_memory();
}


///
/// Frees memory
///
void GENNAlg::free_memory(){
    if(ga != NULL){
        delete ga;
        ga = NULL;
    }
}

///
/// Set current Dataset for running algorithm
/// @param new_set Dataset
/// 
void GENNAlg::set_dataset(Dataset* new_set){
   set = new_set;
   GEObjective::setDataset(new_set);
}


///
/// Sets parameters within the algorithm for the 
/// analysis.
/// @param alg_param AlgorithmParams
/// @param numExchanges total number of times the best model will be passed to other algorithms
/// when more than one algorithm
/// @param numGenos number of Genotypes in set
/// @param numContin number of Continuous variables in set
/// 
void GENNAlg::set_params(AlgorithmParams& alg_param, int numExchanges, int numGenos, int numContin){
    
    map<string, string>::iterator mapIter;
    
    for(mapIter = alg_param.params.begin(); mapIter != alg_param.params.end(); 
      mapIter++){     
        switch(param_map[mapIter->first]){
            case noMatchParam:
                throw AthenaExcept("No match for parameter " + mapIter->first +
                        "in Algorithm GENN");
                break;
            case minSizeParam:
                minSize = Stringmanip::stouint(mapIter->second);
                break;
            case maxSizeParam:
                maxSize = Stringmanip::stouint(mapIter->second);
                #ifdef PARALLEL
                  if(maxSize > MAX_GENOME_SIZE)
                    throw AthenaExcept(mapIter->first + " can be no more than " + 
                      Stringmanip::itos(MAX_GENOME_SIZE));
                #endif
                break;
            case tailRatioParam:
                tailRatio = Stringmanip::stodouble(mapIter->second);
                break;
            case growRateParam:
                growRate = Stringmanip::stodouble(mapIter->second);
                break;
            case maxDepthParam:
                maxDepth = Stringmanip::stouint(mapIter->second);
                break;
            case tailSizeParam:
                tailSize = Stringmanip::stouint(mapIter->second);
                break;
            case sensibleInitParam:
                sensibleInit = Stringmanip::check_true_false(mapIter->second);
                break;
            case popSizeParam:
                pop_size = Stringmanip::stouint(mapIter->second);
                break;
            case probCrossParam:
                prob_cross = Stringmanip::stodouble(mapIter->second);
                break;
            case probMutParam:
                prob_mut = Stringmanip::stodouble(mapIter->second);
                break;
            case gramFileParam:
                grammarFile = mapIter->second;
                break;
            case stepSize:
                step_size = Stringmanip::stouint(mapIter->second);
                break;
            case calcType:
                calculatorName = Stringmanip::to_upper(mapIter->second);
                break;        
            case useEffectiveXO:
                effectiveXO = Stringmanip::check_true_false(mapIter->second);
                break;
            case useAllSnps:
                useAllVars = Stringmanip::check_true_false(mapIter->second);
                break;
            case useAllCovariates:
                useAllCovars = Stringmanip::check_true_false(mapIter->second);
                break;
            case requireAll:
                requireAllVars = Stringmanip::check_true_false(mapIter->second);
                break;
            case requireAllOnce:
                requireAllVarsOnce = Stringmanip::check_true_false(mapIter->second);
                break;
            case bioInitFract:
                init_bio_fract = Stringmanip::stodouble(mapIter->second);     
                break;
            case restrictVarGens:
                ngens_var_restrict = Stringmanip::stouint(mapIter->second);
                break;
            case bioModelSelection:
                if(BioModelSelectionMap.find(Stringmanip::to_upper(mapIter->second)) == 
                    BioModelSelectionMap.end())
                  throw AthenaExcept("No match for bio model selection type " + mapIter->second);
                else
                  biofilter_selector_type = BioModelSelectionMap[Stringmanip::to_upper(mapIter->second)];
                break;
            case gaSelection:
                if(GASelectorMap.find(Stringmanip::to_upper(mapIter->second)) ==
                  GASelectorMap.end())
                  throw AthenaExcept("No match for GA selection type " + mapIter->second);
                else
                  gaSelector = GASelectorMap[Stringmanip::to_upper(mapIter->second)];
                  break;
#ifdef ATHENA_BLOAT_CONTROL
            case doubleTournF:
                doubletourneyF = Stringmanip::stouint(mapIter->second);
                break;
            case doubleTournD:
                doubletourneyD = Stringmanip::stodouble(mapIter->second);
                break;
            case doubleTournFitFirst:
                fitfirst = Stringmanip::check_true_false(mapIter->second);
                break;
#endif
            case blockCrossGens:
                ngens_block_cross = Stringmanip::stouint(mapIter->second);
                break;
            case resetVarsAtMigration:
                reset_restricted_at_migration = Stringmanip::check_true_false(mapIter->second);
                break;
            case bpstart:
                bp_first_gen = Stringmanip::stoi(mapIter->second);
                break;
            case bpfreq:    
                bp_freq_gen = Stringmanip::stoi(mapIter->second);
                break;
#ifdef ATHENA_BLOAT_CONTROL
            case prunePlantFract:
#endif
            default:
                throw AthenaExcept("No match for parameter " + mapIter->first +
                        " in Algorithm GENN");               
        }
    }
    
    if(step_size > num_generations){
        step_size = num_generations;
    }
    
    num_generations = step_size * numExchanges;
    num_genotypes = numGenos;
    num_continuous = numContin;
    
    // first optimization of backpropagation 
    bp_next_opt = bp_first_gen;
       
    set_ga_params();
    
}


///
/// Initializes parameters for basic run of algorithm.  These
/// can be modified using the set_params function
///
void GENNAlg::initialize_params(){
   
   myRank = 0;
   minSize = 50; // Minimum size for Random Initialization
   maxSize = 200; // Maximum size for Random Initialization

   tailRatio = 0.0;
   growRate = 0.5;

   sensibleInit = false;
   maxDepth = 10;
   tailSize = 0;

   grammarFile = "";
   wrapEvents = 0;

   pop_size = 100;
   num_generations = 100;
   prob_cross = 0.9;
   prob_mut = 0.01;
   init_bio_fract = 0.0;
   
   step_size = 100;

   effectiveXO = false;
   randSeed = 7;
 
   calculatorName = "BALANCEDACC";
   
   // establish map for parameters
   param_map["MINSIZE"] = minSizeParam;
   param_map["MAXSIZE"] = maxSizeParam;
   param_map["TAILRATIO"] = tailRatioParam;
   param_map["GROWRATE"] = growRateParam;
   param_map["MAXDEPTH"] = maxDepthParam;
   param_map["TAILSIZE"] = tailSizeParam;
   param_map["SENSIBLEINIT"] = sensibleInitParam;
   param_map["POPSIZE"] = popSizeParam;
//    param_map["NUMGEN"] = numGenParam;
   param_map["PROBCROSS"] = probCrossParam;
   param_map["PROBMUT"] = probMutParam;
   param_map["GRAMMARFILE"] = gramFileParam;
   param_map["GENSPERSTEP"] = stepSize;
   param_map["CALCTYPE"] = calcType;
   param_map["EFFECTIVEXO"] = useEffectiveXO;
   param_map["INCLUDEALLSNPS"] = useAllSnps;
   param_map["REQUIREALLVARS"] = requireAll;
   param_map["REQUIREALLONCE"] = requireAllOnce;
   param_map["BIOFILTERFRACT"] = bioInitFract;
   param_map["NUMGENSRESTRICTVARS"] = restrictVarGens;
   param_map["BIOMODELSELECTION"] = bioModelSelection;
   param_map["BLOCKCROSSGENS"] = blockCrossGens;
   param_map["RESETVARSMIGRATION"] = resetVarsAtMigration;
   param_map["BACKPROPFREQ"] = bpfreq;
   param_map["BACKPROPSTART"] = bpstart;
   param_map["GASELECTION"] = gaSelection;
#ifdef ATHENA_BLOAT_CONTROL
   param_map["DOUBTOURNF"] = doubleTournF;
   param_map["DOUBTOURND"] = doubleTournD;
   param_map["DOUBTOURNFITFIRST"] = doubleTournFitFirst;
   param_map["PRUNEPLANT"] = prunePlantFract;
#endif
   
   BioModelSelectionMap["ROULETTE"] = rouletteSelect;
   BioModelSelectionMap["ORDERED"] = orderedSelect;
   
#ifdef ATHENA_BLOAT_CONTROL
   GASelectorMap["DOUBLE"] = DoubleTournamentSelection;
#endif
   GASelectorMap["ROULETTE"] = RouletteWheelSelection;
   
   biofilter_selector_type = orderedSelect;
#ifdef ATHENA_BLOAT_CONTROL
   pruneAndPlantFract = 0.0;
#endif

   parameters_set = false;
   useAllVars = false;
   requireAllVars = false;
   requireAllVarsOnce = false;
   num_genotypes = 0;
   num_continuous = 0;
   ngens_var_restrict = 0;
   // default is to add any additional variables to migration instead of replacing
   reset_restricted_at_migration = false;
   ngens_block_cross = 0;
   
   ga = NULL;
   gelog = NULL;
   
   maxbest = true;
   
   // these control when first backpropagation occurs and how often thereafter
   bp_first_gen = -1;
   bp_freq_gen = 0;
   bp_next_opt = -1;
   
   logTypeSelected = LogNone;
   gaSelector = RouletteWheelSelection;
#ifdef ATHENA_BLOAT_CONTROL
   fitfirst = false;
   doubletourneyF = 7;
   doubletourneyD = 1.4;
#endif
   
   #ifdef PARALLEL
     genomeInfo = 112;
     genomeArray=113;
   #endif
   
}


///
/// Sets parameters for use with GAlib
/// @param alg_params AlgorithmParams
/// @throws AthenaExcept on error
///
void GENNAlg::set_ga_params(){   
    GARandomSeed(randSeed);
    srand(randSeed);   
    // first free ga memory if already run
    free_memory();
    
    //Initialize the GE mapper
	 //Set maximum number of wrapping events per mapping
	 mapper.setMaxWraps(wrapEvents);
	 
    // adjust grammar when ott dummy encoding has been used
    // or for shorthand ways of specifying variables
    adjuster.read_grammar_file(grammarFile);
    if(useAllVars){
      adjuster.include_all_vars(num_genotypes, num_continuous);
    }
    else{
      adjuster.expand_variables();
      if(dummy_encoded)
        adjuster.double_genotype_grammar();
     }
    adjuster.set_mapper(mapper);

    set_mapper_prefs(mapper);
	  
	  if(!requireAllVars && !requireAllVarsOnce)
      GEObjective::setSolutionType("NN", calculatorName);
    else if(requireAllVars){
      vector<string> var_strings = adjuster.get_variables();
      GEObjective::setSolutionType("NNALL", calculatorName, var_strings);
    }
    else if(requireAllVarsOnce){
      vector<string> var_strings = adjuster.get_variables();
      GEObjective::setSolutionType("NNONCE", calculatorName, var_strings);
    }

   if(calculatorName.compare("RSQUARED")==0){
     pop.setConvertScores(true);
   }

   // add mapper to objective function
   GEObjective::setMapper(&mapper);
   if(logTypeSelected == LogNone)
     GEObjective::addLogging(false);
   else
     GEObjective::addLogging(true);

   InitGEgenome::setMinSize(minSize);
   InitGEgenome::setMaxSize(maxSize);
   InitGEgenome::setMapper(&mapper);
}


///
/// sets mapper preferences
/// @param hemannMapper AthenaGrammarSI
///
void GENNAlg::set_mapper_prefs(AthenaGrammarSI& hemannMapper){
    // Mapper settings.
	  hemannMapper.setGenotypeMaxCodonValue(INT_MAX);
	  hemannMapper.setPopSize(pop_size);
	  hemannMapper.setGrow(growRate);
	  hemannMapper.setMaxDepth(maxDepth);
	  
	  if(tailSize){
		  hemannMapper.setTailSize(tailSize);
	  }
	  else{
		  hemannMapper.setTailRatio(tailRatio);
	  }
	  hemannMapper.setRestrictRule("<v>");
	  
	  // set up reverse grammar that can be used with resetting genome based
	  // on optimization (such as backpropagation)
    hemannMapper.constructReverseGrammar();
}


///
/// Establishes the algorithm for the run based on the parameters
/// set in this algorithm.
///
void GENNAlg::run(){
    while(!ga->done()){
        ga->step();
    }
}


///
/// Runs an indicated interval for the algorithm. This interval is set
/// in the step_size variable.  
///
void GENNAlg::step(){

GE1DArrayGenome::myrank = myRank;
    for(unsigned int i=0; i < step_size; i++){
        if(!ga->done()){
            if(logTypeSelected!=LogNone){
                gelog->add_generation();
            }
         ga->step();

         restrict_steps_done++;
          if(ngens_var_restrict && restrict_steps_done == ngens_var_restrict){
           // have to convert the networks back to original mapping
           GEObjective::setMapper(&mapper);
           GE1DArrayGenome::setMapper(&mapper);
           convertNetworks(restrictMapper, mapper);
         }
         
         // check for need to change crossovers
         if(ngens_block_cross && restrict_steps_done == ngens_block_cross){
            ResetCrossover();
         }
 
         // check to see if the back propagation should be done at this generation
         if(int(restrict_steps_done) == bp_next_opt && bp_first_gen >= 0){
          runBackPropagation();
          bp_next_opt += bp_freq_gen;
         }
         fillLog();
        }
    }

    #ifdef PARALLEL
      // if restricted variables has been used and are still in effect
      // need to convert all networks back to original grammar and exchange
      // then construct new grammar restricted to only variables in the 
      // new population and continue to process
      if(ngens_var_restrict && restrict_steps_done < ngens_var_restrict){
      // have to convert the networks back to original mapping
         GEObjective::setMapper(&mapper);
         GE1DArrayGenome::setMapper(&mapper);
         convertNetworks(restrictMapper, mapper);
      } 

      SendAndReceiveStruct(totalNodes, myRank);

      // when running in parallel transfer around populations
      if(ngens_var_restrict && restrict_steps_done < ngens_var_restrict){
        // after transfer construct new grammar and set for use
        setRestrictedGrammar(reset_restricted_at_migration);
      }
    #endif
    
    // only need to fill population at this point not at end of each generation
    fill_population();
}



///
/// Fills the algorithm log with data from the population
///
void GENNAlg::fillLog(){
    
    // don't fill log if not being kept
    if(logTypeSelected==LogNone){
      return;
    }

    unsigned int num_inds = ga->population().size();
    
    float worst_score = GEObjective::get_worst_score();
    for(unsigned int curr_ind =0; curr_ind < num_inds; curr_ind++){
  
      GE1DArrayGenome& genome = (GE1DArrayGenome&)(ga->population().individual(curr_ind));

      if(worst_score != genome.score()){
        gelog->add_network();
        gelog->add_nn_size(genome.getEffectiveSize());
        gelog->add_nn_depth(genome.getDepth());
        if(!pop.getConvertScores()){
          gelog->add_fitness(genome.score(), genome.get_genos());
        }
        else{  // add converted score to log file
          gelog->add_fitness(pop[0]->adjust_score_out(genome.score(), genome.getNumIndsEvaluated(),
            genome.getSSTotal()), genome.get_genos());
        }
        gelog->add_num_genos(genome.getNumGenes());
        gelog->add_num_covars(genome.getNumCovars());
        gelog->add_epochs(genome.getNumEpochsTrained());
        // zero out epochs
        genome.setNumEpochsTrained(0);
      } 
    }
    gelog->complete_gen();
    
    #ifdef PARALLEL
      gelog->SendReceiveLogs(totalNodes, myRank);
//     if(myRank==0)
//     
//       gelog->receiveLogs(totalNodes);
//     else
//       gelog->sendLog();
    #endif
    
    
    // output log information -- appended to existing log files
    if(myRank==0){
        writeLog();
// exit(1);
    }
}


///
/// Gets log and saves in vector
///
void GENNAlg::saveLog(){

  // need to get all slaves logs when running parallel
  #ifdef PARALLEL
    if(myRank==0)
      gelog->receiveLogs(totalNodes);
    else
      gelog->sendLog();
  #endif

  logs.push_back(gelog);
}

///
/// Starts fresh log 
///
void GENNAlg::startLog(int num_snps){
  if(gelog != NULL){
    delete gelog;
    gelog=NULL;
  }
  gelog = new NNLog(num_snps);
  gelog->set_max_best(maxbest);
}


///
/// Outputs current population of genome 
/// 
void GENNAlg::output_alg_inds(){

   unsigned int num_inds = ga->population().size();
   
    for(unsigned int curr_ind = 0; curr_ind < num_inds; curr_ind++){
      GE1DArrayGenome& genome = (GE1DArrayGenome&)(ga->population().individual(curr_ind));
      mapper.setGenotype(genome);  
      Phenotype const *phenotype=mapper.getPhenotype();
      unsigned int pheno_size=(*phenotype).size();
      vector<string> symbols(pheno_size, "");
     
      for(unsigned int i=0; i<pheno_size; ++i){
         symbols[i] = *((*phenotype)[i]);
      }
      cout << "IND in ga: " << curr_ind+1 << " score: " << genome.score() << endl;
           
    }
   
}


///
/// Use in analyzing algorithm.  Checks population
/// and reports information on solutions that contain
/// the desired variables.
///
void GENNAlg::AnalyzePopulation(){

  // place all models into population
  fill_population();
  
  std::set<string> target5;
  target5.insert("G9");
  target5.insert("G10");

  std::set<string> target10;
  target10.insert("G19");
  target10.insert("G20");
  
  int modrank=1;
  bool contains5, contains10;
  vector<string> symbols;
  // examine each model and output rank and score when matching desired variables
  for(Solution* sol=pop.best(); sol != NULL; sol=pop.GetNext()){
    symbols = sol->get_symbols();
    contains5 = contains10 = false;
    
    vector<string> vars;
    for(vector<string>::iterator iter = symbols.begin(); iter != symbols.end();
      ++iter){
      std::size_t loc = iter->find_first_of("G");
      if(loc == 0){
        vars.push_back(*iter);
      }
      if(target5.find(*iter) != target5.end()){
        contains5 = true;
      }
      else if(target10.find(*iter) != target10.end()){
        contains10 = true;
      }
    }
    
    // print out all variables in model along with score and rank in the population
    // when solution contains any of the indicated variables
    if(contains5 || contains10){
      string variables;
      for(vector<string>::iterator iter = vars.begin(); iter != vars.end();
        ++iter){
        variables += *iter + " ";
      }
    }    
    modrank++;
  }

}



///
/// Fills the population after a step to allow for exchange of models with 
/// other algorithms.
/// 
void GENNAlg::fill_population(){

    unsigned int num_inds = ga->population().size();
    
    // clear solution population
    pop.clear();
    
    for(unsigned int curr_ind = 0; curr_ind < num_inds; curr_ind++){
        Solution* sol = GEObjective::getBlankSolution();
      GE1DArrayGenome& genome = (GE1DArrayGenome&)(ga->population().individual(curr_ind));
      mapper.setGenotype(genome);
      Phenotype const *phenotype=mapper.getPhenotype();
      
      unsigned int pheno_size=(*phenotype).size();
      vector<string> symbols(pheno_size, "");
      
      for(unsigned int i=0; i<pheno_size; ++i){
         symbols[i] = *((*phenotype)[i]);
      }
      
      sol->set_symbols(symbols);
      sol->fitness(genome.score());
      sol->testval(genome.getTestValue());
      pop.insert(sol);
    }
}


///
/// Establishes test data set values for the population
/// @param test_set Dataset containing individuals for the test set
///
void GENNAlg::test_solution(Dataset* test_set, int nproc){
   
   // use first dataset
   set = test_set;
   GEObjective::setDataset(test_set);
 
  for(int i=0; i < nproc; i++){
    GE1DArrayGenome bestgenome = (GE1DArrayGenome&)ga->population().best(i);
    float test_score = GEObjective::GEObjectiveFunc(bestgenome);
    ((GE1DArrayGenome&)ga->population().best(i)).setTestValue(test_score);
    pop[i]->testval(test_score);
  }
}

///
/// Evaluates and outputs dataset to indicated output stream.
/// @param set Dataset
/// @param os ostream to write to
/// @param model index of model to output
/// 
void GENNAlg::output_ind_evals(Dataset* set, ostream& os, int model){
  GEObjective::setDataset(set);
  GE1DArrayGenome bestgenome = (GE1DArrayGenome&)ga->population().best(model);
  GEObjective::GEObjectiveFuncOut(bestgenome, os);
}

///
/// Returns the snps and covariates present in the network
///
vector<string> GENNAlg::getBestVariables(){
  vector<string> symbols;
  return symbols;
}


///
/// Sets all genomes in a population to either effective crossover or one point
///
void GENNAlg::ResetCrossover(){
    if(effectiveXO){
      ga->crossover(GE1DArrayGenome::effCrossover);
    }
    else{
      ga->crossover(GE1DArrayGenome::OnePointCrossover);
    }
}

///
/// Initializes the algorithm
///
void GENNAlg::initialize(){
    restrict_steps_done=0;
    
    // free ga if already established
    free_memory();
    
    // establish genome type and assign intialization and objective functions
    GE1DArrayGenome genome(50);

    // need to add objective function
    genome.evaluator(GEObjective::GEObjectiveFunc);
    
    if(sensibleInit){
        genome.initializer(InitGEgenome::initFuncSI);
    }
    else{
        genome.initializer(InitGEgenome::initFuncRandom);
    }
    
    // assign crossover type
    if(ngens_block_cross > 0)
      genome.crossover(GE1DArrayGenome::blockCrossover);
    else if(!effectiveXO){
        genome.crossover(GE1DArrayGenome::OnePointCrossover);
    }
    else{
        genome.crossover(GE1DArrayGenome::effCrossover);
    }
    
    // use point mutator for GEListGenome
    genome.mutator(GE1DArrayGenome::codonMutator);
    
    // controls allowable sizes
    genome.resizeBehaviour(1, maxSize);
    GEObjective::setMaxGenomeSize(maxSize);
    
    // Set up the algorithm
    ga = new GASimpleGA(genome);
    
#ifdef ATHENA_BLOAT_CONTROL
    ga->setPrunePlant(pruneAndPlantFract);
    ga->pruneplant(GE1DArrayGenome::prune_and_plant);
#endif

    // Set selector type
    GASelectionScheme* selector;
    switch(gaSelector){
      case NoMatchSelector:
      case RouletteWheelSelection:
        selector = new GARouletteWheelSelector;
        break;
#ifdef ATHENA_BLOAT_CONTROL
      case DoubleTournamentSelection:
        selector = new GADoubleTournamentSelector(doubletourneyD, doubletourneyF, fitfirst);
        break;
#endif
    };
    
    ga->selector(*selector);
    // safe to delete selector as algorithm clones it
    delete selector; 
    // scaling
    ga->scaling(GANoScaling());
    // individuals in population
    ga->populationSize(pop_size);
    ga->pMutation(prob_mut);        
    ga->pCrossover(prob_cross);
    ga->nGenerations(num_generations);

    if(GEObjective::max_best()){
      ga->maximize();
      maxbest = true;
    }
    else{
      ga->minimize();
      pop.sort_ascending();
      maxbest = false;
      gelog->set_max_best(maxbest);
    }
   
    mapper.resetGrammarModels();
    
    mapper.setVariableCodonMap();
    GE1DArrayGenome::setMapper(&mapper);
    ga->initialize();

    // when restricted variables are requested 
    // fill new mapper with grammar with variables from original mapper
    if(ngens_var_restrict > 0){
      setRestrictedGrammar(true);
    }
    
    // first optimization of backpropagation 
    bp_next_opt = bp_first_gen;

    // run optimization after initialization when indicated
    if(bp_next_opt == 0){
      runBackPropagation();
      bp_next_opt += bp_freq_gen;
    }
    fill_population();
    fillLog();
}


///
/// Runs backpropagation on current population of genomes
///
void GENNAlg::runBackPropagation(){
  
  unsigned int num_inds = ga->population().size();
  
  for(unsigned int curr_ind = 0; curr_ind < num_inds; curr_ind++){
    GEObjective::optimizeSolution(ga->population().individual(curr_ind));
  }
  ga->evaluatePop();
}


///
/// Establish restricted grammar based on current population
/// @param clearVariables when true only variables in current populations will
/// be included in the new restricted grammar.  When false it will keep existing
/// variables and add new ones from population.
///
void GENNAlg::setRestrictedGrammar(bool clearVariables){
  if(clearVariables)
    adjuster.clear_variables();
  fill_population();
  for(Solution* sol=pop.best(); sol != NULL; sol=pop.GetNext()){
    adjuster.add_variables(sol->get_symbols());
  }
  adjuster.edit_only_var_included();
  adjuster.set_mapper(restrictMapper, myRank);
  set_mapper_prefs(restrictMapper);
  restrictMapper.setVariableCodonMap();
  GEObjective::setMapper(&restrictMapper);
  GE1DArrayGenome::setMapper(&restrictMapper);
  convertNetworks(mapper, restrictMapper);
}


///
/// Converts networks from one mapping to other
/// @param currentMapper
/// @param newMapper
///
void GENNAlg::convertNetworks(AthenaGrammarSI& currentMapper, AthenaGrammarSI& newMapper){
    unsigned int num_inds = ga->population().size();
    for(unsigned int curr_ind = 0; curr_ind < num_inds; curr_ind++){
      GE1DArrayGenome& genome = (GE1DArrayGenome&)(ga->population().individual(curr_ind));            
      currentMapper.convertGenomeVariables(newMapper, genome);
      // reset the genome codons to the new values
      unsigned int genome_size = currentMapper.getGenotype()->size();
      genome.resize(genome_size);
      int i=0;

 	    // Now copy genotype onto genome
    	Genotype::const_iterator genIt=(currentMapper.getGenotype())->begin();
    	while(genIt!=(currentMapper.getGenotype())->end()){
	      genome.gene(i, *genIt);
    	  genIt++;
	      i++;
    	}  	
    }
    ga->statistics().setBestIndividual(ga->population());
}


void GENNAlg::outputGenome(GAGenome& g){
  GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);
  cout << "Score: " << genome.score() << " ";
  
  int length = genome.length();
  for(int i=0; i<length; i++){
    cout << (unsigned int)(genome.gene(i)) << " ";
  }
  cout << endl;
}

///
/// Prepares for logging
/// @param outname
/// @param cv
///
void GENNAlg::prepareLog(string basename, int cv){
    int totalPopSize = pop_size;
    
    main_log_filename = basename + ".cv." + Stringmanip::itos(cv) + ".log";
    fitness_log_filename =  basename + ".cv." + Stringmanip::itos(cv) +
        ".fitness.log";    
    snpname_log_filename = basename + ".cv." + Stringmanip::itos(cv) + 
        ".snpsize.log";    
  #ifdef PARALLEL
    if(myRank==0){
  #endif
  
  ofstream outfile;
  switch(logTypeSelected){
    case LogDetailed:
        outfile.open(fitness_log_filename.c_str(), ios::out);
        if(!outfile.is_open()){
            throw AthenaExcept(fitness_log_filename + " unable to open for writing log file");
          }   
      gelog->output_fitness_headers(outfile);
      outfile.close();  
      outfile.open(snpname_log_filename.c_str(), ios::out);
      if(!outfile.is_open()){
            throw AthenaExcept(snpname_log_filename + " unable to open for writing log file");
       }   
       gelog->output_snp_headers(outfile);
      outfile.close();      
    case LogSummary:
    outfile.open(main_log_filename.c_str(), ios::out);
        if(!outfile.is_open()){
            throw AthenaExcept(main_log_filename + " unable to open for writing log file");
        }  
      gelog->output_main_headers(outfile); 
      outfile.close();        
    case LogNone:
    break;
  }
  
  
  #ifdef PARALLEL
    }
  #endif
}


///
/// Writes current log to output
/// @param outname
///
void GENNAlg::writeLog(){

    ofstream outfile;
    outfile.open(main_log_filename.c_str(), ios::app);
    if(!outfile.is_open()){
        throw AthenaExcept(main_log_filename + " unable to open for writing log file");
     }
     gelog->output_log(outfile);
     outfile.close();

//  int totalPopSize = pop_size;
//   #ifdef PARALLEL
//     if(logTypeSelected != LogNone){
//       if(myRank==0)
//         gelog->receiveLogs(totalNodes);
//       else
//         gelog->sendLog();      
//       if(logTypeSelected == LogDetailed){
//         if(myRank==0)
//           gelog->receiveDetailedLogs(totalNodes);
//         else
//           gelog->sendDetailedLog();
//       }
//     }
//     totalPopSize *= totalNodes;
//     
//     if(myRank==0){
//   #endif
// 
//   if(logTypeSelected != LogNone){
//     string filename = basename + ".cv." + Stringmanip::itos(cv) + ".log";
//     ofstream outfile;
//     outfile.open(filename.c_str(), ios::out);
//       if(!outfile.is_open()){
//         throw AthenaExcept(filename + " unable to open for writing log file");
//       }   
//       gelog->output_log(outfile);
//       outfile.close();
//     if(logTypeSelected == LogDetailed){
//       string fitnessname = basename + ".cv." + Stringmanip::itos(cv) +
//         ".fitness.log";
//       ofstream fitout;
//       fitout.open(fitnessname.c_str(), ios::out);
//       if(!fitout.is_open()){
//         throw AthenaExcept(fitnessname + " unable to open for writing log file");
//       }
//       gelog->output_fitness(fitout, totalPopSize);
//       fitout.close();
//       
//       string snpname = basename + ".cv." + Stringmanip::itos(cv) + 
//         ".snpsize.log";
//       ofstream sizeout;
//       sizeout.open(snpname.c_str(), ios::out);
//       if(!sizeout.is_open()){
//         throw AthenaExcept(snpname + " unable to open for writing log file");
//       }
//       gelog->output_snp_sizes(sizeout, totalPopSize);
//       sizeout.close();
//     }
//   }
//   
//   #ifdef PARALLEL
//   } // end write for master
//   #endif
  
//   delete gelog;
//   gelog=NULL;
}


///
/// Clears logs from the algorithm
///
void GENNAlg::clearLogs(){
  for(unsigned int i=0; i<logs.size(); i++){
    delete logs[i];
  }
  logs.clear();
}


///
/// Writes graphical representation to stream
/// @param os
/// @param sol Solution to output to stream
///
void GENNAlg::writeGraphical(ostream& os, Solution* sol, data_manage::Dataholder* holder,
  bool map_used, bool ott_dummy){
  GEObjective::outputModel(os, sol, holder, map_used, ott_dummy);
}

///
/// Returns graphical file extension appropriate to solution
/// @return extension
///
std::string GENNAlg::getGraphicalFileExt(){
  return GEObjective::getGraphicalExt();
}



///
/// Retrieves the models from BioFilter and stores the information in the algorithm
/// @param filename File with biofilter models
/// @param biofiletype Type of biological filter file (BINARY or TEXT)
/// @param holder Dataholder that contains map file info to convert model names
/// into indexes in the 
///
void GENNAlg::getBioModels(std::string filename, std::string biofiletype, data_manage::Dataholder* holder){
  
  BioFilterModelCollection collection(filename, 100000, biofiletype);
  setBioModels(collection, holder);
}



///
/// Fills biomodel collection from archive files 
/// @param genegeneFile genegene filename
/// @param archiveFile arcchive filename
/// @param holder Dataholder that contains map file info to convert model names
/// into indexes
///
void GENNAlg::getBioModelsArchive(string genegeneFile, string archiveFile, data_manage::Dataholder* holder){
  BioFilterModelCollection collection(genegeneFile, archiveFile, 10000);
  setBioModels(collection, holder);
}

///
/// Sets the bio models based on collection passed
///
void GENNAlg::setBioModels(BioFilterModelCollection& collection, data_manage::Dataholder* holder){

  vector<BioModel> biomodels;

  int num_models_needed = (unsigned int)(init_bio_fract * pop_size);

  if(biofilter_selector_type == rouletteSelect){  
    for(int i=0; i<num_models_needed; i++){
      biomodels.push_back(collection.GetRandomModel());
    }
  }
  else{
    // set starting point for this algorithm -- in parallel mode this will give each 
    // node a different set
    collection.SetStartModel(myRank * num_models_needed);
    for(int i=0; i<num_models_needed; i++){
      biomodels.push_back(collection.GetNextModel());
    }
  }

  mapper.clearModels();

  // from models find the variables to use and store in mapper for use in initializing dataset
  for(vector<BioModel>::iterator iter=biomodels.begin(); iter != biomodels.end(); iter++){
    vector<int> indexes;
    int index;
    for(unsigned int i=0; i<iter->idstring.size(); i++){
      index = holder->get_geno_index(iter->idstring[i]);
      indexes.push_back(index);
    }
    mapper.addModel(indexes);
  }
	  mapper.setModelCodons(dummy_encoded);
}



#ifdef PARALLEL

void GENNAlg::setRank(int rank){
        Algorithm::setRank(rank); 
        InitGEgenome::setrank(rank);
        GEObjective::setrank(rank);
}


///
/// Alternate MPI messaging taking advantage of MPI_Allgather functionality using struct
/// @param totalNodes
/// @param myRank
///
void GENNAlg::SendAndReceiveStruct(int totalNodes, int myRank){

  struct_mpi * send = new struct_mpi;
  GE1DArrayGenome& genome = (GE1DArrayGenome&)ga->statistics().bestIndividual();
  // package genome Info
  send->genomeParams[0] = genome.length();
  send->genomeParams[1] = genome.score();;
  send->genomeParams[2] = genome.getEffectiveSize();
  send->genomeParams[3] = genome.getNumGenes();
  send->genomeParams[4] = genome.getNumCovars();
  send->genomeParams[5] = genome.getNumIndsEvaluated();
  send->genomeParams[6] = genome.getSSTotal();
  
  // package codons
  for(int i=0; i<genome.length(); i++){
    send->codons[i] = genome.gene(i);
  }
  
  // prepare receiving array
  struct_mpi * recv = new struct_mpi[totalNodes];
  
  MPI_Allgather (send, sizeof(*send), MPI_BYTE, recv, sizeof(*send), MPI_BYTE, MPI_COMM_WORLD);
  
  updateWithMigration(recv, totalNodes, myRank);
  
  delete send;
  delete [] recv;
  
}


///
/// Alternate MPI messaging taking advantage of MPI_Allgather functionality
/// @param totalNodes
/// @param myRank
///
void GENNAlg::SendAndReceive(int totalNodes, int myRank){

  // each sends the basic info for best model
  int nParams = 7;
  float * genomeParams = new float[nParams];
  int totalParams = nParams * totalNodes;
  float * allGenomeParams = new float[totalParams];

  GE1DArrayGenome& genome = (GE1DArrayGenome&)ga->statistics().bestIndividual();
  // package genome Info
  genomeParams[0] = genome.length();
  genomeParams[1] = genome.score();;
  genomeParams[2] = genome.getEffectiveSize();
  genomeParams[3] = genome.getNumGenes();
  genomeParams[4] = genome.getNumCovars();
  genomeParams[5] = genome.getNumIndsEvaluated();
  genomeParams[6] = genome.getSSTotal();

  MPI_Allgather(genomeParams, nParams, MPI_FLOAT, allGenomeParams, nParams, MPI_FLOAT, MPI_COMM_WORLD);
  // determine largest genome
  int max_genome=0;
  for(int i=0; i<totalNodes; i++){
    if(allGenomeParams[i*nParams] > max_genome)
      max_genome = allGenomeParams[i*nParams];
  }
  
  int * sendGenome = new int[max_genome];
  int total_genome_size = max_genome * totalNodes;
  int * totalGenomes = new int[total_genome_size];

  for(int i=0; i<genome.length(); i++){
    sendGenome[i]=genome.gene(i);
  }

  // then each sends and receives the actual codons and
  // incorporates those into its own population
  // needs to be large enough to hold the largest genome
  MPI_Allgather(sendGenome, max_genome, MPI_INT, totalGenomes, max_genome, MPI_INT, MPI_COMM_WORLD);
  updateWithMigration(allGenomeParams, totalGenomes, totalNodes, myRank, max_genome);
  
  delete [] genomeParams;
  delete [] allGenomeParams;
  delete [] sendGenome;
  delete [] totalGenomes;
}


///
/// Incorporates migration into population
/// @param genomes 
/// @param totalNodes
/// @param myRank
///
void GENNAlg::updateWithMigration(struct_mpi* mpi_genomes, int totalNodes, int myRank){
  GAPopulation pop(ga->population());
  
  for(int node=0; node < totalNodes; node++){
    if(myRank==node){
      continue;
    }
    
    GAGenome *tmpind = ga->population().individual(0).clone();
    GE1DArrayGenome& genome = (GE1DArrayGenome&)*tmpind;
    int len = mpi_genomes[node].genomeParams[0];
    genome.length(len);
    genome.setEffectiveSize(mpi_genomes[node].genomeParams[2]);
    genome.setNumGenes(mpi_genomes[node].genomeParams[3]);
    genome.setNumCovars(mpi_genomes[node].genomeParams[4]);
    genome.setNumIndsEvaluated(mpi_genomes[node].genomeParams[5]);
    genome.setSSTotal(mpi_genomes[node].genomeParams[6]);
    for(int i=0; i<len; i++){
      genome.gene(i, mpi_genomes[node].codons[i]);
    }
    genome.score(mpi_genomes[node].genomeParams[1]);

    pop.add(genome);
    
    delete tmpind;
  }

  // remove worst individuals from population
  for(int i=0; i < totalNodes-1; i++)
    pop.destroy();
    
  ga->population(pop);

}


///
/// Incorporates migration into population
/// @param genomes 
/// @param totalNodes
/// @param myRank
///
// void GENNAlg::updateWithMigration(struct_mpi* mpi_genomes, int totalNodes, int myRank){
//   GAPopulation pop(ga->population());
//   
//   for(int node=0; node < totalNodes; node++){
//     if(myRank==node){
//       continue;
//     }
//     
//     GAGenome *tmpind = ga->population().individual(0).clone();
//     GE1DArrayGenome& genome = (GE1DArrayGenome&)*tmpind;
//     int len = mpi_genomes[node].genomeParams[0];
//     genome.length(len);
//     genome.setEffectiveSize(mpi_genomes[node].genomeParams[2]);
//     genome.setNumGenes(mpi_genomes[node].genomeParams[3]);
//     genome.setNumCovars(mpi_genomes[node].genomeParams[4]);
//     genome.setNumIndsEvaluated(mpi_genomes[node].genomeParams[5]);
//     genome.setSSTotal(mpi_genomes[node].genomeParams[6]);
// if(myRank==1){
//   cout << "in update " << myRank << " len=" << len << " effsize=" << mpi_genomes[node].genomeParams[2] << " numgenes=" << mpi_genomes[node].genomeParams[3] << endl;
// }
//     for(int i=0; i<len; i++){
//  if(myRank==1){
//    cout << "in update " << myRank << " adding for i=" << i << " codon " << mpi_genomes[node].codons[i] << endl;
//  }
//       genome.gene(i, mpi_genomes[node].codons[i]);
//     }
//     genome.score(mpi_genomes[node].genomeParams[1]);
// 
//     pop.add(genome);
//     
//     delete tmpind;
//   }
// 
//   // remove worst individuals from population
//   for(int i=0; i < totalNodes-1; i++)
//     pop.destroy();
//     
//   ga->population(pop);
// 
// }

///
/// Incorporates migration into population
/// @param stats
/// @param codons
/// @param totalNodes
/// @param myRank
/// @param max_length Maximum length of a genome
///
void GENNAlg::updateWithMigration(float* stats, int* codons, int totalNodes, int myRank, int max_length){
  int currCodon=0, len;

  GAPopulation pop(ga->population());
//if(myRank==1){
//  cout << "myRank=" << myRank << " totalNodes=" << totalNodes << endl;
//  cout << "myRank=" << myRank << " max_length=" << max_length << endl;
//}
  
  for(int node=0; node < totalNodes; node++){
    if(myRank == node){ // skip when same node 
      if(!max_length){
        currCodon += (unsigned int)(stats[node*7]);
      }
      else{
        currCodon += max_length;
      }
      continue;
    }
        
    GAGenome *tmpind = ga->population().individual(0).clone();
   
    GE1DArrayGenome& genome = (GE1DArrayGenome&)*tmpind;

    int number_params = 7;
    len = stats[node*number_params];
    genome.length(len);
    genome.setEffectiveSize(stats[node*number_params+2]);
    genome.setNumGenes(stats[node*number_params+3]);
    genome.setNumCovars(stats[node*number_params+4]);
    genome.setNumIndsEvaluated(stats[node*number_params+5]);
    genome.setSSTotal(stats[node*number_params+6]);
    
    int final_codon = max_length + currCodon;
    
    for(int i=0; i<len; i++){
//if(myRank==1){
//  cout << "in update " << myRank << " adding for i=" << i << " codon " << codons[currCodon] << endl;
//}
      genome.gene(i, codons[currCodon++]);
    }

    if(max_length && currCodon < final_codon){
      currCodon = final_codon;
    }

    genome.score(stats[node*number_params+1]);

    pop.add(genome);
    
    delete tmpind;
  }

  // remove worst individuals from population
  for(int i=0; i < totalNodes-1; i++)
    pop.destroy();
    
  ga->population(pop);
}

#endif
