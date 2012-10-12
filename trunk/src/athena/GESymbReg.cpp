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
 * File:   GESymbReg.cpp
 * Author: dudeksm
 *
 * Created on August 10, 2010, 2:03 PM
 */
#include "GESymbReg.h"
#include "GEObjective.h"
#include "Terminals.h"
#include "GENNGrammarAdjuster.h"
#include <ga/ga.h>
#include <ctime>
 
 ///
/// Sets parameters within the algorithm for the 
/// analysis.
/// @param alg_param AlgorithmParams
/// @param numExchanges total number of times the best model will be passed to other algorithms
/// when more than one algorithm
/// @param numGenos number of Genotypes in set
/// @param numContin number of Continuous variables in set
/// 
void GESymbReg::set_params(AlgorithmParams& alg_param, int numExchanges, int numGenos, int numContin){
    
// cout << "set params called in GESymbReg" << endl;   
    map<string, string>::iterator mapIter;
    
    for(mapIter = alg_param.params.begin(); mapIter != alg_param.params.end(); 
      mapIter++){
//cout << "param=" << mapIter->first << " value=" << mapIter->second << endl;        
        switch(param_map[mapIter->first]){
            case noMatchParam:
                throw AthenaExcept("No match for parameter " + mapIter->first +
                        "in Algorithm GE Symbolic Regression");
                break;
            case minSizeParam:
                minSize = Stringmanip::stouint(mapIter->second);
                break;
            case maxSizeParam:
                maxSize = Stringmanip::stouint(mapIter->second);
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
            case doubleTournF:
                doubletourneyF = Stringmanip::stouint(mapIter->second);
                break;
            case doubleTournD:
                doubletourneyD = Stringmanip::stodouble(mapIter->second);
                break;
            case doubleTournFitFirst:
                fitfirst = Stringmanip::check_true_false(mapIter->second);
                break;
            case blockCrossGens:
                ngens_block_cross = Stringmanip::stouint(mapIter->second);
                break;
            case resetVarsAtMigration:
                reset_restricted_at_migration = Stringmanip::check_true_false(mapIter->second);
                break;
#ifdef ATHENA_BLOAT_CONTROL
            case prunePlantFract:
                pruneAndPlantFract = Stringmanip::stodouble(mapIter->second);
                break;
#endif
            case bpstart:
                throw AthenaExcept("No match for parameter " + mapIter->first +
                        " in Algorithm GE Symbolic Regression");
                break;
            case bpfreq:    
                throw AthenaExcept("No match for parameter " + mapIter->first +
                        " in Algorithm GE Symbolic Regression");
                break;
            default:
                throw AthenaExcept("No match for parameter " + mapIter->first +
                        " in Algorithm GE Symbolic Regression");               
        }
    }
    
    if(step_size > num_generations){
        step_size = num_generations;
    }
    
    num_generations = step_size * numExchanges;
    num_genotypes = numGenos;
    num_continuous = numContin;

    set_ga_params();
// cout << "done with param setting" << endl;
// exit(1);    
}


///
/// Sets parameters for use with GAlib
/// @param alg_params AlgorithmParams
/// @throws AthenaExcept on error
///
void GESymbReg::set_ga_params(){   

// cout << randSeed << " Randseed " << endl;


    GARandomSeed(randSeed);
    srand(randSeed);
 
//    GAParameterList params;
//    GASimpleGA
//cout << "in set_ga_params" << endl;    
    // first free ga memory if already run
    free_memory();
    
    //Initialize the GE mapper
	 //Set maximum number of wrapping events per mapping
	 mapper.setMaxWraps(wrapEvents);
	 
	 //Load grammar
    
    // cout << "read file " << grammarFile << " total # of rules = " << mapper.size() << endl;
    
    // adjust grammar when ott dummy encoding has been used
    // or for shorthand ways of specifying variables
//     GENNGrammarAdjuster adjuster;
    adjuster.read_grammar_file(grammarFile);
    if(useAllVars){
// cout << "num_genotypes=" << num_genotypes << " num_continuous=" << num_continuous << endl;
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
      GEObjective::setSolutionType("SYMBREG", calculatorName);
    else if(requireAllVars){
      vector<string> var_strings = adjuster.get_variables();
      GEObjective::setSolutionType("SYMBREGALL", calculatorName, var_strings);
    }
    else if(requireAllVarsOnce){
      vector<string> var_strings = adjuster.get_variables();
      GEObjective::setSolutionType("SYMBREGONCE", calculatorName, var_strings);
    }

   if(calculatorName.compare("RSQUARED")==0){
     pop.setConvertScores(true);
   }

   // add mapper to objective function
   GEObjective::setMapper(&mapper);
   InitGEgenome::setMinSize(minSize);
   InitGEgenome::setMaxSize(maxSize);
   InitGEgenome::setMapper(&mapper);
//    ga->initialize(7);
}


///
/// Runs an indicated interval for the algorithm. This interval is set
/// in the step_size variable.  
///
void GESymbReg::step(){

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

      if(ngens_var_restrict && restrict_steps_done < ngens_var_restrict){
        // after transfer construct new grammar and set for use
        setRestrictedGrammar(reset_restricted_at_migration);
      }

    #endif
    
    // only need to fill population at this point not at end of each generation
    fill_population();
}



