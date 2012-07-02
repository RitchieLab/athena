#include "GEObjective.h"
#include <math.h>
#include "Structs.h"

AthenaGrammarSI* GEObjective::mapper = NULL;
data_manage::Dataset* GEObjective::set = NULL;
SolutionCreator* GEObjective::sol_creator = NULL;
unsigned int GEObjective::maxGenSize = 250;
bool GEObjective::additional_logging = false;

int GEObjective::rank=0;

///
/// Objective function
/// @param g GAGenome which is being evaluated
/// @return fitness
///
float GEObjective::GEObjectiveFunc(GAGenome& g){

   GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);

	//Assign genotype to mapper
	mapper->setGenotype(genome); 
   Phenotype const *phenotype=mapper->getPhenotype();
   
   float fitness;
   genome.setValid(phenotype->getValid());

   if(phenotype->getValid()){
       unsigned int pheno_size=(*phenotype).size();
       vector<string> symbols(pheno_size, "");     
       
      for(unsigned int i=0; i<pheno_size; ++i){
          symbols[i] = *((*phenotype)[i]);
      }

      sol_creator->establish_solution(symbols, set);

      fitness = sol_creator->evaluate(set);
     
      if(additional_logging){
        sol_creator->detailed_logging();
        genome.setDepth(sol_creator->get_detailed_log());
        genome.setGramDepth(mapper->buildDerivationTree());
//cout << "added depth=" << genome.getDepth() << endl;
      }
      
      sol_creator->free_solution();

      genome.setEffectiveSize(mapper->getGenotype()->getEffectiveSize());   
      genome.setNumCovars(sol_creator->get_num_covars());
      genome.setNumGenes(sol_creator->get_num_genes());
      genome.add_genos(sol_creator->getGeneIndexes());
      genome.add_covars(sol_creator->getCovarIndexes());
      genome.setNumIndsEvaluated(sol_creator->getNumIndsEvaluated());
      
      // when set 
      if(genome.getNumIndsEvaluated() != int(set->num_inds())){
        genome.setSSTotal(sol_creator->get_calculator_constant());
      }
      else{
        genome.setSSTotal(set->get_sstotal());
      }
   }
   else{
        // set fitness to worst score initially
     fitness = sol_creator->get_worst();
   }
  return fitness;
}

///
/// Outputs symbols from genome to output
///
void GEObjective::OutputSymbols(GAGenome& g, ostream& os){
  GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);

	//Assign genotype to mapper
	mapper->setGenotype(genome);
   Phenotype const *phenotype=mapper->getPhenotype();  
       unsigned int pheno_size=(*phenotype).size();
       vector<string> symbols(pheno_size, "");     

      for(unsigned int i=0; i<pheno_size; ++i){
          symbols[i] = *((*phenotype)[i]);        
            os << symbols[i] << " ";
      }   
    os << endl;
}


///
/// Objective function which pass output stream to evaluator so that individual evaluations
/// can be displayed.
/// @param g GAGenome which is being evaluated
/// @return fitness
///
float GEObjective::GEObjectiveFuncOut(GAGenome& g, ostream& os){

   GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);
   
	//Assign genotype to mapper
	mapper->setGenotype(genome);

   Phenotype const *phenotype=mapper->getPhenotype();
   
   float fitness;

   if(phenotype->getValid()){
 
       unsigned int pheno_size=(*phenotype).size();
       vector<string> symbols(pheno_size, "");     

      for(unsigned int i=0; i<pheno_size; ++i){
          symbols[i] = *((*phenotype)[i]);
      }

      sol_creator->establish_solution(symbols, set);
      fitness = sol_creator->evaluate_with_output(set, os);
      sol_creator->free_solution();
      
      genome.setEffectiveSize(mapper->getGenotype()->getEffectiveSize());   

   }
   else{
        // set fitness to worst score initially
     fitness = sol_creator->get_worst();
   }
   return fitness;
}


///
/// Replaces blocks
/// @param g GAGenome
/// @param blocks codonBlocks defining the blocks that will be replaced
///
void GEObjective::insertBlocks(GE1DArrayGenome& genome, vector<AthenaGrammarSI::codonBlocks>& blocks){

    // go through original list and copy to new genome 
    // use vector to create new genome list
    vector<int> new_codons;
    unsigned int curr_block = 0;
 
    vector<int> tempcodons;
    tempcodons.push_back(0);
    tempcodons.push_back(1);
    tempcodons.push_back(0);
    tempcodons.push_back(0);
    tempcodons.push_back(1);
    tempcodons.push_back(0);
    tempcodons.push_back(1);
    tempcodons.push_back(4);
    tempcodons.push_back(9);
    tempcodons.push_back(9);
    tempcodons.push_back(9);
    tempcodons.push_back(9);
    tempcodons.push_back(2);
    
    tempcodons.push_back(1);
    tempcodons.push_back(3);
    tempcodons.push_back(1);
    tempcodons.push_back(4);
    tempcodons.push_back(7);
    
    // replace existing operator chunks with new values and set that to be the genome
    // copy entire genome including possibly unused portion at end after effective size
    for(int i=0; i<genome.size();){
      if(curr_block < blocks.size() && i==blocks[curr_block].start){ //  && curr_block < new_weights.size()){
        // set i to be equal to end of original block so that copying will continue
        // from the correct position
        i=blocks[curr_block].end;
        // in mapper have function that returns codon list for a specified value
        new_codons.insert(new_codons.end(), tempcodons.begin(), tempcodons.end());
        curr_block++;
      }
      else{
        new_codons.push_back(genome.gene(i));
        i++;
      }
    }
    // resize and copy new codons into the genome
    genome.resize(new_codons.size());
    int i=0;
    for(vector<int>::iterator iter=new_codons.begin(); iter != new_codons.end(); ++iter){
      genome.gene(i++, *iter);
    }
}


///
/// Optimizes current model using process provided by SolutionCreator
/// @param GAGenome to optimize
///
void GEObjective::optimizeSolution(GAGenome& g){

  float old_score, opt_score;
  
  GE1DArrayGenome& genome = static_cast<GE1DArrayGenome&>(g);
	//Assign genotype to mapper
	vector<AthenaGrammarSI::codonBlocks> blocks = mapper->setGenotypeOpt(genome, 
	  sol_creator->getOptIncluded(), sol_creator->getStartOptSymbol());

  Phenotype const *phenotype=mapper->getPhenotype();
  
  // run optimization on the model if it is a valid solution
  if(phenotype->getValid()){
    unsigned int pheno_size=(*phenotype).size();
    vector<string> symbols(pheno_size, "");     

    old_score = genome.score();

    for(unsigned int i=0; i<pheno_size; ++i){
     symbols[i] = *((*phenotype)[i]);
    }

    int numEpochsTrained = sol_creator->optimizeSolution(symbols, set);
   
    // after optimization, get the new constant list
    vector<symbVector> new_weights = sol_creator->getOptimizedSymbols();
    opt_score = sol_creator->getOptimizedScore();

    // skip the optimization if it isn't improving 
    if(opt_score < genome.score()){

    // go through original list and copy to new genome 
    // use vector to create new genome list
      vector<int> new_codons;
      unsigned int curr_block = 0;    
    // replace existing operator chunks with new values and set that to be the genome
    // copy entire genome including possibly unused portion at end after effective size
      for(int i=0; i<genome.size();){
        if(curr_block < blocks.size() && i==blocks[curr_block].start){
        // set i to be equal to end of original block so that copying will continue
        // from the correct position
          i=blocks[curr_block].end;
          vector<int> tempcodons = mapper->translateOptValue(new_weights[curr_block]);
        // in mapper have function that returns codon list for a specified value
          new_codons.insert(new_codons.end(), tempcodons.begin(), tempcodons.end());
          curr_block++;
        }
        else{
          new_codons.push_back(genome.gene(i));
          i++;
        }
      }

    // resize if new codon size is able to fit within the max genome
      if(new_codons.size() <= maxGenSize){
        genome.resize(new_codons.size());
        int i=0;
        for(vector<int>::iterator iter=new_codons.begin(); iter != new_codons.end(); ++iter){
          genome.gene(i++, *iter);
        }
      }
    
      genome.setNumEpochsTrained(numEpochsTrained);

    // as last step evaluate new genome and set score in it to be new value
      genome.score(GEObjective::GEObjectiveFunc(genome));
    }
  }  	

}



