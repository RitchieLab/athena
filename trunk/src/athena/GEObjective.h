/* 
 * File:   GEObjective.h
 * Author: dudeksm
 *
 * Created on November 13, 2008, 3:50 PM
 */

#ifndef _GEOBJECTIVE_H
#define	_GEOBJECTIVE_H

#include <GE/ge.h>
#include <ga/ga.h>
#include <Dataset.h>
#include "SolutionFactory.h"
#include "GE1DArrayGenome.h"
#include "Population.h"
#include "HemannGrammarSI.h"

///
/// Provides objective function for use with GALIb and GE library
///
class GEObjective{
    
public:
    
    /// Used to assign * fitness values in GE
    static float GEObjectiveFunc(GAGenome& g);
    
    /// Used to output individual evaluations
    static float GEObjectiveFuncOut(GAGenome& g, ostream& os);
    
    /// Outputs symbols from genome to output
    static void OutputSymbols(GAGenome& g, ostream& os);
    
    /// sets the mapper to use
    static void setMapper(HemannGrammarSI* m){
      mapper = m;
      mapper->setLeftOptBound(sol_creator->getLeftOptBound());
      mapper->setRightOptBound(sol_creator->getRightOptBound());
      mapper->setArgSymbols(sol_creator->getOptArgSymbols());
    }
    
    /// sets the Dataset for objective function to work with
    static void setDataset(data_manage::Dataset* ds){
      set = ds; 
      sol_creator->set_calculator_constant(ds->get_sstotal());
    }
    
    /// sets Solution type for objective function
    static void setSolutionType(std::string solution_name, std::string calculator_name){
        sol_creator = SolutionFactory::create_solution(solution_name);
        sol_creator->set_calculator(calculator_name);
    }
    
    /// Alternative method for creating Solution
    static void setSolutionType(std::string solution_name, std::string calculator_name,
      vector<string>& vars){
      sol_creator = SolutionFactory::create_solution(solution_name, vars);
      sol_creator->set_calculator(calculator_name);
    }
    
    static Solution* getBlankSolution(){return sol_creator->create_new_solution();}
    
    static bool max_best(){return sol_creator->max_best();}
    
    static std::string getGraphicalExt(){return sol_creator->graphicExt();}
    
    static void outputModel(ostream& os, Solution* sol, data_manage::Dataholder* data,
      bool map_used, bool ott_dummy){
      sol_creator->establish_solution(sol->get_symbols());
      sol_creator->graphical_output(os, data, map_used, ott_dummy);
    }

    /// Optimizes current model using process provided by SolutionCreator
    static void optimizeSolution(GAGenome& g);
    
    /// Sets maximum genome size
    static void setMaxGenomeSize(unsigned int maxSize){maxGenSize = maxSize;}
   
    // returns worst score
    static float get_worst_score(){return sol_creator->get_worst();}

    static void setrank(int r){rank=r;}
    static int rank;
private:

    static void insertBlocks(GE1DArrayGenome& g, vector<HemannGrammarSI::codonBlocks>& blocks);
    static HemannGrammarSI* mapper;
    static data_manage::Dataset* set;
    static SolutionCreator* sol_creator;
    static unsigned int maxGenSize;  
};

#endif	/* _GEOBJECTIVE_H */

