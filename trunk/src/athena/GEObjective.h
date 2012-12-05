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
#include "AthenaGrammarSI.h"

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
    static void setMapper(AthenaGrammarSI* m){
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
    
    static std::string calculator_name(){return sol_creator->calculator_name();}
    
    static Solution* getBlankSolution(){return sol_creator->create_new_solution();}
    
    static bool max_best(){return sol_creator->max_best();}
    
    static std::string getGraphicalExt(){return sol_creator->graphicExt();}
    
    static void outputModel(ostream& os, Solution* sol, data_manage::Dataholder* data,
      bool map_used, bool ott_dummy, bool continmap_used){
      sol_creator->establish_solution(sol->get_symbols());
      sol_creator->graphical_output(os, data, map_used, ott_dummy, continmap_used);
    }

    /// Optimizes current model using process provided by SolutionCreator
    static void optimizeSolution(GAGenome& g);
    
    /// Sets maximum genome size
    static void setMaxGenomeSize(unsigned int maxSize){maxGenSize = maxSize;}
   
    // returns worst score
    static float get_worst_score(){return sol_creator->get_worst();}
    
    static void addLogging(bool val){additional_logging=val;}

    static void setrank(int r){rank=r;}
    static int rank;
private:

    static void insertBlocks(GE1DArrayGenome& g, vector<AthenaGrammarSI::codonBlocks>& blocks);
    static AthenaGrammarSI* mapper;
    static data_manage::Dataset* set;
    static SolutionCreator* sol_creator;
    static unsigned int maxGenSize;
    static bool additional_logging;  
};

#endif	/* _GEOBJECTIVE_H */

