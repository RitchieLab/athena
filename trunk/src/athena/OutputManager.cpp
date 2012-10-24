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
#include "OutputManager.h"
#include <Stringmanip.h>
#include <iomanip>
#include <sstream>

using namespace std;


///
/// Sets up new files with appropriate headers for output
///
void OutputManager::setFiles(bool mapfile_used, string fitness_name){

	string summaryName = basename + ".athena.sum";
	int width = 30;
    if(!mapfile_used){
      width = 20;
    }
	
	ofstream outfile;
	outfile.open(summaryName.c_str(), ios::out);
    outfile << setw(5) << left << "CV" << setw(width) << "Variables" << " " << setw(20) << fitness_name + " Training"
            << " " << setw(10) << "Testing" << endl;	
	outfile.close();
	
}


///
/// outputs summary of best models
/// @param pop Population containing best model
/// @param data Dataholder that can translate the snps back to original IDs
/// @param mapfile_used true when an actual map file was used and there are original
/// names to output
/// @param dummy_encoded when true genotypes need to be adjusted back to reflect
/// original genotype positions
///
void OutputManager::outputSummary(Population& pop, int currPop,
  data_manage::Dataholder& data,  bool mapfile_used, bool dummy_encoded,
  bool continmap_used, std::string fitness_name){

    string summaryName = basename + ".athena.sum";
    
    ofstream outfile;
    outfile.open(summaryName.c_str(), ios::app);
   
    Solution* bestSolution;
    
    string prefix, continprefix;
    int width = 30;
    if(!mapfile_used){
      prefix = "G";
      width = 20;
    }
    if(!continmap_used){
      continprefix = "C";
    }
    
//     outfile << setw(5) << left << "CV" << setw(width) << "Variables" << " " << setw(20) << fitness_name + " Training"
//             << " " << setw(10) << "Testing" << endl;    
     
//     for(unsigned int currPop=0; currPop < pops.size(); currPop++){
//         bestSolution = pops[currPop].best();
        bestSolution = pop.best();
        vector<int> genos = bestSolution->get_genotypes(dummy_encoded);
        vector<int> covars = bestSolution->get_covariates();
        
        outfile << setw(5) << left << currPop+1;
        stringstream ss;
        for(unsigned int g=0; g < genos.size(); g++){
            ss << prefix << data.get_geno_name(genos[g]-1) << " ";
        }
        stringstream cs;
        for(unsigned int c=0; c < covars.size(); c++){
            cs << continprefix << data.get_covar_name(covars[c]-1) << " ";
        }
        
        outfile << setw(width) << ss.str() + cs.str() << " ";
        
        outfile << setw(20) << bestSolution->fitness() << " ";
        outfile << setw(10) << bestSolution->testval();
        outfile << endl;
//     }
    
    outfile.close();
    
}



///
/// outputs a file for each best model
/// @param pops Population vector
/// @param nmodels Number of models to output
/// @param currPop population number matching the cross-validation 
/// @param scaleInfo Informtion on the scaling done in model
/// @param data
/// @param map_used
///
void OutputManager::outputBestModels(Population& pop, int nmodels, int currPop,
  string scaleInfo, data_manage::Dataholder& data, bool map_used, bool ott_dummy, bool continmap_used){
  
    Solution* bestSolution;   
//     for(unsigned int currPop=0; currPop < pops.size(); currPop++){
      for(int mod=0; mod < nmodels; mod++){
        ofstream outfile;
        
        string currFileName = basename + ".cv" + Stringmanip::itos(currPop+1) + "." + 
          Stringmanip::itos(mod+1) + ".best";
          
        cout << "Writing best model file: " << currFileName << endl;
        outfile.open(currFileName.c_str(), ios::out);
        if(!outfile.is_open()){
            throw AthenaExcept(currFileName + " unable to open for writing best model");
        }
        bestSolution = pop[mod];
        
        outfile << "CV: " << currPop+1 << endl;
        outfile << "Model Rank: " << mod + 1 << endl;
        outfile << "Training result: " << bestSolution->fitness() << endl;
        outfile << "Testing result: " << bestSolution->testval() << endl;
        outfile << "Model:" << endl;
        bestSolution->output_clean(outfile, data, map_used, ott_dummy, continmap_used);
//         outfile << endl << endl <<"Grammar-compatible version:" << endl;
//         outfile << scaleInfo;
//         bestSolution->output_solution(outfile);
        
        outfile.close();
      }
//     }   
}


///
/// returns a stream for writing
/// @param filename
/// @return ostream to write to 
///
std::ostream& OutputManager::getStream(std::string filename){
  log_stream.open(filename.c_str(), ios::out);
  if(!log_stream.is_open()){
    throw AthenaExcept(filename + " unable to open for writing results");
  }
  
  return log_stream;
  
}


///
/// returns a stream for writing
/// @param filename
/// @return ostream to write to 
///
void OutputManager::outputGraphic(Algorithm* alg, Population& pop, int currPop, std::string basename,
  int nmodels, data_manage::Dataholder& data, bool map_used, bool ott_dummy, bool continmap_used){

  Solution* currSolution;

  // to create a png file use 
  // dot -Tpng <dot file> -o <outfile name>

  string ext = alg->getGraphicalFileExt();
//   for(unsigned int currPop=0; currPop < pops.size(); currPop++){
    for(int mod=0; mod < nmodels; mod++){
      ofstream outfile;
      string currFileName = basename + ".cv" + Stringmanip::itos(currPop+1) + "." + 
          Stringmanip::itos(mod+1) + ext;
      cout << "Writing file " << currFileName << endl;         
      outfile.open(currFileName.c_str(), ios::out);
      if(!outfile.is_open()){
          throw AthenaExcept(currFileName + " unable to open for writing best model");
      }  
      
      currSolution = pop[mod];
      alg->writeGraphical(outfile, currSolution, &data, map_used, ott_dummy, continmap_used);
      outfile.close();
    }
//   }
}


