#include "OutputManager.h"
#include <Stringmanip.h>
#include <iomanip>
#include <sstream>

using namespace std;


///
/// outputs summary of best models
/// @param pops vector containing Populations of solutions to summarize
/// @param data Dataholder that can translate the snps back to original IDs
/// @param mapfile_used true when an actual map file was used and there are original
/// names to output
/// @param dummy_encoded when true genotypes need to be adjusted back to reflect
/// original genotype positions
///
void OutputManager::outputSummary(vector<Population>& pops,
  data_manage::Dataholder& data,  bool mapfile_used, bool dummy_encoded){
    
    string summaryName = basename + ".athena.sum";
    
    ofstream outfile;
    outfile.open(summaryName.c_str(), ios::out);
   
    Solution* bestSolution;
    
    string prefix;
    int width = 30;
    if(!mapfile_used){
      prefix = "G";
      width = 20;
    }
    
    outfile << setw(5) << left << "CV" << setw(width) << "Variables" << " " << setw(10) << "Training"
            << " " << setw(10) << "Testing" << endl;    
     
    for(unsigned int currPop=0; currPop < pops.size(); currPop++){
        bestSolution = pops[currPop].best();
        
        vector<int> genos = bestSolution->get_genotypes(dummy_encoded);
        vector<int> covars = bestSolution->get_covariates();
        
        outfile << setw(5) << currPop+1;
        stringstream ss;
        for(unsigned int g=0; g < genos.size(); g++){
            ss << prefix << data.get_geno_name(genos[g]-1) << " ";
        }
        stringstream cs;    
        for(unsigned int c=0; c < covars.size(); c++){
            cs << "C"  << data.get_covar_name(covars[c]-1) << " ";
        }
        
        outfile << setw(width) << ss.str() + cs.str() << " ";
        
        outfile << setw(10) << bestSolution->fitness() << " ";
        outfile << setw(10) << bestSolution->testval();
        outfile << endl;
    }
    
    outfile.close();
    
}



///
/// outputs a file for each best model
/// @param pops Population vector
/// @param nmodels Number of models to output
/// @param scaleInfo Informtion on the scaling done in model
/// @param data
/// @param map_used
///
void OutputManager::outputBestModels(vector<Population>& pops, int nmodels,
  string scaleInfo, data_manage::Dataholder& data, bool map_used, bool ott_dummy){
    
    Solution* bestSolution;   
    for(unsigned int currPop=0; currPop < pops.size(); currPop++){
      for(int mod=0; mod < nmodels; mod++){
        ofstream outfile;
        
        string currFileName = basename + ".cv" + Stringmanip::itos(currPop+1) + "." + 
          Stringmanip::itos(mod+1) + ".best";
          
        cout << "Writing best model file: " << currFileName << endl;
        outfile.open(currFileName.c_str(), ios::out);
        if(!outfile.is_open()){
            throw HemannExcept(currFileName + " unable to open for writing best model");
        }
        bestSolution = pops[currPop][mod];
        
        outfile << "CV: " << currPop+1 << endl;
        outfile << "Model Rank: " << mod + 1 << endl;
        outfile << "Training result: " << bestSolution->fitness() << endl;
        outfile << "Testing result: " << bestSolution->testval() << endl;
        outfile << "Model:" << endl;
        bestSolution->output_clean(outfile, data, map_used, ott_dummy);
        outfile << endl << endl <<"Grammar-compatible version:" << endl;
        outfile << scaleInfo;
        bestSolution->output_solution(outfile);
        
        outfile.close();
      }
    }   
}


///
/// returns a stream for writing
/// @param filename
/// @return ostream to write to 
///
std::ostream& OutputManager::getStream(std::string filename){
  log_stream.open(filename.c_str(), ios::out);
  if(!log_stream.is_open()){
    throw HemannExcept(filename + " unable to open for writing results");
  }
  
  return log_stream;
  
}


///
/// returns a stream for writing
/// @param filename
/// @return ostream to write to 
///
void OutputManager::outputGraphic(Algorithm* alg, std::vector<Population>& pops, std::string basename,
  int nmodels, data_manage::Dataholder& data, bool map_used, bool ott_dummy){

  Solution* currSolution;

  // to create a png file use 
  // dot -Tpng <dot file> -o <outfile name>

  string ext = alg->getGraphicalFileExt();
  for(unsigned int currPop=0; currPop < pops.size(); currPop++){
    for(int mod=0; mod < nmodels; mod++){
      ofstream outfile;
      string currFileName = basename + ".cv" + Stringmanip::itos(currPop+1) + "." + 
          Stringmanip::itos(mod+1) + ext;
      cout << "Writing file " << currFileName << endl;         
      outfile.open(currFileName.c_str(), ios::out);
      if(!outfile.is_open()){
          throw HemannExcept(currFileName + " unable to open for writing best model");
      }  
      
      currSolution = pops[currPop][mod];
      alg->writeGraphical(outfile, currSolution, &data, map_used, ott_dummy);
      outfile.close();
    }
  }
}


