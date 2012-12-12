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
//ModelLogParser.cpp
#include "ModelLogParser.h"
#include <fstream>
#include <algorithm>
#include <sstream>

using namespace std;

bool compareLogModels(logModel first, logModel second){
  if(first.fitness < second.fitness)
    return true;
  else
    return false;
}

logModel ModelLogParser::get_model(string& line){

    logModel mod;    
    stringstream is(line);
    is >> mod.gen >> mod.rank >> mod.fitness >> mod.gram_depth >> mod.nn_depth 
        >> mod.n_g >> mod.n_c;
    getline(is, mod.model);
    return mod;
}

void ModelLogParser::compile_files(vector<string>& filenames, string outfilename, float not_valid){

    // when only a single file 
    // only need to change name of the filename to match the output name
    if(filenames.size()==1){
        string command = "mv " + filenames[0] + " " + outfilename;
        system(command.c_str());
        return;
    }
    else
    	return;
    
    vector<vector<logModel> > all_models;

    // alternatively pull together all models at each generation,
    // sort them and then write them out to the new file
    vector<string>::iterator file_iter;
    for(file_iter=filenames.begin(); file_iter != filenames.end(); ++file_iter){
        // add an additional generation when needed
        parse_file(*file_iter, all_models);
    }
    
    // sort each generation by fitness
    for(vector<vector<logModel> >::iterator iter=all_models.begin(); iter != all_models.end(); ++iter){
        sort(iter->begin(), iter->end(), compareLogModels);
    }
    
    // output the sorted models to the combined file
    ofstream out_stream(outfilename.c_str(), ios::out);
    
    write_output(out_stream, all_models, not_valid);
    
}


void ModelLogParser::parse_file(string filename, vector<vector<logModel> >& models){

     ifstream log_stream(filename.c_str(), ios::in);
     if(!log_stream.is_open()){
        throw AthenaExcept("Error:  Unable to open " + filename + "\n");
     }
     // skip the header line
     string header;
     getline(log_stream, header);
     vector<logModel> temp;
     
     // get all the models
     while(!log_stream.eof()){
        string line;
        getline(log_stream, line);
        logModel log_model = get_model(line);
        if(int(models.size()) < log_model.gen+1){
            models.push_back(temp);
        }
        models[log_model.gen].push_back(log_model);
     }
}


void ModelLogParser::write_output(ostream & os, vector<vector<logModel> >& models,
    float not_valid){
    os << "GEN\tRANK\tFITNESS\tGRAM_DEPTH\tNN_DEPTH\tNUM_G\tNUM_C\tMODEL\n";
    int gen=0;

//    os >> gen >> rank >> fitness >> gram_depth >> nn_depth >> n_g >> n_c;
//    getline(os, mod.model);
    
    for(vector<vector<logModel> >::iterator iter=models.begin(); iter != models.end(); ++iter){
        int rank = 1;
        for(vector<logModel>::iterator mod_iter=iter->begin(); mod_iter != iter->end(); ++mod_iter){
            os << gen << "\t" << rank <<  "\t";
            if(mod_iter->fitness == not_valid)
                os << "NA";
            else
                os << mod_iter->fitness;
            os << "\t" << mod_iter->gram_depth <<
              "\t" << mod_iter->nn_depth << "\t" << mod_iter->n_g << "\t" << mod_iter->n_c << "\t" 
              << mod_iter->model << endl;
            rank++;
        }
        gen++;
    }
    

}



