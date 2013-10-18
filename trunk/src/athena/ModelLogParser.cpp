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
#include <deque>
#include <stdio.h>

using namespace std;

bool compareLogModels(logModel first, logModel second){
	if(first.fitness < second.fitness)
		return true;
	else
		return false;
}



logModel ModelLogParser::getModel(string& line){

		logModel mod;    
		stringstream is(line);
		is >> mod.gen >> mod.rank >> mod.fitness >> mod.gramDepth >> mod.nnDepth 
				>> mod.n_g >> mod.n_c;
		getline(is, mod.model);
		return mod;
}


void ModelLogParser::compileVariableFiles(std::vector<std::string>& filenames, std::string outFilename){
		// when only a single file 
		// only need to change name of the filename to match the output name
		if(filenames.size()==1){
				string command = "mv " + filenames[0] + " " + outFilename;
				system(command.c_str());
				return;
		}
		
		// compile all files into a single file
		// open filehandle for each and then read each 
		deque<ifstream*> inputFH;
		
		vector<string>::iterator iter;
		for(iter=filenames.begin(); iter != filenames.end(); iter++){
			ifstream * newStream = new ifstream;
			newStream->open(iter->c_str(), ios::in);
			inputFH.push_back(newStream);
		}
		
		ofstream outStream(outFilename.c_str(), ios::out);
		
		bool done=false;
		vector<string> lines;
		string nextGen, gen = "";
		while(!done){
			lines.clear();
			for(deque<ifstream*>::iterator iter=inputFH.begin(); iter != inputFH.end(); iter++){
				done=getNextGen(*iter, lines, nextGen);
			}
			if(gen.length() > 0)
			outputVariables(outStream, lines, gen);	
			gen = nextGen;
		}

		for(deque<ifstream*>::iterator iter=inputFH.begin(); iter != inputFH.end(); iter++){
			(*iter)->close();
			delete *iter;
		}
		outStream.close();
		
		
		//delete temporary files
		for(iter=filenames.begin(); iter != filenames.end(); iter++){
			remove(iter->c_str());
		}
		
}


void ModelLogParser::outputVariables(ofstream& outStream, vector<string>& lines, string& gen){
	
	outStream << gen << "\n";
	for(vector<string>::iterator iter=lines.begin(); iter != lines.end(); iter++){
		outStream << *iter << "\n";
	}
	
}


bool ModelLogParser::getNextGen(ifstream* in, vector<string>& lines, string& nextGen){
	
	string line;
	bool endGen = false;
	
	while(!endGen){
		getline(*in, line);
		if(in->eof() || line.find("Gen")!=string::npos){
			nextGen=line;
			endGen=true;
		}
		else{
			lines.push_back(line);
		}
		
	}
	
	if(in->eof())
		return true;
	else
		return false;
}


void ModelLogParser::compileFiles(vector<string>& filenames, string outFilename, float notValid){

		// when only a single file 
		// only need to change name of the filename to match the output name
		if(filenames.size()==1){
				string command = "mv " + filenames[0] + " " + outFilename;
				system(command.c_str());
				return;
		}
		else
			return;
		
		vector<vector<logModel> > allModels;

		// alternatively pull together all models at each generation,
		// sort them and then write them out to the new file
		vector<string>::iterator fileIter;
		for(fileIter=filenames.begin(); fileIter != filenames.end(); ++fileIter){
				// add an additional generation when needed
				parseFile(*fileIter, allModels);
		}
		
		// sort each generation by fitness
		for(vector<vector<logModel> >::iterator iter=allModels.begin(); iter != allModels.end(); ++iter){
				sort(iter->begin(), iter->end(), compareLogModels);
		}
		
		// output the sorted models to the combined file
		ofstream outStream(outFilename.c_str(), ios::out);
		
		writeOutput(outStream, allModels, notValid);   
}



void ModelLogParser::parseFile(string filename, vector<vector<logModel> >& models){

		 ifstream logStream(filename.c_str(), ios::in);
		 if(!logStream.is_open()){
				throw AthenaExcept("Error:  Unable to open " + filename + "\n");
		 }
		 // skip the header line
		 string header;
		 getline(logStream, header);
		 vector<logModel> temp;
		 
		 // get all the models
		 while(!logStream.eof()){
				string line;
				getline(logStream, line);
				logModel logModel = getModel(line);
				if(int(models.size()) < logModel.gen+1){
						models.push_back(temp);
				}
				models[logModel.gen].push_back(logModel);
		 }
}


void ModelLogParser::writeOutput(ostream & os, vector<vector<logModel> >& models,
		float notValid){
		os << "GEN\tRANK\tFITNESS\tGRAM_DEPTH\tNN_DEPTH\tNUM_G\tNUM_C\tMODEL\n";
		int gen=0;

		for(vector<vector<logModel> >::iterator iter=models.begin(); iter != models.end(); ++iter){
				int rank = 1;
				for(vector<logModel>::iterator modIter=iter->begin(); modIter != iter->end(); ++modIter){
						os << gen << "\t" << rank <<  "\t";
						if(modIter->fitness == notValid)
								os << "NA";
						else
								os << modIter->fitness;
						os << "\t" << modIter->gramDepth <<
							"\t" << modIter->nnDepth << "\t" << modIter->n_g << "\t" << modIter->n_c << "\t" 
							<< modIter->model << endl;
						rank++;
				}
				gen++;
		}
}



