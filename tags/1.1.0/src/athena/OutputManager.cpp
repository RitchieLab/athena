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
#include <set>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <algorithm>
// #include <deque>

using namespace std;


///
/// Sets up new files with appropriate headers for output
///
void OutputManager::setFiles(bool mapFileUsed, string fitnessName,
	std::vector<std::string> additionalHeaders){

	addHeaders = additionalHeaders;
	string summaryName = getSummaryFileName();
	string progressFileName = getProgressFileName();

	// create blank progress file as this is a run from start
	ofstream progOut;
	progOut.open(progressFileName.c_str(), ios::out);
	progOut.close();

	int width = 30;
		if(!mapFileUsed){
			width = 20;
		}

	ofstream outfile;
	outfile.open(summaryName.c_str(), ios::out);
		outfile << "CV\tVariables\t"  << fitnessName + " Training\tTesting";
	for(unsigned int i=0; i<additionalHeaders.size(); i++){
		outfile << "\tTraining-" << additionalHeaders[i];
	}
	for(unsigned int i=0; i<additionalHeaders.size(); i++){
		outfile <<  "\tTesting-" << additionalHeaders[i];
	}

	outfile << endl;
	outfile.close();
}

///
/// outputs summary of best models
/// @param pop Population containing best model
/// @param data Dataholder that can translate the snps back to original IDs
/// @param mapFileUsed true when an actual map file was used and there are original
/// names to output
/// @param dummyEncoded when true genotypes need to be adjusted back to reflect
/// original genotype positions
///
void OutputManager::outputSummary(Population& pop, int currPop,
	data_manage::Dataholder& data,  Algorithm* alg, bool mapFileUsed,
	bool dummyEncoded, bool continMapUsed,
	std::string fitnessName){

		string summaryName = getSummaryFileName();
		string progressFileName = getProgressFileName();

		ofstream outfile;
		outfile.open(summaryName.c_str(), ios::app);

		Solution* bestSolution;

		string prefix, continPrefix;
		if(!mapFileUsed){
			prefix = "G";
		}
		if(!continMapUsed){
			continPrefix = "C";
		}

				bestSolution = pop.best();
				vector<int> genos = bestSolution->getGenotypes(dummyEncoded);
				vector<int> covars = bestSolution->getCovariates();

				outfile << setw(5) << left << currPop+1;
				stringstream ss;
				for(unsigned int g=0; g < genos.size(); g++){
						ss << prefix << data.getGenoName(genos[g]-1) << " ";
				}
				stringstream cs;
				for(unsigned int c=0; c < covars.size(); c++){
						cs << continPrefix << data.getCovarName(covars[c]-1) << " ";
				}

		outfile << "\t" << ss.str() + cs.str() << "\t" << bestSolution->fitness() <<
			"\t" << bestSolution->testVal();

		for(std::vector<string>::iterator iter=bestSolution->getAdditionalOutput().begin();
			iter != bestSolution->getAdditionalOutput().end();
			++iter){
			outfile << "\t" << *iter;
		}
		outfile << endl;
		outfile.close();

		// write to progress file
		ofstream progFile;
		progFile.open(progressFileName.c_str(), ios::app);

		// each model will take up 2 lines with the first being the equation and the
		// second being the internal representation
		progFile << currPop+1 << "\t";
		alg->writeEquation(progFile, bestSolution, &data,
			mapFileUsed, dummyEncoded, continMapUsed);
		progFile << endl;
		progFile << currPop+1 << "\t";
		bestSolution->outputSolution(progFile);
		progFile.close();
}

///
/// Reads progress file for model information and inclusion in final summary file
///
void OutputManager::fillProgress(){
	string progressFileName = getProgressFileName();
	ifstream progFile;
	progFile.open(progressFileName.c_str(), ios::in);

	equationLines.clear();
	modelLines.clear();
	string line;
	while(getline(progFile, line)){
		equationLines.push_back(line);
		getline(progFile, line);
		modelLines.push_back(line);
	}

	progFile.close();
	remove(progressFileName.c_str());
}



void OutputManager::outputBest(Solution* bestSolution, data_manage::Dataholder& data,
	bool mapFileUsed, bool dummyEncoded, bool continMapUsed,
	std::string fitnessName){

		string filename = basename + ".overall.best";

		string prefix, continPrefix;
	 	if(!mapFileUsed){
			prefix = "G";
		}
		if(!continMapUsed){
			continPrefix = "C";
		}

		ofstream outfile;
		outfile.open(filename.c_str(), ios::out);
		vector<int> genos = bestSolution->getGenotypes(dummyEncoded);
		vector<int> covars = bestSolution->getCovariates();

		stringstream ss;
		for(unsigned int g=0; g < genos.size(); g++){
			ss << prefix << data.getGenoName(genos[g]-1) << " ";
		}
		stringstream cs;
		for(unsigned int c=0; c < covars.size(); c++){
			cs << continPrefix << data.getCovarName(covars[c]-1) << " ";
		}

		outfile << "Variables:\t" << ss.str() + cs.str() << endl;
		outfile << fitnessName << ":\t" << bestSolution->fitness() << endl;
		if(!addHeaders.empty()){
			for(unsigned int i=0; i<addHeaders.size(); i++){
				outfile << addHeaders[i] << ":\t" << bestSolution->getAdditionalOutput()[i] << endl;
			}
		}

		outfile << "\nModel:" << endl;
		bestSolution->outputClean(outfile, data, mapFileUsed, dummyEncoded, continMapUsed);
		outfile.close();

}



///
/// outputs a file for each best model
/// @param pops Population vector
/// @param nmodels Number of models to output
/// @param currPop population number matching the cross-validation
/// @param scaleInfo Information on the scaling done in model
/// @param data
/// @param mapUsed
///
void OutputManager::outputBestModels(Population& pop, int nmodels, int currPop,
	string scaleInfo, data_manage::Dataholder& data, bool mapUsed, bool ottDummy,
	bool continMapUsed){

		Solution* bestSolution;
			for(int mod=0; mod < nmodels; mod++){
				ofstream outfile;

				string currFileName = basename + ".cv" + Stringmanip::numberToString(currPop+1) + "." +
					Stringmanip::numberToString(mod+1) + ".best";

				cout << "Writing best model file: " << currFileName << endl;
				outfile.open(currFileName.c_str(), ios::out);
				if(!outfile.is_open()){
						throw AthenaExcept(currFileName + " unable to open for writing best model");
				}
				bestSolution = pop[mod];

				outfile << "CV: " << currPop+1 << endl;
				outfile << "Model Rank: " << mod + 1 << endl;
				outfile << "Training result: " << bestSolution->fitness() << endl;
				outfile << "Testing result: " << bestSolution->testVal() << endl;
				if(!addHeaders.empty()){
					for(unsigned int i=0; i<addHeaders.size(); i++){
						outfile << "Training-" << addHeaders[i] << ":\t" << bestSolution->getAdditionalOutput()[i] << endl;
					}
				}
				outfile << "Model:" << endl;
				bestSolution->outputClean(outfile, data, mapUsed, ottDummy, continMapUsed);

				outfile.close();
			}
}

///
/// Output all models for this population in one file
///
void OutputManager::outputAllModels(Algorithm* alg, Population& pop, int rank, int currPop,
	string scaleInfo, data_manage::Dataholder& data, bool mapUsed, bool ottDummy,
	bool continMapUsed, bool testingDone){
	string currFileName = basename + ".cv" + Stringmanip::numberToString(currPop+1) + "." +
		Stringmanip::numberToString(rank+1) + ".all";
	ofstream outfile;
	outfile.open(currFileName.c_str(), ios::out);
	if(!outfile.is_open()){
		throw AthenaExcept(currFileName + " unable to open for writing all models");
	}
	outfile << "All Models From Final Generation\n";
	outfile << "Model\tEquation\tTraining\tTesting";
	for(size_t i=0; i < addHeaders.size(); i++){
		outfile << "\tTraining-" << addHeaders[i];
	}
	if(testingDone){
		for(size_t i=0; i < addHeaders.size(); i++){
			outfile << "\tTesting-" << addHeaders[i];
		}
	}
	outfile << "\n";
	Solution* bestSolution;
	for(int mod=0; mod < pop.getPopSize(); mod++){
		bestSolution = pop[mod];
		bestSolution->outputClean(outfile, data, mapUsed, ottDummy, continMapUsed);
		outfile << "\t";
		try{
			alg->writeEquation(outfile, pop[mod], &data, mapUsed, ottDummy, continMapUsed);
		}catch(AthenaExcept& ae){
			outfile << "incomplete";
		}
		outfile << "\t" << bestSolution->fitness() << "\t" << bestSolution->testVal();
		for(size_t i=0; i < addHeaders.size(); i++){
			outfile << "\t";
			if(bestSolution->getAdditionalOutput().size() > i)
				outfile << bestSolution->getAdditionalOutput()[i];
		}
		if(testingDone){
			for(size_t i=0; i < addHeaders.size(); i++){
				outfile << "\t";
				if(bestSolution->getAdditionalOutput().size() > i + addHeaders.size())
					outfile << bestSolution->getAdditionalOutput()[i + addHeaders.size()];
			}
		}
	outfile << "\n";
	}
}


///
/// Renames all models file (for single) or combines all models for
/// multiple files
///
void OutputManager::combineAllModels(int nProcs, int currCV, Algorithm* alg){

	string finalName = basename + ".cv" + Stringmanip::numberToString(currCV+1) + ".all";
	string currFilename;
	int sortColumn = alg->allModelSortColumn();

	if(nProcs==1){
		currFilename =  basename + ".cv" + Stringmanip::numberToString(currCV+1) + ".1.all";
		string command = "mv " + currFilename + " " + finalName;
		system(command.c_str());
	}
	else{
		// parse each one
		string line, header1, header2;
		vector<string> lines;
		vector<ModelInfo> modInfo;
		ModelInfo m;
		std::vector<std::string> linePcs;
		for(int p = 1; p <= nProcs; p++){
			currFilename = basename + ".cv" + Stringmanip::numberToString(currCV+1) + "." +
				Stringmanip::numberToString(p) + ".all";
			ifstream infile;
			infile.open(currFilename.c_str());
			getline(infile, header1);
			getline(infile, header2);
			while(getline(infile, line)){
				lines.push_back(line);
				linePcs = Stringmanip::split(line, '\t');
				m.lineNo = lines.size()-1;
				m.score = Stringmanip::stringToNumber<float>(linePcs[sortColumn]);
				modInfo.push_back(m);
			}
			infile.close();
			remove(currFilename.c_str());
		}
		// sort the lines
		sort(modInfo.begin(), modInfo.end(), mySorter);

		ofstream outfile;
		outfile.open(finalName.c_str(), ios::out);
		outfile << header1 << "\n" << header2 << "\n";
		for(vector<ModelInfo>::iterator infoIter = modInfo.begin(); infoIter != modInfo.end();
			++infoIter){
			outfile << lines[infoIter->lineNo] << "\n";
		}
		outfile.close();

	}

}

///
/// returns a stream for writing
/// @param filename
/// @return ostream to write to
///
std::ostream& OutputManager::getStream(std::string filename){
	logStream.open(filename.c_str(), ios::out);
	if(!logStream.is_open()){
		throw AthenaExcept(filename + " unable to open for writing results");
	}

	return logStream;

}


///
/// Output validation result files
///
void OutputManager::writeValidation(string fitnessName, std::vector<std::string> additionalHeaders,
	vector<Solution*> models, data_manage::Dataholder& data, bool mapUsed, bool dummyEncoded,
	bool continMapUsed, Algorithm* alg){

	int width = 30;
	if(!mapUsed){
		width = 20;
	}

	string filename = basename + ".validation";
	ofstream outfile;
	outfile.open(filename.c_str(), ios::out);

	outfile << "Model\tVariables\t"  << fitnessName + " Training\tTesting";
	for(unsigned int i=0; i<additionalHeaders.size(); i++){
		outfile << "\tTraining-" << additionalHeaders[i];
	}
	for(unsigned int i=0; i<additionalHeaders.size(); i++){
		outfile <<  "\tTesting-" << additionalHeaders[i];
	}
	outfile << "\n";

	string prefix, continPrefix;
	if(!mapUsed){
		prefix = "G";
	}
	if(!continMapUsed){
		continPrefix = "C";
	}

	Solution* bestSolution;
	for(size_t i=0; i<models.size(); i++){
		bestSolution = models[i];

		vector<int> genos = bestSolution->getGenotypes(dummyEncoded);
		vector<int> covars = bestSolution->getCovariates();

		outfile << setw(5) << left << i+1;
		stringstream ss;
		for(unsigned int g=0; g < genos.size(); g++){
			ss << prefix << data.getGenoName(genos[g]-1) << " ";
		}
		stringstream cs;
		for(unsigned int c=0; c < covars.size(); c++){
			cs << continPrefix << data.getCovarName(covars[c]-1) << " ";
		}

		outfile << "\t" << ss.str() + cs.str() << "\t" << bestSolution->fitness() <<
			"\t" << bestSolution->testVal();

		for(std::vector<string>::iterator iter=bestSolution->getAdditionalOutput().begin();
			iter != bestSolution->getAdditionalOutput().end();
			++iter){
			outfile << "\t" << *iter;
		}
		outfile << endl;
	}

	outfile << "\nModel\tEquation\n";

	for(size_t i=0; i < models.size(); i++){
		outfile << i+1 << "\t";
		alg->writeEquation(outfile, models[i], &data,
			mapUsed, dummyEncoded, continMapUsed);
		outfile << endl;
	}

	outfile << "\n\n**** For use by ATHENA when running models with independent datasets ****\n";
	outfile << "\nModel\tInternal ATHENA representation\n";
	for(size_t i=0; i<models.size(); i++){
			outfile << i+1 << "\t";
			models[i]->outputSolution(outfile);
	}

	outfile.close();
}


///
/// returns a stream for writing
/// @param filename
/// @return ostream to write to
///
void OutputManager::outputEquations(Algorithm* alg, vector<Solution*>& bestSolutions,
			data_manage::Dataholder& data, bool mapUsed, bool ottDummy,
			bool continMapUsed){

	string summaryName = getSummaryFileName();
	fillProgress();
	ofstream outfile;
	outfile.open(summaryName.c_str(), ios::app);
	outfile << "\nCV\tModel\n";

	vector<string>::iterator strIter;
	for(strIter=equationLines.begin(); strIter != equationLines.end(); ++strIter){
		outfile << *strIter << "\n";
	}

	outfile << "\n\n**** For use by ATHENA when running models with independent datasets ****\n";
	outfile << "\nCV\tInternal ATHENA representation\n";
	for(strIter=modelLines.begin(); strIter != modelLines.end(); ++strIter){
		outfile << *strIter << "\n";
	}

	outfile.close();

}


///
/// Outputs graphical representation of model
/// @param filename
/// @return ostream to write to
///
void OutputManager::outputGraphic(Algorithm* alg, Population& pop, int currPop, std::string basename,
	int nmodels, data_manage::Dataholder& data, bool mapUsed, bool ottDummy,
	bool continMapUsed, std::string imgWriter){

	Solution* currSolution;

	string ext = alg->getGraphicalFileExt();
	if(ext.length() > 0){
		for(int mod=0; mod < nmodels; mod++){
			ofstream outfile;
			string imgFileBaseName = basename + ".cv" + Stringmanip::numberToString(currPop+1) + "." +
					Stringmanip::numberToString(mod+1);
			string currFileName = imgFileBaseName + ext;
			cout << "Writing file " << currFileName << endl;
			outfile.open(currFileName.c_str(), ios::out);
			if(!outfile.is_open()){
					throw AthenaExcept(currFileName + " unable to open for writing " + ext + " file.");
			}

			currSolution = pop[mod];
			alg->writeGraphical(outfile, currSolution, &data, mapUsed, ottDummy, continMapUsed);
			outfile.close();
			if(imgWriter.length() > 1)
				alg->produceGraphic(currFileName, imgFileBaseName, imgWriter);
		}
	}
}


///
/// Converts scores for results using balanced accuracy
/// @param score original output score
/// @return status as 1 or 0
///
int OutputManager::scoreConversion(float score){
	return score > 0.5?1:0;
}

///
/// outputs list of individuals with scores for best model
/// @param is istream
/// @param base base portion of output file name
///
void OutputManager::validationIndOutput(vector<std::stringstream*>& ss, std::string base){

	for(size_t i=0; i<ss.size(); i++){
		string currFileName = base + ".indscores." + Stringmanip::numberToString(i+1) + ".validation.txt";
		ofstream of(currFileName.c_str());
		of << "IND ID\tScore\tObserved\n";
		of << ss[i]->str() << endl;
		of.close();
	}
}


///
/// outputs list of individuals with scores for best model
/// divided into training and testing sets
/// @param is istream
/// @param base base portion of output file name
///
void OutputManager::outputInds(std::istream &is, std::string base, string fitnessName){

	std::multiset<indOutScores,scoreComp> empty;

	vector<multiset<indOutScores,scoreComp> > trainingSet;
	vector<multiset<indOutScores,scoreComp> > testingSet;
	multiset<indOutScores,scoreComp> bestSet;
	int currTrain=0, currTest=0;

	bool balacc=false;
	if(fitnessName.find("BALANCEDACC") != string::npos)
		balacc = true;

	string line;
	string temp, id, scoreStr, finalScore, diffStr, predStatus="";
	int cv;
	bool useTrain=true, useTest=false;
	double score, obs, diff, missDiffValue=-1000.0;

	while(getline(is, line)){
		stringstream ss(line);
		if(line.find("CV") != string::npos){
			if(line.find("Best") == string::npos){
				ss >> temp >> cv;
				currTrain = cv-1;
				currTest=cv-1;
			}
			else{
				useTrain=false;
				useTest=false;
			}
		}
		else if(line.find("Train") != string::npos){
			useTrain=true;
			useTest=false;
			trainingSet.push_back(empty);
		}
		else if(line.find("Test") != string::npos){
			useTest=true;
			useTrain=false;
			testingSet.push_back(empty);
		}
		else{
			ss >> id >> scoreStr >> obs;
			if(scoreStr.find("Missing") != string::npos){
				diffStr = "NA";
				diff = missDiffValue;
				if(balacc)
					predStatus="\tMissing";
			}
			else{
				stringstream ss(scoreStr);
				ss >> score;
				if(!balacc)
					diff = abs(score-obs);
				else
					diff = int(abs(scoreConversion(score)-obs));
				diffStr = Stringmanip::numberToString(diff);
				if(balacc)
					predStatus = "\t" + Stringmanip::numberToString(scoreConversion(score));
			}

			indOutScores newOutput;
			newOutput.output = id + "\t" + scoreStr + predStatus + "\t" + Stringmanip::numberToString(obs);

			newOutput.diff = diff;

			if(useTrain){
					trainingSet[currTrain].insert(newOutput);
			}
			else if(useTest){
					testingSet[currTest].insert(newOutput);
			}
			else{
					bestSet.insert(newOutput);
			}
		}
	}

	vector<string> best, emptyVec;
	vector<vector<string> > training, testing;
	std::multiset<indOutScores,scoreComp>::iterator setIter;
	for(size_t i=0; i<trainingSet.size(); i++){
		training.push_back(emptyVec);
		for(setIter = trainingSet[i].begin(); setIter != trainingSet[i].end(); ++setIter){
			training[i].push_back(setIter->output);
		}
	}
	for(size_t i=0; i<testingSet.size(); i++){
		testing.push_back(emptyVec);
		for(setIter = testingSet[i].begin(); setIter != testingSet[i].end(); ++setIter){
			testing[i].push_back(setIter->output);
		}
	}
	for(setIter=bestSet.begin(); setIter!=bestSet.end(); ++setIter){
		best.push_back(setIter->output);
	}

	string addHeader="", trainTestHeader="";
	if(balacc){
		addHeader="\tPred-Status";
		trainTestHeader="\t";
	}

	// write to output file
	string currFileName = base + ".indscores.txt";
	ofstream of(currFileName.c_str());

	cout << "Writing individual output file: " << currFileName << endl;

	if(best.empty()){
		of << "Training\t\t\t" + trainTestHeader + "Testing";
		for(int j=1; j<cv; j++){
			of << "\t\t\t" + trainTestHeader + "Training\t\t\t" + trainTestHeader+ "Testing";
		}
		of << endl;
		of << "CV#1-ID\tPred" + addHeader + "\tObs\tCV#1-ID\tPred" + addHeader + "\tObs";
		for(int j=1; j<cv; j++){
			of << "\tCV#" << j+1 << "-ID\tPred" + addHeader + "\tObs\tCV#" << j+1 << "-ID\tPred"
				+ addHeader + "\tObs";
		}
		of << endl;

		for(unsigned int i=0; i<training[0].size(); i++){
			for(int j=0; j<cv; j++){
					if(i < training[j].size()){
						of << training[j][i] << "\t";
					}
					else{
						of << "\t\t";
						if(balacc)
							of << "\t";
					}
					if(!testing.empty() && i < testing[j].size()){
						of << testing[j][i] << "\t";
					}
					else{
						of << "\t\t\t";
						if(balacc)
							of << "\t";
					}
				}
			of << endl;
			}
	}
	else{

		of << "Entire Set";
		for(int j=1; j<cv+1; j++){
			of << "\t\t\tTraining\t\t\t" + trainTestHeader + "Testing";
		}
		of << endl;
		of << "Best-ID\tPred" + addHeader + "\tObs";
		for(int j=0; j<cv; j++){
			of << "\tCV#" << j+1 << "-ID\tPred" + addHeader + "\tObs\tCV#" << j+1 << "-ID\tPred" + addHeader + "\tObs";
		}
		of << endl;

		for(unsigned int i=0; i<best.size(); i++){
			of << best[i] << "\t";
			for(int j=0; j<cv; j++){
				if(i < training[j].size()){
					of << training[j][i] << "\t";
				}
				else{
					of << "\t\t\t";
					if(balacc)
						of << "\t";
				}
				if(!testing.empty() && i < testing[j].size()){
					of << testing[j][i] << "\t";
				}
				else{
					of << "\t\t\t";
					if(balacc)
						of << "\t";
				}
			}
			of << endl;
		}
	}
	of.close();
}
