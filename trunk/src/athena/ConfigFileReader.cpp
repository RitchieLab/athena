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

#include "Config.h"

#include "ConfigFileReader.h"
#include <sstream>
#include <iostream>

///
/// Contructor
///
ConfigFileReader::ConfigFileReader(string configFile){
		initializeKeywords();
		readConfig(configFile);
}

///
/// Initialize keyWord map
/// @return 
///
void ConfigFileReader::initializeKeywords(){
	
	keywordMap["ALGORITHM"] = keyAlgorithm;
	keywordMap["END"] = keyEnd;
	keywordMap["DATASET"] = keyDataset;
	keywordMap["CV"] = keyCV;
	keywordMap["RANDSEED"]= keyRandSeed;
	keywordMap["OUT"] = keyOut;
	keywordMap["MISSINGVALUE"] = keyMissingValue;
	keywordMap["CONTINMISS"] = keyContinMiss;
	keywordMap["INPUT"] = keyDatasetType;
	keywordMap["MAPFILE"] = keyMapFile;
	keywordMap["CONTINFILE"] = keyContinFile;
	keywordMap["CONTINMAP"] = keyContinMapFile;
	keywordMap["TESTCONTINFILE"] = keyContinTestFile;
	keywordMap["TRAINCONTINFILE"] = keyContinTrainFile;
	keywordMap["IDINCLUDED"] = keyIDIncluded;
	keywordMap["DUMMYENCODE"] = keyDummyEncode;
	keywordMap["NUMSTEPS"]= keyNumExchanges;
	keywordMap["WRITECV"] = keyCVOutput;
	keywordMap["STATUSADJUST"] = keyStatusChange;
	keywordMap["CONTINADJUST"] = keyContinChange;
	keywordMap["INDOUTPUT"] = keyBestModelIndOutput;
	keywordMap["ALLNODESBEST"] = keyOutputAllNodeBest;
	keywordMap["TRAINFILE"] = keyTrainFile;
	keywordMap["TESTFILE"] = keyTestFile;
	keywordMap["BIOFILTERFILE"] = keyBioFilterFile;
	keywordMap["SUMMARYONLY"] = keySummaryOnly;
	keywordMap["STATUSMISSINGVALUE"] = keyStatusMissingValue;
	keywordMap["BIOFILETYPE"] = keyBioFileType;
	keywordMap["BIOGENEFILE"] = keyBioGeneFile;
	keywordMap["BIOARCHIVEFILE"] = keyBioArchiveFile;
	keywordMap["LOG"] = keyLogType;
	keywordMap["SPLITFILE"] = keySplitFile;
	keywordMap["CVSTART"]=keyCVStart;
	keywordMap["BESTSELECT"]=keySelectBestModel;
	keywordMap["IMAGEWRITER"]=keyImgWriter;
}



///
/// Reads configuration file and stores parameters
/// @param configFile
/// @return Config
///
Config ConfigFileReader::readConfig(string configFile){
	 Config configuration;
		
	 ifstream configStream(configFile.c_str(), ios::in);
	 if(!configStream.is_open()){
		 throw AthenaExcept("Failed in attempt to open file " + configFile 
				+ "\n");
	 }
	 
	 string line, keyWord, dataName,  DataType, outputName, mapName,
					 continFile, idIncluded, value;
	 int missingValue, randSeed, numCV, numExchanges;
	 float continMissValue, statusMissing;
	 while(!configStream.eof()){
			 getline(configStream, line);
			 
			 if(skipLine(line)){
					 continue;
			 }
			 
			 stringstream ss(line);
			 
			 ss >> keyWord;
			 keyWord = Stringmanip::to_upper(keyWord);      
			 switch(keywordMap[keyWord]){
					 case keyNoMatch:
							 throw AthenaExcept(keyWord + " is not a valid keyWord in configuration file");
							 break;
					 case keyDataset:
							 ss >> dataName;             
							 configuration.setDataSetName(dataName);
							 break;
					 case keyOut:
							 ss >> outputName;
							 configuration.setOutputName(outputName);
							 break;
					 case keyMapFile:
							 ss >> mapName;
							 configuration.setMapName(mapName);
							 break;
					 case keyContinMapFile:
							 ss >> mapName;
							 configuration.setContinMapName(mapName);
							 break;
					 case keyMissingValue:
							 ss >> missingValue;
							 configuration.setMissingValue(missingValue, keyWord);
							 break;
					 case keyStatusMissingValue:
							 ss >> statusMissing;
							 configuration.setStatusMissingValue(statusMissing);
							 break;
					 case keyRandSeed:
							 ss >> randSeed;
							 configuration.setRandSeed(randSeed);
							 break;
					 case keyCV:
							 ss >> numCV;
							 configuration.setNumCV(numCV, keyWord);
							 break;
					 case keyEnd:
							 throw(AthenaExcept(keyWord + " is unmatched by ALGORITHM keyWord"));
							 break;
					 case keyDatasetType:
							 ss >> DataType;
							 configuration.setDataSetType(DataType);
							 break;
					 case keyIDIncluded:
							 {
							 ss >> idIncluded;
							 bool id = paramTrue(Stringmanip::to_upper(idIncluded));
							 configuration.setIDinData(id);
							 }
							 break;
					 case keyContinFile:
							 ss >> continFile;
							 configuration.setContinFileName(continFile);
							 break;
					 case keyContinMiss:
							 ss >> continMissValue;
							 configuration.setContinMiss(continMissValue);
							 break;
					 case keyDummyEncode:
							 ss >> value;
							 configuration.setEncodeType(Stringmanip::to_upper(value));
							 break;
					 case keyNumExchanges:
							 ss >> numExchanges;
							 configuration.setNumExchanges(numExchanges, keyWord);
							 break;
					 case keyAlgorithm:
							 {
									AlgorithmParams algParam;
									ss >> algParam.name;
									readAlgParams(algParam, configStream);
									configuration.addAlgorithmParam(algParam);
							 }
							 break;
					 case keyStatusChange:
						 ss >> value;
						 configuration.setStatusAdjust(Stringmanip::to_upper(value));
						 break;
					 case keyContinChange:
					 	 ss >> value;
						 configuration.setContinAdjust(Stringmanip::to_upper(value));
						 break;
					 case keyCVOutput:
						 ss >> value;
						 configuration.setCVOutput(Stringmanip::check_true_false(value));
						 break;
					 case keyBestModelIndOutput:
						 ss >> value;
						 configuration.setIndOutput(Stringmanip::check_true_false(value));
						 break;
					 case keyOutputAllNodeBest:
						 ss >> value;
						 configuration.setOutputAllNodesBest(Stringmanip::check_true_false(value));
						 break;
					 case keySummaryOnly:
						 ss >> value;
						 configuration.setSummaryOnly(value);
						 break;
					 case keyTestFile:
						 ss >> value;
						 configuration.setTestFile(value);
						 break;
					 case keyTrainFile:
						 ss >> value;
						 configuration.setTrainFile(value);
						 break;
					 case keyContinTestFile:
							ss >> value;
							configuration.setContinTestFile(value);
							break;
					 case keyContinTrainFile:
							ss >> value;
							configuration.setContinTrainFile(value);
							break;
					 case keyBioFilterFile:
							ss >> value;
							configuration.setBioFilterFile(value);
							break;
					 case keyBioFileType:
							ss >> value;
							configuration.setBioFileType(value);
							break;
					 case keyBioArchiveFile:
							ss >> value;
							configuration.setBioArchiveFile(value);
							break;
					 case keyBioGeneFile:
							ss >> value;
							configuration.setBioGeneFile(value);
							break;
					 case keyLogType:
							ss >> value;
							configuration.setLogType(Stringmanip::to_upper(value));
							break;
					 case keySplitFile:
							ss >> value;
							configuration.setSplitFile(value);
							break;
					 case keyCVStart:
					 	  ss >> value;
					 	  configuration.setStartCV(Stringmanip::stringToNumber<int>(value));
							break;
					 case keySelectBestModel:
					 		ss >> value;
					 		configuration.setSelectBestModel(Stringmanip::check_true_false(value));
					 		break;
					 case keyImgWriter:
					 		ss >> value;
					 		configuration.setImgWriter(value);
					 		break;
					 default:
							throw AthenaExcept(keyWord + " is not a valid keyWord in configuration file");
							break;
			 }
	 }
	 
	 configStream.close();
	 
	 return configuration;
}


///
/// Checks to see if parameter is set to TRUE or ON
/// @param param string containing text to check
/// @return true when set to TRUE or ON, false otherwise
///
bool ConfigFileReader::paramTrue(string param){
	 if(param.compare("TRUE") == 0 || param.compare("ON")==0)
			 return true;
	 return false;
}



///
/// Reads in algorithm parameters and
/// stores them in map for use by appropriate
/// algorithm
/// @param algParam AlgorithmParams
/// @param configStream 
///
void ConfigFileReader::readAlgParams(AlgorithmParams& algParam, ifstream& configStream){
		string line, keyWord;
		bool done = false;
		
		while(!configStream.eof() && !done){
				getline(configStream, line);
				if(skipLine(line)){
						continue;
				}
				
				// grab keyWord
				stringstream ss(line+ " ");
				ss >> keyWord;
				keyWord = Stringmanip::to_upper(keyWord);
				
				switch(keywordMap[keyWord]){
						case keyEnd:
								return;  // done reading parameters for this alogrithm
						default:
								// add parameter to this algorithm for its own use
						{
								string holder, joined;
								bool start = true;
								ss >> holder;
								do{
										if(start == false){
												joined += " ";
										}
										joined += holder;
										start = false;
										ss >> holder;
								}while(!ss.eof());
								algParam.params[keyWord] = joined;
						} 
				}  
				
		}
}



///
/// Checks whether to skip blank lines or comment lines.
/// @param currLine string to check as comment or blank
/// @return true if skipping
///
bool ConfigFileReader::skipLine(string line){
		// skip blank lines
		if(line.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") 
						== string::npos){
			return true;
		} 
		
		int character = 0;
		// skip comments -- begin line with '#' 
		while(line[character] == ' ' || line[character] == '\t'){
				character++;
		}
		
		if(line[character] == '#'){
				return true;
		}

		return false;
	 
}



