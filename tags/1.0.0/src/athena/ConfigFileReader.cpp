
#include "Config.h"

#include "ConfigFileReader.h"
#include <sstream>
#include <iostream>

///
/// Contructor
///
ConfigFileReader::ConfigFileReader(string configfile){
    initialize_keywords();
    read_config(configfile);
}

///
/// Initialize keyword map
/// @return 
///
void ConfigFileReader::initialize_keywords(){
  
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
  keywordMap["TESTCONTINFILE"] = keyContinTestFile;
  keywordMap["TRAINCONTINFILE"] = keyContinTrainFile;
  keywordMap["IDINCLUDED"] = keyIDIncluded;
  keywordMap["DUMMYENCODE"] = keyDummyEncode;
  keywordMap["NUMSTEPS"]= keyNumExchanges;
  keywordMap["WRITECV"] = keyCVOutput;
  keywordMap["STATUSADJUST"] = keyStatusChange;
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
}


///
/// Reads configuration file and stores parameters
/// @param 
/// @return Config
///
Config ConfigFileReader::read_config(string configfile){
   Config configuration;
    
   ifstream configstream(configfile.c_str(), ios::in);
   if(!configstream.is_open()){
     throw AthenaExcept("Failed in attempt to open file " + configfile 
        + "\n");
   }
   
   
   string line, keyword, dataname, alg_name, data_type, output_name, map_name,
           continfile, id_included, value;
   int missingValue, randseed, num_cv, num_exchanges;
   float continMissValue, statusmissing;
   while(!configstream.eof()){
       getline(configstream, line);
       
       if(skip_line(line)){
           continue;
       }
       
       stringstream ss(line);
       
       ss >> keyword;
       keyword = Stringmanip::to_upper(keyword);      
       switch(keywordMap[keyword]){
           case keyNoMatch:
               throw AthenaExcept(keyword + " is not a valid keyword in configuration file");
               break;
           case keyDataset:
               ss >> dataname;             
               configuration.setDataSetName(dataname);
               break;
           case keyOut:
               ss >> output_name;
               configuration.setOutputName(output_name);
               break;
           case keyMapFile:
               ss >> map_name;
               configuration.setMapName(map_name);
               break;
           case keyMissingValue:
               ss >> missingValue;
               configuration.setMissingValue(missingValue, keyword);
               break;
           case keyStatusMissingValue:
               ss >> statusmissing;
               configuration.setStatusMissingValue(statusmissing);
               break;
           case keyRandSeed:
               ss >> randseed;
               configuration.setRandSeed(randseed);
               break;
           case keyCV:
               ss >> num_cv;
               configuration.setNumCV(num_cv, keyword);
               break;
           case keyEnd:
               throw(AthenaExcept(keyword + " is unmatched by ALGORITHM keyword"));
               break;
           case keyDatasetType:
               ss >> data_type;
               configuration.setDataSetType(data_type);
               break;
           case keyIDIncluded:
               {
               ss >> id_included;
               bool id = param_true(Stringmanip::to_upper(id_included));
               configuration.setIDinData(id);
               }
               break;
           case keyContinFile:
               ss >> continfile;
               configuration.setContinFileName(continfile);
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
               ss >> num_exchanges;
               configuration.setNumExchanges(num_exchanges, keyword);
               break;
           case keyAlgorithm:
               {
                  AlgorithmParams alg_param;
                  ss >> alg_param.name;
                  read_alg_params(alg_param, configstream);
                  configuration.addAlgorithmParam(alg_param);
               }
               break;
           case keyStatusChange:
             ss >> value;
             configuration.setStatusAdjust(Stringmanip::to_upper(value));
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
           default:
              throw AthenaExcept(keyword + " is not a valid keyword in configuration file");
              break;
       }
   }
   
   configstream.close();
   
   return configuration;
}


///
/// Checks to see if parameter is set to TRUE or ON
/// @param param string containing text to check
/// @return true when set to TRUE or ON, false otherwise
///
bool ConfigFileReader::param_true(string param){
   if(param.compare("TRUE") == 0 || param.compare("ON")==0)
       return true;
   return false;
}


///
/// Reads in algorithm parameters and
/// stores them in map for use by appropriate
/// algorithm
/// @param alg_apram AlgorithmParams
/// @param configstream 
///
void ConfigFileReader::read_alg_params(AlgorithmParams& alg_param, ifstream& configstream){
    string line, keyword;
    bool done = false;
    
    while(!configstream.eof() && !done){
        getline(configstream, line);
        if(skip_line(line)){
            continue;
        }
        
        // grab keyword
        stringstream ss(line+ " ");
        ss >> keyword;
        keyword = Stringmanip::to_upper(keyword);
        
        switch(keywordMap[keyword]){
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
                alg_param.params[keyword] = joined;
            } 
        }  
        
    }
}


///
/// Checks whether to skip blank lines or comment lines.
/// @param currLine string to check as comment or blank
/// @return true if skipping
///
bool ConfigFileReader::skip_line(string line){
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



