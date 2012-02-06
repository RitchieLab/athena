/* 
 * File:   Config.h
 * Author: dudeksm
 *
 * Created on November 5, 2008, 4:59 PM
 */

#ifndef _CONFIG_H
#define	_CONFIG_H

#include <map>
#include <string>
#include <vector>
#include "AthenaExcept.h"
#include "Structs.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

/// holds parameters and name of the Algorithm
struct AlgorithmParams{
    std::map<std::string, std::string> params;
    std::string name;
};

///
/// Contains configuration parameters for running HEMANN
///
class Config{
    
public:
    /// Constructor
    Config();
    
    std::string getDataSetType(){return data_type;}
    void setDataSetType(std::string dtype){data_type = dtype;}

    std::string getOutputName(){return out_name;}
    void setOutputName(std::string oname){out_name = oname;}
    
    std::string getMapName(){return map_name;}
    void setMapName(std::string mname){map_name = mname;}
    
    int getMissingValue(){return miss_value;}
    void setMissingValue(int mval, std::string id=""){
      if(mval >= 0 && mval <= 2)
        throw AthenaExcept("The missing value " + id +" cannot be between 0 and 2 inclusive.");
      miss_value = mval;
    }
    
    float getStatusMissingValue(){return stat_miss_value;}
    void setStatusMissingValue(float val){stat_miss_value = val;}
    
    int getRandSeed(){return rand_seed;}
    void setRandSeed(int rseed){rand_seed = rseed;}
    
    int getNumCV(){return ncv;}
    void setNumCV(int numcv, std::string id=""){
      if(numcv < 1)
        throw AthenaExcept("Number of cv " + id + " must be greater than zero");
      ncv = numcv;
    }
    
    std::string getDataSetName(){return datafile;}
    void setDataSetName(std::string dname){datafile = dname;}
    
    void addAlgorithmParam(AlgorithmParams alg){alg_params.push_back(alg);}
    std::vector<AlgorithmParams> getAlgorithmParams(){return alg_params;}
    
    std::string getContinFileName(){return continfile;}
    void setContinFileName(std::string cname){continfile = cname;}
    
    bool getIDinData(){return id_included;}
    void setIDinData(bool id){id_included = id;}
    
    float getContinMiss(){return contin_miss;}
    void setContinMiss(float val){contin_miss = val;}
    
    bool getOttEncoded(){return ott_dummy_encoded;}
    void setOttEncoded(bool ott){ott_dummy_encoded = ott;}
    
    int getNumExchanges(){return num_exchanges;}
    void setNumExchanges(int numEx, std::string id=""){
      if(numEx < 0)
        throw AthenaExcept("Number of exchanges " + id + " must be greater than or equal to zero");
      num_exchanges = numEx;
    }
    
    void setCVOutput(bool cvOut){cv_out = cvOut;}
    bool getCVOutput(){return cv_out;}
    
    void setStatusAdjust(std::string stat_change){status_change=stat_change;}
    std::string getStatusAdjust(){return status_change;}
    
    void setIndOutput(bool io){inds_output = io;}
    bool getIndOutput(){return inds_output;}
    
    void setOutputAllNodesBest(bool val){all_nodes_out = val;}
    bool outputAllNodesBest(){return all_nodes_out;}
    
    std::string getTrainFile(){return train_file;}
    void setTrainFile(std::string filename){train_file = filename;}
    
    std::string getTestFile(){return test_file;}
    void setTestFile(std::string filename){test_file = filename;}
    
    void setContinTestFile(std::string filename){contin_test = filename;}
    std::string getContinTestFile(){return contin_test;}
    
    void setContinTrainFile(std::string filename){contin_train = filename;}
    std::string getContinTrainFile(){return contin_train;}
    
    void setBioFilterFile(std::string filename){biofilter_file = filename;}
    std::string getBioFilterFile(){return biofilter_file;}
    
    void setBioGeneFile(std::string filename){biogene_file = filename;}
    std::string getBioGeneFile(){return biogene_file;}
    
    void setBioArchiveFile(std::string filename){bioarchive_file = filename;}
    std::string getBioArchiveFile(){return bioarchive_file;}

    void setBioFileType(std::string filetype){biofilter_file_type = filetype;}
    std::string getBioFileType(){return biofilter_file_type;}

    /// throws an exception if parameters are in error
    void checkConfig();

    enum DataEncodeType{
      None,
      OttDummy,
      StephenDummy
    };
    
    enum SummaryType{
      True,
      False,
      Best
    };
    
    inline void setSummaryOnly(std::string val){
      std::map<std::string, SummaryType>::iterator iter = summary_map.find(val);
      
      if(iter != summary_map.end()){
        summary_only = iter->second;
      }
      else{
        throw AthenaExcept(val + " is not a valid parameter for summary type");
      }  
    }
    
    SummaryType getSummaryOnly(){return summary_only;}
    
    inline void setEncodeType(std::string encodeType){
      std::map<std::string, DataEncodeType>::iterator iter = data_encode_map.find(encodeType);
      
      if(iter != data_encode_map.end())
        encodeDataType = iter->second;
      else
        throw AthenaExcept(encodeType + " is not a valid parameter for data encoding");
      if(encodeDataType == OttDummy)
        setOttEncoded(true);
      else
        setOttEncoded(false);
    }
    
    inline DataEncodeType getEncodeType(){return encodeDataType;}
    
    
    inline void setLogType(std::string logType){
      std::map<std::string, LogType>::iterator iter = log_type_map.find(logType);
      if(iter != log_type_map.end())
        logTypeSelected = iter->second;
      else
        throw AthenaExcept(logType + " is not a valid parameter for log type selection");
    }
    
    inline LogType getLogType(){return logTypeSelected;}
    
    // add parallel code 
    #ifdef PARALLEL
      void sendConfig();
      void receiveConfig();
    #endif
    
private:
    
    void initialize();
    
    DataEncodeType encodeDataType;
    std::map<std::string, DataEncodeType> data_encode_map;
    std::map<std::string, SummaryType> summary_map;
    LogType logTypeSelected;
    std::map<std::string, LogType> log_type_map;
    
    std::string data_type, out_name, map_name, datafile, continfile, status_change, 
      train_file, test_file, contin_test, contin_train, biofilter_file, biofilter_file_type,
      bioarchive_file, biogene_file;
    int miss_value, ncv, rand_seed, num_exchanges;
    float contin_miss, stat_miss_value;
    bool id_included, ott_dummy_encoded, cv_out, inds_output, all_nodes_out;
    std::vector<AlgorithmParams> alg_params;
    SummaryType summary_only;
    
};


#endif	/* _CONFIG_H */

