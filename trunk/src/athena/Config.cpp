#include "Config.h"

using namespace std;

///
/// Constructor
///
Config::Config(){
    initialize();
}


///
/// Initializes variables
///
void Config::initialize(){
    data_type = "TEXT";
    out_name = "athena";
    map_name = "";
    status_change = "";
    rand_seed = 1;
    datafile = "";
    id_included = false;
    continfile = "";
    contin_miss = -9999;
    ott_dummy_encoded = false;
    num_exchanges = 2;
    cv_out = false;
    status_change = "NONE";
    inds_output = false;
    all_nodes_out = false;
    test_file = "";
    train_file = "";
    contin_test = "";
    contin_train = "";
    miss_value = -1;
    biofilter_file = "";
    summary_only = False;
    stat_miss_value = -1.0;
    biofilter_file_type = "TEXT";
    bioarchive_file = biogene_file = "";
    encodeDataType = None;
    data_encode_map["FALSE"] = None;
    data_encode_map["NONE"] = None;
    data_encode_map["OTT"] = OttDummy;
    data_encode_map["STEPHEN"] = StephenDummy;
    data_encode_map["TRUE"] = OttDummy;
    log_type_map["NONE"] = LogNone;
    log_type_map["SUMMARY"] = LogSummary;
    log_type_map["DETAILED"] = LogDetailed;
    logTypeSelected = LogNone;
    summary_map["TRUE"] = True;
    summary_map["FALSE"] = False;
    summary_map["BEST"] = Best;
}

///
/// Checks configuration parameters for errors
/// @throws HemannExcept
///
void Config::checkConfig(){
  if(num_exchanges < 0)
    throw HemannExcept("NUMEXCHANGES must be greater than zero");
}


#ifdef PARALLEL

///
/// Broadcasts all parameters to slave nodes that the slave nodes require
/// for running.
///
void Config::sendConfig(){
  int nParams = 6; 
  // package each type and send to the other nodes  
  // only include those that the slaves need
  int * intvals = new int[nParams];
  
  intvals[0] = miss_value;
  intvals[1] = ncv;
  intvals[2] = rand_seed;
  intvals[3] = num_exchanges;
  intvals[4] = int(ott_dummy_encoded);
  intvals[5] = int(encodeDataType);
  
  MPI_Bcast(intvals, nParams, MPI_INT, 0, MPI_COMM_WORLD);
  
  delete [] intvals;
  
  // package floats -- again only pass those needed
  float cmiss = contin_miss;
  MPI_Bcast(&cmiss, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  // package up the algorithm parameters and pass
  // expand later to handle multiple algorithms
  // for now assume all get the same parameters
  
  // first send number of parameters for algorithm
  int num_alg_params = int(alg_params[0].params.size());
  
  MPI_Bcast(&num_alg_params, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  // size allocated to each string value
  int string_size = 200;
  
  // first string will be name of algorithm to use
  // each parameter is a map pair
  int total_char_size = num_alg_params * string_size * 2 + string_size;
  
  char * params = new char[total_char_size];
  
  int displacement = 0;
  
  uint i, j;
  
  for(i=0; i < alg_params[0].name.size(); i++){
    params[i] = alg_params[0].name[i];
  }
  // add null to end
  params[i] = '\0';
  
  std::map<std::string, std::string>::iterator mapiter;
  mapiter = alg_params[0].params.begin();
  
  for(i=0; i<alg_params[0].params.size(); i++){
    displacement += string_size;
    
    string key = mapiter->first;
    string value = mapiter->second;
    
    for(j=0; j<key.size(); j++)
      params[displacement+j] = key[j];
    // add null to end
    params[displacement + j] = '\0';
    
    displacement += string_size;
    for(j=0; j<value.size(); j++)
      params[displacement+j] = value[j];
    // add null to end
    params[displacement + j] = '\0';
    mapiter++;
    
  }
  
  MPI_Bcast(params, total_char_size, MPI_CHAR, 0, MPI_COMM_WORLD);
  delete []  params;

}


///
/// Receives broadcasted parameters from master node
///
void Config::receiveConfig(){

  // get integer values from master
  int nParams = 6;
  int * intvals = new int[nParams];
  MPI_Bcast(intvals, nParams, MPI_INT, 0, MPI_COMM_WORLD);
  
  miss_value = intvals[0];
  ncv = intvals[1];
  rand_seed = intvals[2];
  num_exchanges = intvals[3];
  ott_dummy_encoded = bool(intvals[4]);
  encodeDataType = static_cast<DataEncodeType>(intvals[5]);
  
  delete [] intvals;
  
  // get float values from master
  MPI_Bcast(&contin_miss, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  // next value should be the number of parameters for the algorithm
  int num_alg_params;
  MPI_Bcast(&num_alg_params, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // size allocated to each string
  int string_size = 200;
  
 // first string will be name of algorithm to use
  int total_char_size = num_alg_params * string_size *2 + string_size;
  
  char * params = new char[total_char_size];
  
  MPI_Bcast(params, total_char_size, MPI_CHAR, 0, MPI_COMM_WORLD);
  
  int displacement = 0;
  
  string algparam;
  
  int i,j;
  
  // get name first
  for(i=0; i<string_size; i++){
    // end when null is found
    if(params[i] == '\0')
      break;
    
    algparam += params[i];
    
  }
  
  AlgorithmParams currentAlg;
  currentAlg.name = algparam;

  algparam.clear();
 
  string key, value;
  
  for(i=0; i<num_alg_params; i++){
    displacement += string_size;
    key.clear();
    value.clear();
    for(j=0; j<string_size; j++){
      if(params[j+displacement] == '\0')
        break;
      key += params[j+displacement];
    }
    displacement += string_size;
    for(j=0; j<string_size; j++){
      if(params[j+displacement] == '\0')
        break;
      value += params[j+displacement];
    }
    currentAlg.params[key] = value;
  }
  
  alg_params.push_back(currentAlg);
  
  delete []  params;  
  
}

#endif
