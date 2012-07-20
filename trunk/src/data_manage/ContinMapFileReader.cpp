#include "ContinMapFileReader.h"
#include <sstream>
#include <fstream>

namespace data_manage
{

///
/// Parses map file and stores loci names in the Dataholder
/// @param mapfile
/// @param dataholder Dataholder
///
void ContinMapFileReader::parse_map_file(string mapfile, Dataholder* dataholder){
  std::ifstream mapStream(mapfile.c_str(), ios::in);

  if(!mapStream.is_open()){
    throw DataExcept("ERROR: Unable to open " + mapfile + "\n");
  }

  string line;

  string continID;

  while(!mapStream.eof()){
    getline(mapStream, line);

    if(line.find_first_of("0123456789") == string::npos){
      continue;
    }

    stringstream ss(line);

    ss >> continID;
    dataholder->add_covar_name(continID);
  }

  mapStream.close();
}


}

