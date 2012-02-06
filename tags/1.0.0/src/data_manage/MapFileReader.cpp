#include "MapFileReader.h"
#include <sstream>
#include <fstream>

namespace data_manage
{

///
/// Parses map file and stores loci names in the Dataholder
/// @param mapfile
/// @param dataholder Dataholder
///
void MapFileReader::parse_map_file(string mapfile, Dataholder* dataholder){
  std::ifstream mapStream(mapfile.c_str(), ios::in);

  if(!mapStream.is_open()){
    throw DataExcept("ERROR: Unable to open " + mapfile + "\n");
  }

  string line;

  string snpID;
  unsigned int chrom;
  unsigned int pos;

  while(!mapStream.eof()){
    getline(mapStream, line);

    if(line.find_first_of("0123456789") == string::npos){
      continue;
    }

    stringstream ss(line);

    ss >> chrom >> snpID >> pos;
    dataholder->add_geno_name(snpID);
  }

  mapStream.close();
}


}

