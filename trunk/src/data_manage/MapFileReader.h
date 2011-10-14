#ifndef MAPFILEREADER_H_
#define MAPFILEREADER_H_

#include "Dataholder.h"
#include "DataExcept.h"

namespace data_manage
{

class MapFileReader
{
public:

  /// parses the map file and stores the names in the Dataholder
  void parse_map_file(std::string mapfile, Dataholder* holder);

};

}

#endif /*MAPFILEREADER_H_*/
