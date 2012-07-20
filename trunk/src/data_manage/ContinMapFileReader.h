#ifndef CONTINMAPFILEREADER_H_
#define CONTINMAPFILEREADER_H_

#include "Dataholder.h"
#include "DataExcept.h"

namespace data_manage
{

class ContinMapFileReader
{
public:

  /// parses the map file and stores the names in the Dataholder
  void parse_map_file(std::string mapfile, Dataholder* holder);

};

}

#endif /*CONTINMAPFILEREADER_H_*/
