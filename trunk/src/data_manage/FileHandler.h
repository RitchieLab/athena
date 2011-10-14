#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include "Dataholder.h"
#include "DataExcept.h"

namespace data_manage
{

/// Provides interface for multiple data files
class FileHandler
{
public:
	/// virtual destructor
	virtual ~FileHandler();
	
	/// provides interface for filling dataholder object with genotypes
	virtual void parse_file(std::string filename, Dataholder* holder, int missingValue,
           float statusMissingValue, bool contains_id=false)=0; 
	
	/// provides interface for filling dataholder object with genotypes from a testing
	/// and a training file
	virtual void parse_file(std::string train_file, std::string test_file, Dataholder* holder,
	  int missingValue, float statusMissingValue, bool contains_id=false)=0;
	
	/// provides interface for writing out dataholder object
	virtual void write_file(std::string filename, Dataholder* holder)=0;
	
};

}

#endif /*FILEHANDLER_H_*/
