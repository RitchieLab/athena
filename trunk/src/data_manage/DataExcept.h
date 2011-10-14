#ifndef DATAEXCEPT_H_
#define DATAEXCEPT_H_

#include <exception>
#include <string>

namespace data_manage
{

/// Exception for data
class DataExcept: public std::exception
{

public:
  DataExcept() throw();
  DataExcept(std::string message);
  ~DataExcept()throw(){};

  virtual const char* what() const throw();
private:
  std::string error;
};

}

#endif /*DATAEXCEPT_H_*/
