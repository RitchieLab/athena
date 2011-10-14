#include "DataExcept.h"

using namespace std;

namespace data_manage
{

DataExcept::DataExcept(string message)
{
  error = message;
}

const char * DataExcept::what() const throw()
{
  return error.c_str();
}

}
