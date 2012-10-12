/*
Copyright Marylyn Ritchie 2011

This file is part of ATHENA.

ATHENA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ATHENA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ATHENA.  If not, see <http://www.gnu.org/licenses/>.
*/

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
