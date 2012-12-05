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
/* 
 * File:   AthenaExcept.h
 * Author: dudeksm
 *
 * Created on November 6, 2008, 12:43 PM
 */

#ifndef _ATHENAEXCEPT_H
#define	_ATHENAEXCEPT_H

#include <exception>
#include <string>
using namespace std;

///
/// Class thrown for exception in PLATO system <br>
/// Error messages are set by the creating class
///

/// Exception for plato
class AthenaExcept: public exception{
        public:
          AthenaExcept() throw();
          AthenaExcept(string message){error=message;}
          ~AthenaExcept()throw(){};
          virtual const char* what() const throw(){return error.c_str();}

        private:
          string error;
};

#endif	/* _ATHENAEXCEPT_H */

