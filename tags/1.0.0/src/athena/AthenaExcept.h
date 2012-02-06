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

