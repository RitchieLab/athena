/* 
 * File:   HemannExcept.h
 * Author: dudeksm
 *
 * Created on November 6, 2008, 12:43 PM
 */

#ifndef _HEMANNEXCEPT_H
#define	_HEMANNEXCEPT_H

#include <exception>
#include <string>
using namespace std;

///
/// Class thrown for exception in PLATO system <br>
/// Error messages are set by the creating class
///

/// Exception for plato
class HemannExcept: public exception{
        public:
          HemannExcept() throw();
          HemannExcept(string message){error=message;}
          ~HemannExcept()throw(){};
          virtual const char* what() const throw(){return error.c_str();}

        private:
          string error;
};

#endif	/* _HEMANNEXCEPT_H */

