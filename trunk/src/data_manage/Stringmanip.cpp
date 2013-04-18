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
#include "Stringmanip.h"
#include <sstream>

using namespace std;

namespace data_manage
{

///
/// checks that the string passed in is a number
/// @param str string to check as number
/// @return true (1)  or false (0)
///
int Stringmanip::is_number(string str){
	for(unsigned int i=0; i<str.length(); i++){
		if((str[i]<'0'||str[i]>'9') && str[i] != '.' && str[i] != '-'){
			return 0;
		}
	}
	return 1;
}


///
/// Converts string to uppercase
/// @param convert string to convert to uppercase
/// @return converted string
///
string Stringmanip::to_upper(string convert){
	// convert keyword to uppercase
	unsigned int strpos, length = convert.length();
	for(strpos=0; strpos<length; strpos++){
		convert[strpos] = toupper(convert[strpos]);
	}
	return convert;
}

///
/// Checks whether the string passed is true or false or
/// can be ON for true result
/// @param value 
/// @return true or false
///
bool Stringmanip::check_true_false(string value){
		value = to_upper(value);
	 if(value.compare("TRUE") == 0 || value.compare("ON")==0)
			return true;
		return false;
		
}


}
