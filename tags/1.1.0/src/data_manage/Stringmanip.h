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
#ifndef STRINGMANIP_H_
#define STRINGMANIP_H_

#include <string>
#include <vector>
#include <sstream>

namespace data_manage
{

class Stringmanip
{
public:

		/// check if number contained in string
		static int is_number(std::string str);
		/// convert to uppercase
		static std::string to_upper(std::string convert);
		
		template <typename T> static std::string numberToString(T number){
			std::ostringstream ss;
			ss << number;
			return ss.str();
		}
		
		template <typename T> static T stringToNumber(std::string & text){
			std::istringstream ss(text);
			T result;
			return ss >> result ? result : 0;
		}

		static bool check_true_false(std::string value);
		
	static std::vector<std::string> &split(const std::string &s, char delim, 
		std::vector<std::string> &elems);

	static std::vector<std::string> split(const std::string &s, char delim);
		
};

}

#endif /*STRINGMANIP_H_*/
