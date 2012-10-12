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

#ifndef RANDOM_FUNC_H
#define RANDOM_FUNC_H

#include<string>

namespace data_manage
{

void rand_seed(unsigned int seed);

float rand_float();

unsigned int rand_uint();

void rand_load_state(std::string filename);

void rand_save_state(std::string filename);

}
#endif