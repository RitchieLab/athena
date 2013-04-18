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
#include "RandomNoReplace.h"

using namespace data_manage;


///
/// Algorithm 3.4.2S of Knuth's book Seminumeric Algorithm
///
void RandomNoReplace::sampleWithoutReplacement(int popSize, int sampSize, std::vector<int>& samples){
		// Use Knuth's variable names
		int& n = sampSize;
		int& N = popSize;

		int t = 0; // total input records dealt with
		int m = 0; // number of items selected so far
		double u;

		while (m < n){
				u = double(rand()) / RAND_MAX; // call a uniform(0,1) random number generator

				if ( (N - t)*u >= n - m ){
						t++;
				}
				else{
						samples[m] = t;
						t++; m++;
				}
		}
}
