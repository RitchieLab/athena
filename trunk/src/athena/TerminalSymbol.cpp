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
#include "TerminalSymbol.h"


///
/// Results are scaled from -1 to 1
/// @param x float with result
/// @return scaled result
///
float TerminalSymbol::ActivateSigmoid(float x)
{
				if(x < -709) return 0.0;
				if(x > 709)  return 1.0;
				return(1.0 / (1.0 + exp(-x)));
}

///
/// Adjusts the result if it is infinite 
/// or not a number
/// @param x float representing result
/// @return adjusted result 
///
float TerminalSymbol::AdjustResult(float x)
{

				if(isinf(x) == 1) {
								return 1.0;
				}
				if(isinf(x) == -1) {
								return -1.0;
				}
				if(isnan(x)) {
								return 0.0;
				}

				return x;
}

