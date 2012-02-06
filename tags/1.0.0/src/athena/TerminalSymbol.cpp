#include "TerminalSymbol.h"


///
/// Results are scaled from -1 to 1
/// @param x float with result
/// @return scaled result
///
float TerminalSymbol::ActivateSigmoid(float x)
{
        if(x < -709) return -1.0;
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

