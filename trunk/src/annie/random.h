/**
 * Random generators
 * $Id: random.h,v 1.1 2004/06/16 10:53:30 opx Exp $
 * @author OP
 */

#ifndef RANDOM_H
#define RANDOM_H

#include "defines.h"

namespace annie
{
///returns a number drawn uniformly from min, max
real uniformRandom(real min, real max);

inline Vector uniformRandomVector(real min, real max, uint size)	{ Vector o(size); for(uint i=0; i<size; i++) o[i] = uniformRandom(min, max); return o; }

/// Generates a random real number between -1.0 and 1.0
real random();

/// Generates a random real number between 0.0 and 1.0
real random01();

/** Generates a random integer between given bounds.
  * @param low Lower bound
  * @param high Upper bound (exclusively)
  * @return A random integer i, low<=i<high
  */
int randomInt(int low, int high);

}

#endif
