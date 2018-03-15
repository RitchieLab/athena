#include "defines.h"

namespace annie	{

///TODO: check these functions

real random()
{
	return (real)((rand()-RAND_MAX/2.0)/RAND_MAX);
}

real random01()
{
	real r = (unsigned int)rand();
	r = r/RAND_MAX;
	return r;
}

int randomInt(int low, int high)
{
	return rand()%(high-low)+low;
}

real uniformRandom(real min, real max)	{
	return min + (max-min)*random01();
}


}	//annie
