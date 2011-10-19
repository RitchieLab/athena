#ifndef DEFINES_H
#define DEFINES_H

#include <cstdlib>
#include "config.h"
#include <string>
#include "Exception.h"
#include <limits>

//TODO: autoconfize the switches & version

#define ANNIE_VERSION "0.70"			
#define MAX_NEURONS_IN_LAYER 1000

#define _SIMPLE_NEURON_STRING "SimpleNeuron"
#define _ABSTRACT_NEURON_STRING "AbstractNeuron"
#define _PARAMETRIZED_NEURON_STRING "ParametrizedNeuron"
#define _GAB_NEURON_STRING "GABNeuron"

#define _CENTER_NEURON_STRING "CenterNeuron"
#define _RECURRENT_NEURON_STRING "RecurrentNeuron"
#define _INPUT_NEURON_STRING "InputNeuron"
#define _KOHONEN_NEURON_STRING "KohonenNeuron"

#define VM_ENABLED
namespace annie { void assert_failed(const char *cond, const char *function, const char *file, int line); }

#ifdef DEBUG
	#define ASSERT(cond)	do { if(!(cond)) assert_failed(#cond, __FUNCTION__, __FILE__, __LINE__); } while(0)
#else
	#define ASSERT(cond)	((void)0)
#endif

#ifdef assert
	#undef assert
#endif
#define assert ASSERT

/**
 * Paranoid check (mostly function parameter checks)
 * These are normally asserts, but...  asserts are typically removed on release. 
 * The library must be (at least a bit robust) - PARANOID gives us a chance to possibly to tune robustness OR speed
 *
 */
#define PARANOID(assert_expression, exception_text)	{ if(!(assert_expression))	throw Exception(exception_text);}

//OP
#define CONTROL
#define SIGMOID_APPROX	//should result in faster sigmoid (and derivation) function computation, is not yet throughly tested, hovewer
//

typedef unsigned int uint;	//it's also in standard C types, so we cannot have it inside annie - it causes trouble when "using namespace annie" ..

namespace annie
{
///Use this instead of double/float for real numbers pertaining to annie
typedef double real;
		
const real REAL_MAX=std::numeric_limits<double>::max(), REAL_MIN=std::numeric_limits<double>::min();


/** The TrainingSet can be saved as a binary file or a text file, the latter allowing
  * users to create a training set without using annie. A binary file is referred to
  * as annie::BINARY_FILE and a text file as annie::TEXT_FILE
  * @see TrainingSet
  * @see TrainingSet::save
  * @see TEXT_FILE
  */
const int BINARY_FILE = 0;

/** @see BINARY_FILE */
const int TEXT_FILE = 1;

/*
#ifdef SIGMOID_APPROX	
const real SIGMOID_APROX_THRESHOLD=30;	/// sigmoid(x) for x \not <- [-SIGMOID_APROX_THRESHOLD, SIGMOID_APROX_THRESHOLD] is set to 0 or 1 resp.
#endif		
*/
// smd changed to match current athena implementation
const real SIGMOID_APROX_THRESHOLD=709;

}

#include "Vector.h"
#include "Control.h"
#ifdef VM_ENABLED
	#define VM(a) if(!!defaultControl.get("verbose")) cerr << a;	//verbose MLP - only for heavy debugging, it's too slowwwww for std usage
#else
	#define VM(a) {}
#endif
#endif // define DEFINES_H

